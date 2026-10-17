// test_buffer_tank_solver.cpp
//
// End-to-end buffer-tank tests that drive the real solvers (E-02).
//
// A two-step cycle: "Purge" feeds the top of the bed from a buffer tank that
// holds kPurgePercent % of the moles "Feed" vented through its top pressure
// boundary; then "Feed" pushes air up the bed and vents the top. Because the
// consumer (Purge) runs BEFORE the source (Feed), cycle 1 Purge uses its
// first-cycle assumption (a fixed velocity) and from cycle 2 on it must use
// the summary of Feed's run in the previous cycle.
//
// The unit tests in test_buffer_tank.cpp exercise the summary maths and the
// pointer wiring; these tests check that BOTH solvers (CoupledSolver and the
// chunked CoupledSolverChuncked used by the WASM build) actually refresh the
// summary at the end of each step, agree with each other, and ignore a stale
// summary left over from a previous run or loaded from a saved deck.

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/boundary_value.h"
#include "adsim/boundary_conditions/buffer_tank_source.h"
#include "adsim/boundary_conditions/flow_boundary_value.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/equations/isotherm.h"
#include "adsim/math/general_functions.h"
#include "adsim/objects/user_data.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/solver/log.h"

namespace adsim {

namespace {

constexpr int kNumCycles = 3;
constexpr double kStepTime = 20.0;          // [s] both steps
constexpr double kPurgePercent = 60.0;     // [%] of Feed's vented moles
constexpr double kFirstCyclePurgeVel = 0.01; // [m/s] first-cycle assumption
const std::string kFeedStep = "Feed";
const std::string kPurgeStep = "Purge";
const std::string kPurgeStream = "Purge In";
const std::string kFeedVentLogKey = "Feed Pressure Boundary Top";
constexpr int kPurgeIdx = 0; // step order: consumer first ...
constexpr int kFeedIdx = 1;  // ... then the buffer-tank source

std::unique_ptr<ConstantBoundaryValue> Constant(double v) {
  auto bv = std::make_unique<ConstantBoundaryValue>();
  bv->constant_value_ = v;
  return bv;
}

std::unique_ptr<VelocityBoundaryValue> Velocity(double v) {
  auto vel = std::make_unique<VelocityBoundaryValue>();
  static_cast<ConstantBoundaryValue*>(vel->profile_.get())->constant_value_ = v;
  return vel;
}

std::unique_ptr<PressureBoundaryCondition> PressureOutlet(BoundaryLocation loc) {
  auto pbc = std::make_unique<PressureBoundaryCondition>();
  pbc->direction_ = BoundaryDirection::OUTLET;
  pbc->location_ = loc;
  pbc->pressure_ = Constant(101325.0);
  pbc->mole_fraction_["Oxygen"] = 0.21;
  pbc->mole_fraction_["Nitrogen"] = 0.79;
  return pbc;
}

std::unique_ptr<FlowBoundaryCondition> FlowInlet(BoundaryLocation loc,
                                                 std::unique_ptr<BoundaryStream> stream) {
  auto fbc = std::make_unique<FlowBoundaryCondition>();
  fbc->direction_ = BoundaryDirection::INLET;
  fbc->location_ = loc;
  fbc->streams_.push_back(std::move(stream));
  return fbc;
}

// Builds a small, fast O2/N2 PSA-like deck with an active buffer tank on the
// Purge inlet sourced from Feed's top pressure boundary.
void BuildBufferTankDeck(UserData& data) {
  data.bed.wall_.thickness_ = 0.0175;
  data.bed.wall_.heat_capacity_ = 537;
  data.bed.wall_.density_ = 7800;
  data.bed.wall_.thermal_conductivity_ = 16.2;

  data.fluid.AddComponent("Oxygen");
  data.fluid.AddComponent("Nitrogen");

  auto ads_id = data.adsorbent_directory.Add("13X");
  auto layer_id = data.bed.layers_.Add("13X");
  auto& layer = data.bed.layers_.Get(layer_id);
  layer.momentum_model_ = MomentumModel::KarmenCozeny;  // fixture predates Ergun: keep its Darcy / Kozeny-Carman pressure drop
  layer.domain_.length = 1.0;
  layer.domain_.n_cells = 10;
  layer.particle_diameter_ = 0.0011;
  layer.eb_ = 0.37;
  layer.composition_[ads_id] = 1.0;

  auto& ads = data.adsorbent_directory.Get(ads_id);
  ads.ep_ = 0.37;
  auto langmuir = std::make_unique<LangmuirIsotherm>();
  langmuir->params_["Oxygen"].m1 = 0.5;
  langmuir->params_["Oxygen"].b1 = 3e-11;
  langmuir->params_["Oxygen"].q1 = 30000;
  langmuir->params_["Nitrogen"].m1 = 1.5;
  langmuir->params_["Nitrogen"].b1 = 3e-11;
  langmuir->params_["Nitrogen"].q1 = 35000;
  ads.adsorption_source_.AddIsotherm({"Oxygen", "Nitrogen"}, std::move(langmuir));
  for (const char* comp : {"Oxygen", "Nitrogen"}) {
    auto ldf = std::make_unique<LinearDrivingForce>();
    auto mtc = std::make_unique<MassTransferCoefficientConstant>();
    mtc->constant_mtc_ = 0.1;
    ldf->mtc_model_ = std::move(mtc);
    ads.adsorption_source_.mass_transfer_models_[comp] = std::move(ldf);
  }

  auto ihtc = std::make_unique<ConstantInternalHeatTransferModel>();
  ihtc->heat_transfer_coefficient_ = 0.0;
  layer.equations_.internal_heat_transfer_model_ = std::move(ihtc);
  auto ehtc = std::make_unique<ConstantExternalHeatTransferModel>();
  ehtc->heat_transfer_coefficient_ = 0.0;
  layer.equations_.external_heat_transfer_model_ = std::move(ehtc);

  data.ic.pressure = 101325.0;
  data.ic.temperature = 298.0;
  data.ic.wall_temperature = 298.0;
  data.ic.velocity = 0.0;
  data.ic.mole_fraction["Oxygen"] = 0.21;
  data.ic.mole_fraction["Nitrogen"] = 0.79;

  // ---- Step 1, Purge: buffer-tank inlet at the top, vent at the bottom ----
  data.cycle.AddStep(kPurgeStep);
  {
    Step& purge = data.cycle.GetStep(kPurgeStep);
    purge.step_time_ = kStepTime;
    auto purge_in = std::make_unique<BoundaryStream>();
    purge_in->name_ = kPurgeStream;
    purge_in->boundary_flow_ = Velocity(kFirstCyclePurgeVel); // first-cycle assumption
    purge_in->mole_fraction["Oxygen"] = 1.0;
    purge_in->mole_fraction["Nitrogen"] = 0.0;
    purge_in->temperature_ = 298.0;
    purge_in->buffer_tank_.Activate(kFeedStep, kPressureBoundaryEndpointName);
    purge_in->buffer_tank_.percentage_ = kPurgePercent;
    purge.top_boundary_condition_ = FlowInlet(BoundaryLocation::TOP, std::move(purge_in));
    purge.bottom_boundary_condition_ = PressureOutlet(BoundaryLocation::BOTTOM);
  }

  // ---- Step 2, Feed: air in at the bottom, vent at the top ----
  data.cycle.AddStep(kFeedStep);
  {
    Step& feed = data.cycle.GetStep(kFeedStep);
    feed.step_time_ = kStepTime;
    auto air = std::make_unique<BoundaryStream>();
    air->name_ = "Feed In";
    air->boundary_flow_ = Velocity(0.1);
    air->mole_fraction["Oxygen"] = 0.21;
    air->mole_fraction["Nitrogen"] = 0.79;
    air->temperature_ = 298.0;
    feed.bottom_boundary_condition_ = FlowInlet(BoundaryLocation::BOTTOM, std::move(air));
    feed.top_boundary_condition_ = PressureOutlet(BoundaryLocation::TOP);
  }

  data.cycle.number_of_cycles_ = kNumCycles;
  data.settings.max_step = 1.0;
}

// E-18: replaces Purge's direct-flow buffer-tank inlet with a VALVE inlet on
// the same buffer tank. The valve's flow stays pressure-driven (upstream
// pressure kPurgeValveUpstreamP); the tank supplies its temperature and
// composition from cycle 2 on.
constexpr double kPurgeValveUpstreamP = 111325.0; // [Pa]
constexpr double kPurgeValveCv = 0.03;
constexpr double kFirstCyclePurgeT = 310.0;       // [K] distinct from the Feed vent
void MakePurgeInletAValve(UserData& data) {
  Step& purge = data.cycle.GetStep(kPurgeStep);
  auto* fbc = dynamic_cast<FlowBoundaryCondition*>(purge.top_boundary_condition_.get());
  ASSERT_NE(fbc, nullptr);
  auto valve = std::make_unique<ValveBoundaryStream>();
  valve->name_ = kPurgeStream;
  valve->valve_.cv_ = kPurgeValveCv;
  valve->boundary_value_ = Constant(kPurgeValveUpstreamP);
  valve->mole_fraction["Oxygen"] = 1.0;    // first-cycle assumption
  valve->mole_fraction["Nitrogen"] = 0.0;
  valve->temperature_ = kFirstCyclePurgeT; // first-cycle assumption
  valve->buffer_tank_.Activate(kFeedStep, kPressureBoundaryEndpointName);
  valve->buffer_tank_.percentage_ = kPurgePercent;
  fbc->streams_[0] = std::move(valve);
}

CycleEndpointSummary* FeedVentSummary(UserData& data) {
  return data.cycle.GetStep(kFeedStep).FindSummaryForEndpoint(kPressureBoundaryEndpointName);
}

// Per-cycle aggregates read back from a solver's boundary log.
struct RunResult {
  std::vector<double> feed_vent_moles;   // [mol] signed moles out of Feed's top, per cycle
  std::vector<double> purge_in_moles;    // [mol] |moles| into Purge's top, per cycle
  std::vector<double> purge_in_flow_min; // [mol/s] min |molar flow| of Purge In, per cycle
  std::vector<double> purge_in_flow_max; // [mol/s] max |molar flow| of Purge In, per cycle
};

RunResult Extract(const BoundaryLog& log) {
  RunResult r;
  for (int c = 0; c < kNumCycles; ++c) {
    auto feed = log.GetVariableData(c, kFeedIdx, kFeedVentLogKey, BoundaryVariableType::MOLAR_FLOW_RATE);
    auto purge = log.GetVariableData(c, kPurgeIdx, kPurgeStream, BoundaryVariableType::MOLAR_FLOW_RATE);
    EXPECT_GE(feed.time.size(), 2u) << "no Feed vent data in cycle " << c;
    EXPECT_GE(purge.time.size(), 2u) << "no Purge inlet data in cycle " << c;
    r.feed_vent_moles.push_back(IntegrateVariableInTime(feed.variable_timeseries, feed.time));
    std::vector<double> abs_flow;
    for (double v : purge.variable_timeseries) abs_flow.push_back(std::abs(v));
    r.purge_in_moles.push_back(IntegrateVariableInTime(abs_flow, purge.time));
    double lo = abs_flow.empty() ? 0.0 : abs_flow.front();
    double hi = lo;
    for (double v : abs_flow) { lo = std::min(lo, v); hi = std::max(hi, v); }
    r.purge_in_flow_min.push_back(lo);
    r.purge_in_flow_max.push_back(hi);
  }
  return r;
}

RunResult RunNonChunked(UserData& data) {
  CoupledSolver solver{data};
  int retval = solver.IntegrateCycle(data.cycle, false);
  EXPECT_GE(retval, 0) << "CoupledSolver did not complete";
  return Extract(solver.GetBoundaryLog());
}

RunResult RunChunked(UserData& data) {
  CoupledSolverChuncked solver{data};
  solver.Start();
  ChunkStatus status;
  do {
    solver.Integrate(50);
    status = solver.GetStatus();
    if (status.GetStatus() == SolverStatus::FAILED) {
      ADD_FAILURE() << "Chunked solver failed: " << status.GetMessage();
      break;
    }
  } while (solver.IsActive());
  EXPECT_EQ(status.GetStatus(), SolverStatus::COMPLETED) << status.GetMessage();
  return Extract(solver.GetBoundaryLog());
}

void ExpectRel(double a, double b, double rel, const std::string& what) {
  const double scale = std::max({std::abs(a), std::abs(b), 1e-30});
  EXPECT_LE(std::abs(a - b) / scale, rel) << what << ": " << a << " vs " << b;
}

// From cycle 2 on, Purge's inlet must be the constant buffer-tank flow
// kPurgePercent % * (moles Feed vented in the previous cycle) / step time.
void ExpectBufferTankDrivesPurge(const RunResult& r, const std::string& label) {
  for (int c = 1; c < kNumCycles; ++c) {
    const double expected_flow =
        kPurgePercent / 100.0 * std::abs(r.feed_vent_moles[c - 1]) / kStepTime;
    ExpectRel(r.purge_in_flow_min[c], expected_flow, 1e-6,
              label + " Purge inlet min flow, cycle " + std::to_string(c + 1));
    ExpectRel(r.purge_in_flow_max[c], expected_flow, 1e-6,
              label + " Purge inlet max flow, cycle " + std::to_string(c + 1));
  }
}

} // namespace

// E-02: the chunked solver (used by the web build) must refresh the buffer-tank
// summary at the end of every step exactly like CoupledSolver does, so both
// solvers produce the same cycle-2+ purge flows.
TEST(BufferTankSolverTest, ChunkedSolverUpdatesSummaryAndMatchesNonChunked) {
  UserData data_ref;
  BuildBufferTankDeck(data_ref);
  std::string err;
  ASSERT_TRUE(data_ref.PerformChecks(err)) << err;
  const RunResult ref = RunNonChunked(data_ref);
  ASSERT_FALSE(HasFailure());

  UserData data_chunk;
  BuildBufferTankDeck(data_chunk);
  ASSERT_TRUE(data_chunk.PerformChecks(err)) << err;
  const RunResult chunk = RunChunked(data_chunk);
  ASSERT_FALSE(HasFailure());

  // The source summary is valid after the chunked run and reflects the last cycle.
  const CycleEndpointSummary* sum = FeedVentSummary(data_chunk);
  ASSERT_NE(sum, nullptr);
  EXPECT_TRUE(sum->is_valid) << "chunked solver never updated the buffer-tank summary";
  EXPECT_EQ(sum->source_cycle_index, kNumCycles - 1);

  // The buffer tank actually changed the purge flow (otherwise this test
  // could not tell a working tank from the first-cycle assumption).
  ASSERT_GT(std::abs(ref.purge_in_moles[1] - ref.purge_in_moles[0]),
            0.05 * ref.purge_in_moles[0])
      << "deck does not distinguish buffer-tank flow from the first-cycle assumption";

  ExpectBufferTankDrivesPurge(ref, "CoupledSolver");
  ExpectBufferTankDrivesPurge(chunk, "CoupledSolverChuncked");

  for (int c = 0; c < kNumCycles; ++c) {
    const std::string cyc = ", cycle " + std::to_string(c + 1);
    ExpectRel(chunk.feed_vent_moles[c], ref.feed_vent_moles[c], 1e-3, "Feed vent moles" + cyc);
    ExpectRel(chunk.purge_in_moles[c], ref.purge_in_moles[c], 1e-3, "Purge inlet moles" + cyc);
  }
}

// A summary marked valid before the run starts (left over from a previous run
// in the same session, or deserialised from a deck saved after a run) must be
// discarded: cycle 1 always uses the first-cycle assumption, in both solvers.
TEST(BufferTankSolverTest, StaleSummaryIsResetAtRunStart) {
  UserData data_fresh;
  BuildBufferTankDeck(data_fresh);
  std::string err;
  ASSERT_TRUE(data_fresh.PerformChecks(err)) << err;
  const RunResult fresh = RunChunked(data_fresh);
  ASSERT_FALSE(HasFailure());

  auto poison = [](UserData& data) {
    CycleEndpointSummary* sum = FeedVentSummary(data);
    ASSERT_NE(sum, nullptr);
    sum->is_valid = true;
    sum->total_moles = 1000.0; // absurd: 30 mol/s of purge if consumed
    sum->avg_temperature = 298.0;
    sum->avg_mole_fraction = {{"Oxygen", 1.0}, {"Nitrogen", 0.0}};
    sum->source_cycle_index = 99;
  };

  {
    UserData data;
    BuildBufferTankDeck(data);
    ASSERT_TRUE(data.PerformChecks(err)) << err;
    poison(data);
    const RunResult r = RunChunked(data);
    ASSERT_FALSE(HasFailure());
    ExpectRel(r.purge_in_moles[0], fresh.purge_in_moles[0], 1e-9,
              "chunked: cycle-1 purge moles with a stale summary");
    EXPECT_EQ(FeedVentSummary(data)->source_cycle_index, kNumCycles - 1);
  }
  {
    UserData data;
    BuildBufferTankDeck(data);
    ASSERT_TRUE(data.PerformChecks(err)) << err;
    poison(data);
    const RunResult r = RunNonChunked(data);
    ASSERT_FALSE(HasFailure());
    ExpectRel(r.purge_in_moles[0], fresh.purge_in_moles[0], 1e-3,
              "non-chunked: cycle-1 purge moles with a stale summary");
    EXPECT_EQ(FeedVentSummary(data)->source_cycle_index, kNumCycles - 1);
  }
}

// E-18: a VALVE inlet with an active buffer tank must take the tank's
// temperature and composition (the molar-flow-weighted averages of Feed's
// vent in the previous cycle) from cycle 2 on, not the first-cycle
// assumption forever.
TEST(BufferTankSolverTest, ValveInletUsesTankTemperatureAndCompositionFromCycle2) {
  UserData data;
  BuildBufferTankDeck(data);
  MakePurgeInletAValve(data);
  ASSERT_FALSE(HasFailure());
  std::string err;
  ASSERT_TRUE(data.PerformChecks(err)) << err;

  CoupledSolver solver{data};
  ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "CoupledSolver did not complete";
  const BoundaryLog& log = solver.GetBoundaryLog();

  for (int c = 0; c < kNumCycles; ++c) {
    const std::string cyc = "cycle " + std::to_string(c + 1);
    double expected_T = kFirstCyclePurgeT;
    double expected_y_o2 = 1.0;
    if (c > 0) {
      // Recompute the tank summary from Feed's vent in the previous cycle.
      const auto n = log.GetVariableData(c - 1, kFeedIdx, kFeedVentLogKey, BoundaryVariableType::MOLAR_FLOW_RATE);
      const auto T = log.GetVariableData(c - 1, kFeedIdx, kFeedVentLogKey, BoundaryVariableType::TEMPERATURE);
      const auto n_o2 = log.GetComponentVariableData(c - 1, kFeedIdx, kFeedVentLogKey, BoundaryComponentVariableType::MOLAR_FLOW, "Oxygen");
      const auto n_n2 = log.GetComponentVariableData(c - 1, kFeedIdx, kFeedVentLogKey, BoundaryComponentVariableType::MOLAR_FLOW, "Nitrogen");
      ASSERT_GE(n.time.size(), 2u) << cyc;
      std::vector<double> abs_n, abs_n_T;
      for (size_t k = 0; k < n.variable_timeseries.size(); ++k) {
        abs_n.push_back(std::abs(n.variable_timeseries[k]));
        abs_n_T.push_back(abs_n.back() * T.variable_timeseries[k]);
      }
      expected_T = IntegrateVariableInTime(abs_n_T, n.time) / IntegrateVariableInTime(abs_n, n.time);
      const double mol_o2 = IntegrateVariableInTime(n_o2.variable_timeseries, n_o2.time);
      const double mol_n2 = IntegrateVariableInTime(n_n2.variable_timeseries, n_n2.time);
      expected_y_o2 = mol_o2 / (mol_o2 + mol_n2);
      // The deck must make the tank distinguishable from the assumption.
      ASSERT_GT(std::abs(expected_y_o2 - 1.0), 0.05) << cyc;
      ASSERT_GT(std::abs(expected_T - kFirstCyclePurgeT), 1.0) << cyc;
    }

    const auto purge_n = log.GetVariableData(c, kPurgeIdx, kPurgeStream, BoundaryVariableType::MOLAR_FLOW_RATE);
    const auto purge_T = log.GetVariableData(c, kPurgeIdx, kPurgeStream, BoundaryVariableType::TEMPERATURE);
    const auto purge_y = log.GetComponentVariableData(c, kPurgeIdx, kPurgeStream, BoundaryComponentVariableType::MOLE_FRACTION, "Oxygen");
    ASSERT_GE(purge_T.time.size(), 2u) << cyc;
    ASSERT_EQ(purge_y.variable_timeseries.size(), purge_T.variable_timeseries.size()) << cyc;
    double max_flow = 0.0;
    for (size_t k = 0; k < purge_T.variable_timeseries.size(); ++k) {
      max_flow = std::max(max_flow, std::abs(purge_n.variable_timeseries[k]));
      EXPECT_NEAR(purge_T.variable_timeseries[k], expected_T, 1e-6)
          << "Purge valve inlet temperature, " << cyc << ", sample " << k;
      EXPECT_NEAR(purge_y.variable_timeseries[k], expected_y_o2, 1e-8)
          << "Purge valve inlet O2 fraction, " << cyc << ", sample " << k;
    }
    EXPECT_GT(max_flow, 1e-5) << "valve inlet carried no flow, " << cyc;
  }
}

} // namespace adsim
