// test_cycle_state_roundtrip.cpp
//
// The cycle-state API (adsim/solver/cycle_state.h, Simulation::FinalState /
// SetInitialState / RunCycle):
//
//  * Round trip: N cycles, take the state, run one cycle from it on a fresh
//    Simulation of the same deck, and compare with N + 1 cycles run straight.
//    Checked on a deck with inter-bed connections (FiveStepProductPurge,
//    coarsened) and on a deck with a buffer tank whose consumer runs before its
//    source (so it reads the PREVIOUS cycle's summary). Dropping the histories
//    from the state must break the match on both, which shows they are needed.
//  * Warm start from a cyclic-steady-state state reaches CSS in a few cycles.
//  * A state that does not fit the deck is refused with a reason.
//
// Tolerance: the restarted cycle differs from the straight one only by the
// integrator's path (a fresh IDA start instead of the restart at the cycle
// boundary, a band probe at another state), so the end states agree to the
// integrator tolerance: max |dx| / max(RMS(field), abstol) over each field of
// the differential state is below kRoundTripTol (measured on Linux/GCC:
// 3e-12 buffer tank, 1.3e-5 Vysyaraju2024, 4e-7 FiveStepProductPurge; the
// bound leaves room for platform differences). Without the histories the
// first two miss by 0.23 and 0.25.

#include <algorithm>
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
#include "adsim/objects/user_data.h"
#include "adsim/simulation.h"
#include "adsim/solver/cycle_state.h"

namespace adsim {

namespace {

constexpr double kRoundTripTol = 1e-4;
// A state without its histories must miss by far more than the tolerance.
constexpr double kMissingHistoryMin = 1e-2;

// Largest relative difference between two states over the differential
// entries of every field: max |a - b| / max(RMS(b), floor) per field.
double StateDifference(const CycleState& a, const CycleState& b) {
  EXPECT_EQ(a.y.size(), b.y.size());
  if (a.y.size() != b.y.size()) return 1e300;
  double worst = 0.0;
  for (const auto& f : b.fields) {
    if (!f.differential) continue;
    double ss = 0.0;
    size_t n = 0;
    for (size_t g : f.index) { ss += b.y[g] * b.y[g]; ++n; }
    const double scale = std::max(std::sqrt(ss / std::max<size_t>(n, 1)), 1e-8);
    for (size_t g : f.index) worst = std::max(worst, std::abs(a.y[g] - b.y[g]) / scale);
  }
  return worst;
}

// ---- a buffer-tank deck (consumer before source) ---------------------------
// Purge feeds the top of the bed from a tank holding 60 % of the moles Feed
// vented through its top in its last run; Feed runs after Purge, so from cycle
// 2 on Purge replays the previous cycle's summary.

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

std::unique_ptr<FlowBoundaryCondition> FlowInlet(BoundaryLocation loc, std::unique_ptr<BoundaryStream> stream) {
  auto fbc = std::make_unique<FlowBoundaryCondition>();
  fbc->direction_ = BoundaryDirection::INLET;
  fbc->location_ = loc;
  fbc->streams_.push_back(std::move(stream));
  return fbc;
}

UserData BufferTankDeck() {
  UserData data;
  data.bed.wall_.thickness_ = 0.0175;
  data.bed.wall_.heat_capacity_ = 537;
  data.bed.wall_.density_ = 7800;
  data.bed.wall_.thermal_conductivity_ = 16.2;

  data.fluid.AddComponent("Oxygen");
  data.fluid.AddComponent("Nitrogen");

  auto ads_id = data.adsorbent_directory.Add("13X");
  auto layer_id = data.bed.layers_.Add("13X");
  auto& layer = data.bed.layers_.Get(layer_id);
  layer.momentum_model_ = MomentumModel::KarmenCozeny;
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

  data.cycle.AddStep("Purge");
  {
    Step& purge = data.cycle.GetStep("Purge");
    purge.step_time_ = 20.0;
    auto in = std::make_unique<BoundaryStream>();
    in->name_ = "Purge In";
    in->boundary_flow_ = Velocity(0.01);  // first-cycle assumption
    in->mole_fraction["Oxygen"] = 1.0;
    in->mole_fraction["Nitrogen"] = 0.0;
    in->temperature_ = 298.0;
    in->buffer_tank_.Activate("Feed", kPressureBoundaryEndpointName);
    in->buffer_tank_.percentage_ = 60.0;
    purge.top_boundary_condition_ = FlowInlet(BoundaryLocation::TOP, std::move(in));
    purge.bottom_boundary_condition_ = PressureOutlet(BoundaryLocation::BOTTOM);
  }
  data.cycle.AddStep("Feed");
  {
    Step& feed = data.cycle.GetStep("Feed");
    feed.step_time_ = 20.0;
    auto air = std::make_unique<BoundaryStream>();
    air->name_ = "Feed In";
    air->boundary_flow_ = Velocity(0.1);
    air->mole_fraction["Oxygen"] = 0.21;
    air->mole_fraction["Nitrogen"] = 0.79;
    air->temperature_ = 298.0;
    feed.bottom_boundary_condition_ = FlowInlet(BoundaryLocation::BOTTOM, std::move(air));
    feed.top_boundary_condition_ = PressureOutlet(BoundaryLocation::TOP);
  }
  data.cycle.EnsureIds();
  data.cycle.number_of_cycles_ = 3;
  data.settings.max_step = 1.0;
  return data;
}

Simulation BufferTankSimulation() {
  Simulation sim;
  sim.SetData(BufferTankDeck());
  return sim;
}

// An example deck with inter-bed connections, on a coarse grid with the CSS
// check off so every run does its cycles.
Simulation ExampleSimulation(const std::string& deck) {
  Simulation sim = Simulation::FromFile("../examples/decks/" + deck + ".json");
  for (auto& layer : sim.Data().bed.layers_.layers_) layer.domain_.n_cells = 20;
  sim.Data().settings.css.method_ = CSSMethod::NONE;
  return sim;
}

// check_history: also check that the state without its histories misses (on
// decks where a connection or tank reads a step that runs later in the cycle,
// so the previous cycle's history matters).
template <class Make>
void ExpectRoundTrip(Make make, int n, const char* label, bool check_history = true) {
  Simulation straight = make();
  ASSERT_TRUE(straight.Run({.cycles = n + 1}).Ok()) << label << ": " << straight.LastResult().message;
  const CycleState reference = straight.FinalState();
  EXPECT_EQ(reference.cycles, n + 1);

  Simulation first = make();
  ASSERT_TRUE(first.Run({.cycles = n}).Ok()) << label << ": " << first.LastResult().message;
  const CycleState x = first.FinalState();
  EXPECT_EQ(x.cycles, n);
  EXPECT_TRUE(x.HasHistory()) << label;

  Simulation restarted = make();
  const RunResult r = restarted.RunCycle(x);
  ASSERT_TRUE(r.Ok()) << label << ": " << r.message;
  EXPECT_EQ(r.cycles_run, 1);
  EXPECT_EQ(restarted.Data().cycle.number_of_cycles_, make().Data().cycle.number_of_cycles_)
      << "RunCycle must not change the deck's number of cycles";
  const CycleState gx = restarted.FinalState();
  EXPECT_EQ(gx.cycles, n + 1);
  const double diff = StateDifference(gx, reference);
  std::cout << "[ cycle state ] " << label << ": round trip difference " << diff << "\n";
  EXPECT_LT(diff, kRoundTripTol) << label;

  if (!check_history) return;
  // Without the histories the restart replays nothing (connections) or the
  // first-cycle assumption (buffer tank): the cycle differs clearly.
  CycleState bare = x;
  for (auto& e : bare.endpoints) { e.history.Clear(); e.summary.Clear(); }
  Simulation no_history = make();
  ASSERT_TRUE(no_history.RunCycle(bare).Ok()) << label;
  const double miss = StateDifference(no_history.FinalState(), reference);
  std::cout << "[ cycle state ] " << label << ": without histories " << miss << "\n";
  EXPECT_GT(miss, kMissingHistoryMin) << label;
}

}  // namespace

TEST(CycleStateRoundTrip, BufferTankDeck) {
  ExpectRoundTrip(BufferTankSimulation, 2, "buffer tank");
}

TEST(CycleStateRoundTrip, ConnectedBedsDeck) {
  // Vysyaraju2024: a pressure-equalisation step replays a partner step that
  // runs later in the cycle, i.e. the previous cycle's logged outlet.
  ExpectRoundTrip([] { return ExampleSimulation("Vysyaraju2024_PSA_exp1"); }, 2, "Vysyaraju2024_PSA_exp1");
  // Four connections; each reads a step that ran earlier in the same cycle.
  ExpectRoundTrip([] { return ExampleSimulation("FiveStepProductPurge"); }, 2, "FiveStepProductPurge",
                  /*check_history=*/false);
}

TEST(CycleStateRoundTrip, RepeatedRunCycleStepsTheCycleMap) {
  // RunCycle() without a state continues from the last final state, so three
  // calls are three cycles of one run.
  Simulation straight = BufferTankSimulation();
  ASSERT_TRUE(straight.Run({.cycles = 3}).Ok());
  Simulation stepped = BufferTankSimulation();
  for (int k = 0; k < 3; ++k) ASSERT_TRUE(stepped.RunCycle().Ok()) << "cycle " << k;
  EXPECT_EQ(stepped.FinalState().cycles, 3);
  EXPECT_LT(StateDifference(stepped.FinalState(), straight.FinalState()), kRoundTripTol);
}

TEST(CycleStateRoundTrip, StateLayoutAndDifferentialVector) {
  Simulation sim = BufferTankSimulation();
  ASSERT_TRUE(sim.RunCycle().Ok());
  CycleState x = sim.FinalState();
  EXPECT_EQ(x.components, (std::vector<std::string>{"Oxygen", "Nitrogen"}));
  ASSERT_EQ(x.layers.size(), 1u);
  EXPECT_EQ(x.layers[0].cells, 10);
  ASSERT_EQ(x.steps.size(), 2u);
  EXPECT_EQ(x.steps[0].name, "Purge");
  EXPECT_EQ(x.y.size(), x.yp.size());
  size_t entries = 0;
  for (const auto& f : x.fields) entries += f.index.size();
  EXPECT_EQ(entries, x.y.size());

  const auto* u = x.FindField("Velocity");
  ASSERT_NE(u, nullptr);
  EXPECT_FALSE(u->differential);
  EXPECT_NE(x.FindField("Pressure"), nullptr);
  const auto y_o2 = x.FieldValues("Mole Fraction", "Oxygen");
  EXPECT_FALSE(y_o2.empty());
  EXPECT_NE(x.FindField("Adsorbed Concentration", "Nitrogen"), nullptr);
  EXPECT_THROW(x.FieldValues("No Such Field"), std::out_of_range);

  // Differential round trip, and a size check.
  std::vector<double> d = x.Differential();
  ASSERT_FALSE(d.empty());
  EXPECT_LT(d.size(), x.y.size());  // velocity and the ghost cells are algebraic
  CycleState copy = x;
  for (double& v : d) v *= 1.0;
  copy.SetDifferential(d);
  EXPECT_EQ(copy.y, x.y);
  EXPECT_THROW(copy.SetDifferential(std::vector<double>(d.size() + 1, 0.0)), std::invalid_argument);
}

TEST(CycleStateRoundTrip, IncompatibleStatesAreRefusedWithAReason) {
  Simulation sim = BufferTankSimulation();
  ASSERT_TRUE(sim.Run({.cycles = 2}).Ok());
  const CycleState x = sim.FinalState();
  EXPECT_EQ(sim.CheckState(x), "");

  {  // a different grid
    Simulation other = BufferTankSimulation();
    other.Data().bed.layers_.layers_[0].domain_.n_cells = 12;
    const std::string why = other.CheckState(x);
    EXPECT_NE(why.find("cell counts"), std::string::npos) << why;
    EXPECT_THROW(other.SetInitialState(x), std::invalid_argument);
  }
  {  // a different component list
    Simulation other = BufferTankSimulation();
    other.Data().fluid.AddComponent("Argon");
    EXPECT_NE(other.CheckState(x).find("components"), std::string::npos);
  }
  {  // a changed step time on a deck whose buffer tank reads the histories
    Simulation other = BufferTankSimulation();
    other.Data().cycle.step_list_[1].step_time_ = 25.0;
    EXPECT_NE(other.CheckState(x).find("changed its time"), std::string::npos) << other.CheckState(x);
    // ... is fine when nothing reads a history.
    CycleState bare = x;
    for (auto& e : bare.endpoints) { e.history.Clear(); e.summary.Clear(); }
    EXPECT_EQ(other.CheckState(bare), "");
  }
  {  // a removed step
    Simulation other = BufferTankSimulation();
    other.Data().cycle.step_list_.pop_back();
    EXPECT_NE(other.CheckState(x).find("step structure"), std::string::npos);
  }
  {  // a different DAE layout (isothermal): found when the run starts
    Simulation other = BufferTankSimulation();
    other.Data().bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
    ASSERT_EQ(other.CheckState(x), "");
    const RunResult r = other.RunCycle(x);
    EXPECT_EQ(r.status, RunStatus::Failed);
    EXPECT_NE(r.message.find("state layout"), std::string::npos) << r.message;
    EXPECT_FALSE(other.HasFinalState());
  }
}

TEST(CycleStateRoundTrip, WarmStartFromCssReachesCssInAFewCycles) {
  // Buffer-tank deck with the Per-Component mass-balance CSS check (min 5
  // cycles, 3 consecutive passes): a cold run needs many cycles; a warm start
  // from its CSS state needs only the consecutive-cycles rule's minimum (the
  // first warm cycle has no previous cycle in the log to compare with).
  auto make = [] {
    Simulation sim = BufferTankSimulation();
    auto& css = sim.Data().settings.css;
    css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
    css.min_cycles_ = 5;
    css.consecutive_cycles_ = 3;
    sim.Data().cycle.number_of_cycles_ = 200;
    return sim;
  };
  Simulation cold = make();
  const RunResult rc = cold.Run();
  ASSERT_EQ(rc.status, RunStatus::ReachedCSS) << rc.message;

  Simulation warm = make();
  warm.SetInitialState(cold.FinalState());
  const RunResult rw = warm.Run();
  std::cout << "[ cycle state ] buffer tank: cold " << rc.cycles_run << " cycles to CSS, warm "
            << rw.cycles_run << "\n";
  ASSERT_EQ(rw.status, RunStatus::ReachedCSS) << rw.message;
  EXPECT_LE(rw.cycles_run, 4);
  EXPECT_LT(rw.cycles_run, rc.cycles_run);
  EXPECT_EQ(warm.FinalState().cycles, rc.cycles_run + rw.cycles_run);
  EXPECT_EQ(warm.Data().settings.css.min_cycles_, 5) << "the CSS floor must be restored after the run";
  EXPECT_TRUE(warm.HasInitialState());
  warm.ClearInitialState();
  EXPECT_FALSE(warm.HasInitialState());
}

}  // namespace adsim
