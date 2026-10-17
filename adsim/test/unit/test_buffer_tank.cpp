// STL includes
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// adsim includes
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/buffer_tank_source.h"
#include "adsim/boundary_conditions/cycle_endpoint_summary.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/equations/equation.h"
#include "adsim/objects/compressor.h"
#include "adsim/objects/fluid.h"
#include "adsim/solver/cycle.h"

// Gtest includes
#include <gtest/gtest.h>


namespace adsim {

  namespace {

    // Builds a StreamData time series sampled at `n_samples` points uniformly
    // over [0, duration]. molar_flow / temperature / mole_fractions per
    // component are supplied as values constant in time (the simplest case).
    StreamData MakeFlatStream(double molar_flow,
                              double temperature,
                              const std::unordered_map<std::string, double>& mole_fraction,
                              double duration,
                              int n_samples = 5) {
      StreamData s;
      const double dt = duration / (n_samples - 1);
      for (int i = 0; i < n_samples; ++i) {
        const double t = i * dt;
        std::unordered_map<std::string, double> component_flow;
        for (const auto& [c, y] : mole_fraction) {
          component_flow[c] = molar_flow * y;
        }
        s.Log(molar_flow,
              molar_flow * temperature, // enthalpy_flow = n_dot * T (consistent with PressureBC::Log)
              temperature,
              101325.0,
              0.0,
              t,
              component_flow,
              mole_fraction);
      }
      return s;
    }

  } // anonymous namespace

  TEST(BufferTankSummaryTest, FlatStreamYieldsExactAverages) {
    std::unordered_map<std::string, double> y = {{"N2", 0.7}, {"O2", 0.3}};
    const double n_dot = 2.5;
    const double T = 305.0;
    const double duration = 4.0;

    StreamData s = MakeFlatStream(n_dot, T, y, duration);
    CycleEndpointSummary summary = ComputeSummary(s, /*cycle_index=*/0);

    ASSERT_TRUE(summary.is_valid);
    EXPECT_NEAR(summary.total_moles, n_dot * duration, 1e-9);
    EXPECT_NEAR(summary.avg_temperature, T, 1e-9);
    EXPECT_NEAR(summary.avg_mole_fraction.at("N2"), 0.7, 1e-9);
    EXPECT_NEAR(summary.avg_mole_fraction.at("O2"), 0.3, 1e-9);
    EXPECT_EQ(summary.source_cycle_index, 0);
  }

  TEST(BufferTankSummaryTest, EmptyStreamIsInvalid) {
    StreamData s; // no Log() calls
    CycleEndpointSummary summary = ComputeSummary(s, 0);
    EXPECT_FALSE(summary.is_valid);
    EXPECT_EQ(summary.total_moles, 0.0);
  }

  TEST(BufferTankSummaryTest, ZeroFlowLeavesTotalMolesZero) {
    std::unordered_map<std::string, double> y = {{"N2", 1.0}};
    StreamData s = MakeFlatStream(/*molar_flow=*/0.0, /*T=*/300.0, y, /*duration=*/2.0);
    CycleEndpointSummary summary = ComputeSummary(s, 0);

    EXPECT_TRUE(summary.is_valid);
    EXPECT_NEAR(summary.total_moles, 0.0, 1e-12);
  }

  TEST(BufferTankSourceTest, InactiveByDefault) {
    BufferTankSource bt;
    EXPECT_FALSE(bt.IsActive());
    EXPECT_FALSE(bt.IsReady());
  }

  TEST(BufferTankSourceTest, ActivateThenDeactivate) {
    BufferTankSource bt;
    bt.Activate("SourceStep", "Stream 1");
    EXPECT_TRUE(bt.IsActive());
    EXPECT_EQ(bt.source_step_name_, "SourceStep");
    EXPECT_EQ(bt.source_endpoint_name_, "Stream 1");

    bt.Deactivate();
    EXPECT_FALSE(bt.IsActive());
    EXPECT_EQ(bt.summary_, nullptr);
  }

  TEST(BufferTankSourceTest, MolarFlowScalesByPercentageAndStepTime) {
    CycleEndpointSummary summary;
    summary.is_valid = true;
    summary.total_moles = 10.0;
    summary.avg_temperature = 310.0;
    summary.avg_mole_fraction = {{"N2", 1.0}};

    BufferTankSource bt;
    bt.Activate("Src", "Stream 1");
    bt.summary_ = &summary;
    bt.percentage_ = 50.0;

    // 50% of 10 mol / 2 s = 2.5 mol/s
    EXPECT_NEAR(bt.GetMolarFlow(2.0), 2.5, 1e-12);
    EXPECT_NEAR(bt.GetTemperature(), 310.0, 1e-12);
    EXPECT_NEAR(bt.GetMoleFraction("N2"), 1.0, 1e-12);
  }

  TEST(BufferTankSourceTest, ZeroSummaryYieldsZeroFlow) {
    CycleEndpointSummary summary;
    summary.is_valid = true;
    summary.total_moles = 0.0;
    BufferTankSource bt;
    bt.Activate("Src", "Stream 1");
    bt.summary_ = &summary;
    EXPECT_NEAR(bt.GetMolarFlow(1.0), 0.0, 1e-15);
  }

  TEST(BufferTankSourceTest, PerformChecksRejectsOutOfRangePercentage) {
    BufferTankSource bt;
    bt.Activate("Src", "Stream 1");
    bt.percentage_ = 110.0;
    std::string err;
    EXPECT_FALSE(bt.PerformChecks(err));
    EXPECT_FALSE(err.empty());
  }

  // Builds a two-step cycle: "Source" with a flow OUTLET and "Consumer"
  // with a flow INLET that pulls from Source's stream via BufferTankSource.
  // Verifies the wiring goes through Step::FindSummaryForEndpoint and
  // Cycle::ResolveBufferTankSources.
  TEST(BufferTankCycleTest, ResolveBindsSummaryPointer) {
    Cycle cycle;
    cycle.AddStep("Source");
    cycle.AddStep("Consumer");

    auto& source = cycle.GetStep("Source");
    auto& consumer = cycle.GetStep("Consumer");

    // Reach into Consumer's default bottom flow boundary (Stream 1) and
    // activate a buffer tank pointing at Source's bottom Stream 1.
    auto* consumer_flow = dynamic_cast<FlowBoundaryCondition*>(consumer.bottom_boundary_condition_.get());
    ASSERT_NE(consumer_flow, nullptr);
    ASSERT_FALSE(consumer_flow->streams_.empty());
    auto& consumer_stream = consumer_flow->streams_[0];

    auto* source_flow = dynamic_cast<FlowBoundaryCondition*>(source.bottom_boundary_condition_.get());
    ASSERT_NE(source_flow, nullptr);
    ASSERT_FALSE(source_flow->streams_.empty());
    auto& source_stream = source_flow->streams_[0];

    consumer_stream->buffer_tank_.Activate("Source", source_stream->name_);
    consumer_stream->buffer_tank_.percentage_ = 80.0;

    // Before resolution there's no pointer.
    EXPECT_EQ(consumer_stream->buffer_tank_.summary_, nullptr);

    cycle.ResolveBufferTankSources();

    EXPECT_EQ(consumer_stream->buffer_tank_.summary_, &source_stream->cycle_summary_);
    // Source has not produced any data yet, so summary is not yet valid.
    EXPECT_FALSE(consumer_stream->buffer_tank_.IsReady());
  }

  TEST(BufferTankCycleTest, UpdateCycleSummaryPopulatesAggregates) {
    Cycle cycle;
    cycle.AddStep("Source");
    auto& source = cycle.GetStep("Source");

    auto* source_flow = dynamic_cast<FlowBoundaryCondition*>(source.bottom_boundary_condition_.get());
    ASSERT_NE(source_flow, nullptr);
    auto& source_stream = source_flow->streams_[0];

    // Fabricate a stream_data_ as if integration had filled it.
    source_stream->stream_data_ = MakeFlatStream(
      /*molar_flow=*/1.5,
      /*T=*/295.0,
      /*mole_fraction=*/{{"N2", 1.0}},
      /*duration=*/3.0);

    source.UpdateCycleSummary(/*cycle_index=*/2);

    const auto* sum = source.FindSummaryForEndpoint(source_stream->name_);
    ASSERT_NE(sum, nullptr);
    EXPECT_TRUE(sum->is_valid);
    EXPECT_NEAR(sum->total_moles, 1.5 * 3.0, 1e-9);
    EXPECT_NEAR(sum->avg_temperature, 295.0, 1e-9);
    EXPECT_EQ(sum->source_cycle_index, 2);
  }

  TEST(BufferTankCycleTest, PerformChecksRejectsUnknownSourceStep) {
    Cycle cycle;
    cycle.AddStep("Consumer");
    auto& consumer = cycle.GetStep("Consumer");

    auto* consumer_flow = dynamic_cast<FlowBoundaryCondition*>(consumer.bottom_boundary_condition_.get());
    ASSERT_NE(consumer_flow, nullptr);
    consumer_flow->streams_[0]->buffer_tank_.Activate("NoSuchStep", "Stream 1");
    // Provide a mole fraction so the unrelated checks pass.
    consumer_flow->streams_[0]->mole_fraction["N2"] = 1.0;
    // Make the top BC an OUTLET so we don't trip the "both INLET" check.
    consumer.top_boundary_condition_->direction_ = BoundaryDirection::OUTLET;
    consumer.bottom_boundary_condition_->direction_ = BoundaryDirection::INLET;

    Fluid fluid;
    fluid.components_ = {"N2"};

    std::string err;
    EXPECT_FALSE(cycle.PerformChecks(fluid, err));
    EXPECT_NE(err.find("unknown source step"), std::string::npos)
        << "Expected unknown-source-step error, got: " << err;
  }

  // --------------------------------------------------------------------------
  // E-18: every inlet stream kind takes its temperature and composition from
  // an active buffer tank (not only the direct-flow BoundaryStream), and a
  // pressure boundary logs the same inlet source that its residuals use.
  // --------------------------------------------------------------------------
  namespace {

    constexpr int kCells = 6; // includes the two ghost cells

    // Uniform O2/N2 bed state with everything the stream flow functions,
    // PressureBoundaryCondition::Apply and the Log() functions read.
    Equations MakeBedState(double pressure, double temperature, double y_o2) {
      Equations state;
      state.pressure.var.assign(kCells, pressure);
      state.pressure.res_var.assign(kCells, 0.0);
      state.fluid_temperature.var.assign(kCells, temperature);
      state.fluid_temperature.res_var.assign(kCells, 0.0);
      state.concentration.var.assign(kCells, pressure / (8.314 * temperature));
      state.concentration.res_var.assign(kCells, 0.0);
      state.heat_capacity.var.assign(kCells, 1000.0);
      state.molar_mass.var.assign(kCells, 0.029);
      state.molar_flow_rate.var.assign(kCells + 1, 0.0);
      const std::vector<std::pair<std::string, double>> gas = {{"O2", y_o2}, {"N2", 1.0 - y_o2}};
      for (const auto& [comp, y_val] : gas) {
        Equation y{ "Mole Fraction", CellType::CellCentre, EquationType::Algebraic, nullptr };
        y.var.assign(kCells, y_val);
        y.res_var.assign(kCells, 0.0);
        state.mole_fraction.emplace(comp, std::move(y));
        IndependentEquation c{ "Molar Concentration", CellType::CellCentre, nullptr };
        c.var.assign(kCells, y_val * state.concentration.var[0]);
        state.molar_concentration.emplace(comp, std::move(c));
      }
      return state;
    }

    // A source summary whose temperature / composition differ from every
    // first-cycle assumption used below. Starts not-yet-valid.
    CycleEndpointSummary MakeTankSummary() {
      CycleEndpointSummary summary;
      summary.is_valid = false;
      summary.total_moles = 5.0;
      summary.avg_temperature = 340.0;
      summary.avg_mole_fraction = {{"O2", 0.35}, {"N2", 0.65}};
      summary.source_cycle_index = 0;
      return summary;
    }

    const std::unordered_map<std::string, double> kAssumedY = {{"O2", 1.0}, {"N2", 0.0}};
    constexpr double kAssumedT = 293.0;

    // Inlet stream-side temperature / O2 fraction implied by the flow functions.
    struct InletView { double n, T, y_o2, y_n2; };
    InletView ViewInlet(IBoundaryStream& s, const Equations& state, BoundaryLocation loc) {
      const double n = s.GetMolarFlow(state, 2.0, loc);
      return { n,
               s.GetEnthalpyFlow(state, 2.0, loc) / n,
               s.GetComponentMolarFlow(state, 2.0, loc, "O2") / n,
               s.GetComponentMolarFlow(state, 2.0, loc, "N2") / n };
    }

  } // anonymous namespace

  TEST(BufferTankStreamKindTest, ValveInletTakesTankTemperatureAndComposition) {
    ValveBoundaryStream valve;
    valve.SetDirection(BoundaryDirection::INLET);
    valve.temperature_ = kAssumedT;
    valve.mole_fraction = kAssumedY;
    valve.valve_.cv_ = 0.2;
    const double P_up = 1.5e5;
    valve.boundary_value_ = std::make_unique<ConstantBoundaryValue>(P_up);

    CycleEndpointSummary summary = MakeTankSummary();
    valve.buffer_tank_.Activate("Src", "Stream 1");
    valve.buffer_tank_.summary_ = &summary;

    const double P_bed = 1.0e5;
    const Equations state = MakeBedState(P_bed, 300.0, 0.21);

    for (BoundaryLocation loc : {BoundaryLocation::BOTTOM, BoundaryLocation::TOP}) {
      const double sign = (loc == BoundaryLocation::BOTTOM) ? 1.0 : -1.0;

      // Source not run yet: first-cycle assumption.
      summary.is_valid = false;
      const double n_assumed = valve.valve_.GetMolarFlowRate(P_up, P_bed, 0.029, kAssumedT, 1000.0);
      ASSERT_GT(n_assumed, 0.0);
      InletView v = ViewInlet(valve, state, loc);
      EXPECT_NEAR(v.n, sign * n_assumed, 1e-12);
      EXPECT_NEAR(v.T, kAssumedT, 1e-9);
      EXPECT_NEAR(v.y_o2, 1.0, 1e-12);
      EXPECT_NEAR(v.y_n2, 0.0, 1e-12);

      // Source summary ready: the tank's temperature and composition. The
      // flow stays pressure-driven (the tank's moles do not set it).
      summary.is_valid = true;
      v = ViewInlet(valve, state, loc);
      EXPECT_NEAR(v.n, sign * valve.valve_.GetMolarFlowRate(P_up, P_bed, 0.029, 340.0, 1000.0), 1e-12);
      EXPECT_NEAR(v.T, 340.0, 1e-9);
      EXPECT_NEAR(v.y_o2, 0.35, 1e-12);
      EXPECT_NEAR(v.y_n2, 0.65, 1e-12);

      // The logged inlet stream reports the same values.
      valve.stream_data_.Clear();
      valve.Log(state, 2.0, loc);
      EXPECT_NEAR(valve.stream_data_.temperature_boundary_.back(), 340.0, 1e-9);
      EXPECT_NEAR(valve.stream_data_.mole_fraction_boundary_.at("O2").back(), 0.35, 1e-12);
    }
  }

  TEST(BufferTankStreamKindTest, PumpInletCompositionAndSuctionTemperatureFollowTank) {
    PumpBoundaryStream pump;
    pump.SetDirection(BoundaryDirection::INLET);
    pump.temperature_ = kAssumedT;
    pump.mole_fraction = kAssumedY;
    auto* power = dynamic_cast<ConstantPumpPowerModel*>(pump.power_.get());
    ASSERT_NE(power, nullptr);
    power->power_ = 2000.0;
    auto* comp = dynamic_cast<Compressor*>(pump.pump_.get());
    ASSERT_NE(comp, nullptr);
    comp->suction_pressure_ = 1e5;

    CycleEndpointSummary summary = MakeTankSummary();
    pump.buffer_tank_.Activate("Src", "Stream 1");
    pump.buffer_tank_.summary_ = &summary;

    const double P_bed = 3e5;
    const Equations state = MakeBedState(P_bed, 320.0, 0.21);

    for (BoundaryLocation loc : {BoundaryLocation::BOTTOM, BoundaryLocation::TOP}) {
      const double sign = (loc == BoundaryLocation::BOTTOM) ? 1.0 : -1.0;

      summary.is_valid = false;
      const double n_assumed = comp->CalculateMolarFlow(P_bed, 2000.0, kAssumedT);
      ASSERT_GT(n_assumed, 0.0);
      InletView v = ViewInlet(pump, state, loc);
      EXPECT_NEAR(v.n, sign * n_assumed, 1e-10 * n_assumed);
      EXPECT_NEAR(v.T, kAssumedT, 1e-9); // isothermal (default) compressor
      EXPECT_NEAR(v.y_o2, 1.0, 1e-12);

      // Ready: composition AND suction temperature come from the tank, so the
      // molar flow and the enthalpy flow are evaluated at the same T.
      summary.is_valid = true;
      const double n_tank = comp->CalculateMolarFlow(P_bed, 2000.0, 340.0);
      ASSERT_GT(std::abs(n_tank - n_assumed), 1e-6 * n_assumed);
      v = ViewInlet(pump, state, loc);
      EXPECT_NEAR(v.n, sign * n_tank, 1e-10 * n_tank);
      EXPECT_NEAR(v.T, 340.0, 1e-9);
      EXPECT_NEAR(v.y_o2, 0.35, 1e-12);
      EXPECT_NEAR(v.y_n2, 0.65, 1e-12);
    }
  }

  // Apply() and Log() must read the inlet from the same source. With a tank
  // that is active but not yet ready, Log() used to skip the first-cycle
  // assumption and log the connected stream instead. (Cycle::PerformChecks
  // rejects a tank together with a connection, but nothing stops the
  // boundary itself from being in that state.)
  TEST(BufferTankStreamKindTest, PressureInletLogMatchesApplyWhileTankNotReady) {
    PressureBoundaryCondition bc;
    bc.location_ = BoundaryLocation::BOTTOM;
    bc.direction_ = BoundaryDirection::INLET;
    bc.temperature_ = kAssumedT;
    bc.mole_fraction_ = kAssumedY;
    bc.pressure_ = std::make_unique<ConstantBoundaryValue>(1.0e5);

    CycleEndpointSummary summary = MakeTankSummary();
    bc.buffer_tank_.Activate("Src", "Stream 1");
    bc.buffer_tank_.summary_ = &summary;

    BoundaryStream source;
    for (double t : {0.0, 5.0, 10.0}) {
      source.stream_data_.Log(0.1, 0.1 * 360.0, 360.0, 1.0e5, 0.0, t,
                              {{"O2", 0.05}, {"N2", 0.05}}, {{"O2", 0.5}, {"N2", 0.5}});
    }
    bc.connection_.MakeConnection(1, source);

    for (bool ready : {false, true}) {
      summary.is_valid = ready;
      Equations state = MakeBedState(1.0e5, 300.0, 0.21);
      state.molar_flow_rate.var.assign(kCells + 1, 0.1);
      bc.Apply(state, 2.0);
      const double T_apply = 0.5 * (state.fluid_temperature.var[0] + state.fluid_temperature.var[1]
                                    - state.fluid_temperature.res_var[0]);
      const double y_apply = state.mole_fraction.at("O2").var[0] - state.mole_fraction.at("O2").res_var[0];

      EXPECT_NEAR(T_apply, ready ? 340.0 : kAssumedT, 1e-9) << "ready=" << ready;
      EXPECT_NEAR(y_apply, ready ? 0.35 : 1.0, 1e-12) << "ready=" << ready;

      bc.stream_data_.Clear();
      bc.Log(state, 2.0);
      EXPECT_NEAR(bc.stream_data_.temperature_boundary_.back(), T_apply, 1e-9) << "ready=" << ready;
      EXPECT_NEAR(bc.stream_data_.mole_fraction_boundary_.at("O2").back(), y_apply, 1e-12) << "ready=" << ready;
    }
  }

} // end adsim namespace
