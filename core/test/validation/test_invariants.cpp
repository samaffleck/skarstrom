/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Invariants: properties the solution must have everywhere, always.
//
// Unlike the analytical tests these do not need a closed form. They assert
// things that are true of the exact solution for structural reasons -- a
// maximum principle, non-negativity, closure of a composition -- and that a
// well-behaved discretisation must inherit. Their value is coverage: they run
// over every cell and every logged instant rather than one scalar functional,
// so they catch localised misbehaviour (an overshoot at a front, a negative
// loading in a nearly-empty cell) that an integral comparison averages away.

#include "breakthrough_case.h"

namespace adsim {

  using namespace adsim::validation;

  // ===========================================================================
  // 1. Discrete maximum principle for temperature
  // ===========================================================================
  //
  // Convection and diffusion of a scalar with no source obey a maximum
  // principle: the solution can never leave the envelope set by its initial and
  // boundary values. Feeding hot gas into a cold bed, no cell may ever exceed
  // the feed temperature or fall below the initial bed temperature.
  //
  // Upwind differencing preserves this; a central or higher-order scheme without
  // a limiter does not, and the symptom is a small overshoot right at the front.
  // Because that overshoot is local and transient, an integral energy balance
  // will not see it -- this test will.

  TEST(Invariants, TemperatureObeysDiscreteMaximumPrinciple) {
    for (bool heating : {true, false}) {
      BreakthroughCase c;
      c.adsorbing = false;   // no heat of adsorption: a pure thermal front
      // The thermal wave is ~1000x slower than the gas, so the run has to be
      // long for the front to reach the middle of the column.
      c.run_time = 20000.0;
      c.max_step = 20.0;
      c.bed_temperature = heating ? 300.0 : 360.0;
      c.feed_temperature = heating ? 360.0 : 300.0;

      UserData data;
      Build(data, c);

      CoupledSolver solver{data};
      std::string error;
      ASSERT_TRUE(data.PerformChecks(error)) << error;
      ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

      const double lower = std::min(c.bed_temperature, c.feed_temperature);
      const double upper = std::max(c.bed_temperature, c.feed_temperature);
      // The threshold is set by the integrator, not by the scheme: the
      // temperature equation carries an absolute tolerance of 1e-4 K, so
      // anything at that level is IDA's error control rather than a spatial
      // overshoot. A non-monotone discretisation overshoots by a few percent of
      // the 60 K jump -- four orders of magnitude above this.
      const double tolerance = 1e-3;  // [K]

      const auto T = solver.GetLog().GetVariable(0, 0, VariableType::TEMPERATURE);
      ASSERT_FALSE(T.values.empty());

      double worst_overshoot = 0.0;
      double worst_undershoot = 0.0;
      for (size_t n = 0; n < T.values.size(); ++n) {
        for (size_t i = 1; i + 1 < T.values[n].size(); ++i) {
          worst_overshoot = std::max(worst_overshoot, T.values[n][i] - upper);
          worst_undershoot = std::max(worst_undershoot, lower - T.values[n][i]);
        }
      }

      EXPECT_LT(worst_overshoot, tolerance)
        << (heating ? "Heating" : "Cooling") << " front overshoots the envelope by "
        << worst_overshoot << " K";
      EXPECT_LT(worst_undershoot, tolerance)
        << (heating ? "Heating" : "Cooling") << " front undershoots the envelope by "
        << worst_undershoot << " K";

      // The front must actually have moved into the bed, or the bounds hold
      // trivially because nothing changed.
      const auto& final_profile = T.values.back();
      const double mid = final_profile[final_profile.size() / 2];
      EXPECT_GT(std::abs(mid - c.bed_temperature), 0.1 * std::abs(upper - lower))
        << "Thermal front did not develop; the maximum principle holds trivially";
    }
  }

  // ===========================================================================
  // 2. Composition and loading bounds
  // ===========================================================================
  //
  // Mole fractions must lie in [0, 1] and sum to one; adsorbed loadings must be
  // non-negative and must not exceed what the isotherm allows at the highest
  // partial pressure the cell ever sees. The existing IsMassConserved helper
  // checks the sum only; the individual bounds are what catch a component going
  // slightly negative in a cell the front has not reached yet, which is a
  // classic symptom of a species equation that is not positivity-preserving.

  TEST(Invariants, MoleFractionsAndLoadingsStayBounded) {
    BreakthroughCase c;
    c.adsorbing = true;
    c.run_time = 400.0;

    UserData data;
    Build(data, c);

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    const auto& log = solver.GetLog();
    const std::vector<std::string> components{"Nitrogen", "Oxygen"};

    // Mole fraction bounds, per component.
    std::vector<VariableLogEntry> y;
    for (const auto& component : components) {
      y.push_back(log.GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION,
                                           component));
      ASSERT_FALSE(y.back().values.empty()) << "No data logged for " << component;
    }

    double worst_negative = 0.0;
    double worst_above_one = 0.0;
    double worst_sum_error = 0.0;
    for (size_t n = 0; n < y[0].values.size(); ++n) {
      for (size_t i = 1; i + 1 < y[0].values[n].size(); ++i) {
        double sum = 0.0;
        for (size_t k = 0; k < components.size(); ++k) {
          const double value = y[k].values[n][i];
          worst_negative = std::max(worst_negative, -value);
          worst_above_one = std::max(worst_above_one, value - 1.0);
          sum += value;
        }
        worst_sum_error = std::max(worst_sum_error, std::abs(sum - 1.0));
      }
    }

    // As above, the thresholds sit just clear of the integrator's own tolerance
    // (the mole-fraction equations carry abstol 1e-8) and far below the
    // percent-level excursions a positivity-violating scheme would produce.
    EXPECT_LT(worst_negative, 1e-6) << "A mole fraction went negative by " << worst_negative;
    EXPECT_LT(worst_above_one, 1e-6) << "A mole fraction exceeded 1 by " << worst_above_one;
    EXPECT_LT(worst_sum_error, 1e-6) << "Mole fractions summed to 1 +/- " << worst_sum_error;

    // Loading bounds. The Henry isotherm is unbounded in pressure, so the
    // analytic ceiling here is its value at the feed partial pressure and the
    // coldest temperature the bed reaches.
    const auto q = log.GetComponentVariable(0, 0, ComponentVariableType::ADSORBED_CONCENTRATION,
                                            "Oxygen");
    const auto T = log.GetVariable(0, 0, VariableType::TEMPERATURE);
    const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
    ASSERT_FALSE(q.values.empty());

    double worst_negative_loading = 0.0;
    double worst_excess_loading = 0.0;
    double peak_loading = 0.0;
    for (size_t n = 0; n < q.values.size(); ++n) {
      for (size_t i = 1; i + 1 < q.values[n].size(); ++i) {
        const double loading = q.values[n][i];
        worst_negative_loading = std::max(worst_negative_loading, -loading);
        peak_loading = std::max(peak_loading, loading);

        // Equilibrium loading at this cell's feed-composition ceiling.
        const double q_max = kCaseHenryK * P.values[n][i] * c.y_feed *
                             std::exp(kCaseHenryQ / (kR * T.values[n][i]));
        worst_excess_loading = std::max(worst_excess_loading, loading - q_max);
      }
    }

    EXPECT_LT(worst_negative_loading, 1e-12)
      << "Adsorbed loading went negative by " << worst_negative_loading << " mol/kg";
    EXPECT_LT(worst_excess_loading, 1e-9)
      << "Adsorbed loading exceeded the equilibrium value at the feed partial "
         "pressure by " << worst_excess_loading << " mol/kg";
    EXPECT_GT(peak_loading, 0.0) << "Nothing adsorbed; the bounds hold trivially";
  }

  // ===========================================================================
  // 3. Bounds hold through a multi-step pressure swing
  // ===========================================================================
  //
  // The single-step test above never reverses the flow. A pressurisation /
  // feed / blowdown sequence does, and reversal is where boundedness is most
  // likely to break: the upwind branch flips, the boundary conditions swap
  // roles, and loadings are driven back down towards zero.

  TEST(Invariants, BoundsHoldThroughAPressureSwingCycle) {
    UserData data;
    data.bed.wall_.thickness_ = 0.0175;
    data.bed.wall_.heat_capacity_ = 537;
    data.bed.wall_.density_ = 7800;
    data.bed.wall_.thermal_conductivity_ = 16.2;

    data.fluid.AddComponent("Nitrogen");
    data.fluid.AddComponent("Oxygen");

    auto ads = data.adsorbent_directory.Add("13X");
    auto lid = data.bed.layers_.Add("13X");
    auto& layer = data.bed.layers_.Get(lid);
    layer.momentum_model_ = MomentumModel::KarmenCozeny;  // fixture predates Ergun: keep its Darcy / Kozeny-Carman pressure drop
    layer.domain_.length = 1.0;
    layer.domain_.n_cells = 40;
    layer.eb_ = 0.4;
    layer.particle_diameter_ = 0.002;
    layer.composition_[ads] = 1.0;
    std::get<CylinderBed>(data.bed.shape_).diameter = 0.1;

    auto& adsorbent = data.adsorbent_directory.Get(ads);
    adsorbent.ep_ = 0.35;
    adsorbent.solid_density_ = 1130.0;
    adsorbent.heat_capacity_ = 1070.0;

    // Both components adsorb, with different affinities: a real separation.
    auto isotherm = std::make_unique<LangmuirIsotherm>();
    isotherm->params_["Oxygen"].m1 = 3.0;
    isotherm->params_["Oxygen"].b1 = 3.0e-11;
    isotherm->params_["Oxygen"].q1 = 20000.0;
    isotherm->params_["Nitrogen"].m1 = 3.0;
    isotherm->params_["Nitrogen"].b1 = 5.0e-11;
    isotherm->params_["Nitrogen"].q1 = 24000.0;
    adsorbent.adsorption_source_.AddIsotherm({"Oxygen", "Nitrogen"}, std::move(isotherm));
    for (const auto& component : {"Oxygen", "Nitrogen"}) {
      auto ldf = std::make_unique<LinearDrivingForce>();
      auto mtc = std::make_unique<MassTransferCoefficientConstant>();
      mtc->constant_mtc_ = 0.5;
      ldf->mtc_model_ = std::move(mtc);
      adsorbent.adsorption_source_.mass_transfer_models_[component] = std::move(ldf);
      adsorbent.adsorption_source_.heat_of_adsorption_[component] =
        (std::string(component) == "Nitrogen") ? 24000.0 : 20000.0;
    }

    auto ihtc = std::make_unique<ConstantInternalHeatTransferModel>();
    ihtc->heat_transfer_coefficient_ = 20.0;
    layer.equations_.internal_heat_transfer_model_ = std::move(ihtc);
    auto ehtc = std::make_unique<ConstantExternalHeatTransferModel>();
    ehtc->heat_transfer_coefficient_ = 10.0;
    layer.equations_.external_heat_transfer_model_ = std::move(ehtc);

    data.bed.ambient_temperature_ = 298.0;
    data.ic.pressure = 1.0e5;
    data.ic.temperature = 298.0;
    data.ic.wall_temperature = 298.0;
    data.ic.velocity = 0.0;
    data.ic.mole_fraction["Nitrogen"] = 0.79;
    data.ic.mole_fraction["Oxygen"] = 0.21;

    auto make_flow_inlet = [&](double molar_flow, BoundaryLocation location) {
      auto stream = std::make_unique<BoundaryStream>();
      auto flow = std::make_unique<MolarFlowBoundaryValue>();
      static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = molar_flow;
      stream->boundary_flow_ = std::move(flow);
      stream->mole_fraction["Nitrogen"] = 0.79;
      stream->mole_fraction["Oxygen"] = 0.21;
      stream->temperature_ = 298.0;
      auto bc = std::make_unique<FlowBoundaryCondition>();
      bc->direction_ = BoundaryDirection::INLET;
      bc->location_ = location;
      bc->streams_.push_back(std::move(stream));
      return bc;
    };
    auto make_pressure_outlet = [&](double pressure, BoundaryLocation location) {
      auto bc = std::make_unique<PressureBoundaryCondition>();
      bc->direction_ = BoundaryDirection::OUTLET;
      bc->location_ = location;
      auto value = std::make_unique<ConstantBoundaryValue>();
      value->constant_value_ = pressure;
      bc->pressure_ = std::move(value);
      bc->mole_fraction_["Nitrogen"] = 0.79;
      bc->mole_fraction_["Oxygen"] = 0.21;
      return bc;
    };
    auto make_closed = [&](BoundaryLocation location) {
      auto bc = std::make_unique<ClosedBoundaryCondition>();
      bc->location_ = location;
      return bc;
    };

    // Pressurise (closed top), feed (open top), blow down (reverse flow out the
    // bottom).
    data.cycle.AddStep("Pressurise");
    Step& pressurise = data.cycle.GetStep("Pressurise");
    pressurise.step_time_ = 30.0;
    pressurise.bottom_boundary_condition_ = make_flow_inlet(0.02, BoundaryLocation::BOTTOM);
    pressurise.top_boundary_condition_ = make_closed(BoundaryLocation::TOP);

    data.cycle.AddStep("Feed");
    Step& feed = data.cycle.GetStep("Feed");
    feed.step_time_ = 120.0;
    feed.bottom_boundary_condition_ = make_flow_inlet(0.02, BoundaryLocation::BOTTOM);
    feed.top_boundary_condition_ = make_pressure_outlet(3.0e5, BoundaryLocation::TOP);

    data.cycle.AddStep("Blowdown");
    Step& blowdown = data.cycle.GetStep("Blowdown");
    blowdown.step_time_ = 120.0;
    blowdown.bottom_boundary_condition_ = make_pressure_outlet(1.0e5, BoundaryLocation::BOTTOM);
    blowdown.top_boundary_condition_ = make_closed(BoundaryLocation::TOP);

    data.cycle.number_of_cycles_ = 3;
    data.settings.reltol = 1e-8;
    data.settings.max_step = 1.0;

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    const auto& log = solver.GetLog();
    const std::vector<std::string> components{"Nitrogen", "Oxygen"};
    const double saturation = 3.0;  // m1, shared by both components

    double worst_negative = 0.0;
    double worst_above_one = 0.0;
    double worst_sum_error = 0.0;
    double worst_negative_loading = 0.0;
    double worst_excess_loading = 0.0;
    bool flow_reversed = false;

    for (int cycle : log.sys_data.cycle_index_set_) {
      for (int step : log.sys_data.step_index_set_) {
        std::vector<VariableLogEntry> y, q;
        for (const auto& component : components) {
          y.push_back(log.GetComponentVariable(cycle, step, ComponentVariableType::MOLE_FRACTION,
                                               component));
          q.push_back(log.GetComponentVariable(
            cycle, step, ComponentVariableType::ADSORBED_CONCENTRATION, component));
        }
        if (y[0].values.empty()) continue;

        for (size_t n = 0; n < y[0].values.size(); ++n) {
          for (size_t i = 1; i + 1 < y[0].values[n].size(); ++i) {
            double sum = 0.0;
            for (size_t k = 0; k < components.size(); ++k) {
              const double value = y[k].values[n][i];
              worst_negative = std::max(worst_negative, -value);
              worst_above_one = std::max(worst_above_one, value - 1.0);
              sum += value;

              const double loading = q[k].values[n][i];
              worst_negative_loading = std::max(worst_negative_loading, -loading);
              worst_excess_loading = std::max(worst_excess_loading, loading - saturation);
            }
            worst_sum_error = std::max(worst_sum_error, std::abs(sum - 1.0));
          }
        }

        const auto u = log.GetFaceVariable(cycle, step, VariableType::VELOCITY);
        for (const auto& snapshot : u.values) {
          for (double velocity : snapshot) {
            if (velocity < -1e-6) flow_reversed = true;
          }
        }
      }
    }

    // The cycle must genuinely reverse, or this is just a longer version of the
    // single-step test.
    EXPECT_TRUE(flow_reversed) << "Flow never reversed; the cycle does not exercise "
                                  "the upwind branch switch";

    EXPECT_LT(worst_negative, 1e-6) << "A mole fraction went negative by " << worst_negative;
    EXPECT_LT(worst_above_one, 1e-6) << "A mole fraction exceeded 1 by " << worst_above_one;
    EXPECT_LT(worst_sum_error, 1e-4) << "Mole fractions summed to 1 +/- " << worst_sum_error;
    EXPECT_LT(worst_negative_loading, 1e-9)
      << "Adsorbed loading went negative by " << worst_negative_loading << " mol/kg";
    EXPECT_LT(worst_excess_loading, 1e-9)
      << "Adsorbed loading exceeded the Langmuir saturation capacity by "
      << worst_excess_loading << " mol/kg";
  }

} // namespace adsim
