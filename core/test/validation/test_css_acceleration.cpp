/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// CSS acceleration: does it actually reduce the number of cycles to reach the
// cyclic steady state?
//
// A PSA process reaches CSS only asymptotically: the end-of-cycle bed state maps
// to the next end-of-cycle state through a contraction whose fixed point is the
// CSS. Unaccelerated, the error decays geometrically, so hundreds of cycles can
// be needed for a slowly-contracting cycle. AccelerateCyclicSteadyState fits the
// slowest error mode from three successive start-of-cycle profiles and
// extrapolates the adsorbed-phase loading to the fixed point (Aitken Delta^2),
// re-seeding the solver so the *next* cycle starts much closer to CSS.
//
// This test pins the claim that the acceleration EARNS its keep: run the SAME
// cycle twice -- once with acceleration off, once on -- and require the
// accelerated run to reach the (independently established) CSS profile in
// strictly, and substantially, fewer cycles, while landing on the SAME CSS (the
// extrapolation must not bias the fixed point).
//
// The convergence metric is the L2 norm (over interior cells) of the difference
// between each end-of-cycle nitrogen adsorbed-loading profile and a well-converged
// reference CSS profile. "Cycles to CSS" is the first cycle whose distance falls
// below a fixed tolerance.

#include "breakthrough_case.h"

namespace adsim {

  using namespace adsim::validation;

  namespace {
    double InteriorL2Diff(const std::vector<double>& a, const std::vector<double>& b) {
      double s = 0.0;
      for (size_t i = 1; i + 1 < a.size() && i + 1 < b.size(); ++i) {
        const double d = a[i] - b[i];
        s += d * d;
      }
      return std::sqrt(s);
    }
    double InteriorL2Norm(const std::vector<double>& a) {
      double s = 0.0;
      for (size_t i = 1; i + 1 < a.size(); ++i) s += a[i] * a[i];
      return std::sqrt(s);
    }

    // Build a minimal single-adsorbent Skarstrom-like PSA cycle. Single adsorbent
    // is required: the acceleration only re-seeds adsorbed loadings for a
    // single-adsorbent bed (the logged aggregate loading maps back to the lone
    // adsorbent's loading).
    void BuildPsaCycle(UserData& data, double mtc = 0.5) {
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
      layer.domain_.n_cells = 20;
      layer.eb_ = 0.4;
      layer.particle_diameter_ = 0.002;
      layer.composition_[ads] = 1.0;
      std::get<CylinderBed>(data.bed.shape_).diameter = 0.1;

      auto& adsorbent = data.adsorbent_directory.Get(ads);
      adsorbent.ep_ = 0.35;
      adsorbent.solid_density_ = 1130.0;
      adsorbent.heat_capacity_ = 1070.0;

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
        auto mtc_model = std::make_unique<MassTransferCoefficientConstant>();
        mtc_model->constant_mtc_ = mtc;
        ldf->mtc_model_ = std::move(mtc_model);
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

      data.settings.reltol = 1e-8;
      data.settings.max_step = 2.0;
      data.settings.max_cycles_to_keep_ = 500;  // keep the whole trajectory
      // CSS break is handled by this test, not the solver; acceleration must run
      // every cycle regardless of the CSS method.
      data.settings.css.method_ = CSSMethod::NONE;
    }

    // Run the cycle for n_cycles and return the end-of-cycle (end of the last
    // step) O2 gas-phase loading profile for each completed cycle.
    std::vector<std::vector<double>> RunAndCollect(bool accelerate, int n_cycles, double mtc = 0.5) {
      UserData data;
      BuildPsaCycle(data, mtc);
      data.cycle.number_of_cycles_ = n_cycles;
      data.settings.acc_css.method_ =
        accelerate ? AccelerationMethod::AITKEN : AccelerationMethod::NONE;

      CoupledSolver solver{data};
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

      const auto& log = solver.GetLog();
      const int last_step = *log.sys_data.step_index_set_.rbegin();
      std::vector<std::vector<double>> end_of_cycle;
      for (int cycle_index : log.sys_data.cycle_index_set_) {
        // The adsorbed loading is the slowly-converging quantity that drives CSS
        // (and the one the acceleration re-seeds). Nitrogen is the stronger
        // adsorber, so its loading settles slowest.
        const auto c = log.GetComponentVariable(cycle_index, last_step,
            ComponentVariableType::ADSORBED_CONCENTRATION, "Nitrogen");
        if (c.values.empty()) continue;
        end_of_cycle.push_back(c.values.back());
      }
      return end_of_cycle;
    }

    // Replace the Pressurise step's flow inlet with a pressure inlet whose value
    // is a linear ramp: it starts at the bed's boundary pressure when the step is
    // initialised (InitialiseStartValue) and reaches end_pressure at the end of
    // the step (InitialiseStepTime). The value is a function of the time relative
    // to the step start, so it exercises step_start_time_.
    void UseLinearPressurisationRamp(UserData& data, double end_pressure) {
      Step& pressurise = data.cycle.GetStep("Pressurise");
      auto bc = std::make_unique<PressureBoundaryCondition>();
      bc->direction_ = BoundaryDirection::INLET;
      bc->location_ = BoundaryLocation::BOTTOM;
      bc->temperature_ = 298.0;
      bc->mole_fraction_["Nitrogen"] = 0.79;
      bc->mole_fraction_["Oxygen"] = 0.21;
      auto ramp = std::make_unique<LinearBoundaryValue>();
      ramp->end_value_ = end_pressure;
      bc->pressure_ = std::move(ramp);
      pressurise.bottom_boundary_condition_ = std::move(bc);
    }

    struct RampRun {
      // Bottom-face pressure at the first logged instant (t_rel = 0) of the
      // Pressurise step, one entry per cycle.
      std::vector<double> ramp_start_pressure;
      // Bottom-face pressure at the end of the cycle's last step (Blowdown).
      std::vector<double> blowdown_end_pressure;
      int applied = 0;
      int skipped = 0;
      bool completed = false;
    };

    template <typename Solver>
    RampRun CollectRampRun(const Solver& solver) {
      RampRun run;
      const auto& log = solver.GetLog();
      const int last_step = *log.sys_data.step_index_set_.rbegin();
      for (int cycle_index : log.sys_data.cycle_index_set_) {
        const auto p0 = log.GetFaceVariable(cycle_index, 0, VariableType::PRESSURE);
        const auto pN = log.GetFaceVariable(cycle_index, last_step, VariableType::PRESSURE);
        if (p0.values.empty() || pN.values.empty()) continue;
        run.ramp_start_pressure.push_back(p0.values.front().front());
        run.blowdown_end_pressure.push_back(pN.values.back().front());
      }
      return run;
    }

    // Run the PSA cycle with a linear pressurisation ramp as its first step and
    // record the ramp's starting boundary pressure per cycle plus the CSS
    // acceleration log messages. chunked selects CoupledSolverChuncked (the
    // GUI/web path) instead of CoupledSolver.
    RampRun RunRamp(bool accelerate, bool chunked, int n_cycles) {
      UserData data;
      BuildPsaCycle(data);
      UseLinearPressurisationRamp(data, 3.0e5);
      data.cycle.number_of_cycles_ = n_cycles;
      data.settings.acc_css.method_ =
        accelerate ? AccelerationMethod::AITKEN : AccelerationMethod::NONE;

      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;

      int applied = 0, skipped = 0;
      auto on_message = [&](LogLevel, const std::string& msg) {
        if (msg.find("CSS acceleration applied") != std::string::npos) ++applied;
        if (msg.find("CSS acceleration skipped") != std::string::npos) ++skipped;
      };

      RampRun run;
      if (chunked) {
        CoupledSolverChuncked solver{data};
        solver.SetMessageCallback(on_message);
        solver.SetLogLevel(LogLevel::MEDIUM);
        solver.Start();
        EXPECT_TRUE(solver.IsActive()) << solver.GetStatus().GetMessage();
        for (int chunks = 0; solver.IsActive() && chunks < 200000; ++chunks) {
          solver.Integrate(50);
        }
        run = CollectRampRun(solver);
        run.completed = !solver.IsActive()
          && static_cast<int>(run.ramp_start_pressure.size()) == n_cycles;
      } else {
        CoupledSolver solver{data};
        solver.SetMessageCallback(on_message);
        solver.SetLogLevel(LogLevel::MEDIUM);
        const int flag = solver.IntegrateCycle(data.cycle, false);
        EXPECT_GE(flag, 0) << "Solver did not complete";
        run = CollectRampRun(solver);
        run.completed = flag >= 0
          && static_cast<int>(run.ramp_start_pressure.size()) == n_cycles;
      }
      run.applied = applied;
      run.skipped = skipped;
      return run;
    }

    // E-14: a CSS re-seed must be validated and committed with the first step's
    // time-varying boundary values evaluated at t_rel = 0. With a stale
    // step_start_time_ (the start of the previous cycle's last step) the ramp is
    // evaluated past its end, so the consistent-IC solve pins the boundary to the
    // ramp's END value; the next cycle then starts its "ramp" at full pressure
    // (or the re-seed is rejected outright and acceleration never applies).
    void CheckReseedUsesRampStartValue(bool chunked) {
      const int n_cycles = 8;
      const RampRun base = RunRamp(/*accelerate=*/false, chunked, n_cycles);
      const RampRun acc  = RunRamp(/*accelerate=*/true,  chunked, n_cycles);
      ASSERT_TRUE(base.completed) << "Unaccelerated ramp run did not complete";
      ASSERT_TRUE(acc.completed) << "Accelerated ramp run did not complete";

      // The acceleration must actually be applied (the probe accepts the re-seed)
      // and never rejected for this well-behaved cycle.
      EXPECT_GT(acc.applied, 0) << "CSS acceleration was never applied";
      EXPECT_EQ(acc.skipped, 0) << "CSS re-seed rejected " << acc.skipped
                                << " time(s) (re-seed not integrable)";
      EXPECT_EQ(base.applied + base.skipped, 0);

      // GUARD: the blowdown must leave the bed well below the ramp's end value,
      // or the start-vs-end distinction below is vacuous.
      for (int k = 0; k < n_cycles; ++k) {
        ASSERT_LT(acc.blowdown_end_pressure[k], 1.2e5) << "cycle " << k;
      }

      // Every cycle after the first starts its pressurisation ramp at the bed's
      // end-of-blowdown pressure (the ramp's start value), not at its 3 bar end
      // value -- in the unaccelerated run and, crucially, in the cycles that
      // begin from a committed re-seed.
      for (int k = 1; k < n_cycles; ++k) {
        const double expected = acc.blowdown_end_pressure[k - 1];
        EXPECT_NEAR(acc.ramp_start_pressure[k], expected, 0.02e5)
          << "Accelerated run, cycle " << k << ": pressurisation ramp started at "
          << acc.ramp_start_pressure[k] << " Pa instead of the end-of-blowdown "
          << "pressure " << expected << " Pa";
        EXPECT_NEAR(base.ramp_start_pressure[k], base.blowdown_end_pressure[k - 1], 0.02e5)
          << "Unaccelerated run, cycle " << k;
        // And consistent with the unaccelerated run.
        EXPECT_NEAR(acc.ramp_start_pressure[k], base.ramp_start_pressure[k], 0.05e5)
          << "cycle " << k;
      }
    }
  } // namespace

  TEST(CSSAcceleration, ReseedEvaluatesFirstStepRampAtStepStart) {
    CheckReseedUsesRampStartValue(/*chunked=*/false);
  }

  TEST(CSSAcceleration, ReseedEvaluatesFirstStepRampAtStepStartChunked) {
    CheckReseedUsesRampStartValue(/*chunked=*/true);
  }

  // ===========================================================================
  // Acceleration reduces the number of cycles to CSS
  // ===========================================================================

  TEST(CSSAcceleration, ReducesCyclesToCSS) {
    // Establish the true CSS from a long, unaccelerated reference run. The
    // end-of-cycle map here contracts slowly (rho ~ 0.9), so ~90 cycles drives the
    // reference far below the convergence tolerance used later and lets it serve
    // as ground truth for BOTH runs.
    std::vector<std::vector<double>> reference = RunAndCollect(/*accelerate=*/false, 90);
    ASSERT_GE(reference.size(), 80u) << "Reference run did not complete";
    const std::vector<double>& css_ref = reference.back();
    const double scale = InteriorL2Norm(css_ref);
    ASSERT_GT(scale, 0.0);

    // The two runs under test: same cycle, acceleration off vs on.
    const int n_cycles = 55;
    std::vector<std::vector<double>> base = RunAndCollect(/*accelerate=*/false, n_cycles);
    std::vector<std::vector<double>> acc  = RunAndCollect(/*accelerate=*/true,  n_cycles);
    ASSERT_EQ(base.size(), static_cast<size_t>(n_cycles)) << "Unaccelerated run aborted early";
    ASSERT_EQ(acc.size(), static_cast<size_t>(n_cycles))
      << "Accelerated run aborted early -- acceleration must never fail the run";

    // GUARD: the cycle must genuinely evolve the bed, or "reaching CSS" is
    // vacuous. The first end-of-cycle profile must sit clearly away from CSS.
    const double move = InteriorL2Diff(base.front(), css_ref);
    ASSERT_GT(move, 0.04 * scale)
      << "First end-of-cycle profile is already at CSS; the cycle is not "
         "exercising the bed";

    // "Cycles to CSS": first cycle whose end-of-cycle loading is within tol of the
    // reference CSS. tol is far above the reference's own residual error, so the
    // reference is legitimate ground truth for both runs.
    const double tol = 0.005 * scale;
    auto cycles_to_css = [&](const std::vector<std::vector<double>>& eoc) -> int {
      for (size_t k = 0; k < eoc.size(); ++k) {
        if (InteriorL2Diff(eoc[k], css_ref) <= tol) return static_cast<int>(k);
      }
      return -1;
    };

    const int cycles_base = cycles_to_css(base);
    const int cycles_acc  = cycles_to_css(acc);

    ASSERT_GE(cycles_base, 0) << "Unaccelerated run never reached CSS tolerance";
    ASSERT_GE(cycles_acc, 0) << "Accelerated run never reached CSS tolerance";

    // 1. Acceleration must land on the SAME cyclic steady state: the final
    //    accelerated end-of-cycle profile must match the reference CSS. This
    //    guards against the extrapolation biasing the fixed point (i.e. "fewer
    //    cycles" achieved by converging to the wrong state).
    const double final_gap = InteriorL2Diff(acc.back(), css_ref);
    EXPECT_LT(final_gap, tol)
      << "Accelerated run converged to a different state than the unaccelerated "
      << "CSS (gap " << final_gap << " vs profile norm " << scale << ")";

    // 2. The core claim: acceleration reaches CSS in strictly fewer cycles.
    EXPECT_LT(cycles_acc, cycles_base)
      << "CSS acceleration did not reduce the cycle count (accelerated reached "
      << "CSS at cycle " << cycles_acc << ", unaccelerated at cycle "
      << cycles_base << ")";

    // 3. And the reduction must be worthwhile, not a one-cycle rounding win: at
    //    least a 20% reduction in cycles to CSS.
    EXPECT_LE(cycles_acc, static_cast<int>(0.80 * cycles_base))
      << "CSS acceleration reduced the cycle count only marginally (accelerated "
      << cycles_acc << " vs unaccelerated " << cycles_base << ")";
  }

} // namespace adsim
