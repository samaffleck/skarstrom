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
// below a fixed tolerance. The claim is checked with each advection scheme
// (upwind, WENO, van Leer).
//
// ExtrapolatesGeometricSequences checks the extrapolation itself on synthetic
// end-of-cycle states that follow exact geometric sequences: the jump lands on
// the fixed point of a fast mode, is limited (not dropped) for a slow one, and
// is skipped, with the reason reported, for growing steps and for steps that
// are not one geometric mode. SlowGeometricModeReachesCssSooner checks a deck
// (Nguyen 2023 X, 20 cells) whose approach to CSS is one slow mode, where the
// fitted jump is far beyond the limit.

#include <cmath>
#include <cstdio>
#include <string>

#include "psa_cycle_case.h"
#include "adsim/simulation.h"

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

    // The shared PSA cycle on a coarse grid (a single adsorbent, as the
    // acceleration requires), logging the whole trajectory.
    void BuildAccelerationCase(UserData& data, double mtc = 0.5) {
      BuildPsaCycle(data, 20, mtc);
      data.settings.reltol = 1e-8;
      data.settings.max_step = 2.0;
      data.settings.max_cycles_to_keep_ = 500;  // keep the whole trajectory
      // CSS break is handled by this test, not the solver; acceleration must run
      // every cycle regardless of the CSS method.
      data.settings.css.method_ = CSSMethod::NONE;
    }

    struct AccelerationRun {
      // End-of-cycle (end of the last step) nitrogen loading profile, one per
      // completed cycle.
      std::vector<std::vector<double>> end_of_cycle;
      int applied = 0;  // committed re-seeds
      int skipped = 0;  // rejected windows and rejected re-seeds
      std::vector<std::string> skip_messages;
    };

    // Run the cycle for n_cycles with the given advection scheme and collect the
    // end-of-cycle loading profiles and the CSS acceleration log messages.
    AccelerationRun RunAndCollect(bool accelerate, int n_cycles,
                                  AdvectionScheme scheme = AdvectionScheme::Upwind,
                                  double mtc = 0.5) {
      UserData data;
      BuildAccelerationCase(data, mtc);
      data.settings.advection_scheme_ = scheme;
      data.cycle.number_of_cycles_ = n_cycles;
      data.settings.acc_css.method_ =
        accelerate ? AccelerationMethod::AITKEN : AccelerationMethod::NONE;

      AccelerationRun run;
      CoupledSolver solver{data};
      solver.SetMessageCallback([&](LogLevel, const std::string& msg) {
        if (msg.find("CSS acceleration applied") != std::string::npos) ++run.applied;
        if (msg.find("CSS acceleration skipped") != std::string::npos) {
          ++run.skipped;
          run.skip_messages.push_back(msg);
        }
      });
      solver.SetLogLevel(LogLevel::MEDIUM);
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

      const auto& log = solver.GetLog();
      const int last_step = *log.sys_data.step_index_set_.rbegin();
      for (int cycle_index : log.sys_data.cycle_index_set_) {
        // The adsorbed loading is the slowly-converging quantity that drives CSS
        // (and the one the acceleration re-seeds). Nitrogen is the stronger
        // adsorber, so its loading settles slowest.
        const auto c = log.GetComponentVariable(cycle_index, last_step,
            ComponentVariableType::ADSORBED_CONCENTRATION, "Nitrogen");
        if (c.values.empty()) continue;
        run.end_of_cycle.push_back(c.values.back());
      }
      return run;
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
      BuildAccelerationCase(data);
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

  // ===========================================================================
  // The extrapolation on exact geometric sequences
  // ===========================================================================
  //
  // Feed AccelerateCyclicSteadyState three synthetic end-of-cycle states and
  // check the jump against the closed form. For a single geometric mode,
  // x_k = x* + rho^k e, the steps are r_k = rho^(k-1) (rho - 1) e and the Aitken
  // coefficient is mu = rho / (rho - 1), so x2 - mu r2 = x* exactly. A jump
  // limited to m natural steps (mu = -m) lands at
  //     x2 + m r2 = x* + rho^2 e (1 - m (1 - rho) / rho).
  namespace {
    struct SyntheticField {
      Equation* eq;
      double base;  // x*, uniform over the cells
      double amp;   // amplitude of e relative to base, signed (mole fractions pair up)
    };

    // The extrapolable fields of the PSA case with a base value and an error
    // amplitude each: temperatures around 300 K, loadings around 1 mol/kg and
    // the two mole fractions with opposite errors, so they still sum to one.
    std::vector<SyntheticField> SyntheticFields(UserData& data) {
      Equations& state = data.state;
      std::vector<SyntheticField> fields;
      for (Equation* eq : state.GetEquations()) {
        if (eq->GetEquationType() != EquationType::Differential) continue;
        if (eq == &state.concentration || eq == &state.fluid_temperature) continue;
        double base = 1.0, amp = 0.01;
        if (eq == &state.solid_temperature || eq == &state.wall_temperature) base = 300.0;
        for (auto& [component, y] : state.mole_fraction) {
          if (eq != &y) continue;
          base = component == "Nitrogen" ? 0.79 : 0.21;
          amp = component == "Nitrogen" ? 0.01 / 0.79 : -0.01 / 0.21;
        }
        fields.push_back({eq, base, amp});
      }
      return fields;
    }

    // Shape of the error over the cells (the first sine mode on the cell
    // centres) and a second, orthogonal one (the second mode).
    double Mode(int k, size_t c, size_t n) {
      const double pi = std::acos(-1.0);
      return std::sin(pi * k * (static_cast<double>(c) + 0.5) / static_cast<double>(n));
    }
  } // namespace

  TEST(CSSAcceleration, ExtrapolatesGeometricSequences) {
    UserData data;
    BuildAccelerationCase(data);
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    data.Initialise();
    const std::vector<SyntheticField> fields = SyntheticFields(data);
    ASSERT_GE(fields.size(), 4u);  // two loadings, two mole fractions at least
    bool has_mole_fraction = false;
    for (const auto& f : fields) has_mole_fraction |= (f.base < 1.0);
    ASSERT_TRUE(has_mole_fraction);

    // x_k = x* + rho^k e1 + kick_k e2 in every field (e2 orthogonal to e1).
    auto set_sample = [&](double rho, int k, double kick) {
      for (const auto& f : fields) {
        const size_t n = f.eq->var.size();
        for (size_t c = 0; c < n; ++c) {
          f.eq->var[c] = f.base * (1.0 + f.amp * (std::pow(rho, k) * Mode(1, c, n) + kick * Mode(2, c, n)));
        }
      }
    };
    // Feeds samples k = 0, 1, 2 at cycles start, start+1, start+2.
    struct Outcome { bool moved; AccelerationNote note; };
    auto run = [&](double rho, int start, const double kicks[3]) {
      AccelerateCyclicSteadyState acc;
      acc.method_ = AccelerationMethod::AITKEN;
      Outcome out{false, {}};
      for (int k = 0; k < 3; ++k) {
        set_sample(rho, k, kicks[k]);
        out.note = AccelerationNote{};
        out.moved = acc.Accelerate(data.state, start + k, &out.note);
        if (k < 2) EXPECT_FALSE(out.moved) << "jumped before three samples";
      }
      return out;
    };
    // Largest relative deviation of the state from x* + factor rho^2 e1.
    auto deviation = [&](double rho, double factor) {
      double worst = 0.0;
      for (const auto& f : fields) {
        const size_t n = f.eq->var.size();
        for (size_t c = 0; c < n; ++c) {
          const double expected = f.base * (1.0 + f.amp * factor * rho * rho * Mode(1, c, n));
          worst = std::max(worst, std::abs(f.eq->var[c] - expected) / f.base);
        }
      }
      return worst;
    };
    const double no_kick[3] = {0.0, 0.0, 0.0};
    const double m = AccelerateCyclicSteadyState::kMaxStepMultiple;

    // 1. A fast mode (rho = 0.6, mu = -1.5): the jump lands on x* exactly.
    {
      const Outcome out = run(0.6, 2, no_kick);
      EXPECT_TRUE(out.moved) << out.note.skipped;
      EXPECT_TRUE(out.note.skipped.empty());
      EXPECT_LT(deviation(0.6, 0.0), 1e-12);
    }

    // 2. A slow mode (rho = 0.965, mu = -27.6, beyond the limit m = 20): the
    //    jump is limited to m natural steps and closes the fraction
    //    m (1 - rho) / rho (73 %) of the gap, undershooting x*. (The earlier
    //    safeguard skipped this window.)
    {
      const double rho = 0.965;
      ASSERT_GT(rho / (1.0 - rho), m);  // GUARD: the fit asks for more than the limit
      const Outcome out = run(rho, 2, no_kick);
      EXPECT_TRUE(out.moved) << out.note.skipped;
      EXPECT_NE(out.note.fit.find("limited"), std::string::npos) << out.note.fit;
      const double factor = 1.0 - m * (1.0 - rho) / rho;
      ASSERT_GT(factor, 0.0);
      ASSERT_LT(factor, 0.9);  // GUARD: the jump is a sizeable part of the gap
      EXPECT_LT(deviation(rho, factor), 1e-12);
    }

    // 3. Growing steps (rho = 1.02): nothing to extrapolate to; the state is
    //    left as it is and the reason is reported.
    {
      const double rho = 1.02;
      const Outcome out = run(rho, 2, no_kick);
      EXPECT_FALSE(out.moved);
      EXPECT_NE(out.note.skipped.find("not shrinking"), std::string::npos) << out.note.skipped;
      EXPECT_LT(deviation(rho, 1.0), 1e-15);
    }

    // 3b. An oscillating mode (rho = -0.5, mu = +1/3): the jump lands on x*
    //     exactly. A growing oscillation (rho = -1.2, mu = +0.55), which no
    //     contracting mode explains, is skipped and never limited or applied.
    //     (A positive mu of 1/2 or more always comes with |r2| >= |r1|, so the
    //     step-ratio test catches it.)
    {
      const Outcome out = run(-0.5, 2, no_kick);
      EXPECT_TRUE(out.moved) << out.note.skipped;
      EXPECT_EQ(out.note.fit.find("limited"), std::string::npos) << out.note.fit;
      EXPECT_LT(deviation(-0.5, 0.0), 1e-12);
    }
    {
      const double rho = -1.2;
      const Outcome out = run(rho, 2, no_kick);
      EXPECT_FALSE(out.moved);
      EXPECT_FALSE(out.note.skipped.empty());
      EXPECT_LT(deviation(rho, 1.0), 1e-15);
    }

    // 4. Steps that shrink slowly (rho = 0.99) while turning: a component
    //    along the second mode appears in the last sample only, 12 % of the
    //    step. The step ratio is sqrt(0.99^2 + 0.12^2) = 0.997 < 1 and the steps
    //    point the same way (single-mode rate 0.99), but the Aitken coefficient
    //    comes out POSITIVE, (0.12^2 - 0.0099) / (0.01^2 + 0.12^2) = +0.31:
    //    jumping would move the state back along its own trajectory. The window
    //    is skipped and reported.
    {
      const double rho = 0.99;
      const double step = rho - 1.0;  // r1 along e1, per unit e1
      const double kicks[3] = {0.0, 0.0, 0.12 * std::abs(step)};
      const Outcome out = run(rho, 2, kicks);
      EXPECT_FALSE(out.moved);
      EXPECT_NE(out.note.skipped.find("geometric"), std::string::npos)
        << out.note.skipped << " / " << out.note.fit;
    }
  }

  // The same cycle, with acceleration off and on, for one advection scheme.
  // The claim must hold for every scheme: the end-of-cycle map has the same
  // slow mode (rate ~0.73 early, ~0.966 late) whichever scheme computes the
  // fluxes. An earlier version of the safeguard dropped every window whose
  // fitted jump exceeded five natural steps; with WENO or van Leer the third
  // window landed just past that limit (mu = -5.01), no later window came near
  // it again (mu -> -28 on the slow mode), and acceleration stopped after two
  // jumps: the accelerated run reached CSS at the same cycle as the plain one.
  void CheckReducesCyclesToCSS(AdvectionScheme scheme) {
    // Establish the true CSS from a long, unaccelerated reference run. The
    // end-of-cycle map here contracts slowly (rho ~ 0.966 late in the run), so
    // ~90 cycles drives the reference well below the convergence tolerance used
    // later (its own residual is ~0.2 tol) and lets it serve as ground truth for
    // BOTH runs.
    const std::vector<std::vector<double>> reference =
      RunAndCollect(/*accelerate=*/false, 90, scheme).end_of_cycle;
    ASSERT_GE(reference.size(), 80u) << "Reference run did not complete";
    const std::vector<double>& css_ref = reference.back();
    const double scale = InteriorL2Norm(css_ref);
    ASSERT_GT(scale, 0.0);

    // The two runs under test: same cycle, acceleration off vs on.
    const int n_cycles = 55;
    const AccelerationRun base_run = RunAndCollect(/*accelerate=*/false, n_cycles, scheme);
    const AccelerationRun acc_run  = RunAndCollect(/*accelerate=*/true,  n_cycles, scheme);
    const auto& base = base_run.end_of_cycle;
    const auto& acc  = acc_run.end_of_cycle;
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
    std::printf("  cycles to CSS: %d unaccelerated, %d accelerated (%d jumps applied, %d skipped)\n",
                cycles_base, cycles_acc, acc_run.applied, acc_run.skipped);

    ASSERT_GE(cycles_base, 0) << "Unaccelerated run never reached CSS tolerance";
    ASSERT_GE(cycles_acc, 0) << "Accelerated run never reached CSS tolerance";

    // 1. Acceleration must land on the SAME cyclic steady state: the final
    //    accelerated end-of-cycle profile must match the reference CSS. This
    //    guards against the extrapolation biasing the fixed point (i.e. "fewer
    //    cycles" achieved by converging to the wrong state). And it must end at
    //    least as close to it as the plain run of the same length does.
    const double final_gap = InteriorL2Diff(acc.back(), css_ref);
    const double base_final_gap = InteriorL2Diff(base.back(), css_ref);
    EXPECT_LT(final_gap, tol)
      << "Accelerated run converged to a different state than the unaccelerated "
      << "CSS (gap " << final_gap << " vs profile norm " << scale << ")";
    EXPECT_LE(final_gap, base_final_gap)
      << "Accelerated run ended further from the CSS than the plain run (gap "
      << final_gap << " vs " << base_final_gap << ")";

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

    // 4. Acceleration keeps working for the whole run: this smooth cycle offers
    //    a window every three cycles from the fifth on (17 in 55 cycles), and
    //    none should be rejected. (The safeguard that stopped after two jumps
    //    applied 2 here with WENO and van Leer.)
    EXPECT_GE(acc_run.applied, 15) << "CSS acceleration stopped being applied";
    EXPECT_EQ(acc_run.skipped, 0) << "first skip: "
      << (acc_run.skip_messages.empty() ? std::string() : acc_run.skip_messages.front());
    EXPECT_EQ(base_run.applied + base_run.skipped, 0);
  }

  TEST(CSSAcceleration, ReducesCyclesToCSSUpwind) {
    CheckReducesCyclesToCSS(AdvectionScheme::Upwind);
  }

  TEST(CSSAcceleration, ReducesCyclesToCSSWeno) {
    CheckReducesCyclesToCSS(AdvectionScheme::WENO);
  }

  TEST(CSSAcceleration, ReducesCyclesToCSSVanLeer) {
    CheckReducesCyclesToCSS(AdvectionScheme::VanLeer);
  }

  // ===========================================================================
  // A deck whose approach to CSS is one slow geometric mode
  // ===========================================================================
  //
  // Nguyen 2023 experiment X on a coarse 20-cell grid (upwind): after ~100
  // cycles the end-of-cycle state change decays as one clean geometric mode of
  // rate ~0.985 (the CO2 heel building up), so the fitted coefficient is
  // mu ~ -65 to -170, far beyond any fixed cap. (On 80 cells the same deck
  // decays at 0.995.) The earlier safeguard dropped those windows and applied 3
  // jumps in 435 cycles: the accelerated run stopped at the same cycle as the
  // plain one. The deck is the cheapest one found that shows the regime.
  //
  // Checks: the accelerated run stops in at most half the plain run's cycles,
  // and on KPIs that agree with the plain run's within twice the error the
  // stopping rule leaves on a plain run, g eps rho / (1 - rho) (see
  // test_css_convergence.cpp), with rho and g measured on the plain run's last
  // cycles.
  namespace {
    struct DeckCssRun {
      RunResult result;
      std::vector<CSSRecord> css;   // the logged cycles
      std::vector<KPIRecord> kpis;  // every cycle
      int applied = 0;
    };

    DeckCssRun RunSlowModeDeck(bool accelerate) {
      Simulation sim = Simulation::FromFile("../examples/decks/Nguyen2023_Exp_X.json");
      sim.SetParameter("settings.advection_scheme", static_cast<int>(AdvectionScheme::Upwind));
      sim.SetParameter("layer[1].cells", 20);
      sim.SetParameter("settings.max_cycles_to_keep", 25);
      sim.SetParameter("settings.acceleration_method",
                       static_cast<int>(accelerate ? AccelerationMethod::AITKEN : AccelerationMethod::NONE));
      DeckCssRun run;
      RunOptions options;
      options.log_level = LogLevel::MEDIUM;
      options.on_message = [&](LogLevel, const std::string& msg) {
        if (msg.find("CSS acceleration applied") != std::string::npos) ++run.applied;
      };
      run.result = sim.Run(options);
      run.css = sim.CSSHistory();
      run.kpis = sim.KPIs();
      return run;
    }
  } // namespace

  TEST(CSSAcceleration, SlowGeometricModeReachesCssSooner) {
    const DeckCssRun plain = RunSlowModeDeck(false);
    const DeckCssRun acc = RunSlowModeDeck(true);
    ASSERT_TRUE(plain.result.css_reached) << plain.result.message;
    ASSERT_TRUE(acc.result.css_reached) << acc.result.message;
    ASSERT_GE(plain.css.size(), 21u);
    const double tol = plain.css.back().state_tolerance;
    ASSERT_GT(tol, 0.0);

    // The plain run's rate over its last ten logged cycles, and the KPI change
    // per unit state change over the same cycles.
    const size_t n = plain.css.size();
    const double rho = std::pow(plain.css[n - 1].state_change / plain.css[n - 11].state_change, 0.1);
    double sum_change = 0.0;
    for (size_t k = n - 10; k < n; ++k) sum_change += plain.css[k].state_change;
    const int last = plain.css[n - 1].cycle, first = plain.css[n - 11].cycle;
    auto kpi_at = [](const std::vector<KPIRecord>& kpis, int cycle) -> const KPIRecord& {
      for (const auto& r : kpis) if (r.cycle == cycle && r.valid) return r;
      ADD_FAILURE() << "no KPIs at cycle " << cycle;
      return kpis.front();
    };
    const KPIRecord& plain_end = kpi_at(plain.kpis, last);
    const KPIRecord& plain_first = kpi_at(plain.kpis, first);
    const KPIRecord& acc_end = kpi_at(acc.kpis, acc.result.css_cycle);

    std::printf("  plain: CSS at cycle %d (rate %.4f); accelerated: cycle %d, %d jumps\n",
                plain.result.css_cycle + 1, rho, acc.result.css_cycle + 1, acc.applied);

    // GUARD: the deck is in the slow-mode regime, and the plain run is long.
    ASSERT_GT(rho, 0.97) << "the approach to CSS is not slow; the deck no longer tests the limit";
    ASSERT_GE(plain.result.css_cycle, 200);

    EXPECT_LE(acc.result.css_cycle, plain.result.css_cycle / 2)
      << "acceleration did not halve the cycles to CSS";
    EXPECT_GE(acc.applied, 10);

    for (auto [name, field] : {std::pair{"purity", &KPIRecord::purity},
                               std::pair{"recovery", &KPIRecord::recovery}}) {
      const double g = std::abs(plain_end.*field - plain_first.*field) / sum_change;
      const double allowed = g * tol * rho / (1.0 - rho);
      const double diff = std::abs(acc_end.*field - plain_end.*field);
      std::printf("  %s: plain %.7f, accelerated %.7f, difference %.2e (plain stop error %.2e)\n",
                  name, plain_end.*field, acc_end.*field, diff, allowed);
      ASSERT_GT(allowed, 0.0);
      EXPECT_LT(diff, 2.0 * allowed) << name;
    }
  }

} // namespace adsim
