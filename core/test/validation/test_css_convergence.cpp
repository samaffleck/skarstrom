/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Convergence to cyclic steady state, measured on the KPIs.
//
// test_cyclic_steady_state checks that the end-of-cycle profile settles and
// test_css_acceleration that acceleration needs fewer cycles. Neither says how
// far the KPIs a run reports are from their CSS values when the stopping rule
// fires, nor checks the model behind the rule's tolerance. That model is the
// linearisation of the end-of-cycle map x_{k+1} = M(x_k) at its fixed point
// x*: x_{k+1} - x* = J (x_k - x*), so once the slowest mode (eigenvalue rho)
// dominates,
//
//   (1) the state change Delta_k = |x_k - x_{k-1}| decays geometrically,
//       Delta_{k+1} / Delta_k -> rho;
//   (2) every smooth functional of the state -- a KPI -- converges at the same
//       rate: f_k = f* + A rho^k, so f* can be extrapolated from three cycles
//       m apart (Aitken's delta-squared, exact for a single geometric mode) and
//       the KPI of any later cycle predicted, f_n = f* + (f_k - f*) rho^(n-k);
//   (3) the KPI change per cycle is proportional to the state change,
//       |f_{k+1} - f_k| = g Delta_k, so the error left when the rule stops at a
//       state change Delta_K <= eps is
//         |f_K - f*| = sum_{j >= K} g Delta_j rho = g Delta_K rho / (1 - rho)
//       -- proportional to the state tolerance, and to 1/(1 - rho) (the claim
//       of docs/cyclic-steady-state.md, "What the defaults deliver").
//
// Two variants of the shared PSA cycle (psa_cycle_case.h) span the range: the
// isothermal cycle with a slow LDF (k = 0.02 1/s) contracts by rho ~ 0.1 per
// cycle (a fast, clean geometric decay through five decades), and the cycle
// with its energy balance on contracts by rho ~ 0.966 (a slow thermal mode,
// the regime in which the stopping rule matters). The KPIs are those of the O2
// product (AddO2ProductKpis).
//
// The KPIs carry cycle-to-cycle noise at the integrator tolerance (about 1e-9
// for the fast case and 3e-6 for the slow one, at reltol 1e-8), so every
// comparison is restricted to cycles whose signal is far above it.

#include "psa_cycle_case.h"

#include <cstdio>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    struct CycleHistory {
      std::vector<double> state_change;  // Delta_k (k >= 1; entry 0 is unused)
      std::vector<O2Kpis> kpis;
      std::vector<int> css_stop;         // cycles at which a stopping rule held (filled by callers)
    };

    void BuildCssCase(UserData& data, bool isothermal, int n_cells, double mtc, int cycles,
                      bool accelerate) {
      BuildPsaCycle(data, n_cells, mtc);
      if (isothermal) data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
      AddO2ProductKpis(data);
      data.cycle.number_of_cycles_ = cycles;
      data.settings.reltol = 1e-8;
      data.settings.max_step = 0.0;
      data.settings.max_cycles_to_keep_ = cycles;  // the whole history
      data.settings.css.method_ = CSSMethod::NONE;
      data.settings.acc_css.method_ = accelerate ? AccelerationMethod::AITKEN : AccelerationMethod::NONE;
    }

    CycleHistory Collect(const UserData& data, const CoupledSolver& solver) {
      CycleHistory h;
      const auto& log = solver.GetLog();
      for (int c : log.sys_data.cycle_index_set_) {
        h.kpis.push_back(O2KpisAt(data, log, solver.GetBoundaryLog(), c));
        h.state_change.push_back(c > 0 ? CyclicSteadyState::EndOfCycleStateChange(log, c) : 0.0);
      }
      return h;
    }

    double Recovery(const O2Kpis& k) { return k.recovery; }
    double Purity(const O2Kpis& k) { return k.purity; }

    // Least-squares slope of ln(Delta_k) over cycles [k0, k1]: ln(rho).
    double FittedRho(const std::vector<double>& delta, int k0, int k1) {
      double sx = 0, sy = 0, sxx = 0, sxy = 0;
      const int n = k1 - k0 + 1;
      for (int k = k0; k <= k1; ++k) {
        const double y = std::log(delta[k]);
        sx += k; sy += y; sxx += double(k) * k; sxy += k * y;
      }
      return std::exp((n * sxy - sx * sy) / (n * sxx - sx * sx));
    }

    // Aitken's delta-squared limit of f at cycles k - 2m, k - m, k.
    double AitkenLimit(double f0, double f1, double f2) {
      const double d1 = f1 - f0, d2 = f2 - f1;
      return f2 - d2 * d2 / (d2 - d1);
    }

    // First cycle at which the deck's stopping rule (with `state_tolerance`)
    // holds, evaluated on the logged history exactly as the solver does after
    // every cycle; -1 if it never does.
    int StopCycle(const UserData& data, const CoupledSolver& solver, double state_tolerance) {
      CyclicSteadyState rule;
      rule.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
      rule.state_tolerance_ = state_tolerance;
      const auto& log = solver.GetLog();
      for (int c : log.sys_data.cycle_index_set_)
        if (rule.CheckCSS(log, c).converged) return c;
      (void)data;
      return -1;
    }

  } // namespace

  // ===========================================================================
  // Fast contraction: a clean geometric decay, and KPIs that follow it
  // ===========================================================================

  TEST(CSSConvergence, StateChangeDecaysGeometricallyAndTheKpisFollow) {
    constexpr int kCycles = 12;
    UserData data;
    BuildCssCase(data, /*isothermal=*/true, 20, /*mtc=*/0.02, kCycles, false);
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    CoupledSolver solver{data};
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);
    const CycleHistory h = Collect(data, solver);
    ASSERT_EQ(h.kpis.size(), static_cast<size_t>(kCycles));
    for (const auto& k : h.kpis) ASSERT_TRUE(k.defined);

    // The window: cycles whose state change is at least 30x the noise floor
    // (the state change of the last cycles, ~5e-9), leaving out the first
    // cycle (still far from the fixed point, where the map is not linear).
    const double floor = std::max(h.state_change[kCycles - 1], h.state_change[kCycles - 2]);
    std::vector<int> window;
    for (int k = 1; k < kCycles; ++k)
      if (h.state_change[k] > 30.0 * floor && h.state_change[k] < 1e-2) window.push_back(k);
    // GUARD: the decay spans at least four decades in the window.
    ASSERT_GE(window.size(), 4u);
    ASSERT_GT(h.state_change[window.front()] / h.state_change[window.back()], 1e3);

    // (1) Geometric decay: every successive ratio equals the fitted rate. The
    // ratios approach rho as (lambda_2 / rho)^k; this cycle's second mode is
    // already gone at the start of the window (measured ratios 0.108, 0.106,
    // 0.105, 0.105), so 10 % is a wide band around a real geometric sequence and
    // far from anything else.
    const double rho = FittedRho(h.state_change, window.front(), window.back());
    std::printf("  fast cycle: rho = %.4f over cycles %d-%d\n", rho, window.front(), window.back());
    EXPECT_LT(rho, 0.5) << "contraction rate " << rho;
    for (size_t i = 1; i < window.size(); ++i) {
      const double ratio = h.state_change[window[i]] / h.state_change[window[i - 1]];
      EXPECT_NEAR(ratio / rho, 1.0, 0.1) << "cycle " << window[i] << ": ratio " << ratio;
    }

    // (2) + (3) The KPIs converge at the same rate: their remaining error
    // e_k = |f_k - f*| (f* = the last cycle, at the noise floor) is a fixed
    // multiple of the state change, e_k = g rho / (1 - rho) Delta_k.
    for (auto kpi : {&Recovery, &Purity}) {
      const double f_star = kpi(h.kpis.back());
      std::vector<double> ratio;
      for (int k : window) {
        const double e = std::abs(kpi(h.kpis[k]) - f_star);
        ratio.push_back(e / h.state_change[k]);
      }
      // Leave out the last cycle of the window: its error is only ~100x the
      // floor, and f* itself carries the floor.
      for (size_t i = 1; i + 1 < ratio.size(); ++i)
        EXPECT_NEAR(ratio[i] / ratio[0], 1.0, 0.1)
          << "KPI error / state change at cycle " << window[i] << " is " << ratio[i]
          << ", at cycle " << window[0] << " " << ratio[0];
    }

    // The stopping rule: tightening the state tolerance by a decade costs one
    // cycle (rho ~ 0.1) and buys a decade of KPI accuracy -- the error left at
    // the stop is proportional to the tolerance. (At 1e-3 the rule's minimum
    // of five cycles decides the stop.)
    double prev_error = 0.0;
    int prev_stop = -1;
    for (double eps : {1e-3, 1e-4, 1e-5, 1e-6}) {
      const int stop = StopCycle(data, solver, eps);
      ASSERT_GE(stop, 0) << "the rule never held at state tolerance " << eps;
      ASSERT_LT(stop, kCycles - 2);
      const double e = std::abs(h.kpis[stop].recovery - h.kpis.back().recovery);
      std::printf("  state tolerance %.0e: stop at cycle %d, recovery error %.2e\n", eps, stop, e);
      // At the stop the last three state changes are below eps, and the error
      // is g rho / (1 - rho) times the last of them: below that factor times eps.
      EXPECT_LE(h.state_change[stop], eps);
      if (prev_stop >= 0) {
        // One more cycle per decade of tolerance (rho ~ 0.1), and the error
        // falls by rho per extra cycle: measured 0.105, 0.105, 0.106.
        EXPECT_GT(stop, prev_stop);
        EXPECT_NEAR(e / prev_error, std::pow(rho, stop - prev_stop), 0.1 * std::pow(rho, stop - prev_stop))
          << "the KPI error left at the stop did not fall with the tolerance as the geometric model predicts";
      }
      prev_error = e;
      prev_stop = stop;
    }
  }

  // ===========================================================================
  // Slow contraction: KPIs predicted from the geometric model
  // ===========================================================================
  //
  // With the energy balance on, the cycle has a slow thermal mode. From the
  // KPIs at cycles 39, 59 and 79 alone, (2) extrapolates the CSS value f* and
  // predicts the KPI at every later cycle; the run must follow the prediction
  // to within its noise. Then (3) gives the error the default stopping rule
  // leaves, which must match the measured one.

  TEST(CSSConvergence, SlowCycleFollowsTheGeometricModel) {
    constexpr int kCycles = 140;
    UserData data;
    BuildCssCase(data, /*isothermal=*/false, 10, /*mtc=*/0.5, kCycles, false);
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    CoupledSolver solver{data};
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);
    const CycleHistory h = Collect(data, solver);
    ASSERT_EQ(h.kpis.size(), static_cast<size_t>(kCycles));

    // (1) The slow mode: the rate over cycles 40-100 (state change 2e-4 to
    // 2.5e-5, far above the ~1e-6 noise floor).
    const double rho = FittedRho(h.state_change, 40, 100);
    std::printf("  slow cycle: rho = %.4f\n", rho);
    // GUARD: genuinely slow (so the stopping rule's error is amplified by
    // 1/(1 - rho) ~ 30) and a contraction.
    ASSERT_GT(rho, 0.9);
    ASSERT_LT(rho, 1.0);
    for (int k = 41; k <= 100; ++k)
      EXPECT_NEAR(h.state_change[k] / h.state_change[k - 1] / rho, 1.0, 0.02) << "cycle " << k;

    // Noise of the KPI: the RMS second difference over the last 20 cycles,
    // where the geometric signal is smallest. A second difference of white
    // noise of standard deviation s has standard deviation sqrt(6) s.
    auto noise_of = [&](double (*kpi)(const O2Kpis&)) {
      double s = 0.0;
      int n = 0;
      for (int k = kCycles - 20; k < kCycles; ++k, ++n) {
        const double d2 = kpi(h.kpis[k]) - 2.0 * kpi(h.kpis[k - 1]) + kpi(h.kpis[k - 2]);
        s += d2 * d2;
      }
      return std::sqrt(s / n / 6.0);
    };

    for (auto kpi : {&Recovery, &Purity}) {
      const double f39 = kpi(h.kpis[39]), f59 = kpi(h.kpis[59]), f79 = kpi(h.kpis[79]);
      const double f_star = AitkenLimit(f39, f59, f79);
      const double r20 = (f79 - f59) / (f59 - f39);  // rho^20 seen by the KPI
      const double noise = noise_of(kpi);
      std::printf("  KPI: f* = %.8f, rho from the KPI = %.4f, noise %.1e\n", f_star,
                  std::pow(r20, 1.0 / 20.0), noise);
      // The KPI converges at the state's rate.
      EXPECT_NEAR(std::pow(r20, 1.0 / 20.0), rho, 0.005);
      // GUARD: the remaining error at cycle 79 is far above the noise.
      ASSERT_GT(std::abs(f79 - f_star), 30.0 * noise);
      for (int n = 99; n < kCycles; n += 10) {
        const double predicted = f_star + (f79 - f_star) * std::pow(r20, (n - 79) / 20.0);
        EXPECT_LT(std::abs(kpi(h.kpis[n]) - predicted), 4.0 * noise + 0.01 * std::abs(f79 - f_star))
          << "cycle " << n << ": KPI " << kpi(h.kpis[n]) << ", geometric prediction " << predicted;
      }

      // (3) The error the default rule (state tolerance 2e-5, three cycles)
      // leaves, against g rho / (1 - rho) Delta_K, with g the KPI change per
      // unit state change measured on cycles 60-100.
      const int stop = StopCycle(data, solver, CyclicSteadyState::kDefaultStateTolerance);
      ASSERT_GT(stop, 0);
      ASSERT_LT(stop, kCycles);
      const double g = std::abs(kpi(h.kpis[100]) - kpi(h.kpis[60]))
                     / [&] { double s = 0; for (int k = 61; k <= 100; ++k) s += h.state_change[k]; return s; }();
      const double predicted_error = g * h.state_change[stop] * rho / (1.0 - rho);
      const double measured_error = std::abs(kpi(h.kpis[stop]) - f_star);
      std::printf("  default rule stops at cycle %d: error %.2e, predicted %.2e\n", stop,
                  measured_error, predicted_error);
      EXPECT_NEAR(measured_error / predicted_error, 1.0, 0.1);
    }
  }

  // ===========================================================================
  // Aitken acceleration lands on the same KPIs
  // ===========================================================================
  //
  // The accelerated run, stopped by the default rule, must report the KPIs of
  // the plain run's CSS: within the error the rule leaves on a plain run
  // ((3) above), so acceleration does not move the fixed point by more than
  // the tolerance already allows.

  TEST(CSSConvergence, AitkenAccelerationReachesTheSameKpis) {
    constexpr int kPlainCycles = 80;
    UserData plain_data;
    BuildCssCase(plain_data, false, 10, 0.5, kPlainCycles, false);
    std::string error;
    ASSERT_TRUE(plain_data.PerformChecks(error)) << error;
    CoupledSolver plain{plain_data};
    ASSERT_GE(plain.IntegrateCycle(plain_data.cycle, false), 0);
    const CycleHistory hp = Collect(plain_data, plain);
    const double rho = FittedRho(hp.state_change, 40, 79);

    UserData acc_data;
    BuildCssCase(acc_data, false, 10, 0.5, 200, true);
    acc_data.settings.css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;  // default tolerances
    ASSERT_TRUE(acc_data.PerformChecks(error)) << error;
    CoupledSolver acc{acc_data};
    ASSERT_GE(acc.IntegrateCycle(acc_data.cycle, false), 0);
    ASSERT_TRUE(acc.ReachedCSS()) << "accelerated run did not reach CSS in 200 cycles";
    const CycleHistory ha = Collect(acc_data, acc);
    const int acc_stop = static_cast<int>(ha.kpis.size()) - 1;

    for (auto kpi : {&Recovery, &Purity}) {
      const double f_star = AitkenLimit(kpi(hp.kpis[39]), kpi(hp.kpis[59]), kpi(hp.kpis[79]));
      // The plain-run stop error at the default tolerance, from (3).
      const double g = std::abs(kpi(hp.kpis[79]) - kpi(hp.kpis[59]))
                     / [&] { double s = 0; for (int k = 60; k <= 79; ++k) s += hp.state_change[k]; return s; }();
      const double allowed = g * CyclicSteadyState::kDefaultStateTolerance * rho / (1.0 - rho);
      const double e = std::abs(kpi(ha.kpis.back()) - f_star);
      std::printf("  accelerated stop at cycle %d: KPI %.8f, f* %.8f, error %.2e (plain-rule error %.2e)\n",
                  acc_stop, kpi(ha.kpis.back()), f_star, e, allowed);
      EXPECT_LT(e, allowed) << "accelerated CSS KPI is further from the CSS than the rule allows";
    }
  }

} // namespace adsim
