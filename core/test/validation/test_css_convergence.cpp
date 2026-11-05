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
//       rate: f_k = f* + A rho^k, so with rho known from the state, f* and A
//       follow from a linear least-squares fit over many cycles and the KPI of
//       any later cycle can be predicted;
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
// The KPIs carry cycle-to-cycle noise that depends on the integrator's step
// sequence (they integrate the logged boundary flows, and the flow transient
// after a pressure switch is sampled only at IDA's own steps; see
// test_integrator_tolerance.cpp). Its size differs between platforms: at
// reltol 1e-8 about 1e-9 for the fast case and 2e-6 for the slow one on
// Linux, and it is not monotone in reltol (1e-6 at reltol 1e-9). So the tests
// measure it in the run itself and compare only quantities whose signal is
// far above it, and every CSS value f* is a least-squares fit over many
// cycles rather than an extrapolation from three (Aitken's delta-squared
// amplifies the noise of its three samples about five-fold: on macOS that
// alone moved f* by 15 % of the error under test).

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

    // The log-linear least-squares fit of Delta over cycles [k0, k1],
    // evaluated at cycle k.
    double FittedDelta(const std::vector<double>& delta, int k0, int k1, int k) {
      double sx = 0, sy = 0, sxx = 0, sxy = 0;
      const int n = k1 - k0 + 1;
      for (int j = k0; j <= k1; ++j) {
        const double y = std::log(delta[j]);
        sx += j; sy += y; sxx += double(j) * j; sxy += j * y;
      }
      const double slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
      return std::exp((sy - slope * sx) / n + slope * k);
    }

    // Least-squares fit of f_k = f* + A rho^(k - k0) over cycles [k0, k1] with
    // rho given (linear in f* and A): f*, the RMS residual (the KPI noise,
    // with 2 degrees of freedom removed), the standard error of f*, and the
    // fitted value and its standard error at any cycle.
    struct GeometricTail {
      double f_star = 0.0, amplitude = 0.0, rho = 0.0, noise = 0.0;
      int k0 = 0;
      double inv[2][2] = {{0, 0}, {0, 0}};  // (X^T X)^-1
      double At(int k) const { return f_star + amplitude * std::pow(rho, k - k0); }
      double StdErrAt(int k) const {
        const double x = std::pow(rho, k - k0);
        return noise * std::sqrt(inv[0][0] + 2.0 * x * inv[0][1] + x * x * inv[1][1]);
      }
      double StdErrFStar() const { return noise * std::sqrt(inv[0][0]); }
    };

    GeometricTail FitTail(const std::vector<O2Kpis>& kpis, double (*kpi)(const O2Kpis&), double rho,
                          int k0, int k1) {
      double s1 = 0, sx = 0, sxx = 0, sy = 0, sxy = 0;
      for (int k = k0; k <= k1; ++k) {
        const double x = std::pow(rho, k - k0), y = kpi(kpis[k]);
        s1 += 1; sx += x; sxx += x * x; sy += y; sxy += x * y;
      }
      const double det = s1 * sxx - sx * sx;
      GeometricTail t;
      t.rho = rho;
      t.k0 = k0;
      t.amplitude = (s1 * sxy - sx * sy) / det;
      t.f_star = (sy - t.amplitude * sx) / s1;
      t.inv[0][0] = sxx / det; t.inv[0][1] = t.inv[1][0] = -sx / det; t.inv[1][1] = s1 / det;
      double ss = 0.0;
      for (int k = k0; k <= k1; ++k) {
        const double r = kpi(kpis[k]) - t.At(k);
        ss += r * r;
      }
      t.noise = std::sqrt(ss / (k1 - k0 + 1 - 2));
      return t;
    }

    // First cycle at which the deck's stopping rule (with `state_tolerance`)
    // holds, evaluated on the logged history exactly as the solver does after
    // every cycle; -1 if it never does.
    // min_cycles: the rule's minimum number of cycles (default 5).
    int StopCycle(const UserData& data, const CoupledSolver& solver, double state_tolerance,
                  int min_cycles = CyclicSteadyState{}.min_cycles_) {
      CyclicSteadyState rule;
      rule.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
      rule.state_tolerance_ = state_tolerance;
      rule.min_cycles_ = min_cycles;
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
    constexpr int kCycles = 14;
    constexpr int kTail = 5;  // the last cycles, all at the noise floor
    UserData data;
    BuildCssCase(data, /*isothermal=*/true, 20, /*mtc=*/0.02, kCycles, false);
    // Measured with upwind. With the default WENO the state change settles on
    // a higher noise floor (the nonlinear weights respond to tolerance-level
    // differences), so the window above 30x the floor spans 2.9 decades
    // (876x) instead of the guard's 3.
    data.settings.advection_scheme_ = AdvectionScheme::Upwind;
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
    double floor = 0.0;
    for (int k = kCycles - kTail; k < kCycles; ++k) floor = std::max(floor, h.state_change[k]);
    std::vector<int> window;
    for (int k = 1; k < kCycles - kTail; ++k)
      if (h.state_change[k] > 30.0 * floor && h.state_change[k] < 1e-2) window.push_back(k);
    // GUARD: the decay spans at least three decades in the window.
    ASSERT_GE(window.size(), 4u);
    ASSERT_GT(h.state_change[window.front()] / h.state_change[window.back()], 1e3);

    // (1) Geometric decay: every successive ratio equals the fitted rate. The
    // ratios approach rho as (lambda_2 / rho)^k; this cycle's second mode is
    // already gone at the start of the window (measured ratios 0.108, 0.106,
    // 0.105, 0.105 on Linux; rho = 0.106 on Linux, 0.104 on macOS), so 10 % is
    // a wide band around a real geometric sequence and far from anything else.
    const double rho = FittedRho(h.state_change, window.front(), window.back());
    std::printf("  fast cycle: rho = %.4f over cycles %d-%d\n", rho, window.front(), window.back());
    EXPECT_LT(rho, 0.5) << "contraction rate " << rho;
    for (size_t i = 1; i < window.size(); ++i) {
      const double ratio = h.state_change[window[i]] / h.state_change[window[i - 1]];
      EXPECT_NEAR(ratio / rho, 1.0, 0.1) << "cycle " << window[i] << ": ratio " << ratio;
    }

    // The KPI's CSS value and noise from the last cycles: f* their mean, the
    // noise their largest deviation from it (about 2 standard deviations).
    // An error e measured against it carries a relative noise of about
    // noise / e, so a ratio of two such errors is good to sqrt(2) noise / e:
    // within 5 % (half the 10 % bands below) for e > 30 noise. Only errors
    // above that are compared.
    auto tail_of = [&](double (*kpi)(const O2Kpis&), double& noise) {
      double mean = 0.0;
      for (int k = kCycles - kTail; k < kCycles; ++k) mean += kpi(h.kpis[k]) / kTail;
      noise = 0.0;
      for (int k = kCycles - kTail; k < kCycles; ++k) noise = std::max(noise, std::abs(kpi(h.kpis[k]) - mean));
      return mean;
    };

    // (2) + (3) The KPIs converge at the same rate: their remaining error
    // e_k = |f_k - f*| is a fixed multiple of the state change,
    // e_k = g rho / (1 - rho) Delta_k.
    for (auto kpi : {&Recovery, &Purity}) {
      double noise = 0.0;
      const double f_star = tail_of(kpi, noise);
      std::vector<double> ratio;
      std::vector<int> cycles;
      for (int k : window) {
        const double e = std::abs(kpi(h.kpis[k]) - f_star);
        if (e < 30.0 * noise) continue;
        ratio.push_back(e / h.state_change[k]);
        cycles.push_back(k);
      }
      // GUARD: at least three cycles above the noise (measured: cycles 2-6,
      // errors 9e-4 to 1e-7 against a noise of ~2e-9 in recovery).
      ASSERT_GE(ratio.size(), 3u) << "KPI noise " << noise << " leaves too few cycles to compare";
      for (size_t i = 1; i < ratio.size(); ++i)
        EXPECT_NEAR(ratio[i] / ratio[0], 1.0, 0.1)
          << "KPI error / state change at cycle " << cycles[i] << " is " << ratio[i]
          << ", at cycle " << cycles[0] << " " << ratio[0];
    }

    // The stopping rule (minimum cycle count lowered to 1, so that the state
    // tolerance alone decides): tightening the state tolerance by a decade
    // costs one cycle (rho ~ 0.1) and buys a decade of KPI accuracy -- the
    // error left at the stop is proportional to the tolerance. Measured on
    // Linux: stops at cycles 4 to 8 for 1e-2 to 1e-6, recovery errors 9.9e-6,
    // 1.0e-6, 1.1e-7, then 1.2e-8 and 2e-9, which are at the KPI noise (and
    // on macOS the noise is several times larger), so only the steps whose
    // error is above 30x the measured noise compare errors; every step checks
    // the state change and the cycle count.
    double noise = 0.0;
    const double f_star = tail_of(&Recovery, noise);
    double prev_error = 0.0;
    int prev_stop = -1, compared = 0;
    for (double eps : {1e-2, 1e-3, 1e-4, 1e-5, 1e-6}) {
      const int stop = StopCycle(data, solver, eps, /*min_cycles=*/1);
      ASSERT_GE(stop, 0) << "the rule never held at state tolerance " << eps;
      ASSERT_LT(stop, kCycles - kTail);
      const double e = std::abs(h.kpis[stop].recovery - f_star);
      std::printf("  state tolerance %.0e: stop at cycle %d, recovery error %.2e (noise %.1e)\n", eps, stop, e, noise);
      EXPECT_LE(h.state_change[stop], eps);
      if (prev_stop >= 0) {
        EXPECT_GT(stop, prev_stop);
        if (e > 30.0 * noise) {
          // The error falls by rho per extra cycle (measured 0.105).
          EXPECT_NEAR(e / prev_error, std::pow(rho, stop - prev_stop), 0.1 * std::pow(rho, stop - prev_stop))
            << "the KPI error left at the stop did not fall with the tolerance as the geometric model predicts";
          ++compared;
        }
      }
      prev_error = e;
      prev_stop = stop;
    }
    EXPECT_GE(compared, 1) << "no tolerance step had a KPI error above the noise";
  }

  // ===========================================================================
  // Slow contraction: KPIs predicted from the geometric model
  // ===========================================================================
  //
  // With the energy balance on, the cycle has a slow thermal mode. From the
  // KPIs of cycles 40-79 alone, (2) fits the CSS value f* and predicts the KPI
  // at every later cycle; the run must follow the prediction to within its
  // noise. Then (3) gives the error the default stopping rule leaves, which
  // must match the measured one to within the uncertainty of the measurement.

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
    std::printf("  slow cycle: rho = %.5f\n", rho);
    // GUARD: genuinely slow (so the stopping rule's error is amplified by
    // 1/(1 - rho) ~ 30) and a contraction.
    ASSERT_GT(rho, 0.9);
    ASSERT_LT(rho, 1.0);
    for (int k = 41; k <= 100; ++k)
      EXPECT_NEAR(h.state_change[k] / h.state_change[k - 1] / rho, 1.0, 0.02) << "cycle " << k;

    const int stop = StopCycle(data, solver, CyclicSteadyState::kDefaultStateTolerance);
    ASSERT_GT(stop, 100);
    ASSERT_LT(stop, kCycles - 1);

    for (auto kpi : {&Recovery, &Purity}) {
      const char* name = kpi == &Recovery ? "recovery" : "purity";
      // The KPI converges at the state's rate: the ratio of its changes over
      // cycles 39-59 and 59-79 is rho^20. Each change (~1.7e-4 in recovery)
      // carries sqrt(2) x the noise (~2e-6), so the rate is good to about
      // 0.001; 0.005 is five times that.
      const double r20 = (kpi(h.kpis[79]) - kpi(h.kpis[59])) / (kpi(h.kpis[59]) - kpi(h.kpis[39]));
      const double rho_kpi = std::pow(r20, 1.0 / 20.0);
      EXPECT_NEAR(rho_kpi, rho, 0.005) << name;

      // (2) Prediction: f* and A from cycles 40-79, with the state's rate.
      // GUARD: the remaining error at cycle 79 is far above the noise.
      const GeometricTail early = FitTail(h.kpis, kpi, rho, 40, 79);
      std::printf("  %s: f* from cycles 40-79 = %.8f, noise %.1e\n", name, early.f_star, early.noise);
      ASSERT_GT(std::abs(kpi(h.kpis[79]) - early.f_star), 30.0 * early.noise);
      for (int n = 99; n < kCycles; n += 10) {
        // Noise of the cycle plus the extrapolation error of the fit, 4 sigma,
        // plus the effect of the rate's own uncertainty (1 % of the remaining
        // error at cycle 79, for an error of 2.5e-4 in rho over 60 cycles).
        const double tol = 4.0 * std::hypot(early.noise, early.StdErrAt(n))
                         + 0.01 * std::abs(kpi(h.kpis[79]) - early.f_star);
        EXPECT_LT(std::abs(kpi(h.kpis[n]) - early.At(n)), tol)
          << name << ", cycle " << n << ": KPI " << kpi(h.kpis[n]) << ", geometric prediction " << early.At(n);
      }

      // (3) The error the default rule (state tolerance 2e-5, three cycles)
      // leaves, against g Delta_K rho / (1 - rho). Each factor is measured so
      // that its noise is small and known:
      //   f*      least-squares fit over cycles 40-139 (standard error from the
      //           fit, ~1e-6 in recovery);
      //   Delta_K the state change at the stop from a log-linear fit over the
      //           30 cycles before it (the raw value carries ~1e-6 noise);
      //   g       the KPI change over cycles 40-130 per unit summed state
      //           change (the change, ~4e-4, carries sqrt(2) x the noise).
      // The measured error |f_K - f*| itself carries the noise of f_K. The
      // tolerance is 3 standard deviations of the ratio from these, plus the
      // ambiguity of rho in the model (state or KPI rate: |rho_kpi - rho| /
      // (1 - rho), 0 on Linux and 2 % on macOS). Measured on Linux: recovery
      // 6.11e-5 against 6.21e-5 predicted (ratio 0.985, tolerance 0.108),
      // purity 7.20e-5 against 7.29e-5 (0.988, tolerance 0.038). With f* from
      // an Aitken extrapolation of cycles 39, 59 and 79 the ratios were 0.94
      // and 0.92 on Linux and 0.885 on macOS: the extrapolation's noise, not
      // the model.
      const GeometricTail full = FitTail(h.kpis, kpi, rho, 40, kCycles - 1);
      const double delta_stop = FittedDelta(h.state_change, stop - 30, stop, stop);
      double sum_delta = 0.0;
      for (int k = 41; k <= 130; ++k) sum_delta += h.state_change[k];
      const double df = kpi(h.kpis[130]) - kpi(h.kpis[40]);
      const double g = std::abs(df) / sum_delta;
      const double predicted_error = g * delta_stop * rho / (1.0 - rho);
      const double measured_error = std::abs(kpi(h.kpis[stop]) - full.f_star);
      const double rel_sigma = std::sqrt(std::pow(full.noise / measured_error, 2)
                                         + std::pow(full.StdErrFStar() / measured_error, 2)
                                         + std::pow(std::sqrt(2.0) * full.noise / std::abs(df), 2));
      const double tol = 3.0 * rel_sigma + std::abs(rho_kpi - rho) / (1.0 - rho);
      std::printf("  %s: default rule stops at cycle %d: error %.3e, predicted %.3e (ratio %.3f, "
                  "tolerance %.3f; noise %.1e, f* %.8f +- %.1e)\n", name, stop, measured_error,
                  predicted_error, measured_error / predicted_error, tol, full.noise, full.f_star,
                  full.StdErrFStar());
      EXPECT_NEAR(measured_error / predicted_error, 1.0, tol) << name;
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
      const double f_star = FitTail(hp.kpis, kpi, rho, 40, kPlainCycles - 1).f_star;
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
