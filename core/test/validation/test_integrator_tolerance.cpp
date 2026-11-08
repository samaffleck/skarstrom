/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Time-integration error: tolerance proportionality, the absolute-tolerance
// floor, and independence of the maximum step size.
//
// NumericalVerification.SolutionIsIndependentOfIntegratorTolerance checks one
// pair of tolerances. These tests measure how the error falls over several
// decades of reltol, where it stops falling and why, and that max_step only
// changes which steps IDA takes, never the answer.
//
// The theory. IDA accepts a step when the weighted RMS norm of its local error
// estimate is below 1, with weights 1/(reltol |y_i| + abstol_i). This is
// error-per-step control of a variable-order BDF method; for order k the
// global error then scales as tol^(k/(k+1)) (Shampine, "Numerical Solution of
// Ordinary Differential Equations", 1994, sec. 6.3), between tol^0.5 (k = 1)
// and tol^0.83 (k = 5), and as tol^1 for error-per-unit-step control. So the
// observed slope of log(error) against log(reltol) must lie in [0.5, 1], and
// once reltol |y| falls below abstol the weights stop changing and so does the
// error: it floors at a small multiple of abstol. The absolute tolerances are
// fixed in the code (docs/numerics.md): 1e-8 for mole fractions [-] and
// loadings [mol/kg].

#include "breakthrough_case.h"
#include "psa_cycle_case.h"

#include <cstdio>

namespace adsim {

  using namespace adsim::validation;

  // kUpwindNote: these runs set Upwind, the scheme every number below was
  // measured with. With the default WENO the 20 -> 40 cell change of the
  // cycle recovery is 1.8e-5, below its error at the default reltol (5.2e-5),
  // so the "tolerance is not the dominant error" claim holds for upwind only;
  // and WENO takes more free steps (2153 on the breakthrough, 4414 on the
  // cycle), so a 0.1 s cap no longer doubles them and the guards that prove
  // the cap active would need a smaller cap.
  namespace {

    constexpr double kAbsTolMoleFraction = 1e-8;  // docs/numerics.md
    constexpr double kAbsTolLoading = 1e-8;       // [mol/kg], docs/numerics.md

    // Least-squares slope of log10(e) against log10(x).
    double LogLogSlope(const std::vector<double>& x, const std::vector<double>& e) {
      double sx = 0, sy = 0, sxx = 0, sxy = 0;
      const double n = static_cast<double>(x.size());
      for (size_t i = 0; i < x.size(); ++i) {
        const double a = std::log10(x[i]), b = std::log10(e[i]);
        sx += a; sy += b; sxx += a * a; sxy += a * b;
      }
      return (n * sxy - sx * sy) / (n * sxx - sx * sx);
    }

    // ------------------------------------------------------------------------
    // Breakthrough: final profiles at three times (front at 1/3, 1/2 and 2/3
    // of the bed), compared cell by cell at identical times (the step end).
    // ------------------------------------------------------------------------
    struct BreakthroughProfiles {
      std::vector<std::vector<double>> y;  // O2 mole fraction
      std::vector<std::vector<double>> q;  // O2 loading [mol/kg]
      long steps = 0;
    };

    BreakthroughProfiles RunBreakthrough(double reltol, double y_feed, double max_step = 0.0) {
      BreakthroughProfiles p;
      for (double t_end : {50.0, 75.0, 100.0}) {
        BreakthroughCase c;
        c.adsorbing = true;
        c.n_cells = 40;
        c.reltol = reltol;
        c.y_feed = y_feed;
        c.run_time = t_end;
        c.max_step = max_step;
        UserData data;
        Build(data, c);
        data.settings.advection_scheme_ = AdvectionScheme::Upwind;  // measured with upwind (kUpwindNote)
        std::string error;
        EXPECT_TRUE(data.PerformChecks(error)) << error;
        CoupledSolver solver{data};
        EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "reltol " << reltol;
        const auto& log = solver.GetLog();
        p.y.push_back(log.GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION, "Oxygen").values.back());
        p.q.push_back(log.GetComponentVariable(0, 0, ComponentVariableType::ADSORBED_CONCENTRATION, "Oxygen").values.back());
        p.steps += solver.GetStatistics().num_steps;
      }
      return p;
    }

    double MaxDiff(const std::vector<std::vector<double>>& a, const std::vector<std::vector<double>>& b) {
      double m = 0.0;
      for (size_t k = 0; k < a.size(); ++k)
        for (size_t i = 1; i + 1 < a[k].size(); ++i) m = std::max(m, std::abs(a[k][i] - b[k][i]));
      return m;
    }

    double MaxAbs(const std::vector<std::vector<double>>& a) {
      double m = 0.0;
      for (const auto& v : a)
        for (size_t i = 1; i + 1 < v.size(); ++i) m = std::max(m, std::abs(v[i]));
      return m;
    }

    // ------------------------------------------------------------------------
    // Cycle: the shared PSA cycle (energy balance on), two cycles on 20 cells.
    // ------------------------------------------------------------------------
    struct CycleResult {
      O2Kpis kpis;                // of the second cycle
      std::vector<double> q_n2;   // N2 loading profile at the end of the second cycle
      long steps = 0;
    };

    CycleResult RunCycle(double reltol, double max_step = 0.0) {
      UserData data;
      BuildPsaCycle(data, 20);
      data.settings.advection_scheme_ = AdvectionScheme::Upwind;  // measured with upwind (kUpwindNote)
      AddO2ProductKpis(data);
      data.cycle.number_of_cycles_ = 2;
      data.settings.reltol = reltol;
      data.settings.max_step = max_step;
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      CoupledSolver solver{data};
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "reltol " << reltol << ", max_step " << max_step;
      CycleResult r;
      const auto& log = solver.GetLog();
      r.kpis = O2KpisAt(data, log, solver.GetBoundaryLog(), 1);
      EXPECT_TRUE(r.kpis.defined);
      r.q_n2 = log.GetComponentVariable(1, 2, ComponentVariableType::ADSORBED_CONCENTRATION, "Nitrogen").values.back();
      r.steps = solver.GetStatistics().num_steps;
      return r;
    }

    double MaxRelDiff(const std::vector<double>& a, const std::vector<double>& b) {
      double m = 0.0, s = 0.0;
      for (size_t i = 1; i + 1 < a.size(); ++i) {
        m = std::max(m, std::abs(a[i] - b[i]));
        s = std::max(s, std::abs(b[i]));
      }
      return m / s;
    }

  } // namespace

  // ===========================================================================
  // 1. Breakthrough: error proportional to tol^(k/(k+1)), then the abstol floor
  // ===========================================================================

  TEST(IntegratorTolerance, BreakthroughErrorFallsWithReltolUntilTheAbstolFloor) {
    // y_feed = 0.5: reltol |y| exceeds abstol_y down to reltol = 2e-8, so the
    // range 1e-3 .. 1e-7 is reltol-controlled.
    const BreakthroughProfiles ref = RunBreakthrough(1e-12, 0.5);
    std::vector<double> tols = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7};
    std::vector<double> err;
    for (double rt : tols) {
      const BreakthroughProfiles p = RunBreakthrough(rt, 0.5);
      err.push_back(MaxDiff(p.y, ref.y));
      std::printf("  y_feed 0.5, reltol %.0e: max |y - y_ref| = %.2e (%ld steps)\n", rt, err.back(), p.steps);
    }
    // GUARD: the error at the loosest tolerance is far above the reference's
    // own (the floor, ~1e-8), so the slope measures the tolerance.
    ASSERT_GT(err.front(), 1e3 * kAbsTolMoleFraction);
    // The error falls with every decade of reltol down to 1e-6, by uneven
    // factors where IDA changes BDF order. The last point is compared over two
    // decades: its single-decade factor is noise from the order switch and
    // moves with the step sequence. Measured (1e-3 .. 1e-7):
    //   Linux GCC     2.57e-4, 2.94e-5, 2.88e-6, 1.30e-6, 3.12e-7
    //   Windows MSVC  2.57e-4, 2.94e-5, 2.88e-6, 3.39e-7, 3.60e-7
    //   before inert loadings left the unknowns (a different step sequence):
    //                 2.7e-4, 3.9e-5, 2.5e-6, 2.1e-6, 4.6e-7
    // The slope bounds below carry the convergence claim over the whole range.
    for (size_t i = 1; i + 1 < err.size(); ++i)
      EXPECT_LT(err[i], err[i - 1]) << "reltol " << tols[i];
    EXPECT_LT(err.back(), err[err.size() - 3]) << "reltol " << tols.back() << " vs " << tols[tols.size() - 3];
    const double slope = LogLogSlope(tols, err);
    std::printf("  observed slope d log(error) / d log(reltol) = %.2f\n", slope);
    EXPECT_GE(slope, 0.5) << "global error falls slower than error-per-step BDF control allows";
    EXPECT_LE(slope, 1.0) << "global error falls faster than tolerance proportionality";

    // The floor. With y_feed = 0.05 the loading is ~6e-3 mol/kg, so
    // reltol |q| < abstol_q for reltol < 1.6e-6: below that the weights, and
    // the error, stop changing. Between 1e-8 and 1e-10 the loading error must
    // stay at a small multiple of abstol_q (measured 0.75e-8 .. 1.4e-8 mol/kg)
    // instead of falling by the factor of 100 a proportional error would.
    const BreakthroughProfiles ref_dilute = RunBreakthrough(1e-12, 0.05);
    const double q_scale = MaxAbs(ref_dilute.q);
    ASSERT_LT(1e-8 * q_scale, kAbsTolLoading) << "case is not abstol-limited at reltol 1e-8";
    std::vector<double> floor;
    for (double rt : {1e-8, 1e-9, 1e-10}) {
      floor.push_back(MaxDiff(RunBreakthrough(rt, 0.05).q, ref_dilute.q));
      std::printf("  y_feed 0.05, reltol %.0e: max |q - q_ref| = %.2e mol/kg\n", rt, floor.back());
    }
    EXPECT_GT(floor.back(), 0.1 * floor.front()) << "loading error kept falling below the abstol crossover";
    for (double f : floor) {
      EXPECT_LT(f, 10.0 * kAbsTolLoading) << "loading error floor far above abstol";
      EXPECT_GT(f, 0.1 * kAbsTolLoading) << "loading error floor far below abstol";
    }
  }

  // ===========================================================================
  // 2. Cycle: the bed state converges with reltol, and at the default
  //    tolerance the KPI error is below the grid error
  // ===========================================================================
  //
  // The KPIs integrate the boundary flows over the logged points by the
  // trapezoid rule. A step that opens a pressure boundary drives a flow
  // transient whose time scale is the pressure diffusion time of one cell,
  // mu eb dz^2 / (K P) ~ 2e-5 s here, resolved only by IDA's own steps, so the
  // quadrature error of that transient falls with the step sizes IDA picks
  // there rather than with the local error: measured, the recovery converges
  // with reltol at a slope of 0.42 (1.4e-3 at reltol 1e-3, 6.2e-5 at the
  // default 1e-6, 3.0e-5 at 1e-7), and the end-of-cycle loadings (which do not
  // depend on the quadrature) at 0.77, inside the BDF range. The test asserts
  // the state rate from the theory above, that the KPI error falls
  // monotonically, and that at the default reltol it is below the change from
  // 20 to 40 cells.

  TEST(IntegratorTolerance, CycleStateConvergesAndDefaultToleranceKpiErrorIsBelowGridError) {
    const CycleResult ref = RunCycle(1e-11);
    std::vector<double> tols = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7};
    std::vector<double> err, err_rec;
    for (double rt : tols) {
      const CycleResult r = RunCycle(rt);
      err.push_back(MaxRelDiff(r.q_n2, ref.q_n2));
      err_rec.push_back(std::abs(r.kpis.recovery - ref.kpis.recovery));
      std::printf("  reltol %.0e: end-of-cycle N2 loading error %.2e (relative), recovery error %.2e\n",
                  rt, err.back(), err_rec.back());
    }
    ASSERT_GT(err.front(), 1e-5);  // GUARD: far above the reference's own error
    for (size_t i = 1; i < err.size(); ++i) EXPECT_LT(err[i], err[i - 1]) << "reltol " << tols[i];
    const double slope = LogLogSlope(tols, err);
    const double slope_kpi = LogLogSlope(tols, err_rec);
    std::printf("  slope: state %.2f, recovery %.2f\n", slope, slope_kpi);
    EXPECT_GE(slope, 0.5);
    EXPECT_LE(slope, 1.0);
    // The KPI error still falls with reltol, monotonically.
    for (size_t i = 1; i < err_rec.size(); ++i) EXPECT_LT(err_rec[i], err_rec[i - 1]) << "reltol " << tols[i];

    // The grid error of the same KPI: 20 against 40 cells at reltol 1e-8.
    UserData fine;
    BuildPsaCycle(fine, 40);
    fine.settings.advection_scheme_ = AdvectionScheme::Upwind;  // as RunCycle
    AddO2ProductKpis(fine);
    fine.cycle.number_of_cycles_ = 2;
    fine.settings.reltol = 1e-8;
    fine.settings.max_step = 0.0;
    CoupledSolver fine_solver{fine};
    ASSERT_GE(fine_solver.IntegrateCycle(fine.cycle, false), 0);
    const double rec_40 = O2KpisAt(fine, fine_solver.GetLog(), fine_solver.GetBoundaryLog(), 1).recovery;
    const double grid_change = std::abs(rec_40 - RunCycle(1e-8).kpis.recovery);
    const double default_tol_error = err_rec[3];  // reltol 1e-6
    std::printf("  recovery: grid change 20 -> 40 cells %.2e, error at the default reltol %.2e\n",
                grid_change, default_tol_error);
    // At the default reltol the tolerance is not the dominant error even on
    // the coarse grid (measured: 6.2e-5 against a 20 -> 40 cell change of
    // 2.4e-4 after two cycles).
    EXPECT_LT(default_tol_error, grid_change);
  }

  // ===========================================================================
  // 3. max_step: the answer does not depend on it
  // ===========================================================================
  //
  // max_step caps the step IDA may take; it never relaxes the error test. So a
  // cap above the steps IDA picks anyway changes nothing at all, and a cap
  // below them changes the result only by re-sampling the integration error:
  // by less than the tolerance error of a looser reltol. Decks that set
  // max_step (Vysyaraju 2024 uses 0.5 s) do so for robustness, not accuracy:
  // with a 1 s cap IDA stalled at its minimum step in one equalisation step of
  // that deck, and with 0.5 s (or reltol 1e-7 at 1 s) the run completes and the
  // two agree to 0.01 points in purity and recovery
  // (core/test/papers/vysyarajua_2024).

  TEST(IntegratorTolerance, BreakthroughDoesNotDependOnMaxStep) {
    const BreakthroughProfiles free_run = RunBreakthrough(1e-8, 0.5, 0.0);
    // A cap far above any step IDA takes on these 50-100 s runs is inactive:
    // the run is bit-for-bit the same.
    const BreakthroughProfiles inactive = RunBreakthrough(1e-8, 0.5, 10.0);
    EXPECT_EQ(inactive.steps, free_run.steps);
    EXPECT_EQ(MaxDiff(inactive.y, free_run.y), 0.0);
    // Tighter caps (measured: 2 s and 0.5 s leave this run unchanged too, 0.1 s
    // more than doubles the steps): the same answer within the tolerance error
    // of a 10x looser reltol.
    const double tol_error = MaxDiff(RunBreakthrough(1e-7, 0.5).y, free_run.y);
    for (double cap : {2.0, 0.5, 0.1}) {
      const BreakthroughProfiles capped = RunBreakthrough(1e-8, 0.5, cap);
      const double cap_change = MaxDiff(capped.y, free_run.y);
      std::printf("  max_step %.1f s: %ld steps (free %ld), profile change %.2e, tolerance error %.2e\n", cap,
                  capped.steps, free_run.steps, cap_change, tol_error);
      if (cap == 0.1) ASSERT_GT(capped.steps, 2 * free_run.steps) << "the 0.1 s cap is not active";  // GUARD
      EXPECT_LT(cap_change, tol_error) << "max_step " << cap;
    }
  }

  TEST(IntegratorTolerance, CycleKpisAndStateDoNotDependOnMaxStep) {
    const CycleResult free_run = RunCycle(1e-8, 0.0);
    const CycleResult loose = RunCycle(1e-6, 0.0);
    const double tol_rec = std::abs(loose.kpis.recovery - free_run.kpis.recovery);
    const double tol_pur = std::abs(loose.kpis.purity - free_run.kpis.purity);
    const double tol_e = std::abs(loose.kpis.specific_energy - free_run.kpis.specific_energy);
    const double tol_q = MaxRelDiff(loose.q_n2, free_run.q_n2);
    for (double cap : {10.0, 2.0, 0.5, 0.1}) {
      const CycleResult r = RunCycle(1e-8, cap);
      const double d_rec = std::abs(r.kpis.recovery - free_run.kpis.recovery);
      const double d_pur = std::abs(r.kpis.purity - free_run.kpis.purity);
      const double d_e = std::abs(r.kpis.specific_energy - free_run.kpis.specific_energy);
      const double d_q = MaxRelDiff(r.q_n2, free_run.q_n2);
      std::printf("  max_step %4.1f s: %ld steps, recovery %+.1e, purity %+.1e, energy %+.1e, loading %.1e\n",
                  cap, r.steps, d_rec, d_pur, d_e / free_run.kpis.specific_energy, d_q);
      if (cap == 0.1) ASSERT_GT(r.steps, 2 * free_run.steps) << "the 0.1 s cap is not active";  // GUARD
      // Bounded by the tolerance error of a 100x looser reltol.
      EXPECT_LT(d_rec, tol_rec) << "max_step " << cap;
      EXPECT_LT(d_pur, tol_pur) << "max_step " << cap;
      EXPECT_LT(d_e, tol_e) << "max_step " << cap;
      EXPECT_LT(d_q, tol_q) << "max_step " << cap;
    }
  }

} // namespace adsim
