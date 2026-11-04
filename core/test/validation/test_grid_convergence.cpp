/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Grid convergence of the cycle KPIs at cyclic steady state.
//
// test_numerical_verification and test_advection_schemes measure the order of
// accuracy on single breakthrough curves. What a user reports, though, is a
// cycle KPI -- purity, recovery, productivity, specific energy -- at cyclic
// steady state, and nothing so far said how far those are from their
// grid-converged values, or at what rate they get there. These tests run the
// shared PSA cycle (psa_cycle_case.h) to CSS on four grids, each twice as fine
// as the last (10, 20, 40, 80 cells), with each advection scheme, and apply the
// standard verification analysis (Roache, "Verification and Validation in
// Computational Science and Engineering", 1998; ASME V&V 20):
//
//   * Monotone convergence: the successive differences d_k = f(2N) - f(N) keep
//     one sign and shrink.
//   * Observed order from the three finest grids (Richardson):
//         p = log2( d(20 -> 40) / d(40 -> 80) ).
//     The KPI error of a consistent scheme is E(h) = a h^q + b h^(q+1) + ...
//     with q its formal order. When the two leading terms have the same sign,
//     the difference ratio d(h)/d(h/2) = 2^q (a + c1 b h)/(a + c2 b h) with
//     c1 > c2 > 0 lies between 2^q and 2^(q+1), so the observed order on a
//     pre-asymptotic grid lies in [q, q + 1] and falls toward q as the grid is
//     refined. For upwind q = 1 (the leading term is the numerical diffusion
//     u dx / 2), so p must lie in [1, 2].
//   * Grid convergence index of the finest grid, safety factor Fs = 1.25:
//         GCI = Fs |d(40 -> 80) / f(80)| / (2^p - 1),
//     a band around f(80) meant to contain the grid-converged value. The test
//     checks that claim directly: the band of the upwind KPIs must contain the
//     WENO and van Leer KPIs at 80 cells, which are an order of magnitude closer
//     to the limit (an independent estimate of it from a different scheme).
//
// The higher-order schemes. The WENO and van Leer interior faces are third and
// (limited) second order, but the two bed-end faces take the upwind value in
// every scheme (the ghost cell holds a face value; see test_advection_schemes),
// so a quantity read at a bed end is formally first order for all three. The
// KPIs here are read at the feed end (where the face carries the feed
// composition exactly) and the product end, where the gradient is small, so the
// first-order boundary term has a small coefficient: on these grids the
// higher-order schemes converge at p = 2.4 to 2.5, against 1.55 for upwind
// (measured, below). The test requires them to converge faster than upwind (a
// larger observed order and a smaller finest-grid difference) rather than a
// specific order, since their asymptotic order is set by the boundary.
//
// Measured (Release build, reltol 1e-8; O2 product purity and recovery, and the
// specific energy [J/kg]; p from 20/40/80 cells, GCI relative to f(80)):
//
//   KPI       scheme    f(10)    f(20)    f(40)    f(80)     p     GCI      Richardson
//   purity    upwind    0.38801  0.39304  0.39503  0.39572   1.54  1.2e-3   0.39608
//             WENO      0.39339  0.39555  0.39602  0.39611   2.50  5.7e-5   0.39613
//             van Leer  0.39420  0.39577  0.39606  0.39611   2.47  3.6e-5   0.39612
//   recovery  upwind    0.67868  0.68348  0.68533  0.68597   1.55  6.0e-4   0.68630
//             WENO      0.68377  0.68581  0.68625  0.68633   2.48  3.1e-5   0.68634
//             van Leer  0.68453  0.68601  0.68628  0.68633   2.40  2.2e-5   0.68634
//   energy    upwind    854701   848717   846428   845650    1.56  5.9e-4   845250
//             WENO      848354   845841   845308   845213    2.48  3.1e-5   845192
//             van Leer  847419   845599   845271   845209    2.41  2.1e-5   845194
//
// and, offline, upwind on 80/160/320 and 160/320/640 cells gives p = 1.37 and 1.26: the
// observed order falls toward 1 as expected. The upwind GCI at 80 cells is
// 0.12 % (purity) and 0.06 % (recovery), i.e. about 4e-4 absolute. The
// example decks need far more cells than this cycle: with upwind on their
// 30 cells, Perez2019_A and Ward2024 are 0.6 and 0.8 points of recovery from
// their grid-converged values (scripts/verification/grid_convergence.py;
// docs/numerics.md, "Convergence of the results").
//
// Why isothermal. With the energy balance on, this cycle approaches CSS with a
// slow thermal mode (the wall), contracting by only 0.966 per cycle, so its KPIs
// move by ~1e-5 per cycle even after 40 cycles. Isothermal, the end-of-cycle
// map contracts by about 0.04 per cycle and the state change is below 1e-7
// within five cycles, so nine cycles leave a CSS error far below the grid
// differences (GUARD 2 checks it) and the whole study runs in seconds.
//
// KPI noise. At CSS the KPIs still vary from cycle to cycle, by an amount set
// by the integrator's step sequence (they integrate the logged boundary
// flows, and the flow transient after a pressure switch is sampled only at
// IDA's own steps; see test_integrator_tolerance.cpp), so it differs between
// platforms and between Jacobian builds: up to 5.6e-6 in the van Leer
// recovery on macOS once the coloured Jacobian of the higher-order schemes
// matched the band one exactly. Each KPI is therefore the mean of the four
// CSS cycles 5-8, and its noise the standard error of that mean.
//
// The productivity of this cycle is its recovery times a constant (the feed is
// a fixed molar flow), so it is checked through the recovery; the specific
// energy (feed compressor work per kg of O2 produced) depends on the inlet
// pressure history as well and is checked on its own.
//
// GUARDS: (1) the grid effect is present and large, far above the noise of
// the KPIs at CSS: the coarsest difference above 100 and the finest above 5
// standard deviations of a difference (sqrt(2) x the standard error). Then
// the signs are certain and the observed order is good to about 0.3
// (1/ln 2 x 1/5), a third of the smallest margin asserted on it (WENO and
// van Leer exceed upwind's order by 0.85 or more; measured noise on Linux:
// standard error up to 1.7e-6 in recovery, finest differences 4.8e-5 to
// 6.4e-4); (2) the bed is at CSS (end-of-cycle state change below 1e-6 from
// the first averaged cycle on).

#include "psa_cycle_case.h"

#include <array>
#include <cstdio>
#include <map>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    constexpr std::array<int, 4> kGrids = {10, 20, 40, 80};
    constexpr int kCycles = 9;
    constexpr int kFirstAveraged = 5;  // cycles 5..8 are averaged
    constexpr double kSafetyFactor = 1.25;  // Roache's factor for a three-grid study

    enum Kpi { kPurity = 0, kRecovery = 1, kEnergy = 2, kNumKpis = 3 };
    const char* KpiName(int k) {
      static const char* names[] = {"purity", "recovery", "specific energy"};
      return names[k];
    }

    struct GridRun {
      std::array<double, kNumKpis> value{};  // mean over the averaged CSS cycles
      std::array<double, kNumKpis> noise{};  // standard error of that mean
      double state_change = 0.0;             // largest end-of-cycle state change of those cycles
    };

    // The shared PSA cycle, isothermal, with the O2 product KPIs.
    void BuildKpiCycle(UserData& data, int n_cells, AdvectionScheme scheme) {
      BuildPsaCycle(data, n_cells);
      data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
      AddO2ProductKpis(data);
      data.cycle.number_of_cycles_ = kCycles;
      data.settings.reltol = 1e-8;
      data.settings.max_step = 0.0;
      data.settings.max_cycles_to_keep_ = kCycles;
      data.settings.css.method_ = CSSMethod::NONE;
      data.settings.advection_scheme_ = scheme;
    }

    std::array<double, kNumKpis> KpisAt(const UserData& data, const CoupledSolver& solver, int cycle) {
      const O2Kpis k = O2KpisAt(data, solver.GetLog(), solver.GetBoundaryLog(), cycle);
      EXPECT_TRUE(k.defined) << "KPI undefined at cycle " << cycle;
      return {k.purity, k.recovery, k.specific_energy};
    }

    GridRun RunGrid(int n_cells, AdvectionScheme scheme) {
      GridRun out;
      UserData data;
      BuildKpiCycle(data, n_cells, scheme);
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      CoupledSolver solver{data};
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0)
        << "Solver did not complete (scheme " << static_cast<int>(scheme) << ", " << n_cells << " cells)";
      constexpr int n = kCycles - kFirstAveraged;
      std::array<std::array<double, kNumKpis>, n> f{};
      for (int c = kFirstAveraged; c < kCycles; ++c) {
        f[c - kFirstAveraged] = KpisAt(data, solver, c);
        out.state_change = std::max(out.state_change,
                                    CyclicSteadyState::EndOfCycleStateChange(solver.GetLog(), c));
      }
      for (int k = 0; k < kNumKpis; ++k) {
        double mean = 0.0, var = 0.0;
        for (const auto& row : f) mean += row[k] / n;
        for (const auto& row : f) var += (row[k] - mean) * (row[k] - mean) / (n - 1);
        out.value[k] = mean;
        out.noise[k] = std::sqrt(var / n);
      }
      return out;
    }

    // The four-grid study of one scheme, run once per test binary.
    const std::array<GridRun, kGrids.size()>& Study(AdvectionScheme scheme) {
      static std::map<int, std::array<GridRun, kGrids.size()>> cache;
      auto it = cache.find(static_cast<int>(scheme));
      if (it == cache.end()) {
        std::array<GridRun, kGrids.size()> runs;
        for (size_t i = 0; i < kGrids.size(); ++i) runs[i] = RunGrid(kGrids[i], scheme);
        it = cache.emplace(static_cast<int>(scheme), runs).first;
      }
      return it->second;
    }

    struct Richardson {
      std::array<double, kGrids.size() - 1> d{};  // f(2N) - f(N), coarse to fine
      double order = 0.0;         // from the three finest grids
      double gci = 0.0;           // relative, of the finest grid
      double extrapolated = 0.0;  // f_fine + d_fine / (2^p - 1)
      double fine = 0.0;
    };

    Richardson Analyse(const std::array<GridRun, kGrids.size()>& runs, int k) {
      Richardson r;
      for (size_t i = 0; i + 1 < runs.size(); ++i) r.d[i] = runs[i + 1].value[k] - runs[i].value[k];
      const double d_mid = r.d[r.d.size() - 2], d_fine = r.d.back();
      r.fine = runs.back().value[k];
      r.order = std::log2(d_mid / d_fine);
      const double denom = std::pow(2.0, r.order) - 1.0;
      r.gci = kSafetyFactor * std::abs(d_fine / r.fine) / denom;
      r.extrapolated = r.fine + d_fine / denom;
      return r;
    }

    // Guards and monotone convergence, common to every scheme.
    void CheckMonotoneConvergence(AdvectionScheme scheme, const char* name) {
      const auto& runs = Study(scheme);
      for (const auto& run : runs) {
        // GUARD 2: at CSS.
        EXPECT_LT(run.state_change, 1e-6) << name << ": not at cyclic steady state";
      }
      for (int k = 0; k < kNumKpis; ++k) {
        const Richardson r = Analyse(runs, k);
        // Standard deviation of a grid difference: sqrt(2) x the largest
        // standard error of the KPI means.
        double noise = 0.0;
        for (const auto& run : runs) noise = std::max(noise, run.noise[k]);
        const double sigma_d = std::sqrt(2.0) * noise;
        // GUARD 1: the grid effect is far above the noise, so the differences
        // below measure the discretisation and nothing else.
        EXPECT_GT(std::abs(r.d.front()), 100.0 * sigma_d)
          << name << " " << KpiName(k) << ": grid effect " << r.d.front() << " vs noise " << sigma_d;
        EXPECT_GT(std::abs(r.d.back()), 5.0 * sigma_d)
          << name << " " << KpiName(k) << ": finest difference " << r.d.back() << " vs noise " << sigma_d;
        std::printf("  %-8s %-15s finest difference %.2e, sd of a difference %.1e\n", name, KpiName(k),
                    std::abs(r.d.back()), sigma_d);
        for (size_t i = 1; i < r.d.size(); ++i) {
          EXPECT_GT(r.d[i] * r.d[0], 0.0)
            << name << " " << KpiName(k) << ": oscillatory convergence (differences "
            << r.d[i - 1] << ", " << r.d[i] << ")";
          EXPECT_LT(std::abs(r.d[i]), std::abs(r.d[i - 1]))
            << name << " " << KpiName(k) << ": grid difference did not shrink ("
            << r.d[i - 1] << " -> " << r.d[i] << ")";
        }
      }
    }

  } // namespace

  // ===========================================================================
  // Upwind: monotone convergence at an observed order in [1, 2]
  // ===========================================================================

  TEST(CycleKpiGridConvergence, UpwindConvergesMonotonicallyAtFirstOrder) {
    CheckMonotoneConvergence(AdvectionScheme::Upwind, "upwind");
    const auto& runs = Study(AdvectionScheme::Upwind);
    for (int k = 0; k < kNumKpis; ++k) {
      const Richardson r = Analyse(runs, k);
      // Formal order 1, next term second order: see the file comment.
      EXPECT_GE(r.order, 1.0) << "upwind " << KpiName(k) << ": observed order " << r.order;
      EXPECT_LE(r.order, 2.0) << "upwind " << KpiName(k) << ": observed order " << r.order;
      // The coarsest grid is far from converged: the 10-cell KPI is off by
      // several times the finest-grid difference (the guard that the study
      // spans a real change, not rounding).
      EXPECT_GT(std::abs(runs.front().value[k] - runs.back().value[k]), 5.0 * std::abs(r.d.back()));
    }
  }

  // ===========================================================================
  // WENO and van Leer: monotone, faster than upwind, and inside its GCI band
  // ===========================================================================

  void CheckFasterThanUpwind(AdvectionScheme scheme, const char* name) {
    CheckMonotoneConvergence(scheme, name);
    const auto& up = Study(AdvectionScheme::Upwind);
    const auto& sc = Study(scheme);
    for (int k = 0; k < kNumKpis; ++k) {
      const Richardson ru = Analyse(up, k);
      const Richardson rs = Analyse(sc, k);
      EXPECT_GT(rs.order, ru.order)
        << name << " " << KpiName(k) << ": observed order " << rs.order
        << " is not above upwind's " << ru.order;
      // A higher-order interior reconstruction must also give a smaller error
      // at every grid: each finest-grid difference is below upwind's.
      for (size_t i = 0; i < rs.d.size(); ++i)
        EXPECT_LT(std::abs(rs.d[i]), std::abs(ru.d[i]))
          << name << " " << KpiName(k) << ": grid difference " << i << " is not below upwind's";
    }
  }

  // ===========================================================================
  // The GCI is a true error band: every scheme converges to the same KPIs
  // ===========================================================================
  //
  // Three schemes with different truncation errors must converge to the same
  // grid-independent KPIs. The upwind GCI band around its 80-cell value is
  // Roache's estimate of where that limit is; the higher-order 80-cell values,
  // which are about ten times closer to it, must fall inside the band. And the
  // two higher-order schemes must agree with each other within the sum of
  // their own bands.

  //
  // (One test, so that ctest, which runs every test in its own process, runs
  // each four-grid study once.)

  TEST(CycleKpiGridConvergence, HigherOrderSchemesConvergeFasterToTheSameKpis) {
    CheckFasterThanUpwind(AdvectionScheme::WENO, "WENO");
    CheckFasterThanUpwind(AdvectionScheme::VanLeer, "van Leer");
    const auto& up = Study(AdvectionScheme::Upwind);
    const auto& weno = Study(AdvectionScheme::WENO);
    const auto& vl = Study(AdvectionScheme::VanLeer);
    for (int k = 0; k < kNumKpis; ++k) {
      const Richardson ru = Analyse(up, k);
      const Richardson rw = Analyse(weno, k);
      const Richardson rv = Analyse(vl, k);
      for (const auto* r : {&ru, &rw, &rv}) {
        const auto& runs = (r == &ru) ? up : (r == &rw) ? weno : vl;
        std::printf("  %-15s %-8s f = %.6g %.6g %.6g %.6g  p = %.2f  GCI = %.2e  f_ext = %.6g\n", KpiName(k),
                    r == &ru ? "upwind" : r == &rw ? "WENO" : "vanLeer", runs[0].value[k], runs[1].value[k],
                    runs[2].value[k], runs[3].value[k], r->order, r->gci, r->extrapolated);
      }
      const double band_up = ru.gci * std::abs(ru.fine);
      EXPECT_LT(std::abs(rw.fine - ru.fine), band_up)
        << KpiName(k) << ": WENO at 80 cells (" << rw.fine << ") is outside the upwind GCI band "
        << ru.fine << " +- " << band_up;
      EXPECT_LT(std::abs(rv.fine - ru.fine), band_up)
        << KpiName(k) << ": van Leer at 80 cells (" << rv.fine << ") is outside the upwind GCI band "
        << ru.fine << " +- " << band_up;
      // The Richardson-extrapolated upwind value is closer to the higher-order
      // solution than the 80-cell upwind value itself.
      EXPECT_LT(std::abs(rw.fine - ru.extrapolated), std::abs(rw.fine - ru.fine))
        << KpiName(k) << ": Richardson extrapolation of upwind does not move toward WENO";
      const double band_hi = rw.gci * std::abs(rw.fine) + rv.gci * std::abs(rv.fine);
      EXPECT_LT(std::abs(rw.fine - rv.fine), band_hi)
        << KpiName(k) << ": WENO (" << rw.fine << ") and van Leer (" << rv.fine
        << ") disagree by more than their GCI bands (" << band_hi << ")";
    }
  }

} // namespace adsim
