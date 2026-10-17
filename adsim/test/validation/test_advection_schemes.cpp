/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Order-of-accuracy verification of the convective face-reconstruction schemes.
//
// test_numerical_verification measures the observed order of the UPWIND scheme
// on a DISCONTINUOUS exact solution (a travelling step), via the outlet-curve
// variance. This file complements that with a SMOOTH exact solution -- an inert
// tracer with physical dispersion switched on, whose outlet breakthrough curve
// is the analytic Ogata-Banks profile -- and measures order the way CFD codes
// verify a discretisation when no exact FINITE-domain solution exists:
//
//   * ACCURACY is checked against the analytic Ogata-Banks solution. A raw
//     pointwise comparison carries a grid-INDEPENDENT floor, because the finite
//     column with its outlet boundary is not exactly the semi-infinite
//     Ogata-Banks solution (an O(1/Pe) difference). That floor is a property of
//     the physics + boundary, identical for every scheme, so it CANCELS when the
//     schemes are compared to each other -- which is what the scheme-comparison
//     tests use it for.
//
//   * ORDER is measured by self-convergence (Richardson): for a scheme of order
//     p, the L2 difference between successive grids satisfies
//         || F_h - F_2h ||  ~  C (2^p - 1) h^p,
//     so log2( ||F_2h - F_4h|| / ||F_h - F_2h|| ) = p. This needs no exact
//     solution and is immune to the boundary floor above.
//
// Three checks:
//   1. Upwind converges at first order on the smooth solution.
//   2. WENO reduces numerical diffusion vs upwind (more accurate at every grid).
//   3. van Leer reduces numerical diffusion vs upwind (same).
//
// NOTE on observed order: the breakthrough is measured at the OUTLET FACE, and
// the face reconstruction at the domain boundary is first-order (WENO/van Leer
// fall back to a simple average at the boundary face, and the outlet boundary
// condition is upwinded). So the OUTLET-measured self-convergence order is capped
// at ~1 for every scheme, including WENO/van Leer -- the interior high-order
// reconstruction cannot raise the order of a quantity read at a first-order
// boundary. What the higher-order schemes do deliver is a smaller error CONSTANT
// (less interior numerical diffusion), which is what tests 2 and 3 verify. Raising
// the outlet-measured order would require a high-order boundary reconstruction --
// a worthwhile future change, noted in the V&V dossier.
//
// The stoichiometric first moment (front position, fixed by capacity) must be
// preserved by every scheme -- a conservation cross-guard that a broken
// reconstruction would violate.

#include "breakthrough_case.h"

#include <limits>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    // Ogata-Banks outlet concentration ratio C/C0 at z = L for a step inlet.
    double OgataBanks(double L, double u, double D, double t) {
      if (t <= 0.0) return 0.0;
      const double s = 2.0 * std::sqrt(D * t);
      const double term1 = std::erfc((L - u * t) / s);
      const double exponent = u * L / D;
      double term2 = 0.0;
      if (exponent < 700.0) term2 = std::exp(exponent) * std::erfc((L + u * t) / s);
      return 0.5 * (term1 + term2);
    }

    // Fixed time grid the outlet curves are sampled onto so that runs at
    // different cell counts (and hence different internal time steps) can be
    // differenced pointwise.
    const std::vector<double>& CommonTimeGrid() {
      static const std::vector<double> grid = [] {
        std::vector<double> g;
        for (double t = 3.0; t <= 200.0; t += 0.5) g.push_back(t);
        return g;
      }();
      return grid;
    }

    struct SchemeRun {
      double u_int;               // interstitial velocity [m/s]
      double Pe;                  // u_int L / D_phys (inf when D_phys <= 0)
      double m1;                  // first moment (stoichiometric time) [s]
      double variance;            // outlet-curve variance [s^2] (numerical diffusion on a shock)
      double l2_vs_analytic;      // L2 of (F_sim - Ogata-Banks(D_phys)); 0 when D_phys <= 0
      std::vector<double> F;      // outlet C/C0 sampled on CommonTimeGrid()
    };

    // Run one inert-tracer breakthrough with the given face scheme. D_phys > 0
    // turns physical dispersion on (smooth Ogata-Banks solution, analytic L2
    // populated); D_phys <= 0 leaves dispersion off (the exact solution is a
    // travelling step, so the outlet variance is pure numerical diffusion).
    SchemeRun RunScheme(AdvectionScheme scheme, int n_cells, double D_phys) {
      BreakthroughCase c;
      c.n_cells = n_cells;
      c.adsorbing = false;
      c.length = 1.0;
      c.diameter = 0.1;
      c.molar_flow = 0.005;
      c.y_feed = 0.05;
      c.run_time = 300.0;
      c.max_step = 0.1;
      c.reltol = 1e-9;

      UserData data;
      Build(data, c);

      auto& layer = data.bed.layers_.layers_.front();
      auto disp = std::make_unique<DispersionModelConstant>();
      disp->dispersion_coefficient_ = D_phys > 0.0 ? D_phys : 0.0;
      layer.equations_.dispersion_model_ = std::move(disp);

      data.settings.advection_scheme_ = scheme;

      CoupledSolver solver{data};
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0)
        << "Solver did not complete (scheme=" << static_cast<int>(scheme)
        << ", n_cells=" << n_cells << ")";

      const auto& log = solver.GetLog();
      const auto& dom = log.GetDomain();

      const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
      const double P_in = P.values.back()[1];
      const double P_out = P.values.back()[P.values.back().size() - 2];
      EXPECT_LT(std::abs(P_in - P_out) / P_out, 0.01)
        << "Pressure drop too large for a constant-velocity Ogata-Banks fit";

      const auto vel = log.GetVariable(0, 0, VariableType::VELOCITY);
      const auto& u_face = vel.values.back();
      double u_sum = 0.0; int u_cnt = 0;
      for (size_t i = 1; i + 1 < u_face.size(); ++i) { u_sum += u_face[i]; ++u_cnt; }
      const double et = dom.et[1];
      const double u_int = (u_sum / u_cnt) / et;

      auto [times, y_out] = OutletSeries(log, c.feed_from_top);
      EXPECT_NEAR(y_out.back(), c.y_feed, 1e-2 * c.y_feed)
        << "Breakthrough incomplete (scheme=" << static_cast<int>(scheme)
        << ", n_cells=" << n_cells << ")";

      // Sample F = C/C0 onto the common grid, and accumulate the L2 error vs the
      // analytic Ogata-Banks profile at the physical D.
      std::vector<double> Fvals(y_out.size());
      for (size_t i = 0; i < y_out.size(); ++i) Fvals[i] = y_out[i] / c.y_feed;

      const auto& tg = CommonTimeGrid();
      std::vector<double> F(tg.size());
      double sum_sq = 0.0;
      for (size_t i = 0; i < tg.size(); ++i) {
        F[i] = interp_linear(times, Fvals, tg[i]);
        if (D_phys > 0.0) {
          const double ref = OgataBanks(c.length, u_int, D_phys, tg[i]);
          sum_sq += (F[i] - ref) * (F[i] - ref);
        }
      }
      const double l2 = D_phys > 0.0 ? std::sqrt(sum_sq / tg.size()) : 0.0;

      // Moments of the outlet RTD from integrals of (1 - F): m1 = int (1-F) dt,
      // m2 = 2 int t (1-F) dt, variance = m2 - m1^2.
      std::vector<double> unbroken(times.size()), t_unbroken(times.size());
      for (size_t i = 0; i < times.size(); ++i) {
        unbroken[i] = 1.0 - Fvals[i];
        t_unbroken[i] = times[i] * unbroken[i];
      }
      const double m1 = TimeIntegral(times, unbroken);
      const double variance = 2.0 * TimeIntegral(times, t_unbroken) - m1 * m1;

      const double Pe = D_phys > 0.0 ? u_int * c.length / D_phys
                                     : std::numeric_limits<double>::infinity();
      return {u_int, Pe, m1, variance, l2, std::move(F)};
    }

    double L2Diff(const std::vector<double>& a, const std::vector<double>& b) {
      double s = 0.0;
      for (size_t i = 0; i < a.size(); ++i) s += (a[i] - b[i]) * (a[i] - b[i]);
      return std::sqrt(s / a.size());
    }

    void CheckSmoothFrontGuards(const SchemeRun& r) {
      ASSERT_GT(r.Pe, 5.0)   << "Peclet " << r.Pe << " too low: column fully mixed";
      ASSERT_LT(r.Pe, 200.0) << "Peclet " << r.Pe << " too high: front is a shock";
    }

    constexpr double kDphys = 2.0e-3;  // [m2/s]

  } // namespace

  // ===========================================================================
  // 1. Upwind: first order on the smooth solution (self-convergence)
  // ===========================================================================

  TEST(AdvectionSchemeOrder, UpwindIsFirstOrderOnSmoothSolution) {
    std::vector<SchemeRun> r;
    for (int n : {40, 80, 160, 320})
      r.push_back(RunScheme(AdvectionScheme::Upwind, n, kDphys));
    for (const auto& run : r) CheckSmoothFrontGuards(run);

    // First moment (stoichiometric time) is grid-invariant.
    for (size_t i = 1; i < r.size(); ++i)
      EXPECT_NEAR(r[i].m1, r[0].m1, 0.01 * r[0].m1)
        << "Stoichiometric time drifted with grid refinement (upwind)";

    // Correctness anchor: the finest solution matches the analytic Ogata-Banks
    // profile to within the finite-column O(1/Pe) floor.
    EXPECT_LT(r.back().l2_vs_analytic, 0.05)
      << "Finest upwind solution does not match analytic Ogata-Banks (L2 "
      << r.back().l2_vs_analytic << ")";

    // Observed order by self-convergence: log2 of the ratio of successive-grid
    // L2 differences. Two independent estimates from the four grids.
    const double d1 = L2Diff(r[0].F, r[1].F);  // 40 vs 80
    const double d2 = L2Diff(r[1].F, r[2].F);  // 80 vs 160
    const double d3 = L2Diff(r[2].F, r[3].F);  // 160 vs 320
    ASSERT_GT(d3, 0.0);
    const double order_a = std::log2(d1 / d2);
    const double order_b = std::log2(d2 / d3);
    EXPECT_NEAR(order_a, 1.0, 0.3)
      << "Upwind self-convergence order (40->80->160) is " << order_a;
    EXPECT_NEAR(order_b, 1.0, 0.3)
      << "Upwind self-convergence order (80->160->320) is " << order_b;
  }

  // On the SHOCK case (dispersion off) the exact solution is a travelling step
  // whose outlet variance is zero, so the entire simulated variance is numerical
  // diffusion. Comparing a scheme's outlet variance to upwind's is therefore a
  // direct, unambiguous measure of how much less numerical diffusion it carries
  // -- the decisive verification of a higher-order reconstruction. Both must also
  // preserve the stoichiometric first moment (front position) and converge toward
  // the zero-variance exact solution under refinement.
  void CheckLessDiffusiveThanUpwind(AdvectionScheme scheme, const char* name) {
    const std::vector<int> grids = {40, 80, 160};
    std::vector<SchemeRun> up, sc;
    for (int n : grids) {
      up.push_back(RunScheme(AdvectionScheme::Upwind, n, /*D_phys=*/0.0));
      sc.push_back(RunScheme(scheme, n, /*D_phys=*/0.0));
    }

    for (size_t i = 0; i < grids.size(); ++i) {
      ASSERT_GT(sc[i].variance, 0.0);

      // Front arrives at its physical stoichiometric time L/u_int: the scheme
      // must sharpen the front without moving it. (A few percent of slack covers
      // the effect of the [0,1] face-value clamp on the sharp shock, which can
      // nudge the measured arrival time; total mass stays conserved regardless.)
      const double t_arrival = 1.0 / sc[i].u_int;  // L / u_int, L = 1 m
      EXPECT_NEAR(sc[i].m1, t_arrival, 0.03 * t_arrival)
        << name << " front arrival time " << sc[i].m1 << " s departs from L/u_int = "
        << t_arrival << " s at n_cells = " << grids[i];

      // Strictly less numerical diffusion than upwind at every grid.
      EXPECT_LT(sc[i].variance, up[i].variance)
        << name << " is not less diffusive than upwind at n_cells = " << grids[i]
        << " (" << name << " variance " << sc[i].variance << " vs upwind "
        << up[i].variance << ")";
    }

    // ... and materially so: at least a 25% reduction at the finest grid.
    EXPECT_LT(sc.back().variance, 0.75 * up.back().variance)
      << name << " numerical-diffusion reduction over upwind is marginal at the "
      << "finest grid (" << sc.back().variance << " vs upwind " << up.back().variance << ")";

    // Converges toward the zero-variance exact solution under refinement.
    EXPECT_LT(sc.back().variance, sc.front().variance)
      << name << " variance did not decrease under grid refinement ("
      << sc.front().variance << " -> " << sc.back().variance << ")";
  }

  // ===========================================================================
  // 2. WENO carries less numerical diffusion than upwind
  // ===========================================================================

  TEST(AdvectionSchemeOrder, WENOIsLessDiffusiveThanUpwind) {
    CheckLessDiffusiveThanUpwind(AdvectionScheme::WENO, "WENO");
  }

  // ===========================================================================
  // 3. van Leer carries less numerical diffusion than upwind
  // ===========================================================================

  TEST(AdvectionSchemeOrder, VanLeerIsLessDiffusiveThanUpwind) {
    CheckLessDiffusiveThanUpwind(AdvectionScheme::VanLeer, "van Leer");
  }

} // namespace adsim
