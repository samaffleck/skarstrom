/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Analytical validation case C: retarded advection-dispersion for a LINEAR
// isotherm.
//
// Case A (test_analytic_dispersion) pins the physical dispersion term for an
// INERT tracer against Ogata-Banks. The stoichiometric-time test
// (test_analytic_transport) pins the front POSITION for an adsorbing Henry
// species via a dispersion-invariant first moment. Neither exercises what this
// test does: how a linear isotherm RETARDS a dispersion-broadened front -- the
// coupling of the isotherm slope with the advection AND dispersion terms.
//
// For a linear (Henry) isotherm, q* = [k exp(Q/RT) R_g T] C, the species balance
// with dispersion reduces to the advection-dispersion equation with a
// retardation factor R:
//
//     R dC/dt = D d2C/dz2 - u_int dC/dz,
//     R = 1 + (rho_b/et) dq*/dC = 1 + (rho_b/et) k exp(Q/RT) R_g T,
//
// i.e. Ogata-Banks with velocity u* = u_int/R and dispersion D* = D/R. Two exact
// consequences follow for the outlet breakthrough curve, and they are the sharp
// signatures this test checks:
//
//   * POSITION:  mean_ads / mean_inert = R.
//   * WIDTH:     var_ads  / var_inert  = R^2.
//
// (The second follows from var = 2 D* L / u*^3 = R^2 * (2 D L / u^3).)
//
// Crucially these are checked as RATIOS against an inert front run at the SAME
// dispersion, flow and grid. The retarded and inert fronts share the same
// Peclet number (u*L/D* = u_int L/D), so their finite-column O(1/Pe) corrections
// and the upwind scheme's numerical diffusion cancel in the ratio, leaving R and
// R^2 as clean targets that do not depend on the accuracy of any absolute
// semi-infinite closed form. The isotherm-derived R (>~5 here) is large, so the
// test genuinely exercises the coupling and cannot be satisfied by inert physics.
//
// The mass-transfer coefficient is set fast so LDF resistance does not add a
// second, competing source of spreading; what is left is dispersion acting on a
// retarded front.

#include "breakthrough_case.h"

namespace adsim {

  using namespace adsim::validation;

  namespace {

    // First two central moments of the breakthrough curve F (0 -> 1), treating
    // its increments dF as the residence-time distribution. Returns {mean,
    // variance} in time units, normalised by the total rise.
    std::pair<double, double> CurveMoments(const std::vector<double>& t,
                                           const std::vector<double>& F) {
      double norm = 0.0, m1 = 0.0;
      for (size_t i = 1; i < t.size(); ++i) {
        const double dF = F[i] - F[i - 1];
        const double tm = 0.5 * (t[i] + t[i - 1]);
        norm += dF;
        m1 += dF * tm;
      }
      if (norm <= 0.0) return {0.0, 0.0};
      m1 /= norm;
      double m2 = 0.0;
      for (size_t i = 1; i < t.size(); ++i) {
        const double dF = F[i] - F[i - 1];
        const double tm = 0.5 * (t[i] + t[i - 1]);
        m2 += dF * (tm - m1) * (tm - m1);
      }
      return {m1, m2 / norm};
    }

    struct FrontMoments {
      double u_int;   // interstitial velocity [m/s]
      double mean;    // breakthrough-time first moment [s]
      double var;     // breakthrough-time variance [s^2]
      double Pe;      // u_int L / D_phys (== u* L / D* for the retarded case)
      double D_num;   // upwind numerical diffusion u_int*dx/2 [m2/s]
    };

    // Build, run and reduce one breakthrough case (inert or adsorbing) at a fixed
    // physical dispersion. Adsorbing cases get a fast LDF so mass-transfer
    // resistance is negligible against dispersion.
    FrontMoments RunFront(int n_cells, bool adsorbing, double D_phys) {
      BreakthroughCase c;
      c.n_cells = n_cells;
      c.adsorbing = adsorbing;
      c.length = 1.0;        // [m]
      c.diameter = 0.1;      // [m]
      c.molar_flow = 0.005;  // [mol/s]
      c.y_feed = 0.02;       // dilute -> isotherm stays linear
      c.run_time = 900.0;    // long enough for a retarded (R~5) breakthrough
      c.max_step = 0.2;
      c.reltol = 1e-9;

      UserData data;
      Build(data, c);

      auto& layer = data.bed.layers_.layers_.front();
      auto disp = std::make_unique<DispersionModelConstant>();
      disp->dispersion_coefficient_ = D_phys;
      layer.equations_.dispersion_model_ = std::move(disp);

      if (adsorbing) {
        auto& adsorbent = data.adsorbent_directory.Get(layer.composition_.begin()->first);
        auto ldf = std::make_unique<LinearDrivingForce>();
        auto mtc = std::make_unique<MassTransferCoefficientConstant>();
        mtc->constant_mtc_ = 50.0;  // [1/s], fast: LDF spread << dispersion
        ldf->mtc_model_ = std::move(mtc);
        adsorbent.adsorption_source_.mass_transfer_models_["Oxygen"] = std::move(ldf);
      }

      CoupledSolver solver{data};
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0)
        << "Solver did not complete (adsorbing=" << adsorbing
        << ", n_cells=" << n_cells << ")";

      const auto& log = solver.GetLog();
      const auto& dom = log.GetDomain();

      // Essentially isobaric so u is constant along the column.
      const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
      const double P_in = P.values.back()[1];
      const double P_out = P.values.back()[P.values.back().size() - 2];
      EXPECT_LT(std::abs(P_in - P_out) / P_out, 0.01)
        << "Pressure drop too large for a constant-velocity fit";

      const auto vel = log.GetVariable(0, 0, VariableType::VELOCITY);
      const auto& u_face = vel.values.back();
      double u_sum = 0.0; int u_cnt = 0;
      for (size_t i = 1; i + 1 < u_face.size(); ++i) { u_sum += u_face[i]; ++u_cnt; }
      const double et = dom.et[1];
      const double u_int = (u_sum / u_cnt) / et;

      auto [times, y_out] = OutletSeries(log, c.feed_from_top);
      EXPECT_NEAR(y_out.back(), c.y_feed, 1e-2 * c.y_feed)
        << "Breakthrough incomplete (adsorbing=" << adsorbing
        << ", n_cells=" << n_cells << ")";
      std::vector<double> F(times.size());
      for (size_t i = 0; i < times.size(); ++i) F[i] = y_out[i] / c.y_feed;
      auto [mean, var] = CurveMoments(times, F);

      const double dx = dom.cell_width[1];
      return {u_int, mean, var, u_int * c.length / D_phys, u_int * dx / 2.0};
    }

    // Retardation factor from the Henry isotherm installed by Build().
    double RetardationFactor(double T) {
      const double rho_b = (1.0 - kCaseVoidFraction) * kCaseSolidDensity;
      const double henry_slope = kCaseHenryK * std::exp(kCaseHenryQ / (kR * T)) * kR * T; // [m3/kg]
      return 1.0 + (rho_b / kCaseVoidFraction) * henry_slope;
    }

  } // namespace

  // ===========================================================================
  // Retarded advection-dispersion: position retards by R, width by R^2
  // ===========================================================================

  TEST(AnalyticRetardedDispersion, LinearIsothermRetardsFrontPositionAndWidth) {
    const double D_phys = 3.0e-3;  // [m2/s]
    const double R = RetardationFactor(298.0);

    // GUARD: retardation must dominate, or this degenerates into the inert test.
    ASSERT_GT(R, 3.0) << "Retardation factor " << R << " too small";

    struct GridRatios {
      int n_cells;
      double mean_ratio;  // mean_ads / mean_inert   (expect R)
      double var_ratio;   // var_ads  / var_inert    (expect R^2)
    };
    std::vector<GridRatios> grids;

    for (int n_cells : {80, 160}) {
      const FrontMoments inert = RunFront(n_cells, /*adsorbing=*/false, D_phys);
      const FrontMoments ads   = RunFront(n_cells, /*adsorbing=*/true,  D_phys);

      // GUARDS on the inert front (identical geometry/flow, so they carry over):
      // dispersion dominates numerical diffusion, and the front is
      // dispersion-broadened rather than a shock or fully mixed.
      ASSERT_GT(D_phys, 3.0 * inert.D_num)
        << "Physical dispersion does not dominate numerical diffusion at n_cells = "
        << n_cells << " (D_phys=" << D_phys << ", D_num=" << inert.D_num << ")";
      ASSERT_GT(inert.Pe, 5.0)   << "Peclet " << inert.Pe << " too low";
      ASSERT_LT(inert.Pe, 200.0) << "Peclet " << inert.Pe << " too high";

      // The two runs must see the same interstitial velocity, else the ratio
      // targets R / R^2 would not hold.
      ASSERT_NEAR(ads.u_int, inert.u_int, 0.01 * inert.u_int)
        << "Adsorbing and inert runs disagree on interstitial velocity";

      grids.push_back({n_cells, ads.mean / inert.mean, ads.var / inert.var});
    }

    // -----------------------------------------------------------------------
    // 1. POSITION retards by R. Both fronts carry the same O(1/Pe) finite-column
    //    offset, so it cancels and the mean ratio is R to a tight tolerance.
    for (const auto& g : grids) {
      EXPECT_NEAR(g.mean_ratio, R, 0.03 * R)
        << "mean_ads/mean_inert = " << g.mean_ratio << " != R = " << R
        << " at n_cells = " << g.n_cells << " (position is not retarded by the "
        << "isotherm-derived factor)";
    }

    // -----------------------------------------------------------------------
    // 2. WIDTH retards by R^2. This is the dispersion<->isotherm coupling: the
    //    dispersive spread of the retarded front is R^2 times the inert spread.
    //    Numerical diffusion and finite-Pe corrections cancel in the ratio, so
    //    the target is the clean R^2 with a modest tolerance for the residual
    //    O(1/Pe^2) term and moment-integration error.
    for (const auto& g : grids) {
      EXPECT_NEAR(g.var_ratio, R * R, 0.08 * R * R)
        << "var_ads/var_inert = " << g.var_ratio << " != R^2 = " << R * R
        << " at n_cells = " << g.n_cells << " (dispersive width is not retarded "
        << "by R^2; the isotherm<->dispersion coupling is wrong)";
    }

    // 3. The ratios must be grid-stable (they are exact, so refinement must not
    //    move them appreciably).
    EXPECT_NEAR(grids.back().var_ratio, grids.front().var_ratio,
                0.05 * grids.front().var_ratio)
      << "Variance ratio drifts with grid refinement (coarse "
      << grids.front().var_ratio << ", fine " << grids.back().var_ratio << ")";
  }

} // namespace adsim
