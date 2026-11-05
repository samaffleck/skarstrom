/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Analytical validation of a MASS-TRANSFER-LIMITED front: a linear (Henry)
// isotherm with linear-driving-force kinetics, no axial dispersion, isothermal,
// against the Anzelius / Thomas closed form for the WHOLE outlet curve.
//
// The other transport tests pin the front position (first moment, which is
// independent of the kinetics) or run with a fast LDF so the front is set by
// dispersion. Here the LDF rate is what shapes the front, and the shape is
// compared point by point.
//
// Model. For a dilute adsorbate in an inert carrier at constant pressure and
// temperature, with interstitial velocity u, total voidage eps (= eps_b, the
// particles are non-porous here), bulk density rho_b and Henry slope
// K = dq*/dc [m3/kg]:
//
//     dc/dt + u dc/dz + (rho_b/eps) dq/dt = 0,     dq/dt = k (K c - q),
//
// with c(0, t) = c0 and a clean bed. In the frame tau = k (t - z/u) and with the
// dimensionless distance x = k (R - 1) z / u, where R = 1 + (rho_b/eps) K is
// the retardation factor, this is dc/dx = -dQ/dtau, dQ/dtau = c - Q (Q = q/K),
// whose solution (Anzelius 1926; Thomas 1944; Ruthven 1984, ch. 8) is
//
//     c/c0 = J(x, tau) = 1 - int_0^x exp(-tau - s) I0(2 sqrt(tau s)) ds,
//
// or equivalently, since its Laplace transform in tau is exp(-x s / (1 + s)) / s,
//
//     J(x, tau) = exp(-x) + int_0^tau h(x, sigma) dsigma,
//     h(x, tau) = exp(-x - tau) sqrt(x / tau) I1(2 sqrt(x tau)).
//
// J is evaluated here by Simpson quadrature of h, with I0 and I1 summed from
// their ascending series (nothing from core/src). Its cumulants in tau are
// x, 2x and 6x, so the outlet curve has mean L R / u, variance
// 2 (R - 1) L / (u k) and skewness 2.12 / sqrt(x): a clearly asymmetric curve.
//
// The discretisation. The solver's first-order upwind scheme turns the column
// into n tanks in series, each with its own LDF. The Laplace transfer function
// of that semi-discrete system is exactly
//
//     G_n(s) = (1 + (dx/u) g(s))^(-L/dx),   g(s) = s (1 + (R - 1) / (1 + s/k)),
//
// so ln G_n = -(L/u) g + (L dx / 2u^2) g^2 + O(dx^2). Two exact consequences
// are checked:
//
//   * VARIANCE. The second cumulant of G_n is EXACTLY
//         var_n = 2 (R - 1) L / (u k) + L dx R^2 / u^2
//     (the second term is the u dx / 2 numerical dispersion acting on a front
//     retarded by R) on every grid, with no higher-order term.
//
//   * SHAPE. Since g^2 F^(s) = u^2 d2F^/dL2 for the exact curve, the leading
//     error of the whole outlet curve is known in closed form:
//         F_n(t) = J + (L dx / 2) d2F/dz2 |_(z = L, t) + O(dx^2),
//     with d2F/dz2 computed analytically from J. So max |F_n - J| must fall at
//     FIRST order, and the residual after subtracting this correction at
//     SECOND order. That is a far sharper statement than any single tolerance:
//     it holds only if the LDF rate, the Henry capacity, the velocity and the
//     upwind flux all enter exactly as the model says.
//
// Neglected terms (each estimated, and well below the second-order residual):
// the sorption effect (the carrier velocity drops as the adsorbate leaves the
// gas), an O(y_feed) distortion of about 0.3 y_feed in F, so y_feed = 1e-4 makes
// it about 3e-5; the pressure drop (about 70 Pa on 1 bar, so u varies by
// < 0.1 %); and time integration at reltol = 1e-9.
//
// GUARDS: the curve is far from the equilibrium step and from a symmetric
// (Gaussian) curve of the same mean and variance (max difference about 0.04,
// some 250x the finest-grid residual), so neither the kinetics nor the
// asymmetry can be missing; and the LDF variance dominates the numerical one.

#include "breakthrough_case.h"

#include <algorithm>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    // -------------------------------------------------------------------------
    // Reference solution, written from the textbook forms above
    // -------------------------------------------------------------------------

    // exp(-z) I0(z) and exp(-z) I1(z) from the ascending series
    //   I0(z) = sum (z/2)^(2m) / (m!)^2,  I1(z) = sum (z/2)^(2m+1) / (m! (m+1)!).
    // All terms are positive, so the sum is accurate to rounding for the
    // arguments used here (z < 60, where the largest term is ~1e25).
    double ScaledI0(double z) {
      double term = 1.0, sum = 1.0;
      for (int m = 1; m < 500; ++m) {
        term *= (0.25 * z * z) / (double(m) * double(m));
        sum += term;
        if (term < 1e-17 * sum) break;
      }
      return sum * std::exp(-z);
    }

    double ScaledI1(double z) {
      double term = 0.5 * z, sum = term;
      for (int m = 1; m < 500; ++m) {
        term *= (0.25 * z * z) / (double(m) * double(m + 1));
        sum += term;
        if (term < 1e-17 * sum) break;
      }
      return sum * std::exp(-z);
    }

    // exp(-x - tau) I_n(2 sqrt(x tau)), written as
    // exp(-(sqrt x - sqrt tau)^2) * [exp(-z) I_n(z)] so nothing overflows.
    double DecayedI0(double x, double tau) {
      const double z = 2.0 * std::sqrt(x * tau);
      const double d = std::sqrt(x) - std::sqrt(tau);
      return std::exp(-d * d) * ScaledI0(z);
    }
    double DecayedI1(double x, double tau) {
      const double z = 2.0 * std::sqrt(x * tau);
      const double d = std::sqrt(x) - std::sqrt(tau);
      return std::exp(-d * d) * ScaledI1(z);
    }

    // Outlet residence-time density in tau, dJ/dtau for tau > 0.
    double ThomasH(double x, double tau) {
      if (tau <= 0.0) return x * std::exp(-x);  // the tau -> 0+ limit
      return std::sqrt(x / tau) * DecayedI1(x, tau);
    }

    // J(x, tau) = exp(-x) + int_0^tau h(x, s) ds (Simpson, 2000 panels; h is
    // smooth, so the quadrature error is far below 1e-10). J = 0 for tau < 0:
    // the feed has not reached z yet.
    double ThomasJ(double x, double tau) {
      if (tau < 0.0) return 0.0;
      if (tau == 0.0) return std::exp(-x);
      const int n = 2000;
      const double hs = tau / n;
      double sum = ThomasH(x, 0.0) + ThomasH(x, tau);
      for (int i = 1; i < n; ++i) sum += (i % 2 ? 4.0 : 2.0) * ThomasH(x, i * hs);
      return std::exp(-x) + sum * hs / 3.0;
    }

    // The other textbook form, 1 - int_0^x exp(-tau - s) I0(2 sqrt(tau s)) ds,
    // used only to cross-check ThomasJ.
    double ThomasJByDistance(double x, double tau) {
      const int n = 2000;
      const double hs = x / n;
      double sum = DecayedI0(tau, 0.0) + DecayedI0(tau, x);
      for (int i = 1; i < n; ++i) sum += (i % 2 ? 4.0 : 2.0) * DecayedI0(tau, i * hs);
      return 1.0 - sum * hs / 3.0;
    }

    // Second partial derivatives of J for tau > 0, from J_x = -E I0 and
    // J_tau = E sqrt(x/tau) I1 (E = exp(-x - tau), z = 2 sqrt(x tau)) with
    // I0' = I1 and I1'(z) = I0 - I1/z:
    //   J_xx     = E [ I0 - sqrt(tau/x) I1 ]
    //   J_xtau   = E [ I0 - sqrt(x/tau) I1 ]
    //   J_tautau = E [ (x/tau) I0 - sqrt(x/tau) I1 - sqrt(x) tau^(-3/2) I1 ]
    struct SecondDerivatives { double xx, xt, tt; };
    SecondDerivatives ThomasSecondDerivatives(double x, double tau) {
      const double EI0 = DecayedI0(x, tau);
      const double EI1 = DecayedI1(x, tau);
      return {EI0 - std::sqrt(tau / x) * EI1,
              EI0 - std::sqrt(x / tau) * EI1,
              (x / tau) * EI0 - std::sqrt(x / tau) * EI1
                - std::sqrt(x) * std::pow(tau, -1.5) * EI1};
    }

    // The exact outlet curve F(z, t) = J(a z, k (t - z/u)), a = k (R - 1) / u,
    // and its second derivative in z at fixed t:
    //   F_zz = a^2 J_xx - 2 a (k/u) J_xtau + (k/u)^2 J_tautau.
    struct ThomasFront {
      double k;      // LDF rate [1/s]
      double R;      // retardation factor [-]
      double u;      // interstitial velocity [m/s]
      double a() const { return k * (R - 1.0) / u; }  // [1/m]
      double F(double z, double t) const { return ThomasJ(a() * z, k * (t - z / u)); }
      double Fzz(double z, double t) const {
        const double tau = k * (t - z / u);
        if (tau <= 1e-9) return 0.0;  // the feed has not arrived (the e^-x jump is ~1e-5 here)
        const auto d = ThomasSecondDerivatives(a() * z, tau);
        const double b = k / u;
        return a() * a() * d.xx - 2.0 * a() * b * d.xt + b * b * d.tt;
      }
    };

    // Retardation factor of the Henry isotherm Build() installs (ep = 0, so the
    // total voidage is the bed voidage).
    double HenryRetardation(double T) {
      const double rho_b = (1.0 - kCaseVoidFraction) * kCaseSolidDensity;
      const double K = kCaseHenryK * std::exp(kCaseHenryQ / (kR * T)) * kR * T;  // dq*/dc [m3/kg]
      return 1.0 + (rho_b / kCaseVoidFraction) * K;
    }

    // -------------------------------------------------------------------------
    // One run of the plain breakthrough case with a slow LDF
    // -------------------------------------------------------------------------

    constexpr double kLdfRate = 0.1;     // [1/s]: x = k (R-1) L/u ~ 12 transfer units
    constexpr double kYFeed = 1.0e-4;    // dilute: sorption effect ~3e-5 in F
    constexpr double kRunTime = 450.0;   // [s]: mean + 6 sigma of the outlet curve

    struct LdfRun {
      std::vector<double> t, F;  // outlet y / y_feed
      double u = 0.0;            // interstitial velocity [m/s]
      double dx = 0.0;           // [m]
      double dP_rel = 0.0;       // pressure drop / P
    };

    LdfRun RunLdfFront(int n_cells) {
      BreakthroughCase c;
      c.n_cells = n_cells;
      c.adsorbing = true;   // Henry isotherm on Oxygen
      c.y_feed = kYFeed;
      c.run_time = kRunTime;
      c.max_step = 0.5;
      c.reltol = 1e-9;

      UserData data;
      Build(data, c);
      // The reference carries the exact upwind error (the default is WENO).
      data.settings.advection_scheme_ = AdvectionScheme::Upwind;
      data.settings.max_log_points_per_step_ = 450;  // one snapshot per second

      auto& layer = data.bed.layers_.layers_.front();
      auto& adsorbent = data.adsorbent_directory.Get(layer.composition_.begin()->first);
      auto ldf = std::make_unique<LinearDrivingForce>();
      auto mtc = std::make_unique<MassTransferCoefficientConstant>();
      mtc->constant_mtc_ = kLdfRate;
      ldf->mtc_model_ = std::move(mtc);
      adsorbent.adsorption_source_.mass_transfer_models_["Oxygen"] = std::move(ldf);

      CoupledSolver solver{data};
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0)
        << "Solver did not complete (n_cells = " << n_cells << ")";

      const auto& log = solver.GetLog();
      const auto& dom = log.GetDomain();

      LdfRun run;
      const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
      const auto& Pf = P.values.back();
      run.dP_rel = std::abs(Pf[1] - Pf[Pf.size() - 2]) / Pf[Pf.size() - 2];

      // Interstitial velocity: mean superficial velocity over the interior faces
      // at the end of the run, divided by the voidage.
      const auto vel = log.GetVariable(0, 0, VariableType::VELOCITY);
      const auto& uf = vel.values.back();
      double u_sum = 0.0;
      int u_cnt = 0;
      for (size_t i = 1; i + 1 < uf.size(); ++i) { u_sum += uf[i]; ++u_cnt; }
      run.u = (u_sum / u_cnt) / dom.et[1];
      run.dx = dom.cell_width[1];

      auto [times, y_out] = OutletSeries(log, c.feed_from_top);
      run.t = times;
      run.F.resize(times.size());
      for (size_t i = 0; i < times.size(); ++i) run.F[i] = y_out[i] / kYFeed;
      return run;
    }

    // Mean and variance of a complete breakthrough curve from its survival
    // function: mean = int (1 - F) dt, E[t^2] = int 2 t (1 - F) dt.
    std::pair<double, double> CurveMoments(const std::vector<double>& t,
                                           const std::vector<double>& F) {
      std::vector<double> s1(t.size()), s2(t.size());
      for (size_t i = 0; i < t.size(); ++i) {
        s1[i] = 1.0 - F[i];
        s2[i] = 2.0 * t[i] * (1.0 - F[i]);
      }
      const double m1 = TimeIntegral(t, s1);
      return {m1, TimeIntegral(t, s2) - m1 * m1};
    }

    double ObservedOrder(double e_coarse, double e_fine) {
      return std::log(e_coarse / e_fine) / std::log(2.0);
    }

  } // namespace

  // ===========================================================================
  // The reference solution is right
  // ===========================================================================

  // Checks the quadrature of J against the other textbook form and against two
  // exact identities, and the analytic second derivatives against finite
  // differences of J, so a slip in the reference cannot pass for a solver
  // defect (or hide one). No solver run.
  TEST(LdfThomasFront, ReferenceSolutionSatisfiesItsIdentities) {
    for (const auto& [x, tau] : std::vector<std::pair<double, double>>{
           {12.0, 3.0}, {12.0, 12.0}, {12.0, 25.0}, {6.0, 9.0}, {3.0, 0.5}}) {
      // The two integral forms agree.
      EXPECT_NEAR(ThomasJ(x, tau), ThomasJByDistance(x, tau), 1e-10) << "x=" << x << " tau=" << tau;
      // J(a, b) + J(b, a) = 1 + exp(-a - b) I0(2 sqrt(a b))  (Ruthven 1984, eq. 8.21).
      EXPECT_NEAR(ThomasJ(x, tau) + ThomasJ(tau, x), 1.0 + DecayedI0(x, tau), 1e-10)
        << "x=" << x << " tau=" << tau;

      const double e = 1e-3;
      const auto d = ThomasSecondDerivatives(x, tau);
      const double fd_xx = (ThomasJ(x + e, tau) - 2.0 * ThomasJ(x, tau) + ThomasJ(x - e, tau)) / (e * e);
      const double fd_tt = (ThomasJ(x, tau + e) - 2.0 * ThomasJ(x, tau) + ThomasJ(x, tau - e)) / (e * e);
      const double fd_xt = (ThomasJ(x + e, tau + e) - ThomasJ(x + e, tau - e)
                            - ThomasJ(x - e, tau + e) + ThomasJ(x - e, tau - e)) / (4.0 * e * e);
      EXPECT_NEAR(d.xx, fd_xx, 1e-5) << "x=" << x << " tau=" << tau;
      EXPECT_NEAR(d.tt, fd_tt, 1e-5) << "x=" << x << " tau=" << tau;
      EXPECT_NEAR(d.xt, fd_xt, 1e-5) << "x=" << x << " tau=" << tau;
    }
    // At the moment the feed arrives the bed has had no time to adsorb, so the
    // gas is only depleted by the instantaneous uptake: J(x, 0) = exp(-x).
    EXPECT_NEAR(ThomasJByDistance(5.0, 0.0), std::exp(-5.0), 1e-12);
    // And the cumulants of the tau-distribution are x and 2x.
    const double x = 12.0;
    std::vector<double> tau, J;
    for (int i = 0; i <= 3000; ++i) { tau.push_back(0.02 * i); J.push_back(ThomasJ(x, 0.02 * i)); }
    const auto [mean, var] = CurveMoments(tau, J);
    EXPECT_NEAR(mean, x, 1e-4);
    EXPECT_NEAR(var, 2.0 * x, 1e-3);
  }

  // ===========================================================================
  // Linear isotherm + LDF: the whole outlet curve against Thomas
  // ===========================================================================

  TEST(LdfThomasFront, LinearIsothermOutletCurveMatchesThomasSolution) {
    const double T = 298.0;
    const double R = HenryRetardation(T);
    const double L = 1.0;

    struct Grid {
      int n_cells;
      double E;         // max |F_n - J|
      double residual;  // max |F_n - J - (L dx / 2) F_zz|
      double var_excess_rel;  // (var_n - L dx R^2/u^2) / var_LDF - 1
    };
    std::vector<Grid> grids;
    double guard_gauss = 0.0;  // max |F_finest - Gaussian with Thomas moments|

    for (int n_cells : {50, 100, 200}) {
      const LdfRun run = RunLdfFront(n_cells);
      ASSERT_FALSE(run.t.empty());
      const ThomasFront ref{kLdfRate, R, run.u};

      const double x_L = ref.a() * L;                          // transfer units
      const double var_ldf = 2.0 * (R - 1.0) * L / (run.u * kLdfRate);
      const double var_num = L * run.dx * R * R / (run.u * run.u);
      const double mean_ref = L * R / run.u;

      // --- the derivation's assumptions ------------------------------------
      ASSERT_LT(run.dP_rel, 2e-3) << "Pressure drop too large for a constant velocity";
      ASSERT_GT(run.F.back(), 1.0 - 1e-4) << "Breakthrough incomplete; the moments would be truncated";
      // GUARD: a genuinely kinetic front of a strongly retained species.
      ASSERT_GT(R, 3.0);
      ASSERT_GT(x_L, 5.0) << "Too few transfer units";
      ASSERT_LT(x_L, 30.0) << "Too many transfer units: the curve would be nearly Gaussian";
      // GUARD: the LDF spreading dominates the numerical spreading.
      ASSERT_LT(var_num, 0.2 * var_ldf) << "Numerical dispersion dominates at n = " << n_cells;

      // --- 1. moments: mean and the EXACT semi-discrete variance -----------
      const auto [mean, var] = CurveMoments(run.t, run.F);
      EXPECT_NEAR(mean, mean_ref, 2e-3 * mean_ref)
        << "First moment " << mean << " s != L R / u = " << mean_ref << " s at n = " << n_cells;
      const double var_excess_rel = (var - var_num) / var_ldf - 1.0;
      // Tolerance: the time quadrature of the moments (1 s sampling, ~1e-4) and
      // the O(y_feed) sorption effect (~1e-4); 0.2 % covers both with margin.
      EXPECT_NEAR(var - var_num, var_ldf, 2e-3 * var_ldf)
        << "Outlet variance minus the upwind contribution L dx R^2/u^2 = " << var - var_num
        << " s^2 != Thomas variance 2 (R-1) L/(u k) = " << var_ldf << " s^2 at n = " << n_cells;

      // --- 2. the curve, point by point ------------------------------------
      double E = 0.0, residual = 0.0, gauss = 0.0;
      const double sd = std::sqrt(var_ldf);
      for (size_t i = 0; i < run.t.size(); ++i) {
        const double t = run.t[i];
        const double J = ref.F(L, t);
        const double correction = 0.5 * L * run.dx * ref.Fzz(L, t);
        E = std::max(E, std::abs(run.F[i] - J));
        residual = std::max(residual, std::abs(run.F[i] - J - correction));
        const double G = 0.5 * std::erfc(-(t - mean_ref) / (sd * std::sqrt(2.0)));
        gauss = std::max(gauss, std::abs(run.F[i] - G));
      }
      grids.push_back({n_cells, E, residual, var_excess_rel});
      guard_gauss = gauss;
      std::cout << "[LdfThomas] n=" << n_cells << " x_L=" << x_L << " R=" << R
                << " u=" << run.u << " mean=" << mean << " (ref " << mean_ref << ")"
                << " var-var_num=" << var - var_num << " (ref " << var_ldf << ")"
                << " E=" << E << " residual=" << residual << " max|F-Gauss|=" << gauss << "\n";
    }

    // The leading error is the numerical dispersion: first order ...
    for (size_t g = 1; g < grids.size(); ++g) {
      const double p = ObservedOrder(grids[g - 1].E, grids[g].E);
      EXPECT_NEAR(p, 1.0, 0.15)
        << "max |F - J| converges at order " << p << " between n = " << grids[g - 1].n_cells
        << " and " << grids[g].n_cells << ", expected 1 (upwind)";
    }
    // ... and it is EXACTLY the closed-form correction (L dx / 2) F_zz: what is
    // left converges at second order. The finest pair is where the O(y_feed)
    // sorption effect (~3e-5) starts to show, so it gets a lower bound only.
    const double p_res_coarse = ObservedOrder(grids[0].residual, grids[1].residual);
    const double p_res_fine = ObservedOrder(grids[1].residual, grids[2].residual);
    EXPECT_NEAR(p_res_coarse, 2.0, 0.3)
      << "Residual after the O(dx) correction converges at order " << p_res_coarse
      << " (n = 50 -> 100), expected 2";
    EXPECT_GT(p_res_fine, 1.5)
      << "Residual after the O(dx) correction converges at order " << p_res_fine
      << " (n = 100 -> 200), expected 2";
    EXPECT_LT(grids.back().residual, 0.05 * grids.back().E)
      << "The closed-form O(dx) correction does not account for the grid error";

    // GUARD: the shape is resolved. A symmetric curve with the right mean and
    // variance misses the simulated one by far more than the residual, so the
    // agreement with J is a statement about the asymmetric Thomas shape.
    EXPECT_GT(guard_gauss, 0.03) << "Outlet curve is too close to a Gaussian to test the shape";
    EXPECT_GT(guard_gauss, 20.0 * grids.back().residual);
  }

} // namespace adsim
