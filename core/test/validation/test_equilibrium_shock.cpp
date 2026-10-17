/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Analytical validation against equilibrium (local-equilibrium) theory.
//
// The linear-isotherm breakthrough-time test in test_analytic_transport.cpp
// pins the front speed for a HENRY isotherm, where the wave speed is the same at
// every concentration and the moment is trivially the column capacity. This test
// does the harder, nonlinear case: a favourable LANGMUIR isotherm, where the
// self-sharpening front travels as a genuine shock and its speed is set by the
// SECANT of the isotherm across the jump, not its slope at any single point.
//
// Equilibrium theory (Rhee, Aris & Amundson; LeVan & Carta in Perry's) gives the
// shock velocity from a mass balance across the discontinuity. For a single
// adsorbing species stepped from a clean bed (q = 0, c = 0) to the feed state,
// with negligible dispersion and fast mass transfer (q -> q*), the shock speed is
//
//     u_sh = u / ( 1 + (rho_b / eps_t) * (dq / dc) )
//
// where dq = q*_feed - 0 and dc = c_feed - 0 are the jumps ACROSS the shock, u is
// the INTERSTITIAL velocity, rho_b the bed bulk density and eps_t the total void
// fraction. (The solver logs the superficial velocity; a concentration front
// travels at u_int = u_super / eps_t, so we convert -- see the Ogata-Banks
// dispersion test for the same conversion and its justification.) The
// stoichiometric breakthrough time at the outlet is then
//
//     t_st = L / u_sh = (L / u) * ( 1 + (rho_b / eps_t) * (dq / dc) ).
//
// This is exact for the shock POSITION even when the numerical front has some
// finite width: the stoichiometric time is the first moment of the breakthrough
// curve, and a symmetric numerical smear does not move the centre of mass. We
// therefore verify t_st via the flux-weighted first moment, exactly as the Henry
// test does, but with dq computed from the Langmuir secant.
//
// The Langmuir isotherm as implemented is
//
//     q*(p) = m1 * K(T) * p / ( 1 + K(T) * p ),   K(T) = b1 * exp(q1 / (R T)),
//
// with p = P * y the partial pressure of the adsorbing species.
//
// GUARDS: (1) the isotherm must be genuinely nonlinear across the step, i.e.
// K*p_feed is NOT << 1, so the secant dq/dc differs materially from the Henry
// slope; if it were linear this would collapse into the existing Henry test.
// (2) the adsorbed capacity must dominate the gas capacity, so the front speed is
// actually set by the solid.

#include "breakthrough_case.h"

namespace adsim {

  using namespace adsim::validation;

  namespace {
    // Langmuir loading [mol/kg] at partial pressure p [Pa] and temperature T [K],
    // written straight from the isotherm's algebraic form so the reference does
    // not touch the solver's own evaluation.
    double LangmuirLoading(double m1, double b1, double q1, double T, double p) {
      const double K = b1 * std::exp(q1 / (kR * T));
      return m1 * K * p / (1.0 + K * p);
    }
  } // namespace

  // ===========================================================================
  // Langmuir shock-front speed
  // ===========================================================================

  TEST(EquilibriumShock, LangmuirShockFrontSpeedMatchesSecantBalance) {
    // Langmuir parameters. Chosen so that at the feed partial pressure the bed is
    // well up the isotherm's knee (K * p_feed of order unity), which is where the
    // secant and the local slope diverge and the wave is a true shock.
    const double m1 = 3.0;      // [mol/kg] saturation capacity
    const double b1 = 2.0e-8;   // [1/Pa] pre-exponential affinity
    const double q1 = 20000.0;  // [J/mol] adsorption energy
    const double HOA = 0.0;     // isothermal by construction

    const double y_feed = 0.5;  // NOT dilute: pushes p_feed well up the isotherm
    const double T = 298.0;

    BreakthroughCase c;
    c.n_cells = 200;            // fine grid: sharpens the front, tightens the moment
    c.adsorbing = false;        // we install our own Langmuir isotherm below
    c.length = 1.0;
    c.diameter = 0.1;
    c.molar_flow = 0.02;        // [mol/s]
    c.y_feed = y_feed;
    c.feed_temperature = T;
    c.bed_temperature = T;
    c.run_time = 2500.0;        // ~2.3 x t_st: full breakthrough of the strong adsorbate
    c.max_step = 1.0;
    c.reltol = 1e-9;

    UserData data;
    Build(data, c);

    auto& layer = data.bed.layers_.layers_.front();
    // Find the adsorbent this layer uses.
    AdsorbentID ads_id = layer.composition_.begin()->first;
    auto& adsorbent = data.adsorbent_directory.Get(ads_id);

    // Install a single-site Langmuir isotherm on the adsorbing species (Oxygen),
    // with fast LDF kinetics so the local-equilibrium assumption q -> q* holds.
    auto isotherm = std::make_unique<LangmuirIsotherm>();
    isotherm->params_["Oxygen"].m1 = m1;
    isotherm->params_["Oxygen"].b1 = b1;
    isotherm->params_["Oxygen"].q1 = q1;
    adsorbent.adsorption_source_.AddIsotherm({"Oxygen"}, std::move(isotherm));

    auto ldf = std::make_unique<LinearDrivingForce>();
    auto mtc = std::make_unique<MassTransferCoefficientConstant>();
    mtc->constant_mtc_ = 5.0;   // fast: front is equilibrium-controlled, not kinetic
    ldf->mtc_model_ = std::move(mtc);
    adsorbent.adsorption_source_.mass_transfer_models_["Oxygen"] = std::move(ldf);
    adsorbent.adsorption_source_.heat_of_adsorption_["Oxygen"] = HOA;

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    const auto& log = solver.GetLog();

    // --- confirm the derivation's assumptions -------------------------------
    // Isobaric (so P is a single number in the isotherm) ...
    const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
    const double P_col = P.values.back()[1];
    ASSERT_LT(std::abs(P_col - kCasePressure) / kCasePressure, 0.01)
      << "Column is not isobaric; the single-pressure isotherm secant does not apply";

    // ... and isothermal (HOA = 0), so K(T) is a single value.
    const auto Tvar = log.GetVariable(0, 0, VariableType::TEMPERATURE);
    for (const auto& snap : {Tvar.values.front(), Tvar.values.back()}) {
      for (size_t i = 1; i + 1 < snap.size(); ++i) {
        ASSERT_NEAR(snap[i], T, 1.0) << "Run is not isothermal; K(T) is not constant";
      }
    }

    // Interstitial velocity. The solver logs the SUPERFICIAL (Darcy) velocity;
    // the continuity equation et*dC/dt + d(C u_super)/dz = 0 makes concentration
    // fronts travel at the interstitial velocity u_int = u_super / et. Equilibrium
    // theory is written in that interstitial frame, so we convert. The velocity is
    // taken well upstream of the front at the final time, where the feed-side gas
    // has reached its steady superficial velocity.
    const auto vel = log.GetVariable(0, 0, VariableType::VELOCITY);
    const auto& u_final = vel.values.back();
    const double eps_t = kCaseVoidFraction;                              // ep = 0 => et = eb
    const double u_super = u_final[1];                                   // first interior face (feed side)
    const double u = u_super / eps_t;                                   // interstitial velocity
    ASSERT_GT(u, 0.0);

    // --- analytical shock speed via the secant balance ----------------------
    const double rho_b = (1.0 - kCaseVoidFraction) * kCaseSolidDensity;  // [kg/m3]
    const double p_feed = P_col * y_feed;                                // [Pa]

    const double dq = LangmuirLoading(m1, b1, q1, T, p_feed) - 0.0;      // [mol/kg]
    const double dc = p_feed / (kR * T) - 0.0;                           // [mol/m3]
    const double secant = dq / dc;                                      // [m3/kg]

    const double t_st_expected = (c.length / u) * (1.0 + (rho_b / eps_t) * secant);

    // GUARD 1: the isotherm must be in its nonlinear region. If K*p_feed were
    // small the secant would equal the Henry slope m1*K and this would degenerate
    // into the linear test. Require the secant to be well below the Henry slope.
    const double K_feed = b1 * std::exp(q1 / (kR * T));
    const double henry_slope_per_pressure = m1 * K_feed;                // dq/dp at p->0
    const double secant_per_pressure = dq / p_feed;                     // dq/dp across step
    ASSERT_GT(K_feed * p_feed, 0.3)
      << "K*p_feed = " << K_feed * p_feed << " is too small; isotherm is effectively linear";
    ASSERT_LT(secant_per_pressure, 0.7 * henry_slope_per_pressure)
      << "Secant slope is too close to the Henry slope; the nonlinearity is not exercised "
      << "(secant/pressure " << secant_per_pressure << ", Henry/pressure "
      << henry_slope_per_pressure << ")";

    // GUARD 2: the solid must dominate. Compare the adsorbed retardation term to
    // the gas term (which is 1).
    const double solid_retardation = (rho_b / eps_t) * secant;
    ASSERT_GT(solid_retardation, 5.0)
      << "Solid capacity does not dominate (retardation factor " << solid_retardation
      << "); the front speed is not set by adsorption";
    ASSERT_LT(t_st_expected, 0.6 * c.run_time) << "Run too short to see full breakthrough";

    // --- measured flux-weighted first moment (stoichiometric time) ----------
    auto [times, y_out] = OutletSeries(log, c.feed_from_top);
    const auto n_dot = log.GetFaceVariable(0, 0, VariableType::MOLAR_FLOW_RATE);
    ASSERT_FALSE(times.empty());
    ASSERT_NEAR(y_out.back(), y_feed, 5e-3 * y_feed)
      << "Adsorbate has not fully broken through; the first moment would be truncated";

    // Weight by outlet molar flux of the adsorbing species: adsorption removes
    // moles during the transient so the outlet flow is not the feed flow.
    const double feed_tracer_flow = c.molar_flow * y_feed;
    std::vector<double> unbroken(times.size(), 0.0);
    for (size_t n = 0; n < times.size(); ++n) {
      const double tracer_out = n_dot.values[n].back() * y_out[n];
      unbroken[n] = 1.0 - tracer_out / feed_tracer_flow;
    }
    const double t_st_measured = TimeIntegral(times, unbroken);

    // Tolerance: the first moment is exact for the shock position; the residual
    // is the finite-grid smear of the (self-sharpening) front plus the small
    // amount of gas the front leaves behind it. On a 200-cell grid this is a low
    // single-digit percent.
    EXPECT_NEAR(t_st_measured, t_st_expected, 0.03 * t_st_expected)
      << "Langmuir shock breakthrough time disagrees with the secant balance "
      << "(measured " << t_st_measured << " s, analytic " << t_st_expected
      << " s, retardation " << solid_retardation << ")";
  }

} // namespace adsim
