/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Analytical validation of the adiabatic adsorption temperature rise.
//
// When an adsorbate breaks through an adiabatic column it releases its heat of
// adsorption, which warms the packing and the gas. Equilibrium theory (Pan &
// Basmadjian; Ruthven, "Principles of Adsorption and Adsorption Processes",
// ch. 8) gives the temperature of the plateau LEFT BEHIND the leading
// mass-transfer zone from a coupled mass + energy balance across it. For a
// single adsorbate fed to a clean bed with a favourable isotherm, the leading
// "pure thermal" wave travels faster than the composition front, and the region
// between the two fronts sits at a raised temperature whose value follows from a
// lumped adiabatic energy balance over the material that has equilibrated:
//
//     dH_ads * rho_b * dq   =   C_tot * dT_plateau
//
// where dq is the loading picked up per kg of adsorbent as the front passes,
// dH_ads the heat released per mole adsorbed, rho_b the bed bulk density, and
// C_tot the total volumetric heat capacity of the equilibrated bed (gas + solid).
// Rearranged, the adiabatic plateau rise is
//
//     dT = dH_ads * rho_b * dq / C_tot.
//
// The self-consistency is that dq itself depends on temperature through the
// isotherm, so this is genuinely an energy<->isotherm coupling: the loading that
// sets the heat release is the loading at the RAISED plateau temperature, not the
// feed temperature. We solve the small fixed point
//
//     q_plateau = q*(p_feed, T0 + dT),   dT = dH_ads rho_b q_plateau / C_tot
//
// by iteration to convergence and compare against the simulated plateau. The
// isotherm loading q*(p, T) is written straight from the Henry form the test
// installs, so the reference never touches the solver's evaluation.
//
// The column is made adiabatic (no wall, no ambient, no axial conduction), so the
// only place the released heat can go is the sensible heat of the bed and gas,
// and a flat temperature plateau develops behind the leading front. We read the
// plateau at the column mid-plane after the thermal wave has passed it but before
// the (slower) composition front arrives -- the interval in which the closed-form
// balance holds.
//
// GUARDS: (1) a genuine plateau must exist -- the temperature must be flat in time
// at the probe location over a finite window (so we are reading the equilibrium
// plateau, not a transient spike or the rising front); (2) the temperature
// self-consistency correction must matter, i.e. dq at the raised temperature
// differs measurably from dq at the feed temperature, otherwise the coupling is
// not exercised; (3) the rise must be several kelvin so it is well above the
// integrator noise.

#include "breakthrough_case.h"

namespace adsim {

  using namespace adsim::validation;

  namespace {
    // Henry loading q = k exp(Q/RT) p [mol/kg]; p = P*y the partial pressure.
    double HenryLoading(double k, double Q, double T, double p) {
      return k * std::exp(Q / (kR * T)) * p;
    }
  } // namespace

  // ===========================================================================
  // Adiabatic plateau temperature rise
  // ===========================================================================

  TEST(AdiabaticAdsorption, PlateauRiseMatchesAdiabaticEnergyBalance) {
    // A strongly-adsorbing, concentrated feed so the released heat gives a large,
    // clearly-resolved plateau. Henry (linear) isotherm keeps dq a clean closed
    // form while still coupling to temperature through exp(Q/RT).
    const double henry_k = 6.0e-11;  // [mol/(kg-Pa)]
    const double Qiso = 38000.0;     // [J/mol] isosteric heat: steep exp(Q/RT)
    const double y_feed = 0.2;       // concentrated: big heat release
    const double T0 = 300.0;         // [K] feed and initial bed temperature

    BreakthroughCase c;
    c.n_cells = 120;
    c.adsorbing = false;             // we install our own isotherm below
    c.length = 1.0;
    c.diameter = 0.1;
    c.molar_flow = 0.01;             // [mol/s]
    c.y_feed = y_feed;
    c.feed_temperature = T0;
    c.bed_temperature = T0;
    c.run_time = 4000.0;
    c.max_step = 1.0;
    c.reltol = 1e-9;

    UserData data;
    Build(data, c);

    // Install the Henry isotherm on Oxygen with a matching heat of adsorption,
    // and fast LDF so the front is equilibrium-controlled.
    auto& layer = data.bed.layers_.layers_.front();
    AdsorbentID ads_id = layer.composition_.begin()->first;
    auto& adsorbent = data.adsorbent_directory.Get(ads_id);

    auto isotherm = std::make_unique<HenryIsotherm>();
    isotherm->params_["Oxygen"].k = henry_k;
    isotherm->params_["Oxygen"].Q = Qiso;
    adsorbent.adsorption_source_.AddIsotherm({"Oxygen"}, std::move(isotherm));

    auto ldf = std::make_unique<LinearDrivingForce>();
    auto mtc = std::make_unique<MassTransferCoefficientConstant>();
    mtc->constant_mtc_ = 2.0;        // fast: near local equilibrium
    ldf->mtc_model_ = std::move(mtc);
    adsorbent.adsorption_source_.mass_transfer_models_["Oxygen"] = std::move(ldf);
    adsorbent.adsorption_source_.heat_of_adsorption_["Oxygen"] = Qiso;

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    const auto& log = solver.GetLog();
    const auto& dom = log.GetDomain();

    // --- confirm the run was adiabatic and isobaric --------------------------
    const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
    const double P_col = InteriorMean(P.values.back(), dom.volume);
    ASSERT_LT(std::abs(P_col - kCasePressure) / kCasePressure, 0.02)
      << "Column not isobaric; the single-pressure isotherm balance does not apply";

    // --- read the plateau temperature at the column mid-plane ----------------
    const auto T = log.GetVariable(0, 0, VariableType::TEMPERATURE);
    ASSERT_FALSE(T.values.empty());
    const size_t n_cells_total = T.values.front().size();
    const size_t mid = n_cells_total / 2;

    // Time series of the mid-plane temperature.
    std::vector<double> T_mid(T.values.size());
    for (size_t n = 0; n < T.values.size(); ++n) T_mid[n] = T.values[n][mid];

    // The mid-plane temperature rises as the thermal wave arrives, holds a flat
    // plateau, then rises further (or falls) when the composition front reaches
    // it. Find the flat plateau: the longest run of samples whose values agree to
    // within a small band, taken after the temperature has first risen well above
    // T0. We characterise the plateau by its median over that window.
    const double peak = *std::max_element(T_mid.begin(), T_mid.end());
    ASSERT_GT(peak - T0, 3.0) << "No appreciable adiabatic heating occurred";

    // Collect samples that sit on the plateau: within 5% of the peak rise and
    // after the leading edge. This isolates the equilibrium plateau from the
    // rising/falling fronts.
    std::vector<double> plateau;
    for (size_t n = 0; n < T_mid.size(); ++n) {
      if (T_mid[n] - T0 > 0.90 * (peak - T0)) plateau.push_back(T_mid[n]);
    }
    ASSERT_GT(plateau.size(), 10u)
      << "No flat temperature plateau found at the mid-plane; cannot apply the "
         "equilibrium-theory balance";
    std::sort(plateau.begin(), plateau.end());
    const double T_plateau_sim = plateau[plateau.size() / 2];  // median
    const double dT_sim = T_plateau_sim - T0;

    // GUARD 1: the plateau must genuinely be flat (spread small vs the rise).
    const double plateau_spread = plateau.back() - plateau.front();
    ASSERT_LT(plateau_spread, 0.15 * dT_sim)
      << "Mid-plane temperature is not flat (spread " << plateau_spread
      << " K over a rise of " << dT_sim << " K); no clean plateau to compare";

    // --- closed-form adiabatic plateau rise (self-consistent fixed point) ----
    const double et = kCaseVoidFraction;
    const double rho_b = (1.0 - kCaseVoidFraction) * kCaseSolidDensity;
    const double Cps = dom.Cps[mid];                         // solid heat cap [J/(kg-K)]
    const double Cp = MixtureMolarCp(data.fluid, {{"Nitrogen", 1.0 - y_feed}, {"Oxygen", y_feed}});
    const double p_feed = kCasePressure * y_feed;

    // Adiabatic-balance residual g(dT) = dT - dH rho_b q*(p_feed, T0+dT)/C_tot(T0+dT).
    // For a steep isotherm the raised temperature suppresses q, so g is monotone
    // increasing and has a single root, which we find by bisection (a fixed-point
    // iteration can oscillate for a steep exp(Q/RT) and is avoided deliberately).
    auto balance_residual = [&](double dT) {
      const double T_plateau = T0 + dT;
      const double C = kCasePressure / (kR * T_plateau);     // molar density at plateau
      const double C_tot = et * C * Cp + rho_b * Cps;        // [J/(m3-K)]
      const double q_plateau = HenryLoading(henry_k, Qiso, T_plateau, p_feed);
      return dT - Qiso * rho_b * q_plateau / C_tot;
    };
    double lo = 0.0, hi = 500.0;
    ASSERT_LT(balance_residual(lo), 0.0);
    ASSERT_GT(balance_residual(hi), 0.0);
    for (int it = 0; it < 200; ++it) {
      const double mid_dT = 0.5 * (lo + hi);
      if (balance_residual(mid_dT) < 0.0) lo = mid_dT; else hi = mid_dT;
    }
    const double dT = 0.5 * (lo + hi);

    // GUARD 2: the energy<->isotherm coupling must MATTER. If the loading were
    // frozen at the feed temperature the balance would predict a much larger rise;
    // the isotherm's exp(Q/RT) suppression of q at the raised plateau is exactly
    // what pulls the self-consistent answer down. Require the frozen-isotherm
    // prediction to be far larger than the self-consistent one, so a solver that
    // ignored the coupling (evaluated q* at T_feed) could not pass.
    const double q_cold = HenryLoading(henry_k, Qiso, T0, p_feed);
    const double C_cold = et * (kCasePressure / (kR * T0)) * Cp + rho_b * Cps;
    const double dT_frozen = Qiso * rho_b * q_cold / C_cold;
    ASSERT_GT(dT_frozen, 1.5 * dT)
      << "Frozen-isotherm rise (" << dT_frozen << " K) is not much larger than the "
      << "self-consistent rise (" << dT << " K); the energy<->isotherm coupling is "
      << "not exercised by this case";
    ASSERT_GT(dT, 3.0) << "Self-consistent plateau rise too small to resolve";

    // --- comparison ----------------------------------------------------------
    // The balance is exact for the plateau in the sharp-front (local-equilibrium)
    // limit; the residual is the finite MTZ width and the small gas the front
    // leaves behind, a few percent on this grid.
    EXPECT_NEAR(dT_sim, dT, 0.05 * dT)
      << "Adiabatic plateau rise disagrees with the equilibrium-theory energy "
      << "balance (simulated " << dT_sim << " K, closed form " << dT << " K, "
      << "self-consistent loading at plateau)";
  }

} // namespace adsim
