/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Analytical validation of the linear-driving-force (LDF) mass-transfer kinetics.
//
// The solver writes the adsorbed-concentration equation as
//
//     dq/dt = k ( q* - q ),
//
// so if the equilibrium loading q* is held constant while q starts away from it,
// the exact trajectory is a single exponential relaxation
//
//     q(t) = q* + (q0 - q*) exp(-k t).
//
// The whole point is the RATE: the time constant is exactly 1/k, independent of
// geometry, gas inventory or spatial discretisation.
//
// Realising a constant-q* relaxation is subtler than it looks, because the solver
// ALWAYS initialises the loading at equilibrium with the gas-phase initial
// condition (Equations::SetAdsorbedPhaseInitialConditions). A closed bed that
// starts uniform therefore starts AT equilibrium and never relaxes. To create a
// genuine disequilibrium we instead start the bed CLEAN (no adsorbate in the gas,
// so q0 = 0 by that same equilibrium initialisation) and then step the adsorbate
// in at a HIGH molar flow. We arrange a strong separation of timescales:
//
//     residence time  L/u  <<  LDF time constant  1/k,
//
// so the gas partial pressure fills the whole column to its feed value almost
// instantly and near-uniformly, presenting every cell with an effectively
// constant q* long before the (slow) loading has relaxed. The loading then obeys
// q(t) = q* (1 - exp(-k t)) from t ~ 0, and we recover k from the trajectory.
//
// This keeps the essential feature the task asks for -- a fixed q* and a pure
// exponential in q -- while respecting how the solver seeds the adsorbed phase.
//
// GUARDS: (1) the gas fill must be fast compared with 1/k (asserted from the
// residence time), so q* really is constant over the relaxation; (2) the loading
// must actually climb most of the way to q*, i.e. there is appreciable
// disequilibrium to relax.

#include "validation_helpers.h"

namespace adsim {

  using namespace adsim::validation;

  class LdfKineticsTest : public UserDataTest {
  public:
    LayerID layer_id{};
    AdsorbentID ads_id{};

    static constexpr double kLength = 1.0;
    static constexpr double kDiameter = 0.1;
    static constexpr double kVoidFraction = 0.4;
    static constexpr double kTemperature = 298.0;
    static constexpr double kPressure = 1.0e5;

    void BuildBed(int n_cells) {
      AddComponents({"Nitrogen", "Oxygen"});
      auto ids = AddLayer("Packing");
      layer_id = ids.first;
      ads_id = ids.second;

      auto& layer = GetLayer(layer_id);
      layer.domain_.length = kLength;
      layer.domain_.n_cells = n_cells;
      layer.eb_ = kVoidFraction;
      layer.particle_diameter_ = 0.002;
      std::get<CylinderBed>(data.bed.shape_).diameter = kDiameter;

      auto& adsorbent = GetAdsorbent(ads_id);
      adsorbent.ep_ = 0.0;               // total void == bed void
      adsorbent.solid_density_ = 1000.0;
      adsorbent.heat_capacity_ = 900.0;
      adsorbent.thermal_conductivity_ = 0.0;

      // Isothermal, no wall coupling. Dispersion is turned UP (below) to keep the
      // gas front from smearing the fill, but that only affects how fast/uniform
      // the gas becomes -- it does not enter the loading ODE.
      AddConstantInternalHeatTransferCoefficientModel(0.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);
    }
  };

  // ===========================================================================
  // LDF exponential step response
  // ===========================================================================

  TEST_F(LdfKineticsTest, LoadingRelaxesExponentiallyAtRateK) {
    data.info.title = "Validation: LDF exponential relaxation";

    BuildBed(/*n_cells=*/40);

    // Small Henry capacity so adsorption barely perturbs the gas, and a slow LDF
    // rate so the loading relaxes on a much longer timescale than the gas fill.
    const double henry_k = 1.0e-11;    // [mol/(kg-Pa)] small on purpose
    const double Qiso = 15000.0;       // [J/mol]
    const double k_ldf = 0.02;         // [1/s] LDF rate constant -- the target
    const double y_ads = 0.5;          // feed adsorbate mole fraction

    auto isotherm = std::make_unique<HenryIsotherm>();
    isotherm->params_["Oxygen"].k = henry_k;
    isotherm->params_["Oxygen"].Q = Qiso;
    AddIsotherm({"Oxygen"}, ads_id, std::move(isotherm));
    AddLinearDrivingForceModel(ads_id, "Oxygen", k_ldf);
    SetHOA(ads_id, "Oxygen", 0.0);     // isothermal: q* = k_H exp(Q/RT) p is fixed

    // Strong axial dispersion so the gas front is smeared into a near-uniform,
    // rapid fill rather than a travelling shock; this makes q* spatially uniform
    // during the relaxation. (Dispersion does not appear in the loading ODE.)
    auto disp = std::make_unique<DispersionModelConstant>();
    disp->dispersion_coefficient_ = 0.05;  // [m2/s], very large: fills fast
    GetLayer(layer_id).equations_.dispersion_model_ = std::move(disp);

    // Bed starts CLEAN: no adsorbate in the gas, so the equilibrium
    // initialisation seeds q0 = 0 -- the maximum disequilibrium once the feed
    // arrives.
    AddInitialConditions(kPressure, kTemperature, kTemperature,
                         {{"Nitrogen", 1.0}, {"Oxygen", 0.0}}, 0.0);

    // High molar flow: the residence time is a few seconds, << 1/k = 50 s.
    const double molar_flow = 0.2;     // [mol/s]
    const double run_time = 400.0;     // ~8 time constants at k = 0.02 /s
    AddStep("Fill", run_time);
    Step& step = data.cycle.GetStep("Fill");

    auto inlet = CreateNoValveInletStream(CreateMolarFlowBC(molar_flow),
                                          {{"Nitrogen", 1.0 - y_ads}, {"Oxygen", y_ads}},
                                          kTemperature);
    auto inlet_bc = std::make_unique<FlowBoundaryCondition>();
    inlet_bc->direction_ = BoundaryDirection::INLET;
    inlet_bc->location_ = BoundaryLocation::BOTTOM;
    inlet_bc->streams_.push_back(std::move(inlet));
    step.bottom_boundary_condition_ = std::move(inlet_bc);

    auto outlet_bc = std::make_unique<PressureBoundaryCondition>();
    outlet_bc->direction_ = BoundaryDirection::OUTLET;
    outlet_bc->location_ = BoundaryLocation::TOP;
    outlet_bc->pressure_ = CreateConstantBoundaryValue(kPressure);
    outlet_bc->mole_fraction_["Nitrogen"] = 1.0 - y_ads;
    outlet_bc->mole_fraction_["Oxygen"] = y_ads;
    step.top_boundary_condition_ = std::move(outlet_bc);

    ConfigureCycle(1);
    data.settings.reltol = 1e-10;
    data.settings.max_step = 0.2;

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    const auto& log = solver.GetLog();
    const auto& dom = log.GetDomain();

    // --- analytic equilibrium loading ---------------------------------------
    // q* = k_H exp(Q/RT) * p, p = P * y (the concentration-independent Henry form).
    const double p_ads = kPressure * y_ads;
    const double q_star = henry_k * std::exp(Qiso / (kR * kTemperature)) * p_ads;
    ASSERT_GT(q_star, 0.0);

    // --- confirm the timescale separation (fast, uniform gas fill) ----------
    // Interstitial velocity from the steady field, and the residence time it
    // implies. The gas must fill in a small fraction of 1/k.
    const auto vel = log.GetVariable(0, 0, VariableType::VELOCITY);
    const auto& u_final = vel.values.back();
    double u_sum = 0.0; int u_cnt = 0;
    for (size_t i = 1; i + 1 < u_final.size(); ++i) { u_sum += u_final[i]; ++u_cnt; }
    const double u = u_sum / u_cnt;
    const double t_residence = kLength / u;
    ASSERT_LT(t_residence, 0.1 / k_ldf)
      << "Gas residence time " << t_residence << " s is not << 1/k = " << 1.0 / k_ldf
      << " s; q* is not constant during the relaxation";

    // Isothermal check: k(T) constant.
    const auto Tvar = log.GetVariable(0, 0, VariableType::TEMPERATURE);
    const double Tf = InteriorMean(Tvar.values.back(), dom.volume);
    ASSERT_NEAR(Tf, kTemperature, 1.0) << "Run is not isothermal";

    // --- measured loading trajectory (spatial mean over the bed) ------------
    const auto q = log.GetComponentVariable(0, 0, ComponentVariableType::ADSORBED_CONCENTRATION,
                                            "Oxygen");
    ASSERT_FALSE(q.values.empty());
    std::vector<double> q_mean(q.values.size());
    for (size_t n = 0; n < q.values.size(); ++n) {
      q_mean[n] = InteriorMean(q.values[n], dom.volume);
    }

    // GUARDS on disequilibrium: starts empty, ends at equilibrium.
    ASSERT_LT(q_mean.front(), 0.02 * q_star)
      << "Bed did not start clean (initial mean loading " << q_mean.front()
      << ", q* " << q_star << "); there is no disequilibrium to relax";
    ASSERT_GT(q_mean.back(), 0.95 * q_star)
      << "Loading did not approach equilibrium (final " << q_mean.back()
      << ", q* " << q_star << "); the relaxation is not resolved";

    // --- recover the time constant ------------------------------------------
    // For q(t) = q* (1 - exp(-k t)), ln(1 - q/q*) = -k t is a straight line
    // through the origin with slope -k. Fit over the window 20%..90% of q*: below
    // 20% the gas is still filling (q* not yet fully established), above 90% the
    // saturated tail is dominated by log-of-small-number noise.
    double num = 0.0, den = 0.0;
    int fit_points = 0;
    for (size_t n = 0; n < q.times.size(); ++n) {
      const double frac = q_mean[n] / q_star;
      if (frac < 0.20 || frac > 0.90) continue;
      const double lhs = std::log(1.0 - frac);  // = -k t
      const double t = q.times[n];
      num += t * lhs;
      den += t * t;
      ++fit_points;
    }
    ASSERT_GT(fit_points, 10) << "Too few points in the relaxation window to fit k";
    const double k_measured = -num / den;

    // Tolerance: the loading ODE is exercised with the spatial operator's only
    // job being to keep the gas uniform (which the guard confirms it does to a
    // few percent). The recovered rate matches the input to a few percent; the
    // floor is the residual gas-fill transient bleeding into the fit window,
    // whose speed depends on the gas viscosity / pressure drop (~3.7% with the
    // real component-viscosity model, vs ~1% under the old constant-mu stub).
    EXPECT_NEAR(k_measured, k_ldf, 0.05 * k_ldf)
      << "Recovered LDF rate constant " << k_measured << " /s disagrees with the input "
      << k_ldf << " /s; the loading equation is not dq/dt = k (q* - q)";

    // Independent shape cross-check: at one time constant, q = q*(1 - 1/e).
    const double t_tau = 1.0 / k_ldf;
    const double q_at_tau = interp_linear(q.times, q_mean, t_tau);
    EXPECT_NEAR(q_at_tau / q_star, 1.0 - std::exp(-1.0), 0.02)
      << "Loading at t = 1/k is " << q_at_tau / q_star << " q*, expected "
      << 1.0 - std::exp(-1.0);
  }

  // ===========================================================================
  // Second-order LDF (LDF2): dq/dt = k (q* - q) |q* - q|
  // ===========================================================================
  //
  // With q* held constant the exact solution is
  //
  //     q* - q(t) = (q* - q0) / (1 + k |q* - q0| t),
  //
  // i.e. 1/|q* - q| = 1/|q* - q0| + k t is a straight line of slope k. That is
  // hyperbolic, not exponential: at the "time constant" t = 1/(k |q* - q0|) the
  // remaining disequilibrium is 1/2 (LDF would give 1/e). The same bed set-up as
  // the first-order test realises a constant q*: adsorption steps the adsorbate
  // in onto a clean bed; desorption starts in equilibrium with the adsorbate and
  // purges it out with the carrier, so q* -> 0 and q decays as q0/(1 + k q0 t).

  TEST(LdfRateModel, SecondOrderRateIsSignPreservingAndDiffersFromFirstOrder) {
    LinearDrivingForce ldf;
    LinearDrivingForceSecondOrder ldf2;
    // Adsorption (q* > q): both positive.
    EXPECT_DOUBLE_EQ(ldf.Rate(1.0, 0.5, 2.0), 1.0);
    EXPECT_DOUBLE_EQ(ldf2.Rate(1.0, 0.5, 2.0), 0.5);
    // Desorption (q* < q): both NEGATIVE -- the literal k(q*-q)^2 would be +0.5.
    EXPECT_DOUBLE_EQ(ldf.Rate(0.0, 0.5, 2.0), -1.0);
    EXPECT_DOUBLE_EQ(ldf2.Rate(0.0, 0.5, 2.0), -0.5);
    // Equilibrium: no uptake.
    EXPECT_DOUBLE_EQ(ldf2.Rate(0.7, 0.7, 2.0), 0.0);
    // Quadratic in the driving force.
    EXPECT_DOUBLE_EQ(ldf2.Rate(3.0, 1.0, 0.5), 2.0);
    EXPECT_DOUBLE_EQ(ldf2.Rate(1.0, 3.0, 0.5), -2.0);
  }

  class LdfSecondOrderKineticsTest : public LdfKineticsTest {
  public:
    static constexpr double kHenry = 1.0e-11;   // [mol/(kg-Pa)] small on purpose
    static constexpr double kQiso = 15000.0;    // [J/mol]
    static constexpr double kRateScale = 0.04;  // k |q* - q0| [1/s]

    static double QStar(double y_ads) {
      return kHenry * std::exp(kQiso / (kR * kTemperature)) * kPressure * y_ads;
    }

    struct Trajectory {
      std::vector<double> t;
      std::vector<double> q_mean;
      double t_residence = 0.0;
      double T_final = 0.0;
      double q_star_final = 0.0;  // bed-mean equilibrium loading at the end
    };

    // Fill or purge a bed whose Oxygen kinetics are LDF2 with coefficient k2.
    Trajectory Run(double y_initial, double y_feed, double k2, double run_time) {
      BuildBed(/*n_cells=*/40);

      auto isotherm = std::make_unique<HenryIsotherm>();
      isotherm->params_["Oxygen"].k = kHenry;
      isotherm->params_["Oxygen"].Q = kQiso;
      AddIsotherm({"Oxygen"}, ads_id, std::move(isotherm));
      AddLinearDrivingForceModelSecondOrder(ads_id, "Oxygen", k2);
      SetHOA(ads_id, "Oxygen", 0.0);

      auto disp = std::make_unique<DispersionModelConstant>();
      disp->dispersion_coefficient_ = 0.05;
      GetLayer(layer_id).equations_.dispersion_model_ = std::move(disp);

      AddInitialConditions(kPressure, kTemperature, kTemperature,
                           {{"Nitrogen", 1.0 - y_initial}, {"Oxygen", y_initial}}, 0.0);

      AddStep("Feed", run_time);
      Step& step = data.cycle.GetStep("Feed");
      auto inlet = CreateNoValveInletStream(CreateMolarFlowBC(0.2),
                                            {{"Nitrogen", 1.0 - y_feed}, {"Oxygen", y_feed}},
                                            kTemperature);
      auto inlet_bc = std::make_unique<FlowBoundaryCondition>();
      inlet_bc->direction_ = BoundaryDirection::INLET;
      inlet_bc->location_ = BoundaryLocation::BOTTOM;
      inlet_bc->streams_.push_back(std::move(inlet));
      step.bottom_boundary_condition_ = std::move(inlet_bc);

      auto outlet_bc = std::make_unique<PressureBoundaryCondition>();
      outlet_bc->direction_ = BoundaryDirection::OUTLET;
      outlet_bc->location_ = BoundaryLocation::TOP;
      outlet_bc->pressure_ = CreateConstantBoundaryValue(kPressure);
      outlet_bc->mole_fraction_["Nitrogen"] = 1.0 - y_feed;
      outlet_bc->mole_fraction_["Oxygen"] = y_feed;
      step.top_boundary_condition_ = std::move(outlet_bc);

      ConfigureCycle(1);
      data.settings.reltol = 1e-10;
      data.settings.max_step = 0.2;

      CoupledSolver solver{data};
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

      const auto& log = solver.GetLog();
      const auto& dom = log.GetDomain();
      Trajectory out;

      const auto vel = log.GetVariable(0, 0, VariableType::VELOCITY);
      const auto& u_final = vel.values.back();
      double u_sum = 0.0; int u_cnt = 0;
      for (size_t i = 1; i + 1 < u_final.size(); ++i) { u_sum += u_final[i]; ++u_cnt; }
      out.t_residence = kLength / (u_sum / u_cnt);

      const auto Tvar = log.GetVariable(0, 0, VariableType::TEMPERATURE);
      out.T_final = InteriorMean(Tvar.values.back(), dom.volume);

      const auto q = log.GetComponentVariable(0, 0, ComponentVariableType::ADSORBED_CONCENTRATION,
                                              "Oxygen");
      out.t = q.times;
      out.q_mean.resize(q.values.size());
      for (size_t n = 0; n < q.values.size(); ++n) out.q_mean[n] = InteriorMean(q.values[n], dom.volume);

      const auto qs = log.GetComponentVariable(0, 0, ComponentVariableType::EQUILIBRIUM_ADSORBED_CONCENTRATION,
                                               "Oxygen");
      if (!qs.values.empty()) out.q_star_final = InteriorMean(qs.values.back(), dom.volume);
      return out;
    }

    // Least-squares slope of 1/|q* - q| against t over the window where the
    // remaining disequilibrium fraction |q* - q|/|q* - q0| lies in [lo, hi].
    static double FitSecondOrderK(const Trajectory& tr, double q_star, double q0,
                                  double lo, double hi, int& n_points) {
      const double d0 = std::abs(q_star - q0);
      double St = 0, Sy = 0, Stt = 0, Sty = 0; n_points = 0;
      for (size_t n = 0; n < tr.t.size(); ++n) {
        const double d = std::abs(q_star - tr.q_mean[n]);
        const double frac = d / d0;
        if (frac < lo || frac > hi) continue;
        const double y = 1.0 / d;
        St += tr.t[n]; Sy += y; Stt += tr.t[n] * tr.t[n]; Sty += tr.t[n] * y;
        ++n_points;
      }
      return (n_points * Sty - St * Sy) / (n_points * Stt - St * St);
    }
  };

  TEST_F(LdfSecondOrderKineticsTest, AdsorptionFollowsHyperbolicSecondOrderLaw) {
    data.info.title = "Validation: LDF2 adsorption";
    const double y_ads = 0.5;
    const double q0 = 0.0;                          // clean bed
    const double k2 = kRateScale / QStar(y_ads);    // [1/(s mol/kg)]
    const auto tr = Run(/*y_initial=*/0.0, /*y_feed=*/y_ads, k2, /*run_time=*/400.0);
    ASSERT_FALSE(tr.q_mean.empty());

    // 1/(q* - q) is very sensitive to q* near the tail, so use the solver's own
    // bed-mean equilibrium loading rather than the outlet-pressure estimate: the
    // pressure drop raises the mean partial pressure (and q*) by a few percent.
    const double q_star = tr.q_star_final;
    ASSERT_NEAR(q_star, QStar(y_ads), 0.05 * QStar(y_ads)) << "Unexpected equilibrium loading";

    ASSERT_LT(tr.t_residence, 0.1 / kRateScale) << "q* is not constant during the relaxation";
    ASSERT_NEAR(tr.T_final, kTemperature, 1.0) << "Run is not isothermal";
    ASSERT_LT(tr.q_mean.front(), 0.02 * q_star) << "Bed did not start clean";
    // Hyperbolic tail: at t = 400 s, 1 - 1/(1 + 16) = 0.94 of q*.
    ASSERT_GT(tr.q_mean.back(), 0.90 * q_star) << "Loading did not approach equilibrium";

    int n = 0;
    const double k_measured = FitSecondOrderK(tr, q_star, q0, 0.10, 0.80, n);
    ASSERT_GT(n, 10) << "Too few points in the relaxation window to fit k";
    EXPECT_NEAR(k_measured, k2, 0.05 * k2)
      << "Recovered LDF2 coefficient " << k_measured << " disagrees with input " << k2
      << "; the loading equation is not dq/dt = k (q* - q)|q* - q|";

    // At t = 1/(k q*) the analytic LDF2 loading is q*/2; first-order LDF with
    // the same k q* would give q*(1 - 1/e) = 0.632 q*.
    const double t_half = 1.0 / (k2 * q_star);
    const double frac = interp_linear(tr.t, tr.q_mean, t_half) / q_star;
    EXPECT_NEAR(frac, 0.5, 0.02) << "LDF2 loading at t = 1/(k q*) should be q*/2";
    EXPECT_GT(std::abs(frac - (1.0 - std::exp(-1.0))), 0.1)
      << "LDF2 trajectory is indistinguishable from first-order LDF";
  }

  TEST_F(LdfSecondOrderKineticsTest, DesorptionFollowsHyperbolicSecondOrderLaw) {
    data.info.title = "Validation: LDF2 desorption";
    const double y_init = 0.5;
    const double q0 = QStar(y_init);         // starts in equilibrium with the gas
    const double q_star = 0.0;               // purged with pure carrier
    const double k2 = kRateScale / q0;       // [1/(s mol/kg)]
    const auto tr = Run(/*y_initial=*/y_init, /*y_feed=*/0.0, k2, /*run_time=*/400.0);
    ASSERT_FALSE(tr.q_mean.empty());

    ASSERT_LT(tr.t_residence, 0.1 / kRateScale) << "q* is not constant during the relaxation";
    ASSERT_NEAR(tr.T_final, kTemperature, 1.0) << "Run is not isothermal";
    ASSERT_NEAR(tr.q_mean.front(), q0, 0.02 * q0) << "Bed did not start at equilibrium";
    // The literal k (q* - q)^2 could never desorb: the loading must fall.
    ASSERT_LT(tr.q_mean.back(), 0.10 * q0) << "Loading did not desorb";

    int n = 0;
    const double k_measured = FitSecondOrderK(tr, q_star, q0, 0.10, 0.80, n);
    ASSERT_GT(n, 10) << "Too few points in the relaxation window to fit k";
    EXPECT_NEAR(k_measured, k2, 0.05 * k2)
      << "Recovered LDF2 coefficient " << k_measured << " disagrees with input " << k2
      << "; the loading equation is not dq/dt = k (q* - q)|q* - q|";

    // At t = 1/(k q0): analytic LDF2 gives q0/2; first-order LDF gives q0/e.
    const double t_half = 1.0 / (k2 * q0);
    const double frac = interp_linear(tr.t, tr.q_mean, t_half) / q0;
    EXPECT_NEAR(frac, 0.5, 0.02) << "LDF2 loading at t = 1/(k q0) should be q0/2";
    EXPECT_GT(std::abs(frac - std::exp(-1.0)), 0.1)
      << "LDF2 trajectory is indistinguishable from first-order LDF";
  }

} // namespace adsim
