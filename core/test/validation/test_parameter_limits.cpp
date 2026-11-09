/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Parameter limits and degenerate mixtures.
//
// As a parameter goes to zero or infinity the model must reach a known
// limit: an equilibrium front, an inert tracer, an isothermal bed, a Henry
// isotherm, no separation. These tests do not stop at "the limit is close";
// they measure the RATE at which it is approached as the parameter is varied
// over a geometric sequence, and compare it with the rate theory gives (error
// ~ 1/k, ~ k, ~ dH, ~ 1/C, ~ b p, ~ swing^2, ~ y). A defect that leaves the
// limit approximately right but approaches it the wrong way (a wrong power, a
// missing factor, a term that does not vanish) shows up as the wrong order.
//
// The observed order between two members of a sequence that differ by a
// factor r is log(e1/e2)/log(r). Each tolerance on an order is the size of
// the next term of the expansion at the members used, stated with the test.
//
// Degenerate mixtures: two components that are physically the same are not
// separated (they keep their feed ratio in every cell, at all times, through
// a pressure-swing cycle); a component with zero affinity is an inert.

#include "breakthrough_case.h"
#include "psa_cycle_case.h"
#include "adsim/objects/deck_edit.h"

#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    struct Curve {
      std::vector<double> t, y, tracer_flux;  // outlet time, mole fraction, tracer molar flow
      double max_dT = 0.0;     // max |T - T0| over interior cells and logged instants [K]
      double max_dTw = 0.0;    // the same for the wall
      std::vector<double> final_y;  // end-of-run tracer profile (interior cells)
    };

    Adsorbent& PackingOf(UserData& data) {
      return data.adsorbent_directory.Get(data.bed.layers_.layers_.front().composition_.begin()->first);
    }

    void SetLdfCoefficient(UserData& data, double k) {
      auto& model = PackingOf(data).adsorption_source_.mass_transfer_models_.at("Oxygen");
      auto* ldf = dynamic_cast<LinearDrivingForce*>(model.get());
      ASSERT_NE(ldf, nullptr);
      auto* mtc = dynamic_cast<MassTransferCoefficientConstant*>(ldf->mtc_model_.get());
      ASSERT_NE(mtc, nullptr);
      mtc->constant_mtc_ = k;
    }

    Curve RunCase(UserData& data, const BreakthroughCase& c) {
      Curve r;
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      CoupledSolver solver{data};
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
      const auto& log = solver.GetLog();
      std::tie(r.t, r.y) = OutletSeries(log, c.feed_from_top);
      const auto n = log.GetFaceVariable(0, 0, VariableType::MOLAR_FLOW_RATE);
      for (size_t k = 0; k < r.t.size(); ++k) r.tracer_flux.push_back(n.values[k].back() * r.y[k]);
      for (const auto& s : log.GetVariable(0, 0, VariableType::TEMPERATURE).values)
        for (size_t i = 1; i + 1 < s.size(); ++i) r.max_dT = std::max(r.max_dT, std::abs(s[i] - c.bed_temperature));
      for (const auto& s : log.GetVariable(0, 0, VariableType::WALL_TEMPERATURE).values)
        for (size_t i = 1; i + 1 < s.size(); ++i) r.max_dTw = std::max(r.max_dTw, std::abs(s[i] - c.bed_temperature));
      const auto y = log.GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION, "Oxygen").values.back();
      r.final_y.assign(y.begin() + 1, y.end() - 1);
      return r;
    }

    // max_t |y_a(t) - y_b(t)| on a fine common grid.
    double MaxDifference(const Curve& a, const Curve& b, double t_end) {
      double worst = 0.0;
      for (int k = 0; k <= 1000; ++k) {
        const double t = t_end * k / 1000.0;
        worst = std::max(worst, std::abs(interp_linear(a.t, a.y, t) - interp_linear(b.t, b.y, t)));
      }
      return worst;
    }

    double Order(double e_coarse, double e_fine, double ratio) {
      return std::log(e_coarse / e_fine) / std::log(ratio);
    }

    // The breakthrough case's column: superficial velocity u_s [m/s] and the
    // dimensionless Henry capacity H = rho_b K R T exp(Q/RT) (adsorbed moles
    // per bed volume over gas concentration), from the case constants.
    double CaseSuperficialVelocity(const BreakthroughCase& c) {
      return c.molar_flow * kR * c.bed_temperature / (0.25 * M_PI * c.diameter * c.diameter * kCasePressure);
    }
    double CaseHenryCapacity(const BreakthroughCase& c) {
      const double rho_b = (1.0 - kCaseVoidFraction) * kCaseSolidDensity;
      return rho_b * kCaseHenryK * std::exp(kCaseHenryQ / (kR * c.bed_temperature)) * kR * c.bed_temperature;
    }

  } // namespace

  // ===========================================================================
  // 1. k_LDF -> infinity: the equilibrium-theory front, at rate 1/k
  // ===========================================================================
  //
  // Linear isotherm, no dispersion. Equilibrium theory (k = infinity) gives a
  // step at the stoichiometric time t_st = (L/u_s)(et + H); the LDF
  // resistance spreads it. In Laplace space the outlet transfer function of
  // the semi-discrete column (N upwind cells, each a stirred tank of
  // residence dx/u_s) is G(s) = [1 + (dx/u_s)(et s + H s k/(s + k))]^-N, and
  // its cumulants are exact:
  //
  //     mean      = (L/u_s)(et + H)                       (independent of k)
  //     variance  = 2 (L/u_s) H / k  +  (L dx/u_s^2)(et + H)^2
  //
  // The LDF part and the upwind part ADD, so the grid drops out of
  // differences: var(k) - var(2k) = (L/u_s) H / k exactly. The test checks
  // that the mean does not move with k, that those differences halve as k
  // doubles (order 1 in 1/k), and that their size is (L/u_s) H / k.
  //
  // The derivation is for a trace adsorbate. At a finite feed fraction y the
  // uptake also changes the gas velocity, which sharpens the front by a
  // relative O(y) amount (measured: -2.2 y to -2.6 y, the same on 80 and 160
  // cells).
  // So the size is checked in the trace limit: the relative error must halve
  // with y (order 1), and its Richardson extrapolation to y = 0 must remove
  // at least 90 % of it.

  TEST(ParameterLimits, LdfFrontApproachesEquilibriumTheoryAsOneOverK) {
    const std::vector<double> ks{1.0, 2.0, 4.0, 8.0};
    const std::vector<double> ys{0.005, 0.0025};
    std::map<double, std::vector<double>> variance, mean;  // keyed by y_feed

    BreakthroughCase c;
    c.adsorbing = true;
    c.run_time = 500.0;
    c.max_step = 0.5;
    for (double y : ys) {
      for (double k : ks) {
        c.y_feed = y;
        UserData data;
        Build(data, c);
        // The cumulants above are those of N upwind cells (the default is WENO).
        data.settings.advection_scheme_ = AdvectionScheme::Upwind;
        SetLdfCoefficient(data, k);
        const Curve r = RunCase(data, c);
        ASSERT_NEAR(r.y.back(), y, 1e-4 * y) << "Breakthrough incomplete";
        std::vector<double> unbroken(r.t.size()), t_unbroken(r.t.size());
        for (size_t n = 0; n < r.t.size(); ++n) {
          unbroken[n] = 1.0 - r.tracer_flux[n] / (c.molar_flow * y);
          t_unbroken[n] = r.t[n] * unbroken[n];
        }
        const double m1 = TimeIntegral(r.t, unbroken);
        mean[y].push_back(m1);
        variance[y].push_back(2.0 * TimeIntegral(r.t, t_unbroken) - m1 * m1);
      }
    }

    const double u_s = CaseSuperficialVelocity(c);
    const double H = CaseHenryCapacity(c);
    const double t_st = (c.length / u_s) * (kCaseVoidFraction + H);
    std::map<double, std::vector<double>> relative_error;
    for (double y : ys) {
      // The mean is the column capacity over the feed rate, whatever k is.
      for (double m1 : mean[y]) {
        EXPECT_NEAR(m1, mean[y].front(), 1e-3 * t_st) << "The front moved with k (y = " << y << ")";
        EXPECT_NEAR(m1, t_st, 2e-3 * t_st) << "Mean is not the stoichiometric time";
      }
      for (size_t i = 0; i + 1 < ks.size(); ++i) {
        const double dvar = variance[y][i] - variance[y][i + 1];
        const double expected = (c.length / u_s) * H / ks[i];
        relative_error[y].push_back(dvar / expected - 1.0);
        std::cout << "[LDF k->inf] y " << y << ", k " << ks[i] << " -> " << ks[i + 1] << ": var drop "
                  << dvar << " s^2, trace theory " << expected << " s^2\n";
        if (i + 2 < ks.size()) {
          const double next = variance[y][i + 1] - variance[y][i + 2];
          // Exactly 1 for the linear problem; the O(y) velocity effect moves
          // it by < y.
          EXPECT_NEAR(Order(dvar, next, 2.0), 1.0, 0.01) << "Spreading does not decay as 1/k";
        }
      }
    }
    // Guard: at k = 1 the LDF spreading is comparable to the numerical one.
    ASSERT_GT(variance[ys[0]][0] - variance[ys[0]][1], 0.25 * variance[ys[0]][1]);

    for (size_t i = 0; i + 1 < ks.size(); ++i) {
      const double r1 = relative_error[ys[0]][i], r2 = relative_error[ys[1]][i];
      const double r0 = 2.0 * r2 - r1;  // Richardson in y (order 1)
      std::cout << "[LDF k->inf] relative error of the var drop at y = " << ys[0] << ", " << ys[1]
                << ": " << r1 << ", " << r2 << "; extrapolated to y = 0: " << r0 << "\n";
      EXPECT_NEAR(Order(std::abs(r1), std::abs(r2), 2.0), 1.0, 0.25) << "Non-trace error is not O(y)";
      EXPECT_LT(std::abs(r0), 0.1 * std::abs(r1)) << "The trace limit does not reach (L/u_s) H / k";
    }
  }

  // ===========================================================================
  // 2. k_LDF -> 0: an inert tracer, at rate k
  // ===========================================================================
  //
  // With k -> 0 nothing is taken up while the gas passes and the outlet curve
  // tends to the inert tracer's. To first order in k the bed starts clean
  // (q = 0) and each element of gas loses tracer at the rate rho_b k q* =
  // k H c per unit bed volume during its residence et L/u_s, so behind the
  // front the outlet is short of the feed by the fraction
  //
  //     deficit = 1 - exp(-k H L / u_s)  ~  k H L/u_s  (1 - k H L/(2 u_s) - ...),
  //
  // (H L/u_s = 121 s here). The test checks the outlet curve's distance from
  // the inert curve is O(k) (order 1 within 0.05, the size of the next term at
  // the largest k), and that the deficit per unit k, Richardson-extrapolated
  // to k = 0, is H L / u_s, with the extrapolation removing at least 90 % of
  // the raw error.

  TEST(ParameterLimits, VanishingLdfCoefficientGivesAnInertTracer) {
    BreakthroughCase c;
    c.y_feed = 0.01;
    c.run_time = 150.0;
    c.adsorbing = false;
    UserData inert;
    Build(inert, c);
    const Curve r_inert = RunCase(inert, c);

    const std::vector<double> ks{1e-4, 2e-4, 4e-4};
    const double t_eval = 100.0;  // behind the inert front (residence time ~25 s)
    std::vector<double> distance, deficit_per_k;
    c.adsorbing = true;
    for (double k : ks) {
      UserData data;
      Build(data, c);
      SetLdfCoefficient(data, k);
      const Curve r = RunCase(data, c);
      distance.push_back(MaxDifference(r, r_inert, c.run_time));
      deficit_per_k.push_back((1.0 - interp_linear(r.t, r.tracer_flux, t_eval) / (c.molar_flow * c.y_feed)) / k);
    }
    const double expected = CaseHenryCapacity(c) * c.length / CaseSuperficialVelocity(c);
    for (size_t i = 0; i < ks.size(); ++i) {
      std::cout << "[LDF k->0] k " << ks[i] << ": max |y - y_inert| " << distance[i] << ", deficit/k "
                << deficit_per_k[i] << " s (H L/u_s = " << expected << " s)\n";
      if (i + 1 < ks.size())
        EXPECT_NEAR(Order(distance[i + 1], distance[i], 2.0), 1.0, 0.05) << "Not an O(k) approach";
    }
    ASSERT_GT(distance[0], 10.0 * 1e-8) << "Effect below the integrator tolerance";
    const double extrapolated = 2.0 * deficit_per_k[0] - deficit_per_k[1];
    std::cout << "[LDF k->0] extrapolated deficit/k " << extrapolated << " s\n";
    EXPECT_LT(std::abs(extrapolated - expected), 0.1 * std::abs(deficit_per_k[0] - expected));
  }

  // ===========================================================================
  // 3. dH -> 0: an adiabatic run approaches the isothermal one, linearly
  // ===========================================================================
  //
  // The breakthrough case is adiabatic (no wall exchange), so the heat of
  // adsorption heats the bed and lowers its capacity. The source is
  // rho_b dH dq/dt, linear in dH, so the temperature rise and the change in
  // the outlet curve are O(dH): halving dH must halve both (order 1; the next
  // term is O(dH^2), measured order 1.000 +/- 0.001 for dH <= 10 kJ/mol). At
  // dH = 0 what is left is the compression heating of the gas, and the run
  // must then match the explicit isothermal mode far more closely than the
  // smallest dH effect probed.

  TEST(ParameterLimits, AdiabaticRunApproachesIsothermalAsHeatOfAdsorptionVanishes) {
    BreakthroughCase c;
    c.adsorbing = true;
    c.run_time = 400.0;
    c.reltol = 1e-8;
    UserData iso;
    Build(iso, c);
    iso.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
    const Curve r_iso = RunCase(iso, c);

    const double dH0 = 40000.0;  // [J/mol]
    const std::vector<double> lambdas{0.25, 0.125, 0.0625, 0.0};
    std::vector<Curve> runs;
    for (double lambda : lambdas) {
      UserData data;
      Build(data, c);
      PackingOf(data).adsorption_source_.heat_of_adsorption_["Oxygen"] = lambda * dH0;
      runs.push_back(RunCase(data, c));
    }
    const Curve& r0 = runs.back();
    std::vector<double> shift;
    for (size_t i = 0; i + 1 < runs.size(); ++i) {
      shift.push_back(MaxDifference(runs[i], r0, c.run_time));
      std::cout << "[dH->0] dH " << lambdas[i] * dH0 << " J/mol: max |dT| " << runs[i].max_dT
                << " K, max |y - y(dH=0)| / y_feed " << shift.back() / c.y_feed << "\n";
    }
    const double residual = MaxDifference(r0, r_iso, c.run_time);
    std::cout << "[dH->0] dH 0: max |y - y_isothermal| / y_feed " << residual / c.y_feed
              << ", max |dT| " << r0.max_dT << " K\n";

    ASSERT_GT(runs.front().max_dT, 0.05) << "The heat of adsorption barely heats the bed";
    for (size_t i = 0; i + 1 < shift.size(); ++i) {
      EXPECT_NEAR(Order(shift[i], shift[i + 1], 2.0), 1.0, 0.02) << "Curve shift not O(dH)";
      EXPECT_NEAR(Order(runs[i].max_dT, runs[i + 1].max_dT, 2.0), 1.0, 0.02) << "Temperature rise not O(dH)";
    }
    EXPECT_LT(residual, 0.05 * shift.back()) << "dH = 0 does not reach the isothermal run";
  }

  // ===========================================================================
  // 4. Wall heat capacity -> infinity / bed-wall coefficient -> infinity
  // ===========================================================================
  //
  // The wall exchanges heat with the bed (h_i = 20 W/m2K, the bed heated by
  // dH = 40 kJ/mol) and nothing else (no ambient loss, zero-flux wall ends).
  // Its temperature then changes by the heat it receives over its heat
  // capacity, so as rho_w Cp_w -> infinity the wall stays at T0 with an
  // excursion ~ 1/(rho_w Cp_w): a tenfold heat capacity must give a tenfold
  // smaller excursion (order 1; the bed temperature it responds to changes by
  // O(1/C) too, so the next term is O(1/C^2): order within 0.02 from 10x).
  //
  // With that isothermal wall, the bed itself approaches the wall temperature
  // as h_i -> infinity: once the bed-wall time constant rho_b Cp/(h_i a) is
  // well below the 2 s of the LDF uptake that releases the heat, the bed's
  // excursion is quasi-steady, source/(h_i a) ~ 1/h_i (order within 0.05 from
  // h_i = 1e5 W/m2K, where that time constant is 0.14 s).

  TEST(ParameterLimits, WallAndBedBecomeIsothermalAsHeatCapacityAndExchangeGrow) {
    BreakthroughCase c;
    c.adsorbing = true;
    c.run_time = 400.0;
    c.reltol = 1e-8;
    auto make = [&](UserData& data, double wall_capacity_scale, double h_i) {
      Build(data, c);
      PackingOf(data).adsorption_source_.heat_of_adsorption_["Oxygen"] = 40000.0;
      auto ihtc = std::make_unique<ConstantInternalHeatTransferModel>();
      ihtc->heat_transfer_coefficient_ = h_i;
      data.bed.layers_.layers_.front().equations_.internal_heat_transfer_model_ = std::move(ihtc);
      data.bed.equations_.wall_temperature_boundary_condition_ =
        CreateWallTemperatureBoundaryCondition(WallTemperatureBoundaryConditionType::ZERO_FLUX);
      data.bed.ambient_temperature_ = c.bed_temperature;
      data.bed.wall_.heat_capacity_ = 537.0 * wall_capacity_scale;
    };

    std::vector<double> wall_excursion;
    double bed_excursion = 0.0;
    for (double scale : {1.0, 10.0, 100.0, 1000.0}) {
      UserData data;
      make(data, scale, 20.0);
      const Curve r = RunCase(data, c);
      wall_excursion.push_back(r.max_dTw);
      bed_excursion = std::max(bed_excursion, r.max_dT);
      std::cout << "[wall C->inf] wall heat capacity x" << scale << ": max |Tw - T0| " << r.max_dTw
                << " K, max |T - T0| " << r.max_dT << " K\n";
    }
    ASSERT_GT(wall_excursion.front(), 0.01) << "The wall barely heats";
    ASSERT_GT(bed_excursion, 0.2);
    for (size_t i = 1; i + 1 < wall_excursion.size(); ++i)
      EXPECT_NEAR(Order(wall_excursion[i], wall_excursion[i + 1], 10.0), 1.0, 0.02)
        << "Wall excursion not O(1/C)";

    std::vector<double> bed;
    for (double h_i : {1e5, 1e6, 1e7}) {
      UserData data;
      make(data, 1e4, h_i);
      const Curve r = RunCase(data, c);
      bed.push_back(r.max_dT);
      std::cout << "[h_i->inf] h_i " << h_i << ": max |T - T0| " << r.max_dT << " K, max |Tw - T0| "
                << r.max_dTw << " K\n";
    }
    for (size_t i = 0; i + 1 < bed.size(); ++i)
      EXPECT_NEAR(Order(bed[i], bed[i + 1], 10.0), 1.0, 0.05) << "Bed excursion not O(1/h_i)";
  }

  // ===========================================================================
  // 5. b -> 0: Langmuir reduces to Henry, at rate b p
  // ===========================================================================
  //
  // q = m b' p / (1 + b' p) with m b' held at the Henry constant: as b' p ->
  // 0 the isotherm is H p (1 - b' p + ...), and the simulated outlet curve
  // must approach the Henry run's with an error O(b' p_feed) (order 1; next
  // term O((b'p)^2), so within 0.05 for b'p_feed <= 0.1). This is the
  // simulation-level statement; the isotherm function's own limit is in
  // test_thermodynamic_consistency.cpp.

  TEST(ParameterLimits, LangmuirRunReducesToHenryRunAsAffinityVanishes) {
    BreakthroughCase c;
    c.adsorbing = true;
    c.run_time = 400.0;
    c.reltol = 1e-8;
    UserData henry;
    Build(henry, c);
    const Curve r_henry = RunCase(henry, c);

    const double e = std::exp(kCaseHenryQ / (kR * c.bed_temperature));
    const double p_feed = c.y_feed * kCasePressure;
    const std::vector<double> bps{0.1, 0.05, 0.025};
    std::vector<double> distance;
    for (double bp : bps) {
      UserData data;
      Build(data, c);
      auto isotherm = std::make_unique<LangmuirIsotherm>();
      const double b1 = bp / (e * p_feed);
      isotherm->params_["Oxygen"].b1 = b1;
      isotherm->params_["Oxygen"].m1 = kCaseHenryK / b1;  // same Henry constant
      isotherm->params_["Oxygen"].q1 = kCaseHenryQ;
      PackingOf(data).adsorption_source_.AddIsotherm({"Oxygen"}, std::move(isotherm));
      distance.push_back(MaxDifference(RunCase(data, c), r_henry, c.run_time));
      std::cout << "[b->0] b' p_feed " << bp << ": max |y_L - y_H| / y_feed " << distance.back() / c.y_feed << "\n";
    }
    ASSERT_GT(distance.front(), 0.05 * c.y_feed) << "The Langmuir curvature barely matters";
    for (size_t i = 0; i + 1 < distance.size(); ++i)
      EXPECT_NEAR(Order(distance[i], distance[i + 1], 2.0), 1.0, 0.05) << "Not an O(b p) approach";
  }

  // ===========================================================================
  // 6. Pressure ratio -> 1: no separation
  // ===========================================================================
  //
  // A pressure-swing cycle separates only through the swing: what blowdown
  // removes is the gas the pressure change desorbs. With P_high = (1 + s)
  // P_low the moles swung are O(s) and their composition differs from the
  // feed by O(s), so at CSS the heavy (blowdown) product's deviation from the
  // feed composition is O(s) and the light product's enrichment, which is
  // their product spread over the feed step, is O(s^2). Both orders are
  // checked (within 0.15: the feed step's own pressure drop adds a swing of
  // ~0.002 to every member, which moves the orders by ~2 x 0.002/0.025), and
  // with no swing the product must equal the feed to within the enrichment
  // that residual pressure drop explains (< 1 % of the s = 0.025 case).
  //
  // Isothermal, pressure-ramp pressurisation to P_high, 0.02 mol/s feed for
  // 120 s with the product at P_high, blowdown at the feed end to P_low.

  namespace {
    void BuildSwing(UserData& data, double P_high, double P_low) {
      BuildPsaCycle(data, 30);
      data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
      auto bc = std::make_unique<PressureBoundaryCondition>();
      bc->direction_ = BoundaryDirection::INLET;
      bc->location_ = BoundaryLocation::BOTTOM;
      auto ramp = std::make_unique<LinearBoundaryValue>();  // from the bed pressure, over the step
      ramp->end_value_ = P_high;
      bc->pressure_ = std::move(ramp);
      bc->temperature_ = 298.0;
      bc->mole_fraction_["Nitrogen"] = 0.79;
      bc->mole_fraction_["Oxygen"] = 0.21;
      data.cycle.GetStep("Pressurise").bottom_boundary_condition_ = std::move(bc);
      auto* top = dynamic_cast<PressureBoundaryCondition*>(data.cycle.GetStep("Feed").top_boundary_condition_.get());
      static_cast<ConstantBoundaryValue*>(top->pressure_.get())->constant_value_ = P_high;
      auto* bottom = dynamic_cast<PressureBoundaryCondition*>(data.cycle.GetStep("Blowdown").bottom_boundary_condition_.get());
      static_cast<ConstantBoundaryValue*>(bottom->pressure_.get())->constant_value_ = P_low;
      data.ic.pressure = P_low;
    }

    // Cycle-mean O2 mole fraction of a step's pressure-boundary stream.
    double MeanOxygenFraction(const StreamData& sd) {
      double n_o2 = 0.0, n = 0.0;
      const auto& o2 = sd.component_molar_flow_.at("Oxygen");
      for (size_t k = 1; k < sd.time_series_.size(); ++k) {
        const double dt = sd.time_series_[k] - sd.time_series_[k - 1];
        n_o2 += 0.5 * (o2[k] + o2[k - 1]) * dt;
        n += 0.5 * (sd.molar_flow_[k] + sd.molar_flow_[k - 1]) * dt;
      }
      return n != 0.0 ? n_o2 / n : 0.0;
    }
  } // namespace

  TEST(ParameterLimits, PressureSwingSeparationVanishesAsThePressureRatioGoesToOne) {
    const std::vector<double> swings{0.1, 0.05, 0.025, 0.0};
    std::vector<double> enrichment, heavy_deviation;
    for (double s : swings) {
      UserData data;
      BuildSwing(data, 1.0e5 * (1.0 + s), 1.0e5);
      data.cycle.number_of_cycles_ = 12;
      data.settings.reltol = 1e-8;
      data.settings.max_step = 2.0;
      std::string error;
      ASSERT_TRUE(data.PerformChecks(error)) << error;
      CoupledSolver solver{data};
      ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);
      auto* light = dynamic_cast<PressureBoundaryCondition*>(data.cycle.GetStep("Feed").top_boundary_condition_.get());
      auto* heavy = dynamic_cast<PressureBoundaryCondition*>(data.cycle.GetStep("Blowdown").bottom_boundary_condition_.get());
      enrichment.push_back(MeanOxygenFraction(light->stream_data_) - 0.21);
      heavy_deviation.push_back(MeanOxygenFraction(heavy->stream_data_) - 0.21);
      std::cout << "[ratio->1] P_high/P_low " << 1.0 + s << ": light-product O2 enrichment "
                << enrichment.back() << ", heavy-product O2 deviation " << heavy_deviation.back() << "\n";
    }
    ASSERT_GT(enrichment.front(), 1e-4) << "The 10 % swing barely separates";
    ASSERT_LT(heavy_deviation.front(), -1e-2);
    for (size_t i = 0; i + 2 < swings.size(); ++i) {
      EXPECT_NEAR(Order(enrichment[i], enrichment[i + 1], 2.0), 2.0, 0.15) << "Enrichment not O(s^2)";
      EXPECT_NEAR(Order(heavy_deviation[i], heavy_deviation[i + 1], 2.0), 1.0, 0.15) << "Heavy product not O(s)";
    }
    EXPECT_LT(std::abs(enrichment.back()), 0.01 * enrichment[swings.size() - 2])
      << "Without a pressure swing the cycle still separates";
  }

  // ===========================================================================
  // 7. Trace feed: the single-component (linear chromatography) limit
  // ===========================================================================
  //
  // O2 fed at a mole fraction y into a bed of N2, both adsorbing on a
  // competitive (extended) Langmuir isotherm. As y -> 0 (a) the carrier's
  // fields must tend to those of the single-component run (pure N2 feed),
  // with an error O(y); (b) the O2 outlet curve normalised by y must tend to
  // a limit curve, O(y); and (c) that limit is linear chromatography with the
  // trace's Henry constant at the carrier's loading,
  //
  //     H_eff = rho_b m_O2 b'_O2 R T / (1 + b'_N2 p_N2),
  //
  // so its first moment is t_st = (L/u_s)(et + H_eff). Orders within 0.1
  // (next term O(y) relative, y <= 0.04); the moment, Richardson-extrapolated
  // to y = 0, within 0.2 % (the trapezoidal quadrature of the logged curve).

  TEST(ParameterLimits, TraceFeedReducesToTheSingleComponentLimit) {
    auto make = [](UserData& data, double y) {
      BreakthroughCase c;
      c.adsorbing = true;
      c.y_feed = y;
      c.run_time = 400.0;
      c.reltol = 1e-8;
      Build(data, c);
      // Measured with upwind, whose O(y) front sharpening is the same on 80 and
      // 160 cells; WENO's differs slightly (order 1.10 against the 1 +- 0.1
      // asserted), so the run is upwind (the default is WENO).
      data.settings.advection_scheme_ = AdvectionScheme::Upwind;
      auto& source = PackingOf(data).adsorption_source_;
      auto isotherm = std::make_unique<ExtendedLangmuirIsotherm>();
      isotherm->params_["Oxygen"].m1 = 3.0;
      isotherm->params_["Oxygen"].b1 = 3.0e-11;
      isotherm->params_["Oxygen"].q1 = 20000.0;
      isotherm->params_["Nitrogen"].m1 = 3.0;
      isotherm->params_["Nitrogen"].b1 = 5.0e-11;
      isotherm->params_["Nitrogen"].q1 = 24000.0;
      source.AddIsotherm({"Oxygen", "Nitrogen"}, std::move(isotherm));
      auto ldf = std::make_unique<LinearDrivingForce>();
      auto mtc = std::make_unique<MassTransferCoefficientConstant>();
      mtc->constant_mtc_ = 0.5;
      ldf->mtc_model_ = ToResistances(std::move(mtc));  // generic model, as the example decks store it
      source.mass_transfer_models_["Nitrogen"] = std::move(ldf);
      source.heat_of_adsorption_["Nitrogen"] = 0.0;
      return c;
    };
    struct Run { Curve curve; std::vector<double> q_n2; double m1 = 0.0; };
    auto run = [&](double y) {
      UserData data;
      const BreakthroughCase c = make(data, y);
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      CoupledSolver solver{data};
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0);
      const auto& log = solver.GetLog();
      Run r;
      std::tie(r.curve.t, r.curve.y) = OutletSeries(log, false);
      const auto n = log.GetFaceVariable(0, 0, VariableType::MOLAR_FLOW_RATE);
      std::vector<double> unbroken(r.curve.t.size(), 0.0);
      if (y > 0.0) {
        for (size_t k = 0; k < r.curve.t.size(); ++k) {
          unbroken[k] = 1.0 - n.values[k].back() * r.curve.y[k] / (c.molar_flow * y);
          r.curve.y[k] /= y;
        }
        r.m1 = TimeIntegral(r.curve.t, unbroken);
      }
      const auto q = log.GetComponentVariable(0, 0, ComponentVariableType::ADSORBED_CONCENTRATION, "Nitrogen").values.back();
      r.q_n2.assign(q.begin() + 1, q.end() - 1);
      return r;
    };

    const Run pure = run(0.0);  // the single-component run
    const std::vector<double> ys{0.04, 0.02, 0.01};
    std::vector<Run> runs;
    for (double y : ys) runs.push_back(run(y));

    std::vector<double> carrier_error, curve_step;
    for (size_t i = 0; i < runs.size(); ++i) {
      double worst = 0.0, scale = 0.0;
      for (size_t j = 0; j < pure.q_n2.size(); ++j) {
        worst = std::max(worst, std::abs(runs[i].q_n2[j] - pure.q_n2[j]));
        scale = std::max(scale, pure.q_n2[j]);
      }
      carrier_error.push_back(worst / scale);
      if (i + 1 < runs.size()) curve_step.push_back(MaxDifference(runs[i].curve, runs[i + 1].curve, 400.0));
      std::cout << "[trace] y " << ys[i] << ": N2 loading vs pure N2 (relative) " << carrier_error.back()
                << ", O2 first moment " << runs[i].m1 << " s\n";
    }
    ASSERT_GT(pure.q_n2.front(), 0.01) << "The carrier does not adsorb";
    for (size_t i = 0; i + 1 < carrier_error.size(); ++i)
      EXPECT_NEAR(Order(carrier_error[i], carrier_error[i + 1], 2.0), 1.0, 0.1) << "Carrier not O(y)";
    ASSERT_GT(curve_step.front(), 1e-3);
    EXPECT_NEAR(Order(curve_step[0], curve_step[1], 2.0), 1.0, 0.1) << "Normalised trace curve not O(y)";

    // (c) The limit curve's first moment, extrapolated to y = 0 (it is O(y)).
    const double T = 298.0;
    const double H_eff = (1.0 - kCaseVoidFraction) * kCaseSolidDensity * 3.0 * 3.0e-11 *
                         std::exp(20000.0 / (kR * T)) * kR * T /
                         (1.0 + 5.0e-11 * std::exp(24000.0 / (kR * T)) * kCasePressure);
    const BreakthroughCase c;
    const double t_st = (c.length / CaseSuperficialVelocity(c)) * (kCaseVoidFraction + H_eff);
    const double m1_limit = 2.0 * runs[2].m1 - runs[1].m1;
    std::cout << "[trace] first moment extrapolated to y = 0: " << m1_limit << " s; linear chromatography "
              << t_st << " s\n";
    ASSERT_GT(H_eff, 0.5 * kCaseVoidFraction) << "The trace is barely retained";
    EXPECT_NEAR(m1_limit, t_st, 2e-3 * t_st);
  }

  // ===========================================================================
  // 8. Degenerate mixture: identical components are not separated
  // ===========================================================================
  //
  // Split the shared PSA cycle's O2 into two components with identical
  // isotherm (competitive Langmuir over all three), kinetics and heat of
  // adsorption: O2 and Ar, fed and initialised in the ratio 2:1. Their
  // species equations are then the same linear equation (the coefficients
  // depend only on the shared P, T, u and the competitive denominator), so
  // y_O2 = 2 y_Ar and q_O2 = 2 q_Ar must hold in every cell at every instant,
  // through pressurisation, feed and reversed-flow blowdown, while the pair
  // is separated from N2. (Their gas properties, molar mass, Cp, viscosity,
  // differ, but they enter only the shared mixture fields.) Tolerance: one
  // error-test weight, reltol |x| + atol (atol 1e-8 for both, docs/numerics.md).

  namespace {
    void SetAllMoleFractions(UserData& data, const std::map<std::string, double>& y) {
      for (const auto& [c, v] : y) data.ic.mole_fraction[c] = v;
      for (auto& step : data.cycle.step_list_) {
        for (auto* bc : {step.top_boundary_condition_.get(), step.bottom_boundary_condition_.get()}) {
          if (auto* f = dynamic_cast<FlowBoundaryCondition*>(bc))
            for (auto& s : f->streams_) for (const auto& [c, v] : y) s->mole_fraction[c] = v;
          if (auto* p = dynamic_cast<PressureBoundaryCondition*>(bc))
            for (const auto& [c, v] : y) p->mole_fraction_[c] = v;
        }
      }
    }
  } // namespace

  TEST(DegenerateMixtures, IdenticalComponentsKeepTheirFeedRatio) {
    UserData data;
    BuildPsaCycle(data, 30);
    // Upwind faces are linear in the cell values, so identical components
    // keep their ratio to round-off. WENO's smoothness weights use an
    // absolute epsilon (1e-12), so two identical components at different
    // mole fractions get slightly different weights and their ratio drifts by
    // 5.4e-6 in the gas (3.3e-7 in the loadings) here: a property of the
    // scheme, documented in docs/numerics.md.
    data.settings.advection_scheme_ = AdvectionScheme::Upwind;
    std::string error;
    ASSERT_TRUE(AddComponentToDeck(data, "Argon", error)) << error;
    auto& source = PackingOf(data).adsorption_source_;
    auto isotherm = std::make_unique<ExtendedLangmuirIsotherm>();
    for (const char* c : {"Oxygen", "Argon"}) {
      isotherm->params_[c].m1 = 3.0;
      isotherm->params_[c].b1 = 3.0e-11;
      isotherm->params_[c].q1 = 20000.0;
    }
    isotherm->params_["Nitrogen"].m1 = 3.0;
    isotherm->params_["Nitrogen"].b1 = 5.0e-11;
    isotherm->params_["Nitrogen"].q1 = 24000.0;
    source.AddIsotherm({"Oxygen", "Nitrogen", "Argon"}, std::move(isotherm));
    auto ldf = std::make_unique<LinearDrivingForce>();
    auto mtc = std::make_unique<MassTransferCoefficientConstant>();
    mtc->constant_mtc_ = 0.5;  // = O2's in the PSA case
    ldf->mtc_model_ = ToResistances(std::move(mtc));  // generic model, as the example decks store it
    source.mass_transfer_models_["Argon"] = std::move(ldf);
    source.heat_of_adsorption_["Argon"] = 20000.0;  // = O2's
    SetAllMoleFractions(data, {{"Nitrogen", 0.79}, {"Oxygen", 0.14}, {"Argon", 0.07}});
    const double reltol = 1e-8;
    data.cycle.number_of_cycles_ = 3;
    data.settings.reltol = reltol;
    data.settings.max_step = 2.0;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    CoupledSolver solver{data};
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);

    const auto& log = solver.GetLog();
    double worst_y = 0.0, worst_q = 0.0, y_min = 1.0, y_max = 0.0;
    bool y_ok = true, q_ok = true;
    for (int cycle : log.sys_data.cycle_index_set_) {
      for (int step : log.sys_data.step_index_set_) {
        const auto yo = log.GetComponentVariable(cycle, step, ComponentVariableType::MOLE_FRACTION, "Oxygen");
        const auto ya = log.GetComponentVariable(cycle, step, ComponentVariableType::MOLE_FRACTION, "Argon");
        const auto qo = log.GetComponentVariable(cycle, step, ComponentVariableType::ADSORBED_CONCENTRATION, "Oxygen");
        const auto qa = log.GetComponentVariable(cycle, step, ComponentVariableType::ADSORBED_CONCENTRATION, "Argon");
        for (size_t n = 0; n < yo.values.size(); ++n) {
          for (size_t i = 1; i + 1 < yo.values[n].size(); ++i) {
            const double dy = std::abs(yo.values[n][i] - 2.0 * ya.values[n][i]);
            const double dq = std::abs(qo.values[n][i] - 2.0 * qa.values[n][i]);
            worst_y = std::max(worst_y, dy);
            worst_q = std::max(worst_q, dq);
            y_ok = y_ok && dy <= reltol * std::abs(yo.values[n][i]) + 1e-8;
            q_ok = q_ok && dq <= reltol * std::abs(qo.values[n][i]) + 1e-8;
            y_min = std::min(y_min, yo.values[n][i]);
            y_max = std::max(y_max, yo.values[n][i]);
          }
        }
      }
    }
    std::cout << "[identical components] max |y_O2 - 2 y_Ar| " << worst_y << ", max |q_O2 - 2 q_Ar| "
              << worst_q << " mol/kg; y_O2 spans " << y_min << " .. " << y_max << "\n";
    // Guard: the pair is genuinely separated from N2 (y_O2 = 0.14 in the feed).
    ASSERT_LT(y_min, 0.10);
    ASSERT_GT(y_max, 0.20);
    EXPECT_TRUE(y_ok) << "Identical components left their feed ratio (worst " << worst_y << ")";
    EXPECT_TRUE(q_ok) << "Identical components' loadings left their feed ratio (worst " << worst_q << ")";
  }

  // ===========================================================================
  // 9. Degenerate mixture: a component with zero affinity is an inert
  // ===========================================================================
  //
  // A Henry isotherm with k = 0 for the tracer, and a large heat of adsorption
  // (which must have no effect, since nothing adsorbs), must reproduce the
  // inert-tracer run: the loading stays at zero (to its absolute tolerance)
  // and the end-of-run profile matches to within E_inert + E_zero + w, as in
  // test_metamorphic_relations.cpp (E: each run's error, measured by
  // tightening its reltol 100x; w: one error-test weight). The two runs do
  // not share their unknowns (an Inert loading is not a DAE unknown, a Henry
  // k = 0 loading is), so their integration errors are not correlated and
  // each is measured; the bound is the same metamorphic argument as 2E + w
  // for two runs with one layout, not a looser constant.

  TEST(DegenerateMixtures, ZeroAffinityComponentBehavesAsAnInert) {
    const double reltol = 1e-6;
    BreakthroughCase c;
    c.run_time = 15.0;  // the inert front mid-column (residence time ~25 s)
    c.reltol = reltol;
    auto inert_run = [&](double rt) {
      BreakthroughCase ci = c;
      ci.adsorbing = false;
      ci.reltol = rt;
      UserData data;
      Build(data, ci);
      return RunCase(data, ci);
    };
    const Curve inert = inert_run(reltol);
    const Curve inert_fine = inert_run(reltol * 1e-2);

    // The zero-affinity run at a relative tolerance: its final Oxygen
    // profile (interior cells) and its largest |q| over the run.
    auto zero_run = [&](double rt, double& q_max) {
      BreakthroughCase cz = c;
      cz.adsorbing = true;
      cz.reltol = rt;
      UserData zero;
      Build(zero, cz);
      auto& source = PackingOf(zero).adsorption_source_;
      auto henry = std::make_unique<HenryIsotherm>();
      henry->params_["Oxygen"].k = 0.0;
      henry->params_["Oxygen"].Q = kCaseHenryQ;
      source.AddIsotherm({"Oxygen"}, std::move(henry));
      source.heat_of_adsorption_["Oxygen"] = 40000.0;
      std::string error;
      EXPECT_TRUE(zero.PerformChecks(error)) << error;
      CoupledSolver solver{zero};
      EXPECT_GE(solver.IntegrateCycle(zero.cycle, false), 0);
      q_max = 0.0;
      for (const auto& s : solver.GetLog().GetComponentVariable(0, 0, ComponentVariableType::ADSORBED_CONCENTRATION, "Oxygen").values)
        for (size_t i = 1; i + 1 < s.size(); ++i) q_max = std::max(q_max, std::abs(s[i]));
      const auto y = solver.GetLog().GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION, "Oxygen").values.back();
      return std::vector<double>(y.begin() + 1, y.end() - 1);
    };
    double q_max = 0.0, q_max_fine = 0.0;
    const std::vector<double> y_zero = zero_run(reltol, q_max);
    const std::vector<double> y_zero_fine = zero_run(reltol * 1e-2, q_max_fine);

    ASSERT_EQ(y_zero.size(), inert.final_y.size());
    ASSERT_EQ(y_zero_fine.size(), inert.final_y.size());
    double discrepancy = 0.0, E = 0.0, E_zero = 0.0, y_hi = 0.0, y_lo = 1.0;
    for (size_t i = 0; i < y_zero.size(); ++i) {
      discrepancy = std::max(discrepancy, std::abs(y_zero[i] - inert.final_y[i]));
      E = std::max(E, std::abs(inert.final_y[i] - inert_fine.final_y[i]));
      E_zero = std::max(E_zero, std::abs(y_zero[i] - y_zero_fine[i]));
      y_hi = std::max(y_hi, inert.final_y[i]);
      y_lo = std::min(y_lo, inert.final_y[i]);
    }
    const double w = reltol * y_hi + 1e-8;
    std::cout << "[zero affinity] max |q| " << q_max << " mol/kg; profile discrepancy " << discrepancy
              << ", integration errors: inert run " << E << ", zero-affinity run " << E_zero << ", w " << w << "\n";
    ASSERT_GT(y_hi, 0.9 * c.y_feed) << "Front has not entered the bed";
    ASSERT_LT(y_lo, 0.1 * c.y_feed) << "Front has left the bed";
    EXPECT_LE(q_max, 1e-8) << "A component with no affinity was adsorbed";
    EXPECT_LE(discrepancy, E + E_zero + w) << "A zero-affinity component does not behave as an inert";
  }

} // namespace adsim
