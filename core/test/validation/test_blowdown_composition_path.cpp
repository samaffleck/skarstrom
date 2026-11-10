/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Analytical validation of a PRESSURE-CHANGE step: the equilibrium-theory
// composition path of a binary mixture during blowdown.
//
// test_analytic_transport's valve blowdown pins the pressure-time history of an
// inert gas. With adsorption the more interesting quantity is the COMPOSITION:
// as the pressure falls, the strongly adsorbed component desorbs in larger
// proportion and the gas left in the column is enriched in it. That
// enrichment is what a PSA blowdown step delivers, and its closed form is the
// basis of the Knaebel-Hill and Chan-Hassan-Hill cycle analyses.
//
// Derivation (isothermal, local equilibrium, linear isotherms q_i = K_i p_i,
// negligible pressure drop). With one end closed and the bed initially uniform,
// every gas element sees the same pressure history P(t), so the composition
// stays uniform in z and the gas leaving has the bed's composition y (of the
// strong component A). Per unit bed volume each component's inventory is
// a_i P y_i, with
//
//     a_i = eps / (R T) + rho_b K_i    [mol/(m3 Pa)],
//
// and the outflow removes the components in the ratio y : 1 - y, so
//
//     d(a_A P y) / y = d(a_B P (1 - y)) / (1 - y).
//
// Expanding, a_A (dP/P + dy/y) = a_B (dP/P - dy/(1 - y)), i.e. with
// beta = a_B / a_A < 1,
//
//     d ln P = - [ 1/y + beta/(1 - y) ] dy / (1 - beta),
//
//     P / P0 = (y0 / y)^(1/(1-beta)) * ((1 - y) / (1 - y0))^(beta/(1-beta)).
//
// (Knaebel & Hill 1985, Chem. Eng. Sci. 40, 2351, eq. for the blowdown step, in
// their notation beta_i = 1 / (1 + rho_b K_i R T / eps).) The path depends only on
// the isotherm slopes, not on the valve, the flow or the time.
//
// What the engine has to get right: the gas and adsorbed-phase accumulation of
// both components at a changing pressure, and that the outflow carries the
// local composition. Neglected terms, each estimated and checked:
//
//   * LDF lag: q trails q* by (dq*/dt)/k, so ln P is off the path by
//     O(|d ln P/dt| / k). The test runs two LDF rates and asserts the deviation
//     falls as 1/k, and that the extrapolated (k -> inf) deviation is small;
//   * pressure drop: < 1e-4 of P here (Darcy, about 60 Pa at the peak flow);
//   * temperature: the isothermal mode holds T exactly.
//
// GUARDS: the enrichment is large (y from 0.3 to 0.76), and the ideal-gas
// (no adsorption, beta = 1) path would keep y constant, so the test cannot pass
// with the adsorbed phase missing or with the wrong isotherm slope.

#include "validation_helpers.h"

#include <algorithm>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    constexpr double kT = 298.0;
    constexpr double kEps = 0.4;
    constexpr double kRhoS = 1000.0;
    constexpr double kKA = 1.6e-6;    // [mol/(kg Pa)] strong: CO2
    constexpr double kKB = 2.0e-7;    // [mol/(kg Pa)] weak: N2
    constexpr double kP0 = 5.0e5;     // [Pa]
    constexpr double kPdown = 1.0e5;  // [Pa]
    constexpr double kY0 = 0.3;       // initial mole fraction of A
    constexpr double kCv = 5.0e-3;    // valve coefficient
    constexpr double kRunTime = 600.0;

    double RhoB() { return (1.0 - kEps) * kRhoS; }
    double Capacity(double K) { return kEps / (kR * kT) + RhoB() * K; }  // a_i
    double Beta() { return Capacity(kKB) / Capacity(kKA); }

    // ln(P/P0) on the equilibrium-theory path at composition y.
    double PathLogPressure(double y) {
      const double b = Beta();
      return -(std::log(y / kY0) - b * std::log((1.0 - y) / (1.0 - kY0))) / (1.0 - b);
    }

    struct Sample { double t, P, y, y_spread, dlnP_dt; };

    std::vector<Sample> RunBlowdown(double ldf_rate) {
      UserData data;
      data.bed.wall_.thickness_ = 0.005;
      data.fluid.AddComponent("Nitrogen");
      data.fluid.AddComponent("CarbonDioxide");
      auto ads = data.adsorbent_directory.Add("Packing");
      auto lid = data.bed.layers_.Add("Packing");
      auto& layer = data.bed.layers_.Get(lid);
      layer.momentum_model_ = MomentumModel::KarmenCozeny;
      layer.domain_.length = 1.0;
      layer.domain_.n_cells = 20;
      layer.eb_ = kEps;
      layer.particle_diameter_ = 0.002;
      layer.composition_[ads] = 1.0;
      std::get<CylinderBed>(data.bed.shape_).diameter = 0.1;
      auto& adsorbent = data.adsorbent_directory.Get(ads);
      adsorbent.ep_ = 0.0;
      adsorbent.solid_density_ = kRhoS;
      adsorbent.heat_capacity_ = 900.0;

      auto iso = std::make_unique<HenryIsotherm>();
      iso->params_["CarbonDioxide"].k = kKA;
      iso->params_["CarbonDioxide"].Q = 0.0;
      iso->params_["Nitrogen"].k = kKB;
      iso->params_["Nitrogen"].Q = 0.0;
      adsorbent.adsorption_source_.AddIsotherm({"CarbonDioxide", "Nitrogen"}, std::move(iso));
      for (const char* c : {"CarbonDioxide", "Nitrogen"}) {
        auto ldf = std::make_unique<LinearDrivingForce>();
        auto mtc = std::make_unique<MassTransferCoefficientConstant>();
        mtc->constant_mtc_ = ldf_rate;
        ldf->mtc_model_ = ToResistances(std::move(mtc));  // generic model, as the example decks store it
        adsorbent.adsorption_source_.mass_transfer_models_[c] = std::move(ldf);
        adsorbent.adsorption_source_.heat_of_adsorption_[c] = 0.0;
      }
      data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;

      data.ic.pressure = kP0;
      data.ic.temperature = kT;
      data.ic.wall_temperature = kT;
      data.ic.velocity = 0.0;
      data.ic.mole_fraction["Nitrogen"] = 1.0 - kY0;
      data.ic.mole_fraction["CarbonDioxide"] = kY0;

      data.cycle.AddStep("Blowdown");
      Step& step = data.cycle.GetStep("Blowdown");
      step.step_time_ = kRunTime;
      auto closed = std::make_unique<ClosedBoundaryCondition>();
      closed->location_ = BoundaryLocation::BOTTOM;
      step.bottom_boundary_condition_ = std::move(closed);
      auto valve = std::make_unique<ValveBoundaryStream>();
      auto p_down = std::make_unique<ConstantBoundaryValue>();
      p_down->constant_value_ = kPdown;
      valve->boundary_value_ = std::move(p_down);
      valve->valve_.cv_ = kCv;
      auto outlet_bc = std::make_unique<FlowBoundaryCondition>();
      outlet_bc->direction_ = BoundaryDirection::OUTLET;
      outlet_bc->location_ = BoundaryLocation::TOP;
      outlet_bc->streams_.push_back(std::move(valve));
      step.top_boundary_condition_ = std::move(outlet_bc);

      data.cycle.number_of_cycles_ = 1;
      data.settings.reltol = 1e-9;
      data.settings.max_step = 1.0;
      data.settings.max_log_points_per_step_ = 300;

      CoupledSolver solver{data};
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

      const auto& log = solver.GetLog();
      const auto& dom = log.GetDomain();
      const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
      const auto y = log.GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION, "CarbonDioxide");
      std::vector<Sample> out;
      for (size_t n = 0; n < P.values.size(); ++n) {
        const auto& yn = y.values[n];
        double lo = 1.0, hi = 0.0;
        for (size_t i = 1; i + 1 < yn.size(); ++i) { lo = std::min(lo, yn[i]); hi = std::max(hi, yn[i]); }
        out.push_back({P.times[n], InteriorMean(P.values[n], dom.volume),
                       InteriorMean(yn, dom.volume), hi - lo, 0.0});
      }
      for (size_t n = 1; n + 1 < out.size(); ++n)
        out[n].dlnP_dt = std::log(out[n + 1].P / out[n - 1].P) / (out[n + 1].t - out[n - 1].t);
      return out;
    }

  } // namespace

  TEST(BlowdownCompositionPath, FollowsEquilibriumTheoryPath) {
    // GUARD: a strongly selective pair, so the path is far from y = const.
    ASSERT_LT(Beta(), 0.3);

    // Two LDF rates: the path deviation must fall as 1/k.
    struct Result { double k, worst_dev, worst_lag, y_end, P_end, worst_spread; };
    std::vector<Result> results;
    for (double k : {5.0, 10.0}) {
      const auto samples = RunBlowdown(k);
      ASSERT_GT(samples.size(), 10u);
      Result r{k, 0.0, 0.0, samples.back().y, samples.back().P, 0.0};
      for (const auto& s : samples) {
        const double dev = std::log(s.P / kP0) - PathLogPressure(s.y);
        r.worst_dev = std::max(r.worst_dev, std::abs(dev));
        r.worst_lag = std::max(r.worst_lag, std::abs(s.dlnP_dt) / k);
        r.worst_spread = std::max(r.worst_spread, s.y_spread);
      }
      results.push_back(r);
      std::cout << "[BlowdownPath] k=" << k << " max|lnP - path(y)|=" << r.worst_dev
                << " max|dlnP/dt|/k=" << r.worst_lag << " y_end=" << r.y_end
                << " P_end=" << r.P_end << " max y spread=" << r.worst_spread << "\n";
    }

    for (const auto& r : results) {
      // GUARD: the blowdown went most of the way and enriched A strongly.
      ASSERT_LT(r.P_end, 1.1 * kPdown);
      ASSERT_GT(r.y_end, kY0 + 0.4) << "Enrichment too small to test the path";
      // The composition stayed uniform (the premise of the lumped derivation).
      EXPECT_LT(r.worst_spread, 0.01) << "Composition not uniform along the bed";
      // The deviation from the path is bounded by the LDF lag.
      EXPECT_LT(r.worst_dev, r.worst_lag)
        << "ln P departs from the equilibrium-theory path by " << r.worst_dev
        << " at k = " << r.k << ", more than the LDF lag bound " << r.worst_lag;
    }
    // ... and falls as 1/k: halving the lag halves the deviation.
    const double ratio = results[0].worst_dev / results[1].worst_dev;
    EXPECT_NEAR(ratio, 2.0, 0.3)
      << "Path deviation does not scale with the LDF lag 1/k (ratio " << ratio << ")";
    // Richardson-extrapolated to k -> inf, the deviation vanishes.
    const double extrapolated = 2.0 * results[1].worst_dev - results[0].worst_dev;
    EXPECT_LT(std::abs(extrapolated), 0.2 * results[1].worst_dev);

    // End-point check in mole-fraction units: y at the final pressure.
    const auto& r = results.back();
    double lo = kY0, hi = 0.999;  // solve PathLogPressure(y) = ln(P_end/P0) by bisection
    for (int it = 0; it < 100; ++it) {
      const double mid = 0.5 * (lo + hi);
      (PathLogPressure(mid) > std::log(r.P_end / kP0) ? lo : hi) = mid;
    }
    EXPECT_NEAR(r.y_end, 0.5 * (lo + hi), 2e-3)
      << "Final composition " << r.y_end << " != equilibrium-theory " << 0.5 * (lo + hi);
  }

} // namespace adsim
