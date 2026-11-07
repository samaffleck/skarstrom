/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Analytical validation of THERMAL FRONTS in a flowing bed.
//
// The energy tests in test_energy_balance seal the bed, so convection of heat
// never enters; the adiabatic-plateau test has convection, but its front speed
// is set by the adsorbate. These two tests pin the convective heat transport
// itself, and the bed -> wall -> ambient loss along a flowing column.
//
// 1. Inert-bed thermal front speed. A step in feed temperature T0 -> T1 into an
//    inert, adiabatic bed with no axial conduction travels as a front. The
//    total energy per bed volume is eps C Cv T + rho_b Cps T, and the gas part
//    eps C Cv T = eps Cv P / R is CONSTANT at constant pressure, so across the
//    front the energy and mole balances (molar flux G, interstitial gas
//    inventory eps C) read
//
//        rho_b Cps [T] w = Cp [G T],    eps [C] w = [G],
//
//    and with C T = P / R on both sides (Ruthven 1984, ch. 10; the same
//    algebra as a concentration shock) they give
//
//        w = G1 Cp / (rho_b Cps + eps C1 Cp),
//
//    with G1 and C1 = P / (R T1) on the feed side. The front reaches the outlet
//    at t = L / w, which is the first moment of the outlet temperature curve
//    (a mass/energy balance, so independent of the numerical smearing). The
//    eps C1 Cp term is the gas heat capacity: the test makes it ~7 % of the
//    total, so leaving it out (or using Cv) is caught.
//
//    The front is a contact discontinuity: with G2 from the mole balance the
//    characteristic speeds lambda = G Cp / (rho_b Cps + eps C Cp) on the two
//    sides are both exactly w (G2 (rho_b Cps + eps C1 Cp) = G1 (rho_b Cps +
//    eps C2 Cp)), so a numerically smeared front still travels at w and its
//    outlet moment is L / w. The test computes the two speeds and uses their
//    bracket |L/lambda1 - L/lambda2| (zero to rounding) as the physics part of
//    the tolerance; the rest, 1e-3, covers the pressure drop (C varies by
//    < 0.1 % of a 7 % term) and the moment quadrature.
//
// 2. Steady flow with wall losses (constant gas Cp: the CONSTANT mode). At
//    steady state, with the wall's axial
//    conduction off, the gas loses heat through the bed -> wall film (h_i, on
//    the inner surface pi d_i) and the wall -> ambient film (h_o, on pi d_o) in
//    series, so with an overall U, 1/(U pi d_i) = 1/(h_i pi d_i) + 1/(h_o pi d_o),
//
//        n_dot Cp dT/dz = -U pi d_i (T - T_amb)
//        T(z) = T_amb + (T_in - T_amb) exp(-lambda z),  lambda = U pi d_i / (n_dot Cp).
//
//    (n_dot Cp is constant along the bed at steady state: no adsorption, so
//    the molar flow is constant, and Cp is the cached molar value.) The
//    first-order upwind scheme turns the exponential into a geometric series
//    (1 + lambda dx)^(-i), so the error is O(lambda^2 dx); the test asserts
//    first-order convergence of the whole profile toward exp(-lambda z) and
//    that the series of films (not either film alone) sets lambda. The run
//    uses the upwind scheme for that reason: WENO and van Leer reconstruct
//    the face temperature too, and with WENO the outlet error converges at
//    second order (checked as well; the profile's max error stays first
//    order, set by the upwind inlet face where the profile has a slope).

#include "breakthrough_case.h"

#include <algorithm>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    double ObservedOrder(double e_coarse, double e_fine) {
      return std::log(e_coarse / e_fine) / std::log(2.0);
    }

    // Pure nitrogen everywhere (y_feed = 0): one gas, constant Cp.
    BreakthroughCase InertCase(int n_cells, double run_time) {
      BreakthroughCase c;
      c.n_cells = n_cells;
      c.adsorbing = false;
      c.y_feed = 0.0;
      c.run_time = run_time;
      c.max_step = 1.0;
      c.reltol = 1e-9;
      return c;
    }

  } // namespace

  // ===========================================================================
  // 1. Inert-bed thermal front speed
  // ===========================================================================

  TEST(ThermalFronts, InertBedFrontTravelsAtEnergyBalanceSpeed) {
    const double T0 = 298.0, T1 = 308.0;
    const double Cps = 10.0;  // [J/(kg-K)]: small, so the gas term is ~7 % of the total
    BreakthroughCase c = InertCase(100, 900.0);
    c.bed_temperature = T0;
    c.feed_temperature = T1;

    UserData data;
    Build(data, c);
    data.settings.max_log_points_per_step_ = 900;
    auto& layer = data.bed.layers_.layers_.front();
    auto& adsorbent = data.adsorbent_directory.Get(layer.composition_.begin()->first);
    adsorbent.heat_capacity_ = Cps;
    data.bed.equations_.wall_temperature_boundary_condition_ =
      CreateWallTemperatureBoundaryCondition(WallTemperatureBoundaryConditionType::ZERO_FLUX);

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
    const auto& log = solver.GetLog();
    const auto& dom = log.GetDomain();

    const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
    const double P_mean = InteriorMean(P.values.back(), dom.volume);
    const auto& Pf = P.values.back();
    ASSERT_LT(std::abs(Pf[1] - Pf[Pf.size() - 2]) / P_mean, 2e-3) << "Not isobaric";

    // Reference, from the inputs and the cached gas Cp.
    const double area = 0.25 * M_PI * c.diameter * c.diameter;
    const double G1 = c.molar_flow / area;                         // [mol/(m2 s)]
    const double Cp = MixtureMolarCp(data.fluid, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}});
    const double rho_b = (1.0 - kCaseVoidFraction) * kCaseSolidDensity;
    const double solid = rho_b * Cps;
    const double C1 = P_mean / (kR * T1), C2 = P_mean / (kR * T0);
    const double gas1 = kCaseVoidFraction * C1 * Cp;
    const double w = G1 * Cp / (solid + gas1);
    const double t_front = c.length / w;
    // Characteristic speeds on either side of the front; G2 from eps [C] w = [G].
    const double G2 = G1 + w * kCaseVoidFraction * (C2 - C1);
    const double lambda1 = G1 * Cp / (solid + gas1);
    const double lambda2 = G2 * Cp / (solid + kCaseVoidFraction * C2 * Cp);
    const double bracket = std::abs(c.length / lambda1 - c.length / lambda2);

    // GUARD: the gas heat capacity matters, and the run sees the whole front.
    const double t_no_gas = c.length * solid / (G1 * Cp);
    ASSERT_GT(std::abs(t_front - t_no_gas) / t_front, 0.05);
    ASSERT_LT(t_front, 0.5 * c.run_time);

    // Outlet temperature curve and its first moment.
    const auto T = log.GetFaceVariable(0, 0, VariableType::TEMPERATURE);
    std::vector<double> survival;
    for (const auto& snap : T.values) survival.push_back(1.0 - (snap.back() - T0) / (T1 - T0));
    ASSERT_LT(std::abs(survival.back()), 1e-4) << "Thermal front has not fully broken through";
    const double t_moment = TimeIntegral(T.times, survival);

    // Tolerance: the shock-vs-characteristic bracket (physics), plus 1e-3 for
    // the pressure-drop variation of C and the moment quadrature (1 s samples).
    const double tol = bracket + 1e-3 * t_front;
    std::cout << "[ThermalFronts] t_moment=" << t_moment << " s, L/w=" << t_front
              << " s, no-gas-term " << t_no_gas << " s, bracket " << bracket << " s\n";
    EXPECT_NEAR(t_moment, t_front, tol)
      << "Thermal front reaches the outlet at " << t_moment << " s, energy balance gives "
      << t_front << " s (without the gas heat capacity: " << t_no_gas << " s)";
  }

  // ===========================================================================
  // 2. Steady flow with wall losses: exponential axial temperature profile
  // ===========================================================================

  TEST(ThermalFronts, SteadyFlowWithWallLossGivesExponentialProfile) {
    const double T_amb = 298.0, T_in = 348.0;
    const double hi = 1.0, ho = 1.5;   // [W/(m2-K)]
    const double t_wall = 0.005;       // [m]

    struct Grid { int n; double err_profile, err_outlet; };
    std::vector<Grid> grids, weno;
    double lambda_L = 0.0, theta_out_ref = 0.0;

    for (AdvectionScheme scheme : {AdvectionScheme::Upwind, AdvectionScheme::WENO})
    for (int n : {25, 50, 100}) {
      BreakthroughCase c = InertCase(n, 3000.0);
      c.bed_temperature = T_amb;
      c.feed_temperature = T_in;
      UserData data;
      Build(data, c);
      // The exponential profile solves a linear ODE: it needs a constant gas
      // Cp. With the default ideal-gas C_p0(T) (Fluid) the profile is not
      // exactly exponential (N2's C_p0 changes by 0.3 % from 298 to 348 K), so
      // the gas uses the CONSTANT mode at N2's C_p0 of the mean temperature.
      data.fluid.gas_heat_capacity_mode_ = GasHeatCapacityMode::CONSTANT;
      for (const char* gas : {"Nitrogen", "Oxygen"})
        data.fluid.gas_heat_capacity_[gas] =
          Fluid::EvalIdealGasHeatCapacity(*Fluid::IdealGasCoeffs(gas), 0.5 * (T_amb + T_in));
      data.settings.advection_scheme_ = scheme;
      data.bed.ambient_temperature_ = T_amb;
      data.bed.wall_.thickness_ = t_wall;
      data.bed.wall_.heat_capacity_ = 5.0;          // light wall: fast to reach steady state
      data.bed.wall_.thermal_conductivity_ = 0.0;   // no axial conduction in the wall
      data.bed.equations_.wall_temperature_boundary_condition_ =
        CreateWallTemperatureBoundaryCondition(WallTemperatureBoundaryConditionType::ZERO_FLUX);
      auto& layer = data.bed.layers_.layers_.front();
      auto& adsorbent = data.adsorbent_directory.Get(layer.composition_.begin()->first);
      adsorbent.heat_capacity_ = 10.0;              // light packing: fast to reach steady state
      auto ihtc = std::make_unique<ConstantInternalHeatTransferModel>();
      ihtc->heat_transfer_coefficient_ = hi;
      layer.equations_.internal_heat_transfer_model_ = std::move(ihtc);
      auto ehtc = std::make_unique<ConstantExternalHeatTransferModel>();
      ehtc->heat_transfer_coefficient_ = ho;
      layer.equations_.external_heat_transfer_model_ = std::move(ehtc);

      CoupledSolver solver{data};
      std::string error;
      ASSERT_TRUE(data.PerformChecks(error)) << error;
      ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete (n = " << n << ")";
      const auto& log = solver.GetLog();
      const auto& dom = log.GetDomain();

      const double d_i = c.diameter, d_o = c.diameter + 2.0 * t_wall;
      const double UA = 1.0 / (1.0 / (hi * M_PI * d_i) + 1.0 / (ho * M_PI * d_o));  // [W/(m K)]
      const double Cp = MixtureMolarCp(data.fluid, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}});
      const double lambda = UA / (c.molar_flow * Cp);
      lambda_L = lambda * c.length;
      theta_out_ref = std::exp(-lambda_L);

      // Steady state reached: the profile no longer moves.
      const auto T = log.GetVariable(0, 0, VariableType::TEMPERATURE);
      const auto& Tf = T.values.back();
      const auto& Tp = T.values[T.values.size() - 2];
      double drift = 0.0;
      for (size_t i = 1; i + 1 < Tf.size(); ++i) drift = std::max(drift, std::abs(Tf[i] - Tp[i]));
      ASSERT_LT(drift, 1e-4 * (T_in - T_amb)) << "Not at steady state (n = " << n << ")";

      // Profile at the cell centres, and the outlet face.
      double err = 0.0;
      double z = 0.0;
      for (size_t i = 1; i + 1 < Tf.size(); ++i) {
        const double zc = z + 0.5 * dom.cell_width[i];
        z += dom.cell_width[i];
        const double theta = (Tf[i] - T_amb) / (T_in - T_amb);
        err = std::max(err, std::abs(theta - std::exp(-lambda * zc)));
      }
      const auto Tface = log.GetFaceVariable(0, 0, VariableType::TEMPERATURE);
      const double theta_out = (Tface.values.back().back() - T_amb) / (T_in - T_amb);
      (scheme == AdvectionScheme::Upwind ? grids : weno).push_back({n, err, std::abs(theta_out - theta_out_ref)});
      std::cout << "[ThermalFronts] " << (scheme == AdvectionScheme::Upwind ? "upwind" : "WENO") << " n=" << n << " lambda L=" << lambda_L << " theta_out=" << theta_out
                << " (exp " << theta_out_ref << ") max profile error " << err << "\n";
    }

    // GUARD: a strong decay, set by BOTH films (each alone would give a lambda
    // at least 1.6x off).
    ASSERT_GT(lambda_L, 0.5);
    ASSERT_LT(lambda_L, 3.0);

    for (size_t g = 1; g < grids.size(); ++g) {
      const double p_prof = ObservedOrder(grids[g - 1].err_profile, grids[g].err_profile);
      const double p_out = ObservedOrder(grids[g - 1].err_outlet, grids[g].err_outlet);
      EXPECT_NEAR(p_prof, 1.0, 0.15) << "Profile error order " << p_prof << " at n = " << grids[g].n;
      EXPECT_NEAR(p_out, 1.0, 0.15) << "Outlet error order " << p_out << " at n = " << grids[g].n;
    }
    // The leading upwind error at the outlet is theta (lambda L)(lambda dx)/2.
    const double predicted = theta_out_ref * lambda_L * (lambda_L / grids.back().n) / 2.0;
    EXPECT_NEAR(grids.back().err_outlet, predicted, 0.15 * predicted)
      << "Outlet error " << grids.back().err_outlet << " != leading upwind error " << predicted;
    // WENO's reconstructed face temperature: second-order outlet error
    // (measured 2.01 and 2.00), far below upwind's on the same grid.
    for (size_t g = 1; g < weno.size(); ++g) {
      const double p_out = ObservedOrder(weno[g - 1].err_outlet, weno[g].err_outlet);
      EXPECT_NEAR(p_out, 2.0, 0.2) << "WENO outlet error order " << p_out << " at n = " << weno[g].n;
    }
    EXPECT_LT(weno.back().err_outlet, 0.1 * grids.back().err_outlet);
  }

} // namespace adsim
