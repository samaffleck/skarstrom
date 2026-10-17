/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Analytical validation of the energy balance.
//
// Each test drives the solver into a regime where the coupled PDE system
// collapses to a problem with a closed-form (or exactly-solvable lumped)
// solution, and compares the simulated state against that solution. The
// expected values are derived from thermodynamics, never from a previous run.
//
// The recurring trick is a CLOSED, SPATIALLY UNIFORM bed: both ends shut, so
// convection and axial conduction vanish identically and every cell follows the
// same ordinary differential equation. Whatever the discretisation does in
// space is then irrelevant, and any discrepancy is a genuine defect in the
// energy equation itself.

#include "validation_helpers.h"

#include <array>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    // Geometry / packing shared by the closed-bed energy tests.
    constexpr double kLength = 1.0;    // [m]
    constexpr double kDiameter = 0.1;  // [m]
    constexpr double kBedVoidFraction = 0.4;
    constexpr double kSolidDensity = 1000.0;  // [kg/m3]
    constexpr int kNumberOfCells = 20;

  } // namespace

  // Fixture that builds a closed, adiabatic, non-adsorbing bed. Derived tests
  // switch individual physics back on so that exactly one mechanism is active.
  class ClosedBedEnergyTest : public UserDataTest {
  public:
    LayerID layer_id{};
    AdsorbentID ads_id{};

    // Builds the bed. `solid_heat_capacity` sets the thermal ballast of the
    // packing: small values make the gas-phase (and hence the P-V work) term a
    // large share of the total, which is what makes the Cp-vs-Cv distinction
    // observable.
    void BuildBed(double solid_heat_capacity) {
      AddComponents({"Nitrogen", "Oxygen"});

      auto ids = AddLayer("Inert Packing");
      layer_id = ids.first;
      ads_id = ids.second;

      auto& layer = GetLayer(layer_id);
      layer.domain_.length = kLength;
      layer.domain_.n_cells = kNumberOfCells;
      layer.eb_ = kBedVoidFraction;
      layer.particle_diameter_ = 0.002;
      std::get<CylinderBed>(data.bed.shape_).diameter = kDiameter;

      auto& adsorbent = GetAdsorbent(ads_id);
      adsorbent.solid_density_ = kSolidDensity;
      adsorbent.heat_capacity_ = solid_heat_capacity;
      adsorbent.thermal_conductivity_ = 0.0;
      // No macropore void: the total void fraction is then exactly eb, so the
      // analytical gas inventory needs no extra bookkeeping.
      adsorbent.ep_ = 0.0;
    }

    // Bulk density of the packing [kg adsorbent / m3 bed].
    static double BulkDensity() { return (1.0 - kBedVoidFraction) * kSolidDensity; }

    // Total bed volume [m3].
    double BedVolume() const {
      return 0.25 * M_PI * kDiameter * kDiameter * kLength;
    }

    // Make the bed adiabatic: no bed<->wall and no wall<->ambient coupling.
    void MakeAdiabatic() {
      AddConstantInternalHeatTransferCoefficientModel(0.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);
    }

    void SealBothEnds(const std::string& step_name) {
      AddClosedBC(BoundaryLocation::BOTTOM, step_name);
      AddClosedBC(BoundaryLocation::TOP, step_name);
    }

    // Tighten the integrator so that the comparison is against the analytical
    // solution and not against IDA's default error control.
    void UseTightTolerances(double max_step) {
      data.settings.reltol = 1e-10;
      data.settings.max_step = max_step;
    }
  };

  // ===========================================================================
  // 1. Quiescent state
  // ===========================================================================
  //
  // A closed, adiabatic bed with no source and no adsorption has zero right hand
  // side: the exact solution is the initial condition, held forever. Any drift
  // is a spurious source term somewhere in the residuals.

  TEST_F(ClosedBedEnergyTest, QuiescentClosedBedIsExactlyPreserved) {
    data.info.title = "Validation: quiescent closed bed";

    BuildBed(900.0);
    MakeAdiabatic();

    const double P0 = 5.0e5;  // [Pa]
    const double T0 = 320.0;  // [K]
    AddInitialConditions(P0, T0, T0, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}}, 0.0);

    AddStep("Hold", 500.0);
    SealBothEnds("Hold");
    ConfigureCycle(1);
    UseTightTolerances(5.0);

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    const auto& log = solver.GetLog();
    const auto T = log.GetVariable(0, 0, VariableType::TEMPERATURE);
    const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
    ASSERT_FALSE(T.values.empty());

    for (size_t i = 1; i + 1 < T.values.back().size(); ++i) {
      EXPECT_NEAR(T.values.back()[i], T0, 1e-6)
        << "Temperature drifted in cell " << i << " with no driving force";
      EXPECT_NEAR(P.values.back()[i], P0, 1e-6 * P0)
        << "Pressure drifted in cell " << i << " with no driving force";
    }
  }

  // ===========================================================================
  // 2. Electrical heating of a closed bed
  // ===========================================================================
  //
  // Closed, adiabatic, non-adsorbing bed with a uniform volumetric heat input Q
  // [W/m3]. Per unit bed volume the energy equation is
  //
  //     [ et * C * Cp + rho_b * Cps ] dT/dt  =  et * dP/dt + Q
  //
  // The bed is rigid and sealed, so the molar density C is a constant of the
  // motion and the ideal-gas law gives dP/dt = C * R * dT/dt. Substituting,
  //
  //     [ et * C * (Cp - R) + rho_b * Cps ] dT/dt = Q
  //     [ et * C * Cv       + rho_b * Cps ] dT/dt = Q
  //
  // Every coefficient is constant, so the exact solution is a straight line:
  //
  //     T(t) = T0 + Q * t / ( et * C0 * Cv + rho_b * Cps )
  //     P(t) = C0 * R * T(t)
  //
  // This pins down three things at once: the accumulation term, the P-V work
  // term (drop it and the ramp is ~Cv/Cp too shallow), and the ideal-gas
  // coupling between the pressure and temperature equations. The packing is
  // given a deliberately low heat capacity and the gas a high pressure so the
  // gas carries a large share of the thermal mass — otherwise the solid ballast
  // would mask any error in the gas terms.

  TEST_F(ClosedBedEnergyTest, ElectricHeatingGivesAnalyticTemperatureRamp) {
    data.info.title = "Validation: closed-bed electrical heating";

    const double Cps = 50.0;  // [J/(kg-K)] light packing so the gas term matters
    BuildBed(Cps);
    MakeAdiabatic();

    const double P0 = 20.0e5;  // [Pa]
    const double T0 = 300.0;   // [K]
    AddInitialConditions(P0, T0, T0, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}}, 0.0);

    const double heater_power = 100.0;  // [W]
    const double heating_time = 100.0;  // [s]
    AddStep("Heat", heating_time);
    SealBothEnds("Heat");

    EnergySource heater;
    heater.enabled_ = true;
    heater.power_ = heater_power;
    data.cycle.GetStep("Heat").energy_sources_[layer_id] = heater;

    ConfigureCycle(1);
    UseTightTolerances(1.0);

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    // --- analytical solution -------------------------------------------------
    const double Cv = MixtureMolarCv(data.fluid, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}});
    const double C0 = P0 / (kR * T0);
    const double volumetric_heat_capacity =
      kBedVoidFraction * C0 * Cv + BulkDensity() * Cps;   // [J/(m3-K)]
    const double q_volumetric = heater_power / BedVolume();  // [W/m3]
    const double dTdt = q_volumetric / volumetric_heat_capacity;  // [K/s]

    // Sanity check on the test design itself: the gas must carry enough of the
    // thermal mass for the work term to be observable.
    const double gas_share = kBedVoidFraction * C0 * Cv / volumetric_heat_capacity;
    ASSERT_GT(gas_share, 0.1) << "Test is not sensitive to the gas-phase energy terms";

    const auto& log = solver.GetLog();
    const auto T = log.GetVariable(0, 0, VariableType::TEMPERATURE);
    const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
    const auto& dom = log.GetDomain();
    ASSERT_FALSE(T.values.empty());

    std::vector<double> T_mean;
    std::vector<double> P_mean;
    T_mean.reserve(T.values.size());
    P_mean.reserve(P.values.size());
    for (size_t n = 0; n < T.values.size(); ++n) {
      T_mean.push_back(InteriorMean(T.values[n], dom.volume));
      P_mean.push_back(InteriorMean(P.values[n], dom.volume));
    }

    for (double t : {25.0, 50.0, 75.0, 100.0}) {
      const double expected = T0 + dTdt * t;
      const double actual = InterpolateSeries(T.times, T_mean, t);
      EXPECT_NEAR(actual, expected, 0.005 * (expected - T0))
        << "Temperature ramp deviates from Q*t/(et*C*Cv + rho_b*Cps) at t = " << t
        << " s (rise: simulated " << actual - T0 << " K, analytic " << expected - T0 << " K)";
    }

    // The rise must also be *inconsistent* with the Cp-based (no work term)
    // prediction, otherwise the test above would pass for the wrong reason.
    const double Cp = MixtureMolarCp(data.fluid, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}});
    const double dTdt_no_work =
      q_volumetric / (kBedVoidFraction * C0 * Cp + BulkDensity() * Cps);
    ASSERT_GT(std::abs(dTdt - dTdt_no_work) * heating_time, 0.5)
      << "Test cannot distinguish the P-V work term";

    // Rigid sealed vessel: the pressure must track C0 * R * T exactly.
    const double T_final = InterpolateSeries(T.times, T_mean, heating_time);
    const double P_final = InterpolateSeries(P.times, P_mean, heating_time);
    EXPECT_NEAR(P_final, C0 * kR * T_final, 1e-4 * P0)
      << "Sealed bed violates P = C0*R*T";
  }

  // ===========================================================================
  // 3. Heat-driven desorption in a closed bed
  // ===========================================================================
  //
  // Same closed adiabatic bed, but now the packing is a real adsorbent holding a
  // Henry-law loading, and the heater drives desorption. Because the bed stays
  // spatially uniform, the PDE system reduces exactly to three coupled ODEs in
  // (T, C, q) which this test integrates independently with a fine RK4 step.
  // That reduced model is written straight from the governing equations, not
  // from the solver's discretisation, so agreement exercises the sign and
  // magnitude of the heat-of-adsorption source, the loss of gas to the adsorbed
  // phase, and the LDF kinetics all at once.
  //
  //     et dC/dt = -rho_b dq/dt                          (moles are conserved)
  //     dq/dt    = k_ldf ( q*(P,T) - q )                 (linear driving force)
  //     [et C Cp + rho_b Cps] dT/dt = et dP/dt + dH rho_b dq/dt + Q
  //     P = C R T,  q* = k P exp(Qiso / R T)

  TEST_F(ClosedBedEnergyTest, HeatDrivenDesorptionMatchesLumpedModel) {
    data.info.title = "Validation: closed-bed heat-driven desorption";

    // Light packing: the latent heat of desorption is then a large share of the
    // energy budget rather than a rounding error on the solid's sensible heat.
    const double Cps = 50.0;  // [J/(kg-K)]
    BuildBed(Cps);
    MakeAdiabatic();

    const double P0 = 1.0e5;  // [Pa]
    const double T0 = 300.0;  // [K]

    // Henry isotherm q = k * p * exp(Qiso / RT) [mol/kg], with the isosteric
    // heat matching the heat released per mole adsorbed (van 't Hoff exact for
    // this form -- see test_thermodynamic_consistency.cpp).
    const double henry_k = 4.5e-11;  // [mol/(kg-Pa)] -> q0 ~ 0.1 mol/kg
    const double Qiso = 25000.0;     // [J/mol]
    const double k_ldf = 0.5;        // [1/s]

    auto isotherm = std::make_unique<HenryIsotherm>();
    isotherm->params_["Nitrogen"].k = henry_k;
    isotherm->params_["Nitrogen"].Q = Qiso;
    AddIsotherm({"Nitrogen"}, ads_id, std::move(isotherm));
    AddLinearDrivingForceModel(ads_id, "Nitrogen", k_ldf);
    SetHOA(ads_id, "Nitrogen", Qiso);

    AddInitialConditions(P0, T0, T0, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}}, 0.0);

    const double heater_power = 70.0;  // [W]
    const double heating_time = 200.0; // [s]
    AddStep("Regenerate", heating_time);
    SealBothEnds("Regenerate");

    EnergySource heater;
    heater.enabled_ = true;
    heater.power_ = heater_power;
    data.cycle.GetStep("Regenerate").energy_sources_[layer_id] = heater;

    ConfigureCycle(1);
    UseTightTolerances(1.0);

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    // --- independent lumped reference solution -------------------------------
    const double Cp = MixtureMolarCp(data.fluid, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}});
    const double rho_b = BulkDensity();
    const double et = kBedVoidFraction;
    const double q_volumetric = heater_power / BedVolume();

    auto equilibrium_loading = [&](double P, double T) {
      return henry_k * P * std::exp(Qiso / (kR * T));
    };

    // State: y = (T, C, q).
    auto derivatives = [&](const std::array<double, 3>& y) {
      const double T = y[0], C = y[1], q = y[2];
      const double P = C * kR * T;
      const double dqdt = k_ldf * (equilibrium_loading(P, T) - q);
      const double dCdt = -rho_b * dqdt / et;
      // dP/dt = R (C dT/dt + T dC/dt); substitute into the energy equation and
      // solve for dT/dt.
      //   [et C Cp + rho_b Cps] dT/dt - et R C dT/dt = et R T dC/dt + dH rho_b dq/dt + Q
      const double lhs_coeff = et * C * Cp + rho_b * Cps - et * kR * C;
      const double rhs = et * kR * T * dCdt + Qiso * rho_b * dqdt + q_volumetric;
      return std::array<double, 3>{rhs / lhs_coeff, dCdt, dqdt};
    };

    std::array<double, 3> y{T0, P0 / (kR * T0), equilibrium_loading(P0, T0)};
    const double dt = 1e-3;
    const int n_steps = static_cast<int>(heating_time / dt);
    for (int n = 0; n < n_steps; ++n) {
      const auto k1 = derivatives(y);
      std::array<double, 3> y2{}, y3{}, y4{};
      for (int j = 0; j < 3; ++j) y2[j] = y[j] + 0.5 * dt * k1[j];
      const auto k2 = derivatives(y2);
      for (int j = 0; j < 3; ++j) y3[j] = y[j] + 0.5 * dt * k2[j];
      const auto k3 = derivatives(y3);
      for (int j = 0; j < 3; ++j) y4[j] = y[j] + dt * k3[j];
      const auto k4 = derivatives(y4);
      for (int j = 0; j < 3; ++j) {
        y[j] += dt / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
      }
    }
    const double T_expected = y[0];
    const double q_expected = y[2];
    const double P_expected = y[1] * kR * y[0];

    // --- comparison ----------------------------------------------------------
    const auto& log = solver.GetLog();
    const auto& dom = log.GetDomain();
    const auto T = log.GetVariable(0, 0, VariableType::TEMPERATURE);
    const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
    const auto q = log.GetComponentVariable(0, 0, ComponentVariableType::ADSORBED_CONCENTRATION,
                                            "Nitrogen");
    ASSERT_FALSE(T.values.empty());
    ASSERT_FALSE(q.values.empty());

    const double T_final = InteriorMean(T.values.back(), dom.volume);
    const double P_final = InteriorMean(P.values.back(), dom.volume);
    const double q_final = InteriorMean(q.values.back(), dom.volume);

    // Desorption must actually have happened, and the latent heat must be a
    // meaningful share of the energy budget -- otherwise the test would pass
    // even if the heat-of-adsorption source term were dropped entirely.
    const double q_initial = InteriorMean(q.values.front(), dom.volume);
    ASSERT_LT(q_final, 0.8 * q_initial) << "No appreciable desorption occurred";
    const double latent_heat = Qiso * rho_b * (q_initial - q_expected) * BedVolume();
    ASSERT_GT(latent_heat, 0.2 * heater_power * heating_time)
      << "Latent heat is too small a share of the energy input for this test to "
         "discriminate the heat-of-adsorption source term";

    EXPECT_NEAR(T_final, T_expected, 0.01 * (T_expected - T0))
      << "Final temperature disagrees with the lumped energy balance "
      << "(simulated " << T_final << " K, reference " << T_expected << " K)";
    EXPECT_NEAR(q_final, q_expected, 0.01 * std::abs(q_initial - q_expected))
      << "Final loading disagrees with the lumped mass balance "
      << "(simulated " << q_final << " mol/kg, reference " << q_expected << " mol/kg)";
    EXPECT_NEAR(P_final, P_expected, 0.01 * std::abs(P_expected - P0))
      << "Final pressure disagrees with the lumped model";
  }

  // ===========================================================================
  // 4. Global energy balance with wall and ambient losses
  // ===========================================================================
  //
  // The bed<->wall<->ambient chain is the one part of the energy model with no
  // closed-form answer, but it must still obey conservation. For a closed bed
  // with no source and no adsorption, first-law bookkeeping over the whole
  // column (bed + wall) gives, with no approximation:
  //
  //     dU_bed + dU_wall + Q_ambient = 0
  //
  // and the internal exchange must be a pure transfer:
  //
  //     Q_bed->wall (as seen by the bed) == Q_bed->wall (as seen by the wall)
  //
  // The reference fluxes here use the textbook cylinder areas -- pi*di*dx on the
  // inside, pi*do*dx on the outside -- which is what the heat transfer
  // coefficients are defined against.
  //
  // This test originally failed by ~60 %: UpdateWallTemperature normalised its
  // fluxes by ((d_i + t_w)^2 - d_i^2) rather than the annulus area factor
  // (d_o^2 - d_i^2), and wrote flux_i in diameter form but flux_o in radius
  // form, so the wall absorbed ~2x the heat the bed gave up and rejected ~0.5x
  // the physical Newton-cooling flux. Fixed in adsorbent_equations.cpp.

  TEST_F(ClosedBedEnergyTest, ClosedBedWallAndAmbientEnergyBalanceCloses) {
    data.info.title = "Validation: bed/wall/ambient energy balance";

    const double Cps = 900.0;  // [J/(kg-K)]
    BuildBed(Cps);

    const double hi = 60.0;  // [W/(m2-K)] bed -> wall
    const double ho = 12.0;  // [W/(m2-K)] wall -> ambient
    AddConstantInternalHeatTransferCoefficientModel(hi, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(ho, layer_id);

    const double P0 = 1.0e5;
    const double T0 = 400.0;   // hot bed
    const double T_ambient = 300.0;
    data.bed.ambient_temperature_ = T_ambient;
    data.bed.wall_.thickness_ = 0.005;
    data.bed.wall_.density_ = 7900.0;
    data.bed.wall_.heat_capacity_ = 500.0;
    data.bed.wall_.thermal_conductivity_ = 16.2;

    // Insulate the wall ends. The default wall boundary condition pins the ghost
    // cells to a FIXED temperature (293 K), which drains energy axially out of
    // the end cells and would leave the control volume open.
    data.bed.equations_.wall_temperature_boundary_condition_ =
      CreateWallTemperatureBoundaryCondition(WallTemperatureBoundaryConditionType::ZERO_FLUX);

    AddInitialConditions(P0, T0, T0, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}}, 0.0);

    AddStep("Cool", 600.0);
    SealBothEnds("Cool");
    ConfigureCycle(1);
    UseTightTolerances(2.0);

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    const auto& log = solver.GetLog();
    const double Cv = MixtureMolarCv(data.fluid, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}});

    const double dU_bed = ClosedBedInternalEnergyChange(log, 0, 0, Cv);
    const double dU_wall = WallInternalEnergyChange(log, 0, 0);
    const double Q_ambient = AmbientHeatLoss(log, 0, 0, T_ambient);
    const double Q_bed_to_wall = BedToWallHeatTransfer(log, 0, 0);

    // The bed must actually cool appreciably, or the balance is trivially met.
    ASSERT_LT(dU_bed, -100.0) << "Bed did not cool; test proves nothing";

    // (a) The bed's loss must equal what the wall receives.
    EXPECT_NEAR(-dU_bed, Q_bed_to_wall, 0.02 * std::abs(dU_bed))
      << "Bed internal-energy loss does not match the integrated bed->wall flux "
      << "(dU_bed = " << dU_bed << " J, Q_bed->wall = " << Q_bed_to_wall << " J)";

    // (b) Closed system: everything the bed and wall lose must leave to ambient.
    const double imbalance = dU_bed + dU_wall + Q_ambient;
    const double scale = std::max({std::abs(dU_bed), std::abs(dU_wall), std::abs(Q_ambient)});
    EXPECT_NEAR(imbalance, 0.0, 0.02 * scale)
      << "Global energy balance does not close: dU_bed = " << dU_bed
      << " J, dU_wall = " << dU_wall << " J, Q_ambient = " << Q_ambient
      << " J, imbalance = " << imbalance << " J ("
      << 100.0 * imbalance / scale << "% of the largest term)";
  }

  // ===========================================================================
  // 5. Two-body lumped cooling: closed-form transient
  // ===========================================================================
  //
  // Conservation says the books balance; it does not say the wall responds at
  // the right rate. With the bed spatially uniform and the wall axially
  // insulated, the bed+wall pair is a two-state linear system whose solution is
  // a sum of two exponentials with analytically known rates:
  //
  //     d/dt [x1]   [ -a       a     ] [x1]        x1 = T  - T_amb
  //          [x2] = [  b   -(b + c)  ] [x2]        x2 = Tw - T_amb
  //
  //     a = (pi di h_i) / C_bed,   b = (pi di h_i) / C_wall,
  //     c = (pi do h_o) / C_wall            (all per unit column length)
  //
  // The eigenvalues of that 2x2 matrix set the fast (bed/wall equilibration) and
  // slow (loss to ambient) time constants. Getting either heat-transfer
  // coefficient, either heat capacity or either area wrong moves them.

  TEST_F(ClosedBedEnergyTest, WallAmbientCoolingMatchesTwoBodyExponentialSolution) {
    data.info.title = "Validation: two-body lumped cooling";

    const double Cps = 900.0;  // [J/(kg-K)]
    BuildBed(Cps);

    const double hi = 60.0;  // [W/(m2-K)]
    const double ho = 12.0;  // [W/(m2-K)]
    AddConstantInternalHeatTransferCoefficientModel(hi, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(ho, layer_id);

    const double P0 = 1.0e5;
    const double T0 = 400.0;
    const double T_ambient = 300.0;
    const double wall_thickness = 0.005;
    const double wall_density = 7900.0;
    const double wall_cp = 500.0;

    data.bed.ambient_temperature_ = T_ambient;
    data.bed.wall_.thickness_ = wall_thickness;
    data.bed.wall_.density_ = wall_density;
    data.bed.wall_.heat_capacity_ = wall_cp;
    data.bed.wall_.thermal_conductivity_ = 16.2;
    data.bed.equations_.wall_temperature_boundary_condition_ =
      CreateWallTemperatureBoundaryCondition(WallTemperatureBoundaryConditionType::ZERO_FLUX);

    AddInitialConditions(P0, T0, T0, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}}, 0.0);

    const double run_time = 1200.0;
    AddStep("Cool", run_time);
    SealBothEnds("Cool");
    ConfigureCycle(1);
    UseTightTolerances(2.0);

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    // --- closed-form two-body solution --------------------------------------
    const double Cv = MixtureMolarCv(data.fluid, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}});
    const double C0 = P0 / (kR * T0);
    const double d_i = kDiameter;
    const double d_o = kDiameter + 2.0 * wall_thickness;

    // Heat capacities per unit column length [J/(K-m)].
    const double C_bed = (kBedVoidFraction * C0 * Cv + BulkDensity() * Cps) *
                         0.25 * M_PI * d_i * d_i;
    const double C_wall = wall_density * wall_cp * 0.25 * M_PI * (d_o * d_o - d_i * d_i);

    const double a = M_PI * d_i * hi / C_bed;
    const double b = M_PI * d_i * hi / C_wall;
    const double c = M_PI * d_o * ho / C_wall;

    // Eigenvalues of [[-a, a], [b, -(b+c)]].
    const double trace = -(a + b + c);
    const double determinant = a * c;  // (-a)(-(b+c)) - (a)(b) = a*c
    const double disc = std::sqrt(trace * trace - 4.0 * determinant);
    const double lambda1 = 0.5 * (trace + disc);
    const double lambda2 = 0.5 * (trace - disc);
    ASSERT_LT(lambda1, 0.0);
    ASSERT_LT(lambda2, 0.0);

    // Eigenvectors: for eigenvalue L, (-a - L) v1 + a v2 = 0  =>  v = (a, a + L).
    // Initial condition x1(0) = x2(0) = T0 - T_amb fixes the coefficients.
    const double x0 = T0 - T_ambient;
    // Solve  k1*(a) + k2*(a)         = x0
    //        k1*(a + lambda1) + k2*(a + lambda2) = x0
    const double det_v = a * (a + lambda2) - a * (a + lambda1);
    ASSERT_GT(std::abs(det_v), 1e-12);
    const double k1 = (x0 * (a + lambda2) - a * x0) / det_v;
    const double k2 = (a * x0 - x0 * (a + lambda1)) / det_v;

    auto bed_temperature = [&](double t) {
      return T_ambient + k1 * a * std::exp(lambda1 * t) + k2 * a * std::exp(lambda2 * t);
    };
    auto wall_temperature = [&](double t) {
      return T_ambient + k1 * (a + lambda1) * std::exp(lambda1 * t) +
                         k2 * (a + lambda2) * std::exp(lambda2 * t);
    };

    // The two modes must be genuinely separated, or the test degenerates into a
    // single-exponential check that would not see a wrong wall heat capacity.
    // lambda1 is the slow (loss-to-ambient) mode, lambda2 the fast
    // (bed/wall equilibration) mode; they must be well separated or the test
    // degenerates into a single-exponential check that would not see a wrong
    // wall heat capacity.
    ASSERT_GT(lambda2 / lambda1, 5.0)
      << "Fast and slow modes are not separated; test is not discriminating";

    // --- comparison ----------------------------------------------------------
    const auto& log = solver.GetLog();
    const auto& dom = log.GetDomain();
    const auto T = log.GetVariable(0, 0, VariableType::TEMPERATURE);
    const auto Tw = log.GetVariable(0, 0, VariableType::WALL_TEMPERATURE);
    ASSERT_FALSE(T.values.empty());

    std::vector<double> T_mean, Tw_mean;
    for (size_t n = 0; n < T.values.size(); ++n) {
      T_mean.push_back(InteriorMean(T.values[n], dom.volume));
      Tw_mean.push_back(InteriorMean(Tw.values[n], dom.volume));
    }

    for (double t : {30.0, 120.0, 400.0, 800.0, 1200.0}) {
      const double T_expected = bed_temperature(t);
      const double Tw_expected = wall_temperature(t);
      EXPECT_NEAR(InterpolateSeries(T.times, T_mean, t), T_expected,
                  0.01 * (T_expected - T_ambient))
        << "Bed temperature deviates from the two-body solution at t = " << t << " s";
      EXPECT_NEAR(InterpolateSeries(Tw.times, Tw_mean, t), Tw_expected,
                  0.01 * (Tw_expected - T_ambient))
        << "Wall temperature deviates from the two-body solution at t = " << t << " s";
    }
  }

  // ===========================================================================
  // 6. Open-system energy balance
  // ===========================================================================
  //
  // With flow the control volume is open and the first law reads
  //
  //     dU  =  int ( n_dot_in h_in - n_dot_out h_out ) dt
  //
  // with h = Cp (T - T_ref) and u = Cv (T - T_ref) for an ideal gas. Unlike the
  // closed-bed tests this one is NOT expected to hold exactly on a coarse grid:
  // the temperature equation discretises convection in non-conservative form
  // (rho Cp u_avg (T_e - T_w) rather than a difference of face fluxes), so
  // neighbouring cells' fluxes do not telescope and a small residual survives.
  //
  // For a single-component, isobaric feed the two forms differ by
  // T d(rho Cp u)/dz, which vanishes identically because the molar flux is
  // constant -- so the imbalance is pure discretisation error and must converge
  // to zero under refinement. That convergence is the assertion.

  TEST_F(ClosedBedEnergyTest, OpenSystemEnergyBalanceConvergesUnderGridRefinement) {
    data.info.title = "Validation: open-system energy balance";

    const double P_outlet = 1.0e5;
    const double T_bed = 300.0;
    const double T_feed = 360.0;   // hot feed drives a thermal front
    const double molar_flow = 0.01;
    const double run_time = 20000.0;  // long enough for the thermal front to break through

    std::vector<double> relative_imbalance;
    for (int n_cells : {20, 40, 80, 160}) {
      UserData local;
      local.bed.wall_.thickness_ = 0.0175;
      local.bed.wall_.heat_capacity_ = 537;
      local.bed.wall_.density_ = 7800;
      local.bed.wall_.thermal_conductivity_ = 16.2;

      local.fluid.AddComponent("Nitrogen");
      local.fluid.AddComponent("Oxygen");

      auto ads = local.adsorbent_directory.Add("Inert Packing");
      auto lid = local.bed.layers_.Add("Inert Packing");
      auto& layer = local.bed.layers_.Get(lid);
      layer.momentum_model_ = MomentumModel::KarmenCozeny;  // fixture predates Ergun: keep its Darcy / Kozeny-Carman pressure drop
      layer.domain_.length = kLength;
      layer.domain_.n_cells = n_cells;
      layer.eb_ = kBedVoidFraction;
      layer.particle_diameter_ = 0.003;  // low resistance: nearly isobaric
      layer.composition_[ads] = 1.0;
      std::get<CylinderBed>(local.bed.shape_).diameter = kDiameter;

      auto& adsorbent = local.adsorbent_directory.Get(ads);
      adsorbent.ep_ = 0.0;
      adsorbent.solid_density_ = kSolidDensity;
      adsorbent.heat_capacity_ = 900.0;
      adsorbent.thermal_conductivity_ = 0.0;

      // Adiabatic: the wall must not participate, so the only energy paths are
      // the two boundary enthalpy flows.
      auto ihtc = std::make_unique<ConstantInternalHeatTransferModel>();
      ihtc->heat_transfer_coefficient_ = 0.0;
      layer.equations_.internal_heat_transfer_model_ = std::move(ihtc);
      auto ehtc = std::make_unique<ConstantExternalHeatTransferModel>();
      ehtc->heat_transfer_coefficient_ = 0.0;
      layer.equations_.external_heat_transfer_model_ = std::move(ehtc);

      // No axial conduction. The inlet ghost cell sits at the feed temperature,
      // so a non-zero effective conductivity would pump heat across the bottom
      // face on top of the convected enthalpy -- a second energy path that the
      // two-term balance below does not model.
      auto heat_disp = std::make_unique<HeatDispersionModelConstant>();
      heat_disp->dispersion_coefficient_ = 0.0;
      layer.equations_.heat_dispersion_model_ = std::move(heat_disp);

      local.ic.pressure = P_outlet;
      local.ic.temperature = T_bed;
      local.ic.wall_temperature = T_bed;
      local.ic.velocity = 0.0;
      local.ic.mole_fraction["Nitrogen"] = 1.0;
      local.ic.mole_fraction["Oxygen"] = 0.0;

      local.cycle.AddStep("Feed");
      Step& step = local.cycle.GetStep("Feed");
      step.step_time_ = run_time;

      auto inlet = std::make_unique<BoundaryStream>();
      auto flow = std::make_unique<MolarFlowBoundaryValue>();
      static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = molar_flow;
      inlet->boundary_flow_ = std::move(flow);
      inlet->mole_fraction["Nitrogen"] = 1.0;
      inlet->mole_fraction["Oxygen"] = 0.0;
      inlet->temperature_ = T_feed;
      auto inlet_bc = std::make_unique<FlowBoundaryCondition>();
      inlet_bc->direction_ = BoundaryDirection::INLET;
      inlet_bc->location_ = BoundaryLocation::BOTTOM;
      inlet_bc->streams_.push_back(std::move(inlet));
      step.bottom_boundary_condition_ = std::move(inlet_bc);

      auto outlet_bc = std::make_unique<PressureBoundaryCondition>();
      outlet_bc->direction_ = BoundaryDirection::OUTLET;
      outlet_bc->location_ = BoundaryLocation::TOP;
      auto outlet_p = std::make_unique<ConstantBoundaryValue>();
      outlet_p->constant_value_ = P_outlet;
      outlet_bc->pressure_ = std::move(outlet_p);
      outlet_bc->mole_fraction_["Nitrogen"] = 1.0;
      outlet_bc->mole_fraction_["Oxygen"] = 0.0;
      step.top_boundary_condition_ = std::move(outlet_bc);

      local.cycle.number_of_cycles_ = 1;
      local.settings.reltol = 1e-10;
      local.settings.max_step = 20.0;

      CoupledSolver solver{local};
      std::string error;
      ASSERT_TRUE(local.PerformChecks(error)) << error;
      ASSERT_GE(solver.IntegrateCycle(local.cycle, false), 0)
        << "Solver did not complete at n_cells = " << n_cells;

      const auto& log = solver.GetLog();
      const auto& dom = log.GetDomain();
      const double Cp = MixtureMolarCp(local.fluid, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}});
      const double Cv = Cp - kR;

      const auto T = log.GetVariable(0, 0, VariableType::TEMPERATURE);
      const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
      const auto n_dot = log.GetFaceVariable(0, 0, VariableType::MOLAR_FLOW_RATE);

      // Reference state at the bed's initial temperature, so U(t=0) = 0.
      const auto& T_end = T.values.back();
      const auto& P_end = P.values.back();
      double dU = 0.0;
      for (size_t i = 1; i + 1 < T_end.size(); ++i) {
        const double C = P_end[i] / (kR * T_end[i]);
        dU += dom.volume[i] * (dom.et[i] * C * Cv + dom.bulk_density[i] * dom.Cps[i]) *
              (T_end[i] - T_bed);
      }

      // Enthalpy carried in at the bottom face and out at the top face.
      const double H_in = molar_flow * Cp * (T_feed - T_bed) * run_time;
      std::vector<double> outflow_power(T.values.size(), 0.0);
      for (size_t n = 0; n < T.values.size(); ++n) {
        const double n_out = n_dot.values[n].back();          // [mol/s], + upwards
        const double T_out = T.values[n][T.values[n].size() - 2];
        outflow_power[n] = n_out * Cp * (T_out - T_bed);
      }
      const double H_out = TimeIntegral(T.times, outflow_power);

      ASSERT_GT(H_in, 0.0);
      // The outlet term must carry real weight, or only the inlet side is tested.
      EXPECT_GT(H_out, 0.1 * H_in)
        << "Thermal front has not broken through; the outlet enthalpy term is untested";

      relative_imbalance.push_back(std::abs(dU - (H_in - H_out)) / H_in);
    }

    ASSERT_EQ(relative_imbalance.size(), 4u);

    // Halving the cell width must halve the closure error: that is first-order
    // convergence, which is the expected order for the upwind discretisation.
    // Asserting the observed *order* is far stronger than asserting any single
    // tolerance -- it says the scheme is consistent, not merely accurate enough
    // on one grid.
    for (size_t i = 1; i < relative_imbalance.size(); ++i) {
      const double order_ratio = relative_imbalance[i - 1] / relative_imbalance[i];
      EXPECT_NEAR(order_ratio, 2.0, 0.3)
        << "Energy-balance closure error is not converging at first order: "
        << relative_imbalance[i - 1] << " -> " << relative_imbalance[i]
        << " (ratio " << order_ratio << ")";
    }

    EXPECT_LT(relative_imbalance.back(), 0.005)
      << "Open-system energy balance does not close on the finest grid (error = "
      << 100.0 * relative_imbalance.back() << "% of the enthalpy input)";
  }

  // ===========================================================================
  // 7. Cyclic steady state
  // ===========================================================================
  //
  // At cyclic steady state the state returns to where it started after a full
  // cycle, so the net energy accumulation over one cycle must be zero: whatever
  // the heater puts in must leave through the wall. The mass balance already has
  // this as a convergence criterion (CyclicSteadyState::CheckMassBalance...);
  // the energy balance had no equivalent.
  //
  // A sealed bed is used so the cycle has no flow enthalpy terms and the balance
  // is exact rather than discretisation-limited: heat in from the electric
  // heater, out through the wall to ambient.

  TEST_F(ClosedBedEnergyTest, SealedThermalCycleEnergyBalanceClosesAtCyclicSteadyState) {
    data.info.title = "Validation: cyclic energy balance";

    const double Cps = 900.0;
    BuildBed(Cps);

    const double hi = 80.0;
    const double ho = 15.0;
    AddConstantInternalHeatTransferCoefficientModel(hi, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(ho, layer_id);

    const double P0 = 1.0e5;
    const double T_ambient = 300.0;
    data.bed.ambient_temperature_ = T_ambient;
    data.bed.wall_.thickness_ = 0.005;
    data.bed.wall_.density_ = 7900.0;
    data.bed.wall_.heat_capacity_ = 500.0;
    data.bed.wall_.thermal_conductivity_ = 16.2;
    data.bed.equations_.wall_temperature_boundary_condition_ =
      CreateWallTemperatureBoundaryCondition(WallTemperatureBoundaryConditionType::ZERO_FLUX);

    AddInitialConditions(P0, T_ambient, T_ambient, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}}, 0.0);

    const double heater_power = 120.0;  // [W]
    const double heat_time = 300.0;
    const double cool_time = 900.0;

    AddStep("Heat", heat_time);
    SealBothEnds("Heat");
    EnergySource heater;
    heater.enabled_ = true;
    heater.power_ = heater_power;
    data.cycle.GetStep("Heat").energy_sources_[layer_id] = heater;

    AddStep("Cool", cool_time);
    SealBothEnds("Cool");

    const int n_cycles = 15;
    ConfigureCycle(n_cycles);
    data.settings.reltol = 1e-9;
    data.settings.max_step = 5.0;

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    const auto& log = solver.GetLog();
    const int last_cycle = *log.sys_data.cycle_index_set_.rbegin();
    const double Cv = MixtureMolarCv(data.fluid, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}});

    // Energy in over the cycle is the heater duty; energy out is the integrated
    // ambient loss over both steps. The state change over the cycle is the sum
    // of the two steps' accumulations.
    double dU = 0.0;
    double Q_ambient = 0.0;
    for (int step_index : log.sys_data.step_index_set_) {
      dU += ClosedBedInternalEnergyChange(log, last_cycle, step_index, Cv);
      dU += WallInternalEnergyChange(log, last_cycle, step_index);
      Q_ambient += AmbientHeatLoss(log, last_cycle, step_index, T_ambient);
    }
    const double Q_heater = heater_power * heat_time;

    // Cyclic steady state: the residual state change must be small compared with
    // the energy that flows through the cycle.
    EXPECT_LT(std::abs(dU), 0.02 * Q_heater)
      << "Cycle " << last_cycle << " has not reached cyclic steady state "
      << "(net state change = " << dU << " J vs " << Q_heater << " J in)";

    // First law over the cycle, including the residual drift.
    const double imbalance = dU - (Q_heater - Q_ambient);
    EXPECT_NEAR(imbalance, 0.0, 0.02 * Q_heater)
      << "Cyclic energy balance does not close: heater in = " << Q_heater
      << " J, ambient out = " << Q_ambient << " J, state change = " << dU
      << " J, imbalance = " << imbalance << " J";
  }

  // ===========================================================================
  // 8. Live per-step energy-closure assertion
  // ===========================================================================
  //
  // The mass balance has a live, per-step conservation assertion applied to every
  // step of every run (IsMolesConserved / IsComponentMolesConserved in
  // user_data_test.h). Until now the energy balance had no equivalent live check:
  // the tests above each construct a bespoke closed form. This test adds the
  // energy analogue -- SealedBedEnergyImbalance in validation_helpers.h -- and
  // exercises it over a multi-step sealed run with a heater, wall and ambient
  // losses, exactly the ingredients the mole balance never sees.
  //
  // For a SEALED bed the first law over bed + wall is exact (no boundary
  // enthalpy): dU_bed + dU_wall + Q_ambient - Q_source = 0 every step. The check
  // reads the heater duty back from the logged ENERGY_SOURCE field, so it is
  // self-contained and would catch a sign error or missing term in any of the
  // accumulation, wall-coupling, ambient-loss or source contributions.

  TEST_F(ClosedBedEnergyTest, SealedBedPerStepEnergyBalanceCloses) {
    data.info.title = "Validation: live per-step energy balance";

    const double Cps = 900.0;
    BuildBed(Cps);

    const double hi = 70.0;  // [W/(m2-K)] bed -> wall
    const double ho = 14.0;  // [W/(m2-K)] wall -> ambient
    AddConstantInternalHeatTransferCoefficientModel(hi, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(ho, layer_id);

    const double P0 = 1.0e5;
    const double T_ambient = 300.0;
    data.bed.ambient_temperature_ = T_ambient;
    data.bed.wall_.thickness_ = 0.005;
    data.bed.wall_.density_ = 7900.0;
    data.bed.wall_.heat_capacity_ = 500.0;
    data.bed.wall_.thermal_conductivity_ = 16.2;
    data.bed.equations_.wall_temperature_boundary_condition_ =
      CreateWallTemperatureBoundaryCondition(WallTemperatureBoundaryConditionType::ZERO_FLUX);

    AddInitialConditions(P0, T_ambient, T_ambient, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}}, 0.0);

    // A heated step (source + wall loss active) and an unheated cooling step
    // (wall loss only), so the check must close with the source both present and
    // absent and with the internal energy rising and falling.
    const double heater_power = 110.0;  // [W]
    AddStep("Heat", 300.0);
    SealBothEnds("Heat");
    EnergySource heater;
    heater.enabled_ = true;
    heater.power_ = heater_power;
    data.cycle.GetStep("Heat").energy_sources_[layer_id] = heater;

    AddStep("Cool", 600.0);
    SealBothEnds("Cool");

    ConfigureCycle(1);
    UseTightTolerances(3.0);

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    const auto& log = solver.GetLog();
    const double Cv = MixtureMolarCv(data.fluid, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}});

    // The effect under test must be large: the heater duty read back from the log
    // must be a substantial number of joules, and the bed must actually warm and
    // then cool, or a passing balance proves nothing.
    const double Q_heat = EnergySourceHeat(log, 0, 0);
    ASSERT_GT(Q_heat, 1.0e4) << "Heater duty is negligible; the source term is untested";
    const double dU_heat = ClosedBedInternalEnergyChange(log, 0, 0, Cv);
    const double dU_cool = ClosedBedInternalEnergyChange(log, 0, 1, Cv);
    ASSERT_GT(dU_heat, 1.0e3) << "Bed did not warm during the heated step";
    ASSERT_LT(dU_cool, -1.0e3) << "Bed did not cool during the unheated step";

    // Live per-step energy closure over the whole run.
    const double imbalance = SealedBedEnergyImbalance(log, Cv, T_ambient);
    EXPECT_LT(imbalance, 0.01)
      << "Live per-step energy balance does not close: worst relative imbalance "
      << 100.0 * imbalance << "% of the largest term in some step";
  }

} // namespace adsim
