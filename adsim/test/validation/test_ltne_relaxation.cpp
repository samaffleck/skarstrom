/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Analytical validation of the two-temperature (LTNE) gas-solid energy exchange.
//
// In local-thermal-non-equilibrium mode the solver carries SEPARATE fluid and
// solid temperature equations coupled by a volumetric heat-transfer term
// hp * (Tf - Ts):
//
//     et rho Cpg dV dTf/dt = -hp dV (Tf - Ts) + et dV dP/dt + Q      (gas)
//     rho_b Cps  dV dTs/dt = +hp dV (Tf - Ts)                        (solid)
//
// (adsorbent_equations.cpp, UpdateFluidTemperature / UpdateSolidTemperature; the
// gas keeps the P-V work term, the solid does not.) With the bed CLOSED, no
// adsorption, no wall coupling and spatially uniform, every cell obeys this same
// pair of ODEs and the discretisation cannot influence the answer.
//
// The test proceeds in two phases within one sealed run:
//
//   * A short heating phase deposits a known duty into the GAS only (the electric
//     energy source enters the fluid equation), driving Tf above Ts and creating
//     a genuine gas-solid disequilibrium -- the solver seeds Tf = Ts at t = 0, so
//     without this the pair would start already equilibrated and never relax.
//
//   * A relaxation phase (heater off, still sealed and adiabatic) in which the
//     two temperatures approach a common equilibrium exponentially. Because the
//     P-V work term makes dTf/dt use Cv rather than Cp (exactly as in the sealed
//     electric-heating energy test), the relaxation is the two-body linear system
//
//         d/dt (Tf - Ts) = -lambda (Tf - Ts),
//         lambda = hp V ( 1/C_gas + 1/C_solid ),
//
//     with C_gas = et C Cv V and C_solid = rho_b Cps V (per bed volume the V
//     cancels: lambda = hp ( 1/(et C Cv) + 1/(rho_b Cps) )). The equilibrium
//     temperature is the heat-capacity-weighted mean, fixed by conservation.
//
//   * Energy conservation is exact: the heater duty must equal the combined
//     sensible-energy rise of gas + solid over the whole run.
//
// The volumetric coefficient hp is built by the model as hp = h * av with
// av = 6 (1 - eb) / dp (ConstantParticleHeatTransferModel::Update), so the RATE
// this test recovers pins h, the specific area, AND the two heat capacities at
// once.
//
// GUARDS: (1) the gas-solid temperature split at the end of heating must be
// large, so there is a real disequilibrium to relax; (2) the relaxation must
// resolve most of the way to equilibrium; (3) the two heat capacities must be
// comparable, so lambda genuinely depends on both (a wrong solid Cps could not
// hide behind a dominant gas term or vice versa).

#include "validation_helpers.h"

namespace adsim {

  using namespace adsim::validation;

  class LtneRelaxationTest : public UserDataTest {
  public:
    LayerID layer_id{};
    AdsorbentID ads_id{};

    static constexpr double kLength = 1.0;      // [m]
    static constexpr double kDiameter = 0.1;    // [m]
    static constexpr double kVoidFraction = 0.4;
    static constexpr double kParticleDiameter = 0.002;  // [m]
    static constexpr int kNumberOfCells = 20;

    // Chosen so the gas and solid heat capacities are of the same order (see the
    // guard): high pressure fattens the gas term, a light packing thins the solid
    // term.
    static constexpr double kPressure = 30.0e5;     // [Pa]
    static constexpr double kSolidDensity = 1000.0; // [kg/m3]
    static constexpr double kSolidCp = 60.0;        // [J/(kg-K)] light packing
    static constexpr double kParticleHtc = 8.0;     // [W/(m2-K)] gas-to-particle

    void BuildBed() {
      AddComponents({"Nitrogen", "Oxygen"});
      auto ids = AddLayer("Inert Packing");
      layer_id = ids.first;
      ads_id = ids.second;

      auto& layer = GetLayer(layer_id);
      layer.domain_.length = kLength;
      layer.domain_.n_cells = kNumberOfCells;
      layer.eb_ = kVoidFraction;
      layer.particle_diameter_ = kParticleDiameter;
      std::get<CylinderBed>(data.bed.shape_).diameter = kDiameter;

      auto& adsorbent = GetAdsorbent(ads_id);
      adsorbent.ep_ = 0.0;  // total void == bed void
      adsorbent.solid_density_ = kSolidDensity;
      adsorbent.heat_capacity_ = kSolidCp;
      adsorbent.thermal_conductivity_ = 0.0;

      // Adiabatic: no bed<->wall, no wall<->ambient. The solid can then only
      // exchange energy with the gas, and the gas only with the solid.
      AddConstantInternalHeatTransferCoefficientModel(0.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);

      // Two-temperature (LTNE) mode with a constant gas-to-particle coefficient.
      data.bed.equations_.thermal_equilibrium_mode_ =
        ThermalEquilibriumMode::NON_LOCAL;
      auto php = std::make_unique<ConstantParticleHeatTransferModel>();
      php->heat_transfer_coefficient_ = kParticleHtc;
      data.bed.equations_.particle_heat_transfer_model_ = std::move(php);
    }

    double BedVolume() const {
      return 0.25 * M_PI * kDiameter * kDiameter * kLength;
    }
    static double BulkDensity() { return (1.0 - kVoidFraction) * kSolidDensity; }
  };

  // ===========================================================================
  // Gas-solid two-temperature relaxation
  // ===========================================================================

  TEST_F(LtneRelaxationTest, GasAndSolidRelaxToCommonTemperatureAtKnownRate) {
    data.info.title = "Validation: LTNE gas-solid relaxation";

    BuildBed();

    const double T0 = 300.0;  // [K]
    AddInitialConditions(kPressure, T0, T0, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}}, 0.0);

    // Phase 1: heat the gas only, briefly, to open a Tf - Ts gap. The steady gap
    // during heating is ~ Q_volumetric / hp; the power is chosen to make it large
    // (order 10 K) so there is a clear disequilibrium to relax afterwards.
    const double heater_power = 1500.0;  // [W]
    const double heat_time = 60.0;       // [s]
    AddStep("Heat", heat_time);
    AddClosedBC(BoundaryLocation::BOTTOM, "Heat");
    AddClosedBC(BoundaryLocation::TOP, "Heat");
    EnergySource heater;
    heater.enabled_ = true;
    heater.power_ = heater_power;
    data.cycle.GetStep("Heat").energy_sources_[layer_id] = heater;

    // Phase 2: heater off, sealed and adiabatic. The two temperatures relax.
    const double relax_time = 1200.0;   // [s]
    AddStep("Relax", relax_time);
    AddClosedBC(BoundaryLocation::BOTTOM, "Relax");
    AddClosedBC(BoundaryLocation::TOP, "Relax");

    ConfigureCycle(1);
    data.settings.reltol = 1e-10;
    data.settings.max_step = 2.0;

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    const auto& log = solver.GetLog();
    const auto& dom = log.GetDomain();

    // --- closed-form coefficients (written from the governing ODEs) ----------
    const double Cv = MixtureMolarCv(data.fluid, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}});
    const double C0 = kPressure / (kR * T0);
    const double et = kVoidFraction;
    const double rho_b = BulkDensity();

    // Volumetric heat capacities [J/(m3-K)]. The gas uses Cv because the sealed
    // rigid vessel converts the Cp in the temperature equation to Cv, exactly as
    // in the closed-bed electric-heating test.
    const double c_gas = et * C0 * Cv;          // gas, per bed volume
    const double c_solid = rho_b * kSolidCp;    // solid, per bed volume

    // Volumetric gas-to-particle coefficient hp = h * av, av = 6 (1 - eb) / dp.
    const double av = 6.0 * (1.0 - kVoidFraction) / kParticleDiameter;
    const double hp = kParticleHtc * av;        // [W/(m3-K)]

    // Relaxation rate of the temperature DIFFERENCE (per bed volume the volume
    // factor cancels): lambda = hp (1/c_gas + 1/c_solid).
    const double lambda = hp * (1.0 / c_gas + 1.0 / c_solid);

    // GUARD: the two thermal masses must be comparable so lambda depends on both.
    const double cap_ratio = c_gas / c_solid;
    ASSERT_GT(cap_ratio, 0.2) << "Gas thermal mass negligible; lambda insensitive to it";
    ASSERT_LT(cap_ratio, 5.0) << "Solid thermal mass negligible; lambda insensitive to it";

    // --- logged trajectories (spatial means over the interior) ---------------
    const auto Tf = log.GetVariable(0, 1, VariableType::TEMPERATURE);        // Relax step
    const auto Ts = log.GetVariable(0, 1, VariableType::SOLID_TEMPERATURE);
    ASSERT_FALSE(Tf.values.empty());
    ASSERT_FALSE(Ts.values.empty());

    std::vector<double> Tf_mean(Tf.values.size()), Ts_mean(Ts.values.size());
    for (size_t n = 0; n < Tf.values.size(); ++n)
      Tf_mean[n] = InteriorMean(Tf.values[n], dom.volume);
    for (size_t n = 0; n < Ts.values.size(); ++n)
      Ts_mean[n] = InteriorMean(Ts.values[n], dom.volume);

    // GUARD 1: heating must have opened a real gas-solid gap.
    const double dT0 = Tf_mean.front() - Ts_mean.front();
    ASSERT_GT(dT0, 5.0)
      << "Gas-solid temperature split at start of relaxation is only " << dT0
      << " K; there is no disequilibrium to relax";

    // The equilibrium temperature is the heat-capacity-weighted mean, and is a
    // constant of the relaxation (no energy leaves the sealed adiabatic bed).
    const double T_eq_expected =
      (c_gas * Tf_mean.front() + c_solid * Ts_mean.front()) / (c_gas + c_solid);

    // GUARD 2: the relaxation must actually resolve.
    const double dT_end = Tf_mean.back() - Ts_mean.back();
    ASSERT_LT(std::abs(dT_end), 0.1 * dT0)
      << "Gas and solid have not equilibrated by the end of the run (gap " << dT_end
      << " K vs initial " << dT0 << " K)";

    // Both must land on the same weighted-mean equilibrium temperature.
    EXPECT_NEAR(Tf_mean.back(), T_eq_expected, 0.02 * (Tf_mean.front() - T_eq_expected))
      << "Gas did not settle at the conserved weighted-mean temperature";
    EXPECT_NEAR(Ts_mean.back(), T_eq_expected, 0.02 * (T_eq_expected - Ts_mean.front()))
      << "Solid did not settle at the conserved weighted-mean temperature";

    // --- recover the relaxation rate from the difference trajectory ----------
    // (Tf - Ts)(t) = dT0 exp(-lambda t): ln(gap) = ln(dT0) - lambda t is a
    // straight line with slope -lambda. Fit over the window 15%..85% of dT0 to
    // stay clear of the tiny initial numerical transient and the noisy tail.
    const auto& t_rel = Tf.times;
    const double t_start = t_rel.front();
    double num = 0.0, den = 0.0;
    int fit_points = 0;
    for (size_t n = 0; n < t_rel.size(); ++n) {
      const double gap = Tf_mean[n] - Ts_mean[n];
      const double frac = gap / dT0;
      if (frac < 0.15 || frac > 0.85) continue;
      const double t = t_rel[n] - t_start;
      const double lhs = std::log(gap) - std::log(dT0);  // = -lambda t
      num += t * lhs;
      den += t * t;
      ++fit_points;
    }
    ASSERT_GT(fit_points, 10) << "Too few points in the relaxation window to fit lambda";
    const double lambda_measured = -num / den;

    EXPECT_NEAR(lambda_measured, lambda, 0.03 * lambda)
      << "Recovered gas-solid relaxation rate " << lambda_measured << " /s disagrees "
      << "with hp (1/c_gas + 1/c_solid) = " << lambda << " /s; the two-temperature "
      << "coupling term is wrong (hp = " << hp << " W/(m3-K), c_gas = " << c_gas
      << ", c_solid = " << c_solid << " J/(m3-K))";

    // --- exact energy conservation over the whole run ------------------------
    // Heater duty in must equal the combined sensible-energy rise of gas + solid.
    const auto Tf_all_heat = log.GetVariable(0, 0, VariableType::TEMPERATURE);
    const auto Ts_all_heat = log.GetVariable(0, 0, VariableType::SOLID_TEMPERATURE);
    const double Tf_initial = InteriorMean(Tf_all_heat.values.front(), dom.volume);
    const double Ts_initial = InteriorMean(Ts_all_heat.values.front(), dom.volume);

    // End of run = end of the relaxation step.
    const double Tf_final = Tf_mean.back();
    const double Ts_final = Ts_mean.back();

    const double dU = BedVolume() * (c_gas * (Tf_final - Tf_initial) +
                                     c_solid * (Ts_final - Ts_initial));
    const double Q_in = heater_power * heat_time;
    EXPECT_NEAR(dU, Q_in, 0.02 * Q_in)
      << "Sealed adiabatic LTNE bed does not conserve energy: heater in = " << Q_in
      << " J, sensible rise (gas + solid) = " << dU << " J";
  }

  // ===========================================================================
  // LTNE bed/wall energy closure
  // ===========================================================================
  //
  // In LTNE the bed touches the wall through BOTH phases. The wall equation
  // (UpdateWallTemperature, NON_LOCAL branch) splits the internal exchange by the
  // bed voidage:
  //
  //     Q_bed->wall = pi di dx hi [ eb (Tf - Tw) + (1 - eb) (Ts - Tw) ]
  //
  // so for the exchange to be a pure transfer the gas must give up exactly the
  // eb share and the solid the (1 - eb) share. Any other split creates or
  // destroys (1 - eb) pi di dx hi (Tf - Ts) per cell whenever the two phases are
  // out of equilibrium.
  //
  // Scenario: sealed, inert (non-adsorbing, so the adsorbed-phase energy is
  // identically zero) LTNE bed, hot (400 K) against a cold (300 K) wall that is
  // insulated from the ambient (ho = 0) and at its ends (ZERO_FLUX). The gas
  // has the smaller thermal mass and cools first; the solid sheds most of its
  // heat through the gas, so Tf != Ts throughout the transient. Bed + wall is a
  // closed system, so
  //
  //     dU_gas + dU_solid + dU_adsorbed + dU_wall = 0      (dU_adsorbed = 0)
  //
  // and the heat the bed gives up must equal what the wall receives.
  //
  // This test originally failed: UpdateFluidTemperature removed the full
  // hi (Tf - Tw) from the gas while the solid had no wall term at all, so the
  // bed and wall disagreed on the exchanged heat (bug E-03).

  TEST_F(LtneRelaxationTest, BedWallExchangeConservesEnergyDuringTransient) {
    data.info.title = "Validation: LTNE bed/wall energy closure";

    BuildBed();

    // Switch the bed<->wall coupling back on; keep the wall insulated outside.
    const double hi = 50.0;  // [W/(m2-K)] bed -> wall
    AddConstantInternalHeatTransferCoefficientModel(hi, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);

    // Weaker gas-to-particle coupling than the relaxation test, so the gas and
    // solid stay well apart while the wall draws heat from both.
    auto php = std::make_unique<ConstantParticleHeatTransferModel>();
    php->heat_transfer_coefficient_ = 1.0;  // [W/(m2-K)]
    data.bed.equations_.particle_heat_transfer_model_ = std::move(php);

    data.bed.ambient_temperature_ = 300.0;
    data.bed.wall_.thickness_ = 0.005;
    data.bed.wall_.density_ = 7900.0;
    data.bed.wall_.heat_capacity_ = 500.0;
    data.bed.wall_.thermal_conductivity_ = 16.2;
    // The default FIXED wall BC would drain energy axially out of the end cells.
    data.bed.equations_.wall_temperature_boundary_condition_ =
      CreateWallTemperatureBoundaryCondition(WallTemperatureBoundaryConditionType::ZERO_FLUX);

    const double T0 = 400.0;   // hot bed (gas and solid)
    const double Tw0 = 300.0;  // cold wall
    AddInitialConditions(kPressure, T0, Tw0, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}}, 0.0);

    AddStep("Cool", 300.0);
    AddClosedBC(BoundaryLocation::BOTTOM, "Cool");
    AddClosedBC(BoundaryLocation::TOP, "Cool");

    ConfigureCycle(1);
    data.settings.reltol = 1e-10;
    data.settings.max_step = 1.0;

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    const auto& log = solver.GetLog();
    const auto& dom = log.GetDomain();

    const auto Tf = log.GetVariable(0, 0, VariableType::TEMPERATURE);
    const auto Ts = log.GetVariable(0, 0, VariableType::SOLID_TEMPERATURE);
    const auto Tw = log.GetVariable(0, 0, VariableType::WALL_TEMPERATURE);
    const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
    const auto hi_log = log.GetVariable(0, 0, VariableType::INTERNAL_HEAT_TRANSFER_COEFFICIENT);
    ASSERT_FALSE(Tf.values.empty());
    ASSERT_EQ(Tf.values.size(), Ts.values.size());
    ASSERT_EQ(Tf.values.size(), Tw.values.size());
    ASSERT_EQ(Tf.values.size(), hi_log.values.size());

    // GUARD 1: there must be a real gas-solid disequilibrium while the bed
    // exchanges heat with the wall, otherwise the (1 - eb) share is untested.
    double max_gap = 0.0;
    for (size_t n = 0; n < Tf.values.size(); ++n) {
      max_gap = std::max(max_gap, std::abs(InteriorMean(Ts.values[n], dom.volume) -
                                           InteriorMean(Tf.values[n], dom.volume)));
    }
    ASSERT_GT(max_gap, 5.0)
      << "Gas-solid split during the transient is only " << max_gap
      << " K; the LTNE wall split is not exercised";

    // --- internal-energy changes (closed, rigid, non-adsorbing bed) ----------
    // Gas uses Cv (the P-V work term in the gas equation converts Cp to Cv at
    // fixed C); the gas inventory C0 is a constant of the motion.
    const double Cv = MixtureMolarCv(data.fluid, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}});
    const auto& Tf0 = Tf.values.front();
    const auto& Ts0 = Ts.values.front();
    const auto& Tf1 = Tf.values.back();
    const auto& Ts1 = Ts.values.back();
    const auto& P0 = P.values.front();
    double dU_gas = 0.0, dU_solid = 0.0;
    for (size_t i = 1; i + 1 < Tf0.size(); ++i) {
      const double C0 = P0[i] / (kR * Tf0[i]);
      dU_gas += dom.volume[i] * dom.et[i] * C0 * Cv * (Tf1[i] - Tf0[i]);
      dU_solid += dom.volume[i] * dom.bulk_density[i] * dom.Cps[i] * (Ts1[i] - Ts0[i]);
    }
    const double dU_adsorbed = 0.0;  // inert packing: no adsorbed phase
    const double dU_bed = dU_gas + dU_solid + dU_adsorbed;
    const double dU_wall = WallInternalEnergyChange(log, 0, 0);

    // Heat crossing the inner wall face, integrated from the textbook flux over
    // pi di dx with the LTNE voidage split.
    std::vector<double> power(Tf.values.size(), 0.0);
    for (size_t n = 0; n < Tf.values.size(); ++n) {
      double q = 0.0;
      for (size_t i = 1; i + 1 < Tf.values[n].size(); ++i) {
        const double eb = kVoidFraction;
        q += InternalWallArea(dom, i) * hi_log.values[n][i] *
             (eb * (Tf.values[n][i] - Tw.values[n][i]) +
              (1.0 - eb) * (Ts.values[n][i] - Tw.values[n][i]));
      }
      power[n] = q;
    }
    const double Q_bed_to_wall = TimeIntegral(Tf.times, power);

    // GUARD 2: the bed must cool appreciably, or the balance is trivially met.
    const double C_bed = BedVolume() * (kVoidFraction * kPressure / (kR * T0) * Cv +
                                        BulkDensity() * kSolidCp);
    ASSERT_LT(dU_bed, -0.5 * C_bed * (T0 - Tw0))
      << "Bed gave up less than half of its excess energy; test proves nothing";

    // (a) The wall gains what crosses its inner face.
    EXPECT_NEAR(dU_wall, Q_bed_to_wall, 0.01 * std::abs(dU_wall))
      << "Wall energy gain " << dU_wall << " J does not match the integrated "
      << "LTNE bed->wall flux " << Q_bed_to_wall << " J";

    // (b) The bed loses that same heat.
    EXPECT_NEAR(-dU_bed, Q_bed_to_wall, 0.01 * std::abs(dU_bed))
      << "Bed energy loss (gas " << dU_gas << " J + solid " << dU_solid
      << " J) does not match the integrated LTNE bed->wall flux "
      << Q_bed_to_wall << " J";

    // (c) Closed, insulated system: total energy is conserved. Both sides are
    // end-state differences (no time quadrature), so this holds far tighter
    // than (a)/(b); the pre-fix code created ~35 % of the wall gain here.
    const double imbalance = dU_bed + dU_wall;
    EXPECT_NEAR(imbalance, 0.0, 1e-3 * std::abs(dU_wall))
      << "LTNE bed + wall does not conserve energy: dU_gas = " << dU_gas
      << " J, dU_solid = " << dU_solid << " J, dU_wall = " << dU_wall
      << " J, imbalance = " << imbalance << " J ("
      << 100.0 * imbalance / std::abs(dU_wall) << "% of the wall gain)";
  }

} // namespace adsim
