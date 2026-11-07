/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Heat in structured layers (monolith, laminate), each against a result
// written out here independently of core/src:
//
//  1. Column wall to solid (LTNE): in a structured layer the whole bed-to-wall
//     flux goes into the solid (the gas reaches the wall only through the
//     channel walls). In a closed bed with no gas-solid exchange the gas
//     temperature then stays put exactly, and the solid and the column wall
//     relax as a two-body system,
//         d(Ts - Tw)/dt = -H (1/c_s + 1/c_w) (Ts - Tw),  H = 4 h_i / d_i,
//     with c_s = rho_b Cp_s and c_w = rho_w Cp_w (d_o^2 - d_i^2) / d_i^2 per
//     bed volume. A packed bed, where the gas takes the eb share, is the
//     contrast.
//  2. Gas to channel wall at a known Nusselt number (LTNE): the film
//     coefficient of fully developed laminar flow between parallel plates,
//     Nu_T = 7.541 on D_h = 2 x gap (Shah and London 1978), gives
//     h_p = Nu k_g / D_h a_v, and the gas-solid temperature gap of a closed
//     bed decays as exp(-lambda t), lambda = h_p (1/c_g + 1/c_s), c_g =
//     OFA C Cv (a sealed rigid bed turns the gas's Cp into Cv).
//  3. Steady axial conduction through the solid (Fourier): a heater in the
//     bottom layer and an equal cooler in the top one drive a steady flux Q
//     through the middle layer's solid, whose temperature is then linear with
//     slope -Q / (A lambda_s), lambda_s = (1 - OFA) k_s.
//  4. Energy closure with solid conduction on: heater duty = rise in the
//     internal energy of gas, solid and column wall (sealed, adiabatic
//     outside). The gas term is sum V et Cv dP / R, exact for an ideal gas
//     whatever moves between cells.
//  5. The channel heat-dispersion model's conductivity: OFA k_g under LTNE,
//     OFA k_g + (1 - OFA) k_s under LTE.
//  6. Packed beds take no solid conduction: an LTNE packed bed with a
//     temperature gradient gives the same answer whatever the solid
//     conductivity.

#include "validation_helpers.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace adsim {

  using namespace adsim::validation;

  class StructuredHeat : public UserDataTest {
  public:
    static constexpr double kDiameter = 0.05;  // [m]
    static constexpr double kT0 = 300.0;       // [K]

    // A laminate layer with an inert solid (no isotherm, ep = 0, so et = OFA).
    LayerID AddLaminate(const std::string& name, double length, int cells, double sheet, double gap,
                        double rho_s, double cp_s, double k_s) {
      auto [lid, aid] = AddLayer(name);
      auto& layer = GetLayer(lid);
      layer.kind_ = LayerKind::Laminate;
      layer.structured_.sheet_thickness_ = sheet;
      layer.structured_.gap_ = gap;
      layer.domain_.length = length;
      layer.domain_.n_cells = cells;
      layer.equations_.heat_dispersion_model_ = std::make_unique<HeatDispersionModelChannel>();
      auto& ads = GetAdsorbent(aid);
      ads.ep_ = 0.0;
      ads.solid_density_ = rho_s;
      ads.heat_capacity_ = cp_s;
      ads.thermal_conductivity_ = k_s;
      AddConstantInternalHeatTransferCoefficientModel(0.0, lid);
      AddConstantExternalHeatTransferCoefficientModel(0.0, lid);
      return lid;
    }

    void Ltne(std::unique_ptr<IParticleHeatTransferModel> model) {
      data.bed.equations_.thermal_equilibrium_mode_ = ThermalEquilibriumMode::NON_LOCAL;
      data.bed.equations_.particle_heat_transfer_model_ = std::move(model);
    }

    static std::unique_ptr<IParticleHeatTransferModel> ConstantHp(double h) {
      auto m = std::make_unique<ConstantParticleHeatTransferModel>();
      m->heat_transfer_coefficient_ = h;
      return m;
    }

    void Sealed(const std::string& step, double time) {
      AddStep(step, time);
      AddClosedBC(BoundaryLocation::BOTTOM, step);
      AddClosedBC(BoundaryLocation::TOP, step);
    }

    void Heater(const std::string& step, LayerID layer, double power) {
      EnergySource s;
      s.enabled_ = true;
      s.power_ = power;
      data.cycle.GetStep(step).energy_sources_[layer] = s;
    }

    void Run(CoupledSolver& solver) {
      std::string error;
      ASSERT_TRUE(data.EvaluateAllExpressions(error) && data.PerformChecks(error)) << error;
      ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
    }

    double Cv() const { return MixtureMolarCv(data.fluid, {{"Nitrogen", 1.0}}); }
  };

  // 1. The whole column-wall flux goes into the solid of a structured layer.
  TEST_F(StructuredHeat, ColumnWallHeatsOnlyTheSolidOfAStructuredLayer) {
    data.info.title = "Validation: structured layer, column wall to solid";
    AddComponents({"Nitrogen"});
    std::get<CylinderBed>(data.bed.shape_).diameter = kDiameter;
    const double hi = 100.0, rho_s = 2000.0, cp_s = 2000.0;
    const LayerID lid = AddLaminate("Sheets", 0.1, 10, 2e-4, 2e-4, rho_s, cp_s, 1.0);
    AddConstantInternalHeatTransferCoefficientModel(hi, lid);
    data.bed.wall_.thickness_ = 0.005;
    Ltne(ConstantHp(0.0));  // gas and solid exchange nothing directly
    const double Tw0 = 330.0;
    AddInitialConditions(1e5, kT0, Tw0, {{"Nitrogen", 1.0}});
    Sealed("Relax", 600.0);
    ConfigureCycle(1);
    data.settings.reltol = 1e-10;
    data.settings.max_step = 1.0;

    CoupledSolver solver{data};
    Run(solver);
    const auto& log = solver.GetLog();
    const auto& dom = log.GetDomain();
    const auto Tf = log.GetVariable(0, 0, VariableType::TEMPERATURE);
    const auto Ts = log.GetVariable(0, 0, VariableType::SOLID_TEMPERATURE);
    const auto Tw = log.GetVariable(0, 0, VariableType::WALL_TEMPERATURE);

    const auto& layer = data.bed.layers_.layers_.front();
    const double ofa = layer.GetOpenFraction();
    const double di = kDiameter, d_o = di + 2.0 * data.bed.wall_.thickness_;
    const double c_s = (1.0 - ofa) * rho_s * cp_s;
    const double c_w = data.bed.wall_.density_ * data.bed.wall_.heat_capacity_ * (d_o * d_o - di * di) / (di * di);
    const double H = 4.0 * hi / di;
    const double lambda = H * (1.0 / c_s + 1.0 / c_w);
    ASSERT_GT(c_w / c_s, 0.2);
    ASSERT_LT(c_w / c_s, 5.0);

    double worst_gas = 0.0, worst_gap = 0.0;
    const double gap0 = (InteriorMean(Ts.values.front(), dom.volume) - InteriorMean(Tw.values.front(), dom.volume));
    ASSERT_NEAR(gap0, kT0 - Tw0, 1e-9);
    for (size_t n = 0; n < Ts.values.size(); ++n) {
      const double t = Ts.times[n] - Ts.times.front();
      for (size_t i = 1; i + 1 < Tf.values[n].size(); ++i) worst_gas = std::max(worst_gas, std::abs(Tf.values[n][i] - kT0));
      const double gap = InteriorMean(Ts.values[n], dom.volume) - InteriorMean(Tw.values[n], dom.volume);
      worst_gap = std::max(worst_gap, std::abs(gap - gap0 * std::exp(-lambda * t)) / std::abs(gap0));
    }
    const double t_end = Ts.times.back() - Ts.times.front();
    std::cout << "[wall to solid] lambda " << lambda << " 1/s over " << lambda * t_end << " time constants; "
              << "worst gas drift " << worst_gas << " K, worst gap error " << worst_gap << " of the initial gap\n";
    ASSERT_GT(lambda * t_end, 3.0) << "the relaxation must run most of its course";
    EXPECT_LT(worst_gas, 1e-8) << "the gas must take no share of the wall flux in a structured layer";
    EXPECT_LT(worst_gap, 1e-4);
    // The common temperature from conservation.
    const double T_eq = (c_s * kT0 + c_w * Tw0) / (c_s + c_w);
    EXPECT_NEAR(InteriorMean(Ts.values.back(), dom.volume), T_eq, 1e-3 * (Tw0 - kT0) + std::abs(gap0) * std::exp(-lambda * t_end));
  }

  // 1b. The contrast: a packed bed's gas takes the eb share of the wall flux.
  TEST_F(StructuredHeat, PackedBedGasStillTakesItsShareOfTheWallFlux) {
    AddComponents({"Nitrogen"});
    std::get<CylinderBed>(data.bed.shape_).diameter = kDiameter;
    const LayerID lid = AddLaminate("Beads", 0.1, 10, 2e-4, 2e-4, 2000.0, 2000.0, 1.0);
    auto& layer = GetLayer(lid);
    layer.kind_ = LayerKind::Packed;
    layer.equations_.heat_dispersion_model_ = std::make_unique<HeatDispersionModelConstant>();
    AddConstantInternalHeatTransferCoefficientModel(100.0, lid);
    data.bed.wall_.thickness_ = 0.005;
    Ltne(ConstantHp(0.0));
    AddInitialConditions(1e5, kT0, 330.0, {{"Nitrogen", 1.0}});
    Sealed("Relax", 60.0);
    ConfigureCycle(1);
    CoupledSolver solver{data};
    Run(solver);
    const auto Tf = solver.GetLog().GetVariable(0, 0, VariableType::TEMPERATURE);
    EXPECT_GT(InteriorMean(Tf.values.back(), solver.GetLog().GetDomain().volume) - kT0, 1.0);
  }

  // 2. Gas-solid relaxation at the rate of Shah and London's Nu for plates.
  TEST_F(StructuredHeat, ChannelWallRelaxesAtTheShahAndLondonNusseltRate) {
    data.info.title = "Validation: structured layer, gas to channel wall";
    AddComponents({"Nitrogen"});
    std::get<CylinderBed>(data.bed.shape_).diameter = kDiameter;
    const double rho_s = 2000.0, cp_s = 20.0, gap = 2e-3, sheet = 2e-3, P0 = 30e5;
    const LayerID lid = AddLaminate("Sheets", 0.1, 10, sheet, gap, rho_s, cp_s, 1.0);
    Ltne(CreateParticleHeatTransferModel(ParticleHeatTransferType::FILM));
    AddInitialConditions(P0, kT0, kT0, {{"Nitrogen", 1.0}});
    Sealed("Heat", 2.0);
    Heater("Heat", lid, 25.0);  // into the gas
    Sealed("Relax", 4.0);
    ConfigureCycle(1);
    data.settings.reltol = 1e-10;
    data.settings.max_step = 0.005;

    CoupledSolver solver{data};
    Run(solver);
    const auto& log = solver.GetLog();
    const auto& dom = log.GetDomain();
    const auto Tf = log.GetVariable(0, 1, VariableType::TEMPERATURE);
    const auto Ts = log.GetVariable(0, 1, VariableType::SOLID_TEMPERATURE);
    const auto kg = log.GetVariable(0, 1, VariableType::THERMAL_CONDUCTIVITY);
    const auto hp = log.GetVariable(0, 1, VariableType::PARTICLE_HEAT_TRANSFER_COEFFICIENT);

    const auto& layer = data.bed.layers_.layers_.front();
    const double ofa = layer.GetOpenFraction();
    const double dh = 2.0 * gap, av = 2.0 / (gap + sheet), nu = 7.541;
    ASSERT_DOUBLE_EQ(layer.GetHydraulicDiameter(), dh);
    // The mean temperature hardly moves while the gap relaxes, so k_g is
    // taken at the end; the film coefficient must be Nu k_g / D_h a_v there.
    const size_t mid = 5;
    const double k_g = kg.values.back()[mid];
    const double h_p = nu * k_g / dh * av;
    EXPECT_NEAR(hp.values.back()[mid], h_p, 1e-10 * h_p);

    const double C0 = P0 / (kR * kT0);
    const double c_g = ofa * C0 * Cv();
    const double c_s = (1.0 - ofa) * rho_s * cp_s;
    ASSERT_GT(c_g / c_s, 0.2);
    ASSERT_LT(c_g / c_s, 5.0);
    const double lambda = h_p * (1.0 / c_g + 1.0 / c_s);

    const double gap0 = InteriorMean(Tf.values.front(), dom.volume) - InteriorMean(Ts.values.front(), dom.volume);
    ASSERT_GT(gap0, 2.0) << "the heater must open a real gas-solid gap";
    double worst = 0.0;
    for (size_t n = 0; n < Tf.values.size(); ++n) {
      const double t = Tf.times[n] - Tf.times.front();
      if (lambda * t > 6.0) break;  // the gap is then below the integrator's noise
      const double g = InteriorMean(Tf.values[n], dom.volume) - InteriorMean(Ts.values[n], dom.volume);
      worst = std::max(worst, std::abs(g - gap0 * std::exp(-lambda * t)) / gap0);
    }
    std::cout << "[gas to channel wall] h_p " << h_p << " W/(m3 K), lambda " << lambda << " 1/s, initial gap "
              << gap0 << " K, worst error " << worst << " of the initial gap\n";
    // k_g moves with the mean temperature (~1e-3 over the relaxation).
    EXPECT_LT(worst, 3e-3);
  }

  // 3. Steady conduction through the middle layer's solid.
  TEST_F(StructuredHeat, SteadySolidConductionIsLinearWithFouriersFlux) {
    data.info.title = "Validation: steady axial conduction through the solid";
    AddComponents({"Nitrogen"});
    std::get<CylinderBed>(data.bed.shape_).diameter = kDiameter;
    const double k_s = 20.0, length = 0.05;
    const int cells = 10;
    const LayerID a = AddLaminate("Heated", length, cells, 2e-3, 2e-3, 2000.0, 100.0, k_s);
    const LayerID b = AddLaminate("Conductor", length, cells, 2e-3, 2e-3, 2000.0, 100.0, k_s);
    const LayerID c = AddLaminate("Cooled", length, cells, 2e-3, 2e-3, 2000.0, 100.0, k_s);
    // Isolate the solid path: no gas conduction anywhere.
    for (LayerID id : {a, b, c}) GetLayer(id).equations_.heat_dispersion_model_ = std::make_unique<HeatDispersionModelConstant>();
    Ltne(CreateParticleHeatTransferModel(ParticleHeatTransferType::FILM));
    AddInitialConditions(1e5, kT0, kT0, {{"Nitrogen", 1.0}});
    const double Q = 4.0;  // [W]
    Sealed("Conduct", 400.0);
    Heater("Conduct", a, Q);
    Heater("Conduct", c, -Q);
    ConfigureCycle(1);
    data.settings.reltol = 1e-10;
    data.settings.max_step = 2.0;

    CoupledSolver solver{data};
    Run(solver);
    const auto& log = solver.GetLog();
    const auto Ts = log.GetVariable(0, 0, VariableType::SOLID_TEMPERATURE).values.back();
    const auto& layer = data.bed.layers_.Get(b);
    const double ofa = layer.GetOpenFraction();
    const double lambda_s = (1.0 - ofa) * k_s;
    const double A = M_PI * kDiameter * kDiameter / 4.0;
    const double dx = length / cells;
    const double slope_expected = -Q / (A * lambda_s);

    // Cell centres of the middle layer; skip the cell either side of its
    // faces, which sit next to the heated and cooled layers.
    const int first = layer.domain_.start_index + 1, last = layer.domain_.end_index - 1;
    const double slope = (Ts[last] - Ts[first]) / ((last - first) * dx);
    double worst_line = 0.0;
    for (int i = first; i <= last; ++i) {
      const double line = Ts[first] + slope * (i - first) * dx;
      worst_line = std::max(worst_line, std::abs(Ts[i] - line));
    }
    const double dT = std::abs(slope_expected) * length;
    std::cout << "[steady conduction] slope " << slope << " K/m, Fourier " << slope_expected << " K/m; worst departure"
              << " from a line " << worst_line << " K over " << dT << " K\n";
    ASSERT_GT(dT, 5.0) << "the flux must drive a real gradient";
    EXPECT_NEAR(slope, slope_expected, 1e-5 * std::abs(slope_expected));
    EXPECT_LT(worst_line, 1e-6 * dT);
    // The faces either side of the middle layer carry the same flux.
    EXPECT_NEAR((Ts[layer.domain_.start_index] - Ts[layer.domain_.start_index - 1]) / dx, slope_expected,
                1e-5 * std::abs(slope_expected));
    EXPECT_NEAR((Ts[layer.domain_.end_index + 1] - Ts[layer.domain_.end_index]) / dx, slope_expected,
                1e-5 * std::abs(slope_expected));
  }

  // 4. Energy closes with solid conduction, gas conduction and the wall flux.
  TEST_F(StructuredHeat, EnergyClosesWithSolidConduction) {
    data.info.title = "Validation: energy closure, structured layers";
    AddComponents({"Nitrogen"});
    std::get<CylinderBed>(data.bed.shape_).diameter = kDiameter;
    const double rho_s = 2000.0, cp_s = 300.0;
    const LayerID a = AddLaminate("Heated", 0.05, 10, 2e-3, 2e-3, rho_s, cp_s, 20.0);
    const LayerID b = AddLaminate("Walled", 0.05, 10, 2e-3, 2e-3, rho_s, cp_s, 20.0);
    AddConstantInternalHeatTransferCoefficientModel(50.0, b);  // heat leaves only through b's solid
    data.bed.wall_.thickness_ = 0.003;
    Ltne(CreateParticleHeatTransferModel(ParticleHeatTransferType::FILM));
    const double P0 = 1e5;
    AddInitialConditions(P0, kT0, kT0, {{"Nitrogen", 1.0}});
    Sealed("Heat", 200.0);
    Heater("Heat", a, 5.0);
    Sealed("Rest", 200.0);
    ConfigureCycle(1);
    data.settings.reltol = 1e-10;
    data.settings.max_step = 1.0;

    CoupledSolver solver{data};
    Run(solver);
    const auto& log = solver.GetLog();
    const auto& dom = log.GetDomain();
    const double cv = Cv();
    double worst = 0.0;
    for (int step : {0, 1}) {
      const auto P = log.GetVariable(0, step, VariableType::PRESSURE);
      const auto Ts = log.GetVariable(0, step, VariableType::SOLID_TEMPERATURE);
      double dU_gas = 0.0, dU_solid = 0.0;
      for (size_t i = 1; i + 1 < P.values.front().size(); ++i) {
        dU_gas += dom.volume[i] * dom.et[i] * cv * (P.values.back()[i] - P.values.front()[i]) / kR;
        dU_solid += dom.volume[i] * dom.bulk_density[i] * dom.Cps[i] * (Ts.values.back()[i] - Ts.values.front()[i]);
      }
      const double dU_wall = WallInternalEnergyChange(log, 0, step);
      const double Q_in = EnergySourceHeat(log, 0, step);
      const double imbalance = dU_gas + dU_solid + dU_wall - Q_in;
      const double scale = std::max({std::abs(dU_solid), std::abs(dU_wall), std::abs(Q_in)});
      std::cout << "[energy closure] step " << step << ": gas " << dU_gas << " J, solid " << dU_solid << " J, wall "
                << dU_wall << " J, heater " << Q_in << " J, imbalance " << imbalance << " J\n";
      if (step == 0) ASSERT_GT(Q_in, 900.0);
      ASSERT_GT(std::abs(dU_wall), 0.05 * scale) << "the wall must take a real share through b's solid";
      worst = std::max(worst, std::abs(imbalance) / scale);
    }
    EXPECT_LT(worst, 1e-4);
  }

  // 5. The channel model's axial conductivity under LTNE and LTE.
  TEST_F(StructuredHeat, ChannelConductivityIsGasAloneUnderLtneAndLumpedUnderLte) {
    for (bool ltne : {true, false}) {
      UserData& d = data;
      d = UserData{};
      SetUp();
      AddComponents({"Nitrogen"});
      std::get<CylinderBed>(d.bed.shape_).diameter = kDiameter;
      const double k_s = 3.0;
      const LayerID lid = AddLaminate("Sheets", 0.05, 6, 3e-4, 2e-4, 2000.0, 900.0, k_s);
      AddConstantInternalHeatTransferCoefficientModel(0.0, lid);
      if (ltne) Ltne(ConstantHp(10.0));
      AddInitialConditions(1e5, kT0, kT0, {{"Nitrogen", 1.0}});
      Sealed("Hold", 1.0);
      ConfigureCycle(1);
      CoupledSolver solver{d};
      Run(solver);
      const auto& log = solver.GetLog();
      const auto lam = log.GetVariable(0, 0, VariableType::HEAT_DISPERSION).values.back();
      const auto kg = log.GetVariable(0, 0, VariableType::THERMAL_CONDUCTIVITY).values.back();
      const double ofa = d.bed.layers_.layers_.front().GetOpenFraction();
      EXPECT_NEAR(ofa, 2e-4 / 5e-4, 1e-15);
      const double expected = ofa * kg[3] + (ltne ? 0.0 : (1.0 - ofa) * k_s);
      EXPECT_NEAR(lam[3], expected, 1e-12 * expected) << (ltne ? "LTNE" : "LTE");
      EXPECT_DOUBLE_EQ(log.GetDomain().solid_conductivity[3], (1.0 - ofa) * k_s);
    }
  }

  // 6. A packed bed has no solid conduction path: with a gradient along it,
  //    its LTNE result does not depend on the solid conductivity.
  TEST_F(StructuredHeat, PackedBedTakesNoSolidConduction) {
    std::vector<double> reference;
    for (double k_s : {0.0, 50.0}) {
      UserData& d = data;
      d = UserData{};
      SetUp();
      AddComponents({"Nitrogen"});
      std::get<CylinderBed>(d.bed.shape_).diameter = kDiameter;
      const LayerID a = AddLaminate("Heated", 0.05, 8, 2e-3, 2e-3, 2000.0, 100.0, k_s);
      const LayerID b = AddLaminate("Plain", 0.05, 8, 2e-3, 2e-3, 2000.0, 100.0, k_s);
      for (LayerID id : {a, b}) {
        GetLayer(id).kind_ = LayerKind::Packed;
        GetLayer(id).equations_.heat_dispersion_model_ = std::make_unique<HeatDispersionModelConstant>();
      }
      Ltne(ConstantHp(5.0));
      AddInitialConditions(1e5, kT0, kT0, {{"Nitrogen", 1.0}});
      Sealed("Heat", 30.0);
      Heater("Heat", a, 5.0);
      ConfigureCycle(1);
      CoupledSolver solver{d};
      Run(solver);
      EXPECT_EQ(solver.GetLog().GetDomain().solid_conductivity[3], 0.0);
      const auto Ts = solver.GetLog().GetVariable(0, 0, VariableType::SOLID_TEMPERATURE).values.back();
      if (reference.empty()) {
        reference = Ts;
        ASSERT_GT(Ts[3] - Ts[12], 1.0) << "there must be a gradient for conduction to act on";
      } else {
        ASSERT_EQ(Ts.size(), reference.size());
        for (size_t i = 0; i < Ts.size(); ++i) EXPECT_EQ(Ts[i], reference[i]) << "cell " << i;
      }
    }
  }

}  // namespace adsim
