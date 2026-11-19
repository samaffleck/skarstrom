/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// test_internal_heat_transfer_kinds.cpp
//
// The bed-to-wall (internal) heat transfer kinds (heat_transfer_model.h):
//   * the shared 1-D lumping and the wall Nusselt correlations against hand
//     calculations (Martin-Nilles at Re = 0 and kappa = 1; the Li-Finlayson
//     power law);
//   * dixon, now written as its wall coefficient, its radial conductivity and
//     the shared lumping, gives bitwise the values of the former inline code;
//   * solver level: in a closed, quiescent bed the logged h_i of
//     martin_nilles is the hand value from the logged gas conductivity, and
//     li_finlayson (no stagnant term) gives h_i = 0;
//   * li_finlayson and martin_nilles run on a bundled dixon deck, give a
//     different temperature history, and survive a save and load.
// Short grids and one cycle: these check that the feature works.

#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "../validation/validation_helpers.h"

#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/deck_edit.h"
#include "adsim/objects/model_kinds.h"
#include "adsim/simulation.h"

namespace adsim {

namespace {

const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";

// The Dixon apparent coefficient as InternalHeatTransferModel_Dixon::Update
// wrote it inline before the lumping became a shared step (kept verbatim).
double FormerInlineDixon(double k_g, double k_s, double eb, double dp, double D_c, double re, double vel,
                         double rho_g, double Cpg) {
  double h_g = (k_g / dp) * (0.0835 * std::pow(re, 0.91));
  double theta_w = 0.0024 * std::pow((D_c / dp), 1.58);
  double h_0 = (k_g / dp) * (2 * eb + (1 - eb) / ((k_g / k_s) * (1.0 / 3.0) + theta_w));
  double h_w = h_0 + h_g;
  double phi = 0.22 * eb * eb;
  double k_e0 = k_g * (eb + ((1 - eb) / ((k_g / k_s) * (2.0 / 3.0) + phi)));
  double d = dp / D_c;
  double Pe_H = 8.65 * (1 + 19.4 * d * d);
  double k_eg = std::abs(vel) * rho_g * Cpg * dp / Pe_H;
  double k_e = k_e0 + k_eg;
  double Bi = h_w * D_c / (2.0 * k_e);
  double h_l_invr = (1.0 / h_w) + (D_c / (6.0 * k_e) * ((Bi + 3.0) / (Bi + 4.0)));
  return 1.0 / h_l_invr;
}

Simulation ShortDeck(const std::string& name) {
  Simulation sim = Simulation::FromFile((kExamples / name).string());
  sim.SetParameter("cycle.number_of_cycles", 1);
  for (auto& layer : sim.Data().bed.layers_.layers_) layer.domain_.n_cells = 20;
  return sim;
}

}  // namespace

// ---- correlations -----------------------------------------------------------

TEST(InternalHeatTransferKinds, LumpingMatchesHandCalculation) {
  // h_w = 100, k_e = 0.5, D = 0.03: Bi = 3, 1/h_i = 0.01 + 0.01 (6/7) = 0.13/7.
  EXPECT_NEAR(LumpedWallHeatTransferCoefficient(100.0, 0.5, 0.03), 7.0 / 0.13, 1e-12);
  // The radial resistance only adds: h_i < h_w, and h_i -> h_w as k_e grows.
  EXPECT_LT(LumpedWallHeatTransferCoefficient(100.0, 0.5, 0.03), 100.0);
  EXPECT_NEAR(LumpedWallHeatTransferCoefficient(100.0, 1e9, 0.03), 100.0, 1e-5);
  // No wall coefficient, no heat to the wall.
  EXPECT_EQ(LumpedWallHeatTransferCoefficient(0.0, 0.5, 0.03), 0.0);
}

TEST(InternalHeatTransferKinds, LiFinlaysonPowerLaw) {
  // 0.17 Re^0.79 (spheres), 0.16 Re^0.93 (cylinders) at Re = 100.
  EXPECT_NEAR(LiFinlaysonWallNusselt(100.0, false), 6.463219737449542, 1e-12);
  EXPECT_NEAR(LiFinlaysonWallNusselt(100.0, true), 11.590975361199844, 1e-12);
  // Doubling Re multiplies h_w by 2^0.79 (spheres), 2^0.93 (cylinders).
  for (double re : {25.0, 300.0, 4000.0}) {
    EXPECT_NEAR(LiFinlaysonWallNusselt(2 * re, false) / LiFinlaysonWallNusselt(re, false), std::pow(2.0, 0.79), 1e-12);
    EXPECT_NEAR(LiFinlaysonWallNusselt(2 * re, true) / LiFinlaysonWallNusselt(re, true), std::pow(2.0, 0.93), 1e-12);
  }
  EXPECT_EQ(LiFinlaysonWallNusselt(0.0, false), 0.0);
}

TEST(InternalHeatTransferKinds, MartinNillesHandCalculation) {
  // Re -> 0 with kappa = 1 (k0 = k_g): Nu_w = 1.3 + 5 d_p/D.
  const double k0_unit = ZehnerSchlunderConductivityRatio(0.4, 1.0);
  EXPECT_NEAR(k0_unit, 1.0, 1e-12);
  for (double d : {0.02, 0.1, 0.25})
    EXPECT_NEAR(MartinNillesWallNusselt(0.0, 0.7, k0_unit, d), 1.3 + 5.0 * d, 1e-12);
  // Re = 100, Pr = 0.7, k0/k_g = 3.643, d_p/D = 0.05:
  // 1.55 * 3.643 + 0.19 * 100^0.75 * 0.7^(1/3).
  EXPECT_NEAR(MartinNillesWallNusselt(100.0, 0.7, 3.643, 0.05), 10.981468079260992, 1e-12);
  // The radial conductivity: k_g [k0/k_g + Re Pr / K_r], K_r = 7 [2 - (1 - 2 d)^2] = 8.33 at d = 0.05.
  EXPECT_NEAR(ZehnerSchlunderRadialConductivity(0.025, 3.643, 100.0, 0.7, 0.05), 0.025 * (3.643 + 70.0 / 8.33),
              1e-14);
  EXPECT_NEAR(ZehnerSchlunderRadialConductivity(0.025, 3.643, 0.0, 0.7, 0.05), 0.025 * 3.643, 1e-15);
}

// dixon is the shared lumping of DixonWallCoefficient and
// DixonRadialConductivity; the composition is bitwise the former inline code.
TEST(InternalHeatTransferKinds, DixonSplitIsBitwiseTheFormerFormula) {
  for (double k_g : {0.018, 0.026, 0.15})
    for (double k_s : {0.1, 0.6, 20.0})
      for (double eb : {0.33, 0.4, 0.48})
        for (double dp : {5e-4, 2e-3})
          for (double D_c : {0.02, 0.1})
            for (double re : {0.0, 3.7, 150.0})
              for (double vel : {0.0, -0.05, 0.3}) {
                const double rho = 1.2, cp = 1040.0;
                const double h_w = DixonWallCoefficient(k_g, k_s, eb, dp, D_c, re);
                const double k_e = DixonRadialConductivity(k_g, k_s, eb, dp, D_c, vel, rho, cp);
                const double now = LumpedWallHeatTransferCoefficient(h_w, k_e, D_c);
                const double before = FormerInlineDixon(k_g, k_s, eb, dp, D_c, re, vel, rho, cp);
                ASSERT_EQ(now, before) << k_g << ' ' << k_s << ' ' << eb << ' ' << dp << ' ' << D_c << ' ' << re
                                       << ' ' << vel;
              }
}

// ---- solver level -----------------------------------------------------------

// A closed, quiescent, non-adsorbing bed (Re = 0 everywhere): martin_nilles
// gives h_w = (1.3 + 5 d_p/D) k0 / d_p and k_e = k0, so the logged h_i is the
// hand value from the logged gas conductivity; li_finlayson has no stagnant
// term and gives h_i = 0.
class ClosedBedWallCoefficient : public UserDataTest {
public:
  static constexpr double kDiameter = 0.05;   // [m]
  static constexpr double kDp = 0.002;        // [m]
  static constexpr double kEb = 0.4;
  static constexpr double kSolidConductivity = 0.5;  // [W/(m K)]

  // Returns the logged h_i and k_g in the interior cells at the end of the step.
  void Run(const std::string& kind, std::vector<double>& hi, std::vector<double>& kg) {
    AddComponents({"Nitrogen", "Oxygen"});
    auto [layer_id, ads_id] = AddLayer("Packing");
    auto& layer = GetLayer(layer_id);
    layer.domain_.length = 0.5;
    layer.domain_.n_cells = 10;
    layer.eb_ = kEb;
    layer.particle_diameter_ = kDp;
    std::get<CylinderBed>(data.bed.shape_).diameter = kDiameter;
    auto& ads = GetAdsorbent(ads_id);
    ads.solid_density_ = 1000.0;
    ads.heat_capacity_ = 900.0;
    ads.thermal_conductivity_ = kSolidConductivity;
    ads.ep_ = 0.0;
    AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);
    build::SetInternalHeatTransfer(data, layer_id, kind);
    AddInitialConditions(2e5, 310.0, 300.0, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}}, 0.0);
    AddStep("Hold", 20.0);
    AddClosedBC(BoundaryLocation::BOTTOM, "Hold");
    AddClosedBC(BoundaryLocation::TOP, "Hold");
    ConfigureCycle(1);

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
    const auto& log = solver.GetLog();
    const auto h = log.GetVariable(0, 0, VariableType::INTERNAL_HEAT_TRANSFER_COEFFICIENT);
    const auto k = log.GetVariable(0, 0, VariableType::THERMAL_CONDUCTIVITY);
    ASSERT_FALSE(h.values.empty());
    hi.assign(h.values.back().begin() + 1, h.values.back().end() - 1);
    kg.assign(k.values.back().begin() + 1, k.values.back().end() - 1);
  }
};

TEST_F(ClosedBedWallCoefficient, MartinNillesAtRestIsTheStagnantLimit) {
  std::vector<double> hi, kg;
  Run("martin_nilles", hi, kg);
  ASSERT_EQ(hi.size(), kg.size());
  ASSERT_FALSE(hi.empty());
  for (size_t i = 0; i < hi.size(); ++i) {
    const double k0 = kg[i] * ZehnerSchlunderConductivityRatio(kEb, kSolidConductivity / kg[i]);
    const double h_w = (1.3 + 5.0 * kDp / kDiameter) * k0 / kDp;
    const double Bi = h_w * kDiameter / (2.0 * k0);
    const double expected = 1.0 / (1.0 / h_w + kDiameter / (6.0 * k0) * (Bi + 3.0) / (Bi + 4.0));
    EXPECT_NEAR(hi[i], expected, 1e-9 * expected) << "cell " << i;
    EXPECT_GT(hi[i], 0.0);
  }
}

TEST_F(ClosedBedWallCoefficient, LiFinlaysonAtRestIsZero) {
  std::vector<double> hi, kg;
  Run("li_finlayson", hi, kg);
  ASSERT_FALSE(hi.empty());
  for (size_t i = 0; i < hi.size(); ++i) EXPECT_EQ(hi[i], 0.0) << "cell " << i;
}

// The new kinds run on a bundled dixon deck, change the temperature history,
// and keep their kind through a save and load.
class NewWallKindsOnDixonDeck : public ::testing::TestWithParam<std::string> {};

TEST_P(NewWallKindsOnDixonDeck, RunsDiffersAndRoundTrips) {
  Simulation base = ShortDeck("Streb2021_A1.json");
  ASSERT_NE(dynamic_cast<const InternalHeatTransferModel_Dixon*>(
                base.Data().bed.layers_.layers_.front().equations_.internal_heat_transfer_model_.get()),
            nullptr);
  ASSERT_TRUE(base.Run().Ok());

  Simulation sim = ShortDeck("Streb2021_A1.json");
  for (auto& layer : sim.Data().bed.layers_.layers_) build::SetInternalHeatTransfer(sim.Data(), layer.id_, GetParam());
  const RunResult r = sim.Run();
  ASSERT_TRUE(r.Ok()) << r.message;
  EXPECT_TRUE(r.healthy);
  const auto t0 = base.Series("temperature", "", SeriesLocation::Top).values;
  const auto t1 = sim.Series("temperature", "", SeriesLocation::Top).values;
  ASSERT_FALSE(t1.empty());
  EXPECT_NE(t0, t1);
  const auto hi = sim.Series("internal_heat_transfer_coefficient", "", SeriesLocation::Top).values;
  ASSERT_FALSE(hi.empty());
  for (double h : hi) {
    EXPECT_TRUE(std::isfinite(h));
    EXPECT_GE(h, 0.0);
  }

  Simulation copy = Simulation::FromJson(sim.ToJson());
  const auto* m = copy.Data().bed.layers_.layers_.front().equations_.internal_heat_transfer_model_.get();
  ASSERT_NE(m, nullptr);
  EXPECT_EQ(KindDescriptors(ModelFamily::InternalHeatTransfer)[static_cast<size_t>(
                const_cast<IInternalHeatTransferModel*>(m)->GetModelType())]
                .tag,
            GetParam());
}

INSTANTIATE_TEST_SUITE_P(InternalHeatTransferKinds, NewWallKindsOnDixonDeck,
                         ::testing::Values("li_finlayson", "martin_nilles"));

TEST(InternalHeatTransferKinds, OfferedOnPackedLayersOnly) {
  EXPECT_EQ(build::InternalHeatTransferKinds(),
            (std::vector<std::string>{"constant", "dixon", "leva", "li_finlayson", "martin_nilles"}));
  for (const char* tag : {"li_finlayson", "martin_nilles"}) {
    EXPECT_TRUE(KindAcceptsLayer(ModelFamily::InternalHeatTransfer, tag, "packed")) << tag;
    EXPECT_FALSE(KindAcceptsLayer(ModelFamily::InternalHeatTransfer, tag, "monolith")) << tag;
  }
}

}  // namespace adsim
