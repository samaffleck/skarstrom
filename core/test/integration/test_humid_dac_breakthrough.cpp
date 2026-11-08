// Illustrative humid direct-air-capture breakthrough: 400 ppm CO2 in humid N2
// through a column of the amine sorbent APDES-NFC, non-isothermal, with water
// adsorbing on a GAB isotherm and CO2 on the humid Toth, which reads the
// adsorbed water loading. Writes
// examples/decks/Illustrative_Humid_DAC_Breakthrough.json.
//
// ILLUSTRATIVE ONLY. The sorbent is the library entry "APDES-NFC"
// (data/adsorbents/apdes_nfc_stampi2020.json, loaded unchanged): GAB water,
// humid Toth CO2 (Stampi-Bombelli, Storione and Mazzotti, Adsorption 26 (2020)
// 1183, corrected 2022), LDF k_CO2 = 2e-4 and k_H2O = 2e-3 1/s, heats of
// adsorption 57 and 49 kJ/mol. The feed is the paper's: 400 ppm CO2,
// y_H2O = 0.0115 (49 % RH at 20 degC), balance N2, 20 degC, 1 bar, 50 cm3/s
// through its 8 cm column (50.3 cm2), a superficial velocity of 0.01 m/s.
// Bed-to-wall and wall-to-ambient coefficients are the paper's 3 and
// 26 W/(m2 K). The bed length (0.2 m, 25 cells), voidage, wall thickness and
// the step time are ours: the paper does not run a breakthrough on this
// column. With the slow LDF (k_CO2 = 2e-4 1/s) the rate groups k L / u times
// the capacity ratio are about 10 for both components at this length, so the
// fronts are sharp.
//
// Measured (Linux, release): t50 water 2.67e3 s, t50 CO2 8.98e4 s, t95 CO2
// 1.03e5 s (dry feed: t50 CO2 3.45e4 s); mid-bed temperature peak 300.1 K at
// 681 s; end of step, last cell: q_H2O 2.819 mol/kg (GAB 2.819), q_CO2
// 2.066 mol/kg (humid Toth 2.070; dry run 0.822).
//
// Checks (fast tier, two runs of one step on 25 cells):
//   * the run completes, is healthy, and the step and KPI mole balances close
//     for N2, CO2 and water;
//   * water breaks through first; at the end of the step the water loading at
//     the outlet is the GAB value at the local relative humidity and the CO2
//     loading there is the humid Toth value at the adsorbed water loading;
//   * the CO2 loading is well above the dry value, and the same deck with a
//     dry feed reaches the dry value;
//   * the bed temperature peaks (water sorption heat) before the CO2 front
//     reaches 50 % at the outlet.
// The isotherm spot checks at the feed state are in the last test.

#include "../validation/validation_helpers.h"

#include <algorithm>
#include <filesystem>

#include "adsim/equations/isotherm_kinds.h"
#include "adsim/simulation.h"

namespace adsim {

using namespace adsim::validation;

namespace {

const std::string kN2 = "Nitrogen";
const std::string kCO2 = "CarbonDioxide";
const std::string kH2O = "Water";

constexpr double kFeedCO2 = 400e-6;    // [-]
constexpr double kFeedH2O = 0.0115;    // [-] 49 % RH at 20 degC and 1 bar
constexpr double kTemperature = 293.15;
constexpr double kPressure = 1.0e5;
constexpr double kFeedVelocity = 50e-6 / 5.03e-3;  // [m/s] 50 cm3/s through 50.3 cm2 (0.00994)
constexpr double kStepTime = 1.2e5;  // [s] a first run gave t95(CO2) = 1.03e5 s at the outlet

struct HumidDacDeck : UserDataTest {
  void TestBody() override {}

  void Build(double y_water) {
    SetUp();
    data.info.title = "Illustrative humid DAC breakthrough - 400 ppm CO2 and 49 % RH on APDES-NFC";
    data.info.description =
        "ILLUSTRATIVE, not a validated reproduction. Non-isothermal breakthrough of 400 ppm CO2 in humid N2 "
        "(y_H2O = 0.0115, 49 % RH at 20 degC, 1 bar) on the library adsorbent APDES-NFC: GAB water isotherm, "
        "humid Toth CO2 isotherm that reads the adsorbed water loading, LDF kinetics and heats of adsorption "
        "from the entry (Stampi-Bombelli et al. 2020, corrected 2022). Feed (50 cm3/s through an 8 cm column, "
        "0.01 m/s), column diameter and wall heat transfer as in that paper; the bed length and the step time "
        "are chosen for the example. Water saturates the bed first and releases its heat of sorption; the CO2 front follows with "
        "the humid capacity.";
    data.info.author = "Skarstrom examples";
    AddComponents({kN2, kCO2, kH2O});

    const AdsorbentID ads_id = data.adsorbent_directory.Add("APDES-NFC");
    auto& ads = data.adsorbent_directory.Get(ads_id);
    LoadAdsorbentFromJSON("../data/adsorbents/apdes_nfc_stampi2020.json", ads);
    ads.adsorption_source_.Initialise(data.fluid);

    const LayerID layer_id = data.bed.layers_.Add("APDES-NFC");
    auto& layer = data.bed.layers_.Get(layer_id);
    layer.domain_.length = 0.2;
    layer.domain_.n_cells = 25;
    layer.eb_ = 0.4;
    layer.particle_diameter_ = 7.5e-3;  // the entry's (Table 1 of the paper)
    layer.composition_[ads_id] = 1.0;
    layer.momentum_model_ = MomentumModel::Ergun;
    std::get<CylinderBed>(data.bed.shape_).diameter = 0.08;  // the paper's column
    data.bed.wall_.ApplyMaterial(WallMaterial::SS316L);
    data.bed.wall_.thickness_ = 2.0e-3;
    data.bed.ambient_temperature_ = kTemperature;
    AddConstantInternalHeatTransferCoefficientModel(3.0, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(26.0, layer_id);
    AddLinearDrivingForceModel(ads_id, kN2, 0.1);  // N2 is inert on the entry; a rate model for completeness

    data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::NON_ISOTHERMAL;
    AddInitialConditions(kPressure, kTemperature, kTemperature, {{kN2, 1.0}, {kCO2, 0.0}, {kH2O, 0.0}}, 0.0);
    AddStep("Adsorption", kStepTime);
    std::vector<std::unique_ptr<IBoundaryStream>> in;
    auto s = CreateNoValveInletStream(CreateVelocityBC(kFeedVelocity),
                                      {{kN2, 1.0 - kFeedCO2 - y_water}, {kCO2, kFeedCO2}, {kH2O, y_water}},
                                      kTemperature);
    s->name_ = "Humid air";
    in.push_back(std::move(s));
    AddFlowBC(std::move(in), BoundaryLocation::BOTTOM, BoundaryDirection::INLET, "Adsorption");
    AddPressureOutletBC(CreateConstantBoundaryValue(kPressure), BoundaryLocation::TOP, "Adsorption");
    ConfigureCycle(1);
    for (int c = 0; c < 3; ++c) ConfigureKPI(c, {0}, {BoundaryLocation::BOTTOM}, {0}, {BoundaryLocation::TOP});
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
  }
};

Simulation RunDeck(double y_water, const std::string& name, RunResult& result) {
  HumidDacDeck deck;
  deck.Build(y_water);
  std::filesystem::create_directories(deck.test_path);
  const std::filesystem::path path = deck.test_path / ("HumidDac_" + name + ".json");
  test::SaveToJsonAtomically(deck.data, path);
  if (y_water > 0.0)
    test::SaveToJsonAtomically(deck.data, std::filesystem::current_path() /
                                              "../examples/decks/Illustrative_Humid_DAC_Breakthrough.json");
  Simulation sim = Simulation::FromFile(path.string());
  result = sim.Run();
  return sim;
}

// First time the series reaches `level` (linear interpolation), -1 if never.
double TimeAt(const TimeSeries& s, double level) {
  for (size_t i = 1; i < s.values.size(); ++i)
    if (s.values[i] >= level)
      return s.time[i - 1] + (level - s.values[i - 1]) / (s.values[i] - s.values[i - 1]) * (s.time[i] - s.time[i - 1]);
  return -1.0;
}

double Last(const Simulation& sim, const char* var, const std::string& comp) {
  return sim.Series(var, comp, SeriesLocation::Top).values.back();
}

IIsotherm* IsothermOf(Simulation& sim, const std::string& comp) {
  auto& dir = sim.Data().adsorbent_directory;  // one adsorbent in the deck
  auto& ads = dir.Get(dir.GetAll().begin()->first);
  return ads.adsorption_source_.GetIsotherm(comp).get();
}

// Equilibrium loadings of the library entry at the feed state (20 degC, 1 bar).
struct FeedState {
  double q_water, q_co2_dry, q_co2_humid;
};

FeedState FeedLoadings(Simulation& sim) {
  FeedState f;
  f.q_water = EvalIsotherm(IsothermOf(sim, kH2O), kH2O, kTemperature, kPressure, kFeedH2O);
  f.q_co2_dry = EvalIsotherm(IsothermOf(sim, kCO2), kCO2, kTemperature, kPressure, kFeedCO2, 0.0);
  f.q_co2_humid = EvalIsotherm(IsothermOf(sim, kCO2), kCO2, kTemperature, kPressure, kFeedCO2, f.q_water);
  return f;
}

void ExpectBalancesClose(const Simulation& sim, bool humid) {
  const auto kpis = sim.KPIs();
  ASSERT_EQ(kpis.size(), 3u);  // N2, CO2 and water, one cycle
  for (const auto& k : kpis) {
    if (!humid && k.component == kH2O) continue;  // no water enters: the relative error is undefined (NaN)
    ASSERT_TRUE(std::isfinite(k.quadrature_mass_balance_error)) << k.component;
    EXPECT_LT(std::abs(k.quadrature_mass_balance_error), 1e-4) << k.component;
    EXPECT_LT(std::abs(k.mass_balance_error), 1e-3) << k.component;
  }
  const StepBalance* b = sim.Boundaries().GetStepBalance(sim.LastCycle(), 0);
  ASSERT_NE(b, nullptr);
  ASSERT_TRUE(b->valid);
  EXPECT_LT(std::abs(b->RelativeEndpointDiscrepancy()), 1e-5);
  for (size_t c = 0; c < b->components.size(); ++c) {
    const double scale = std::max(std::abs(b->InventoryChange(c)), b->endpoint_in[c] + b->endpoint_out[c]);
    if (scale < 1e-9) continue;  // water in the dry run: round-off only (~1e-17 mol)
    EXPECT_LT(std::abs(b->EndpointNetIn(c) - b->InventoryChange(c)) / scale, 1e-4) << b->components[c];
  }
}

}  // namespace

TEST(HumidDacBreakthrough, WaterFirstHumidCo2CapacityAndClosedBalances) {
  RunResult r;
  Simulation sim = RunDeck(kFeedH2O, "humid", r);
  ASSERT_TRUE(r.Ok()) << r.message;
  EXPECT_TRUE(r.healthy);
  for (const auto& w : r.health_warnings) ADD_FAILURE() << "run-health warning " << w.code << ": " << w.message;
  ExpectBalancesClose(sim, true);

  // Outlet fronts: water first, then CO2.
  TimeSeries y_w = sim.Series("mole_fraction", kH2O, SeriesLocation::Top);
  TimeSeries y_c = sim.Series("mole_fraction", kCO2, SeriesLocation::Top);
  for (double& v : y_w.values) v /= kFeedH2O;
  for (double& v : y_c.values) v /= kFeedCO2;
  const double tw50 = TimeAt(y_w, 0.5), tc50 = TimeAt(y_c, 0.5), tc95 = TimeAt(y_c, 0.95);
  const TimeSeries T_mid = sim.Series("temperature", "", SeriesLocation::Position, 0.1);
  const size_t i_max = std::max_element(T_mid.values.begin(), T_mid.values.end()) - T_mid.values.begin();
  std::printf("  humid DAC: t50(H2O) = %.0f s, t50(CO2) = %.0f s, t95(CO2) = %.0f s, end c/c0: H2O %.4f CO2 %.4f\n",
              tw50, tc50, tc95, y_w.values.back(), y_c.values.back());
  std::printf("  mid-bed T max %.2f K at %.0f s\n", T_mid.values[i_max], T_mid.time[i_max]);
  ASSERT_GT(tw50, 0.0);
  ASSERT_GT(tc50, 0.0);
  EXPECT_LT(tw50, tc50);
  EXPECT_GT(tc95, 0.0) << "the step should be long enough for the CO2 front to leave the bed";

  // Water sorption heat: the mid-bed temperature peaks before the CO2 front
  // reaches 50 % at the outlet, and the peak is a real rise.
  EXPECT_LT(T_mid.time[i_max], tc50);
  EXPECT_GT(T_mid.values[i_max], kTemperature + 0.5);

  // End of the step, last cell: the loadings against the isotherms at the
  // local state. The bed is back at ambient by then (|T - 293.15 K| < 0.05 K,
  // measured 0.001 K). The relative humidity, and with it the GAB loading,
  // moves by about 6 % per kelvin through p_sat(T), so 0.05 K allows 0.3 %
  // between the local and the feed-state GAB value: 1 % against the feed state.
  // The CO2 LDF (2e-4 1/s) leaves the loading 0.16 % below its equilibrium
  // value at the end of the step: 1 %.
  const AxialProfile T = sim.Profile("temperature", "", -1, 0);
  const AxialProfile P = sim.Profile("pressure", "", -1, 0);
  const AxialProfile yw = sim.Profile("mole_fraction", kH2O, -1, 0);
  const AxialProfile yc = sim.Profile("mole_fraction", kCO2, -1, 0);
  const AxialProfile qw = sim.Profile("adsorbed_concentration", kH2O, -1, 0);
  const AxialProfile qc = sim.Profile("adsorbed_concentration", kCO2, -1, 0);
  const size_t last = (T.rows - 1) * T.cols + (T.cols - 1);
  const double T_end = T.values[last], P_end = P.values[last];
  const double qw_star = EvalIsotherm(IsothermOf(sim, kH2O), kH2O, T_end, P_end, yw.values[last]);
  const double qc_star = EvalIsotherm(IsothermOf(sim, kCO2), kCO2, T_end, P_end, yc.values[last], qw.values[last]);
  const FeedState feed = FeedLoadings(sim);
  std::printf("  end, last cell: T %.3f K, q_H2O %.4f (GAB %.4f, feed %.4f), q_CO2 %.4f (humid Toth %.4f, "
              "feed %.4f, dry %.4f) mol/kg\n",
              T_end, qw.values[last], qw_star, feed.q_water, qc.values[last], qc_star, feed.q_co2_humid,
              feed.q_co2_dry);
  EXPECT_LT(std::abs(T_end - kTemperature), 0.05);
  EXPECT_NEAR(qw.values[last], qw_star, 2e-3 * qw_star);
  EXPECT_NEAR(qw.values[last], feed.q_water, 0.01 * feed.q_water);
  EXPECT_NEAR(qc.values[last], qc_star, 0.01 * qc_star);
  EXPECT_GT(qc.values[last] / feed.q_co2_dry, 2.0);

  // The same deck with a dry feed: the dry capacity, well below the humid one.
  RunResult rd;
  Simulation dry = RunDeck(0.0, "dry", rd);
  ASSERT_TRUE(rd.Ok()) << rd.message;
  EXPECT_TRUE(rd.healthy);
  ExpectBalancesClose(dry, false);
  const AxialProfile qc_dry = dry.Profile("adsorbed_concentration", kCO2, -1, 0);
  const AxialProfile qw_dry = dry.Profile("adsorbed_concentration", kH2O, -1, 0);
  const double q_dry_last = qc_dry.values[(qc_dry.rows - 1) * qc_dry.cols + qc_dry.cols - 1];
  std::printf("  dry run, last cell: q_CO2 %.4f mol/kg\n", q_dry_last);
  EXPECT_LT(*std::max_element(qw_dry.values.begin(), qw_dry.values.end()), 1e-9);
  EXPECT_NEAR(q_dry_last, feed.q_co2_dry, 0.01 * feed.q_co2_dry);
  EXPECT_LT(q_dry_last, qc.values[last]);
  // The dry front leaves the bed earlier: less capacity at the same feed.
  TimeSeries y_c_dry = dry.Series("mole_fraction", kCO2, SeriesLocation::Top);
  for (double& v : y_c_dry.values) v /= kFeedCO2;
  const double tc50_dry = TimeAt(y_c_dry, 0.5);
  std::printf("  dry run: t50(CO2) = %.0f s\n", tc50_dry);
  ASSERT_GT(tc50_dry, 0.0);
  EXPECT_LT(tc50_dry, tc50);
}

// Spot checks of the library entry's isotherms at the paper's feed state.
TEST(HumidDacBreakthrough, LibraryEntryLoadingsAtTheFeedState) {
  HumidDacDeck deck;
  deck.Build(kFeedH2O);
  std::string json;
  SaveToString(deck.data, json);
  Simulation sim = Simulation::FromJson(json);
  const FeedState f = FeedLoadings(sim);

  // GAB in relative humidity (temperature-independent constants), by hand.
  const double x = kFeedH2O * kPressure / 2339.2;  // p_sat(20 degC) = 2339.2 Pa (IAPWS)
  const double C = 0.1489, K = 0.5751, qm = 36.48;
  const double gab = qm * C * K * x / ((1 - K * x) * (1 - K * x + C * K * x));
  EXPECT_NEAR(f.q_water, gab, 2e-3 * gab);

  // Dry Toth 2 at 40 Pa, 20 degC, by hand (Eqs. 7-10 with the corrected b0).
  const double T0 = 296.0, R = 8.314462618;
  const double b = 0.0707 * std::exp(57047.0 / (R * T0) * (T0 / kTemperature - 1.0));
  const double w = 0.4148 - 1.606 * (1.0 - T0 / kTemperature);
  const double p = kFeedCO2 * kPressure;
  const double toth = 2.38 * b * p / std::pow(1.0 + std::pow(b * p, w), 1.0 / w);
  EXPECT_NEAR(f.q_co2_dry, toth, 1e-9 * toth);
  EXPECT_NEAR(f.q_co2_dry, 0.82, 0.01);

  // Humid: the paper prints an equilibrium CO2 loading of about 2.02 mmol/g at
  // the humid feed. That value is read at the bed exit of its simulation, where
  // the temperature and the adsorbed water are not exactly the feed-state
  // values; the corrected parameters give 2.07 at exactly 20 degC and the GAB
  // water loading, within 5 % of the printed value.
  EXPECT_NEAR(f.q_co2_humid, 2.02, 0.05 * 2.02);
  EXPECT_GT(f.q_co2_humid / f.q_co2_dry, 2.0);
}

}  // namespace adsim
