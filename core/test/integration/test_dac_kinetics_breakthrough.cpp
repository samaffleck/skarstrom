// Illustrative direct-air-capture breakthrough with the two DAC kinetics:
// the dual-rate LDF (a fast and a slow site sharing q*) with the
// series-resistance coefficient (film + macropore + Arrhenius micropore, with
// the capacity factor from the isotherm slope), in one end-to-end case.
// Writes examples/decks/Illustrative_DAC_Breakthrough.json.
//
// ILLUSTRATIVE ONLY. The sorbent is the library entry "Lewatit VP OC 1065"
// (data/adsorbents/lewatit_vp_oc_1065_low2023.json, loaded unchanged;
// CO2 isotherm: the dry site of the Young et al. (2021) WADST fit, doi:10.1039/D1EE01272J,
// because the deck has no Water component; N2 from doi:10.1021/acs.jced.3c00401). With
// the earlier dry td-Toth fit of Low et al. the outlet reaches 5/50/95 % of the feed at
// 6409/7393/27311 s; with the Young dry site at 7445/8188/29287 s. The feed is 400 ppm CO2 in N2 at 1 bar and
// 298.15 K, isothermal. The rate parameters (split, tortuosity, particle and
// micropore data) are round numbers chosen to show the two-rate behaviour of
// DAC breakthrough curves (Stampi-Bombelli et al. 2024,
// doi:10.1021/acs.iecr.4c01392), not fitted values.
//
// Checks (fast tier: 20 cells, one feed step):
//   * the run completes and the step and KPI mole balances close;
//   * the outlet shows an early front and a slow tail: the outlet CO2 reaches
//     50 % of the feed well before 95 %;
//   * the same deck with split = 1, K = 1 and the micropore term off is the
//     plain LDF + Haynes deck bit for bit.

#include "../validation/validation_helpers.h"

#include <filesystem>

#include "adsim/objects/deck_builder.h"
#include "adsim/simulation.h"

namespace adsim {

using namespace adsim::validation;

namespace {

constexpr double kFeedCO2 = 400e-6;  // [-] mole fraction
constexpr double kTemperature = 298.15;
constexpr double kPressure = 1.0e5;
constexpr double kFeedVelocity = 0.1;  // [m/s] superficial
constexpr double kStepTime = 6.0e4;    // [s]

enum class Variant {
  DualSeries,         // the example: dual_rate_ldf + series_resistance
  DualSeriesReduced,  // split = 1, K = 1, micropore off
  LdfHaynes,          // the plain reference
};

struct DacDeck : UserDataTest {
  void TestBody() override {}

  void Build(Variant variant) {
    SetUp();
    data.info.title = "Illustrative DAC breakthrough - 400 ppm CO2 on Lewatit VP OC 1065";
    data.info.description =
        "ILLUSTRATIVE, not a validated reproduction. Isothermal breakthrough of 400 ppm CO2 in N2 on the "
        "library adsorbent Lewatit VP OC 1065, with the dual-rate LDF (fast and slow site) and the "
        "series-resistance coefficient (film, macropore and micropore resistances with the capacity factor). "
        "Rate parameters are round numbers chosen to show an early front and a slow tail.";
    data.info.author = "Skarstrom examples";
    AddComponents({"Nitrogen", "CarbonDioxide"});

    const AdsorbentID ads_id = data.adsorbent_directory.Add("Lewatit VP OC 1065");
    auto& ads = data.adsorbent_directory.Get(ads_id);
    LoadAdsorbentFromJSON("../data/adsorbents/lewatit_vp_oc_1065_low2023.json", ads);
    ads.adsorption_source_.Initialise(data.fluid);

    const LayerID layer_id = data.bed.layers_.Add("Lewatit");
    auto& layer = data.bed.layers_.Get(layer_id);
    layer.domain_.length = 0.05;
    layer.domain_.n_cells = 20;
    layer.eb_ = 0.4;
    layer.particle_diameter_ = 6.0e-4;
    layer.composition_[ads_id] = 1.0;
    layer.momentum_model_ = MomentumModel::Ergun;
    std::get<CylinderBed>(data.bed.shape_).diameter = 0.02;
    data.bed.wall_.ApplyMaterial(WallMaterial::SS316L);
    data.bed.wall_.thickness_ = 3.1e-3;
    AddConstantInternalHeatTransferCoefficientModel(0.0, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);

    // Rate models. N2 is a constant LDF (hardly adsorbed). CO2:
    SetRateModels(variant, ads_id);

    data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
    AddInitialConditions(kPressure, kTemperature, kTemperature, {{"Nitrogen", 1.0}, {"CarbonDioxide", 0.0}}, 0.0);
    AddStep("Feed", kStepTime);
    std::vector<std::unique_ptr<IBoundaryStream>> in;
    auto s = CreateNoValveInletStream(CreateVelocityBC(kFeedVelocity),
                                      {{"Nitrogen", 1.0 - kFeedCO2}, {"CarbonDioxide", kFeedCO2}}, kTemperature);
    s->name_ = "Air";
    in.push_back(std::move(s));
    AddFlowBC(std::move(in), BoundaryLocation::BOTTOM, BoundaryDirection::INLET, "Feed");
    AddPressureOutletBC(CreateConstantBoundaryValue(kPressure), BoundaryLocation::TOP, "Feed");
    ConfigureCycle(1);
    ConfigureKPI(1, {0}, {BoundaryLocation::BOTTOM}, {0}, {BoundaryLocation::TOP});
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
  }

  void SetRateModels(Variant variant, AdsorbentID ads_id) {
    auto& models = GetAdsorbent(ads_id).adsorption_source_.mass_transfer_models_;
    AddLinearDrivingForceModel(ads_id, "Nitrogen", 0.1);
    auto series = [](double dc0, double rc, double capacity) {
      auto m = std::make_unique<MassTransferCoefficientSeriesResistance>();
      m->tortuosity_ = 3.0;
      m->macropore_void_fraction_ = 0.338;
      m->crystal_radius_ = rc;
      m->micropore_diffusivity_ = dc0;
      m->micropore_activation_energy_ = 2.0e4;
      m->capacity_factor_ = capacity;
      m->AdoptVersion1Switches();  // terms switched off by zero parameters, as in the version-1 model
      return m;
    };
    if (variant == Variant::LdfHaynes) {
      auto ldf = std::make_unique<LinearDrivingForce>();
      auto h = std::make_unique<MassTransferCoefficient_Haynes>();
      h->tortuosity_ = 3.0;
      h->macropore_void_fraction_ = 0.338;
      ldf->mtc_model_ = std::move(h);
      models["CarbonDioxide"] = std::move(ldf);
      return;
    }
    const bool reduced = variant == Variant::DualSeriesReduced;
    auto dual = std::make_unique<DualRateLinearDrivingForce>();
    dual->mtc_model_ = series(reduced ? 0.0 : kFastDc0, reduced ? 0.0 : kCrystalRadius, reduced ? 1.0 : 0.0);
    dual->slow_mtc_model_ = series(reduced ? 0.0 : kSlowDc0, reduced ? 0.0 : kCrystalRadius, reduced ? 1.0 : 0.0);
    dual->split_ = reduced ? 1.0 : kSplit;
    models["CarbonDioxide"] = std::move(dual);
  }

  static constexpr double kSplit = 0.7;
  static constexpr double kCrystalRadius = 1.0e-5;  // [m]
  static constexpr double kFastDc0 = 2.0e-10;       // [m2/s]
  static constexpr double kSlowDc0 = 1.0e-12;       // [m2/s]
};

// Outlet CO2 over the feed, c/c0 at the top of the bed.
struct Breakthrough {
  std::vector<double> t, y;
  double TimeAt(double fraction) const {
    for (size_t i = 1; i < y.size(); ++i)
      if (y[i] >= fraction) return t[i - 1] + (fraction - y[i - 1]) / (y[i] - y[i - 1]) * (t[i] - t[i - 1]);
    return -1.0;
  }
};

Breakthrough OutletOf(const Simulation& sim) {
  const TimeSeries y = sim.Series("mole_fraction", "CarbonDioxide", SeriesLocation::Top);
  Breakthrough b;
  b.t = y.time;
  for (double v : y.values) b.y.push_back(v / kFeedCO2);
  return b;
}

Simulation RunVariant(Variant variant, const std::string& name, RunResult& result) {
  DacDeck deck;
  deck.Build(variant);
  std::filesystem::create_directories(deck.test_path);
  const std::filesystem::path path = deck.test_path / ("DacKinetics_" + name + ".json");
  test::SaveToJsonAtomically(deck.data, path);
  if (variant == Variant::DualSeries)
    test::SaveToJsonAtomically(deck.data,
                               std::filesystem::current_path() / "../examples/decks/Illustrative_DAC_Breakthrough.json");
  Simulation sim = Simulation::FromFile(path.string());
  result = sim.Run();
  return sim;
}

}  // namespace

TEST(DacKineticsBreakthrough, FastFrontAndSlowTailWithClosedBalances) {
  RunResult r;
  Simulation sim = RunVariant(Variant::DualSeries, "dual_series", r);
  ASSERT_TRUE(r.Ok()) << r.message;

  // Mole balances of the engine's diagnostics.
  const auto kpis = sim.KPIs();
  ASSERT_FALSE(kpis.empty());
  for (const auto& k : kpis) {
    ASSERT_TRUE(std::isfinite(k.mass_balance_error));
    EXPECT_LT(std::abs(k.mass_balance_error), 1e-4);
    EXPECT_LT(std::abs(k.trapezoid_mass_balance_error), 1e-3);
  }
  const StepBalance* b = sim.Boundaries().GetStepBalance(sim.LastCycle(), 0);
  ASSERT_NE(b, nullptr);
  ASSERT_TRUE(b->valid);
  EXPECT_LT(std::abs(b->RelativeEndpointDiscrepancy()), 1e-5);

  // Shape: monotone, an early front and a slow tail.
  const Breakthrough out = OutletOf(sim);
  for (size_t i = 1; i < out.y.size(); ++i) EXPECT_GE(out.y[i], out.y[i - 1] - 1e-4);
  const double t05 = out.TimeAt(0.05), t50 = out.TimeAt(0.5), t95 = out.TimeAt(0.95);
  std::printf("  DAC breakthrough: t5 = %.0f s, t50 = %.0f s, t95 = %.0f s, c/c0 at end %.4f\n", t05, t50, t95,
              out.y.back());
  ASSERT_GT(t50, 0.0);
  ASSERT_GT(t95, 0.0);
  EXPECT_GT(t05, 0.0);
  EXPECT_LT(t50, 0.5 * t95);
}

TEST(DacKineticsBreakthrough, UnitSplitUnitCapacityNoMicroporeIsLdfHaynesBitForBit) {
  RunResult ra, rb;
  Simulation haynes = RunVariant(Variant::LdfHaynes, "ldf_haynes", ra);
  Simulation reduced = RunVariant(Variant::DualSeriesReduced, "dual_series_reduced", rb);
  ASSERT_TRUE(ra.Ok()) << ra.message;
  ASSERT_TRUE(rb.Ok()) << rb.message;
  EXPECT_EQ(ra.stats.num_steps, rb.stats.num_steps);
  EXPECT_EQ(ra.stats.num_res_evals, rb.stats.num_res_evals);
  for (const char* var : {"pressure", "temperature"})
    for (auto where : {SeriesLocation::Top, SeriesLocation::Bottom})
      EXPECT_EQ(haynes.Series(var, "", where).values, reduced.Series(var, "", where).values) << var;
  for (const char* var : {"mole_fraction", "adsorbed_concentration"})
    for (auto where : {SeriesLocation::Top, SeriesLocation::Bottom})
      EXPECT_EQ(haynes.Series(var, "CarbonDioxide", where).values, reduced.Series(var, "CarbonDioxide", where).values)
          << var;
  const auto a = OutletOf(haynes);
  const auto c = OutletOf(reduced);
  EXPECT_EQ(a.t, c.t);
}

}  // namespace adsim
