// test_heat_dispersion_kinds.cpp
//
// The stagnant + dynamic heat-dispersion kinds (dispersion_model.h):
//   * krupiczka at its default dynamic coefficient (0.75) is the yagi_kunii
//     model: on the bundled yagi_kunii decks the two give bitwise equal
//     results;
//   * zehner_schlunder runs and gives a different temperature history;
//   * the dynamic coefficient is a registry parameter (layer[ID].dynamic_coefficient)
//     and keeps its value across a switch between the two kinds;
//   * yagi_kunii still loads and is still accepted by the builder, but is
//     not offered.
// Short grids and one cycle: these check that the feature works.

#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/equations/dispersion_model.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/deck_edit.h"
#include "adsim/objects/model_kinds.h"
#include "adsim/simulation.h"

namespace adsim {

namespace {

const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";

constexpr int kShortCells = 20;

Simulation ShortDeck(const std::string& name) {
  Simulation sim = Simulation::FromFile((kExamples / name).string());
  sim.SetParameter("cycle.number_of_cycles", 1);
  for (auto& layer : sim.Data().bed.layers_.layers_) layer.domain_.n_cells = kShortCells;
  return sim;
}

std::vector<double> TopSeries(const Simulation& sim, const std::string& variable, const std::string& component = {}) {
  return sim.Series(variable, component, SeriesLocation::Top).values;
}

void ExpectBitwiseEqual(const std::vector<double>& a, const std::vector<double>& b, const std::string& what) {
  ASSERT_EQ(a.size(), b.size()) << what;
  ASSERT_FALSE(a.empty()) << what;
  for (size_t i = 0; i < a.size(); ++i) {
    // EXPECT_EQ on doubles is exact; one mismatch is enough to report.
    if (a[i] != b[i]) {
      ADD_FAILURE() << what << " differs at sample " << i << ": " << a[i] << " vs " << b[i];
      return;
    }
  }
}

}  // namespace

class KrupiczkaMatchesYagiKunii : public ::testing::TestWithParam<std::string> {};

TEST_P(KrupiczkaMatchesYagiKunii, Bitwise) {
  Simulation legacy = ShortDeck(GetParam());
  for (const auto& layer : legacy.Data().bed.layers_.layers_)
    ASSERT_NE(dynamic_cast<const HeatDispersionModelYagiandKunii*>(layer.equations_.heat_dispersion_model_.get()),
              nullptr)
        << "the deck is expected to use yagi_kunii (layer " << layer.name_ << ")";
  const RunResult a = legacy.Run();
  ASSERT_TRUE(a.Ok()) << a.message;

  Simulation fresh = ShortDeck(GetParam());
  for (auto& layer : fresh.Data().bed.layers_.layers_) {
    std::string err;
    ASSERT_TRUE(SetHeatDispersionKind(fresh.Data(), layer.id_, AxialHeatDispersionType::KRUPICZKA, err)) << err;
    auto* k = dynamic_cast<HeatDispersionModelKrupiczka*>(layer.equations_.heat_dispersion_model_.get());
    ASSERT_NE(k, nullptr);
    EXPECT_EQ(k->dynamic_coefficient_.Get(), 0.75);
  }
  const RunResult b = fresh.Run();
  ASSERT_TRUE(b.Ok()) << b.message;

  EXPECT_EQ(a.stats.num_steps, b.stats.num_steps);
  ExpectBitwiseEqual(TopSeries(legacy, "temperature"), TopSeries(fresh, "temperature"), "top temperature");
  ExpectBitwiseEqual(TopSeries(legacy, "pressure"), TopSeries(fresh, "pressure"), "top pressure");
  ExpectBitwiseEqual(TopSeries(legacy, "mole_fraction", "CarbonDioxide"),
                     TopSeries(fresh, "mole_fraction", "CarbonDioxide"), "top CO2 mole fraction");
}

// The bundled decks that use yagi_kunii.
INSTANTIATE_TEST_SUITE_P(HeatDispersionKinds, KrupiczkaMatchesYagiKunii,
                         ::testing::Values("Ward2022_Q100ccm.json", "Illustrative_Diluted_Bed.json"));

TEST(HeatDispersionKinds, ZehnerSchlunderRunsAndDiffers) {
  Simulation base = ShortDeck("Ward2022_Q100ccm.json");
  ASSERT_TRUE(base.Run().Ok());

  Simulation zs = ShortDeck("Ward2022_Q100ccm.json");
  for (auto& layer : zs.Data().bed.layers_.layers_) build::SetHeatDispersion(zs.Data(), layer.id_, "zehner_schlunder");
  const RunResult r = zs.Run();
  ASSERT_TRUE(r.Ok()) << r.message;
  EXPECT_TRUE(r.healthy);
  const auto t0 = TopSeries(base, "temperature");
  const auto t1 = TopSeries(zs, "temperature");
  ASSERT_FALSE(t1.empty());
  EXPECT_NE(t0, t1);
}

TEST(HeatDispersionKinds, DynamicCoefficientIsARegistryParameterAndCarriesOver) {
  Simulation sim = ShortDeck("Ward2022_Q100ccm.json");
  const LayerID id = sim.Data().bed.layers_.layers_.front().id_;
  build::SetHeatDispersion(sim.Data(), id, "krupiczka");
  sim.RefreshParameters();
  const std::string path = "layer[" + std::to_string(id) + "].dynamic_coefficient";
  ASSERT_NE(sim.Parameters().Find(path), nullptr) << path;
  EXPECT_EQ(sim.GetParameter(path), 0.75);
  sim.SetParameter(path, 0.6);

  // A switch keeps a parameter of the same name and unit.
  std::string err;
  ASSERT_TRUE(SetHeatDispersionKind(sim.Data(), id, AxialHeatDispersionType::ZEHNER_SCHLUNDER, err)) << err;
  sim.RefreshParameters();
  EXPECT_EQ(sim.GetParameter(path), 0.6);
}

TEST(HeatDispersionKinds, YagiKuniiIsLegacy) {
  const KindDescriptor* yk = FindKind(ModelFamily::HeatDispersion, "yagi_kunii");
  ASSERT_NE(yk, nullptr);
  EXPECT_TRUE(yk->legacy);
  const auto& offered = build::HeatDispersionKinds();
  EXPECT_EQ(std::find(offered.begin(), offered.end(), "yagi_kunii"), offered.end());
  // Still accepted by the builder, so old scripts keep working.
  Simulation sim = ShortDeck("Ward2022_Q100ccm.json");
  const LayerID id = sim.Data().bed.layers_.layers_.front().id_;
  build::SetHeatDispersion(sim.Data(), id, "zehner_schlunder");
  build::SetHeatDispersion(sim.Data(), id, "yagi_kunii");
  EXPECT_NE(dynamic_cast<HeatDispersionModelYagiandKunii*>(
                sim.Data().bed.layers_.layers_.front().equations_.heat_dispersion_model_.get()),
            nullptr);
}

}  // namespace adsim
