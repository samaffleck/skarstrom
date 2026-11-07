// test_support_solids.cpp
//
// The non-adsorbing support solids of the adsorbent library (family
// "support": cordierite, aluminium, stainless steel, ceramic paper) in a
// structured layer. Each support is added to the Perez 2019 A deck from the
// library, mixed by mass with the deck's 13X in the layer's composition, the
// layer is made a monolith (SetLayerKind) and the deck runs one short cycle.
// Checked: the support adsorbs nothing, the layer's solid density is the
// mass-weighted mixture, and the run completes with finite KPIs. Fast tier: a
// coarse grid and one cycle with the CSS check off.

#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/objects/adsorbent_library.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/deck_edit.h"
#include "adsim/simulation.h"

namespace adsim {

  namespace {
    const std::filesystem::path kLibDir = std::filesystem::current_path() / ".." / "data" / "adsorbents";
    const std::string kDeck = "../examples/decks/Perez2019_A.json";
    constexpr double kSupportFraction = 0.6;  // mass fraction of the support in the layer
  }  // namespace

  TEST(SupportSolids, EverySupportRunsInAMonolithLayer) {
    const AdsorbentLibrary lib = LoadAdsorbentLibraryFile((kLibDir / "library.json").string());
    std::vector<std::string> supports;
    for (const auto& m : lib.materials)
      if (m.family == "support") supports.push_back(m.id);
    ASSERT_GE(supports.size(), 4u) << "cordierite, aluminium, stainless steel and ceramic paper";

    for (const auto& material : supports) {
      SCOPED_TRACE(material);
      Simulation sim = Simulation::FromFile(kDeck);
      UserData& data = sim.Data();
      auto& layer = data.bed.layers_.layers_.front();
      ASSERT_EQ(layer.composition_.size(), 1u);
      const AdsorbentID sorbent = layer.composition_.begin()->first;

      AdsorbentID support{};
      std::string err;
      ASSERT_TRUE(AddLibraryAdsorbent(data, kLibDir.string(), material, support, err)) << err;
      Adsorbent& s = data.adsorbent_directory.Get(support);
      EXPECT_EQ(s.ep_.Get(), 0.0);
      EXPECT_GT(s.thermal_conductivity_.Get(), 0.0);
      // Every deck component is inert on the support.
      for (const auto& c : data.fluid.components_) {
        const auto& iso = s.adsorption_source_.GetIsotherm(c);
        ASSERT_TRUE(iso) << c;
        EXPECT_EQ(iso->GetIsothermType(), IsothermType::INERT) << c;
      }

      const LayerID id = data.bed.layers_.layers_.front().id_;
      build::SetLayerComposition(data, id, {{sorbent, 1.0 - kSupportFraction}, {support, kSupportFraction}});
      build::SetLayerKind(data, id, "monolith");
      auto& l = data.bed.layers_.layers_.front();
      l.domain_.n_cells = 20;

      // Mass-weighted mixture: 1 / rho = sum w_i / rho_i.
      const double rho_sorbent = data.adsorbent_directory.Get(sorbent).solid_density_.Get();
      const double rho_mix = 1.0 / ((1.0 - kSupportFraction) / rho_sorbent + kSupportFraction / s.solid_density_.Get());
      EXPECT_NEAR(l.GetMixtureSolidDensity(data.adsorbent_directory), rho_mix, 1e-9 * rho_mix);

      data.settings.css.method_ = CSSMethod::NONE;
      RunOptions options;
      options.cycles = 1;
      options.on_message = [](LogLevel, const std::string&) {};
      const RunResult run = sim.Run(options);
      ASSERT_TRUE(run.Ok()) << run.message;
      EXPECT_EQ(run.cycles_run, 1);
      const auto kpis = sim.KPIs(-1);
      ASSERT_FALSE(kpis.empty());
      for (const auto& k : kpis) {
        EXPECT_TRUE(std::isfinite(k.purity)) << k.component;
        EXPECT_TRUE(std::isfinite(k.recovery)) << k.component;
      }
    }
  }

}  // namespace adsim
