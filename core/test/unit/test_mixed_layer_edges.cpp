// Edges of mixed adsorbent layers (a layer packed with several solids at mass
// fractions w_k): the composition rules (CheckLayerComposition), how the deck
// checks and the run-start warnings report them, and removing a solid that a
// blend uses (RemoveAdsorbentFromDeck). The solved edges (a solid at w = 0, an
// all-inert layer, interfaces between blends) are in
// validation/test_mixed_layers.cpp.

#include <filesystem>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/objects/deck_edit.h"
#include "adsim/objects/porous_media.h"
#include "adsim/objects/user_data.h"
#include "adsim/simulation.h"
#include "adsim/solver/deck_warnings.h"

namespace adsim {
namespace {

  using Severity = LayerCompositionIssue::Severity;

  const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";

  struct Blend {
    AdsorbentDirectory dir;
    MixedAdsorbentLayer layer{"Bottom"};
    AdsorbentID a = dir.Add("Zeolite 13X");
    AdsorbentID b = dir.Add("Alumina");
  };

  size_t Count(const std::vector<LayerCompositionIssue>& issues, Severity s) {
    size_t n = 0;
    for (const auto& i : issues) n += i.severity == s;
    return n;
  }

  TEST(MixedLayerEdges, CompositionSummingToOneIsClean) {
    Blend x;
    x.layer.AddAdsorbent(x.a, 0.7);
    x.layer.AddAdsorbent(x.b, 0.3);
    EXPECT_TRUE(CheckLayerComposition(x.layer, x.dir).empty());
    EXPECT_TRUE(x.layer.ValidateComposition());
    // Within the tolerance (e.g. three fractions typed as 0.3333333).
    x.layer.composition_[x.a] = 0.7 + 0.5 * kCompositionSumTolerance;
    EXPECT_TRUE(CheckLayerComposition(x.layer, x.dir).empty());
    EXPECT_TRUE(x.layer.ValidateComposition());
  }

  TEST(MixedLayerEdges, CompositionNotSummingToOneIsAnError) {
    for (double second : {0.1, 0.5}) {  // sums 0.8 and 1.2
      Blend x;
      x.layer.AddAdsorbent(x.a, 0.7);
      x.layer.AddAdsorbent(x.b, second);
      const auto issues = CheckLayerComposition(x.layer, x.dir);
      ASSERT_EQ(issues.size(), 1u);
      EXPECT_EQ(issues[0].severity, Severity::Error);
      EXPECT_NE(issues[0].message.find("Layer \"Bottom\": mass fractions sum to " +
                                       std::string(second < 0.3 ? "0.8" : "1.2") + "; they must sum to 1"),
                std::string::npos) << issues[0].message;
      // The hint for studies: vary the fractions together.
      EXPECT_NE(issues[0].message.find("set them to x and 1 - x together"), std::string::npos) << issues[0].message;
      EXPECT_FALSE(x.layer.ValidateComposition());
    }
    // A single solid below 1 is the same error (it would scale the solid).
    Blend x;
    x.layer.AddAdsorbent(x.a, 0.5);
    EXPECT_EQ(Count(CheckLayerComposition(x.layer, x.dir), Severity::Error), 1u);
  }

  TEST(MixedLayerEdges, NegativeFractionIsAnError) {
    Blend x;
    x.layer.AddAdsorbent(x.a, 1.2);
    x.layer.AddAdsorbent(x.b, -0.2);  // sums to 1, still invalid
    const auto issues = CheckLayerComposition(x.layer, x.dir);
    ASSERT_EQ(issues.size(), 1u);
    EXPECT_EQ(issues[0].severity, Severity::Error);
    EXPECT_NE(issues[0].message.find("adsorbent \"Alumina\" the mass fraction -0.2"), std::string::npos)
        << issues[0].message;
    EXPECT_FALSE(x.layer.ValidateComposition());
  }

  TEST(MixedLayerEdges, ZeroFractionIsAWarning) {
    Blend x;
    x.layer.AddAdsorbent(x.a, 1.0);
    x.layer.AddAdsorbent(x.b, 0.0);
    const auto issues = CheckLayerComposition(x.layer, x.dir);
    ASSERT_EQ(issues.size(), 1u);
    EXPECT_EQ(issues[0].severity, Severity::Warning);
    EXPECT_NE(issues[0].message.find("\"Alumina\" has mass fraction 0 in layer \"Bottom\""), std::string::npos)
        << issues[0].message;
    EXPECT_TRUE(x.layer.ValidateComposition());

    // Every fraction 0: the warning per solid and the sum error.
    x.layer.composition_[x.a] = 0.0;
    const auto all_zero = CheckLayerComposition(x.layer, x.dir);
    EXPECT_EQ(Count(all_zero, Severity::Warning), 2u);
    EXPECT_EQ(Count(all_zero, Severity::Error), 1u);
  }

  TEST(MixedLayerEdges, EmptyAndUnknownSolidsAreErrors) {
    Blend x;
    auto issues = CheckLayerComposition(x.layer, x.dir);
    ASSERT_EQ(issues.size(), 1u);
    EXPECT_NE(issues[0].message.find("\"Bottom\" has no adsorbent selected"), std::string::npos);
    x.layer.AddAdsorbent(9999, 1.0);
    issues = CheckLayerComposition(x.layer, x.dir);
    ASSERT_EQ(issues.size(), 1u);
    EXPECT_NE(issues[0].message.find("references adsorbent id 9999"), std::string::npos);
  }

  // The deck checks (C++ runs, and Python's sim.check() through
  // Simulation::Check) report every composition error of every layer; the
  // run-start warnings report the warnings; an all-inert solid is clean.
  TEST(MixedLayerEdges, DeckChecksAndWarningsUseTheSameRule) {
    Simulation sim = Simulation::FromFile((kExamples / "Perez2019_A.json").string());
    std::string err;
    ASSERT_TRUE(sim.Check(err)) << err;
    UserData& d = sim.Data();
    ASSERT_FALSE(d.bed.layers_.layers_.empty());
    // By position: adding a layer may move the others.
    auto layer0 = [&]() -> MixedAdsorbentLayer& { return d.bed.layers_.layers_.front(); };
    ASSERT_EQ(layer0().composition_.size(), 1u);
    const AdsorbentID solid = layer0().composition_.begin()->first;
    const std::string name0 = layer0().name_;

    // An inert diluent (every component on the Inert isotherm) at 30 %.
    const AdsorbentID inert = AddAdsorbentToDeck(d, "Glass beads");
    for (const auto& c : d.fluid.components_)
      d.adsorbent_directory.Get(inert).adsorption_source_.SetIsothermType(c, IsothermType::INERT);
    layer0().composition_[solid] = 0.7;
    layer0().composition_[inert] = 0.3;
    EXPECT_TRUE(sim.Check(err)) << err;
    EXPECT_TRUE(LayerCompositionWarnings(d.adsorbent_directory, d.bed).empty());

    // Sum 0.8: refused, with the message the rule gives.
    layer0().composition_[inert] = 0.1;
    err.clear();
    EXPECT_FALSE(sim.Check(err));
    const auto issues = CheckLayerComposition(layer0(), d.adsorbent_directory);
    ASSERT_EQ(issues.size(), 1u);
    EXPECT_NE(err.find(issues[0].message), std::string::npos) << err;

    // A second bad layer is reported too.
    const LayerID top = d.bed.layers_.Add("Top");
    d.bed.layers_.Get(top).composition_[solid] = 0.5;
    err.clear();
    EXPECT_FALSE(sim.Check(err));
    EXPECT_NE(err.find("Layer \"" + name0 + "\": mass fractions sum to 0.8"), std::string::npos) << err;
    EXPECT_NE(err.find("Layer \"Top\": mass fractions sum to 0.5"), std::string::npos) << err;

    // A zero fraction passes the checks and is a run-start warning.
    d.bed.layers_.Remove(top);
    layer0().composition_[solid] = 1.0;
    layer0().composition_[inert] = 0.0;
    err.clear();
    EXPECT_TRUE(sim.Check(err)) << err;
    const auto warnings = LayerCompositionWarnings(d.adsorbent_directory, d.bed);
    ASSERT_EQ(warnings.size(), 1u);
    EXPECT_NE(warnings[0].find("\"Glass beads\" has mass fraction 0"), std::string::npos) << warnings[0];
  }

  // Removing a solid that a blend uses is refused with a message (the user
  // chooses the new composition; the engine never renormalises). Removing
  // the only solid of a layer re-seeds the layer, as before.
  TEST(MixedLayerEdges, RemovingASolidOfABlendIsRefused) {
    UserData d;
    LoadFromJson(d, (kExamples / "Perez2019_A.json").string());
    const AdsorbentID solid = d.bed.layers_.layers_.front().composition_.begin()->first;
    const std::string bottom = d.bed.layers_.layers_.front().name_;
    const AdsorbentID inert = AddAdsorbentToDeck(d, "Glass beads");
    const AdsorbentID spare = AddAdsorbentToDeck(d, "Spare");
    d.bed.layers_.layers_.front().composition_ = {{solid, 0.7}, {inert, 0.3}};
    const LayerID top_id = d.bed.layers_.Add("Top");
    d.bed.layers_.Get(top_id).composition_ = {{solid, 0.5}, {inert, 0.5}};

    std::string err;
    EXPECT_FALSE(RemoveAdsorbentFromDeck(d, inert, err));
    EXPECT_EQ(err, "Glass beads is part of the blend of layers \"" + bottom +
                   "\" (mass fraction 0.3); \"Top\" (mass fraction 0.5); remove it from the layers' compositions first");
    EXPECT_EQ(err, AdsorbentInUseByMixtures(d, inert));
    // Nothing changed.
    EXPECT_TRUE(d.adsorbent_directory.Contains(inert));
    EXPECT_DOUBLE_EQ(d.bed.layers_.layers_.front().composition_.at(inert), 0.3);
    EXPECT_DOUBLE_EQ(d.bed.layers_.Get(top_id).composition_.at(inert), 0.5);

    // An unused solid goes; so does an unknown id (no-op).
    EXPECT_TRUE(AdsorbentInUseByMixtures(d, spare).empty());
    EXPECT_TRUE(RemoveAdsorbentFromDeck(d, spare, err)) << err;
    EXPECT_FALSE(d.adsorbent_directory.Contains(spare));
    EXPECT_TRUE(RemoveAdsorbentFromDeck(d, 9999, err));

    // Once the user has taken it out of both blends, it can be removed.
    d.bed.layers_.layers_.front().composition_ = {{solid, 1.0}};
    d.bed.layers_.Get(top_id).composition_ = {{inert, 1.0}};
    EXPECT_TRUE(AdsorbentInUseByMixtures(d, inert).empty());
    EXPECT_TRUE(RemoveAdsorbentFromDeck(d, inert, err)) << err;
    EXPECT_FALSE(d.adsorbent_directory.Contains(inert));
    // "Top" held only the removed solid: re-seeded with the remaining one.
    const auto& top = d.bed.layers_.Get(top_id).composition_;
    ASSERT_EQ(top.size(), 1u);
    EXPECT_EQ(top.begin()->first, solid);
    EXPECT_DOUBLE_EQ(top.begin()->second, 1.0);
    err.clear();
    EXPECT_TRUE(d.PerformChecks(err)) << err;
  }

  double Sum(const MixedAdsorbentLayer& layer) {
    double s = 0.0;
    for (const auto& [id, w] : layer.composition_) s += w;
    return s;
  }

  // The blend edits keep the fractions summing to 1, so a front end that
  // builds blends through them never makes a deck the check refuses.
  TEST(MixedLayerEdges, BlendEditsKeepTheSumAtOne) {
    UserData d;
    LoadFromJson(d, (kExamples / "Perez2019_A.json").string());
    auto layer = [&]() -> MixedAdsorbentLayer& { return d.bed.layers_.layers_.front(); };
    const AdsorbentID a = layer().composition_.begin()->first;
    const AdsorbentID b = AddAdsorbentToDeck(d, "B");  // no empty layer: not added anywhere
    const AdsorbentID c = AddAdsorbentToDeck(d, "C");
    EXPECT_EQ(layer().composition_.size(), 1u);
    SeedLayerComposition(d, layer());  // a layer that holds a solid is left alone
    EXPECT_EQ(layer().composition_.size(), 1u);
    EXPECT_DOUBLE_EQ(Sum(layer()), 1.0);

    std::string err;
    ASSERT_TRUE(AddSolidToLayer(d, layer(), b, 0.3, err)) << err;  // a: 1 -> 0.7
    EXPECT_DOUBLE_EQ(layer().composition_.at(a), 0.7);
    EXPECT_DOUBLE_EQ(layer().composition_.at(b), 0.3);
    ASSERT_TRUE(AddSolidToLayer(d, layer(), c, 0.5, err)) << err;  // a, b scaled by 0.5
    EXPECT_DOUBLE_EQ(layer().composition_.at(a), 0.35);
    EXPECT_DOUBLE_EQ(layer().composition_.at(b), 0.15);
    EXPECT_NEAR(Sum(layer()), 1.0, 1e-15);
    ASSERT_TRUE(SetLayerMassFraction(d, layer(), a, 0.6, err)) << err;  // b : c stays 0.15 : 0.5
    EXPECT_NEAR(layer().composition_.at(b) / layer().composition_.at(c), 0.3, 1e-14);
    EXPECT_NEAR(Sum(layer()), 1.0, 1e-15);
    ASSERT_TRUE(AddSolidToLayer(d, layer(), a, 0.2, err)) << err;  // already held: a set
    EXPECT_DOUBLE_EQ(layer().composition_.at(a), 0.2);
    EXPECT_NEAR(Sum(layer()), 1.0, 1e-15);
    ASSERT_TRUE(RemoveSolidFromLayer(d, layer(), c, err)) << err;  // a : b kept, sum 1
    EXPECT_NEAR(Sum(layer()), 1.0, 1e-15);
    EXPECT_TRUE(CheckLayerComposition(layer(), d.adsorbent_directory).empty());

    // Others at 0 share the rest equally.
    ASSERT_TRUE(SetLayerMassFraction(d, layer(), a, 1.0, err)) << err;
    EXPECT_DOUBLE_EQ(layer().composition_.at(b), 0.0);
    ASSERT_TRUE(SetLayerMassFraction(d, layer(), a, 0.4, err)) << err;
    EXPECT_DOUBLE_EQ(layer().composition_.at(b), 0.6);

    // Refusals change nothing.
    const auto before = layer().composition_;
    EXPECT_FALSE(AddSolidToLayer(d, layer(), 9999, 0.5, err));
    EXPECT_FALSE(AddSolidToLayer(d, layer(), c, 0.0, err));
    EXPECT_FALSE(AddSolidToLayer(d, layer(), c, 1.5, err));
    EXPECT_FALSE(SetLayerMassFraction(d, layer(), a, -0.1, err));
    EXPECT_FALSE(SetLayerMassFraction(d, layer(), c, 0.5, err));  // not in the layer
    EXPECT_EQ(layer().composition_, before);
    ASSERT_TRUE(RemoveSolidFromLayer(d, layer(), b, err)) << err;
    EXPECT_FALSE(RemoveSolidFromLayer(d, layer(), a, err));  // the only solid
    EXPECT_FALSE(SetLayerMassFraction(d, layer(), a, 0.5, err));
    EXPECT_NE(err.find("only solid"), std::string::npos) << err;
    EXPECT_DOUBLE_EQ(layer().composition_.at(a), 1.0);

    // An empty layer: the first solid gets 1, whatever was asked.
    const LayerID top = d.bed.layers_.Add("Top");
    d.bed.layers_.Get(top).composition_.clear();
    ASSERT_TRUE(AddSolidToLayer(d, d.bed.layers_.Get(top), b, 0.3, err)) << err;
    EXPECT_DOUBLE_EQ(d.bed.layers_.Get(top).composition_.at(b), 1.0);
    err.clear();
    EXPECT_TRUE(d.PerformChecks(err)) << err;
  }

}  // namespace
}  // namespace adsim
