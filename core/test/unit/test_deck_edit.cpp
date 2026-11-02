/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// Tests for the deck I/O and structural-edit helpers that front ends use
// (adsim/objects/deck_io.h, adsim/objects/deck_edit.h).

#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "adsim/objects/deck_builder.h"
#include "adsim/objects/adsorbent_library.h"
#include "adsim/objects/deck_edit.h"
#include "adsim/objects/deck_io.h"
#include "adsim/objects/user_data.h"

namespace adsim {
namespace {

const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";
const std::filesystem::path kLibDir = std::filesystem::current_path() / ".." / "data" / "adsorbents";

std::string ReadText(const std::filesystem::path& p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

UserData LoadPerez() {
  UserData data;
  LoadDeckFlexible(data, ReadText(kExamples / "Perez2019_A.json"));
  return data;
}

}  // namespace

TEST(DeckIo, LoadsBothRootShapes) {
  // File decks use cereal's default "value0" root.
  UserData a = LoadPerez();
  ASSERT_FALSE(a.fluid.components_.empty());
  ASSERT_FALSE(a.cycle.step_list_.empty());

  // String-saved decks use the "simulation" root.
  std::string text;
  SaveToString(a, text);
  ASSERT_NE(text.find("\"simulation\""), std::string::npos);
  UserData b;
  LoadDeckFlexible(b, text);
  EXPECT_EQ(b.fluid.components_, a.fluid.components_);
  EXPECT_EQ(b.cycle.step_list_.size(), a.cycle.step_list_.size());

  UserData c;
  EXPECT_THROW(LoadDeckFlexible(c, ""), std::invalid_argument);
  EXPECT_ANY_THROW(LoadDeckFlexible(c, "{\"not\": \"a deck\"}"));
}

TEST(DeckEdit, AddAndRemoveComponentCascades) {
  UserData data = LoadPerez();
  std::string err;
  EXPECT_FALSE(AddComponentToDeck(data, "Unobtainium", err));
  EXPECT_FALSE(err.empty());

  ASSERT_TRUE(AddComponentToDeck(data, "Argon", err)) << err;
  EXPECT_NE(std::find(data.fluid.components_.begin(), data.fluid.components_.end(), "Argon"),
            data.fluid.components_.end());
  EXPECT_EQ(data.ic.mole_fraction.count("Argon"), 1u);
  // The adsorbed-phase heat capacity gets an editable 0 entry straight away
  // (Fluid::Initialise would only add one at run time).
  ASSERT_EQ(data.fluid.adsorbed_heat_capacity_.count("Argon"), 1u);
  EXPECT_EQ(data.fluid.adsorbed_heat_capacity_.at("Argon").Get(), 0.0);
  for (auto& [id, ads] : data.adsorbent_directory.GetAll()) {
    EXPECT_EQ(ads.adsorption_source_.mass_transfer_models_.count("Argon"), 1u) << ads.name_;
    EXPECT_EQ(ads.adsorption_source_.heat_of_adsorption_.count("Argon"), 1u) << ads.name_;
  }
  // Adding twice is a no-op.
  const auto n = data.fluid.components_.size();
  ASSERT_TRUE(AddComponentToDeck(data, "Argon", err));
  EXPECT_EQ(data.fluid.components_.size(), n);

  RemoveComponentFromDeck(data, "Argon");
  EXPECT_EQ(std::find(data.fluid.components_.begin(), data.fluid.components_.end(), "Argon"),
            data.fluid.components_.end());
  EXPECT_EQ(data.ic.mole_fraction.count("Argon"), 0u);
  EXPECT_EQ(data.fluid.adsorbed_heat_capacity_.count("Argon"), 0u);
  std::string eval_err;
  EXPECT_TRUE(data.EvaluateAllExpressions(eval_err)) << eval_err;
}

TEST(DeckEdit, LibraryAdsorbentAddAndRemove) {
  UserData data = LoadPerez();
  const size_t before = data.adsorbent_directory.GetAll().size();
  AdsorbentID id{};
  std::string err;
  EXPECT_FALSE(AddLibraryAdsorbent(data, kLibDir.string(), "no_such_key", id, err));
  ASSERT_TRUE(AddLibraryAdsorbent(data, kLibDir.string(), "zeolite_13x_perez2019", id, err)) << err;
  EXPECT_EQ(data.adsorbent_directory.GetAll().size(), before + 1);
  // A material id adds that material's recommended entry.
  AdsorbentID by_material{};
  ASSERT_TRUE(AddLibraryAdsorbent(data, kLibDir.string(), "zeolite-13x", by_material, err)) << err;
  const AdsorbentLibrary lib = LoadAdsorbentLibraryFile((kLibDir / "library.json").string());
  const AdsorbentLibraryEntry* rec = lib.Resolve("zeolite-13x");
  ASSERT_NE(rec, nullptr);
  EXPECT_TRUE(rec->recommended);
  Adsorbent rec_blob;
  LoadAdsorbentFromJSON((kLibDir / rec->file).string(), rec_blob);
  EXPECT_EQ(data.adsorbent_directory.Get(by_material).name_, rec_blob.name_);
  RemoveAdsorbentFromDeck(data, by_material);
  EXPECT_FALSE(AddLibraryAdsorbent(data, kLibDir.string(), "no-such-material", by_material, err));
  // Every deck component has physics on the new adsorbent.
  auto& ads = data.adsorbent_directory.Get(id);
  for (const auto& c : data.fluid.components_)
    EXPECT_EQ(ads.adsorption_source_.mass_transfer_models_.count(c), 1u) << c;

  Adsorbent bad;
  EXPECT_FALSE(LoadAdsorbentFromString("not json", bad, err));

  // Removing every adsorbent the layers use re-seeds them with what remains.
  for (auto& layer : data.bed.layers_.layers_) {
    auto comp = layer.composition_;
    for (const auto& [ads_id, frac] : comp)
      if (ads_id != id) RemoveAdsorbentFromDeck(data, ads_id);
  }
  for (auto& layer : data.bed.layers_.layers_) {
    ASSERT_EQ(layer.composition_.size(), 1u);
    EXPECT_EQ(layer.composition_.begin()->first, id);
  }
}

}  // namespace adsim

namespace adsim {

// The edit cascades must leave nothing for run-time initialisation to fill in.
// Fluid::Initialise and Bed::Initialise seed a default for any per-component
// entry a deck lacks (an adsorbed-phase heat capacity, a heat of adsorption, an
// isotherm, a rate model), so a cascade that forgets one is invisible to the
// engine but shows up in an editor that lists the map's entries: the component
// has no field until the deck has been run once. Saving before and after
// Initialise must give the same deck.
namespace {

std::string SavedWithFixedVersion(UserData& data) {
  data.info.sim_version = 0;  // SaveToString increments it
  std::string out;
  SaveToString(data, out);
  return out;
}

void ExpectInitialiseChangesNothing(UserData& data, const char* label) {
  std::string err;
  ASSERT_TRUE(data.EvaluateAllExpressions(err)) << label << ": " << err;
  const std::string before = SavedWithFixedVersion(data);
  ASSERT_NO_THROW(data.fluid.Initialise(data.ic)) << label;
  ASSERT_NO_THROW(data.bed.Initialise(data.adsorbent_directory, data.fluid)) << label;
  for (auto& [id, ads] : data.adsorbent_directory.GetAll())
    const_cast<Adsorbent&>(ads).adsorption_source_.Initialise(data.fluid);
  const std::string after = SavedWithFixedVersion(data);
  if (before == after) return;
  // Report only the lines Initialise added or removed, not both decks.
  auto lines = [](const std::string& text) {
    std::set<std::string> out;
    std::istringstream in(text);
    for (std::string l; std::getline(in, l);) out.insert(l);
    return out;
  };
  const std::set<std::string> a = lines(before), b = lines(after);
  std::string diff;
  for (const auto& l : b) if (!a.count(l)) diff += "  + " + l + "\n";
  for (const auto& l : a) if (!b.count(l)) diff += "  - " + l + "\n";
  ADD_FAILURE() << label << ": Initialise changed the saved deck (a cascade left an entry for run time to seed):\n" << diff;
}

}  // namespace

TEST(DeckEdit, CascadesLeaveNothingForInitialise) {
  // A loaded deck edited through the cascades.
  UserData edited = LoadPerez();
  std::string err;
  ASSERT_TRUE(AddComponentToDeck(edited, "Argon", err)) << err;
  const AdsorbentID fresh = edited.adsorbent_directory.Add("Fresh");
  edited.adsorbent_directory.Get(fresh).adsorption_source_.Initialise(edited.fluid);
  edited.bed.layers_.Add("Fresh layer");
  SeedEmptyLayerCompositions(edited);
  ExpectInitialiseChangesNothing(edited, "Perez2019_A + component + adsorbent + layer");

  // A deck built from nothing (the builder goes through the same cascades).
  UserData built;
  build::AddComponent(built, "CarbonDioxide");
  build::AddComponent(built, "Nitrogen");
  const AdsorbentID zeolite = build::AddAdsorbent(built, "13X");
  build::SetIsotherm(built, zeolite, "extended_dual_site_langmuir_c", {"CarbonDioxide", "Nitrogen"});
  const size_t layer = build::AddLayer(built, "13X");
  build::SetLayerComposition(built, layer, {{zeolite, 1.0}});
  const size_t feed = build::AddStep(built, "Feed");
  build::SetFlowBoundary(built, feed, "bottom", "inlet");
  build::AddStream(built, feed, "bottom", "Feed inlet", "velocity", "constant");
  build::SetPressureBoundary(built, feed, "top", "outlet", "exponential");
  ExpectInitialiseChangesNothing(built, "built from scratch");
}

}  // namespace adsim
