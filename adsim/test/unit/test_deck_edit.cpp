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
#include <sstream>
#include <string>

#include <gtest/gtest.h>

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
