// E-17: removing an adsorbent must not leave dangling ids in layer compositions.
//
// Layer compositions reference adsorbents by AdsorbentID, and every directory
// lookup (AdsorbentDirectory::Get is unordered_map::at) throws std::out_of_range
// on an id that no longer exists. UserData::RemoveAdsorbent cascades the removal
// into every layer, and Bed::PerformChecks rejects any composition id missing
// from the directory with a clear message instead of letting Initialise throw.

// STL includes
#include <filesystem>
#include <string>

// AdSim includes
#include "adsim/objects/user_data.h"

// GTest includes
#include <gtest/gtest.h>


namespace adsim {
namespace {

  const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";

  TEST(AdsorbentRemovalTest, RemoveAdsorbentErasesIdFromEveryLayer) {
    UserData ud;
    const AdsorbentID id_a = ud.adsorbent_directory.Add("A");
    const AdsorbentID id_b = ud.adsorbent_directory.Add("B");

    auto& mixed = ud.bed.layers_.Get(ud.bed.layers_.Add("Mixed"));
    mixed.AddAdsorbent(id_a, 0.5);
    mixed.AddAdsorbent(id_b, 0.5);
    auto& pure_a = ud.bed.layers_.Get(ud.bed.layers_.Add("Pure A"));
    pure_a.AddAdsorbent(id_a, 1.0);
    ud.cost_estimator.adsorbent_cost_per_kg[id_a] = 10.0;
    ud.cost_estimator.adsorbent_cost_per_kg[id_b] = 20.0;

    ud.RemoveAdsorbent(id_a);

    EXPECT_FALSE(ud.adsorbent_directory.Contains(id_a));
    EXPECT_TRUE(ud.adsorbent_directory.Contains(id_b));
    for (const auto& layer : ud.bed.layers_.layers_) {
      EXPECT_EQ(layer.composition_.count(id_a), 0u) << layer.name_;
      for (const auto& [id, frac] : layer.composition_)
        EXPECT_TRUE(ud.adsorbent_directory.Contains(id)) << layer.name_ << " id " << id;
    }
    // Other entries are untouched.
    const auto& mixed_after = ud.bed.layers_.layers_[0].composition_;
    ASSERT_EQ(mixed_after.size(), 1u);
    EXPECT_DOUBLE_EQ(mixed_after.at(id_b), 0.5);
    EXPECT_EQ(ud.cost_estimator.adsorbent_cost_per_kg.count(id_a), 0u);
    EXPECT_DOUBLE_EQ(ud.cost_estimator.adsorbent_cost_per_kg.at(id_b), 20.0);

    // The layer left empty is reported as such, not as a dangling id.
    std::string err;
    EXPECT_FALSE(ud.bed.PerformChecks(ud.adsorbent_directory, err));
    EXPECT_NE(err.find("\"Pure A\" has no adsorbent selected"), std::string::npos) << err;

    // Removing an unknown id is a no-op.
    EXPECT_NO_THROW(ud.RemoveAdsorbent(id_a));
    EXPECT_NO_THROW(ud.RemoveAdsorbent(9999));
  }

  TEST(AdsorbentRemovalTest, BedChecksRejectDanglingCompositionId) {
    AdsorbentDirectory dir;
    const AdsorbentID id_a = dir.Add("A");
    Bed bed;
    auto& layer = bed.layers_.Get(bed.layers_.Add("Bottom"));
    layer.AddAdsorbent(id_a, 1.0);

    std::string err;
    EXPECT_TRUE(bed.PerformChecks(dir, err)) << err;

    dir.Remove(id_a); // bypasses the cascade: the layer now dangles
    err.clear();
    EXPECT_FALSE(bed.PerformChecks(dir, err));
    EXPECT_NE(err.find("Layer \"Bottom\" references adsorbent id " + std::to_string(id_a)), std::string::npos) << err;
  }

  // A full deck whose layer composition references a missing adsorbent fails
  // UserData::PerformChecks with a clear message instead of Initialise throwing
  // std::out_of_range from the directory lookup.
  TEST(AdsorbentRemovalTest, DeckWithDanglingIdFailsChecksInsteadOfThrowing) {
    UserData ud;
    LoadFromJson(ud, (kExamples / "FiveStepProductPurge.json").string());

    std::string err;
    ASSERT_TRUE(ud.PerformChecks(err)) << err;
    ASSERT_FALSE(ud.bed.layers_.layers_.empty());
    const auto& layer = ud.bed.layers_.layers_.front();
    ASSERT_FALSE(layer.composition_.empty());
    const AdsorbentID dangling = layer.composition_.begin()->first;

    ud.adsorbent_directory.Remove(dangling); // what the GUI used to do
    err.clear();
    bool ok = true;
    EXPECT_NO_THROW(ok = ud.PerformChecks(err));
    EXPECT_FALSE(ok);
    EXPECT_NE(err.find("Layer \"" + layer.name_ + "\" references adsorbent id " + std::to_string(dangling)),
              std::string::npos) << err;
  }

} // namespace
} // namespace adsim
