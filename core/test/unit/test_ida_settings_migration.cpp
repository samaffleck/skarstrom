// IDASettings v5 changed the default maximum step from 2 s to no cap (0).
// Decks saved before v5 that still carry the old default migrate on load;
// values the user chose are kept.

#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include <cereal/archives/json.hpp>

#include "adsim/solver/ida_settings.h"

namespace adsim {

  namespace {

    // Serialises default settings, then rewrites them to look like a v4 deck
    // whose Maximum_Step_Size expression is `old_expr`, and loads that back.
    IDASettings LoadAsVersion4(const std::string& old_expr) {
      std::ostringstream out;
      {
        cereal::JSONOutputArchive archive(out);
        IDASettings settings;
        archive(cereal::make_nvp("settings", settings));
      }
      std::string json = out.str();

      // The first class version in the document is IDASettings' own.
      const std::string v5 = "\"cereal_class_version\": 5";
      const auto version_pos = json.find(v5);
      EXPECT_NE(version_pos, std::string::npos);
      json.replace(version_pos, v5.size(), "\"cereal_class_version\": 4");

      // Rewrite the max-step expression (and its cached value) to the old one.
      const auto field = json.find("\"Maximum_Step_Size\"");
      EXPECT_NE(field, std::string::npos);
      const auto value_key = json.find("\"Value\":", field);
      const auto value_end = json.find(',', value_key);
      json.replace(value_key, value_end - value_key, "\"Value\": " + old_expr);
      const auto expr_key = json.find("\"Expression\":", field);
      const auto expr_end = json.find(',', expr_key);
      json.replace(expr_key, expr_end - expr_key, "\"Expression\": \"" + old_expr + "\"");

      std::istringstream in(json);
      cereal::JSONInputArchive archive(in);
      IDASettings loaded;
      archive(cereal::make_nvp("settings", loaded));
      return loaded;
    }

  } // namespace

  TEST(IDASettingsMigration, DefaultIsNoMaxStepCap) {
    IDASettings settings;
    EXPECT_EQ(settings.max_step.Get(), 0.0);
  }

  TEST(IDASettingsMigration, OldDefaultTwoSecondsBecomesNoCap) {
    IDASettings loaded = LoadAsVersion4("2");
    EXPECT_EQ(loaded.max_step.Get(), 0.0);
  }

  TEST(IDASettingsMigration, UserChosenMaxStepIsKept) {
    IDASettings loaded = LoadAsVersion4("0.5");
    EXPECT_EQ(loaded.max_step.Get(), 0.5);
    EXPECT_EQ(loaded.max_step.GetExpr(), "0.5");
  }

} // end adsim namespace
