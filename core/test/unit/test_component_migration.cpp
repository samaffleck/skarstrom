/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// Air was removed as a component; decks, adsorbents and cycles that have it
// load with Nitrogen instead (component_migration.h).

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include <cereal/archives/json.hpp>
#include <cereal/external/rapidjson/document.h>
#include <gtest/gtest.h>

#include "adsim/objects/component_migration.h"
#include "adsim/objects/component_properties.h"
#include "adsim/objects/deck_io.h"
#include "adsim/objects/fluid.h"
#include "adsim/objects/user_data.h"
#include "adsim/simulation.h"

namespace adsim {
namespace {

namespace rj = CEREAL_RAPIDJSON_NAMESPACE;

const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";
const std::filesystem::path kLegacy =
    std::filesystem::current_path() / "test" / "fixtures" / "legacy" / "Ward_2022_Q=50[ccm]_pre_physics_parity.json";

std::string ReadFile(const std::filesystem::path& p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

std::string ReplaceAll(std::string s, const std::string& from, const std::string& to) {
  for (size_t pos = s.find(from); pos != std::string::npos; pos = s.find(from, pos + to.size()))
    s.replace(pos, from.size(), to);
  return s;
}

std::string Snapshot(UserData& data) {
  std::ostringstream oss;
  {
    cereal::JSONOutputArchive ar(oss);
    ar(data);
  }
  return oss.str();
}

rj::Document Parse(const std::string& json) {
  rj::Document d;
  d.Parse(json.c_str());
  EXPECT_FALSE(d.HasParseError());
  return d;
}

// The value of `key` in a cereal map array, or nullptr.
const rj::Value* MapValue(const rj::Value& map, const char* key) {
  for (const auto& e : map.GetArray())
    if (std::string(e["key"].GetString()) == key) return &e["value"];
  return nullptr;
}

TEST(ComponentMigration, AirIsNotAComponent) {
  EXPECT_EQ(std::find(Components::componentNames.begin(), Components::componentNames.end(), "Air"),
            Components::componentNames.end());
  EXPECT_EQ(FindComponentProperties("Air"), nullptr);
  EXPECT_EQ(Fluid::IdealGasCoeffs("Air"), nullptr);
  EXPECT_EQ(Components::renamedComponents.at("Air"), "Nitrogen");
}

// A deck with Air and no Nitrogen: every component reference is renamed;
// stream names are left alone.
TEST(ComponentMigration, RenamesAirToNitrogen) {
  std::string json = R"({"value0": {
    "Fluid": {"Components": ["CarbonDioxide", "Air"],
              "Adsorbed_Heat_Capacity": [{"key": "Air", "value": {"Value": 12.5, "Expression": "12.5"}}]},
    "Initial_Conditions": {"Initial_Composition": [{"key": "CarbonDioxide", "value": 0.0}, {"key": "Air", "value": 1.0}]},
    "Stream": {"Stream_Name": "Air", "Endpoint": "Air"},
    "Sweep": {"Path": "adsorbent[1].isotherm[Air].b"},
    "KPI": {"Product_Component": "Air", "Products": ["Air"]}
  }})";
  ASSERT_TRUE(MigrateRenamedComponents(json));
  const rj::Document d = Parse(json);
  const rj::Value& v = d["value0"];
  ASSERT_EQ(v["Fluid"]["Components"].Size(), 2u);
  EXPECT_STREQ(v["Fluid"]["Components"][1].GetString(), "Nitrogen");
  const rj::Value* cpa = MapValue(v["Fluid"]["Adsorbed_Heat_Capacity"], "Nitrogen");
  ASSERT_NE(cpa, nullptr);
  EXPECT_DOUBLE_EQ((*cpa)["Value"].GetDouble(), 12.5);
  EXPECT_EQ(MapValue(v["Fluid"]["Adsorbed_Heat_Capacity"], "Air"), nullptr);
  const rj::Value* y = MapValue(v["Initial_Conditions"]["Initial_Composition"], "Nitrogen");
  ASSERT_NE(y, nullptr);
  EXPECT_DOUBLE_EQ(y->GetDouble(), 1.0);
  EXPECT_STREQ(v["Stream"]["Stream_Name"].GetString(), "Air");
  EXPECT_STREQ(v["Stream"]["Endpoint"].GetString(), "Air");
  EXPECT_STREQ(v["Sweep"]["Path"].GetString(), "adsorbent[1].isotherm[Nitrogen].b");
  EXPECT_STREQ(v["KPI"]["Product_Component"].GetString(), "Nitrogen");
  EXPECT_STREQ(v["KPI"]["Products"][0].GetString(), "Nitrogen");
}

// A deck with both: Air's mole fractions are added to Nitrogen's, and its
// other entries are dropped so that Nitrogen keeps its own physics.
TEST(ComponentMigration, MergesAirIntoAnExistingNitrogen) {
  std::string json = R"({"value0": {
    "Fluid": {"Components": ["Nitrogen", "Air"],
              "Adsorbed_Heat_Capacity": [{"key": "Air", "value": {"Value": 1.0, "Expression": "1"}},
                                         {"key": "Nitrogen", "value": {"Value": 2.0, "Expression": "2"}}]},
    "Initial_Conditions": {"Initial_Composition": [{"key": "Air", "value": 0.25}, {"key": "Nitrogen", "value": 0.5}]},
    "Step": {"Inlet_Mole_Fraction": [{"key": "Air", "value": {"Value": 0.2, "Expression": "0.2"}},
                                     {"key": "Nitrogen", "value": {"Value": 0.1, "Expression": "x"}}],
             "Boundary_Mole_Fraction": [{"key": "Air", "value": {"Value": 0.25, "Expression": "0.25"}},
                                        {"key": "Nitrogen", "value": {"Value": 0.5, "Expression": "0.5"}}]},
    "Adsorbent": {"Isotherm_Map": [{"key": "Air", "value": 0}, {"key": "Nitrogen", "value": 3}],
                  "Isotherm_Params": [{"key": "Air", "value": {"q": 1}}]}
  }})";
  ASSERT_TRUE(MigrateRenamedComponents(json));
  const rj::Document d = Parse(json);
  const rj::Value& v = d["value0"];
  ASSERT_EQ(v["Fluid"]["Components"].Size(), 1u);
  EXPECT_STREQ(v["Fluid"]["Components"][0].GetString(), "Nitrogen");
  ASSERT_EQ(v["Fluid"]["Adsorbed_Heat_Capacity"].Size(), 1u);
  EXPECT_DOUBLE_EQ((*MapValue(v["Fluid"]["Adsorbed_Heat_Capacity"], "Nitrogen"))["Value"].GetDouble(), 2.0);
  ASSERT_EQ(v["Initial_Conditions"]["Initial_Composition"].Size(), 1u);
  EXPECT_DOUBLE_EQ(MapValue(v["Initial_Conditions"]["Initial_Composition"], "Nitrogen")->GetDouble(), 0.75);
  const rj::Value& inlet = *MapValue(v["Step"]["Inlet_Mole_Fraction"], "Nitrogen");
  EXPECT_NEAR(inlet["Value"].GetDouble(), 0.3, 1e-15);
  EXPECT_STREQ(inlet["Expression"].GetString(), "(x) + (0.2)");
  const rj::Value& boundary = *MapValue(v["Step"]["Boundary_Mole_Fraction"], "Nitrogen");
  EXPECT_STREQ(boundary["Expression"].GetString(), "0.75");
  ASSERT_EQ(v["Adsorbent"]["Isotherm_Map"].Size(), 1u);
  EXPECT_EQ(MapValue(v["Adsorbent"]["Isotherm_Map"], "Nitrogen")->GetInt(), 3);
  EXPECT_EQ(v["Adsorbent"]["Isotherm_Params"].Size(), 0u);
}

// A deck without Air is left as it is, even when a stream is named "Air".
TEST(ComponentMigration, LeavesOtherDecksUnchanged) {
  const std::string original =
      R"({"value0": {"Fluid": {"Components": ["CarbonDioxide", "Nitrogen"]}, "Stream": {"Stream_Name": "Air"}}})";
  std::string json = original;
  EXPECT_FALSE(MigrateRenamedComponents(json));
  EXPECT_EQ(json, original);

  std::string not_json = R"({"Air": )";
  EXPECT_FALSE(MigrateRenamedComponents(not_json));
  EXPECT_EQ(not_json, R"({"Air": )");
}

// A shipped CO2/N2 deck saved with Air in place of Nitrogen loads, through
// every load function, to the same deck as the original.
TEST(ComponentMigration, OldAirDeckLoadsAsTheNitrogenDeck) {
  const std::string n2_json = ReadFile(kExamples / "Ward2024.json");
  ASSERT_FALSE(n2_json.empty());
  const std::string air_json = ReplaceAll(n2_json, "\"Nitrogen\"", "\"Air\"");
  ASSERT_NE(air_json, n2_json);

  UserData expected;
  LoadDeckFlexible(expected, n2_json);
  const std::string want = Snapshot(expected);

  UserData flexible;
  LoadDeckFlexible(flexible, air_json);
  EXPECT_EQ(Snapshot(flexible), want);

  const std::filesystem::path tmp = std::filesystem::temp_directory_path() / "skarstrom_air_deck.json";
  {
    std::ofstream out(tmp, std::ios::binary);
    out << air_json;
  }
  UserData from_file;
  LoadFromJson(from_file, tmp.string());
  std::filesystem::remove(tmp);
  EXPECT_EQ(Snapshot(from_file), want);

  Simulation sim = Simulation::FromJson(air_json);
  const auto& components = sim.Data().fluid.components_;
  EXPECT_NE(std::find(components.begin(), components.end(), "Nitrogen"), components.end());
  EXPECT_EQ(std::find(components.begin(), components.end(), "Air"), components.end());
}

// The legacy fixture (a CO2/Air deck from before the physics-parity options)
// loads with Nitrogen and passes the deck checks.
TEST(ComponentMigration, LegacyAirFixtureLoadsWithNitrogen) {
  UserData data;
  LoadFromJson(data, kLegacy.string());
  EXPECT_EQ(data.fluid.components_, (std::vector<std::string>{"CarbonDioxide", "Nitrogen"}));
  for (const auto& [id, ads] : data.adsorbent_directory.GetAll()) {
    EXPECT_EQ(ads.adsorption_source_.IsothermMap().count("Air"), 0u);
    EXPECT_EQ(ads.adsorption_source_.IsothermMap().count("Nitrogen"), 1u);
  }
  EXPECT_DOUBLE_EQ(static_cast<double>(data.ic.mole_fraction.at("Nitrogen")), 1.0);
  data.EvaluateAllExpressions();
  std::string err;
  EXPECT_TRUE(data.PerformChecks(err)) << err;
}

}  // namespace
}  // namespace adsim
