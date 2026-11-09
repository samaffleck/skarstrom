/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// The deck view (objects/deck_view.h) and the display manifests
// (objects/manifests.h): every example deck's view validates against the
// engine's schema, lists every registry entry exactly once, and setting every
// value back by path gives the same deck; the manifests agree with the
// registry, the KPI list, the warning settings and the example decks.

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <cereal/archives/json.hpp>
#include <cereal/external/rapidjson/document.h>
#include <cereal/external/rapidjson/schema.h>
#include <cereal/external/rapidjson/stringbuffer.h>
#include <cereal/external/rapidjson/writer.h>
#include <gtest/gtest.h>

#include "adsim/math/unit_parser.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/deck_view.h"
#include "adsim/objects/manifests.h"
#include "adsim/objects/model_kinds.h"
#include "adsim/objects/user_data.h"
#include "adsim/parameters/parameter_registry.h"
#include "adsim/solver/operating_limits.h"
#include "adsim/solver/run_health.h"
#include "adsim/solver/sensitivity.h"

namespace adsim {
namespace {

namespace fs = std::filesystem;

const fs::path kExamplesRoot = fs::current_path() / ".." / "examples";
const fs::path kDecks = kExamplesRoot / "decks";

std::vector<fs::path> ExampleDecks() {
  std::vector<fs::path> decks;
  for (const auto& e : fs::directory_iterator(kDecks))
    if (e.path().extension() == ".json") decks.push_back(e.path());
  std::sort(decks.begin(), decks.end());
  return decks;
}

std::string Snapshot(UserData& data) {
  std::ostringstream oss;
  {
    cereal::JSONOutputArchive ar(oss);
    ar(data);
  }
  return oss.str();
}

void Parse(rapidjson::Document& doc, const std::string& text) {
  doc.Parse(text.c_str());
  ASSERT_FALSE(doc.HasParseError()) << "offset " << doc.GetErrorOffset();
}

std::string Write(const rapidjson::Value& v) {
  rapidjson::StringBuffer b;
  rapidjson::Writer<rapidjson::StringBuffer> w(b);
  v.Accept(w);
  return b.GetString();
}

// The first lines where two snapshots differ, with the line before each.
std::string LineDiff(const std::string& a, const std::string& b) {
  std::istringstream sa(a), sb(b);
  std::string la, lb, prev, out;
  int shown = 0;
  while (shown < 8) {
    const bool ga = static_cast<bool>(std::getline(sa, la));
    const bool gb = static_cast<bool>(std::getline(sb, lb));
    if (!ga && !gb) break;
    if (la != lb) {
      out += "  after: " + prev + "\n  - " + la + "\n  + " + lb + "\n";
      ++shown;
    }
    prev = la;
  }
  return out;
}

// Every value object of the view (an object with "path" and "kind").
void CollectValues(const rapidjson::Value& v, std::vector<const rapidjson::Value*>& out) {
  if (v.IsObject()) {
    if (v.HasMember("path") && v.HasMember("kind") && v.HasMember("unit")) {
      out.push_back(&v);
      return;
    }
    for (auto it = v.MemberBegin(); it != v.MemberEnd(); ++it) CollectValues(it->value, out);
  } else if (v.IsArray()) {
    for (const auto& e : v.GetArray()) CollectValues(e, out);
  }
}

class DeckViewTest : public ::testing::Test {
protected:
  static void SetUpTestSuite() {
    rapidjson::Document schema_doc;
    schema_doc.Parse(DeckViewSchema().c_str());
    ASSERT_FALSE(schema_doc.HasParseError());
    schema_ = new rapidjson::SchemaDocument(schema_doc);
  }
  static void TearDownTestSuite() {
    delete schema_;
    schema_ = nullptr;
  }

  // Validation errors of `view` against the schema ("" when valid).
  static std::string Validate(const rapidjson::Document& view) {
    rapidjson::SchemaValidator validator(*schema_);
    if (view.Accept(validator)) return {};
    rapidjson::StringBuffer where, keyword;
    validator.GetInvalidSchemaPointer().StringifyUriFragment(where);
    validator.GetInvalidDocumentPointer().StringifyUriFragment(keyword);
    return std::string("schema ") + where.GetString() + " (" + validator.GetInvalidSchemaKeyword() + ") at " +
           keyword.GetString();
  }

  static rapidjson::SchemaDocument* schema_;
};

rapidjson::SchemaDocument* DeckViewTest::schema_ = nullptr;

TEST_F(DeckViewTest, EmptyDeckViewValidates) {
  UserData data;
  rapidjson::Document view;
  Parse(view, DeckView(data));
  EXPECT_EQ(Validate(view), "");
  EXPECT_EQ(view["view_version"].GetInt(), kDeckViewVersion);
}

TEST_F(DeckViewTest, SchemaEnumeratesTheBuilderKinds) {
  rapidjson::Document schema;
  Parse(schema, DeckViewSchema());
  const auto& iso = schema["definitions"]["AdsorbentComponent"]["properties"]["isotherm"]["enum"];
  std::set<std::string> kinds;
  for (const auto& k : iso.GetArray()) kinds.insert(k.GetString());
  for (const auto& k : build::IsothermKinds()) EXPECT_TRUE(kinds.count(k)) << k;
  const auto& shape = schema["definitions"]["Bed"]["properties"]["shape"]["enum"];
  EXPECT_EQ(shape.Size(), build::BedShapeKinds().size());
}

// Every example deck: the view validates, lists every registry entry exactly
// once (nothing under "other"), shows enumerated values with their choice
// names, and setting every value back by path leaves the deck unchanged.
TEST_F(DeckViewTest, ExampleDecksValidateAndRoundTrip) {
  const auto decks = ExampleDecks();
  ASSERT_GE(decks.size(), 30u);
  for (const auto& deck : decks) {
    SCOPED_TRACE(deck.filename().string());
    // One loaded deck throughout: connection pair ids are numbered at load.
    UserData data;
    ASSERT_NO_THROW(LoadFromJson(data, deck.string()));
    data.EvaluateAllExpressions();
    const std::string before = Snapshot(data);
    rapidjson::Document view;
    Parse(view, DeckView(data));
    EXPECT_EQ(Validate(view), "");
    EXPECT_EQ(view["other"]["values"].Size(), 0u) << Write(view["other"]);
    // The cost estimate and the study settings have their nodes; each KPI
    // entry claims its tracked flag.
    EXPECT_GT(view["cost"]["values"].Size(), 0u);
    EXPECT_GT(view["sensitivity"]["values"].Size(), 0u);
    for (const auto& k : view["kpis"].GetArray()) EXPECT_EQ(k["values"].Size(), 1u) << Write(k);

    std::vector<const rapidjson::Value*> values;
    CollectValues(view, values);
    ParameterRegistry registry(data);
    std::map<std::string, int> seen;
    for (const auto* v : values) ++seen[(*v)["path"].GetString()];
    EXPECT_EQ(values.size(), registry.Size());
    for (const auto& p : registry.All()) EXPECT_EQ(seen[p.Path()], 1) << p.Path();

    for (const auto* v : values) {
      const auto& value = *v;
      const std::string path = value["path"].GetString();
      if (value["choices"].Size() > 0) {
        EXPECT_TRUE(value["choice"].IsString()) << path;
      }
      Parameter* p = registry.Find(path);
      ASSERT_NE(p, nullptr) << path;
      const std::string kind = value["kind"].GetString();
      if (kind == "expression" || kind == "string") {
        p->SetText(value["text"].GetString());
      } else if (kind == "boolean") {
        p->SetValue(value["value"].GetBool() ? 1.0 : 0.0);
      } else if (value["value"].IsNumber()) {
        p->SetValue(value["value"].GetDouble());
      }
    }
    data.EvaluateAllExpressions();
    const std::string after = Snapshot(data);
    EXPECT_TRUE(after == before) << LineDiff(before, after);
  }
}

TEST_F(DeckViewTest, ViewNodesClaimTheirOwnValues) {
  UserData data;
  ASSERT_NO_THROW(LoadFromJson(data, (kDecks / "Perez2019_A.json").string()));
  data.EvaluateAllExpressions();
  rapidjson::Document view;
  Parse(view, DeckView(data));
  ASSERT_GT(view["steps"].Size(), 0u);
  for (const auto& step : view["steps"].GetArray()) {
    const std::string sp = "step[" + std::to_string(step["id"].GetInt()) + "].";
    for (const auto& v : step["values"].GetArray())
      EXPECT_EQ(std::string(v["path"].GetString()).rfind(sp, 0), 0u) << v["path"].GetString();
    for (const char* end : {"bottom", "top"})
      for (const auto& s : step[end]["streams"].GetArray()) {
        const std::string ep = "endpoint[" + std::to_string(s["endpoint"].GetInt()) + "].";
        for (const auto& v : s["values"].GetArray())
          EXPECT_EQ(std::string(v["path"].GetString()).rfind(ep, 0), 0u) << v["path"].GetString();
      }
  }
  // Enumerated layer values carry their names.
  bool found = false;
  for (const auto& v : view["layers"][0]["values"].GetArray()) {
    if (std::string(v["key"].GetString()) != "momentum_model") continue;
    found = true;
    EXPECT_TRUE(v["choice"].IsString());
    EXPECT_GT(v["choices"].Size(), 1u);
  }
  EXPECT_TRUE(found);
}

// A stream's view type is its tag in the stream kind table, the tag that
// set_stream_kind takes, and an end's kind is its tag in the boundary kind
// table (set_boundary_kind). A front end can then show them with the kind lists
// and need not translate (the view once said "stream" for "flow").
TEST_F(DeckViewTest, StreamTypesAndEndKindsAreKindTags) {
  const auto& tags = KindTags(ModelFamily::Stream);
  const auto& end_tags = KindTags(ModelFamily::Boundary);
  rapidjson::Document schema;
  Parse(schema, DeckViewSchema());
  std::vector<std::string> schema_end_kinds;
  for (const auto& t : schema["definitions"]["End"]["properties"]["kind"]["enum"].GetArray())
    schema_end_kinds.push_back(t.GetString());
  EXPECT_EQ(schema_end_kinds, end_tags);
  std::vector<std::string> schema_types;
  for (const auto& t : schema["definitions"]["Stream"]["properties"]["type"]["enum"].GetArray())
    schema_types.push_back(t.GetString());
  std::vector<std::string> expected(tags.begin(), tags.end());
  expected.push_back("unsupported");
  EXPECT_EQ(schema_types, expected);

  std::set<std::string> seen;
  for (const auto& deck : ExampleDecks()) {
    UserData data;
    ASSERT_NO_THROW(LoadFromJson(data, deck.string())) << deck;
    rapidjson::Document view;
    Parse(view, DeckView(data));
    for (const auto& step : view["steps"].GetArray())
      for (const char* end : {"bottom", "top"}) {
        const std::string end_kind = step[end]["kind"].GetString();
        EXPECT_NE(std::find(end_tags.begin(), end_tags.end(), end_kind), end_tags.end()) << deck << ": " << end_kind;
        for (const auto& st : step[end]["streams"].GetArray()) {
          const std::string type = st["type"].GetString();
          EXPECT_NE(std::find(tags.begin(), tags.end(), type), tags.end()) << deck << ": " << type;
          seen.insert(type);
        }
      }
  }
  EXPECT_TRUE(seen.count("flow")) << "no example deck has a set-flow stream";
}

// ---- manifests ------------------------------------------------------------------

TEST(Manifests, KpisFollowTheEngineKpiList) {
  const auto& names = sensitivity::KPINames();
  std::vector<std::string> scalars;
  for (const auto& k : KpiDefinitions()) {
    EXPECT_FALSE(k.label.empty()) << k.key;
    EXPECT_FALSE(k.description.empty()) << k.key;
    EXPECT_TRUE(CheckUnit(k.unit).ok) << k.key << " unit " << k.unit;
    if (!k.per_step && !k.diagnostic && !k.reaction) scalars.push_back(k.key);
    EXPECT_FALSE(k.reaction && k.diagnostic) << k.key;
  }
  EXPECT_EQ(scalars, names);
  // Reported KPIs first, then the diagnostics.
  const auto& defs = KpiDefinitions();
  bool seen_diagnostic = false;
  for (const auto& k : defs) {
    if (k.diagnostic) seen_diagnostic = true;
    else EXPECT_FALSE(seen_diagnostic) << k.key << " follows a diagnostic";
  }
  EXPECT_TRUE(seen_diagnostic);
  rapidjson::Document doc;
  Parse(doc, KpiManifest());
  EXPECT_EQ(doc.Size(), KpiDefinitions().size());
}

TEST(Manifests, WarningsMatchTheSettingsAndOperatingLimits) {
  std::set<std::string> keys, settings;
  for (const auto& w : WarningDefinitions()) {
    EXPECT_TRUE(keys.insert(w.key).second) << w.key;
    EXPECT_TRUE(w.phase == "start" || w.phase == "end") << w.key;
    EXPECT_FALSE(w.message.empty()) << w.key;
    if (!w.setting.empty()) {
      settings.insert(w.setting);
      EXPECT_EQ(w.setting, "settings.warn_" + w.key);
    }
  }
  for (const char* limit : kOperatingLimitList) EXPECT_TRUE(keys.count(limit)) << limit;
  // Every run-health code has its entry, keyed by the code in lower case.
  for (const char* code : kRunHealthCodeList) {
    const auto it = std::find_if(WarningDefinitions().begin(), WarningDefinitions().end(),
                                 [&](const WarningDefinition& w) { return w.code == code; });
    ASSERT_NE(it, WarningDefinitions().end()) << code;
    std::string key = code;
    for (auto& ch : key) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    EXPECT_EQ(it->key, key);
    EXPECT_FALSE(it->setting.empty()) << code;
  }
  // Every warning switch of the registry has its manifest entry, and back.
  std::set<std::string> registry_switches;
  for (const auto* f : StaticFieldsInSection("settings")) {
    const std::string path = std::string("settings.") + f->key;
    if (path.rfind("settings.warn_", 0) == 0) registry_switches.insert(path);
  }
  EXPECT_EQ(settings, registry_switches);
  rapidjson::Document doc;
  Parse(doc, WarningManifest());
  EXPECT_EQ(doc.Size(), WarningDefinitions().size());
}

TEST(Manifests, SuggestedUnitsParseIntoTheirCategory) {
  rapidjson::Document doc;
  Parse(doc, UnitManifest());
  ASSERT_TRUE(doc["categories"].IsArray());
  EXPECT_EQ(doc["categories"].Size(), build::VariableCategoryKinds().size());
  size_t units = 0;
  for (const auto& c : doc["categories"].GetArray()) {
    const std::string key = c["key"].GetString();
    for (const auto& u : c["units"].GetArray()) {
      ++units;
      const std::string symbol = u["symbol"].GetString();
      const UnitCheck check = CheckUnit(symbol, key);
      EXPECT_TRUE(check.ok) << key << ": " << symbol << " " << check.message;
      EXPECT_TRUE(u["factor"].IsNumber()) << symbol;
    }
    // The SI unit of the category is a suggested unit too.
    const std::string si = c["si_unit"].GetString();
    if (!si.empty()) EXPECT_TRUE(CheckUnit(si, key).ok) << key << " SI unit " << si;
  }
  EXPECT_GT(units, 60u);
  for (const auto& o : doc["overrides"].GetArray()) EXPECT_TRUE(CheckUnit(o["symbol"].GetString()).ok);
}

TEST(Manifests, CheckUnitAcceptsWhatTheDecksUseAndRejectsNonsense) {
  // Every unit an example deck writes in a value ("2 [min]", "120 s").
  std::set<std::string> units;
  const std::regex bracket(R"(\[([^\]]+)\]\s*$)");
  for (const auto& deck : ExampleDecks()) {
    UserData data;
    ASSERT_NO_THROW(LoadFromJson(data, deck.string()));
    ParameterRegistry registry(data);
    for (const auto& p : registry.All()) {
      if (p.Info().kind != ParameterKind::Expression) continue;
      const std::string text = p.GetText();
      std::smatch m;
      double number = 0.0;
      std::string unit;
      if (std::regex_search(text, m, bracket)) units.insert(m[1].str());
      else if (SplitQuantity(text, number, unit)) units.insert(unit);
    }
  }
  for (const auto& u : units) EXPECT_TRUE(CheckUnit(u).ok) << u << ": " << CheckUnit(u).message;

  const UnitCheck kpa = CheckUnit("kPa", "pressure");
  EXPECT_TRUE(kpa.ok);
  EXPECT_EQ(kpa.category, "pressure");
  EXPECT_DOUBLE_EQ(kpa.factor, 1000.0);
  const UnitCheck degc = CheckUnit("degC");
  EXPECT_TRUE(degc.ok);
  EXPECT_DOUBLE_EQ(degc.offset, 273.15);
  EXPECT_TRUE(CheckUnit("%", "dimensionless").ok);

  const UnitCheck nonsense = CheckUnit("blargs/fortnightly");
  EXPECT_FALSE(nonsense.ok);
  EXPECT_FALSE(nonsense.message.empty());
  const UnitCheck wrong = CheckUnit("bar", "time");
  EXPECT_FALSE(wrong.ok);
  EXPECT_NE(wrong.message.find("pressure"), std::string::npos) << wrong.message;
}

TEST(Manifests, ExampleDeckManifestListsEveryDeck) {
  const std::string text = ExampleDeckManifest((kExamplesRoot / "example_decks.json").string(), kDecks.string());
  rapidjson::Document doc;
  Parse(doc, text);
  EXPECT_EQ(doc.Size(), ExampleDecks().size());
  for (const auto& d : doc.GetArray()) {
    EXPECT_FALSE(std::string(d["family"].GetString()).empty());
    EXPECT_FALSE(std::string(d["title"].GetString()).empty());
    EXPECT_FALSE(std::string(d["summary"].GetString()).empty()) << d["file"].GetString();
    EXPECT_TRUE(d["hidden"].IsBool());
  }
}

TEST(Manifests, ComponentsHaveFormulaAndMolarMass) {
  rapidjson::Document root;
  Parse(root, ComponentManifest());
  EXPECT_FALSE(std::string(root["cp_reference"].GetString()).empty());
  const auto& doc = root["components"];
  EXPECT_EQ(doc.Size(), Components::componentNames.size());
  for (const auto& c : doc.GetArray()) {
    const std::string name = c["name"].GetString();
    // A component with a Cp polynomial names its species; a held one has none.
    EXPECT_EQ(std::string(c["cp_species"].GetString()).empty(), !c["ideal_gas_heat_capacity"].GetBool()) << name;
    EXPECT_TRUE(c["molar_mass"].IsNumber()) << name;
    if (name != "Air") EXPECT_FALSE(std::string(c["cas"].GetString()).empty()) << name;
  }
  const ComponentInfo co2 = DescribeComponent("CarbonDioxide");
  EXPECT_EQ(co2.formula, "CO2");
  EXPECT_EQ(co2.cas, "124-38-9");
  EXPECT_NEAR(co2.molar_mass, 0.04401, 1e-4);
  EXPECT_TRUE(co2.ideal_gas_heat_capacity);
}

}  // namespace
}  // namespace adsim
