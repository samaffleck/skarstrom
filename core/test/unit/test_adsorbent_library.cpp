// Integrity gate for the curated adsorbent library (data/adsorbents/).
//
// Every catalog entry (library.json) must:
//   * point at a cereal blob that loads with the stock LoadAdsorbentFromJSON;
//   * carry a citation with a DOI (or be explicitly marked a legacy entry whose
//     source is unknown);
//   * describe exactly the blob's adsorbing components, with metadata values
//     (density, Cp, voidage, conductivity, heat of adsorption) equal to the blob;
//   * evaluate its isotherms to finite, non-negative loadings over the stated
//     fit range (and over a generic 273-373 K / 1 Pa-10 bar envelope);
//   * when it claims to be the adsorbent of an example deck (deck_checks), give
//     the same pure-component loadings as that deck's adsorbent.
// The docs page is generated from the same files by
// scripts/docs/gen_adsorbent_library.py, which also validates the catalogue
// (its --check mode fails on drift).

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <algorithm>

#include <gtest/gtest.h>
#include <cereal/external/rapidjson/document.h>

#include "adsim/objects/adsorbent_library.h"
#include "adsim/objects/porous_media.h"
#include "adsim/objects/user_data.h"
#include "adsim/equations/isotherm.h"
#include "adsim/equations/isotherm_kinds.h"  // EvalIsotherm (the evaluator front ends plot with)

namespace adsim {
namespace {

  const std::filesystem::path kLibDir = std::filesystem::current_path() / ".." / "data" / "adsorbents";
  const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";

  std::string ReadText(const std::filesystem::path& p) {
    std::ifstream ifs(p);
    std::stringstream ss;
    ss << ifs.rdbuf();
    return ss.str();
  }

  struct Catalog {
    AdsorbentLibrary lib;
    rapidjson::Document raw;  // for the fields the read-only view does not expose
    Catalog() {
      const std::string text = ReadText(kLibDir / "library.json");
      lib = ParseAdsorbentLibrary(text);
      raw.Parse(text.c_str());
    }
    const rapidjson::Value& RawEntry(const std::string& key) const {
      for (const auto& e : raw["entries"].GetArray())
        if (key == e["key"].GetString()) return e;
      throw std::runtime_error("no raw entry " + key);
    }
  };

  const Catalog& Cat() {
    static const Catalog c;
    return c;
  }

  bool RelEq(double a, double b) {
    return std::abs(a - b) <= 1e-9 * std::max({1.0, std::abs(a), std::abs(b)});
  }

  void LoadBlob(const AdsorbentLibraryEntry& e, Adsorbent& ads) {
    LoadAdsorbentFromJSON((kLibDir / e.file).string(), ads);
  }

  double Loading(Adsorbent& ads, const std::string& comp, double T, double P) {
    return EvalIsotherm(ads.adsorption_source_.GetIsotherm(comp).get(), comp, T, P, 1.0);
  }

}  // namespace

TEST(AdsorbentLibraryTest, CatalogParsesWithUniqueKeysAndFiles) {
  const auto& lib = Cat().lib;
  EXPECT_EQ(lib.schema_version, 3);
  ASSERT_GE(lib.entries.size(), 20u);
  std::set<std::string> keys, files;
  for (const auto& e : lib.entries) {
    EXPECT_TRUE(keys.insert(e.key).second) << "duplicate key " << e.key;
    EXPECT_TRUE(files.insert(e.file).second) << "duplicate file " << e.file;
    // Blob named after the key: a shell-safe asset name (spaces/parentheses
    // break the asset-copy build steps of front ends).
    EXPECT_EQ(e.file, e.key + ".json") << e.key;
    EXPECT_TRUE(std::filesystem::exists(kLibDir / e.file)) << e.key << ": missing " << e.file;
    EXPECT_FALSE(e.label.empty()) << e.key;
    EXPECT_FALSE(e.family.empty()) << e.key;
    EXPECT_FALSE(e.applications.empty()) << e.key;
  }
  // Keys the web build has always exposed stay available (as an entry, or as
  // the alias of a removed entry that resolves to its material).
  for (const char* k : {"activated_carbon", "mg_mof_74"})
    EXPECT_NE(lib.Resolve(k), nullptr) << k;
  // No orphan blobs: every adsorbent file in the folder is catalogued.
  for (const auto& f : std::filesystem::directory_iterator(kLibDir)) {
    const auto name = f.path().filename().string();
    if (f.path().extension() != ".json" || name == "library.json") continue;
    EXPECT_TRUE(files.count(name)) << "orphan adsorbent file " << name;
  }
}

// One recommended entry per material; a material's other entries are named
// variants; material ids, entry keys and aliases resolve without ambiguity.
TEST(AdsorbentLibraryTest, EveryMaterialHasOneRecommendedEntryAndNamedVariants) {
  const auto& lib = Cat().lib;
  ASSERT_FALSE(lib.materials.empty());
  const std::regex id_re(R"(^[a-z0-9]+(-[a-z0-9]+)*$)");
  std::set<std::string> ids, names;
  for (const auto& m : lib.materials) {
    EXPECT_TRUE(std::regex_match(m.id, id_re)) << m.id;
    EXPECT_TRUE(ids.insert(m.id).second) << "duplicate material " << m.id;
    EXPECT_TRUE(names.insert(m.name).second) << "duplicate material name " << m.name;
    EXPECT_EQ(lib.Find(m.id), nullptr) << m.id << " collides with an entry key";
    EXPECT_FALSE(m.recommended_reason.empty()) << m.id;
    const auto variants = lib.Variants(m.id);
    ASSERT_FALSE(variants.empty()) << m.id << " has no entries";
    int n_recommended = 0;
    std::set<std::string> variant_names;
    for (const AdsorbentLibraryEntry* e : variants) {
      n_recommended += e->recommended ? 1 : 0;
      EXPECT_EQ(e->family, m.family) << e->key;
      if (variants.size() == 1) {
        EXPECT_TRUE(e->variant.empty()) << e->key << ": the only entry of " << m.id << " has a variant";
        EXPECT_EQ(e->label, m.name) << e->key;
      } else {
        EXPECT_FALSE(e->variant.empty()) << e->key << ": " << m.id << " has several entries, so each needs a variant";
        EXPECT_TRUE(variant_names.insert(e->variant).second) << e->key << ": duplicate variant";
        EXPECT_EQ(e->label, m.name + " (" + e->variant + ")") << e->key;
      }
    }
    EXPECT_EQ(n_recommended, 1) << m.id;
    ASSERT_NE(variants.front(), nullptr);
    EXPECT_TRUE(variants.front()->recommended) << m.id << ": Variants() lists the recommended entry first";
    EXPECT_EQ(lib.Resolve(m.id), variants.front()) << m.id;
    for (const auto& a : m.aliases) {
      EXPECT_EQ(lib.Find(a), nullptr) << "alias " << a << " is still an entry key";
      EXPECT_EQ(lib.Resolve(a), variants.front()) << "alias " << a;
    }
  }
  for (const auto& e : lib.entries) {
    EXPECT_NE(lib.FindMaterial(e.material_id), nullptr) << e.key << ": unknown material " << e.material_id;
    EXPECT_EQ(lib.Resolve(e.key), &e) << e.key;
  }
  EXPECT_EQ(lib.Resolve("no-such-material"), nullptr);
}

// Each isotherm records the quality of its data; measured_tabulated goes with
// the tabulated source data stored in the entry.
TEST(AdsorbentLibraryTest, EveryIsothermRecordsItsDataQuality) {
  const std::set<std::string> levels{"measured_tabulated", "measured_fit", "refit", "secondary", "not_checked", "unknown"};
  for (const auto& e : Cat().lib.entries) {
    const auto& raw = Cat().RawEntry(e.key);
    ASSERT_TRUE(raw.HasMember("isotherms")) << e.key;
    for (const auto& iso : raw["isotherms"].GetArray()) {
      const std::string comp = iso["component"].GetString();
      ASSERT_TRUE(iso.HasMember("data_quality") && iso["data_quality"].IsString()) << e.key << " " << comp;
      const std::string dq = iso["data_quality"].GetString();
      EXPECT_TRUE(levels.count(dq)) << e.key << " " << comp << ": " << dq;
      EXPECT_EQ(dq == "unknown", e.IsLegacy()) << e.key << " " << comp;
      if (dq == "not_checked") EXPECT_TRUE(iso.HasMember("data_quality_note")) << e.key << " " << comp;
      const bool measured_points = iso.HasMember("source_points") &&
                                   std::string(iso["source_points"]["kind"].GetString()) == "measured";
      EXPECT_EQ(dq == "measured_tabulated", measured_points) << e.key << " " << comp;
    }
  }
}

TEST(AdsorbentLibraryTest, EveryEntryIsCitedOrMarkedLegacy) {
  const auto& lib = Cat().lib;
  const std::regex doi(R"(^10\.\d{4,9}/\S+$)");
  for (const auto& e : lib.entries) {
    if (e.IsLegacy()) {
      EXPECT_EQ(e.legacy_note.rfind("Source unknown", 0), 0u) << e.key << ": legacy note must start 'Source unknown'";
      continue;
    }
    ASSERT_EQ(e.status, "sourced") << e.key;
    ASSERT_FALSE(e.sources.empty()) << e.key << " has no source";
    for (const auto& r : e.sources) {
      const AdsorbentLibrarySource* s = lib.FindSource(r.source_id);
      ASSERT_NE(s, nullptr) << e.key << ": unknown source " << r.source_id;
      EXPECT_FALSE(s->citation.empty()) << r.source_id;
      EXPECT_TRUE(std::regex_match(s->doi, doi)) << r.source_id << ": bad DOI '" << s->doi << "'";
      EXPECT_FALSE(r.refs.empty()) << e.key << ": no table/figure reference for " << r.source_id;
    }
    EXPECT_NE(lib.CitationText(e).find("doi:10."), std::string::npos) << e.key;
  }
}

TEST(AdsorbentLibraryTest, EveryEntryLoadsAndMatchesItsMetadata) {
  const auto& cat = Cat();
  for (const auto& e : cat.lib.entries) {
    SCOPED_TRACE(e.key);
    Adsorbent ads;
    ASSERT_NO_THROW(LoadBlob(e, ads));
    const rapidjson::Value& raw = cat.RawEntry(e.key);

    // Components: exactly the blob's adsorbing (non-inert) isotherm components.
    for (const auto& c : e.components) {
      auto& iso = ads.adsorption_source_.GetIsotherm(c);
      ASSERT_TRUE(iso);
      EXPECT_NE(iso->GetIsothermType(), IsothermType::INERT) << c;
    }

    // Scalar properties equal the blob.
    const auto& props = raw["properties"];
    auto expectProp = [&](const char* name, double blob_value) {
      ASSERT_TRUE(props.HasMember(name)) << name;
      ASSERT_TRUE(props[name]["value"].IsNumber()) << name;
      EXPECT_TRUE(RelEq(props[name]["value"].GetDouble(), blob_value)) << name << ": " << props[name]["value"].GetDouble() << " vs blob " << blob_value;
    };
    expectProp("particle_density", ads.solid_density_.Get());
    expectProp("heat_capacity", ads.heat_capacity_.Get());
    expectProp("particle_voidage", ads.ep_.Get());
    expectProp("thermal_conductivity", ads.thermal_conductivity_.Get());

    // Heat of adsorption: null in the catalog <=> unset (or 0) in the blob.
    const auto& hoa = raw["heat_of_adsorption"];
    const auto& map = ads.adsorption_source_.heat_of_adsorption_;
    for (const auto& c : e.components) {
      ASSERT_TRUE(hoa.HasMember(c.c_str())) << c;
      const auto& v = hoa[c.c_str()]["value"];
      auto it = map.find(c);
      if (v.IsNull()) {
        EXPECT_TRUE(it == map.end() || it->second.Get() == 0.0) << c;
      } else {
        ASSERT_TRUE(it != map.end()) << c;
        EXPECT_TRUE(RelEq(v.GetDouble(), it->second.Get())) << c;
      }
    }
  }
}

TEST(AdsorbentLibraryTest, IsothermsEvaluateFinitelyOverStatedRange) {
  const auto& cat = Cat();
  for (const auto& e : cat.lib.entries) {
    SCOPED_TRACE(e.key);
    Adsorbent ads;
    LoadBlob(e, ads);
    const auto& isos = cat.RawEntry(e.key)["isotherms"];
    for (const auto& im : isos.GetArray()) {
      const std::string c = im["component"].GetString();
      // Stated fit range when given, else the generic envelope.
      double T0 = 273.15, T1 = 373.15, P0 = 1.0, P1 = 1e6;
      if (im["T_range_K"].IsArray()) { T0 = im["T_range_K"][0].GetDouble(); T1 = im["T_range_K"][1].GetDouble(); }
      if (im["P_range_Pa"].IsArray()) { P0 = std::max(1.0, im["P_range_Pa"][0].GetDouble()); P1 = im["P_range_Pa"][1].GetDouble(); }
      for (const auto& range : {std::array<double, 4>{T0, T1, P0, P1}, std::array<double, 4>{273.15, 373.15, 1.0, 1e6}}) {
        for (int i = 0; i <= 4; ++i) {
          const double T = range[0] + (range[1] - range[0]) * i / 4.0;
          for (int j = 0; j <= 12; ++j) {
            const double P = range[2] * std::pow(range[3] / range[2], j / 12.0);
            const double q = Loading(ads, c, T, P);
            ASSERT_TRUE(std::isfinite(q)) << c << " T=" << T << " P=" << P;
            ASSERT_GE(q, 0.0) << c << " T=" << T << " P=" << P;
          }
        }
      }
      // Adsorbs something at the top of the stated range.
      EXPECT_GT(Loading(ads, c, T0, P1), 0.0) << c;
    }
  }
}

TEST(AdsorbentLibraryTest, DeckDerivedEntriesMatchTheirDecks) {
  const auto& cat = Cat();
  int checked = 0;
  for (const auto& e : cat.lib.entries) {
    const auto& raw = cat.RawEntry(e.key);
    if (!raw.HasMember("deck_checks")) continue;
    for (const auto& dc : raw["deck_checks"].GetArray()) {
      SCOPED_TRACE(e.key + " vs " + dc["deck"].GetString());
      UserData deck;
      ASSERT_NO_THROW(LoadFromJson(deck, (kExamples / dc["deck"].GetString()).string()));
      Adsorbent* deck_ads = nullptr;
      for (auto& [id, a] : deck.adsorbent_directory.GetAll())
        if (a.name_ == dc["adsorbent"].GetString()) deck_ads = &deck.adsorbent_directory.Get(id);
      ASSERT_NE(deck_ads, nullptr) << "adsorbent " << dc["adsorbent"].GetString() << " not in deck";
      Adsorbent lib_ads;
      LoadBlob(e, lib_ads);
      for (const auto& c : e.components) {
        if (std::find(deck.fluid.components_.begin(), deck.fluid.components_.end(), c) == deck.fluid.components_.end())
          continue;  // e.g. H2 added to the Streb 13X entry from the paper's table
        for (double T : {273.15, 298.15, 323.15, 373.15})
          for (double P : {1e3, 1e4, 1e5, 1e6}) {
            const double ql = Loading(lib_ads, c, T, P), qd = Loading(*deck_ads, c, T, P);
            EXPECT_NEAR(ql, qd, 1e-9 * std::max(1.0, std::abs(qd))) << c << " T=" << T << " P=" << P;
          }
        ++checked;
      }
    }
  }
  EXPECT_GE(checked, 10);
}

namespace {
  // "Not stated in source (checked: ...)" / "Not stated in source as a single
  // value (checked: ...)" - same rule as scripts/docs/gen_adsorbent_library.py.
  bool NotStated(const rapidjson::Value& v) {
    static const std::regex re(R"(^Not stated in source(?: as a [^(]+)? \(checked: [^)]+\))");
    return v.IsString() && std::regex_search(v.GetString(), re);
  }
  bool Listed(const rapidjson::Value& e, const char* field) {
    if (!e.HasMember("unverified_fields")) return false;
    for (const auto& f : e["unverified_fields"].GetArray())
      if (std::string(f.GetString()) == field) return true;
    return false;
  }
  bool Sourced(const rapidjson::Value& rec) {
    const std::string prov = rec["provenance"].GetString();
    return (prov == "source" || prov == "derived") && !rec["value"].IsNull();
  }
}  // namespace

// Roadmap item 10, box 2: for every sourced entry each required field (fit
// range, heat of adsorption, particle density, heat capacity, particle size,
// example decks; the DOI is checked above) is sourced or explicitly marked
// "Not stated in source (checked: ...)"; a field that could not be checked is
// listed in unverified_fields. Nothing may be silently missing.
TEST(AdsorbentLibraryTest, RequiredFieldsAreSourcedOrExplicitlyNotStated) {
  const auto& cat = Cat();
  int incomplete = 0;
  for (const auto& e : cat.lib.entries) {
    SCOPED_TRACE(e.key);
    const auto& raw = cat.RawEntry(e.key);
    if (e.IsLegacy()) { ++incomplete; continue; }
    bool complete = true;
    auto require = [&](bool ok, const char* field, const std::string& what) {
      if (ok) return;
      complete = false;
      EXPECT_TRUE(Listed(raw, field)) << field << " " << what << " is neither sourced nor 'Not stated in source (checked: ...)'";
    };
    for (const auto& im : raw["isotherms"].GetArray()) {
      const std::string c = im["component"].GetString();
      require(im["T_range_K"].IsArray() || (im.HasMember("T_range_not_stated") && NotStated(im["T_range_not_stated"])), "T_range_K", c);
      require(im["P_range_Pa"].IsArray() || (im.HasMember("P_range_not_stated") && NotStated(im["P_range_not_stated"])), "P_range_Pa", c);
    }
    for (const auto& c : e.components) {
      const auto& h = raw["heat_of_adsorption"][c.c_str()];
      require(Sourced(h) || NotStated(h["ref"]), "heat_of_adsorption", c);
    }
    for (const char* p : {"particle_density", "heat_capacity", "particle_diameter"}) {
      const auto& rec = raw["properties"][p];
      require(Sourced(rec) || NotStated(rec["ref"]), p, "");
    }
    ASSERT_TRUE(raw.HasMember("example_decks") && raw["example_decks"].IsArray());
    EXPECT_TRUE(raw["example_decks"].Size() > 0 || (raw.HasMember("example_deck_note") && raw["example_deck_note"].IsString()))
        << "no example deck and no example_deck_note";
    if (raw.HasMember("unverified_fields")) EXPECT_TRUE(raw.HasMember("verification_note"));
    if (!complete) ++incomplete;
  }
  // Keep this honest: the count only goes down as sources are verified.
  EXPECT_LE(incomplete, 4);
}

// Roadmap item 10, box 3: tabulated source points load, convert to SI and the
// engine isotherm evaluates to a finite goodness of fit at them.
TEST(AdsorbentLibraryTest, SourcePointsLoadConvertAndGiveFiniteFitStatistics) {
  const auto& cat = Cat();
  int blocks = 0, stats = 0;
  for (const auto& e : cat.lib.entries) {
    const auto& raw = cat.RawEntry(e.key);
    Adsorbent ads;
    LoadBlob(e, ads);
    for (const auto& im : raw["isotherms"].GetArray()) {
      if (!im.HasMember("source_points")) continue;
      const std::string c = im["component"].GetString();
      SCOPED_TRACE(e.key + " / " + c);
      const auto& sp = im["source_points"];
      const std::string pu = sp["units"]["p"].GetString(), qu = sp["units"]["q"].GetString(), tu = sp["units"]["T"].GetString();
      const double pf = pu == "Pa" ? 1.0 : pu == "kPa" ? 1e3 : pu == "MPa" ? 1e6 : pu == "bar" ? 1e5 : NAN;
      const double qf = (qu == "mol/kg" || qu == "mmol/g") ? 1.0 : NAN;
      const double toff = tu == "K" ? 0.0 : tu == "C" ? 273.15 : NAN;
      ASSERT_TRUE(std::isfinite(pf) && std::isfinite(qf) && std::isfinite(toff)) << pu << " " << qu << " " << tu;
      const auto& T = sp["T"]; const auto& p = sp["p"]; const auto& q = sp["q"];
      ASSERT_TRUE(T.IsArray() && p.IsArray() && q.IsArray());
      ASSERT_GT(p.Size(), 0u);
      ASSERT_EQ(T.Size(), p.Size());
      ASSERT_EQ(q.Size(), p.Size());
      const std::string quantity = sp["quantity"].GetString();
      double mean = 0.0;
      for (rapidjson::SizeType i = 0; i < q.Size(); ++i) mean += q[i].GetDouble() * qf;
      mean /= q.Size();
      double sse = 0.0, sst = 0.0;
      for (rapidjson::SizeType i = 0; i < p.Size(); ++i) {
        const double Tk = T[i].GetDouble() + toff, pPa = p[i].GetDouble() * pf, qd = q[i].GetDouble() * qf;
        ASSERT_GT(Tk, 0.0);
        ASSERT_GE(pPa, 0.0);
        ASSERT_TRUE(std::isfinite(qd));
        const double qm = Loading(ads, c, Tk, pPa);
        ASSERT_TRUE(std::isfinite(qm)) << "T=" << Tk << " p=" << pPa;
        ASSERT_GE(qm, 0.0);
        sse += (qm - qd) * (qm - qd);
        sst += (qd - mean) * (qd - mean);
      }
      ++blocks;
      if (quantity == "excess") continue;  // model is absolute: no statistic
      const double rmse = std::sqrt(sse / p.Size()), r2 = 1.0 - sse / sst;
      EXPECT_TRUE(std::isfinite(rmse) && std::isfinite(r2)) << "rmse " << rmse << " r2 " << r2;
      std::printf("[source points] %-30s %-14s n=%3u  R2=%7.4f  RMSE=%.4f mol/kg\n", e.key.c_str(), c.c_str(), p.Size(), r2, rmse);
      ++stats;
    }
  }
  EXPECT_GE(blocks, 15);
  EXPECT_GE(stats, 13);
  // Unit conversion spot checks against the source files.
  const auto& lew = cat.RawEntry("lewatit_vp_oc_1065_low2023")["isotherms"][0]["source_points"];
  EXPECT_STREQ(lew["units"]["p"].GetString(), "kPa");  // AIF _units_pressure kPa
  EXPECT_DOUBLE_EQ(lew["p"][0].GetDouble() * 1e3, 9.670951621180774);  // Lewatit_CO2_288K.aif, first point
  const auto& tif = cat.RawEntry("tifsix_3_ni_low2024")["isotherms"][0]["source_points"];
  EXPECT_STREQ(tif["units"]["T"].GetString(), "C");  // AIF _units_temperature degC
  EXPECT_DOUBLE_EQ(tif["T"][0].GetDouble() + 273.15, 288.15);
}

}  // namespace adsim
