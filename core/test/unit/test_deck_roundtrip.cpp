// test_deck_roundtrip.cpp
//
// Decks survive save -> reload exactly, and a run does not depend on the
// order in which a deck lists its map entries.
//
//   * Every number in every bundled example deck is parsed exactly: the
//     reader cereal uses (rapidjson with kParseFullPrecisionFlag) gives the
//     same double as a correctly rounded strtod of the same text.
//   * A plain-number ExpressionParam's text parses (strtod, as Evaluate does)
//     to exactly its stored Value.
//   * What a loaded deck writes reloads to the identical deck: every double
//     bit for bit (rapidjson prints shortest round-trip digits) and every
//     string, including expression texts.
//   * Evaluated, a deck file and its re-saved copy are the identical deck.
//   * Reversing the order of every {key, value} map in a deck leaves a run
//     bit-identical. unordered_map iteration order differs between standard
//     libraries and between a deck and its re-saved copy. Competitive
//     isotherms used to sum over components in that order; with FMA
//     contraction (Apple clang on arm64) even a two-component sum depends on
//     it, which made ObjectIds.OldDeckLoadsDeterministicallyAndRunsIdentically
//     fail on macOS only. This test exposes such a dependence on any platform
//     (with three components it shows without FMA too).

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <cereal/archives/json.hpp>
#include <cereal/external/rapidjson/document.h>
#include <cereal/external/rapidjson/stringbuffer.h>
#include <cereal/external/rapidjson/writer.h>
#include <gtest/gtest.h>

#include "adsim/objects/user_data.h"
#include "adsim/simulation.h"

namespace adsim {
namespace {

const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";

std::string ReadFile(const std::filesystem::path& p) {
  std::ifstream ifs(p);
  std::stringstream ss;
  ss << ifs.rdbuf();
  return ss.str();
}

// cereal's default ("value0") root, without SaveToJson's sim_version bump.
std::string Snapshot(UserData& data) {
  std::ostringstream oss;
  {
    cereal::JSONOutputArchive ar(oss);
    ar(data);
  }
  return oss.str();
}

std::uint64_t Bits(double d) {
  std::uint64_t u;
  std::memcpy(&u, &d, sizeof u);
  return u;
}

bool SameDouble(double a, double b) { return Bits(a) == Bits(b) || (std::isnan(a) && std::isnan(b)); }

std::string Repr(double d) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.17g (0x%016llx)", d, static_cast<unsigned long long>(Bits(d)));
  return buf;
}

// Parses with the same flags cereal's JSONInputArchive uses.
void Parse(rapidjson::Document& doc, const std::string& text) {
  doc.Parse<rapidjson::kParseFullPrecisionFlag | rapidjson::kParseNanAndInfFlag>(text.c_str());
  ASSERT_FALSE(doc.HasParseError()) << "offset " << doc.GetErrorOffset();
}

std::string Join(const std::vector<std::string>& v, size_t max = 25) {
  std::string out;
  for (size_t i = 0; i < v.size() && i < max; ++i) out += "  " + v[i] + "\n";
  if (v.size() > max) out += "  ... and " + std::to_string(v.size() - max) + " more\n";
  return out;
}

std::vector<std::filesystem::path> ExampleDecks() {
  std::vector<std::filesystem::path> decks;
  for (const auto& e : std::filesystem::directory_iterator(kExamples))
    if (e.path().extension() == ".json") decks.push_back(e.path());
  std::sort(decks.begin(), decks.end());
  return decks;
}

bool IsKeyedEntry(const rapidjson::Value& v) { return v.IsObject() && v.HasMember("key"); }

bool SameKey(const rapidjson::Value& a, const rapidjson::Value& b) {
  const auto& ka = a["key"];
  const auto& kb = b["key"];
  if (ka.IsString() && kb.IsString()) return std::string(ka.GetString()) == kb.GetString();
  if (ka.IsNumber() && kb.IsNumber()) return SameDouble(ka.GetDouble(), kb.GetDouble());
  return false;
}

// A {Value, Expression} ExpressionParam whose expression is a plain number must
// parse (with strtod, exactly as ExpressionParam::Evaluate does) to Value.
void CheckExpression(const rapidjson::Value& obj, const std::string& path, std::vector<std::string>& diffs) {
  auto v = obj.FindMember("Value");
  auto e = obj.FindMember("Expression");
  if (v == obj.MemberEnd() || e == obj.MemberEnd() || !v->value.IsNumber() || !e->value.IsString()) return;
  ExpressionParam p;
  p.SetExpr(e->value.GetString());
  if (p.HasExpression()) return;  // a formula: its value is evaluated, not parsed
  const double parsed = std::strtod(e->value.GetString(), nullptr);
  if (!SameDouble(parsed, v->value.GetDouble()))
    diffs.push_back(path + ": Expression \"" + e->value.GetString() + "\" parses to " + Repr(parsed) +
                    " but Value is " + Repr(v->value.GetDouble()));
}

// Walks `from` against `back` and records every number or string that differs
// (and every plain-number expression that disagrees with its Value). Members
// only in `back` are ignored; members only in `from` are reported. Maps
// ({key, value} arrays) are matched by key, not position: an unordered_map's
// iteration order is a standard-library detail.
struct Walker {
  std::vector<std::string> diffs;

  void Walk(const rapidjson::Value& from, const rapidjson::Value& back, const std::string& path) {
    if (from.IsNumber() && back.IsNumber()) {
      if (!SameDouble(from.GetDouble(), back.GetDouble()))
        diffs.push_back(path + ": " + Repr(from.GetDouble()) + " -> " + Repr(back.GetDouble()));
      return;
    }
    if (from.GetType() != back.GetType()) {
      diffs.push_back(path + ": JSON type changed");
      return;
    }
    if (from.IsString()) {
      if (std::string(from.GetString()) != back.GetString())
        diffs.push_back(path + ": \"" + from.GetString() + "\" -> \"" + back.GetString() + "\"");
      return;
    }
    if (from.IsObject()) {
      CheckExpression(from, path, diffs);
      for (auto m = from.MemberBegin(); m != from.MemberEnd(); ++m) {
        const std::string name = m->name.GetString();
        // Archive bookkeeping, not deck data: cereal writes a class version, and
        // a polymorphic type's name with a flagged id, only at the first
        // occurrence, which depends on map iteration order.
        if (name == "cereal_class_version" || name == "polymorphic_id" || name == "polymorphic_name") continue;
        // Connection pair IDs are reassigned on every load from a process-wide
        // counter; they only pair the two ends of a connection.
        if (name == "Pair_ID") continue;
        auto it = back.FindMember(m->name);
        if (it == back.MemberEnd()) {
          diffs.push_back(path + "." + name + ": missing");
          continue;
        }
        Walk(m->value, it->value, path + "." + name);
      }
      return;
    }
    if (from.IsArray()) {
      if (from.Size() != back.Size()) {
        diffs.push_back(path + ": array size " + std::to_string(from.Size()) + " -> " +
                        std::to_string(back.Size()));
        return;
      }
      for (rapidjson::SizeType i = 0; i < from.Size(); ++i) {
        const rapidjson::Value* match = &back[i];
        if (IsKeyedEntry(from[i])) {
          match = nullptr;
          for (rapidjson::SizeType j = 0; j < back.Size(); ++j)
            if (IsKeyedEntry(back[j]) && SameKey(from[i], back[j])) match = &back[j];
          if (!match) {
            diffs.push_back(path + "[" + std::to_string(i) + "]: key missing");
            continue;
          }
        }
        Walk(from[i], *match, path + "[" + std::to_string(i) + "]");
      }
    }
  }
};

// Compares every number in `parsed` with a correctly rounded strtod of its
// text in `as_text` (the same document parsed with kParseNumbersAsStringsFlag),
// and checks every plain-number expression against its Value.
void CheckFileNumbers(const rapidjson::Value& parsed, const rapidjson::Value& as_text, const std::string& path,
                      std::vector<std::string>& diffs) {
  if (parsed.IsNumber()) {
    const char* text = as_text.GetString();
    const double reference = std::strtod(text, nullptr);
    if (!SameDouble(parsed.GetDouble(), reference))
      diffs.push_back(path + ": \"" + text + "\" parsed as " + Repr(parsed.GetDouble()) + ", strtod gives " +
                      Repr(reference));
  } else if (parsed.IsObject()) {
    CheckExpression(parsed, path, diffs);
    for (auto m = parsed.MemberBegin(); m != parsed.MemberEnd(); ++m)
      CheckFileNumbers(m->value, as_text[m->name.GetString()], path + "." + m->name.GetString(), diffs);
  } else if (parsed.IsArray()) {
    for (rapidjson::SizeType i = 0; i < parsed.Size(); ++i)
      CheckFileNumbers(parsed[i], as_text[i], path + "[" + std::to_string(i) + "]", diffs);
  }
}

TEST(DeckRoundTrip, ExampleDecksSurviveSaveAndReloadBitForBit) {
  const auto decks = ExampleDecks();
  ASSERT_GE(decks.size(), 30u);
  for (const auto& deck : decks) {
    SCOPED_TRACE(deck.filename().string());
    const std::string file_text = ReadFile(deck);

    // 1. Every number in the file is parsed exactly; plain-number expression
    //    texts match their values.
    rapidjson::Document file_doc, file_text_doc;
    Parse(file_doc, file_text);
    file_text_doc.Parse<rapidjson::kParseNumbersAsStringsFlag | rapidjson::kParseNanAndInfFlag>(file_text.c_str());
    ASSERT_FALSE(file_text_doc.HasParseError());
    std::vector<std::string> file_diffs;
    CheckFileNumbers(file_doc, file_text_doc, "", file_diffs);
    EXPECT_TRUE(file_diffs.empty()) << "deck file:\n" << Join(file_diffs);

    // 2. What the loaded deck writes reloads to the identical deck.
    UserData a;
    ASSERT_NO_THROW(LoadFromJson(a, deck.string()));
    const std::string saved = Snapshot(a);
    UserData b;
    ASSERT_NO_THROW(LoadFromJsonText(b, saved));
    {
      rapidjson::Document saved_doc, resaved_doc;
      Parse(saved_doc, saved);
      Parse(resaved_doc, Snapshot(b));
      Walker forward, backward;
      forward.Walk(saved_doc["value0"], resaved_doc["value0"], "value0");
      backward.Walk(resaved_doc["value0"], saved_doc["value0"], "value0");
      EXPECT_TRUE(forward.diffs.empty()) << "saved deck -> reloaded deck:\n" << Join(forward.diffs);
      EXPECT_TRUE(backward.diffs.empty()) << "reloaded deck -> saved deck:\n" << Join(backward.diffs);
    }

    // 3. Evaluated, the deck file and its re-saved copy are the same deck
    //    (the pair ObjectIds.OldDeckLoadsDeterministicallyAndRunsIdentically runs).
    a.EvaluateAllExpressions();
    b.EvaluateAllExpressions();
    {
      rapidjson::Document da, db;
      Parse(da, Snapshot(a));
      Parse(db, Snapshot(b));
      Walker evaluated;
      evaluated.Walk(da["value0"], db["value0"], "value0");
      EXPECT_TRUE(evaluated.diffs.empty()) << "file vs re-saved copy after EvaluateAllExpressions:\n"
                                           << Join(evaluated.diffs);
    }
  }
}

// --- Reordering a deck's maps --------------------------------------------

// cereal writes a class's cereal_class_version only at its first occurrence.
// Moves those markers from `from` (the old first map entry) to the same place
// in `to` (the new first entry).
void MoveClassVersions(rapidjson::Value& from, rapidjson::Value& to, rapidjson::Document::AllocatorType& alloc) {
  if (from.IsObject() && to.IsObject()) {
    auto ver = from.FindMember("cereal_class_version");
    if (ver != from.MemberEnd() && !to.HasMember("cereal_class_version")) {
      // cereal reads it before the class's members: put it first.
      rapidjson::Value rebuilt(rapidjson::kObjectType);
      rebuilt.AddMember("cereal_class_version", rapidjson::Value(ver->value, alloc), alloc);
      for (auto m = to.MemberBegin(); m != to.MemberEnd(); ++m)
        rebuilt.AddMember(rapidjson::Value(m->name, alloc), rapidjson::Value(m->value, alloc), alloc);
      to.Swap(rebuilt);
      from.RemoveMember("cereal_class_version");
    }
    for (auto m = from.MemberBegin(); m != from.MemberEnd(); ++m) {
      auto t = to.FindMember(m->name);
      if (t != to.MemberEnd()) MoveClassVersions(m->value, t->value, alloc);
    }
  } else if (from.IsArray() && to.IsArray()) {
    if (from.Size() > 0 && to.Size() > 0) MoveClassVersions(from[0], to[0], alloc);
  }
}

// Reverses every serialised map ({key, value} array) in the deck JSON.
void ReverseMaps(rapidjson::Value& v, rapidjson::Document::AllocatorType& alloc) {
  if (v.IsObject()) {
    for (auto m = v.MemberBegin(); m != v.MemberEnd(); ++m) ReverseMaps(m->value, alloc);
  } else if (v.IsArray()) {
    for (auto& e : v.GetArray()) ReverseMaps(e, alloc);
    if (v.Size() > 1 && IsKeyedEntry(v[0])) {
      for (rapidjson::SizeType i = 0, j = v.Size() - 1; i < j; ++i, --j) v[i].Swap(v[j]);
      MoveClassVersions(v[v.Size() - 1], v[0], alloc);
    }
  }
}

// cereal names a polymorphic type (polymorphic_id with the msb set, plus
// polymorphic_name) at its first occurrence in document order and refers to it
// by id afterwards. Reordering entries breaks that, so renumber in the new order.
constexpr std::uint32_t kPolyMsb = 0x80000000u;

void CollectPolyNames(const rapidjson::Value& v, std::map<std::uint32_t, std::string>& names) {
  if (v.IsObject()) {
    auto id = v.FindMember("polymorphic_id");
    auto nm = v.FindMember("polymorphic_name");
    if (id != v.MemberEnd() && nm != v.MemberEnd()) names[id->value.GetUint() & ~kPolyMsb] = nm->value.GetString();
    for (auto m = v.MemberBegin(); m != v.MemberEnd(); ++m) CollectPolyNames(m->value, names);
  } else if (v.IsArray()) {
    for (const auto& e : v.GetArray()) CollectPolyNames(e, names);
  }
}

void RenumberPoly(rapidjson::Value& v, const std::map<std::uint32_t, std::string>& old_names,
                  std::map<std::string, std::uint32_t>& new_ids, rapidjson::Document::AllocatorType& alloc) {
  if (v.IsObject()) {
    auto id = v.FindMember("polymorphic_id");
    if (id != v.MemberEnd() && id->value.GetUint() != 0) {
      const std::string name = old_names.at(id->value.GetUint() & ~kPolyMsb);
      v.RemoveMember("polymorphic_name");
      id = v.FindMember("polymorphic_id");
      auto it = new_ids.find(name);
      if (it == new_ids.end()) {
        const std::uint32_t nid = static_cast<std::uint32_t>(new_ids.size()) + 1;
        new_ids[name] = nid;
        id->value.SetUint(nid | kPolyMsb);
        // polymorphic_name directly follows polymorphic_id: rebuild the members.
        rapidjson::Value rebuilt(rapidjson::kObjectType);
        for (auto m = v.MemberBegin(); m != v.MemberEnd(); ++m) {
          const bool is_id = std::string(m->name.GetString()) == "polymorphic_id";
          rebuilt.AddMember(rapidjson::Value(m->name, alloc), rapidjson::Value(m->value, alloc), alloc);
          if (is_id) rebuilt.AddMember("polymorphic_name", rapidjson::Value(name.c_str(), alloc), alloc);
        }
        v.Swap(rebuilt);
      } else {
        id->value.SetUint(it->second);
      }
    }
    for (auto m = v.MemberBegin(); m != v.MemberEnd(); ++m) RenumberPoly(m->value, old_names, new_ids, alloc);
  } else if (v.IsArray()) {
    for (auto& e : v.GetArray()) RenumberPoly(e, old_names, new_ids, alloc);
  }
}

// The deck JSON with every serialised map's entries in reverse order.
std::string WithMapsReversed(const std::string& json) {
  rapidjson::Document doc;
  Parse(doc, json);
  std::map<std::uint32_t, std::string> poly_names;
  CollectPolyNames(doc, poly_names);
  ReverseMaps(doc, doc.GetAllocator());
  std::map<std::string, std::uint32_t> new_ids;
  RenumberPoly(doc, poly_names, new_ids, doc.GetAllocator());
  rapidjson::StringBuffer buf;
  rapidjson::Writer<rapidjson::StringBuffer, rapidjson::UTF8<>, rapidjson::UTF8<>, rapidjson::CrtAllocator,
                    rapidjson::kWriteNanAndInfFlag>
      w(buf);
  doc.Accept(w);
  return buf.GetString();
}

TEST(DeckRoundTrip, RunDoesNotDependOnMapOrderInTheDeck) {
  // Two- and three-component decks with competitive isotherms. One cycle is
  // enough: before the fix Perez2019_A differed at ~2e-6 after one cycle in a
  // clang/libc++ build with FMA contraction.
  for (const char* name : {"Perez2019_A.json", "Casas2012_Fig4_50_50_15bar_25C.json", "Streb2021_A1.json"}) {
    SCOPED_TRACE(name);
    const std::string deck = (kExamples / name).string();
    const std::string reversed_json = WithMapsReversed(ReadFile(deck));

    // The reversal must reach the loaded deck (it changes what is saved back).
    UserData a, b;
    LoadFromJson(a, deck);
    ASSERT_NO_THROW(LoadFromJsonText(b, reversed_json));
    ASSERT_NE(Snapshot(a), Snapshot(b)) << "reversing the deck's maps changed nothing";

    Simulation file = Simulation::FromFile(deck);
    Simulation reversed = Simulation::FromJson(reversed_json);
    for (Simulation* s : {&file, &reversed}) {
      s->SetParameter("cycle.number_of_cycles", 1);
      ASSERT_TRUE(s->Run().Ok());
    }
    const auto k1 = file.KPIs();
    const auto k2 = reversed.KPIs();
    ASSERT_EQ(k1.size(), k2.size());
    for (size_t i = 0; i < k1.size(); ++i) {
      EXPECT_TRUE(SameDouble(k1[i].purity, k2[i].purity)) << Repr(k1[i].purity) << " vs " << Repr(k2[i].purity);
      EXPECT_TRUE(SameDouble(k1[i].recovery, k2[i].recovery))
          << Repr(k1[i].recovery) << " vs " << Repr(k2[i].recovery);
    }
    EXPECT_EQ(file.Series("pressure").values, reversed.Series("pressure").values);
    EXPECT_EQ(file.Series("temperature").values, reversed.Series("temperature").values);
    EXPECT_EQ(file.Series("mole_fraction", "CarbonDioxide").values,
              reversed.Series("mole_fraction", "CarbonDioxide").values);
  }
}

}  // namespace
}  // namespace adsim
