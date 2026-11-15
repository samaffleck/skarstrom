/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "adsim/objects/component_migration.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string_view>

#include <cereal/external/rapidjson/document.h>
#include <cereal/external/rapidjson/stringbuffer.h>
#include <cereal/external/rapidjson/writer.h>

#include "adsim/objects/fluid.h"

namespace adsim {

namespace {

namespace rj = CEREAL_RAPIDJSON_NAMESPACE;

// Fields that hold one component name.
constexpr const char* kNameFields[] = {"Component", "Key_Component", "Product_Component", "Key_Reactant"};

bool IsNameField(std::string_view name) {
  for (const char* f : kNameFields)
    if (name == f) return true;
  return false;
}

bool IsMoleFractionMap(std::string_view name) {
  return name.find("Mole_Fraction") != std::string_view::npos || name.find("Composition") != std::string_view::npos;
}

bool IsString(const rj::Value& v, std::string_view s) {
  return v.IsString() && std::string_view(v.GetString(), v.GetStringLength()) == s;
}

// A cereal map: a non-empty array of {"key": string, "value": ...} objects.
bool IsComponentMap(const rj::Value& a) {
  if (!a.IsArray() || a.Empty()) return false;
  for (const auto& e : a.GetArray()) {
    if (!e.IsObject()) return false;
    auto k = e.FindMember("key");
    if (k == e.MemberEnd() || !k->value.IsString() || !e.HasMember("value")) return false;
  }
  return true;
}

bool IsStringArray(const rj::Value& a) {
  if (!a.IsArray() || a.Empty()) return false;
  for (const auto& e : a.GetArray())
    if (!e.IsString()) return false;
  return true;
}

// True (and the number in `out`) if `s` is a plain number.
bool ParseNumber(const char* s, double& out) {
  char* end = nullptr;
  out = std::strtod(s, &end);
  if (end == s) return false;
  while (*end == ' ') ++end;
  return *end == '\0';
}

std::string FormatNumber(double x) {
  char buf[32];
  std::snprintf(buf, sizeof(buf), "%.17g", x);
  return buf;
}

// Add the mole fraction `from` to `into`: a number, or an expression
// parameter {"Value", "Expression", ...}.
void AddMoleFraction(rj::Value& into, const rj::Value& from, rj::Document::AllocatorType& alloc) {
  if (into.IsNumber() && from.IsNumber()) {
    into.SetDouble(into.GetDouble() + from.GetDouble());
    return;
  }
  if (!into.IsObject() || !from.IsObject()) return;
  auto iv = into.FindMember("Value");
  auto fv = from.FindMember("Value");
  if (iv != into.MemberEnd() && fv != from.MemberEnd() && iv->value.IsNumber() && fv->value.IsNumber())
    iv->value.SetDouble(iv->value.GetDouble() + fv->value.GetDouble());
  auto ie = into.FindMember("Expression");
  auto fe = from.FindMember("Expression");
  if (ie != into.MemberEnd() && fe != from.MemberEnd() && ie->value.IsString() && fe->value.IsString()) {
    double a = 0.0, b = 0.0;
    std::string sum;
    if (ParseNumber(ie->value.GetString(), a) && ParseNumber(fe->value.GetString(), b))
      sum = FormatNumber(a + b);
    else
      sum = "(" + std::string(ie->value.GetString()) + ") + (" + fe->value.GetString() + ")";
    ie->value.SetString(sum.c_str(), static_cast<rj::SizeType>(sum.size()), alloc);
  }
}

struct Renamer {
  std::string_view from;
  std::string_view to;
  bool merge;  // the deck already has `to`: drop `from` instead of renaming it
  rj::Document::AllocatorType& alloc;
  bool changed = false;

  void SetString(rj::Value& v, std::string_view s) {
    v.SetString(s.data(), static_cast<rj::SizeType>(s.size()), alloc);
    changed = true;
  }

  void Map(rj::Value& a, std::string_view name) {
    rj::SizeType old_i = a.Size(), new_i = a.Size();
    for (rj::SizeType i = 0; i < a.Size(); ++i) {
      const rj::Value& k = a[i]["key"];
      if (IsString(k, from)) old_i = i;
      else if (IsString(k, to)) new_i = i;
    }
    if (old_i == a.Size()) return;
    const bool mole_fraction = IsMoleFractionMap(name);
    if (new_i != a.Size()) {
      if (mole_fraction) AddMoleFraction(a[new_i]["value"], a[old_i]["value"], alloc);
      a.Erase(a.Begin() + old_i);
      changed = true;
    } else if (merge && !mole_fraction) {
      a.Erase(a.Begin() + old_i);
      changed = true;
    } else {
      SetString(a[old_i]["key"], to);
    }
  }

  void Names(rj::Value& a) {
    bool has_new = false;
    for (const auto& e : a.GetArray()) has_new = has_new || IsString(e, to);
    for (rj::SizeType i = 0; i < a.Size();) {
      if (!IsString(a[i], from)) {
        ++i;
      } else if (has_new) {
        a.Erase(a.Begin() + i);
        changed = true;
      } else {
        SetString(a[i], to);
        has_new = true;
        ++i;
      }
    }
  }

  void Path(rj::Value& v) {
    const std::string old_ref = "[" + std::string(from) + "]";
    const std::string new_ref = "[" + std::string(to) + "]";
    std::string s(v.GetString(), v.GetStringLength());
    bool hit = false;
    for (size_t pos = s.find(old_ref); pos != std::string::npos; pos = s.find(old_ref, pos + new_ref.size())) {
      s.replace(pos, old_ref.size(), new_ref);
      hit = true;
    }
    if (hit) SetString(v, s);
  }

  void Walk(rj::Value& v, std::string_view name) {
    if (v.IsObject()) {
      for (auto& m : v.GetObject()) {
        std::string_view member(m.name.GetString(), m.name.GetStringLength());
        if (m.value.IsString()) {
          if (IsNameField(member) && IsString(m.value, from)) SetString(m.value, to);
          else if (member == "Path") Path(m.value);
        } else {
          Walk(m.value, member);
        }
      }
    } else if (v.IsArray()) {
      if (IsComponentMap(v)) Map(v, name);
      else if (IsStringArray(v)) Names(v);
      for (auto& e : v.GetArray()) Walk(e, name);
    }
  }
};

// The deck's Fluid "Components" array, or nullptr (an adsorbent or a cycle
// file, or a deck too old to have one).
const rj::Value* FluidComponents(const rj::Document& doc) {
  if (!doc.IsObject()) return nullptr;
  for (const char* root : {"value0", "simulation"}) {
    auto r = doc.FindMember(root);
    if (r == doc.MemberEnd() || !r->value.IsObject()) continue;
    auto f = r->value.FindMember("Fluid");
    if (f == r->value.MemberEnd() || !f->value.IsObject()) continue;
    auto c = f->value.FindMember("Components");
    if (c != f->value.MemberEnd() && c->value.IsArray()) return &c->value;
  }
  return nullptr;
}

bool Contains(const rj::Value& names, std::string_view s) {
  for (const auto& e : names.GetArray())
    if (IsString(e, s)) return true;
  return false;
}

constexpr unsigned kParseFlags = rj::kParseFullPrecisionFlag | rj::kParseNanAndInfFlag;

}  // namespace

bool MigrateRenamedComponents(std::string& json) {
  bool any = false;
  for (const auto& r : Components::renamedComponents)
    any = any || json.find("\"" + r.first + "\"") != std::string::npos;
  if (!any) return false;

  rj::Document doc;
  doc.Parse<kParseFlags>(json.c_str(), json.size());
  if (doc.HasParseError()) return false;

  bool changed = false;
  for (const auto& [from, to] : Components::renamedComponents) {
    const rj::Value* components = FluidComponents(doc);
    if (components && !Contains(*components, from)) continue;
    const bool merge = components && Contains(*components, to);
    Renamer renamer{from, to, merge, doc.GetAllocator()};
    renamer.Walk(doc, "");
    changed = changed || renamer.changed;
  }
  if (!changed) return false;

  rj::StringBuffer out;
  rj::Writer<rj::StringBuffer, rj::UTF8<>, rj::UTF8<>, rj::CrtAllocator, rj::kWriteNanAndInfFlag> writer(out);
  doc.Accept(writer);
  json.assign(out.GetString(), out.GetSize());
  return true;
}

}  // namespace adsim
