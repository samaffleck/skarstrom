/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// Every editable value of every example deck has a parameter-registry path.
//
// Method: save the deck, give every registry parameter a value it cannot
// already have (a unique number, the other boolean, the next enum choice, a
// marked string), save again, and compare the two saves leaf by leaf. A leaf
// that did not change is not reachable through the registry. The leaves that
// are not editable values (IDs, versions, names and other structure, run
// results) or that are set through the structural edits of deck_edit.h are
// listed below with the reason; any other unchanged leaf fails the test.

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include <cereal/external/rapidjson/document.h>

#include "adsim/objects/deck_builder.h"
#include "adsim/objects/deck_io.h"
#include "adsim/objects/user_data.h"
#include "adsim/parameters/parameter_registry.h"

namespace adsim {
namespace {

namespace rj = CEREAL_RAPIDJSON_NAMESPACE;
namespace fs = std::filesystem;

const fs::path kExamples = fs::current_path() / ".." / "examples" / "decks";

std::string ReadText(const fs::path& p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

std::string Save(UserData& d) {
  d.info.sim_version = 0;
  std::string out;
  SaveToString(d, out);
  return out;
}

// A leaf's generic path: object keys as they are, array positions as [], and
// the "value" of a {"key", "value"} map entry as [*].
void CollectUnchanged(const rj::Value& a, const rj::Value& b, const std::string& path, std::set<std::string>& out) {
  if (a.IsObject() && b.IsObject()) {
    // An ExpressionParam is one leaf (its value and its text).
    if (a.HasMember("Value") && a.HasMember("Expression")) {
      if (a["Value"] == b["Value"] && a["Expression"] == b["Expression"]) out.insert(path);
      return;
    }
    for (auto it = a.MemberBegin(); it != a.MemberEnd(); ++it) {
      const auto bm = b.FindMember(it->name);
      if (bm == b.MemberEnd()) continue;
      CollectUnchanged(it->value, bm->value, path + "/" + it->name.GetString(), out);
    }
    return;
  }
  if (a.IsArray() && b.IsArray()) {
    const rj::SizeType n = std::min(a.Size(), b.Size());
    for (rj::SizeType i = 0; i < n; ++i) {
      const bool entry = a[i].IsObject() && a[i].HasMember("key") && a[i].HasMember("value");
      if (entry) {
        CollectUnchanged(a[i]["value"], b[i]["value"], path + "[*]", out);
      } else {
        CollectUnchanged(a[i], b[i], path + "[]", out);
      }
    }
    return;
  }
  if (a == b) out.insert(path);
}

// Give every parameter a value it does not have yet.
void MarkEverything(ParameterRegistry& reg) {
  double next = 7.123456789e5;
  int string_n = 0;
  for (auto& p : reg.All()) {
    const auto& info = p.Info();
    switch (info.kind) {
      case ParameterKind::Real:
      case ParameterKind::Expression:
        p.SetValue(next);
        next += 1.0;
        break;
      case ParameterKind::Boolean:
        p.SetValue(p.GetValue() != 0.0 ? 0.0 : 1.0);
        break;
      case ParameterKind::Integer: {
        const double v = p.GetValue();
        if (!info.choices.empty()) {
          if (info.choices.size() > 1)
            p.SetValue(static_cast<double>((static_cast<long long>(v) + 1) % static_cast<long long>(info.choices.size())));
        } else {
          p.SetValue(v + 1.0);
        }
        break;
      }
      case ParameterKind::String:
        p.SetText("marked " + std::to_string(++string_n));
        break;
    }
  }
}

// Leaves that are not values the registry should set, with the reason.
struct Exception {
  const char* pattern;  // a regex over the generic path
  const char* reason;
};

const std::vector<Exception>& Exceptions() {
  static const std::vector<Exception> e = {
      // Bookkeeping.
      {R"(.*/cereal_class_version$)", "archive version"},
      {R"(.*/polymorphic_(id|name)$)", "archive type marker"},
      {R"(.*/ptr_wrapper/(valid|id)$)", "archive pointer marker"},
      {R"(.*/DisplayUnitIndex$)", "display unit of a value (front-end preference)"},
      {R"(.*/(ID|Step_ID|Endpoint_ID|Pair_ID|Next_ID|Next_Layer_ID|Next_Reaction_ID|Layer_ID|Adsorbent_ID|Reaction_ID|Source_Step_ID|Source_Endpoint_ID)$)",
       "object ID"},
      {R"(^/simulation/System_Information/Version$)", "save counter"},
      // Structure, edited through deck_edit.h / deck_builder.h.
      {R"(.*/(Name|Stream_Name|Step_Name|Layer_Name)$)", "name (UserData::RenameStep / RenameStream, or the object's name)"},
      {R"(^/simulation/Fluid/Components\[\]$)", "component list (AddComponentToDeck / RemoveComponentFromDeck)"},
      {R"(.*/Boundary_Condition_(Location|Direction)$)", "which end and which way (SetBoundaryConditionKind, build::Set*Boundary)"},
      {R"(.*/Connection_Type$)", "connection state (ConnectEndpoints / DisconnectPair)"},
      {R"(.*/Use_Pump$)", "pump on or off (SetPumpKind)"},
      {R"(.*/Buffer_Tank_Source/(Active|Source_Step_Name|Source_Endpoint_Name)$)", "buffer tank source (build::SetBufferTank)"},
      {R"(.*/Isotherm_Map\[\*\]$)", "which isotherm a component uses (SetIsothermKind)"},
      {R"(^/simulation/Bed/Layers/Layers\[\]/Kind$)", "packed, monolith or laminate (build::SetLayerKind)"},
      {R"(^/simulation/Bed/Shape/index$)", "bed shape (build::SetBedShape)"},
      {R"(^/simulation/Bed/Wall/Material$)", "wall material, which sets its properties (build::SetWallMaterial)"},
      {R"(^/simulation/KPI_Config/Per_Component\[\]/(Product_Component|(Feed|Product|Recycle)_Endpoints\[\]/(Endpoint|Location|Step))$)",
       "KPI endpoints (build::SetKPI)"},
      {R"(.*/Reactions\[\]/Effectiveness/(Mode|Key_Component)$)", "effectiveness model and key reactant (build::SetEffectiveness)"},
      {R"(.*/Reactions\[\]/Stoichiometry\[\]/Component$)", "which components react (build::SetReactionStoichiometry)"},
      {R"(.*/Reactions\[\]/Stoichiometry\[\]/Phase$)", "phase of a stoichiometry entry: only the gas phase exists"},
      {R"(^/simulation/Experimental_Data/.*)", "measured data series (sim.experimental_data: add, load_csv)"},
      // Values of a variant the deck does not use: the registry lists a
      // field while it takes part in the model, and the structural edit that
      // activates it makes it listed.
      {R"(^/simulation/Bed/Layers/Layers\[\]/Structured_Geometry/.*)", "structured-layer geometry of a packed layer (listed for monoliths and laminates)"},
      {R"(^/simulation/Cycle/Step_List\[\]/Step_Fraction$)", "step fraction in absolute step-time mode (listed in fraction mode)"},
      {R"(^/simulation/Cycle/Total_Cycle_Time$)", "total cycle time in absolute step-time mode (listed in fraction mode)"},
      {R"(.*/Boundary_Pressure/ptr_wrapper/data/Start_Value$)", "start of a pressure ramp: the bed pressure when the step starts, set at run time"},
      {R"(.*/Boundary_Pressure/ptr_wrapper/data/(Ramp_Time|Step_Time)$)", "pressure ramp time while the ramp lasts the step (listed once <kind>_use_step_time is off)"},
      {R"(.*/Pump/ptr_wrapper/data/(Suction_Pressure|Discharge_Pressure|value0/.*)$)", "pump of an endpoint whose pump is off (listed once SetPumpKind turns it on)"},
      {R"(.*/Reactions\[\]/Effectiveness/(Effectiveness|Tortuosity)$)", "the effectiveness value of the mode not in use (eta for constant, tortuosity for Thiele)"},
      {R"(.*/Buffer_Tank_Source/Percentage$)", "share drawn by an inactive buffer tank (listed once the tank is set)"},
      // Legacy.
      {R"(.*/Boundary_Pressure/ptr_wrapper/data/Proportional_Gain$)", "unused field kept for the archive layout of the exponential profile"},
      {R"(^/simulation/Bed/Layers/Layers\[\]/Layer_Equations/Thermal_Equilibrium_Mode$)", "moved to the bed (bed.thermal_equilibrium_mode); discarded on load"},
      // Run results and caches.
      {R"(.*/Cycle_Endpoint_Summary/.*)", "result of the last run (buffer tank source data)"},
  };
  return e;
}

bool Excepted(const std::string& path) {
  for (const auto& e : Exceptions())
    if (std::regex_match(path, std::regex(e.pattern))) return true;
  return false;
}

// The generic paths of `d`'s leaves that no registry parameter reaches and
// no exception explains.
std::set<std::string> Uncovered(UserData& d) {
  const std::string before = Save(d);
  ParameterRegistry reg(d);
  MarkEverything(reg);
  const std::string after = Save(d);
  rj::Document a, b;
  a.Parse(before.c_str());
  b.Parse(after.c_str());
  std::set<std::string> unchanged, out;
  if (a.HasParseError() || b.HasParseError()) return {"(parse error)"};
  CollectUnchanged(a, b, "", unchanged);
  for (const auto& path : unchanged)
    if (!Excepted(path)) out.insert(path);
  return out;
}

}  // namespace

// The example decks have no reactions: a deck with one reaction of every
// rate kind (an LHHW one with an adsorption term, a reversible one, a
// rate expression), both effectiveness modes and both bases.
TEST(RegistryCoverage, EveryReactionValueHasAPath) {
  UserData d;
  for (const char* c : {"Methane", "Water", "CarbonMonoxide", "Hydrogen", "CarbonDioxide"}) build::AddComponent(d, c);
  const AdsorbentID cat = build::AddAdsorbent(d, "Ni catalyst");
  const std::map<std::string, double> smr{{"Methane", -1}, {"Water", -1}, {"CarbonMonoxide", 1}, {"Hydrogen", 3}};
  const std::map<std::string, double> wgs{{"CarbonMonoxide", -1}, {"Water", -1}, {"CarbonDioxide", 1}, {"Hydrogen", 1}};
  build::AddReaction(d, cat, "power", smr, "power_law");
  const ReactionID rev = build::AddReaction(d, cat, "reversible", wgs, "reversible_power_law");
  build::SetReactionBasis(d, cat, rev, "concentration");
  const ReactionID lhhw = build::AddReaction(d, cat, "lhhw", smr, "lhhw");
  build::AddAdsorptionTerm(d, cat, lhhw);
  build::SetEffectiveness(d, cat, lhhw, "thiele", "Methane");
  build::AddReaction(d, cat, "expression", wgs, "expression");
  const std::set<std::string> missing = Uncovered(d);
  std::string report;
  for (const auto& p : missing) report += "  " + p + "\n";
  EXPECT_TRUE(missing.empty()) << missing.size() << " reaction leaves have no registry path:\n" << report;
}

TEST(RegistryCoverage, EveryEditableValueOfEveryExampleDeckHasAPath) {
  std::map<std::string, std::set<std::string>> missing;  // generic path -> decks
  size_t decks = 0;
  for (const auto& entry : fs::directory_iterator(kExamples)) {
    if (entry.path().extension() != ".json") continue;
    ++decks;
    const std::string name = entry.path().stem().string();
    UserData d;
    ASSERT_NO_THROW(LoadDeckFlexible(d, ReadText(entry.path()))) << name;
    for (const auto& path : Uncovered(d)) missing[path].insert(name);
  }
  ASSERT_GE(decks, 30u);
  std::string report;
  for (const auto& [path, where] : missing) {
    report += "  " + path + "  (" + std::to_string(where.size()) + " decks, e.g. " + *where.begin() + ")\n";
  }
  EXPECT_TRUE(missing.empty()) << missing.size() << " deck leaves have no registry path:\n" << report;
}

}  // namespace adsim
