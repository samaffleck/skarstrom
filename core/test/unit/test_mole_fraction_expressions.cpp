// test_mole_fraction_expressions.cpp
//
// Inlet and initial mole fractions are ExpressionParams (an ExpressionMap per
// boundary, stream and the initial conditions), so a composition can be a
// formula of user variables like any other deck value. These tests pin:
//   * user-variable formulas drive the inlet and initial compositions, and the
//     sum-to-1 checks use the evaluated values;
//   * the expressions survive save and reload;
//   * the registry lists the mole fractions as Expression parameters;
//   * a deck written before the change (plain-number maps) still loads.

#include <filesystem>
#include <sstream>
#include <string>
#include <vector>

#include <cereal/archives/json.hpp>
#include <cereal/external/rapidjson/document.h>
#include <cereal/external/rapidjson/stringbuffer.h>
#include <cereal/external/rapidjson/writer.h>
#include <gtest/gtest.h>

#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/objects/deck_io.h"
#include "adsim/objects/user_data.h"
#include "adsim/parameters/parameter_registry.h"

namespace adsim {
namespace {

const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";

// Perez2019_A: CarbonDioxide + Nitrogen, pressure boundaries and a flow stream.
UserData LoadPerez() {
  UserData data;
  LoadFromJson(data, (kExamples / "Perez2019_A.json").string());
  return data;
}

std::string Snapshot(UserData& data) {
  std::ostringstream oss;
  {
    cereal::JSONOutputArchive ar(oss);
    ar(data);
  }
  return oss.str();
}

// Every mole-fraction map in the deck: the initial conditions, then each
// pressure boundary and stream.
std::vector<ExpressionMap*> AllMoleFractionMaps(UserData& data) {
  std::vector<ExpressionMap*> out{&data.ic.mole_fraction};
  for (auto& step : data.cycle.step_list_) {
    for (IBoundaryCondition* bc : {step.top_boundary_condition_.get(), step.bottom_boundary_condition_.get()}) {
      if (auto* p = dynamic_cast<PressureBoundaryCondition*>(bc)) out.push_back(&p->mole_fraction_);
      else if (auto* f = dynamic_cast<FlowBoundaryCondition*>(bc))
        for (auto& s : f->streams_)
          if (s) out.push_back(&s->mole_fraction);
    }
  }
  return out;
}

// Same maps, entry for entry (value and expression text). Whole-deck snapshots
// cannot be compared as text: unordered_map order differs after a reload.
void ExpectSameMoleFractions(UserData& a, UserData& b) {
  auto ma = AllMoleFractionMaps(a);
  auto mb = AllMoleFractionMaps(b);
  ASSERT_EQ(ma.size(), mb.size());
  for (size_t i = 0; i < ma.size(); ++i) {
    ASSERT_EQ(ma[i]->size(), mb[i]->size()) << "map " << i;
    for (const auto& [component, y] : *ma[i]) {
      const ExpressionParam& yb = mb[i]->at(component);
      EXPECT_EQ(yb.Get(), y.Get()) << "map " << i << " " << component;
      EXPECT_EQ(yb.GetExpr(), y.GetExpr()) << "map " << i << " " << component;
    }
  }
}

void AddVariable(UserData& data, const std::string& name, const std::string& expr) {
  UserVariable v;
  v.name = name;
  v.value_expr = expr;
  data.variable_registry.Add(std::move(v));
}

// Rewrites a current deck into the pre-expression form: every mole-fraction
// entry's ExpressionParam becomes its plain Value, and the owning classes'
// archive versions drop back to the last plain-number version.
void DowngradeMoleFractions(rapidjson::Value& v) {
  if (v.IsArray()) {
    for (auto& e : v.GetArray()) DowngradeMoleFractions(e);
    return;
  }
  if (!v.IsObject()) return;
  for (const char* key : {"Initial_Composition", "Inlet_Mole_Fraction", "Boundary_Mole_Fraction"}) {
    auto it = v.FindMember(key);
    if (it == v.MemberEnd()) continue;
    for (auto& entry : it->value.GetArray()) {
      const double value = entry["value"]["Value"].GetDouble();
      entry["value"].SetDouble(value);
    }
    auto ver = v.FindMember("cereal_class_version");
    if (ver != v.MemberEnd()) {
      // InitialConditions 3 -> 2; IBoundaryStream / PressureBoundaryCondition 5 -> 4.
      const unsigned legacy = std::string(key) == "Initial_Composition" ? 2u : 4u;
      ver->value.SetUint(legacy);
    }
  }
  for (auto& m : v.GetObject()) DowngradeMoleFractions(m.value);
}

}  // namespace

TEST(MoleFractionExpressions, UserVariablesDriveInletAndInitialCompositions) {
  UserData data = LoadPerez();
  AddVariable(data, "y_co2", "0.3");
  for (ExpressionMap* m : AllMoleFractionMaps(data)) {
    ASSERT_EQ(m->size(), 2u);
    m->at("CarbonDioxide").SetExpr("y_co2");
    m->at("Nitrogen").SetExpr("1 - y_co2");
  }

  std::string err;
  ASSERT_TRUE(data.EvaluateAllExpressions(err)) << err;
  for (ExpressionMap* m : AllMoleFractionMaps(data)) {
    EXPECT_DOUBLE_EQ(m->at("CarbonDioxide").Get(), 0.3);
    EXPECT_DOUBLE_EQ(m->at("Nitrogen").Get(), 0.7);
  }
  ASSERT_TRUE(data.PerformChecks(err)) << err;

  // Changing the variable moves every composition with it.
  data.variable_registry.Find("y_co2")->value_expr = "0.45";
  ASSERT_TRUE(data.EvaluateAllExpressions(err)) << err;
  for (ExpressionMap* m : AllMoleFractionMaps(data)) {
    EXPECT_DOUBLE_EQ(m->at("CarbonDioxide").Get(), 0.45);
    EXPECT_DOUBLE_EQ(m->at("Nitrogen").Get(), 0.55);
  }
}

TEST(MoleFractionExpressions, FractionUnitsConvert) {
  UserData data = LoadPerez();
  data.ic.mole_fraction.at("CarbonDioxide").SetExpr("400 [ppm]");
  data.ic.mole_fraction.at("Nitrogen").SetExpr("99.96 [%]");
  std::string err;
  ASSERT_TRUE(data.EvaluateAllExpressions(err)) << err;
  EXPECT_DOUBLE_EQ(data.ic.mole_fraction.at("CarbonDioxide").Get(), 4e-4);
  EXPECT_DOUBLE_EQ(data.ic.mole_fraction.at("Nitrogen").Get(), 0.9996);
}

TEST(MoleFractionExpressions, SumCheckUsesTheEvaluatedValues) {
  UserData data = LoadPerez();
  AddVariable(data, "y_co2", "0.3");
  // Every composition sums to 0.3 + 0.3.
  for (ExpressionMap* m : AllMoleFractionMaps(data)) {
    m->at("CarbonDioxide").SetExpr("y_co2");
    m->at("Nitrogen").SetExpr("y_co2");
  }
  std::string err;
  ASSERT_TRUE(data.EvaluateAllExpressions(err)) << err;
  EXPECT_FALSE(data.PerformChecks(err));
  EXPECT_NE(err.find("Sum of mole fractions"), std::string::npos) << err;
}

TEST(MoleFractionExpressions, ExpressionsSurviveSaveAndReload) {
  UserData data = LoadPerez();
  AddVariable(data, "y_co2", "0.3");
  data.ic.mole_fraction.at("CarbonDioxide").SetExpr("y_co2 / 2");
  data.ic.mole_fraction.at("Nitrogen").SetExpr("1 - y_co2 / 2");
  std::string err;
  ASSERT_TRUE(data.EvaluateAllExpressions(err)) << err;

  const std::string saved = Snapshot(data);
  UserData reloaded;
  LoadDeckFlexible(reloaded, saved);
  EXPECT_EQ(reloaded.ic.mole_fraction.at("CarbonDioxide").GetExpr(), "y_co2 / 2");
  EXPECT_EQ(reloaded.ic.mole_fraction.at("Nitrogen").GetExpr(), "1 - y_co2 / 2");
  ASSERT_TRUE(reloaded.EvaluateAllExpressions(err)) << err;
  EXPECT_DOUBLE_EQ(reloaded.ic.mole_fraction.at("CarbonDioxide").Get(), 0.15);
  ExpectSameMoleFractions(data, reloaded);
}

TEST(MoleFractionExpressions, RegistryListsMoleFractionsAsExpressions) {
  UserData data = LoadPerez();
  AddVariable(data, "y_co2", "0.3");
  ParameterRegistry reg(data);

  size_t endpoint_fractions = 0;
  for (auto& p : reg.All()) {
    if (p.Path().find(".y[") == std::string::npos) continue;
    EXPECT_EQ(p.Info().kind, ParameterKind::Expression) << p.Path();
    EXPECT_TRUE(p.Info().tunable) << p.Path();  // a plain number is tunable
    if (p.Path().rfind("endpoint[", 0) == 0) ++endpoint_fractions;
  }
  EXPECT_GT(endpoint_fractions, 0u);

  reg.Get("initial.y[CarbonDioxide]").SetText("y_co2");
  reg.Get("initial.y[Nitrogen]").SetText("1 - y_co2");
  std::string err;
  ASSERT_TRUE(data.EvaluateAllExpressions(err)) << err;
  EXPECT_DOUBLE_EQ(data.ic.mole_fraction.at("CarbonDioxide").Get(), 0.3);
  EXPECT_DOUBLE_EQ(reg.Get("initial.y[Nitrogen]").GetValue(), 0.7);

  // A formula is no longer a tunable number.
  reg.Refresh();
  EXPECT_FALSE(reg.Get("initial.y[CarbonDioxide]").Info().tunable);
}

TEST(MoleFractionExpressions, PlainNumberDecksStillLoad) {
  UserData data = LoadPerez();
  const std::string current = Snapshot(data);

  rapidjson::Document doc;
  doc.Parse<rapidjson::kParseFullPrecisionFlag | rapidjson::kParseNanAndInfFlag>(current.c_str());
  ASSERT_FALSE(doc.HasParseError());
  DowngradeMoleFractions(doc);
  rapidjson::StringBuffer buf;
  rapidjson::Writer<rapidjson::StringBuffer> writer(buf);
  doc.Accept(writer);
  const std::string legacy = buf.GetString();
  ASSERT_EQ(legacy.find("\"Initial_Composition\":[{\"key\":\"CarbonDioxide\",\"value\":{"), std::string::npos);

  UserData loaded;
  LoadDeckFlexible(loaded, legacy);
  // The deck's mole fractions are plain numbers, whose text is the formatted
  // value, so the legacy load reproduces them exactly.
  ExpectSameMoleFractions(data, loaded);
}

}  // namespace adsim
