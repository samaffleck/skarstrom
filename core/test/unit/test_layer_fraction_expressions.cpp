// test_layer_fraction_expressions.cpp
//
// A layer's mass fractions are ExpressionParams (MixedAdsorbentLayer v5), so a
// blend can be written as x and 1 - x over a deck variable x and varied
// through that variable. These tests pin:
//   * a deck written before the change (plain-number fractions) loads with the
//     same numbers, bit for bit, and re-saves as constant expressions;
//   * expression fractions survive save and reload;
//   * the composition check uses the evaluated values and reports a formula
//     that does not evaluate;
//   * the registry lists the fractions as Expression parameters (text works,
//     as from the CLI and Python);
//   * a run with x and 1 - x gives the same result as the same numbers, and a
//     sensitivity sample over variable[x] changes the blend;
//   * the blend edits refuse a layer whose fractions include a formula.

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <cereal/archives/json.hpp>
#include <cereal/external/rapidjson/document.h>
#include <cereal/external/rapidjson/stringbuffer.h>
#include <cereal/external/rapidjson/writer.h>
#include <gtest/gtest.h>

#include "adsim/objects/deck_edit.h"
#include "adsim/objects/deck_io.h"
#include "adsim/objects/user_data.h"
#include "adsim/parameters/parameter_registry.h"
#include "adsim/simulation.h"
#include "adsim/solver/sensitivity.h"

namespace adsim {
namespace {

const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";
const std::filesystem::path kLegacy =
    std::filesystem::current_path() / "test" / "fixtures" / "legacy" / "Ward_2022_Q=50[ccm]_pre_physics_parity.json";

std::string ReadFile(const std::filesystem::path& p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

std::string Snapshot(UserData& data) {
  std::ostringstream oss;
  {
    cereal::JSONOutputArchive ar(oss);
    ar(data);
  }
  return oss.str();
}

// The layers of a saved deck as JSON text (whole-deck snapshots cannot be
// compared as text: the adsorbent directory's order differs after a reload).
std::string LayersJson(const std::string& deck) {
  rapidjson::Document doc;
  doc.Parse<rapidjson::kParseFullPrecisionFlag>(deck.c_str());
  if (doc.HasParseError()) return "";
  rapidjson::StringBuffer buf;
  rapidjson::Writer<rapidjson::StringBuffer> writer(buf);
  doc["value0"]["Bed"]["Layers"].Accept(writer);
  return buf.GetString();
}

// A plain-number deck variable (value_expr stays empty, so it is tunable).
void AddVariable(UserData& data, const std::string& name, double value) {
  UserVariable v;
  v.name = name;
  v.value = value;
  data.variable_registry.Add(std::move(v));
}

// Perez2019_A with its 13X layer blended with an inert filler (default
// properties, inert isotherms); the fractions are left to the caller.
struct BlendDeck {
  UserData data;
  AdsorbentID zeolite = 0;
  AdsorbentID filler = 0;
  MixedAdsorbentLayer& Layer() { return data.bed.layers_.layers_.front(); }

  BlendDeck() {
    LoadFromJson(data, (kExamples / "Perez2019_A.json").string());
    zeolite = data.bed.layers_.layers_.front().composition_.begin()->first;
    filler = AddAdsorbentToDeck(data, "Filler");
    Layer().composition_[filler] = 0.0;
  }
};

// Rewrites a current deck into the pre-expression form: every layer's
// composition entry becomes its plain Value and the layer's archive version
// drops back to 4.
void DowngradeCompositions(rapidjson::Value& layers) {
  for (auto& layer : layers.GetArray()) {
    for (auto& entry : layer["Composition"].GetArray()) {
      const double value = entry["value"]["Value"].GetDouble();
      entry["value"].SetDouble(value);
    }
    auto ver = layer.FindMember("cereal_class_version");
    if (ver != layer.MemberEnd()) ver->value.SetUint(4u);
  }
}

}  // namespace

// The legacy fixture holds plain-number fractions (MixedAdsorbentLayer < v5).
TEST(LayerFractionExpressions, LegacyDeckLoadsTheSameNumbers) {
  const std::string text = ReadFile(kLegacy);
  ASSERT_FALSE(text.empty()) << kLegacy;
  rapidjson::Document doc;
  doc.Parse<rapidjson::kParseFullPrecisionFlag>(text.c_str());
  ASSERT_FALSE(doc.HasParseError());
  const auto& raw_layers = doc["value0"]["Bed"]["Layers"]["Layers"];
  ASSERT_LT(raw_layers[0]["cereal_class_version"].GetUint(), 5u);

  UserData data;
  LoadDeckFlexible(data, text);
  ASSERT_EQ(data.bed.layers_.layers_.size(), raw_layers.Size());
  for (rapidjson::SizeType i = 0; i < raw_layers.Size(); ++i) {
    const auto& raw = raw_layers[i]["Composition"];
    const auto& layer = data.bed.layers_.layers_[i];
    ASSERT_EQ(layer.composition_.size(), raw.Size());
    for (const auto& entry : raw.GetArray()) {
      const auto id = static_cast<AdsorbentID>(entry["key"].GetInt());
      const ExpressionParam& w = layer.composition_.at(id);
      EXPECT_EQ(w.Get(), entry["value"].GetDouble());
      EXPECT_FALSE(w.HasExpression()) << w.GetExpr();
    }
  }
}

// A blend saved by the previous engine (plain numbers) loads bit-identically,
// evaluates to the same numbers and re-saves as constant expressions.
TEST(LayerFractionExpressions, OldBlendDeckLoadsIdentically) {
  BlendDeck b;
  b.Layer().composition_[b.zeolite] = 0.7;
  b.Layer().composition_[b.filler] = 1.0 - 0.7;  // 0.30000000000000004
  rapidjson::Document doc;
  const std::string current = Snapshot(b.data);
  doc.Parse<rapidjson::kParseFullPrecisionFlag>(current.c_str());
  ASSERT_FALSE(doc.HasParseError());
  DowngradeCompositions(doc["value0"]["Bed"]["Layers"]["Layers"]);
  rapidjson::StringBuffer buf;
  rapidjson::Writer<rapidjson::StringBuffer> writer(buf);
  doc.Accept(writer);

  UserData old;
  LoadDeckFlexible(old, buf.GetString());
  std::string err;
  ASSERT_TRUE(old.EvaluateAllExpressions(err)) << err;
  const auto& w = old.bed.layers_.layers_.front().composition_;
  EXPECT_EQ(w.at(b.zeolite).Get(), 0.7);
  EXPECT_EQ(w.at(b.filler).Get(), 1.0 - 0.7);
  EXPECT_FALSE(w.at(b.filler).HasExpression());
  EXPECT_TRUE(old.PerformChecks(err)) << err;
  // Re-saved, the deck is the current one again.
  EXPECT_EQ(LayersJson(Snapshot(old)), LayersJson(current));
}

TEST(LayerFractionExpressions, ExpressionsSurviveSaveAndReload) {
  BlendDeck b;
  AddVariable(b.data, "x", 0.8);
  b.Layer().composition_[b.zeolite].SetExpr("x");
  b.Layer().composition_[b.filler].SetExpr("1 - x");
  std::string err;
  ASSERT_TRUE(b.data.EvaluateAllExpressions(err)) << err;
  EXPECT_EQ(b.Layer().composition_.at(b.zeolite).Get(), 0.8);
  EXPECT_EQ(b.Layer().composition_.at(b.filler).Get(), 1.0 - 0.8);

  const std::string saved = Snapshot(b.data);
  UserData reloaded;
  LoadDeckFlexible(reloaded, saved);
  const auto& w = reloaded.bed.layers_.layers_.front().composition_;
  EXPECT_EQ(w.at(b.zeolite).GetExpr(), "x");
  EXPECT_EQ(w.at(b.filler).GetExpr(), "1 - x");
  EXPECT_EQ(w.at(b.filler).Get(), 1.0 - 0.8);
  EXPECT_EQ(LayersJson(Snapshot(reloaded)), LayersJson(saved));
}

TEST(LayerFractionExpressions, SumCheckUsesTheEvaluatedValues) {
  BlendDeck b;
  AddVariable(b.data, "x", 0.25);
  b.Layer().composition_[b.zeolite].SetExpr("x");
  b.Layer().composition_[b.filler].SetExpr("1 - x");
  std::string err;
  // x and 1 - x pass for any x in [0, 1].
  for (double x : {0.25, 0.5, 0.9, 1.0}) {
    b.data.variable_registry.Find("x")->value = x;
    ASSERT_TRUE(b.data.EvaluateAllExpressions(err)) << err;
    EXPECT_TRUE(b.data.PerformChecks(err)) << x << ": " << err;
    err.clear();
  }

  // A formula that breaks the sum is caught, and the message shows it.
  b.data.variable_registry.Find("x")->value = 0.25;
  b.Layer().composition_[b.filler].SetExpr("1 - 2 * x");
  ASSERT_TRUE(b.data.EvaluateAllExpressions(err)) << err;
  EXPECT_FALSE(b.data.PerformChecks(err));
  EXPECT_NE(err.find("mass fractions (x = 0.25, 1 - 2 * x = 0.5) sum to 0.75"), std::string::npos) << err;
  EXPECT_NE(err.find("set them to x and 1 - x with a deck variable x"), std::string::npos) << err;

  // x outside [0, 1] makes 1 - x negative.
  err.clear();
  b.Layer().composition_[b.filler].SetExpr("1 - x");
  b.data.variable_registry.Find("x")->value = 1.2;
  ASSERT_TRUE(b.data.EvaluateAllExpressions(err)) << err;
  EXPECT_FALSE(b.data.PerformChecks(err));
  EXPECT_NE(err.find("Mass fractions must be between 0 and 1"), std::string::npos) << err;

  // A formula over an unknown variable is reported, not run on a stale value.
  err.clear();
  b.data.variable_registry.Find("x")->value = 0.5;
  b.Layer().composition_[b.filler].SetExpr("1 - no_such_variable");
  ASSERT_TRUE(b.data.EvaluateAllExpressions(err)) << err;
  EXPECT_FALSE(b.data.PerformChecks(err));
  EXPECT_NE(err.find("the mass fraction '1 - no_such_variable'"), std::string::npos) << err;
  EXPECT_NE(err.find("cannot be evaluated"), std::string::npos) << err;
  EXPECT_FALSE(b.Layer().ValidateComposition());
}

TEST(LayerFractionExpressions, RegistrySetsFractionsAsText) {
  BlendDeck b;
  b.Layer().composition_[b.zeolite] = 1.0;
  AddVariable(b.data, "x", 0.6);
  Simulation sim = Simulation::FromJson(Snapshot(b.data));
  const std::string lp = "layer[" + std::to_string(b.Layer().id_) + "].mass_fraction[";
  const std::string pz = lp + std::to_string(b.zeolite) + "]";
  const std::string pf = lp + std::to_string(b.filler) + "]";

  for (const auto& p : sim.Parameters().All())
    if (p.Path().find(".mass_fraction[") != std::string::npos) {
      EXPECT_EQ(p.Info().kind, ParameterKind::Expression) << p.Path();
      EXPECT_TRUE(p.Info().tunable) << p.Path();  // a plain number stays tunable
    }

  sim.SetParameterText(pz, "x");
  sim.SetParameterText(pf, "1 - x");
  EXPECT_EQ(sim.GetParameter(pz), 0.6);
  EXPECT_EQ(sim.GetParameter(pf), 1.0 - 0.6);
  EXPECT_EQ(sim.GetParameterText(pf), "1 - x");
  EXPECT_FALSE(sim.Parameters().Get(pf).Info().tunable);  // vary variable[x] instead
  std::string err;
  EXPECT_TRUE(sim.Check(err)) << err;

  // Sweeping the variable moves both fractions together.
  for (double x : {0.2, 0.5, 0.75}) {
    sim.SetParameter("variable[x]", x);
    EXPECT_EQ(sim.GetParameter(pz), x);
    EXPECT_EQ(sim.GetParameter(pf), 1.0 - x);
    EXPECT_TRUE(sim.Check(err)) << err;
  }

  // An expression over an unknown name is refused and the old text kept.
  EXPECT_THROW(sim.SetParameterText(pf, "1 - zz"), std::invalid_argument);
  EXPECT_EQ(sim.GetParameterText(pf), "1 - x");
}

// x and 1 - x run exactly like the same two numbers; a sensitivity sample over
// variable[x] changes the result.
TEST(LayerFractionExpressions, VariableBlendRunsLikeNumbers) {
  BlendDeck b;
  // Short and coarse (fast tier): one cycle, 20 cells.
  b.Layer().domain_.n_cells = 20;
  b.data.cycle.number_of_cycles_ = 1;

  const double x = 0.8;
  b.Layer().composition_[b.zeolite] = x;
  b.Layer().composition_[b.filler] = 1.0 - x;
  const std::string numeric = Snapshot(b.data);

  AddVariable(b.data, "x", 0.8);
  b.Layer().composition_[b.zeolite].SetExpr("x");
  b.Layer().composition_[b.filler].SetExpr("1 - x");
  const std::string formula = Snapshot(b.data);

  auto outlet = [](const std::string& deck) {
    Simulation sim = Simulation::FromJson(deck);
    const RunResult r = sim.Run();
    EXPECT_TRUE(r.status == RunStatus::Completed || r.status == RunStatus::ReachedCSS) << r.message;
    return sim.Series("mole_fraction", "CarbonDioxide").values;
  };
  const std::vector<double> a = outlet(numeric);
  const std::vector<double> c = outlet(formula);
  ASSERT_FALSE(a.empty());
  ASSERT_EQ(a.size(), c.size());
  for (size_t i = 0; i < a.size(); ++i) ASSERT_EQ(a[i], c[i]) << i;

  // The study path: variable[x] as an input of EvaluateSample.
  using namespace sensitivity;
  const std::vector<Input> inputs = {Uniform("variable[x]", 0.5, 0.9)};
  const std::vector<Output> outputs = {{"CarbonDioxide", "purity"}};
  RunPolicy policy;
  policy.require_css = false;
  policy.cycles = 1;
  const SampleOutcome lo = EvaluateSample(formula, inputs, {0.5}, outputs, policy);
  const SampleOutcome hi = EvaluateSample(formula, inputs, {0.9}, outputs, policy);
  ASSERT_TRUE(lo.accepted) << lo.message;
  ASSERT_TRUE(hi.accepted) << hi.message;
  EXPECT_NE(lo.outputs[0], hi.outputs[0]);
}

TEST(LayerFractionExpressions, BlendEditsRefuseFormulaFractions) {
  BlendDeck b;
  AddVariable(b.data, "x", 0.6);
  b.Layer().composition_[b.zeolite].SetExpr("x");
  b.Layer().composition_[b.filler].SetExpr("1 - x");
  std::string err;
  ASSERT_TRUE(b.data.EvaluateAllExpressions(err)) << err;
  const AdsorbentID third = AddAdsorbentToDeck(b.data, "Third");

  auto unchanged = [&] {
    ASSERT_EQ(b.Layer().composition_.size(), 2u);
    EXPECT_EQ(b.Layer().composition_.at(b.zeolite).GetExpr(), "x");
    EXPECT_EQ(b.Layer().composition_.at(b.filler).GetExpr(), "1 - x");
  };
  EXPECT_FALSE(AddSolidToLayer(b.data, b.Layer(), third, 0.1, err));
  EXPECT_NE(err.find("include expressions (x, 1 - x)"), std::string::npos) << err;
  unchanged();
  err.clear();
  EXPECT_FALSE(SetLayerMassFraction(b.data, b.Layer(), b.zeolite, 0.5, err));
  EXPECT_NE(err.find("cannot be rescaled"), std::string::npos) << err;
  unchanged();
  err.clear();
  EXPECT_FALSE(RemoveSolidFromLayer(b.data, b.Layer(), b.filler, err));
  EXPECT_NE(err.find("cannot be rescaled"), std::string::npos) << err;
  unchanged();

  // Numbers again: the edits rescale as before.
  b.Layer().composition_[b.zeolite] = 0.6;
  b.Layer().composition_[b.filler] = 0.4;
  ASSERT_TRUE(SetLayerMassFraction(b.data, b.Layer(), b.zeolite, 0.5, err)) << err;
  EXPECT_DOUBLE_EQ(b.Layer().composition_.at(b.filler).Get(), 0.5);
  ASSERT_TRUE(AddSolidToLayer(b.data, b.Layer(), third, 0.2, err)) << err;
  EXPECT_DOUBLE_EQ(b.Layer().composition_.at(b.zeolite).Get(), 0.4);
  ASSERT_TRUE(RemoveSolidFromLayer(b.data, b.Layer(), third, err)) << err;
  EXPECT_DOUBLE_EQ(b.Layer().composition_.at(b.zeolite).Get() + b.Layer().composition_.at(b.filler).Get(), 1.0);
}

}  // namespace adsim
