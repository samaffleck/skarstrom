// test_parameter_registry.cpp
//
// The parameter registry (adsim/parameters/parameter_registry.h) is the single
// engine-side list of settable deck values. These tests pin:
//   * every example deck enumerates unique, self-resolving paths whose
//     accessors read and write the live model;
//   * the optimiser's AvailableParameters lists exactly the registry's tunable
//     parameters, backed by the same model values;
//   * renaming an object changes a parameter's alias but not its path.

#include <filesystem>
#include <set>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/objects/user_data.h"
#include "adsim/parameters/parameter_registry.h"
#include "adsim/solver/optimisation.h"


namespace adsim {

namespace {

const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";

std::vector<std::filesystem::path> ExampleDecks() {
  std::vector<std::filesystem::path> out;
  for (const auto& e : std::filesystem::directory_iterator(kExamples))
    if (e.path().extension() == ".json") out.push_back(e.path());
  std::sort(out.begin(), out.end());
  return out;
}

}  // namespace

TEST(ParameterRegistry, ExampleDecksHaveUniqueSelfResolvingPaths) {
  const auto decks = ExampleDecks();
  ASSERT_GT(decks.size(), 10u);
  for (const auto& deck : decks) {
    SCOPED_TRACE(deck.filename().string());
    UserData data;
    LoadFromJson(data, deck.string());
    ParameterRegistry reg(data);
    ASSERT_GT(reg.Size(), StaticFields().size());

    std::set<std::string> paths;
    std::map<std::string, int> alias_count;
    for (const auto& p : reg.All()) {
      EXPECT_TRUE(paths.insert(p.Path()).second) << "duplicate path " << p.Path();
      EXPECT_EQ(p.Path().find(" #"), std::string::npos) << "path needed a uniqueness suffix: " << p.Path();
      alias_count[p.Info().alias]++;
      EXPECT_FALSE(p.Info().name.empty()) << p.Path();
      if (p.Info().tunable) {
        EXPECT_TRUE(p.Info().kind == ParameterKind::Real || p.Info().kind == ParameterKind::Expression);
        EXPECT_LE(p.Info().lower, p.Info().upper) << p.Path();
        EXPECT_NE(p.ValuePointer(), nullptr) << p.Path();
      }
    }
    for (auto& p : reg.All()) {
      EXPECT_EQ(reg.Find(p.Path()), &p) << p.Path();
      if (alias_count[p.Info().alias] == 1) EXPECT_EQ(reg.Find(p.Info().alias), &p) << p.Info().alias;
      else EXPECT_THROW(reg.Find(p.Info().alias), std::invalid_argument) << p.Info().alias;
    }
    EXPECT_EQ(reg.Find("no.such.path"), nullptr);
    EXPECT_THROW(reg.Get("no.such.path"), std::out_of_range);
  }
}

TEST(ParameterRegistry, AccessorsReadAndWriteTheModel) {
  UserData data;
  LoadFromJson(data, (kExamples / "Perez2019_A.json").string());
  ParameterRegistry reg(data);

  Parameter& reltol = reg.Get("settings.reltol");
  EXPECT_EQ(reltol.Info().kind, ParameterKind::Expression);
  EXPECT_FALSE(reltol.Info().tunable);
  reltol.SetValue(3e-7);
  EXPECT_DOUBLE_EQ(data.settings.reltol.Get(), 3e-7);
  EXPECT_EQ(data.settings.reltol.GetExpr(), ExpressionParam::FormatScalar(3e-7));

  Parameter& order = reg.Get("settings.max_order");
  EXPECT_EQ(order.Info().kind, ParameterKind::Integer);
  order.SetValue(3.4);
  EXPECT_EQ(data.settings.max_order, 3);
  order.SetText("4");
  EXPECT_EQ(data.settings.max_order, 4);
  EXPECT_THROW(order.SetText("four"), std::invalid_argument);

  Parameter& title = reg.Get("info.title");
  title.SetText("A title");
  EXPECT_EQ(data.info.title, "A title");
  EXPECT_THROW(title.GetValue(), std::logic_error);

  Parameter& cycles = reg.Get("cycle.number_of_cycles");
  cycles.SetValue(7);
  EXPECT_EQ(data.cycle.number_of_cycles_, 7);

  Parameter& t = reg.Get("step[1].time");
  EXPECT_TRUE(t.Info().tunable);
  EXPECT_EQ(t.Info().unit, "s");
  t.SetValue(42.0);
  EXPECT_DOUBLE_EQ(data.cycle.step_list_[0].step_time_.Get(), 42.0);
  t.SetText("2*21");
  std::string eval_err;
  ASSERT_TRUE(data.EvaluateAllExpressions(eval_err)) << eval_err;
  EXPECT_DOUBLE_EQ(t.GetValue(), 42.0);
  reg.Refresh();
  EXPECT_FALSE(reg.Get("step[1].time").Info().tunable) << "a formula must not be offered to the optimiser";
}

// AvailableParameters (the optimiser's list) is the registry's tunable subset.
TEST(ParameterRegistry, AvailableParametersAreTheTunableParameters) {
  for (const auto& deck : ExampleDecks()) {
    SCOPED_TRACE(deck.filename().string());
    UserData data;
    LoadFromJson(data, deck.string());
    ParameterRegistry reg(data);
    AvailableParameters ap;
    ap.UpdateParameterList(data);
    const auto tunable = reg.Tunable();
    ASSERT_EQ(ap.Size(), tunable.size());
    for (size_t i = 0; i < tunable.size(); ++i) {
      EXPECT_EQ(ap[i].key, tunable[i]->Path());
      EXPECT_EQ(&ap[i].value_ref, tunable[i]->ValuePointer()) << tunable[i]->Path();
      EXPECT_EQ(ap[i].expr_param, tunable[i]->ExpressionPointer()) << tunable[i]->Path();
      EXPECT_EQ(ap[i].name, tunable[i]->Info().name);
      EXPECT_EQ(ap[i].units, tunable[i]->Info().unit);
      EXPECT_DOUBLE_EQ(ap[i].lower_bound, tunable[i]->Info().lower);
      EXPECT_DOUBLE_EQ(ap[i].upper_bound, tunable[i]->Info().upper);
    }
  }
}

// Renaming objects changes aliases, never paths.
TEST(ParameterRegistry, RenamesKeepPaths) {
  UserData data;
  LoadFromJson(data, (kExamples / "Perez2019_A.json").string());
  ParameterRegistry reg(data);
  std::vector<std::string> before;
  for (const auto& p : reg.All()) before.push_back(p.Path());

  data.RenameStep(0, "Adsorb");
  ASSERT_FALSE(data.bed.layers_.layers_.empty());
  data.bed.layers_.layers_.front().name_ = "Renamed layer";
  for (auto& [id, ads] : const_cast<std::unordered_map<AdsorbentID, Adsorbent>&>(data.adsorbent_directory.GetAll()))
    ads.name_ += " (renamed)";
  reg.Refresh();
  std::vector<std::string> after;
  for (const auto& p : reg.All()) after.push_back(p.Path());
  EXPECT_EQ(before, after);
  EXPECT_NE(reg.Find("step[Adsorb].time"), nullptr);
  EXPECT_EQ(reg.Find("step[Feed].time"), nullptr);
  EXPECT_EQ(reg.Find("step[Adsorb].time"), reg.Find("step[1].time"));
}

}  // namespace adsim
