/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// Reactions in a deck (adsim/objects/reaction.h): the builder, serialisation
// (round trip, old v4 adsorbents and decks still load), the parameter
// registry, the component-removal cascade, the deck checks and warnings, and
// the reaction edit cascades of deck_edit.h.

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>

#include <cereal/archives/json.hpp>

#include "adsim/objects/deck_builder.h"
#include "adsim/objects/deck_edit.h"
#include "adsim/objects/reaction.h"
#include "adsim/parameters/parameter_registry.h"
#include "adsim/simulation.h"
#include "adsim/solver/deck_warnings.h"

namespace adsim {
namespace {

const std::string CH4 = "Methane", H2O = "Water", CO = "CarbonMonoxide", H2 = "Hydrogen", CO2 = "CarbonDioxide";

std::string ThrowMessage(const std::function<void()>& f) {
  try {
    f();
  } catch (const std::invalid_argument& e) {
    return e.what();
  }
  return "";
}

struct Reformer {
  Simulation sim;
  AdsorbentID ni = 0, sorbent = 0;
  ReactionID smr = 0, wgs = 0;
  size_t layer = 0;
};

// A steam-reforming catalyst mixed with a CO2 sorbent in one layer.
void Build(Reformer& b) {
  UserData& d = b.sim.Data();
  for (const auto& c : {CH4, H2O, CO, H2, CO2}) build::AddComponent(d, c);
  b.ni = build::AddAdsorbent(d, "Ni catalyst");
  b.sorbent = build::AddAdsorbent(d, "Sorbent");
  b.layer = build::AddLayer(d, "Layer 1");
  build::SetLayerComposition(d, b.layer, {{b.ni, 0.3}, {b.sorbent, 0.7}});
  b.smr = build::AddReaction(d, b.ni, "SMR", {{CH4, -1}, {H2O, -1}, {CO, 1}, {H2, 3}}, "lhhw");
  b.wgs = build::AddReaction(d, b.ni, "WGS", {{CO, -1}, {H2O, -1}, {CO2, 1}, {H2, 1}}, "reversible_power_law");
  EXPECT_EQ(build::AddAdsorptionTerm(d, b.ni, b.smr), 1u);
  EXPECT_EQ(build::AddAdsorptionTerm(d, b.ni, b.smr), 2u);
  b.sim.RefreshParameters();
}

// A reformer deck (three Xu-Froment LHHW rates, one four-term denominator
// squared) saved by the engine before LHHW rates could have several
// denominators. The working directory of the tests is core/.
std::string ReadFixture() {
  std::ifstream in(std::filesystem::current_path() / "test" / "fixtures" / "legacy" / "lhhw_v1_reformer.json");
  std::stringstream ss;
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

TEST(ReactionDeck, BuilderCreatesReactionsAndRejectsBadInput) {
  Reformer b;
  Build(b);
  UserData& d = b.sim.Data();
  const auto& src = d.adsorbent_directory.Get(b.ni).reaction_source_;
  ASSERT_EQ(src.reactions_.size(), 2u);
  EXPECT_EQ(b.smr, 1u);
  EXPECT_EQ(b.wgs, 2u);
  EXPECT_EQ(src.Get(b.smr).rate_->GetType(), ReactionRateType::Lhhw);
  EXPECT_EQ(dynamic_cast<LhhwRate&>(*src.Get(b.smr).rate_).denominators_.size(), 1u);
  EXPECT_EQ(dynamic_cast<LhhwRate&>(*src.Get(b.smr).rate_).denominators_[0].adsorption_.size(), 2u);

  EXPECT_NE(ThrowMessage([&] { build::AddReaction(d, b.ni, "X", {{CH4, -1}}, "michaelis"); }).find("power_law"),
            std::string::npos);
  EXPECT_NE(ThrowMessage([&] { build::AddReaction(d, b.ni, "X", {{"Argon", -1}}); }).find("Argon"),
            std::string::npos);
  EXPECT_NE(ThrowMessage([&] { build::AddAdsorptionTerm(d, b.ni, b.wgs); }).find("lhhw"), std::string::npos);
  EXPECT_NE(ThrowMessage([&] { build::SetReactionKind(d, b.ni, 99, "lhhw"); }).find("99"), std::string::npos);
  EXPECT_NE(ThrowMessage([&] { build::SetEffectiveness(d, b.ni, b.smr, "thiele", H2); }).find("key reactant"),
            std::string::npos);

  build::SetEffectiveness(d, b.ni, b.smr, "thiele", CH4);
  EXPECT_EQ(src.Get(b.smr).effectiveness_.mode_, EffectivenessMode::Thiele);
  build::SetReactionBasis(d, b.ni, b.wgs, "concentration");
  EXPECT_EQ(src.Get(b.wgs).basis_, RateBasis::Concentration);
  build::RemoveAdsorptionTerm(d, b.ni, b.smr, 1);
  EXPECT_EQ(dynamic_cast<LhhwRate&>(*src.Get(b.smr).rate_).denominators_[0].adsorption_.size(), 1u);

  // Denominators: a second one, its terms, and the errors.
  EXPECT_EQ(build::AddDenominator(d, b.ni, b.smr), 2u);
  EXPECT_EQ(build::AddAdsorptionTerm(d, b.ni, b.smr, 2), 1u);
  EXPECT_EQ(build::AddAdsorptionTerm(d, b.ni, b.smr, 2), 2u);
  const auto& dens = dynamic_cast<LhhwRate&>(*src.Get(b.smr).rate_).denominators_;
  ASSERT_EQ(dens.size(), 2u);
  EXPECT_EQ(dens[1].adsorption_.size(), 2u);
  EXPECT_EQ(dens[1].adsorption_[0].exponent_.size(), 5u);  // every component, at 0
  EXPECT_DOUBLE_EQ(dens[1].exponent_.Get(), 1.0);
  EXPECT_NE(ThrowMessage([&] { build::AddAdsorptionTerm(d, b.ni, b.smr, 3); }).find("no denominator 3"),
            std::string::npos);
  EXPECT_NE(ThrowMessage([&] { build::RemoveAdsorptionTerm(d, b.ni, b.smr, 3, 2); }).find("denominator 2"),
            std::string::npos);
  EXPECT_NE(ThrowMessage([&] { build::AddDenominator(d, b.ni, b.wgs); }).find("lhhw"), std::string::npos);
  build::RemoveAdsorptionTerm(d, b.ni, b.smr, 1, 2);
  EXPECT_EQ(dens[1].adsorption_.size(), 1u);
  build::RemoveDenominator(d, b.ni, b.smr, 1);  // the second moves down
  ASSERT_EQ(dens.size(), 1u);
  EXPECT_EQ(dens[0].adsorption_.size(), 1u);
  EXPECT_NE(ThrowMessage([&] { build::RemoveDenominator(d, b.ni, b.smr, 1); }).find("at least one"),
            std::string::npos);
  EXPECT_NE(ThrowMessage([&] { build::RemoveDenominator(d, b.ni, b.smr, 2); }).find("no denominator 2"),
            std::string::npos);
  build::RemoveReaction(d, b.ni, b.wgs);
  EXPECT_EQ(src.reactions_.size(), 1u);
  // IDs are not reused.
  EXPECT_EQ(build::AddReaction(d, b.ni, "WGS2", {{CO, -1}, {H2O, -1}, {CO2, 1}, {H2, 1}}), 3u);

  const auto kinds = build::ReactionRateKinds();
  EXPECT_EQ(kinds, (std::vector<std::string>{"power_law", "reversible_power_law", "lhhw", "expression"}));
}

TEST(ReactionDeck, RegistryPathsAndAliases) {
  Reformer b;
  Build(b);
  auto& reg = b.sim.Parameters();
  const std::string head = "adsorbent[" + std::to_string(b.ni) + "].reaction[1].";
  for (const char* leaf : {"nu[Methane]", "basis", "rate.pre_exp_factor", "rate.activation_energy",
                           "rate.order[Water]", "rate.order[Hydrogen]", "rate.equilibrium_pre_exp_factor",
                           "rate.equilibrium_energy", "rate.denominator[1].exponent", "rate.reversible",
                           "rate.denominator[1].adsorption[1].pre_exp_factor",
                           "rate.denominator[1].adsorption[1].energy",
                           "rate.denominator[1].adsorption[2].exponent[Hydrogen]", "heat_of_reaction", "heat_mode",
                           "heat_reference_temperature", "effectiveness"}) {
    EXPECT_NE(reg.Find(head + leaf), nullptr) << head + leaf;
  }
  EXPECT_EQ(reg.Find(head + "tortuosity"), nullptr);  // constant effectiveness

  Parameter* order = reg.Find("adsorbent[Ni catalyst].reaction[SMR].rate.order[Water]");
  ASSERT_NE(order, nullptr);
  EXPECT_EQ(order->Path(), head + "rate.order[Water]");
  EXPECT_EQ(order->Info().category, "reaction");
  EXPECT_EQ(order->Info().subcategory, "Layer 1, Ni catalyst, SMR");
  EXPECT_TRUE(order->Info().tunable);
  EXPECT_DOUBLE_EQ(order->GetValue(), 1.0);

  const Parameter& nu = reg.Get(head + "nu[Hydrogen]");
  EXPECT_DOUBLE_EQ(nu.GetValue(), 3.0);
  EXPECT_FALSE(nu.Info().tunable);

  const Parameter& k = reg.Get(head + "rate.pre_exp_factor");
  EXPECT_TRUE(k.Info().tunable);
  EXPECT_FALSE(k.Info().unit_checked);
  EXPECT_EQ(k.Info().unit, "mol/(kg-s-bar^2)");  // orders CH4 1 + H2O 1

  const Parameter& basis = reg.Get(head + "basis");
  EXPECT_EQ(basis.Info().choices, (std::vector<std::string>{"partial_pressure", "concentration"}));
  EXPECT_EQ(reg.Get(head + "heat_mode").Info().choices, (std::vector<std::string>{"constant", "kirchhoff"}));
  EXPECT_DOUBLE_EQ(reg.Get(head + "heat_mode").GetValue(), 1.0);  // kirchhoff

  // Values set through the registry reach the model.
  reg.Get(head + "rate.order[Hydrogen]").SetValue(-2.5);
  reg.Get(head + "rate.denominator[1].adsorption[2].exponent[Hydrogen]").SetValue(-1.0);
  const auto& smr = b.sim.Data().adsorbent_directory.Get(b.ni).reaction_source_.Get(b.smr);
  const auto& l = dynamic_cast<const LhhwRate&>(*smr.rate_);
  EXPECT_DOUBLE_EQ(l.order_.at(H2).Get(), -2.5);
  EXPECT_DOUBLE_EQ(l.denominators_[0].adsorption_[1].exponent_.at(H2).Get(), -1.0);

  // The paths of the single-denominator layout still resolve, to denominator 1
  // (by ID and by name), and are not listed.
  for (const std::string& h : {head, std::string("adsorbent[Ni catalyst].reaction[SMR].")}) {
    EXPECT_EQ(reg.Find(h + "rate.denominator_exponent"), reg.Find(head + "rate.denominator[1].exponent"));
    EXPECT_EQ(reg.Find(h + "rate.adsorption[2].exponent[Hydrogen]"),
              reg.Find(head + "rate.denominator[1].adsorption[2].exponent[Hydrogen]"));
    EXPECT_EQ(reg.Find(h + "rate.adsorption[1].energy"), reg.Find(head + "rate.denominator[1].adsorption[1].energy"));
  }
  EXPECT_NE(reg.Find(head + "rate.denominator_exponent"), nullptr);
  EXPECT_EQ(reg.Find(head + "rate.adsorption[3].energy"), nullptr);
  for (const auto& p : reg.All()) {
    EXPECT_EQ(p.Path().find("denominator_exponent"), std::string::npos) << p.Path();
    EXPECT_EQ(p.Path().find("rate.adsorption["), std::string::npos) << p.Path();
  }
  reg.Get(head + "rate.denominator_exponent").SetValue(2.0);
  EXPECT_DOUBLE_EQ(l.denominators_[0].exponent_.Get(), 2.0);
  EXPECT_EQ(reg.Get(head + "rate.denominator[1].adsorption[1].energy").Info().name, "adsorption[1] energy");

  // A second denominator: its exponent and terms under denominator[2].
  build::AddDenominator(b.sim.Data(), b.ni, b.smr);
  build::AddAdsorptionTerm(b.sim.Data(), b.ni, b.smr, 2);
  b.sim.RefreshParameters();
  auto& reg2 = b.sim.Parameters();
  const Parameter& m2 = reg2.Get(head + "rate.denominator[2].exponent");
  EXPECT_DOUBLE_EQ(m2.GetValue(), 1.0);
  EXPECT_EQ(m2.Info().name, "denominator[2] exponent");
  reg2.Get(head + "rate.denominator[2].adsorption[1].exponent[Water]").SetValue(1.0);
  reg2.Get(head + "rate.denominator[2].adsorption[1].exponent[Hydrogen]").SetValue(-1.0);
  b.sim.RefreshParameters();
  const Parameter& K2 = b.sim.Parameters().Get(head + "rate.denominator[2].adsorption[1].pre_exp_factor");
  EXPECT_EQ(K2.Info().name, "denominator[2] adsorption[1] pre-exponential factor");
  EXPECT_EQ(K2.Info().unit, "-");  // exponents sum to 0
  EXPECT_TRUE(K2.Info().tunable);
  EXPECT_NE(b.sim.Parameters().Find("adsorbent[Ni catalyst].reaction[SMR].rate.denominator[2].adsorption[1].energy"),
            nullptr);

  // Thiele mode lists the tortuosity instead.
  build::SetEffectiveness(b.sim.Data(), b.ni, b.smr, "thiele", CH4);
  b.sim.RefreshParameters();
  EXPECT_NE(b.sim.Parameters().Find(head + "tortuosity"), nullptr);
  EXPECT_EQ(b.sim.Parameters().Find(head + "effectiveness"), nullptr);
}

TEST(ReactionDeck, SerialisationRoundTrip) {
  Reformer b;
  Build(b);
  UserData& d = b.sim.Data();
  auto& reg = b.sim.Parameters();
  const std::string head = "adsorbent[" + std::to_string(b.ni) + "].reaction[1].";
  reg.Get(head + "rate.order[Hydrogen]").SetText("-2.5");
  reg.Get(head + "heat_of_reaction").SetValue(206100.0);
  reg.Get(head + "rate.reversible").SetValue(0.0);
  build::AddDenominator(d, b.ni, b.smr);
  build::AddAdsorptionTerm(d, b.ni, b.smr, 2);
  b.sim.RefreshParameters();
  b.sim.Parameters().Get(head + "rate.denominator[2].exponent").SetText("0.5");
  b.sim.Parameters().Get(head + "rate.denominator[2].adsorption[1].energy").SetValue(-6e4);
  build::SetEffectiveness(d, b.ni, b.wgs, "thiele", CO);
  build::SetReactionBasis(d, b.ni, b.wgs, "concentration");
  std::string err;
  ASSERT_TRUE(d.EvaluateAllExpressions(err)) << err;

  // Through this test's own archive and through the library's (ToJson).
  for (const std::string& json : {Snapshot(d), b.sim.ToJson()}) {
    EXPECT_NE(json.find("\"Reaction_Source\""), std::string::npos);
    UserData back;
    ASSERT_NO_THROW(LoadFromJsonText(back, json));
    const auto& src = back.adsorbent_directory.Get(b.ni).reaction_source_;
    ASSERT_EQ(src.reactions_.size(), 2u);
    EXPECT_EQ(src.NextId(), 2u);
    const Reaction& smr = src.Get(b.smr);
    EXPECT_EQ(smr.name_, "SMR");
    EXPECT_EQ(smr.nu_, (std::map<std::string, double>{{CH4, -1}, {H2O, -1}, {CO, 1}, {H2, 3}}));
    const auto& l = dynamic_cast<const LhhwRate&>(*smr.rate_);
    EXPECT_EQ(l.order_.at(H2).GetExpr(), "-2.5");
    EXPECT_FALSE(l.reversible_);
    ASSERT_EQ(l.denominators_.size(), 2u);
    EXPECT_EQ(l.denominators_[0].adsorption_.size(), 2u);
    EXPECT_EQ(l.denominators_[1].adsorption_.size(), 1u);
    EXPECT_EQ(l.denominators_[1].exponent_.GetExpr(), "0.5");
    EXPECT_DOUBLE_EQ(l.denominators_[1].adsorption_[0].energy_.Get(), -6e4);
    EXPECT_DOUBLE_EQ(smr.heat_of_reaction_.Get(), 206100.0);
    EXPECT_EQ(smr.heat_mode_, HeatOfReactionMode::Kirchhoff);
    const Reaction& wgs = src.Get(b.wgs);
    EXPECT_EQ(wgs.rate_->GetType(), ReactionRateType::ReversiblePowerLaw);
    EXPECT_EQ(wgs.basis_, RateBasis::Concentration);
    EXPECT_EQ(wgs.effectiveness_.mode_, EffectivenessMode::Thiele);
    EXPECT_EQ(wgs.effectiveness_.key_component_, CO);
    // The plain sorbent has no reactions.
    EXPECT_TRUE(back.adsorbent_directory.Get(b.sorbent).reaction_source_.Empty());
    // Saving the reloaded reactions reproduces their text (every reaction
    // map is ordered, so the archive is deterministic).
    auto text = [](const ReactionSource& rs) {
      std::ostringstream o;
      {
        cereal::JSONOutputArchive ar(o);
        ar(rs);
      }
      return o.str();
    };
    EXPECT_EQ(text(src), text(d.adsorbent_directory.Get(b.ni).reaction_source_));
  }
}

// An LHHW rate saved before it could have several denominators (LhhwRate
// version 1: Denominator_Exponent and Adsorption_Terms) loads as one
// denominator with the same values; it saves as version 2 (Denominators).
TEST(ReactionDeck, VersionOneLhhwRateLoadsAsOneDenominator) {
  const std::string v1 = ReadFixture();
  ASSERT_NE(v1.find("\"Denominator_Exponent\""), std::string::npos);
  UserData d;
  ASSERT_NO_THROW(LoadFromJsonText(d, v1));
  bool found = false;
  for (auto& [id, ads] : d.adsorbent_directory.GetAll())
    for (auto& r : ads.reaction_source_.reactions_) {
      const auto& l = dynamic_cast<const LhhwRate&>(*r.rate_);
      found = true;
      ASSERT_EQ(l.denominators_.size(), 1u) << r.name_;
      EXPECT_DOUBLE_EQ(l.denominators_[0].exponent_.Get(), 2.0);
      ASSERT_EQ(l.denominators_[0].adsorption_.size(), 4u);
      // Xu-Froment K_H2O p_H2O / p_H2 (the fourth term).
      const auto& t = l.denominators_[0].adsorption_[3];
      EXPECT_DOUBLE_EQ(t.pre_exp_factor_.Get(), 1.77e5);
      EXPECT_DOUBLE_EQ(t.energy_.Get(), 88.68e3);
      EXPECT_DOUBLE_EQ(t.exponent_.at(H2O).Get(), 1.0);
      EXPECT_DOUBLE_EQ(t.exponent_.at(H2).Get(), -1.0);
    }
  EXPECT_TRUE(found);
  const std::string v2 = Snapshot(d);
  EXPECT_EQ(v2.find("\"Denominator_Exponent\""), std::string::npos);
  EXPECT_NE(v2.find("\"Denominators\""), std::string::npos);
  UserData back;
  ASSERT_NO_THROW(LoadFromJsonText(back, v2));
  for (auto& [id, ads] : back.adsorbent_directory.GetAll())
    for (auto& r : ads.reaction_source_.reactions_) {
      const auto& l = dynamic_cast<const LhhwRate&>(*r.rate_);
      ASSERT_EQ(l.denominators_.size(), 1u);
      EXPECT_EQ(l.denominators_[0].adsorption_.size(), 4u);
    }
}

TEST(ReactionDeck, AdsorbentArchiveIsVersionFive) {
  Adsorbent a("x");
  std::ostringstream os;
  {
    cereal::JSONOutputArchive ar(os);
    ar(a);
  }
  const std::string json = os.str();
  EXPECT_NE(json.find("\"cereal_class_version\": 5"), std::string::npos) << json.substr(0, 200);
  EXPECT_NE(json.find("\"Reaction_Source\""), std::string::npos);
}

TEST(ReactionDeck, VersionFourAdsorbentsAndDecksLoadWithoutReactions) {
  namespace fs = std::filesystem;
  // A v4 adsorbent archive (no Reaction_Source node), as every library file
  // and deck written before reactions holds. Loading it into an adsorbent that
  // had reactions leaves none.
  const std::string v4 = R"({
    "value0": {
        "cereal_class_version": 4,
        "Name": "old",
        "Density": {"cereal_class_version": 2, "Value": 1100.0, "Expression": "1100", "DisplayUnitIndex": 0},
        "Specific_Heat_Capacity": {"Value": 900.0, "Expression": "900", "DisplayUnitIndex": 0},
        "Thermal_Conductivity": {"Value": 0.5, "Expression": "0.5", "DisplayUnitIndex": 0},
        "Adsorption_Source": {"cereal_class_version": 1, "Mass_Transfer_Models": [], "Heat_Of_Adsorption": [],
                              "Isotherms": [], "Isotherm_Map": [], "Mixing_Rule": 0},
        "Particle_Void_Fraction": {"Value": 0.35, "Expression": "0.35", "DisplayUnitIndex": 0}
    }
})";
  Adsorbent old;
  old.reaction_source_.Add("stale", {{CH4, -1}, {H2, 2}});
  std::string err;
  ASSERT_TRUE(LoadAdsorbentFromString(v4, old, err)) << err;
  EXPECT_TRUE(old.reaction_source_.Empty());
  EXPECT_EQ(old.name_, "old");
  EXPECT_DOUBLE_EQ(old.solid_density_.Get(), 1100.0);

  // Library files and an example deck still load, with no reactions.
  const fs::path lib = fs::path("..") / "data" / "adsorbents" / "zeolite_13x_perez2019.json";
  if (fs::exists(lib)) {
    std::ifstream in(lib);
    std::stringstream ss;
    ss << in.rdbuf();
    Adsorbent a;
    ASSERT_TRUE(LoadAdsorbentFromString(ss.str(), a, err)) << err;
    EXPECT_TRUE(a.reaction_source_.Empty());
  }
  const fs::path decks = fs::path("..") / "examples" / "decks";
  if (fs::exists(decks)) {
    for (const auto& entry : fs::directory_iterator(decks)) {
      if (entry.path().extension() != ".json") continue;
      UserData d;
      ASSERT_NO_THROW(LoadFromJson(d, entry.path().string())) << entry.path();
      for (const auto& [id, ads] : d.adsorbent_directory.GetAll()) EXPECT_TRUE(ads.reaction_source_.Empty());
      break;  // one deck is enough here (DeckRoundTrip covers them all)
    }
  }
}

TEST(ReactionDeck, RemovingAComponentUsedByAReactionIsRefused) {
  Reformer b;
  Build(b);
  UserData& d = b.sim.Data();
  std::string err;
  EXPECT_FALSE(RemoveComponentFromDeck(d, H2O, err));
  EXPECT_EQ(err, "Water is used by reactions SMR, WGS (on Ni catalyst); remove them first");
  EXPECT_NE(std::find(d.fluid.components_.begin(), d.fluid.components_.end(), H2O), d.fluid.components_.end());

  // A component only in zero-valued entries (an LHHW exponent of 0) is removed
  // and its entries dropped. Add Nitrogen: AddAdsorptionTerm seeds a zero
  // exponent for every deck component.
  build::AddComponent(d, "Nitrogen");
  build::AddAdsorptionTerm(d, b.ni, b.smr);
  const auto& l = dynamic_cast<const LhhwRate&>(*d.adsorbent_directory.Get(b.ni).reaction_source_.Get(b.smr).rate_);
  ASSERT_TRUE(l.denominators_[0].adsorption_.back().exponent_.count("Nitrogen"));
  err.clear();
  EXPECT_TRUE(RemoveComponentFromDeck(d, "Nitrogen", err)) << err;
  EXPECT_FALSE(l.denominators_[0].adsorption_.back().exponent_.count("Nitrogen"));

  // A nonzero exponent counts as use.
  build::AddComponent(d, "Nitrogen");
  EXPECT_EQ(build::AddAdsorptionTerm(d, b.ni, b.smr), 4u);
  b.sim.RefreshParameters();
  b.sim.Parameters().Get("adsorbent[1].reaction[1].rate.adsorption[4].exponent[Nitrogen]").SetValue(1.0);
  EXPECT_FALSE(RemoveComponentFromDeck(d, "Nitrogen", err));
  EXPECT_EQ(err, "Nitrogen is used by reaction SMR (on Ni catalyst); remove it first");

  // After removing the reactions the component goes.
  build::RemoveReaction(d, b.ni, b.smr);
  build::RemoveReaction(d, b.ni, b.wgs);
  EXPECT_TRUE(RemoveComponentFromDeck(d, H2O, err));
  RemoveComponentFromDeck(d, "Nitrogen");  // the unchecked form, nothing uses it now
  EXPECT_EQ(std::find(d.fluid.components_.begin(), d.fluid.components_.end(), "Nitrogen"),
            d.fluid.components_.end());
}

TEST(ReactionDeck, RemovingTheSolidRemovesItsReactions) {
  Reformer b;
  Build(b);
  UserData& d = b.sim.Data();
  // The catalyst is part of the layer's blend: refused until the blend no
  // longer holds it.
  std::string why;
  EXPECT_FALSE(RemoveAdsorbentFromDeck(d, b.ni, why));
  EXPECT_TRUE(d.adsorbent_directory.Contains(b.ni));
  build::SetLayerComposition(d, b.layer, {{b.sorbent, 1.0}});
  EXPECT_TRUE(RemoveAdsorbentFromDeck(d, b.ni, why)) << why;
  EXPECT_FALSE(d.adsorbent_directory.Contains(b.ni));
  std::string err;
  EXPECT_TRUE(RemoveComponentFromDeck(d, H2O, err)) << err;
}

TEST(ReactionDeck, ImportedSolidKeepsReactionsWithMissingComponents) {
  Reformer b;
  Build(b);
  std::ostringstream os;
  {
    cereal::JSONOutputArchive ar(os);
    ar(b.sim.Data().adsorbent_directory.Get(b.ni));
  }
  // A deck with methane and hydrogen only.
  Simulation other;
  UserData& d = other.Data();
  build::AddComponent(d, CH4);
  build::AddComponent(d, H2);
  AdsorbentID id = 0;
  std::string err;
  ASSERT_TRUE(AddAdsorbentFromJson(d, "Ni catalyst", os.str(), id, err)) << err;
  const auto& src = d.adsorbent_directory.Get(id).reaction_source_;
  ASSERT_EQ(src.reactions_.size(), 2u);
  EXPECT_FALSE(src.reactions_[0].IsActive(d.fluid.components_));

  // The deck check reports them once a layer uses the solid.
  err.clear();
  EXPECT_TRUE(ReactionChecks(d.adsorbent_directory, d.bed, d.fluid, err));  // no layer yet
  const size_t layer = build::AddLayer(d, "L");
  build::SetLayerComposition(d, layer, {{id, 1.0}});
  EXPECT_FALSE(ReactionChecks(d.adsorbent_directory, d.bed, d.fluid, err));
  EXPECT_NE(err.find("uses Water, which is not a component of the deck"), std::string::npos) << err;
}

TEST(ReactionDeck, DeckChecks) {
  Reformer b;
  Build(b);
  UserData& d = b.sim.Data();
  std::string err;
  ASSERT_TRUE(d.EvaluateAllExpressions(err)) << err;
  EXPECT_TRUE(ReactionChecks(d.adsorbent_directory, d.bed, d.fluid, err)) << err;

  auto& src = d.adsorbent_directory.Get(b.ni).reaction_source_;
  // Element imbalance.
  build::SetReactionStoichiometry(d, b.ni, b.smr, {{CH4, -1}, {H2O, -1}, {CO, 1}, {H2, 2}});
  err.clear();
  EXPECT_FALSE(ReactionChecks(d.adsorbent_directory, d.bed, d.fluid, err));
  EXPECT_NE(err.find("does not balance the elements: H -2"), std::string::npos) << err;
  build::SetReactionStoichiometry(d, b.ni, b.smr, {{CH4, -1}, {H2O, -1}, {CO, 1}, {H2, 3}});

  // No product.
  build::SetReactionStoichiometry(d, b.ni, b.wgs, {{CO, -1}});
  err.clear();
  EXPECT_FALSE(ReactionChecks(d.adsorbent_directory, d.bed, d.fluid, err));
  EXPECT_NE(err.find("has no product"), std::string::npos) << err;
  build::SetReactionStoichiometry(d, b.ni, b.wgs, {{CO, -1}, {H2O, -1}, {CO2, 1}, {H2, 1}});

  // Negative pre-exponential factor.
  dynamic_cast<ReversiblePowerLawRate&>(*src.Get(b.wgs).rate_).pre_exp_factor_ = -1.0;
  err.clear();
  EXPECT_FALSE(ReactionChecks(d.adsorbent_directory, d.bed, d.fluid, err));
  EXPECT_NE(err.find("pre-exponential factor pre_exp_factor must be finite and not negative"), std::string::npos)
      << err;
  dynamic_cast<ReversiblePowerLawRate&>(*src.Get(b.wgs).rate_).pre_exp_factor_ = 1.0;

  // UserData::PerformChecks runs the reaction checks.
  build::SetReactionStoichiometry(d, b.ni, b.wgs, {{CO, -1}});
  err.clear();
  EXPECT_FALSE(d.PerformChecks(err));
  EXPECT_NE(err.find("has no product"), std::string::npos) << err;
  build::SetReactionStoichiometry(d, b.ni, b.wgs, {{CO, -1}, {H2O, -1}, {CO2, 1}, {H2, 1}});

  // Thiele without a key reactant (the key reactant left the stoichiometry).
  build::SetEffectiveness(d, b.ni, b.wgs, "thiele", CO);
  src.Get(b.wgs).effectiveness_.key_component_ = CO2;
  err.clear();
  EXPECT_FALSE(ReactionChecks(d.adsorbent_directory, d.bed, d.fluid, err));
  EXPECT_NE(err.find("needs a key reactant"), std::string::npos) << err;
}

TEST(ReactionDeck, Warnings) {
  Reformer b;
  Build(b);
  UserData& d = b.sim.Data();
  EXPECT_TRUE(ReactionWarnings(d.adsorbent_directory, d.bed).empty());

  // Nonzero heat of reaction in an isothermal run.
  d.adsorbent_directory.Get(b.ni).reaction_source_.Get(b.smr).heat_of_reaction_ = 206100.0;
  d.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
  auto w = ReactionWarnings(d.adsorbent_directory, d.bed);
  ASSERT_EQ(w.size(), 1u);
  EXPECT_NE(w[0].find("isothermal"), std::string::npos);

  // Non-elementary orders: informational, only on request.
  b.sim.Parameters().Get("adsorbent[1].reaction[1].rate.order[Hydrogen]").SetValue(-2.5);
  EXPECT_EQ(ReactionWarnings(d.adsorbent_directory, d.bed).size(), 1u);
  w = ReactionWarnings(d.adsorbent_directory, d.bed, true);
  ASSERT_EQ(w.size(), 2u);
  EXPECT_NE(w[1].find("orders of Hydrogen"), std::string::npos) << w[1];

  // A catalyst no layer uses.
  build::SetLayerComposition(d, b.layer, {{b.sorbent, 1.0}});
  w = ReactionWarnings(d.adsorbent_directory, d.bed);
  ASSERT_EQ(w.size(), 1u);
  EXPECT_NE(w[0].find("no layer uses it"), std::string::npos) << w[0];
}

// ---- edit cascades (deck_edit.h) --------------------------------------------------------

const LhhwRate& LhhwOf(const Reaction& r) { return dynamic_cast<const LhhwRate&>(*r.rate_); }

const std::map<std::string, ExpressionParam>& OrdersOf(const Reaction& r) {
  if (auto* p = dynamic_cast<const PowerLawRate*>(r.rate_.get())) return p->order_;
  if (auto* p = dynamic_cast<const ReversiblePowerLawRate*>(r.rate_.get())) return p->order_;
  return LhhwOf(r).order_;
}

TEST(ReactionEdit, CascadesEditTheStructure) {
  Reformer b;
  Build(b);
  UserData& d = b.sim.Data();
  const auto& src = d.adsorbent_directory.Get(b.ni).reaction_source_;
  std::string err;

  // New reactions: a default name, a taken name made unique, an empty stoichiometry.
  ReactionID r = 0, r2 = 0;
  ASSERT_TRUE(AddReaction(d, b.ni, "", ReactionRateType::PowerLaw, r, err)) << err;
  EXPECT_EQ(r, 3u);
  EXPECT_EQ(src.Get(r).name_, "Reaction 3");
  EXPECT_TRUE(src.Get(r).nu_.empty());
  EXPECT_EQ(src.Get(r).rate_->GetType(), ReactionRateType::PowerLaw);
  ASSERT_TRUE(AddReaction(d, b.ni, "SMR", ReactionRateType::Expression, r2, err)) << err;
  EXPECT_EQ(src.Get(r2).name_, "SMR 2");
  EXPECT_EQ(src.Get(r2).rate_->GetType(), ReactionRateType::Expression);
  // The same name on another solid is not taken.
  ReactionID other = 0;
  ASSERT_TRUE(AddReaction(d, b.sorbent, "SMR", ReactionRateType::PowerLaw, other, err)) << err;
  EXPECT_EQ(d.adsorbent_directory.Get(b.sorbent).reaction_source_.Get(other).name_, "SMR");

  ASSERT_TRUE(RenameReaction(d, b.ni, r, "Methanation", err)) << err;
  EXPECT_EQ(src.Get(r).name_, "Methanation");
  EXPECT_TRUE(RenameReaction(d, b.ni, r, "Methanation", err));  // its own name

  // Stoichiometry, one entry at a time; new entries take the default orders.
  for (const auto& [c, v] : std::map<std::string, double>{{CO, -1}, {H2, -3}, {CH4, 1}, {H2O, 1}})
    ASSERT_TRUE(SetStoichiometricCoefficient(d, b.ni, r, c, v, err)) << err;
  EXPECT_EQ(src.Get(r).nu_, (std::map<std::string, double>{{CO, -1}, {H2, -3}, {CH4, 1}, {H2O, 1}}));
  EXPECT_DOUBLE_EQ(OrdersOf(src.Get(r)).at(CO).Get(), 1.0);
  EXPECT_DOUBLE_EQ(OrdersOf(src.Get(r)).at(H2).Get(), 3.0);
  EXPECT_DOUBLE_EQ(OrdersOf(src.Get(r)).at(CH4).Get(), 0.0);
  // A changed coefficient keeps the order; 0 removes the entry and its order.
  ASSERT_TRUE(SetStoichiometricCoefficient(d, b.ni, r, CO, -2, err)) << err;
  EXPECT_DOUBLE_EQ(src.Get(r).nu_.at(CO), -2.0);
  EXPECT_DOUBLE_EQ(OrdersOf(src.Get(r)).at(CO).Get(), 1.0);
  ASSERT_TRUE(SetStoichiometricCoefficient(d, b.ni, r, H2O, 0, err)) << err;
  EXPECT_FALSE(src.Get(r).nu_.count(H2O));
  EXPECT_FALSE(OrdersOf(src.Get(r)).count(H2O));
  EXPECT_TRUE(SetStoichiometricCoefficient(d, b.ni, r, H2O, 0, err));  // absent: nothing to do
  // A re-added entry takes the default order again, not the one it had.
  b.sim.RefreshParameters();
  b.sim.Parameters().Get("adsorbent[1].reaction[3].rate.order[Hydrogen]").SetValue(0.5);
  ASSERT_TRUE(SetStoichiometricCoefficient(d, b.ni, r, H2, 0, err)) << err;
  ASSERT_TRUE(SetStoichiometricCoefficient(d, b.ni, r, H2, -3, err)) << err;
  EXPECT_DOUBLE_EQ(OrdersOf(src.Get(r)).at(H2).Get(), 3.0);

  // Kind switch: the orders carry over; a new lhhw rate has one empty denominator.
  ASSERT_TRUE(SetReactionKind(d, b.ni, r, ReactionRateType::Lhhw, err)) << err;
  ASSERT_EQ(src.Get(r).rate_->GetType(), ReactionRateType::Lhhw);
  EXPECT_DOUBLE_EQ(OrdersOf(src.Get(r)).at(CO).Get(), 1.0);
  ASSERT_EQ(LhhwOf(src.Get(r)).denominators_.size(), 1u);
  EXPECT_TRUE(LhhwOf(src.Get(r)).denominators_[0].adsorption_.empty());

  ASSERT_TRUE(SetReactionBasis(d, b.ni, r, RateBasis::Concentration, err)) << err;
  EXPECT_EQ(src.Get(r).basis_, RateBasis::Concentration);

  // Effectiveness: a key reactant, an empty key, back to constant.
  ASSERT_TRUE(SetReactionEffectiveness(d, b.ni, r, EffectivenessMode::Thiele, CO, err)) << err;
  EXPECT_EQ(src.Get(r).effectiveness_.mode_, EffectivenessMode::Thiele);
  EXPECT_EQ(src.Get(r).effectiveness_.key_component_, CO);
  // The key turning into a product clears it (the first reactant is used).
  ASSERT_TRUE(SetStoichiometricCoefficient(d, b.ni, r, CO, 1, err)) << err;
  EXPECT_EQ(src.Get(r).effectiveness_.key_component_, "");
  EXPECT_EQ(src.Get(r).KeyReactant(), H2);
  ASSERT_TRUE(SetStoichiometricCoefficient(d, b.ni, r, CO, -1, err)) << err;
  ASSERT_TRUE(SetReactionEffectiveness(d, b.ni, r, EffectivenessMode::Thiele, H2, err)) << err;
  ASSERT_TRUE(SetStoichiometricCoefficient(d, b.ni, r, H2, 0, err)) << err;  // removed: cleared too
  EXPECT_EQ(src.Get(r).effectiveness_.key_component_, "");
  ASSERT_TRUE(SetReactionEffectiveness(d, b.ni, r, EffectivenessMode::Thiele, "", err)) << err;
  EXPECT_EQ(src.Get(r).effectiveness_.mode_, EffectivenessMode::Thiele);
  ASSERT_TRUE(SetReactionEffectiveness(d, b.ni, r, EffectivenessMode::Constant, CO, err)) << err;
  EXPECT_EQ(src.Get(r).effectiveness_.mode_, EffectivenessMode::Constant);
  EXPECT_EQ(src.Get(r).effectiveness_.key_component_, "");

  // LHHW structure.
  size_t k = 0, j = 0;
  ASSERT_TRUE(AddReactionDenominator(d, b.ni, r, k, err)) << err;
  EXPECT_EQ(k, 2u);
  ASSERT_TRUE(AddReactionAdsorptionTerm(d, b.ni, r, 2, j, err)) << err;
  EXPECT_EQ(j, 1u);
  ASSERT_TRUE(AddReactionAdsorptionTerm(d, b.ni, r, 2, j, err)) << err;
  EXPECT_EQ(j, 2u);
  ASSERT_TRUE(AddReactionAdsorptionTerm(d, b.ni, r, 1, j, err)) << err;
  EXPECT_EQ(j, 1u);
  EXPECT_EQ(LhhwOf(src.Get(r)).denominators_[1].adsorption_[0].exponent_.size(), 5u);  // every component
  ASSERT_TRUE(RemoveReactionAdsorptionTerm(d, b.ni, r, 2, 1, err)) << err;
  EXPECT_EQ(LhhwOf(src.Get(r)).denominators_[1].adsorption_.size(), 1u);
  ASSERT_TRUE(RemoveReactionDenominator(d, b.ni, r, 1, err)) << err;  // the second moves down
  ASSERT_EQ(LhhwOf(src.Get(r)).denominators_.size(), 1u);
  EXPECT_EQ(LhhwOf(src.Get(r)).denominators_[0].adsorption_.size(), 1u);

  // Removal: the others keep their IDs, and IDs are not reused.
  ASSERT_TRUE(RemoveReaction(d, b.ni, r, err)) << err;
  EXPECT_EQ(src.Find(r), nullptr);
  EXPECT_NE(src.Find(b.smr), nullptr);
  EXPECT_NE(src.Find(r2), nullptr);
  ReactionID r3 = 0;
  ASSERT_TRUE(AddReaction(d, b.ni, "Later", ReactionRateType::PowerLaw, r3, err)) << err;
  EXPECT_EQ(r3, 5u);

  // The edited deck reaches the registry.
  b.sim.RefreshParameters();
  EXPECT_NO_THROW(b.sim.Parameters().Get("adsorbent[Ni catalyst].reaction[Later].rate.pre_exp_factor"));
}

TEST(ReactionEdit, RefusalsLeaveTheDeckUnchanged) {
  Reformer b;
  Build(b);
  UserData& d = b.sim.Data();
  build::SetEffectiveness(d, b.ni, b.smr, "thiele", CH4);
  const std::string before = Snapshot(d);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  ReactionID id = 0;
  size_t idx = 0;

  struct Refusal {
    std::string what;
    std::function<bool(std::string&)> call;
    std::string message;  // a part of the expected message
  };
  const std::vector<Refusal> refusals = {
      {"add: unknown solid", [&](std::string& e) { return AddReaction(d, 99, "X", ReactionRateType::PowerLaw, id, e); },
       "no adsorbent with ID 99"},
      {"add: unknown kind",
       [&](std::string& e) { return AddReaction(d, b.ni, "X", static_cast<ReactionRateType>(7), id, e); },
       "unknown reaction rate kind 7"},
      {"remove: unknown reaction", [&](std::string& e) { return RemoveReaction(d, b.ni, 99, e); },
       "no reaction with ID 99 on adsorbent 'Ni catalyst'"},
      {"remove: unknown solid", [&](std::string& e) { return RemoveReaction(d, 99, b.smr, e); },
       "no adsorbent with ID 99"},
      {"remove: another solid's reaction", [&](std::string& e) { return RemoveReaction(d, b.sorbent, b.smr, e); },
       "no reaction with ID 1 on adsorbent 'Sorbent'"},
      {"rename: empty", [&](std::string& e) { return RenameReaction(d, b.ni, b.smr, "", e); }, "needs a name"},
      {"rename: taken", [&](std::string& e) { return RenameReaction(d, b.ni, b.smr, "WGS", e); },
       "already has a reaction named 'WGS'"},
      {"kind: unknown kind",
       [&](std::string& e) { return SetReactionKind(d, b.ni, b.smr, static_cast<ReactionRateType>(-1), e); },
       "unknown reaction rate kind"},
      {"kind: unknown reaction",
       [&](std::string& e) { return SetReactionKind(d, b.ni, 99, ReactionRateType::PowerLaw, e); }, "no reaction"},
      {"coefficient: component not in the deck",
       [&](std::string& e) { return SetStoichiometricCoefficient(d, b.ni, b.smr, "Argon", -1, e); },
       "'Argon' is not a component of this deck"},
      {"coefficient: unknown name",
       [&](std::string& e) { return SetStoichiometricCoefficient(d, b.ni, b.smr, "Unobtainium", -1, e); },
       "not a component"},
      {"coefficient: NaN", [&](std::string& e) { return SetStoichiometricCoefficient(d, b.ni, b.smr, CO, nan, e); },
       "not finite"},
      {"coefficient: infinite",
       [&](std::string& e) { return SetStoichiometricCoefficient(d, b.ni, b.smr, CO, -inf, e); }, "not finite"},
      {"coefficient: unknown reaction",
       [&](std::string& e) { return SetStoichiometricCoefficient(d, b.ni, 99, CO, -1, e); }, "no reaction"},
      {"basis: unknown", [&](std::string& e) { return SetReactionBasis(d, b.ni, b.smr, static_cast<RateBasis>(2), e); },
       "unknown reaction basis kind 2"},
      {"effectiveness: unknown mode",
       [&](std::string& e) {
         return SetReactionEffectiveness(d, b.ni, b.smr, static_cast<EffectivenessMode>(5), "", e);
       },
       "unknown effectiveness kind 5"},
      {"effectiveness: key is a product",
       [&](std::string& e) { return SetReactionEffectiveness(d, b.ni, b.smr, EffectivenessMode::Thiele, H2, e); },
       "'Hydrogen' is not one"},
      {"effectiveness: key not in the reaction",
       [&](std::string& e) { return SetReactionEffectiveness(d, b.ni, b.smr, EffectivenessMode::Thiele, CO2, e); },
       "key reactant"},
      {"denominator: not lhhw", [&](std::string& e) { return AddReactionDenominator(d, b.ni, b.wgs, idx, e); },
       "does not use the lhhw rate kind"},
      {"remove denominator: the last one",
       [&](std::string& e) { return RemoveReactionDenominator(d, b.ni, b.smr, 1, e); }, "at least one denominator"},
      {"remove denominator: absent", [&](std::string& e) { return RemoveReactionDenominator(d, b.ni, b.smr, 2, e); },
       "has no denominator 2"},
      {"remove denominator: 0", [&](std::string& e) { return RemoveReactionDenominator(d, b.ni, b.smr, 0, e); },
       "has no denominator 0"},
      {"remove denominator: not lhhw",
       [&](std::string& e) { return RemoveReactionDenominator(d, b.ni, b.wgs, 1, e); }, "lhhw"},
      {"term: absent denominator",
       [&](std::string& e) { return AddReactionAdsorptionTerm(d, b.ni, b.smr, 2, idx, e); }, "has no denominator 2"},
      {"term: not lhhw", [&](std::string& e) { return AddReactionAdsorptionTerm(d, b.ni, b.wgs, 1, idx, e); }, "lhhw"},
      {"remove term: absent", [&](std::string& e) { return RemoveReactionAdsorptionTerm(d, b.ni, b.smr, 1, 3, e); },
       "has no adsorption term 3"},
      {"remove term: 0", [&](std::string& e) { return RemoveReactionAdsorptionTerm(d, b.ni, b.smr, 1, 0, e); },
       "has no adsorption term 0"},
      {"remove term: absent denominator",
       [&](std::string& e) { return RemoveReactionAdsorptionTerm(d, b.ni, b.smr, 2, 1, e); }, "has no denominator 2"},
      {"remove term: not lhhw",
       [&](std::string& e) { return RemoveReactionAdsorptionTerm(d, b.ni, b.wgs, 1, 1, e); }, "lhhw"},
  };
  for (const auto& r : refusals) {
    std::string err;
    EXPECT_FALSE(r.call(err)) << r.what;
    EXPECT_NE(err.find(r.message), std::string::npos) << r.what << ": " << err;
    EXPECT_EQ(Snapshot(d), before) << r.what;
  }

  // Edits that change nothing succeed and leave the deck as it was too.
  std::string err;
  EXPECT_TRUE(RenameReaction(d, b.ni, b.smr, "SMR", err));
  EXPECT_TRUE(SetReactionKind(d, b.ni, b.smr, ReactionRateType::Lhhw, err));
  EXPECT_TRUE(SetStoichiometricCoefficient(d, b.ni, b.smr, CH4, -1, err));
  EXPECT_TRUE(SetStoichiometricCoefficient(d, b.ni, b.smr, CO2, 0, err));
  EXPECT_EQ(Snapshot(d), before);
}

TEST(ReactionEdit, ReactionKpisNameComponentsAndOutliveReactionEdits) {
  Reformer b;
  Build(b);
  UserData& d = b.sim.Data();
  build::SetReactionKPI(d, CH4, {H2, CO2});
  std::string err;
  // Removing a stoichiometry entry or a whole reaction leaves the components
  // in the deck, so the reaction KPIs stay as they are.
  ASSERT_TRUE(SetStoichiometricCoefficient(d, b.ni, b.wgs, CO2, 0, err)) << err;
  ASSERT_TRUE(RemoveReaction(d, b.ni, b.wgs, err)) << err;
  ASSERT_TRUE(RemoveReaction(d, b.ni, b.smr, err)) << err;
  EXPECT_EQ(d.kpi_config.reaction.key_reactant, CH4);
  EXPECT_EQ(d.kpi_config.reaction.products, (std::vector<std::string>{H2, CO2}));
  // With no reaction using them, the components can go, and the KPIs forget them.
  ASSERT_TRUE(RemoveComponentFromDeck(d, CO2, err)) << err;
  ASSERT_TRUE(RemoveComponentFromDeck(d, CH4, err)) << err;
  EXPECT_FALSE(d.kpi_config.reaction.IsSet());
  EXPECT_EQ(d.kpi_config.reaction.products, (std::vector<std::string>{H2}));
}

// The two-site LHHW rate of Halabi (2011, ch. 8): a metal-site denominator
// and a support-site denominator, each with its adsorption terms, next to a
// reversible water-gas shift. Built once with the builder and once only
// through the edit cascades, the decks are the same.
void TwoSiteBase(Simulation& sim, AdsorbentID& cat) {
  UserData& d = sim.Data();
  for (const auto& c : {CH4, H2O, CO, H2, CO2}) build::AddComponent(d, c);
  cat = build::AddAdsorbent(d, "Rh catalyst");
  const size_t layer = build::AddLayer(d, "Bed");
  build::SetLayerComposition(d, layer, {{cat, 1.0}});
}

TEST(ReactionEdit, CascadesBuildTheBuilderDeck) {
  const std::map<std::string, double> r1 = {{CH4, -1}, {H2O, -1}, {CO, 1}, {H2, 3}};
  const std::map<std::string, double> r2 = {{CO, -1}, {H2O, -1}, {CO2, 1}, {H2, 1}};

  Simulation built;
  AdsorbentID cat = 0;
  TwoSiteBase(built, cat);
  {
    UserData& d = built.Data();
    const ReactionID a = build::AddReaction(d, cat, "R1", r1, "lhhw");
    build::AddAdsorptionTerm(d, cat, a, 1);  // metal sites
    build::AddAdsorptionTerm(d, cat, a, 1);
    EXPECT_EQ(build::AddDenominator(d, cat, a), 2u);
    build::AddAdsorptionTerm(d, cat, a, 2);  // support sites
    build::AddAdsorptionTerm(d, cat, a, 2);
    build::AddAdsorptionTerm(d, cat, a, 2);
    build::SetReactionBasis(d, cat, a, "partial_pressure");
    build::SetEffectiveness(d, cat, a, "thiele", CH4);
    const ReactionID w = build::AddReaction(d, cat, "R2", r2, "reversible_power_law");
    build::SetReactionBasis(d, cat, w, "concentration");
  }

  Simulation edited;
  AdsorbentID cat2 = 0;
  TwoSiteBase(edited, cat2);
  {
    UserData& d = edited.Data();
    std::string err;
    ReactionID a = 0, w = 0;
    size_t k = 0, j = 0;
    ASSERT_TRUE(AddReaction(d, cat2, "R1", ReactionRateType::Lhhw, a, err)) << err;
    for (const auto& [c, v] : r1) ASSERT_TRUE(SetStoichiometricCoefficient(d, cat2, a, c, v, err)) << err;
    ASSERT_TRUE(AddReactionAdsorptionTerm(d, cat2, a, 1, j, err)) << err;
    ASSERT_TRUE(AddReactionAdsorptionTerm(d, cat2, a, 1, j, err)) << err;
    ASSERT_TRUE(AddReactionDenominator(d, cat2, a, k, err)) << err;
    for (int n = 0; n < 3; ++n) ASSERT_TRUE(AddReactionAdsorptionTerm(d, cat2, a, k, j, err)) << err;
    EXPECT_EQ(j, 3u);
    ASSERT_TRUE(SetReactionBasis(d, cat2, a, RateBasis::PartialPressure, err)) << err;
    ASSERT_TRUE(SetReactionEffectiveness(d, cat2, a, EffectivenessMode::Thiele, CH4, err)) << err;
    // The second reaction starts under another name and kind and is set right.
    ASSERT_TRUE(AddReaction(d, cat2, "", ReactionRateType::PowerLaw, w, err)) << err;
    ASSERT_TRUE(SetReactionKind(d, cat2, w, ReactionRateType::ReversiblePowerLaw, err)) << err;
    ASSERT_TRUE(RenameReaction(d, cat2, w, "R2", err)) << err;
    for (const auto& [c, v] : r2) ASSERT_TRUE(SetStoichiometricCoefficient(d, cat2, w, c, v, err)) << err;
    ASSERT_TRUE(SetReactionBasis(d, cat2, w, RateBasis::Concentration, err)) << err;
  }
  EXPECT_EQ(Snapshot(edited.Data()), Snapshot(built.Data()));
}

}  // namespace
}  // namespace adsim
