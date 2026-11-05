/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// Reactions in a deck (adsim/objects/reaction.h): the builder, serialisation
// (round trip, old v4 adsorbents and decks still load), the parameter
// registry, the component-removal cascade, the deck checks and warnings.

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <functional>
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
  EXPECT_EQ(dynamic_cast<LhhwRate&>(*src.Get(b.smr).rate_).adsorption_.size(), 2u);

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
  EXPECT_EQ(dynamic_cast<LhhwRate&>(*src.Get(b.smr).rate_).adsorption_.size(), 1u);
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
                           "rate.equilibrium_energy", "rate.denominator_exponent", "rate.reversible",
                           "rate.adsorption[1].pre_exp_factor", "rate.adsorption[1].energy",
                           "rate.adsorption[2].exponent[Hydrogen]", "heat_of_reaction", "heat_mode",
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
  reg.Get(head + "rate.adsorption[2].exponent[Hydrogen]").SetValue(-1.0);
  const auto& smr = b.sim.Data().adsorbent_directory.Get(b.ni).reaction_source_.Get(b.smr);
  const auto& l = dynamic_cast<const LhhwRate&>(*smr.rate_);
  EXPECT_DOUBLE_EQ(l.order_.at(H2).Get(), -2.5);
  EXPECT_DOUBLE_EQ(l.adsorption_[1].exponent_.at(H2).Get(), -1.0);

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
    EXPECT_EQ(l.adsorption_.size(), 2u);
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
  ASSERT_TRUE(l.adsorption_.back().exponent_.count("Nitrogen"));
  err.clear();
  EXPECT_TRUE(RemoveComponentFromDeck(d, "Nitrogen", err)) << err;
  EXPECT_FALSE(l.adsorption_.back().exponent_.count("Nitrogen"));

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
  RemoveAdsorbentFromDeck(d, b.ni);
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

}  // namespace
}  // namespace adsim
