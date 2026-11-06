/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// Reaction rate kinds (adsim/equations/reaction_kinds.h) against hand-computed
// values: the Xu and Froment (1989) steam-reforming rates at two states,
// elementary power laws, the activity clipping and floor, the concentration
// basis and the expression kind; plus the heat of reaction.

#include <cmath>
#include <map>
#include <string>

#include <gtest/gtest.h>

#include "adsim/equations/reaction_kinds.h"
#include "adsim/math/constants.h"
#include "adsim/objects/reaction.h"

namespace adsim {
namespace {

constexpr double R = constants::kIdealGasConstant;

const std::string CH4 = "Methane", H2O = "Water", CO = "CarbonMonoxide", H2 = "Hydrogen", CO2 = "CarbonDioxide";

// The Xu-Froment denominator: 1 + K_CO p_CO + K_H2 p_H2 + K_CH4 p_CH4 + K_H2O p_H2O / p_H2.
void AddXuFromentDenominator(LhhwRate& r) {
  auto term = [&](double A, double E, std::map<std::string, double> e) {
    LhhwAdsorptionTerm t;
    t.pre_exp_factor_ = A;
    t.energy_ = E;
    for (const auto& [c, v] : e) t.exponent_[c] = v;
    r.adsorption_.push_back(std::move(t));
  };
  term(8.23e-5, -70650.0, {{CO, 1.0}});
  term(6.12e-9, -82900.0, {{H2, 1.0}});
  term(6.65e-4, -38280.0, {{CH4, 1.0}});
  term(1.77e5, 88680.0, {{H2O, 1.0}, {H2, -1.0}});
  r.denominator_exponent_ = 2.0;
}

// r1: CH4 + H2O <-> CO + 3 H2, k1 = 4.225e15 kmol bar^0.5/(kg h), E1 = 240.1 kJ/mol,
// K1 = 1.198e17 exp(-26830/T) bar^2. Rates in mol/(kg s).
Reaction Smr() {
  Reaction r("SMR", 1);
  r.SetStoichiometry({{CH4, -1}, {H2O, -1}, {CO, 1}, {H2, 3}});
  r.SetRateType(ReactionRateType::Lhhw);
  auto& l = dynamic_cast<LhhwRate&>(*r.rate_);
  l.pre_exp_factor_ = 4.225e15 * 1000.0 / 3600.0;
  l.activation_energy_ = 240100.0;
  l.order_[H2] = -2.5;
  l.equilibrium_pre_exp_factor_ = 1.198e17;
  l.equilibrium_energy_ = 26830.0 * R;
  AddXuFromentDenominator(l);
  return r;
}

// r2: CO + H2O <-> CO2 + H2, k2 = 1.955e6 kmol/(kg h bar), E2 = 67.13 kJ/mol,
// K2 = 1.767e-2 exp(4400/T).
Reaction Wgs() {
  Reaction r("WGS", 2);
  r.SetStoichiometry({{CO, -1}, {H2O, -1}, {CO2, 1}, {H2, 1}});
  r.SetRateType(ReactionRateType::Lhhw);
  auto& l = dynamic_cast<LhhwRate&>(*r.rate_);
  l.pre_exp_factor_ = 1.955e6 * 1000.0 / 3600.0;
  l.activation_energy_ = 67130.0;
  l.order_[H2] = -1.0;
  l.equilibrium_pre_exp_factor_ = 1.767e-2;
  l.equilibrium_energy_ = -4400.0 * R;
  AddXuFromentDenominator(l);
  return r;
}

// r3: CH4 + 2 H2O <-> CO2 + 4 H2, k3 = 1.020e15 kmol bar^0.5/(kg h),
// E3 = 243.9 kJ/mol, K3 = K1 K2.
Reaction Overall() {
  Reaction r("Overall", 3);
  r.SetStoichiometry({{CH4, -1}, {H2O, -2}, {CO2, 1}, {H2, 4}});
  r.SetRateType(ReactionRateType::Lhhw);
  auto& l = dynamic_cast<LhhwRate&>(*r.rate_);
  l.pre_exp_factor_ = 1.020e15 * 1000.0 / 3600.0;
  l.activation_energy_ = 243900.0;
  l.order_[H2] = -3.5;
  l.equilibrium_pre_exp_factor_ = 1.198e17 * 1.767e-2;
  l.equilibrium_energy_ = (26830.0 - 4400.0) * R;
  AddXuFromentDenominator(l);
  return r;
}

struct State {
  double T, P;  // K, Pa
  std::map<std::string, double> y;
};

// Reference values computed independently from the Xu-Froment equations as
// printed (pressures in bar, constants as above).
struct Reference {
  State s;
  double den, r1, r2, r3;
};

const Reference kReferences[] = {
    {{773.15, 10e5, {{CH4, 0.2}, {H2O, 0.6}, {H2, 0.1}, {CO, 0.05}, {CO2, 0.05}}},
     5.0394611743493201, 0.033321586424642113, 1.8108131003226489, 0.026736007255304492},
    {{923.15, 5e5, {{CH4, 0.1}, {H2O, 0.3}, {H2, 0.4}, {CO, 0.1}, {CO2, 0.1}}},
     2.7330000080420236, 0.54116020124759323, 1.5516200295459461, 0.059728026846282881},
};

ReactionState At(const State& s) {
  ReactionState st;
  st.T = s.T;
  st.P = s.P;
  st.y = &s.y;
  return st;
}

TEST(ReactionKinds, KindTagsRoundTrip) {
  for (auto t : kAllReactionRateKinds) {
    ReactionRateType back{};
    ASSERT_TRUE(ReactionRateKindFromString(ReactionRateKindString(t), back));
    EXPECT_EQ(back, t);
  }
  EXPECT_STREQ(ReactionRateKindString(ReactionRateType::Lhhw), "lhhw");
  ReactionRateType out{};
  EXPECT_FALSE(ReactionRateKindFromString("langmuir", out));
}

TEST(ReactionKinds, XuFromentSteamReformingAtTwoStates) {
  Reaction r1 = Smr(), r2 = Wgs(), r3 = Overall();
  for (const auto& ref : kReferences) {
    const ReactionState st = At(ref.s);
    const double den = dynamic_cast<LhhwRate&>(*r1.rate_).Denominator(st, RateBasis::PartialPressure);
    EXPECT_NEAR(den, ref.den, 1e-12 * ref.den);
    EXPECT_NEAR(r1.IntrinsicRate(st), ref.r1, 1e-10 * std::abs(ref.r1));
    EXPECT_NEAR(r2.IntrinsicRate(st), ref.r2, 1e-10 * std::abs(ref.r2));
    EXPECT_NEAR(r3.IntrinsicRate(st), ref.r3, 1e-10 * std::abs(ref.r3));
  }
}

TEST(ReactionKinds, LhhwMatchesTheExplicitFormula) {
  // The same r1, written out term by term in this test.
  Reaction r1 = Smr();
  const State& s = kReferences[1].s;
  const double T = s.T;
  auto p = [&](const std::string& c) { return s.y.at(c) * s.P / 1e5; };
  auto K = [&](double A, double E) { return A * std::exp(-E / (R * T)); };
  const double den = 1 + K(8.23e-5, -70650) * p(CO) + K(6.12e-9, -82900) * p(H2) + K(6.65e-4, -38280) * p(CH4) +
                     K(1.77e5, 88680) * p(H2O) / p(H2);
  const double k1 = K(4.225e15 / 3.6, 240100);
  const double K1 = K(1.198e17, 26830 * R);
  const double expected =
      k1 * (std::pow(p(H2), -2.5) * p(CH4) * p(H2O) - std::pow(p(H2), 0.5) * p(CO) / K1) / (den * den);
  EXPECT_NEAR(r1.IntrinsicRate(At(s)), expected, 1e-12 * std::abs(expected));
}

TEST(ReactionKinds, IrreversibleLhhwDropsTheEquilibriumTerm) {
  Reaction r1 = Smr();
  auto& l = dynamic_cast<LhhwRate&>(*r1.rate_);
  l.reversible_ = false;
  const State& s = kReferences[0].s;
  auto p = [&](const std::string& c) { return s.y.at(c) * s.P / 1e5; };
  const double k1 = 4.225e15 / 3.6 * std::exp(-240100 / (R * s.T));
  const double den = kReferences[0].den;
  const double expected = k1 * std::pow(p(H2), -2.5) * p(CH4) * p(H2O) / (den * den);
  EXPECT_NEAR(r1.IntrinsicRate(At(s)), expected, 1e-10 * expected);
}

TEST(ReactionKinds, PowerLawDefaultsToElementaryOrders) {
  Reaction r("A", 1);
  r.SetStoichiometry({{CH4, -1}, {H2O, -2}, {CO2, 1}, {H2, 4}});
  r.SetRateType(ReactionRateType::PowerLaw);
  auto& pl = dynamic_cast<PowerLawRate&>(*r.rate_);
  EXPECT_DOUBLE_EQ(pl.order_.at(CH4).Get(), 1.0);
  EXPECT_DOUBLE_EQ(pl.order_.at(H2O).Get(), 2.0);
  EXPECT_DOUBLE_EQ(pl.order_.at(CO2).Get(), 0.0);
  EXPECT_DOUBLE_EQ(pl.order_.at(H2).Get(), 0.0);

  pl.pre_exp_factor_ = 2.0e3;
  pl.activation_energy_ = 5.0e4;
  const State s{600.0, 3e5, {{CH4, 0.25}, {H2O, 0.5}, {CO2, 0.0}, {H2, 0.25}}};
  const double expected = 2.0e3 * std::exp(-5.0e4 / (R * 600.0)) * 0.75 * 1.5 * 1.5;  // p in bar
  EXPECT_NEAR(r.IntrinsicRate(At(s)), expected, 1e-12 * expected);

  // Concentration basis: a_i = y_i P / (R T) [mol/m3].
  r.basis_ = RateBasis::Concentration;
  const double c = s.P / (R * s.T);
  const double expected_c = 2.0e3 * std::exp(-5.0e4 / (R * 600.0)) * (0.25 * c) * std::pow(0.5 * c, 2);
  EXPECT_NEAR(r.IntrinsicRate(At(s)), expected_c, 1e-12 * expected_c);
}

TEST(ReactionKinds, ReversiblePowerLawUsesTwoMonomials) {
  Reaction r("WGS", 1);
  r.SetStoichiometry({{CO, -1}, {H2O, -1}, {CO2, 1}, {H2, 1}});
  r.SetRateType(ReactionRateType::ReversiblePowerLaw);
  auto& rp = dynamic_cast<ReversiblePowerLawRate&>(*r.rate_);
  rp.pre_exp_factor_ = 10.0;
  rp.activation_energy_ = 0.0;
  rp.equilibrium_pre_exp_factor_ = 4.0;
  rp.equilibrium_energy_ = 0.0;
  // p_CO = 1, p_H2O = 2, p_CO2 = 0.5, p_H2 = 1.5 bar: r = 10 (1 * 2 - 0.5 * 1.5 / 4).
  const State s{500.0, 5e5, {{CO, 0.2}, {H2O, 0.4}, {CO2, 0.1}, {H2, 0.3}}};
  EXPECT_NEAR(r.IntrinsicRate(At(s)), 10.0 * (2.0 - 0.75 / 4.0), 1e-12);
  // No product: the reverse monomial vanishes, the rate stays finite.
  const State fresh{500.0, 5e5, {{CO, 0.5}, {H2O, 0.5}}};
  EXPECT_NEAR(r.IntrinsicRate(At(fresh)), 10.0 * 2.5 * 2.5, 1e-12);
}

TEST(ReactionKinds, ActivitiesAreClippedAndNegativePowersFloored) {
  Reaction r1 = Smr();
  // No hydrogen at all (a reformer inlet): p_H2 is floored at 1e-10 bar in
  // p_H2^-2.5 and in the K_H2O p_H2O / p_H2 term, so the rate is finite.
  const State inlet{773.15, 10e5, {{CH4, 0.25}, {H2O, 0.75}}};
  const double rate = r1.IntrinsicRate(At(inlet));
  EXPECT_TRUE(std::isfinite(rate));
  EXPECT_GT(rate, 0.0);
  // A negative mole fraction (a solver overshoot) is clipped at zero.
  Reaction pl("A", 1);
  pl.SetStoichiometry({{CH4, -1}, {H2, 2}});
  pl.SetRateType(ReactionRateType::PowerLaw);
  const State overshoot{600.0, 1e5, {{CH4, -0.01}, {H2, 1.01}}};
  EXPECT_EQ(pl.IntrinsicRate(At(overshoot)), 0.0);
}

TEST(ReactionKinds, ExpressionKindBindsComponentNames) {
  Reaction r("E", 1);
  r.SetStoichiometry({{CH4, -1}, {H2, 2}});
  r.SetRateType(ReactionRateType::Expression);
  auto& e = dynamic_cast<ExpressionRate&>(*r.rate_);
  e.expression_.SetExpr("2e-3 * exp(-1000/(R*T)) * p_Methane / 1e5 + c_Hydrogen + y_Hydrogen");
  EXPECT_TRUE(e.expression_.UsesStateVariables());
  EXPECT_TRUE(e.UsesComponent(CH4));
  EXPECT_TRUE(e.UsesComponent(H2));
  EXPECT_FALSE(e.UsesComponent(H2O));
  const State s{700.0, 2e5, {{CH4, 0.4}, {H2, 0.6}}};
  const double c = s.P / (R * s.T);
  const double expected = 2e-3 * std::exp(-1000 / (R * 700.0)) * 0.4 * 2 + 0.6 * c + 0.6;
  EXPECT_NEAR(r.IntrinsicRate(At(s)), expected, 1e-12 * expected);

  // Outside a run the expression still parses (placeholder state).
  EvalContext ctx;
  r.EvaluateExpressions(ctx);
  EXPECT_TRUE(e.expression_.LastError().empty()) << e.expression_.LastError();
}

TEST(ReactionKinds, LowerCaseUserVariablesAreNotComponentNames) {
  ExpressionParam p;
  p.SetExpr("2 * p_feed");
  EXPECT_FALSE(p.UsesStateVariables());
  p.SetExpr("2 * y_Water");
  EXPECT_TRUE(p.UsesStateVariables());
}

TEST(ReactionKinds, SetRateTypeKeepsSharedParameters) {
  Reaction r("WGS", 1);
  r.SetStoichiometry({{CO, -1}, {H2O, -1}, {CO2, 1}, {H2, 1}});
  r.SetRateType(ReactionRateType::PowerLaw);
  auto& pl = dynamic_cast<PowerLawRate&>(*r.rate_);
  pl.pre_exp_factor_.SetExpr("12.5");
  pl.activation_energy_ = 4.2e4;
  pl.order_[H2] = -1.0;
  r.SetRateType(ReactionRateType::Lhhw);
  auto& l = dynamic_cast<LhhwRate&>(*r.rate_);
  EXPECT_EQ(l.pre_exp_factor_.GetExpr(), "12.5");
  EXPECT_DOUBLE_EQ(l.activation_energy_.Get(), 4.2e4);
  EXPECT_DOUBLE_EQ(l.order_.at(H2).Get(), -1.0);
  EXPECT_DOUBLE_EQ(l.order_.at(CO).Get(), 1.0);
}

TEST(ReactionKinds, FieldsListEveryLeaf) {
  Reaction r1 = Smr();
  std::vector<std::string> leaves;
  for (auto& [leaf, p] : ReactionRateFields(r1.rate_.get())) leaves.push_back(leaf);
  auto has = [&](const std::string& s) { return std::find(leaves.begin(), leaves.end(), s) != leaves.end(); };
  for (const char* leaf : {"pre_exp_factor", "activation_energy", "order[Water]", "order[Hydrogen]",
                           "equilibrium_pre_exp_factor", "equilibrium_energy", "denominator_exponent",
                           "adsorption[1].pre_exp_factor", "adsorption[1].energy",
                           "adsorption[1].exponent[CarbonMonoxide]", "adsorption[4].exponent[Hydrogen]"})
    EXPECT_TRUE(has(leaf)) << leaf;
}

TEST(ReactionHeat, ConstantAndKirchhoff) {
  Reaction r = Smr();
  r.heat_of_reaction_ = 206.1e3;
  r.heat_reference_temperature_ = 298.15;
  const std::map<std::string, double> cp{{CH4, 35.7}, {H2O, 33.6}, {CO, 29.1}, {H2, 28.8}};
  const double dcp = -35.7 - 33.6 + 29.1 + 3 * 28.8;
  EXPECT_NEAR(r.DeltaHeatCapacity(cp), dcp, 1e-12);
  EXPECT_EQ(r.heat_mode_, HeatOfReactionMode::Kirchhoff);  // the default
  EXPECT_NEAR(r.HeatOfReaction(800.0, cp), 206.1e3 + dcp * (800.0 - 298.15), 1e-6);
  r.heat_mode_ = HeatOfReactionMode::Constant;
  EXPECT_DOUBLE_EQ(r.HeatOfReaction(800.0, cp), 206.1e3);
}

TEST(ReactionElements, BalancedImbalancedAndUnknown) {
  Reaction r = Smr();
  EXPECT_TRUE(r.CheckElements().Balanced());

  r.SetStoichiometry({{CH4, -1}, {H2O, -1}, {CO, 1}, {H2, 2}});  // one H2 short
  const ElementBalance bad = r.CheckElements();
  ASSERT_TRUE(bad.checked);
  ASSERT_EQ(bad.imbalance.size(), 1u);
  EXPECT_DOUBLE_EQ(bad.imbalance.at("H"), -2.0);

  r.SetStoichiometry({{"Air", -1}, {CO2, 1}});
  const ElementBalance unknown = r.CheckElements();
  EXPECT_FALSE(unknown.checked);
  ASSERT_EQ(unknown.unknown.size(), 1u);
  EXPECT_EQ(unknown.unknown[0], "Air");
}

}  // namespace
}  // namespace adsim
