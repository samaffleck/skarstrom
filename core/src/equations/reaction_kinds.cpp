/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/equations/reaction_kinds.h"

#include <algorithm>
#include <cmath>
#include <set>

#include <cereal/types/polymorphic.hpp>

#include "adsim/math/constants.h"

namespace adsim {

  const char* ReactionRateKindString(ReactionRateType t) {
    switch (t) {
      case ReactionRateType::PowerLaw: return "power_law";
      case ReactionRateType::ReversiblePowerLaw: return "reversible_power_law";
      case ReactionRateType::Lhhw: return "lhhw";
      case ReactionRateType::Expression: return "expression";
    }
    return "unknown";
  }

  bool ReactionRateKindFromString(const std::string& kind, ReactionRateType& out) {
    for (auto t : kAllReactionRateKinds) {
      if (kind == ReactionRateKindString(t)) {
        out = t;
        return true;
      }
    }
    return false;
  }

  const char* RateBasisString(RateBasis b) {
    return b == RateBasis::Concentration ? "concentration" : "partial_pressure";
  }

  double Activity(RateBasis basis, double y, double T, double P) {
    if (basis == RateBasis::Concentration) {
      const double Tc = std::max(T, constants::kMinDenomTemperature);
      return y * P / (constants::kIdealGasConstant * Tc);  // [mol/m3]
    }
    return y * P * 1e-5;  // [bar]
  }

  namespace {

    double MoleFraction(const ReactionState& s, const std::string& c) {
      if (!s.y) return 0.0;
      auto it = s.y->find(c);
      return it == s.y->end() ? 0.0 : it->second;
    }

    // a^e with the clipping of the file comment: the activity is clipped at 0,
    // and floored at the basis' minimum when e < 0. e == 0 gives exactly 1.
    double Power(double a, double e, RateBasis basis) {
      if (e == 0.0) return 1.0;
      double base = std::max(a, 0.0);
      if (e < 0.0) {
        const double floor = basis == RateBasis::Concentration ? kMinActivityConcentration : kMinActivityBar;
        base = std::max(base, floor);
      }
      return std::pow(base, e);
    }

    // prod_i a_i^e_i over the entries of `exponents` (cached values).
    double Monomial(const ReactionState& s, RateBasis basis,
                    const std::map<std::string, ExpressionParam>& exponents) {
      double prod = 1.0;
      for (const auto& [c, e] : exponents) prod *= Power(Activity(basis, MoleFraction(s, c), s.T, s.P), e.Get(), basis);
      return prod;
    }

    // prod_i a_i^(n_i + nu_i) over the union of the order and stoichiometry keys.
    double ReverseMonomial(const ReactionState& s, RateBasis basis,
                           const std::map<std::string, ExpressionParam>& order,
                           const std::map<std::string, double>& nu) {
      std::set<std::string> keys;
      for (const auto& [c, n] : order) keys.insert(c);
      for (const auto& [c, v] : nu) keys.insert(c);
      double prod = 1.0;
      for (const auto& c : keys) {
        double e = 0.0;
        if (auto it = order.find(c); it != order.end()) e += it->second.Get();
        if (auto it = nu.find(c); it != nu.end()) e += it->second;
        prod *= Power(Activity(basis, MoleFraction(s, c), s.T, s.P), e, basis);
      }
      return prod;
    }

    double Arrhenius(double A, double E, double T) {
      const double Tc = std::max(T, constants::kMinDenomTemperature);
      return A * std::exp(-E / (constants::kIdealGasConstant * Tc));
    }

    void SeedOrders(std::map<std::string, ExpressionParam>& order, const std::map<std::string, double>& nu) {
      for (const auto& [c, v] : nu) order.try_emplace(c, v < 0.0 ? -v : 0.0);
    }

    bool AnyNonzero(const std::map<std::string, ExpressionParam>& m, const std::string& c) {
      auto it = m.find(c);
      return it != m.end() && it->second.Get() != 0.0;
    }

    void EvaluateMap(std::map<std::string, ExpressionParam>& m, const EvalContext& ctx) {
      for (auto& [c, p] : m) p.Evaluate(ctx);
    }

  }  // namespace

  // ---- power law -------------------------------------------------------------

  double PowerLawRate::Rate(const ReactionState& s, RateBasis basis, const std::map<std::string, double>&) {
    return Arrhenius(pre_exp_factor_.Get(), activation_energy_.Get(), s.T) * Monomial(s, basis, order_);
  }
  void PowerLawRate::SeedComponents(const std::map<std::string, double>& nu) { SeedOrders(order_, nu); }
  bool PowerLawRate::UsesComponent(const std::string& c) const { return AnyNonzero(order_, c); }
  void PowerLawRate::RemoveComponent(const std::string& c) { order_.erase(c); }
  void PowerLawRate::EvaluateExpressions(const EvalContext& ctx) {
    pre_exp_factor_.Evaluate(ctx);
    activation_energy_.Evaluate(ctx);
    EvaluateMap(order_, ctx);
  }

  // ---- reversible power law --------------------------------------------------

  double ReversiblePowerLawRate::Rate(const ReactionState& s, RateBasis basis,
                                      const std::map<std::string, double>& nu) {
    const double k = Arrhenius(pre_exp_factor_.Get(), activation_energy_.Get(), s.T);
    const double K = Arrhenius(equilibrium_pre_exp_factor_.Get(), equilibrium_energy_.Get(), s.T);
    return k * (Monomial(s, basis, order_) - ReverseMonomial(s, basis, order_, nu) / K);
  }
  void ReversiblePowerLawRate::SeedComponents(const std::map<std::string, double>& nu) { SeedOrders(order_, nu); }
  bool ReversiblePowerLawRate::UsesComponent(const std::string& c) const { return AnyNonzero(order_, c); }
  void ReversiblePowerLawRate::RemoveComponent(const std::string& c) { order_.erase(c); }
  void ReversiblePowerLawRate::EvaluateExpressions(const EvalContext& ctx) {
    pre_exp_factor_.Evaluate(ctx);
    activation_energy_.Evaluate(ctx);
    EvaluateMap(order_, ctx);
    equilibrium_pre_exp_factor_.Evaluate(ctx);
    equilibrium_energy_.Evaluate(ctx);
  }

  // ---- LHHW ------------------------------------------------------------------

  double LhhwRate::Denominator(const ReactionState& s, RateBasis basis) const {
    double den = 1.0;
    for (const auto& term : adsorption_)
      den += Arrhenius(term.pre_exp_factor_.Get(), term.energy_.Get(), s.T) * Monomial(s, basis, term.exponent_);
    return den;
  }

  double LhhwRate::Rate(const ReactionState& s, RateBasis basis, const std::map<std::string, double>& nu) {
    const double k = Arrhenius(pre_exp_factor_.Get(), activation_energy_.Get(), s.T);
    double driving = Monomial(s, basis, order_);
    if (reversible_) {
      const double K = Arrhenius(equilibrium_pre_exp_factor_.Get(), equilibrium_energy_.Get(), s.T);
      driving -= ReverseMonomial(s, basis, order_, nu) / K;
    }
    return k * driving / std::pow(Denominator(s, basis), denominator_exponent_.Get());
  }
  void LhhwRate::SeedComponents(const std::map<std::string, double>& nu) { SeedOrders(order_, nu); }
  bool LhhwRate::UsesComponent(const std::string& c) const {
    if (AnyNonzero(order_, c)) return true;
    for (const auto& term : adsorption_)
      if (AnyNonzero(term.exponent_, c)) return true;
    return false;
  }
  void LhhwRate::RemoveComponent(const std::string& c) {
    order_.erase(c);
    for (auto& term : adsorption_) term.exponent_.erase(c);
  }
  void LhhwRate::EvaluateExpressions(const EvalContext& ctx) {
    pre_exp_factor_.Evaluate(ctx);
    activation_energy_.Evaluate(ctx);
    EvaluateMap(order_, ctx);
    equilibrium_pre_exp_factor_.Evaluate(ctx);
    equilibrium_energy_.Evaluate(ctx);
    denominator_exponent_.Evaluate(ctx);
    for (auto& term : adsorption_) {
      term.pre_exp_factor_.Evaluate(ctx);
      term.energy_.Evaluate(ctx);
      EvaluateMap(term.exponent_, ctx);
    }
  }

  // ---- expression ------------------------------------------------------------

  double ExpressionRate::Rate(const ReactionState& s, RateBasis, const std::map<std::string, double>& nu) {
    const double Tc = std::max(s.T, constants::kMinDenomTemperature);
    const double c_total = s.P / (constants::kIdealGasConstant * Tc);
    std::vector<ComponentStateValue> comps;
    auto add = [&](const std::string& name, double y) {
      for (const auto& e : comps)
        if (e.name == name) return;
      comps.push_back({name, y * s.P, y * c_total, y});
    };
    if (s.y)
      for (const auto& [c, y] : *s.y) add(c, y);
    for (const auto& [c, v] : nu) add(c, 0.0);
    EvalContext ctx;
    ctx.registry = s.registry;
    ctx.T = s.T;
    ctx.P = s.P;
    ctx.c = c_total;
    ctx.components = &comps;
    expression_.Evaluate(ctx);
    return expression_.Get();
  }

  bool ExpressionRate::UsesComponent(const std::string& c) const {
    const std::string& e = expression_.GetExpr();
    for (const char* prefix : {"p_", "c_", "y_"}) {
      const std::string name = prefix + c;
      for (size_t pos = e.find(name); pos != std::string::npos; pos = e.find(name, pos + 1)) {
        auto ident = [](char ch) { return std::isalnum(static_cast<unsigned char>(ch)) || ch == '_'; };
        const bool start_ok = pos == 0 || !ident(e[pos - 1]);
        const size_t end = pos + name.size();
        const bool end_ok = end >= e.size() || !ident(e[end]);
        if (start_ok && end_ok) return true;
      }
    }
    return false;
  }

  // ---- factory and fields ----------------------------------------------------

  std::unique_ptr<IReactionRate> CreateReactionRate(ReactionRateType type, const std::map<std::string, double>& nu) {
    std::unique_ptr<IReactionRate> r;
    switch (type) {
      case ReactionRateType::PowerLaw: r = std::make_unique<PowerLawRate>(); break;
      case ReactionRateType::ReversiblePowerLaw: r = std::make_unique<ReversiblePowerLawRate>(); break;
      case ReactionRateType::Lhhw: r = std::make_unique<LhhwRate>(); break;
      case ReactionRateType::Expression: r = std::make_unique<ExpressionRate>(); break;
    }
    if (r) r->SeedComponents(nu);
    return r;
  }

  std::vector<std::pair<std::string, ExpressionParam*>> ReactionRateFields(IReactionRate* rate) {
    std::vector<std::pair<std::string, ExpressionParam*>> out;
    if (!rate) return out;
    auto orders = [&](std::map<std::string, ExpressionParam>& m, const std::string& head) {
      for (auto& [c, p] : m) out.emplace_back(head + "[" + c + "]", &p);
    };
    switch (rate->GetType()) {
      case ReactionRateType::PowerLaw:
        if (auto* r = dynamic_cast<PowerLawRate*>(rate)) {
          out.emplace_back("pre_exp_factor", &r->pre_exp_factor_);
          out.emplace_back("activation_energy", &r->activation_energy_);
          orders(r->order_, "order");
        }
        break;
      case ReactionRateType::ReversiblePowerLaw:
        if (auto* r = dynamic_cast<ReversiblePowerLawRate*>(rate)) {
          out.emplace_back("pre_exp_factor", &r->pre_exp_factor_);
          out.emplace_back("activation_energy", &r->activation_energy_);
          orders(r->order_, "order");
          out.emplace_back("equilibrium_pre_exp_factor", &r->equilibrium_pre_exp_factor_);
          out.emplace_back("equilibrium_energy", &r->equilibrium_energy_);
        }
        break;
      case ReactionRateType::Lhhw:
        if (auto* r = dynamic_cast<LhhwRate*>(rate)) {
          out.emplace_back("pre_exp_factor", &r->pre_exp_factor_);
          out.emplace_back("activation_energy", &r->activation_energy_);
          orders(r->order_, "order");
          out.emplace_back("equilibrium_pre_exp_factor", &r->equilibrium_pre_exp_factor_);
          out.emplace_back("equilibrium_energy", &r->equilibrium_energy_);
          out.emplace_back("denominator_exponent", &r->denominator_exponent_);
          for (size_t j = 0; j < r->adsorption_.size(); ++j) {
            const std::string head = "adsorption[" + std::to_string(j + 1) + "].";
            auto& term = r->adsorption_[j];
            out.emplace_back(head + "pre_exp_factor", &term.pre_exp_factor_);
            out.emplace_back(head + "energy", &term.energy_);
            orders(term.exponent_, head + "exponent");
          }
        }
        break;
      case ReactionRateType::Expression:
        if (auto* r = dynamic_cast<ExpressionRate*>(rate)) out.emplace_back("expression", &r->expression_);
        break;
    }
    return out;
  }

}  // namespace adsim

CEREAL_REGISTER_TYPE(adsim::PowerLawRate)
CEREAL_REGISTER_TYPE(adsim::ReversiblePowerLawRate)
CEREAL_REGISTER_TYPE(adsim::LhhwRate)
CEREAL_REGISTER_TYPE(adsim::ExpressionRate)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IReactionRate, adsim::PowerLawRate)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IReactionRate, adsim::ReversiblePowerLawRate)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IReactionRate, adsim::LhhwRate)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IReactionRate, adsim::ExpressionRate)
