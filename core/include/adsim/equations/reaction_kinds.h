/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Rate laws of heterogeneous reactions on a catalyst solid (objects/reaction.h).
//
// A rate law gives the reaction rate r per kilogram of the solid that carries
// the reaction [mol/(kg s)], in moles of reaction extent: component i is
// produced at nu_i * r. Every kind has a stable lower-case tag (the builder's
// kind name) and a parameter layout listed by ReactionRateFields, the single
// place the registry and the builder read it from.
//
// Activities. A rate is written in the activities a_i of the components:
//   RateBasis::PartialPressure  a_i = p_i = y_i P  in bar
//   RateBasis::Concentration    a_i = C_i = y_i P / (R T)  in mol/m3
// Bar, not Pa, because the kinetics literature (Xu and Froment 1989 and most
// papers after it) tabulates rate and adsorption constants against bar.
//
// Kinds (k(T) = A exp(-E / (R T)), as the Arrhenius mass-transfer coefficient):
//   power_law             r = k(T) prod_i a_i^n_i
//   reversible_power_law  r = k(T) [prod_i a_i^n_i - (1/K(T)) prod_i a_i^(n_i + nu_i)]
//                         = k(T) prod_i a_i^n_i (1 - Q/K(T)),  Q = prod_i a_i^nu_i,
//                         K(T) = A_K exp(-E_K / (R T)) (user-entered; no
//                         thermochemistry table)
//   lhhw                  r = k(T) prod_i a_i^n_i (1 - Q/K(T)) / DEN^m,
//                         DEN = 1 + sum_j K_j(T) prod_i a_i^e_ji,
//                         K_j(T) = A_j exp(-E_j / (R T)); the (1 - Q/K) factor
//                         only when `reversible` is set
//   expression            a user expression of T [K], P [Pa], c [mol/m3] and
//                         each component's p_<C> [Pa], c_<C> [mol/m3] and
//                         y_<C> [-], e.g. "1e-3 * exp(-5e4/(R*T)) * p_Methane/1e5"
// Orders n_i default to -nu_i for reactants and 0 for products (an
// elementary reaction).
//
// The reversible factor is evaluated as the difference of the two monomials,
// never through Q, so a zero product (or reactant) activity gives a finite
// rate. Each activity is clipped at 0 (a solver iterate can overshoot below),
// and an activity raised to a NEGATIVE power is floored at kMinActivityBar
// (partial-pressure basis) or kMinActivityConcentration (concentration basis),
// so a component that is absent (H2 at a steam-reformer inlet, whose order is
// -2.5 in the Xu-Froment rate) gives a large but finite rate.

#include <array>
#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <cereal/access.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include "adsim/math/expression_param.h"
#include "adsim/math/user_variables.h"

namespace adsim {

  // Enumerator values are what the deck stores, so never renumber them.
  enum class ReactionRateType {
    PowerLaw = 0,
    ReversiblePowerLaw = 1,
    Lhhw = 2,
    Expression = 3
  };

  // The activity a rate law is written in (see the file comment). Enumerator
  // values are what the deck stores.
  enum class RateBasis {
    PartialPressure = 0,
    Concentration = 1
  };

  inline constexpr double kMinActivityBar = 1e-10;            // [bar]
  inline constexpr double kMinActivityConcentration = 1e-6;   // [mol/m3]

  // Stable kind tag; "unknown" for an out-of-range value.
  const char* ReactionRateKindString(ReactionRateType t);
  // Reverse of ReactionRateKindString. False if unknown.
  bool ReactionRateKindFromString(const std::string& kind, ReactionRateType& out);
  // "partial_pressure", "concentration".
  const char* RateBasisString(RateBasis b);

  inline constexpr auto kAllReactionRateKinds = std::to_array<ReactionRateType>({
      ReactionRateType::PowerLaw, ReactionRateType::ReversiblePowerLaw,
      ReactionRateType::Lhhw, ReactionRateType::Expression});
  static_assert(kAllReactionRateKinds.size() == static_cast<std::size_t>(ReactionRateType::Expression) + 1,
                "kAllReactionRateKinds must list every ReactionRateType");

  // The local gas state a rate is evaluated at. `y` maps a component to its
  // mole fraction; components absent from it count as 0.
  struct ReactionState {
    double T = 298.15;                                  // [K]
    double P = 1e5;                                     // [Pa]
    const std::map<std::string, double>* y = nullptr;   // mole fractions [-]
    const UserVariableRegistry* registry = nullptr;     // for expression parameters
  };

  // Activity of a component with mole fraction y at (T, P) on a basis.
  double Activity(RateBasis basis, double y, double T, double P);
  // The mole fraction whose activity at (T, P) is the basis' floor
  // (kMinActivityBar or kMinActivityConcentration): below it a negative
  // power of the activity no longer changes.
  double ActivityFloorMoleFraction(RateBasis basis, double T, double P);

  class IReactionRate {
  public:
    virtual ~IReactionRate() = default;
    virtual ReactionRateType GetType() const = 0;

    // Rate per kg of solid [mol/(kg s)] of a reaction with stoichiometry `nu`
    // on `basis`, from the parameters' cached values (call
    // EvaluateExpressions first). The expression kind evaluates its formula
    // at the state.
    virtual double Rate(const ReactionState& state, RateBasis basis,
                        const std::map<std::string, double>& nu) = 0;

    // Seed a default order (and, for LHHW, nothing else) for every component
    // of `nu` that has none: -nu_i for a reactant, 0 for a product.
    virtual void SeedComponents(const std::map<std::string, double>& /*nu*/) {}
    // True when a parameter of this rate depends on `component` (a nonzero
    // order or adsorption exponent; the expression kind: its formula names
    // p_, c_ or y_<component>).
    virtual bool UsesComponent(const std::string& /*component*/) const { return false; }
    // Erase every per-component entry of `component`.
    virtual void RemoveComponent(const std::string& /*component*/) {}
    // The smallest nonzero power to which the rate raises the activity of
    // `component` (cached values): over the orders n_i, the exponents
    // n_i + nu_i of the reverse monomial (reversible kinds) and the exponents
    // e_ji of the LHHW denominator terms; +infinity when the rate does not
    // depend on the component. Below 1 the derivative of the rate in the
    // component's mole fraction is unbounded as it goes to zero (below 0 the
    // rate itself is, up to the activity floor), which the solver's Jacobian
    // treats specially (jacobian.h). The expression kind cannot tell: 1 when
    // its formula names the component, else +infinity.
    virtual double MinExponent(const std::string& component,
                               const std::map<std::string, double>& nu) const;

    virtual void EvaluateExpressions(const EvalContext& ctx) = 0;
    virtual void CompileForRun(const UserVariableRegistry* /*registry*/) {}

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& /*ar*/, std::uint32_t const /*version*/) {}
  };

  // r = k(T) prod a_i^n_i.
  class PowerLawRate : public IReactionRate {
  public:
    ReactionRateType GetType() const override { return ReactionRateType::PowerLaw; }
    double MinExponent(const std::string& component, const std::map<std::string, double>& nu) const override;
    double Rate(const ReactionState& state, RateBasis basis, const std::map<std::string, double>& nu) override;
    void SeedComponents(const std::map<std::string, double>& nu) override;
    bool UsesComponent(const std::string& component) const override;
    void RemoveComponent(const std::string& component) override;
    void EvaluateExpressions(const EvalContext& ctx) override;

    ExpressionParam pre_exp_factor_{1.0};      // A [mol/(kg s activity^n)]
    ExpressionParam activation_energy_{0.0};   // E [J/mol]
    std::map<std::string, ExpressionParam> order_{};  // n_i [-]

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/) {
      ar(cereal::virtual_base_class<IReactionRate>(this),
         cereal::make_nvp("Pre_Exponential_Factor", pre_exp_factor_),
         cereal::make_nvp("Activation_Energy", activation_energy_),
         cereal::make_nvp("Orders", order_));
    }
  };

  // r = k(T) [prod a_i^n_i - (1/K(T)) prod a_i^(n_i + nu_i)].
  class ReversiblePowerLawRate : public IReactionRate {
  public:
    ReactionRateType GetType() const override { return ReactionRateType::ReversiblePowerLaw; }
    double MinExponent(const std::string& component, const std::map<std::string, double>& nu) const override;
    double Rate(const ReactionState& state, RateBasis basis, const std::map<std::string, double>& nu) override;
    void SeedComponents(const std::map<std::string, double>& nu) override;
    bool UsesComponent(const std::string& component) const override;
    void RemoveComponent(const std::string& component) override;
    void EvaluateExpressions(const EvalContext& ctx) override;

    ExpressionParam pre_exp_factor_{1.0};
    ExpressionParam activation_energy_{0.0};
    std::map<std::string, ExpressionParam> order_{};
    ExpressionParam equilibrium_pre_exp_factor_{1.0};  // A_K [activity^sum(nu)]
    ExpressionParam equilibrium_energy_{0.0};          // E_K [J/mol]

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/) {
      ar(cereal::virtual_base_class<IReactionRate>(this),
         cereal::make_nvp("Pre_Exponential_Factor", pre_exp_factor_),
         cereal::make_nvp("Activation_Energy", activation_energy_),
         cereal::make_nvp("Orders", order_),
         cereal::make_nvp("Equilibrium_Pre_Exponential_Factor", equilibrium_pre_exp_factor_),
         cereal::make_nvp("Equilibrium_Energy", equilibrium_energy_));
    }
  };

  // One term K_j(T) prod a_i^e_ji of an LHHW denominator.
  struct LhhwAdsorptionTerm {
    ExpressionParam pre_exp_factor_{1.0};   // A_j [1/activity^sum(e)]
    ExpressionParam energy_{0.0};           // E_j [J/mol] (an adsorption enthalpy, < 0 usually)
    std::map<std::string, ExpressionParam> exponent_{};  // e_ji [-]

    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/) {
      ar(cereal::make_nvp("Pre_Exponential_Factor", pre_exp_factor_),
         cereal::make_nvp("Energy", energy_),
         cereal::make_nvp("Exponents", exponent_));
    }
  };

  // Langmuir-Hinshelwood-Hougen-Watson:
  // r = k(T) prod a_i^n_i (1 - Q/K(T)) / (1 + sum_j K_j(T) prod a_i^e_ji)^m.
  class LhhwRate : public IReactionRate {
  public:
    ReactionRateType GetType() const override { return ReactionRateType::Lhhw; }
    double MinExponent(const std::string& component, const std::map<std::string, double>& nu) const override;
    double Rate(const ReactionState& state, RateBasis basis, const std::map<std::string, double>& nu) override;
    void SeedComponents(const std::map<std::string, double>& nu) override;
    bool UsesComponent(const std::string& component) const override;
    void RemoveComponent(const std::string& component) override;
    void EvaluateExpressions(const EvalContext& ctx) override;

    // The denominator DEN = 1 + sum_j K_j prod a^e (for tests and front ends).
    double Denominator(const ReactionState& state, RateBasis basis) const;

    ExpressionParam pre_exp_factor_{1.0};
    ExpressionParam activation_energy_{0.0};
    std::map<std::string, ExpressionParam> order_{};
    bool reversible_ = true;
    ExpressionParam equilibrium_pre_exp_factor_{1.0};
    ExpressionParam equilibrium_energy_{0.0};
    ExpressionParam denominator_exponent_{1.0};       // m [-]
    std::vector<LhhwAdsorptionTerm> adsorption_{};     // j = 1, 2, ... (index j - 1)

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/) {
      ar(cereal::virtual_base_class<IReactionRate>(this),
         cereal::make_nvp("Pre_Exponential_Factor", pre_exp_factor_),
         cereal::make_nvp("Activation_Energy", activation_energy_),
         cereal::make_nvp("Orders", order_),
         cereal::make_nvp("Reversible", reversible_),
         cereal::make_nvp("Equilibrium_Pre_Exponential_Factor", equilibrium_pre_exp_factor_),
         cereal::make_nvp("Equilibrium_Energy", equilibrium_energy_),
         cereal::make_nvp("Denominator_Exponent", denominator_exponent_),
         cereal::make_nvp("Adsorption_Terms", adsorption_));
    }
  };

  // A user formula of the state (see the file comment). The per-component
  // names are bound for the components of the reaction's stoichiometry.
  class ExpressionRate : public IReactionRate {
  public:
    ReactionRateType GetType() const override { return ReactionRateType::Expression; }
    double Rate(const ReactionState& state, RateBasis basis, const std::map<std::string, double>& nu) override;
    bool UsesComponent(const std::string& component) const override;
    double MinExponent(const std::string& component, const std::map<std::string, double>& nu) const override;
    void EvaluateExpressions(const EvalContext& ctx) override { expression_.Evaluate(ctx); }

    ExpressionParam expression_{0.0};   // r [mol/(kg s)]

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/) {
      ar(cereal::virtual_base_class<IReactionRate>(this),
         cereal::make_nvp("Expression", expression_));
    }
  };

  // A new rate of `type` with default parameters, seeded for `nu`.
  std::unique_ptr<IReactionRate> CreateReactionRate(ReactionRateType type,
                                                    const std::map<std::string, double>& nu = {});

  // The single place that knows each kind's parameter layout: (leaf, field)
  // pairs, leaves as the registry writes them after "rate.":
  //   pre_exp_factor, activation_energy, order[C],
  //   equilibrium_pre_exp_factor, equilibrium_energy,      (reversible kinds)
  //   denominator_exponent, adsorption[j].pre_exp_factor,
  //   adsorption[j].energy, adsorption[j].exponent[C]      (lhhw, j from 1)
  //   expression                                           (expression)
  // The LHHW `reversible` switch is a bool, not listed here.
  std::vector<std::pair<std::string, ExpressionParam*>> ReactionRateFields(IReactionRate* rate);

}  // namespace adsim

// Class versions live next to the inline serialize templates: every
// translation unit that saves or loads a deck must see the same numbers.
CEREAL_CLASS_VERSION(adsim::IReactionRate, 1);
CEREAL_CLASS_VERSION(adsim::PowerLawRate, 1);
CEREAL_CLASS_VERSION(adsim::ReversiblePowerLawRate, 1);
CEREAL_CLASS_VERSION(adsim::LhhwRate, 1);
CEREAL_CLASS_VERSION(adsim::ExpressionRate, 1);
CEREAL_CLASS_VERSION(adsim::LhhwAdsorptionTerm, 1);
