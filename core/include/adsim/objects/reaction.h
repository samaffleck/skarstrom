/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Heterogeneous reactions on a catalyst solid.
//
// Reactions live on a solid (an Adsorbent of the AdsorbentDirectory) beside its
// adsorption physics: Adsorbent::reaction_source_. A catalyst is a solid with
// reactions (usually with inert isotherms); a catalyst mixed into a layer with
// a sorbent by mass fraction makes a sorption-enhanced reactor with no further
// machinery. Rates are per kilogram of the solid that carries them
// [mol/(kg s)], in moles of reaction extent (equations/reaction_kinds.h has the
// rate laws).
//
// This header holds the data model, the checks and the evaluation of a rate
// and a heat of reaction at a given state. The reactions do not yet enter the
// balances of a run.

#include <cstdint>
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

#include "adsim/equations/reaction_kinds.h"
#include "adsim/math/expression_param.h"

namespace adsim {

  class Fluid;
  class AdsorbentDirectory;
  class Bed;

  using ReactionID = std::size_t;

  // Enumerator values are what the deck stores, so never renumber them.
  //   Constant:  dH_r(T) = dH_ref
  //   Kirchhoff: dH_r(T) = dH_ref + sum_i nu_i integral_{T_ref}^{T} Cp_i dT,
  //              with the molar gas heat capacities the energy balance uses
  //              (Fluid::MolarEnthalpyChange; dCp (T - T_ref) for constant
  //              Cp_i). The default.
  enum class HeatOfReactionMode {
    Constant = 0,
    Kirchhoff = 1
  };

  //   Constant: a user effectiveness factor eta (default 1).
  //   Thiele:   eta from a generalised Thiele modulus of the key reactant
  //             (its evaluation is not implemented yet; the fields are kept).
  enum class EffectivenessMode {
    Constant = 0,
    Thiele = 1
  };

  // Phase of a stoichiometry entry. Only the gas phase exists; the field is
  // stored so a later version can add others without a format change.
  enum class Phase {
    Gas = 0
  };

  const char* HeatOfReactionModeString(HeatOfReactionMode m);   // "constant", "kirchhoff"
  const char* EffectivenessModeString(EffectivenessMode m);     // "constant", "thiele"

  class EffectivenessModel {
  public:
    EffectivenessMode mode_ = EffectivenessMode::Constant;
    ExpressionParam effectiveness_{1.0};   // eta [-] (Constant)
    std::string key_component_{};          // key reactant (Thiele)
    ExpressionParam tortuosity_{3.0};      // particle tortuosity [-] (Thiele)

    void EvaluateExpressions(const EvalContext& ctx) {
      effectiveness_.Evaluate(ctx);
      tortuosity_.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/) {
      ar(cereal::make_nvp("Mode", mode_),
         cereal::make_nvp("Effectiveness", effectiveness_),
         cereal::make_nvp("Key_Component", key_component_),
         cereal::make_nvp("Tortuosity", tortuosity_));
    }
  };

  // The element balance of a reaction (Reaction::CheckElements).
  struct ElementBalance {
    // False when some component has no element formula (the balance was not
    // checked; `unknown` lists those components).
    bool checked = false;
    std::vector<std::string> unknown{};
    // sum_i nu_i * (atoms of element e in i), for every element whose sum is
    // not zero within the tolerance. Empty when balanced.
    std::map<std::string, double> imbalance{};
    bool Balanced() const { return checked && imbalance.empty(); }
  };

  class Reaction {
  public:
    Reaction() = default;
    Reaction(std::string name, ReactionID id) : name_(std::move(name)), id_(id) {}
    Reaction(Reaction&&) = default;
    Reaction& operator=(Reaction&&) = default;
    Reaction(const Reaction&) = delete;
    Reaction& operator=(const Reaction&) = delete;

    std::string name_ = "Reaction";
    ReactionID id_ = 0;                               // stable; registry path reaction[ID]
    std::map<std::string, double> nu_{};              // stoichiometry; < 0 reactant, > 0 product
    std::unique_ptr<IReactionRate> rate_{};
    RateBasis basis_ = RateBasis::PartialPressure;
    ExpressionParam heat_of_reaction_{0.0};           // dH_r at T_ref [J/mol of extent], < 0 exothermic
    HeatOfReactionMode heat_mode_ = HeatOfReactionMode::Kirchhoff;
    ExpressionParam heat_reference_temperature_{298.15};  // T_ref [K]
    EffectivenessModel effectiveness_{};

    // Replace the stoichiometry. Exact-zero coefficients are dropped; the rate
    // gets default orders for components it has none for.
    void SetStoichiometry(const std::map<std::string, double>& nu);
    // A new rate of `type`; parameters whose leaf (ReactionRateFields) the
    // new layout shares keep their expression text.
    void SetRateType(ReactionRateType type);

    // Every component the reaction names: the stoichiometry and every
    // per-component rate entry (orders, adsorption exponents).
    std::vector<std::string> Components() const;
    // True when `component` matters to the reaction: a nonzero coefficient,
    // order or adsorption exponent (or the expression kind names it).
    bool UsesComponent(const std::string& component) const;
    // Erase the exact-zero entries of a component the reaction does not use.
    // False (nothing changed) when UsesComponent.
    bool RemoveComponent(const std::string& component);
    // True when every component the reaction uses is in `components`.
    bool IsActive(const std::vector<std::string>& components) const;
    bool HasReactant() const;
    bool HasProduct() const;

    // Element balance with the element formulas of fluid.h; an element sums to
    // zero when |sum| <= tolerance.
    ElementBalance CheckElements(double tolerance = 1e-9) const;

    // Intrinsic rate per kg of solid [mol/(kg s)] at a state (cached values;
    // call EvaluateExpressions first). 0 without a rate.
    double IntrinsicRate(const ReactionState& state);
    // dCp = sum_i nu_i Cp_i [J/(mol K)] from molar heat capacities (a
    // component absent from `molar_cp` counts as 0).
    double DeltaHeatCapacity(const std::map<std::string, double>& molar_cp) const;
    // dH_r(T) [J/mol] by heat_mode_.
    double HeatOfReaction(double T, const std::map<std::string, double>& molar_cp) const;
    // As above with the fluid's heat capacities, integrated from T_ref to T
    // (Fluid::MolarEnthalpyChange; the fluid must be initialised).
    double HeatOfReaction(double T, const Fluid& fluid) const;

    void EvaluateExpressions(const EvalContext& ctx);

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/);
  };

  class ReactionSource {
  public:
    ReactionSource() = default;
    ReactionSource(ReactionSource&&) = default;
    ReactionSource& operator=(ReactionSource&&) = default;
    ReactionSource(const ReactionSource&) = delete;
    ReactionSource& operator=(const ReactionSource&) = delete;

    std::vector<Reaction> reactions_{};

    // A new reaction with a fresh ID (1, 2, ...) and a rate of `type`.
    ReactionID Add(const std::string& name, const std::map<std::string, double>& nu,
                   ReactionRateType type = ReactionRateType::PowerLaw);
    // False when there is no reaction `id`.
    bool Remove(ReactionID id);
    Reaction* Find(ReactionID id);
    const Reaction* Find(ReactionID id) const;
    // Throws std::out_of_range when absent.
    Reaction& Get(ReactionID id);
    const Reaction& Get(ReactionID id) const;
    bool Empty() const { return reactions_.empty(); }
    ReactionID NextId() const { return next_id_; }

    // Names of the reactions that use `component` (Reaction::UsesComponent).
    std::vector<std::string> ReactionsUsing(const std::string& component) const;
    // Erase the exact-zero entries of `component` from every reaction.
    // Reactions that use it are left unchanged (callers refuse that case
    // first; RemoveComponentFromDeck does).
    void RemoveComponent(const std::string& component);
    // Seed default orders for every reaction (idempotent). Reactions with a
    // component the fluid lacks are kept, inactive (Reaction::IsActive).
    void Initialise(const Fluid& fluid);
    void EvaluateExpressions(const EvalContext& ctx);

  private:
    ReactionID next_id_ = 0;

    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/) {
      ar(cereal::make_nvp("Reactions", reactions_),
         cereal::make_nvp("Next_Reaction_ID", next_id_));
    }
  };

  // One stoichiometry entry as stored in the deck.
  struct StoichiometryEntry {
    std::string component;
    double nu = 0.0;
    Phase phase = Phase::Gas;

    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/) {
      ar(cereal::make_nvp("Component", component),
         cereal::make_nvp("Nu", nu),
         cereal::make_nvp("Phase", phase));
    }
  };

  template <class Archive>
  void Reaction::serialize(Archive& ar, std::uint32_t const /*version*/) {
    std::vector<StoichiometryEntry> stoichiometry;
    if constexpr (!Archive::is_loading::value) {
      for (const auto& [c, v] : nu_) stoichiometry.push_back({c, v, Phase::Gas});
    }
    ar(cereal::make_nvp("Name", name_),
       cereal::make_nvp("ID", id_),
       cereal::make_nvp("Stoichiometry", stoichiometry),
       cereal::make_nvp("Rate", rate_),
       cereal::make_nvp("Basis", basis_),
       cereal::make_nvp("Heat_Of_Reaction", heat_of_reaction_),
       cereal::make_nvp("Heat_Mode", heat_mode_),
       cereal::make_nvp("Heat_Reference_Temperature", heat_reference_temperature_),
       cereal::make_nvp("Effectiveness", effectiveness_));
    if constexpr (Archive::is_loading::value) {
      nu_.clear();
      for (const auto& e : stoichiometry) {
        if (e.phase != Phase::Gas)
          throw cereal::Exception("reaction '" + name_ + "': unsupported phase for " + e.component);
        nu_[e.component] = e.nu;
      }
    }
  }

  // ---- deck checks and warnings ------------------------------------------

  // Errors (appended to err as "Error: ...\n" lines) for the reactions of
  // every solid that a layer uses: a component that is not in the fluid, no
  // reactant or no product, an element imbalance with known formulas, a
  // non-finite or negative pre-exponential factor (a non-positive equilibrium
  // one for a reversible rate), a Thiele effectiveness without a key reactant,
  // a missing rate. Solids in no layer are not checked (ReactionWarnings
  // reports them). Expressions must already be evaluated. False on an error.
  bool ReactionChecks(const AdsorbentDirectory& dir, const Bed& bed, const Fluid& fluid, std::string& err);

}  // namespace adsim

CEREAL_CLASS_VERSION(adsim::EffectivenessModel, 1);
CEREAL_CLASS_VERSION(adsim::StoichiometryEntry, 1);
CEREAL_CLASS_VERSION(adsim::Reaction, 1);
CEREAL_CLASS_VERSION(adsim::ReactionSource, 1);
