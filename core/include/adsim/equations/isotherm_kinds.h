/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Isotherm kinds by name, and generic access to an isotherm's parameters.
//
// Every IsothermType has a stable, lower-case string tag ("sips", "toth",
// "dual_site_langmuir_p", ...) that is independent of cereal's polymorphic type
// names. The tags are what the Python package (skarstrom.isotherm_kinds()),
// the isotherm fitter (adsim/solver/isotherm_fit.h) and front ends use to name
// a model. IsothermFields() is the single place that knows each kind's
// parameter layout, and EvalIsotherm() evaluates a single component's loading
// through the same adsim::IsothermFunction free functions the solver uses.

#include <array>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "adsim/equations/isotherm.h"
#include "adsim/math/expression_param.h"

namespace adsim {

// Stable kind tag for every IsothermType. Independent of cereal; "unknown" for
// an out-of-range value.
const char* IsothermKindString(IsothermType t);

// Ordered IsothermType list so the registry enumerates all kinds deterministically.
inline constexpr auto kAllIsothermTypes = std::to_array<IsothermType>({
    IsothermType::INERT, IsothermType::HENRY,
    IsothermType::SINGLE_SITE_LANGMUIR_P, IsothermType::DUAL_SITE_LANGMUIR_P,
    IsothermType::SINGLE_SITE_LANGMUIR_C, IsothermType::DUAL_SITE_LANGMUIR_C,
    IsothermType::EXTENDED_SINGLE_SITE_LANGMUIR_P,
    IsothermType::EXTENDED_DUAL_SITE_LANGMUIR_P,
    IsothermType::EXTENDED_SINGLE_SITE_LANGMUIR_C,
    IsothermType::EXTENDED_DUAL_SITE_LANGMUIR_C,
    IsothermType::ANTI_LANGMUIR, IsothermType::FREUNDLICH,
    IsothermType::SIPS, IsothermType::SIPS_2,
    IsothermType::EXTENDED_SIPS_2, IsothermType::DUAL_SITE_SIPS,
    IsothermType::REDLICH_PETERSON, IsothermType::TOTH,
    IsothermType::TOTH_2, IsothermType::ARRHENIUS_SIPS,
    IsothermType::EXTENDED_LANGMUIR_FREUNDLICH,
    IsothermType::GAB, IsothermType::GAB_YOUNG,
    IsothermType::TOTH_WADST, IsothermType::TOTH_HUMID,
    IsothermType::CUBIC_LANGMUIR});

// A kind missing from this list is rejected by KindFromString. Keep the list a permutation of the
// enum: one entry per IsothermType, each exactly once.
static_assert(kAllIsothermTypes.size() ==
                  static_cast<std::size_t>(IsothermType::CUBIC_LANGMUIR) + 1,
              "kAllIsothermTypes must list every IsothermType");
static_assert(
    [] {
      std::array<bool, kAllIsothermTypes.size()> seen{};
      for (auto t : kAllIsothermTypes) {
        const auto i = static_cast<std::size_t>(t);
        if (i >= seen.size() || seen[i]) return false;
        seen[i] = true;
      }
      return true;
    }(),
    "kAllIsothermTypes must list each IsothermType exactly once");

inline constexpr const auto& AllIsothermTypes() { return kAllIsothermTypes; }

// The single place that knows each variant's parameter layout: returns pointers
// to a component's ExpressionParam fields for the supported kinds. Used for both
// read (Get/GetExpr) and write (SetExpr). Empty vector => kind not yet supported
// for param editing (or the component has no entry). Adding a kind is one case.
std::vector<std::pair<std::string, ExpressionParam*>> IsothermFields(IIsotherm* iso, const std::string& component);

// True when IsothermFields covers this kind (INERT counts: supported, zero params).
bool IsSupportedIsothermKind(IsothermType t);

// True for the competitive kinds, whose loadings couple the components through
// a shared denominator when the adsorbent's mixing rule is none: the extended
// Langmuir kinds (single- and dual-site, both bases), extended Sips 2,
// Arrhenius Sips and extended Langmuir-Freundlich. Under the IAST mixing rule
// (porous_media.h) these are taken in their pure-component limit, so a front
// end offers the non-competitive kinds there instead.
bool IsCompetitiveIsothermKind(IsothermType t);

// True for the kinds whose loading depends on the adsorbed loading of water on
// the same adsorbent (component kCoadsorbedWaterComponent): toth_wadst and
// toth_humid. They have no IAST form (water changes the CO2 sites, not a
// competitive denominator), so the IAST mixing rule rejects them.
bool IsWaterCoupledIsothermKind(IsothermType t);

// True when this component's isotherm reads the adsorbed water loading: a
// water-coupled kind (above), or a dual-site Langmuir kind whose optional
// water factors are set (DslWaterFactors::Active, isotherm.h). A dual-site
// Langmuir with its factors at 0 is the dry kind and returns false.
bool IsothermUsesWaterLoading(IIsotherm* iso, const std::string& component);

// The kind with the same parameter layout in the other family: a competitive
// kind's pure-component counterpart, or a pure kind's competitive form
// (single_site_langmuir_p <-> extended_single_site_langmuir_p, and likewise for
// the dual-site and concentration forms, sips_2 <-> extended_sips_2). Switching
// a component to its paired kind keeps every parameter. False when there is
// none (Arrhenius Sips and extended Langmuir-Freundlich have no pure kind with
// their layout; Toth, Freundlich and the rest have no competitive form).
bool PairedIsothermKind(IsothermType t, IsothermType& out);

// Reverse of IsothermKindString: kind tag -> IsothermType. False if unknown.
bool KindFromString(const std::string& kind, IsothermType& out);

// Guarantee a component has a params entry in `iso` (seeded with the param
// struct's default member initializers). Needed after switching a component to
// an isotherm type whose object already existed (shared by another component),
// because SetIsothermType only seeds params when it CREATES the object.
void EnsureIsothermComponent(IIsotherm* iso, const std::string& component);

// Evaluate the equilibrium loading q [mol/kg] for one component at (T [K], total
// pressure P [Pa], mole fraction yi) using the direct IsothermFunction free
// functions. For a single-component (pure) plot pass yi=1 so partial pressure =
// P. Extended/competitive variants reuse the base function (pure-component limit).
// The GAB kinds take the component's saturation pressure (SaturationPressure). The
// water-coupled kinds (IsWaterCoupledIsothermKind) take the adsorbed water
// loading q_water [mol/kg]; 0 gives their dry isotherm. Other kinds ignore it.
double EvalIsotherm(IIsotherm* iso, const std::string& component, double T, double P, double yi,
                    double q_water = 0.0);

// Isosteric heat of adsorption q_st [J/mol, positive for exothermic
// adsorption] implied by the isotherm's temperature dependence, at the loading
// of the pure component at partial pressure P yi and temperature T:
//   q_st = R T^2 (d ln p / dT)_q = -R T^2 (dq/dT)_p / (dq/d ln p)_T,
// both partials by central differences of EvalIsotherm. For Langmuir with
// b = b0 exp(Q / (R T)) this is Q at every loading; for the GAB kinds it varies
// with loading and includes the heat of condensation through p_sat(T). The
// energy balance does not call this: it takes the adsorbent's constant
// heat_of_adsorption[component]; this gives the value to put there (or checks
// one). NaN when P yi <= 0. q_water as for EvalIsotherm (held fixed).
double IsostericHeat(IIsotherm* iso, const std::string& component, double T, double P, double yi,
                     double q_water = 0.0);

}  // namespace adsim
