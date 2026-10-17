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
    IsothermType::TOTH_2, IsothermType::ARRHENIUS_SIPS});

// A kind missing from this list is rejected by KindFromString. Keep the list a permutation of the
// enum: one entry per IsothermType, each exactly once.
static_assert(kAllIsothermTypes.size() ==
                  static_cast<std::size_t>(IsothermType::ARRHENIUS_SIPS) + 1,
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
double EvalIsotherm(IIsotherm* iso, const std::string& component, double T, double P, double yi);

}  // namespace adsim
