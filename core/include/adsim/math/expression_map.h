/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

#include "adsim/math/expression_param.h"

#include <string>
#include <unordered_map>

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/string.hpp>
#include <cereal/cereal.hpp>

namespace adsim {

// A table of ExpressionParams keyed by name, e.g. an inlet or initial mole
// fraction per component: each entry is a number, a formula of user variables
// or a value with an inline [unit], like any other ExpressionParam.
using ExpressionMap = std::unordered_map<std::string, ExpressionParam>;

inline void EvaluateExpressionMap(ExpressionMap& m, const EvalContext& ctx) {
  for (auto& [key, p] : m) p.Evaluate(ctx);
}

// Serialises `m` under `name`. Archives written before the table became
// expression-valued stored plain numbers under the same name; pass
// `legacy_numbers` (the owning class's archive version is older) to load that
// form. Saving always writes the current class version, so it never takes the
// legacy branch.
template <class Archive>
void SerializeExpressionMap(Archive& ar, const char* name, ExpressionMap& m, bool legacy_numbers) {
  if (legacy_numbers) {
    std::unordered_map<std::string, double> numbers;
    ar(cereal::make_nvp(name, numbers));
    m.clear();
    for (const auto& [key, value] : numbers) m.emplace(key, ExpressionParam(value));
  } else {
    ar(cereal::make_nvp(name, m));
  }
}

} // namespace adsim
