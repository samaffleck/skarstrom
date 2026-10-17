/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

#include "adsim/math/user_variables.h"

#include <string>

namespace adsim {

// Conversion from display unit to SI base: base = ui * mul + add
struct UnitXform {
    double mul = 1.0;
    double add = 0.0;
    double ToBase(double ui)   const { return ui * mul + add; }
    double FromBase(double si) const { return (si - add) / mul; }
};

// Checks whether `expr` ends with a trailing `[symbol]` unit suffix (COMSOL style).
// If found: sets `math_out` to the expression without the bracket, sets `xform_out`
// to the conversion from that unit to SI base, optionally writes the unit's
// dimensional category to `category_out`, and returns true.
// If not found or symbol is unknown: leaves outputs unchanged, returns false.
bool ParseUnitSuffix(const std::string& expr,
                     std::string&       math_out,
                     UnitXform&         xform_out,
                     UnitCategory*      category_out = nullptr);

// Returns the display symbol for a unit given its string key (used in error messages).
// Returns nullptr if the symbol is not in the dictionary.
const char* LookupUnitSymbol(const std::string& key);

} // namespace adsim
