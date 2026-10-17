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
// The symbol is any unit expression (ParseUnitExpression), so "[W/(m^2*K)]"
// and "[kJ/(kg-K)]" work too. The category follows from the dimension
// (Custom when no category has it); %, ppm and mol/mol are MoleFraction.
bool ParseUnitSuffix(const std::string& expr,
                     std::string&       math_out,
                     UnitXform&         xform_out,
                     UnitCategory*      category_out = nullptr);

// Physical dimension as exponents of the SI base quantities the engine uses:
// length, mass, time, temperature and amount of substance.
struct Dimension {
    int length = 0;
    int mass = 0;
    int time = 0;
    int temperature = 0;
    int amount = 0;
    bool operator==(const Dimension&) const = default;
    bool IsDimensionless() const { return *this == Dimension{}; }
};

// A parsed unit: SI value = value * factor + offset. The offset is nonzero
// only for a lone absolute temperature unit (°C, °F).
struct UnitInfo {
    double factor = 1.0;
    double offset = 0.0;
    Dimension dimension;
};

// Parses a unit expression such as "kJ/mol", "W/(m^2*K)", "W/m2-K",
// "J/(kg·K)", "kg m^-3", "1/bar", "mmol/g" or "Nm3/h".
//   * Atoms: any unit symbol the LLNL units library knows, with SI prefixes
//     (kPa, mbar, µm, min, h, atm, psi, degC, kWh...), read strictly so a typo
//     is an error. C, F and A are °C, °F and the ångström; ccm, SCCM, SCFM and
//     Nm3 are gas-flow units at 0 °C and 1 atm. Only units of length, mass,
//     time, temperature and amount of substance are accepted.
//   * Exponents: m2, m^2, m**2, m², s^-1.
//   * Multiplication: '*', '.', '·', '-' or a space. Division: everything after
//     a '/' up to the next '/' is in the denominator, the engineering reading
//     of "W/m2-K" = W/(m²·K); parentheses group.
//   * "-", "1" and "" are dimensionless; "%" is 0.01 and "ppm" 1e-6.
// A temperature unit with an offset counts as absolute only on its own
// ("degC"); inside a compound unit it is a temperature difference.
// Returns false (and, when given, an explanation in `error`) for anything else.
bool ParseUnitExpression(const std::string& text, UnitInfo& out, std::string* error = nullptr);

// "time", "pressure", "heat transfer coefficient"... for the dimensions the
// engine's parameters use; otherwise the base-unit form ("kg·m^-1·s^-2").
std::string DimensionName(const Dimension& dimension);

// Converts `value` from unit `from` to unit `to`. Throws std::invalid_argument
// when either unit does not parse or their dimensions differ.
double ConvertUnit(double value, const std::string& from, const std::string& to);

// Splits a quantity written as "<number> <unit>" ("120 s", "2 min", "25 degC",
// "1.5e5 Pa", "15 %") or "<number> [<unit>]" into the number and the unit text.
// Returns false when `text` is not a single number followed by a parseable
// unit (a plain number, a formula, a unit that does not parse).
bool SplitQuantity(const std::string& text, double& number, std::string& unit);

} // namespace adsim
