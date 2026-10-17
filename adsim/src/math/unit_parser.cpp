#include "adsim/math/unit_parser.h"

#include <regex>
#include <unordered_map>

namespace adsim {

// Each entry pairs the numeric transform with the dimensional category the
// symbol belongs to. Units that don't fit any UnitCategory enum (power,
// inverse pressure, etc.) are tagged Custom — they are only accepted when the
// variable's own category is Custom.
struct UnitEntry {
    UnitXform    xform;
    UnitCategory category;
};

// ---------------------------------------------------------------------------
// Unit dictionary — symbol/alias  →  UnitEntry  (base = ui * mul + add)
// Covers the unit symbols the front ends display, plus ASCII-friendly aliases.
// ---------------------------------------------------------------------------
static const std::unordered_map<std::string, UnitEntry>& UnitMap() {
    static const std::unordered_map<std::string, UnitEntry> kMap = {

        // ---- PRESSURE (base: Pa) ----
        {"Pa",      {{1.0,                0.0}, UnitCategory::Pressure}},
        {"kPa",     {{1e3,                0.0}, UnitCategory::Pressure}},
        {"MPa",     {{1e6,                0.0}, UnitCategory::Pressure}},
        {"bar",     {{1e5,                0.0}, UnitCategory::Pressure}},
        {"atm",     {{101325.0,           0.0}, UnitCategory::Pressure}},
        {"psi",     {{6894.757293168,     0.0}, UnitCategory::Pressure}},
        {"Torr",    {{133.32236842105263, 0.0}, UnitCategory::Pressure}},
        {"torr",    {{133.32236842105263, 0.0}, UnitCategory::Pressure}},

        // ---- TEMPERATURE (base: K) ----
        {"K",       {{1.0,    0.0},                  UnitCategory::Temperature}},
        {"degK",    {{1.0,    0.0},                  UnitCategory::Temperature}},
        {"°C",      {{1.0,    273.15},               UnitCategory::Temperature}},
        {"degC",    {{1.0,    273.15},               UnitCategory::Temperature}},
        {"C",       {{1.0,    273.15},               UnitCategory::Temperature}}, // alias; ambiguous
        {"°F",      {{5.0/9.0, 255.37222222222223},  UnitCategory::Temperature}},
        {"degF",    {{5.0/9.0, 255.37222222222223},  UnitCategory::Temperature}},
        {"F",       {{5.0/9.0, 255.37222222222223},  UnitCategory::Temperature}},

        // ---- LENGTH (base: m) ----
        {"m",       {{1.0,    0.0}, UnitCategory::Length}},
        {"cm",      {{1e-2,   0.0}, UnitCategory::Length}},
        {"mm",      {{1e-3,   0.0}, UnitCategory::Length}},
        {"um",      {{1e-6,   0.0}, UnitCategory::Length}},
        {"µm",      {{1e-6,   0.0}, UnitCategory::Length}},
        {"nm",      {{1e-9,   0.0}, UnitCategory::Length}},
        {"A",       {{1e-10,  0.0}, UnitCategory::Length}},
        {"in",      {{0.0254, 0.0}, UnitCategory::Length}},
        {"ft",      {{0.3048, 0.0}, UnitCategory::Length}},

        // ---- TIME (base: s) ----
        {"s",       {{1.0,     0.0}, UnitCategory::Time}},
        {"ms",      {{1e-3,    0.0}, UnitCategory::Time}},
        {"µs",      {{1e-6,    0.0}, UnitCategory::Time}},
        {"us",      {{1e-6,    0.0}, UnitCategory::Time}},
        {"min",     {{60.0,    0.0}, UnitCategory::Time}},
        {"h",       {{3600.0,  0.0}, UnitCategory::Time}},
        {"hr",      {{3600.0,  0.0}, UnitCategory::Time}},
        {"d",       {{86400.0, 0.0}, UnitCategory::Time}},

        // ---- MASS (base: kg) ----
        {"kg",      {{1.0,        0.0}, UnitCategory::Mass}},
        {"g",       {{1e-3,       0.0}, UnitCategory::Mass}},
        {"mg",      {{1e-6,       0.0}, UnitCategory::Mass}},
        {"t",       {{1000.0,     0.0}, UnitCategory::Mass}},
        {"tonne",   {{1000.0,     0.0}, UnitCategory::Mass}},
        {"lb",      {{0.45359237, 0.0}, UnitCategory::Mass}},

        // ---- ENERGY (base: J) ----
        {"J",       {{1.0,    0.0}, UnitCategory::Energy}},
        {"kJ",      {{1e3,    0.0}, UnitCategory::Energy}},
        {"MJ",      {{1e6,    0.0}, UnitCategory::Energy}},
        {"cal",     {{4.184,  0.0}, UnitCategory::Energy}},
        {"kcal",    {{4184.0, 0.0}, UnitCategory::Energy}},
        {"Btu",     {{1055.06, 0.0}, UnitCategory::Energy}},

        // ---- SPECIFIC MOLAR ENERGY (base: J/mol) ----
        {"J/mol",       {{1.0,     0.0}, UnitCategory::SpecificMolarEnergy}},
        {"kJ/mol",      {{1e3,     0.0}, UnitCategory::SpecificMolarEnergy}},
        {"MJ/mol",      {{1e6,     0.0}, UnitCategory::SpecificMolarEnergy}},
        {"cal/mol",     {{4.184,   0.0}, UnitCategory::SpecificMolarEnergy}},
        {"kcal/mol",    {{4184.0,  0.0}, UnitCategory::SpecificMolarEnergy}},
        {"Btu/mol",     {{1055.06, 0.0}, UnitCategory::SpecificMolarEnergy}},

        // ---- POWER (no matching UnitCategory — usable only in Custom variables) ----
        {"W",       {{1.0,   0.0}, UnitCategory::Custom}},
        {"kW",      {{1e3,   0.0}, UnitCategory::Custom}},
        {"MW",      {{1e6,   0.0}, UnitCategory::Custom}},
        {"GW",      {{1e9,   0.0}, UnitCategory::Custom}},
        {"mW",      {{1e-3,  0.0}, UnitCategory::Custom}},

        // ---- MOLAR AMOUNT (base: mol) ----
        {"mol",     {{1.0,   0.0}, UnitCategory::MolarAmount}},
        {"mmol",    {{1e-3,  0.0}, UnitCategory::MolarAmount}},
        {"kmol",    {{1e3,   0.0}, UnitCategory::MolarAmount}},

        // ---- VOLUME (base: m³) ----
        {"m^3",     {{1.0,           0.0}, UnitCategory::Volume}},
        {"m³",      {{1.0,           0.0}, UnitCategory::Volume}},
        {"L",       {{1e-3,          0.0}, UnitCategory::Volume}},
        {"l",       {{1e-3,          0.0}, UnitCategory::Volume}},
        {"mL",      {{1e-6,          0.0}, UnitCategory::Volume}},
        {"ml",      {{1e-6,          0.0}, UnitCategory::Volume}},
        {"ft^3",    {{0.028316846592, 0.0}, UnitCategory::Volume}},
        {"ft³",     {{0.028316846592, 0.0}, UnitCategory::Volume}},

        // ---- VOLUMETRIC FLOW (base: m³/s) ----
        {"m^3/s",   {{1.0,                   0.0}, UnitCategory::VolumetricFlowRate}},
        {"m³/s",    {{1.0,                   0.0}, UnitCategory::VolumetricFlowRate}},
        {"L/s",     {{1e-3,                  0.0}, UnitCategory::VolumetricFlowRate}},
        {"cm^3/s",  {{1e-6,                  0.0}, UnitCategory::VolumetricFlowRate}},
        {"m^3/h",   {{1.0/3600.0,            0.0}, UnitCategory::VolumetricFlowRate}},
        {"m³/h",    {{1.0/3600.0,            0.0}, UnitCategory::VolumetricFlowRate}},
        {"L/min",   {{1e-3/60.0,             0.0}, UnitCategory::VolumetricFlowRate}},
        {"ccm",     {{1e-6/60.0,             0.0}, UnitCategory::VolumetricFlowRate}},
        {"mL/min",  {{1e-6/60.0,             0.0}, UnitCategory::VolumetricFlowRate}},
        {"ft^3/min",{{0.028316846592/60.0,   0.0}, UnitCategory::VolumetricFlowRate}},

        // ---- MASS FLOW (base: kg/s) ----
        {"kg/s",    {{1.0,           0.0}, UnitCategory::MassFlowRate}},
        {"g/s",     {{1e-3,          0.0}, UnitCategory::MassFlowRate}},
        {"mg/s",    {{1e-6,          0.0}, UnitCategory::MassFlowRate}},
        {"kg/min",  {{1.0/60.0,      0.0}, UnitCategory::MassFlowRate}},
        {"kg/h",    {{1.0/3600.0,    0.0}, UnitCategory::MassFlowRate}},
        {"g/min",   {{1e-3/60.0,     0.0}, UnitCategory::MassFlowRate}},
        {"g/h",     {{1e-3/3600.0,   0.0}, UnitCategory::MassFlowRate}},
        {"t/h",     {{1000.0/3600.0, 0.0}, UnitCategory::MassFlowRate}},

        // ---- MOLAR FLOW (base: mol/s) ----
        // Normal/standard gas units assume 0°C (273.15 K), 1 atm (101325 Pa):
        // molar volume V_m = 22.414e-3 m³/mol, so 1 Nm³ = 1/22.414e-3 ≈ 44.6149 mol.
        {"mol/s",   {{1.0,           0.0}, UnitCategory::MolarFlowRate}},
        {"mol/min", {{1.0/60.0,      0.0}, UnitCategory::MolarFlowRate}},
        {"mol/h",   {{1.0/3600.0,    0.0}, UnitCategory::MolarFlowRate}},
        {"mmol/s",  {{1e-3,          0.0}, UnitCategory::MolarFlowRate}},
        {"umol/s",  {{1e-6,          0.0}, UnitCategory::MolarFlowRate}},
        {"kmol/h",  {{1000.0/3600.0, 0.0}, UnitCategory::MolarFlowRate}},
        {"kmol/s",  {{1000.0,        0.0}, UnitCategory::MolarFlowRate}},
        {"Nm^3/s",  {{1.0/22.414e-3,            0.0}, UnitCategory::MolarFlowRate}},
        {"Nm³/s",   {{1.0/22.414e-3,            0.0}, UnitCategory::MolarFlowRate}},
        {"Nm^3/h",  {{1.0/22.414e-3/3600.0,     0.0}, UnitCategory::MolarFlowRate}},
        {"Nm³/h",   {{1.0/22.414e-3/3600.0,     0.0}, UnitCategory::MolarFlowRate}},
        {"SCCM",    {{1e-6/22.414e-3/60.0,      0.0}, UnitCategory::MolarFlowRate}},
        {"SCFM",    {{0.028316846592/22.414e-3/60.0, 0.0}, UnitCategory::MolarFlowRate}},

        // ---- VELOCITY (base: m/s) ----
        {"m/s",     {{1.0,           0.0}, UnitCategory::Velocity}},
        {"cm/s",    {{1e-2,          0.0}, UnitCategory::Velocity}},
        {"mm/s",    {{1e-3,          0.0}, UnitCategory::Velocity}},
        {"m/min",   {{1.0/60.0,      0.0}, UnitCategory::Velocity}},
        {"km/h",    {{1000.0/3600.0, 0.0}, UnitCategory::Velocity}},
        {"ft/s",    {{0.3048,        0.0}, UnitCategory::Velocity}},
        {"mph",     {{0.44704,       0.0}, UnitCategory::Velocity}},

        // ---- INVERSE TIME (base: 1/s) ----
        {"1/s",     {{1.0,       0.0}, UnitCategory::InverseTime}},
        {"1/min",   {{1.0/60.0,  0.0}, UnitCategory::InverseTime}},
        {"1/h",     {{1.0/3600.0, 0.0}, UnitCategory::InverseTime}},

        // ---- DENSITY (base: kg/m³) ----
        {"kg/m^3",  {{1.0,       0.0}, UnitCategory::Density}},
        {"kg/m³",   {{1.0,       0.0}, UnitCategory::Density}},
        {"g/L",     {{1.0,       0.0}, UnitCategory::Density}},
        {"kg/L",    {{1000.0,    0.0}, UnitCategory::Density}},
        {"g/mL",    {{1000.0,    0.0}, UnitCategory::Density}},
        {"g/cm^3",  {{1000.0,    0.0}, UnitCategory::Density}},
        {"g/cm³",   {{1000.0,    0.0}, UnitCategory::Density}},
        {"lb/ft^3", {{16.018463, 0.0}, UnitCategory::Density}},

        // ---- MOLAR CONCENTRATION (base: mol/m³) ----
        {"mol/m^3", {{1.0,    0.0}, UnitCategory::MolarConcentration}},
        {"mol/m³",  {{1.0,    0.0}, UnitCategory::MolarConcentration}},
        {"mol/L",   {{1000.0, 0.0}, UnitCategory::MolarConcentration}},
        {"mmol/L",  {{1.0,    0.0}, UnitCategory::MolarConcentration}},
        {"umol/L",  {{1e-3,   0.0}, UnitCategory::MolarConcentration}},
        {"kmol/m^3",{{1000.0, 0.0}, UnitCategory::MolarConcentration}},
        {"kmol/m³", {{1000.0, 0.0}, UnitCategory::MolarConcentration}},

        // ---- ADSORBED CONCENTRATION (base: mol/kg) ----
        {"mol/kg",  {{1.0,    0.0}, UnitCategory::AdsorbedConcentration}},
        {"mmol/g",  {{1.0,    0.0}, UnitCategory::AdsorbedConcentration}},
        {"mol/g",   {{1000.0, 0.0}, UnitCategory::AdsorbedConcentration}},
        {"mmol/kg", {{1e-3,   0.0}, UnitCategory::AdsorbedConcentration}},
        {"umol/g",  {{1e-3,   0.0}, UnitCategory::AdsorbedConcentration}},

        // ---- PRESSURE INVERSE (no matching UnitCategory — Custom only) ----
        {"1/Pa",    {{1.0,                       0.0}, UnitCategory::Custom}},
        {"1/kPa",   {{1.0/1e3,                   0.0}, UnitCategory::Custom}},
        {"1/MPa",   {{1.0/1e6,                   0.0}, UnitCategory::Custom}},
        {"1/bar",   {{1.0/1e5,                   0.0}, UnitCategory::Custom}},
        {"1/atm",   {{1.0/101325.0,              0.0}, UnitCategory::Custom}},
        {"1/psi",   {{1.0/6894.757293168,        0.0}, UnitCategory::Custom}},
        {"1/Torr",  {{1.0/133.32236842105263,    0.0}, UnitCategory::Custom}},

        // ---- HEAT CAPACITY (base: J/(kg·K)) ----
        {"J/(kg·K)",    {{1.0,    0.0}, UnitCategory::HeatCapacity}},
        {"kJ/(kg·K)",   {{1e3,    0.0}, UnitCategory::HeatCapacity}},
        {"cal/(g·K)",   {{4184.0, 0.0}, UnitCategory::HeatCapacity}},

        // ---- THERMAL CONDUCTIVITY (base: W/(m·K)) ----
        {"W/(m·K)",     {{1.0,     0.0}, UnitCategory::ThermalConductivity}},
        {"mW/(m·K)",    {{1e-3,    0.0}, UnitCategory::ThermalConductivity}},
        {"W/(cm·K)",    {{100.0,   0.0}, UnitCategory::ThermalConductivity}},

        // ---- MOLAR MASS (base: kg/mol) ----
        {"kg/mol",  {{1.0,   0.0}, UnitCategory::MolarMass}},
        {"g/mol",   {{1e-3,  0.0}, UnitCategory::MolarMass}},
        {"kg/kmol", {{1e-3,  0.0}, UnitCategory::MolarMass}},

        // ---- SPECIFIC ENERGY (base: J/kg) ----
        {"J/kg",        {{1.0,   0.0}, UnitCategory::SpecificEnergy}},
        {"kJ/kg",       {{1e3,   0.0}, UnitCategory::SpecificEnergy}},
        {"MJ/kg",       {{1e6,   0.0}, UnitCategory::SpecificEnergy}},
        {"kWh/kg",      {{3.6e6, 0.0}, UnitCategory::SpecificEnergy}},
        {"kWh/tonne",   {{3.6e3, 0.0}, UnitCategory::SpecificEnergy}},

        // ---- DIFFUSIVITY (base: m²/s) ----
        {"m^2/s",   {{1.0,   0.0}, UnitCategory::Diffusivity}},
        {"m²/s",    {{1.0,   0.0}, UnitCategory::Diffusivity}},
        {"cm^2/s",  {{1e-4,  0.0}, UnitCategory::Diffusivity}},
        {"cm²/s",   {{1e-4,  0.0}, UnitCategory::Diffusivity}},
        {"mm^2/s",  {{1e-6,  0.0}, UnitCategory::Diffusivity}},
        {"mm²/s",   {{1e-6,  0.0}, UnitCategory::Diffusivity}},

        // ---- AREA (base: m²) ----
        {"m^2",     {{1.0,   0.0}, UnitCategory::Area}},
        {"m²",      {{1.0,   0.0}, UnitCategory::Area}},
        {"cm^2",    {{1e-4,  0.0}, UnitCategory::Area}},
        {"cm²",     {{1e-4,  0.0}, UnitCategory::Area}},
        {"mm^2",    {{1e-6,  0.0}, UnitCategory::Area}},
        {"mm²",     {{1e-6,  0.0}, UnitCategory::Area}},

        // ---- PRESSURE GRADIENT (base: Pa/m) ----
        {"Pa/m",    {{1.0,      0.0}, UnitCategory::PressureGradient}},
        {"kPa/m",   {{1e3,      0.0}, UnitCategory::PressureGradient}},
        {"MPa/m",   {{1e6,      0.0}, UnitCategory::PressureGradient}},
        {"bar/m",   {{1e5,      0.0}, UnitCategory::PressureGradient}},
        {"atm/m",   {{101325.0, 0.0}, UnitCategory::PressureGradient}},

        // ---- DIMENSIONLESS / MOLE FRACTION ----
        // % and ppm are tagged as MoleFraction since that's their main use here;
        // bare "-" stays Dimensionless. Empty key is unreachable via the regex
        // but kept for direct lookup parity.
        {"-",       {{1.0, 0.0},   UnitCategory::Dimensionless}},
        {"",        {{1.0, 0.0},   UnitCategory::Dimensionless}},
        {"mol/mol", {{1.0, 0.0},   UnitCategory::MoleFraction}},
        {"%",       {{0.01, 0.0},  UnitCategory::MoleFraction}},
        {"ppm",     {{1e-6, 0.0},  UnitCategory::MoleFraction}},
    };
    return kMap;
}

bool ParseUnitSuffix(const std::string& expr,
                     std::string&       math_out,
                     UnitXform&         xform_out,
                     UnitCategory*      category_out) {
    // Match optional whitespace, then [symbol], at end of string.
    static const std::regex kUnitRegex(R"(\s*\[([^\]]+)\]\s*$)");

    std::smatch m;
    if (!std::regex_search(expr, m, kUnitRegex)) {
        return false;
    }

    const std::string symbol = m[1].str();
    const auto& map = UnitMap();
    auto it = map.find(symbol);
    if (it == map.end()) {
        return false;  // unknown symbol — do not apply conversion
    }

    // math_out is everything before the [symbol] bracket
    math_out  = expr.substr(0, static_cast<size_t>(m.position(0)));
    xform_out = it->second.xform;
    if (category_out) *category_out = it->second.category;
    return true;
}

const char* LookupUnitSymbol(const std::string& key) {
    const auto& map = UnitMap();
    auto it = map.find(key);
    // Return a pointer into the static map's stored key (program lifetime),
    // not into the caller-supplied argument (which may be a temporary).
    return (it != map.end()) ? it->first.c_str() : nullptr;
}

} // namespace adsim
