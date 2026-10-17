#include "adsim/math/unit_parser.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <regex>
#include <stdexcept>
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

// ---------------------------------------------------------------------------
// Unit expressions (ParseUnitExpression)
// ---------------------------------------------------------------------------

namespace {

constexpr Dimension kNone{};
constexpr Dimension kLength{1, 0, 0, 0, 0};
constexpr Dimension kMass{0, 1, 0, 0, 0};
constexpr Dimension kTime{0, 0, 1, 0, 0};
constexpr Dimension kTemperature{0, 0, 0, 1, 0};
constexpr Dimension kAmount{0, 0, 0, 0, 1};
constexpr Dimension kVolume{3, 0, 0, 0, 0};
constexpr Dimension kForce{1, 1, -2, 0, 0};
constexpr Dimension kPressure{-1, 1, -2, 0, 0};
constexpr Dimension kEnergy{2, 1, -2, 0, 0};
constexpr Dimension kPower{2, 1, -3, 0, 0};
constexpr Dimension kVelocity{1, 0, -1, 0, 0};
constexpr Dimension kInverseTime{0, 0, -1, 0, 0};
constexpr Dimension kVolumetricFlow{3, 0, -1, 0, 0};
constexpr Dimension kMolarFlow{0, 0, -1, 0, 1};

Dimension Combine(const Dimension& a, const Dimension& b, int power) {
    return {a.length + power * b.length, a.mass + power * b.mass, a.time + power * b.time,
            a.temperature + power * b.temperature, a.amount + power * b.amount};
}

struct Atom {
    double factor = 1.0;
    double offset = 0.0;
    Dimension dimension;
};

// Units written as a whole word, never with a prefix.
const std::unordered_map<std::string, Atom>& NamedAtoms() {
    static const std::unordered_map<std::string, Atom> kAtoms = {
        {"min",   {60.0, 0.0, kTime}},
        {"h",     {3600.0, 0.0, kTime}},
        {"hr",    {3600.0, 0.0, kTime}},
        {"d",     {86400.0, 0.0, kTime}},
        {"day",   {86400.0, 0.0, kTime}},
        {"atm",   {101325.0, 0.0, kPressure}},
        {"psi",   {6894.757293168, 0.0, kPressure}},
        {"Torr",  {133.32236842105263, 0.0, kPressure}},
        {"torr",  {133.32236842105263, 0.0, kPressure}},
        {"degC",  {1.0, 273.15, kTemperature}},
        {"\xC2\xB0" "C", {1.0, 273.15, kTemperature}},
        {"C",     {1.0, 273.15, kTemperature}},
        {"degF",  {5.0 / 9.0, 255.37222222222223, kTemperature}},
        {"\xC2\xB0" "F", {5.0 / 9.0, 255.37222222222223, kTemperature}},
        {"F",     {5.0 / 9.0, 255.37222222222223, kTemperature}},
        {"degK",  {1.0, 0.0, kTemperature}},
        {"in",    {0.0254, 0.0, kLength}},
        {"ft",    {0.3048, 0.0, kLength}},
        {"A",     {1e-10, 0.0, kLength}},
        {"\xC3\x85", {1e-10, 0.0, kLength}},
        {"lb",    {0.45359237, 0.0, kMass}},
        {"t",     {1000.0, 0.0, kMass}},
        {"tonne", {1000.0, 0.0, kMass}},
        {"Btu",   {1055.06, 0.0, kEnergy}},
        {"Wh",    {3600.0, 0.0, kEnergy}},
        {"kWh",   {3.6e6, 0.0, kEnergy}},
        {"mph",   {0.44704, 0.0, kVelocity}},
        {"ccm",   {1e-6 / 60.0, 0.0, kVolumetricFlow}},
        // Standard gas flows: 0 °C and 1 atm, V_m = 22.414e-3 m³/mol (as the dictionary).
        {"SCCM",  {1e-6 / 22.414e-3 / 60.0, 0.0, kMolarFlow}},
        {"SCFM",  {0.028316846592 / 22.414e-3 / 60.0, 0.0, kMolarFlow}},
        {"%",     {0.01, 0.0, kNone}},
        {"ppm",   {1e-6, 0.0, kNone}},
    };
    return kAtoms;
}

// SI (and SI-style) units that take a prefix: kPa, mmol, cm, kJ, mbar, mL...
const std::unordered_map<std::string, Atom>& PrefixableAtoms() {
    static const std::unordered_map<std::string, Atom> kAtoms = {
        {"m",   {1.0, 0.0, kLength}},
        {"g",   {1e-3, 0.0, kMass}},
        {"s",   {1.0, 0.0, kTime}},
        {"K",   {1.0, 0.0, kTemperature}},
        {"mol", {1.0, 0.0, kAmount}},
        {"Pa",  {1.0, 0.0, kPressure}},
        {"bar", {1e5, 0.0, kPressure}},
        {"J",   {1.0, 0.0, kEnergy}},
        {"cal", {4.184, 0.0, kEnergy}},
        {"W",   {1.0, 0.0, kPower}},
        {"N",   {1.0, 0.0, kForce}},
        {"L",   {1e-3, 0.0, kVolume}},
        {"l",   {1e-3, 0.0, kVolume}},
        {"Hz",  {1.0, 0.0, kInverseTime}},
    };
    return kAtoms;
}

bool LookupAtom(const std::string& name, Atom& out) {
    if (auto it = NamedAtoms().find(name); it != NamedAtoms().end()) {
        out = it->second;
        return true;
    }
    const auto& prefixable = PrefixableAtoms();
    if (auto it = prefixable.find(name); it != prefixable.end()) {
        out = it->second;
        return true;
    }
    static const std::pair<const char*, double> kPrefixes[] = {
        {"G", 1e9}, {"M", 1e6}, {"k", 1e3}, {"h", 1e2}, {"d", 1e-1}, {"c", 1e-2},
        {"m", 1e-3}, {"u", 1e-6}, {"\xC2\xB5", 1e-6}, {"\xCE\xBC", 1e-6}, {"n", 1e-9},
    };
    for (const auto& [prefix, scale] : kPrefixes) {
        const std::string p = prefix;
        if (name.size() <= p.size() || name.compare(0, p.size(), p) != 0) continue;
        if (auto it = prefixable.find(name.substr(p.size())); it != prefixable.end()) {
            out = it->second;
            out.factor *= scale;
            return true;
        }
    }
    return false;
}

std::string Trim(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return {};
    const size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

class UnitExpressionParser {
public:
    explicit UnitExpressionParser(const std::string& text) : s_(text) {}

    bool Parse(UnitInfo& out, std::string& error) {
        Atom result;
        if (!Quotient(result)) {
            error = error_;
            return false;
        }
        SkipSpaces();
        if (i_ != s_.size()) {
            error = "unexpected '" + s_.substr(i_) + "'";
            return false;
        }
        out.factor = result.factor;
        out.dimension = result.dimension;
        // An offset (°C, °F) is absolute only when the unit is that atom alone.
        out.offset = (atoms_ == 1 && lone_offset_) ? offset_ : 0.0;
        return true;
    }

private:
    const std::string& s_;
    size_t i_ = 0;
    std::string error_;
    int atoms_ = 0;
    bool lone_offset_ = false;
    double offset_ = 0.0;

    void SkipSpaces() {
        while (i_ < s_.size() && s_[i_] == ' ') ++i_;
    }

    bool Has(const char* token, size_t at) const { return s_.compare(at, std::char_traits<char>::length(token), token) == 0; }

    // Length of the atom character at `at`: an ASCII letter, '%', or one of
    // the two-byte letters ° µ Å μ. 0 for anything else (², ³ and · included).
    size_t AtomCharLength(size_t at) const {
        if (at >= s_.size()) return 0;
        const unsigned char c = static_cast<unsigned char>(s_[at]);
        if (std::isalpha(c) || c == '%') return 1;
        if (Has("\xC2\xB0", at) || Has("\xC2\xB5", at) || Has("\xC3\x85", at) || Has("\xCE\xBC", at)) return 2;
        return 0;
    }

    bool StartsUnit(size_t at) const { return AtomCharLength(at) > 0 || (at < s_.size() && s_[at] == '('); }

    bool Fail(const std::string& why) {
        if (error_.empty()) error_ = why;
        return false;
    }

    // quotient := product ('/' product)*
    bool Quotient(Atom& out) {
        if (!Product(out)) return false;
        for (;;) {
            SkipSpaces();
            if (i_ >= s_.size() || s_[i_] != '/') return true;
            ++i_;
            Atom denominator;
            if (!Product(denominator)) return false;
            out.factor /= denominator.factor;
            out.dimension = Combine(out.dimension, denominator.dimension, -1);
        }
    }

    // product := factor (separator factor)*, separator = '*' '.' '·' '-' or spaces.
    // Only the first factor of a product may be a number ("1/s"), so "s-1" is
    // rejected rather than read as s·1.
    bool Product(Atom& out) {
        SkipSpaces();
        if (!Factor(out, /*allow_number=*/true)) return false;
        for (;;) {
            const size_t before = i_;
            SkipSpaces();
            bool explicit_separator = false;
            if (i_ < s_.size() && (s_[i_] == '*' || s_[i_] == '.' || s_[i_] == '-')) {
                ++i_;
                explicit_separator = true;
            } else if (Has("\xC2\xB7", i_)) {
                i_ += 2;
                explicit_separator = true;
            }
            if (explicit_separator) {
                SkipSpaces();
                if (!StartsUnit(i_))
                    return Fail("expected a unit after '" + Trim(s_.substr(before, i_ - before)) + "' in '" + s_ + "'");
            } else if (i_ == before || !StartsUnit(i_)) {
                i_ = before;  // the product ends here ('/', ')' or the end)
                return true;
            }
            Atom next;
            if (!Factor(next, /*allow_number=*/false)) return false;
            out.factor *= next.factor;
            out.dimension = Combine(out.dimension, next.dimension, 1);
        }
    }

    // exponent := '^' int | '**' int | digits | '²' | '³'   (optional)
    bool Exponent(int& power) {
        power = 1;
        size_t at = i_;
        if (Has("**", at)) at += 2;
        else if (at < s_.size() && s_[at] == '^') at += 1;
        else if (Has("\xC2\xB2", at)) { i_ += 2; power = 2; return true; }
        else if (Has("\xC2\xB3", at)) { i_ += 2; power = 3; return true; }
        else if (at >= s_.size() || !std::isdigit(static_cast<unsigned char>(s_[at]))) return true;
        const bool explicit_marker = at != i_;
        int sign = 1;
        if (explicit_marker && at < s_.size() && (s_[at] == '-' || s_[at] == '+')) {
            if (s_[at] == '-') sign = -1;
            ++at;
        }
        if (at >= s_.size() || !std::isdigit(static_cast<unsigned char>(s_[at])))
            return Fail("expected an integer exponent in '" + s_ + "'");
        int value = 0;
        while (at < s_.size() && std::isdigit(static_cast<unsigned char>(s_[at]))) {
            value = value * 10 + (s_[at] - '0');
            if (value > 99) return Fail("exponent too large in '" + s_ + "'");
            ++at;
        }
        i_ = at;
        power = sign * value;
        return true;
    }

    // factor := (atom | number | '(' quotient ')') exponent?
    bool Factor(Atom& out, bool allow_number) {
        if (i_ >= s_.size()) return Fail("missing unit at the end of '" + s_ + "'");
        const unsigned char c = static_cast<unsigned char>(s_[i_]);
        if (c == '(') {
            ++i_;
            if (!Quotient(out)) return false;
            SkipSpaces();
            if (i_ >= s_.size() || s_[i_] != ')') return Fail("missing ')' in '" + s_ + "'");
            ++i_;
            atoms_ += 2;  // a group is never a lone offset unit
        } else if (allow_number && (std::isdigit(c) || c == '.')) {
            const char* begin = s_.c_str() + i_;
            char* end = nullptr;
            const double value = std::strtod(begin, &end);
            if (end == begin || !(value > 0.0) || !std::isfinite(value))
                return Fail("bad number in '" + s_ + "'");
            i_ += static_cast<size_t>(end - begin);
            out = Atom{value, 0.0, kNone};
            atoms_ += 2;  // "1 degC" is not an absolute temperature unit either
            return true;  // no exponent on a number
        } else {
            const size_t begin = i_;
            while (size_t n = AtomCharLength(i_)) i_ += n;
            if (i_ == begin) return Fail("unexpected '" + s_.substr(i_) + "' in '" + s_ + "'");
            const std::string name = s_.substr(begin, i_ - begin);
            int power = 1;
            if (!Exponent(power)) return false;
            if (name == "Nm" && power == 3) {
                // Normal cubic metre: the gas amount in 1 m³ at 0 °C, 1 atm.
                out = Atom{1.0 / 22.414e-3, 0.0, kAmount};
                ++atoms_;
                return true;
            }
            Atom atom;
            if (!LookupAtom(name, atom)) return Fail("unknown unit '" + name + "'");
            ++atoms_;
            lone_offset_ = atom.offset != 0.0 && power == 1;
            offset_ = atom.offset;
            out.factor = std::pow(atom.factor, power);
            out.dimension = Combine(kNone, atom.dimension, power);
            return true;
        }
        int power = 1;
        if (!Exponent(power)) return false;
        out.factor = std::pow(out.factor, power);
        out.dimension = Combine(kNone, out.dimension, power);
        return true;
    }
};


// The UnitCategory of a dimension, for units found by ParseUnitExpression
// rather than the dictionary: the lowest-numbered non-Custom category with a
// dictionary symbol of that dimension (Dimensionless for no dimension).
UnitCategory CategoryOf(const Dimension& dimension) {
    if (dimension.IsDimensionless()) return UnitCategory::Dimensionless;
    static const std::vector<std::pair<Dimension, UnitCategory>> kCategories = [] {
        std::vector<std::pair<Dimension, UnitCategory>> v;
        for (const auto& [symbol, entry] : UnitMap()) {
            UnitInfo info;
            if (entry.category == UnitCategory::Custom || !ParseUnitExpression(symbol, info)) continue;
            v.emplace_back(info.dimension, entry.category);
        }
        std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) {
            return static_cast<int>(a.second) < static_cast<int>(b.second);
        });
        return v;
    }();
    for (const auto& [dim, category] : kCategories)
        if (dim == dimension) return category;
    return UnitCategory::Custom;
}

} // namespace

bool ParseUnitExpression(const std::string& text, UnitInfo& out, std::string* error) {
    const std::string t = Trim(text);
    if (t.empty() || t == "-") {
        out = UnitInfo{};
        return true;
    }
    UnitExpressionParser parser(t);
    std::string why;
    UnitInfo parsed;
    if (!parser.Parse(parsed, why)) {
        if (error) *error = why;
        return false;
    }
    out = parsed;
    return true;
}

std::string DimensionName(const Dimension& d) {
    static const std::vector<std::pair<Dimension, const char*>> kNames = {
        {kNone, "dimensionless"},
        {kTime, "time"},
        {kLength, "length"},
        {{2, 0, 0, 0, 0}, "area"},
        {kVolume, "volume"},
        {kMass, "mass"},
        {kTemperature, "temperature"},
        {kAmount, "amount of substance"},
        {kPressure, "pressure"},
        {kEnergy, "energy"},
        {kPower, "power"},
        {kForce, "force"},
        {kVelocity, "velocity"},
        {kInverseTime, "inverse time"},
        {kVolumetricFlow, "volumetric flow"},
        {kMolarFlow, "molar flow"},
        {{0, 1, -1, 0, 0}, "mass flow"},
        {{-3, 1, 0, 0, 0}, "density"},
        {{-3, 0, 0, 0, 1}, "molar concentration"},
        {{0, -1, 0, 0, 1}, "adsorbed amount per mass"},
        {{2, 1, -2, 0, -1}, "molar energy"},
        {{2, 0, -2, 0, 0}, "specific energy"},
        {{2, 0, -2, -1, 0}, "specific heat capacity"},
        {{1, 1, -3, -1, 0}, "thermal conductivity"},
        {{0, 1, -3, -1, 0}, "heat transfer coefficient"},
        {{2, 0, -1, 0, 0}, "diffusivity"},
        {{1, -1, 2, 0, 0}, "inverse pressure"},
        {{3, 0, 0, 0, -1}, "molar volume"},
        {{-1, 1, -1, 0, 0}, "dynamic viscosity"},
        {{0, 1, 0, 0, -1}, "molar mass"},
        {{0, 0, 0, -1, 0}, "inverse temperature"},
        {{-2, 1, -2, 0, 0}, "pressure gradient"},
    };
    for (const auto& [dim, name] : kNames)
        if (dim == d) return name;
    std::string out;
    const std::pair<const char*, int> parts[] = {
        {"kg", d.mass}, {"m", d.length}, {"s", d.time}, {"K", d.temperature}, {"mol", d.amount}};
    for (const auto& [symbol, power] : parts) {
        if (power == 0) continue;
        if (!out.empty()) out += "\xC2\xB7";
        out += symbol;
        if (power != 1) out += "^" + std::to_string(power);
    }
    return out;
}

double ConvertUnit(double value, const std::string& from, const std::string& to) {
    UnitInfo a, b;
    std::string why;
    if (!ParseUnitExpression(from, a, &why))
        throw std::invalid_argument("unknown unit '" + from + "': " + why);
    if (!ParseUnitExpression(to, b, &why))
        throw std::invalid_argument("unknown unit '" + to + "': " + why);
    if (!(a.dimension == b.dimension))
        throw std::invalid_argument("cannot convert '" + from + "' (" + DimensionName(a.dimension) + ") to '" + to +
                                    "' (" + DimensionName(b.dimension) + ")");
    return (value * a.factor + a.offset - b.offset) / b.factor;
}

bool SplitQuantity(const std::string& text, double& number, std::string& unit) {
    const std::string t = Trim(text);
    if (t.empty()) return false;
    const unsigned char first = static_cast<unsigned char>(t[0]);
    if (!(std::isdigit(first) || first == '.' || first == '-' || first == '+')) return false;
    // Hexadecimal, inf and nan are not quantities.
    for (char c : t.substr(0, std::min<size_t>(t.size(), 3)))
        if (c == 'x' || c == 'X') return false;
    const char* begin = t.c_str();
    char* end = nullptr;
    const double value = std::strtod(begin, &end);
    if (end == begin || !std::isfinite(value)) return false;
    std::string rest = Trim(std::string(end));
    if (rest.size() >= 2 && rest.front() == '[' && rest.back() == ']') rest = Trim(rest.substr(1, rest.size() - 2));
    if (rest.empty()) return false;
    // A unit never starts with an operator: "2*x" and "2 + 3" are formulas.
    if (rest[0] == '*' || rest[0] == '+' || rest[0] == '-' || rest[0] == '/' || rest[0] == '^' || rest[0] == '(')
        return false;
    UnitInfo info;
    if (!ParseUnitExpression(rest, info)) return false;
    number = value;
    unit = rest;
    return true;
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
    UnitXform xform;
    UnitCategory category = UnitCategory::Custom;
    if (it != map.end()) {
        xform = it->second.xform;
        category = it->second.category;
    } else {
        UnitInfo info;
        if (!ParseUnitExpression(symbol, info)) {
            return false;  // unknown symbol — do not apply conversion
        }
        xform = {info.factor, info.offset};
        category = CategoryOf(info.dimension);
    }

    // math_out is everything before the [symbol] bracket
    math_out  = expr.substr(0, static_cast<size_t>(m.position(0)));
    xform_out = xform;
    if (category_out) *category_out = category;
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
