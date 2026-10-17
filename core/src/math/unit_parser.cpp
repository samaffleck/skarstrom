#include "adsim/math/unit_parser.h"

#include <units/units.hpp>

namespace llnl = UNITS_NAMESPACE;  // "units": LLNL units, compiled into adsim

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <regex>
#include <stdexcept>
#include <unordered_map>

namespace adsim {

// ---------------------------------------------------------------------------
// Unit expressions (ParseUnitExpression)
//
// The syntax is ours: the separators, exponent forms and the engineering
// reading of '/' (everything after it, up to the next '/', is in the
// denominator) are what decks and the front ends have always written. Each
// unit symbol, its prefixes and the arithmetic come from the LLNL units
// library (https://github.com/LLNL/units), parsed strictly so that a typo is
// an error rather than some obscure unit. A few symbols mean something else
// in process engineering than in LLNL's default reading and are overridden
// (Overrides below): C, F and A are degrees Celsius / Fahrenheit and the
// ångström, not coulomb, farad and ampere; ccm, SCCM, SCFM and Nm³ are the
// usual gas-flow units, with standard conditions of 0 °C and 1 atm
// (V_m = 22.414e-3 m³/mol).
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

constexpr double kStandardMolarVolume = 22.414e-3;  // m³/mol at 0 °C, 1 atm

// Strict parsing of a single symbol: no splitting of unknown words into
// smaller units, no commodities ("kg{CO2}"), locality ("US gal") or
// "per"/"of" phrases.
constexpr std::uint64_t kStrict = llnl::skip_partition_check | llnl::no_commodities |
                                  llnl::no_locality_modifiers | llnl::no_of_operator |
                                  llnl::no_per_operators | llnl::no_addition;

const std::unordered_map<std::string, llnl::precise_unit>& Overrides() {
    using namespace units;
    static const std::unordered_map<std::string, precise_unit> kOverrides = {
        {"C",    precise::degC},
        {"F",    precise::degF},
        {"A",    precise::distance::angstrom},
        {"ccm",  precise::cm.pow(3) / precise::minute},
        {"SCCM", precise_unit(1e-6 / kStandardMolarVolume / 60.0, precise::mol / precise::s)},
        {"sccm", precise_unit(1e-6 / kStandardMolarVolume / 60.0, precise::mol / precise::s)},
        {"SCFM", precise_unit(0.028316846592 / kStandardMolarVolume / 60.0, precise::mol / precise::s)},
        {"ppm",  precise_unit(1e-6, precise::one)},
    };
    return kOverrides;
}

// The engine's five base dimensions; false for a unit outside them (current,
// luminosity, currency, angles, equation units such as dB or pH).
bool ToDimension(const llnl::precise_unit& u, Dimension& out) {
    const auto b = u.base_units();
    if (b.ampere() != 0 || b.candela() != 0 || b.currency() != 0 || b.radian() != 0 || b.is_equation() ||
        b.is_per_unit())
        return false;
    out = {b.meter(), b.kg(), b.second(), b.kelvin(), b.mole()};
    return true;
}

struct Atom {
    llnl::precise_unit unit = llnl::precise::one;
    bool offset_unit = false;  // degC / degF: has an offset when it stands alone
};

bool LookupAtom(const std::string& name, Atom& out) {
    if (auto it = Overrides().find(name); it != Overrides().end()) {
        out.unit = it->second;
    } else {
        out.unit = llnl::unit_from_string(name, kStrict);
        if (!llnl::is_valid(out.unit)) return false;
    }
    Dimension d;
    if (!ToDimension(out.unit, d)) return false;
    out.offset_unit = out.unit.base_units().has_e_flag() && d == kTemperature;
    return true;
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
        llnl::precise_unit result;
        if (!Quotient(result)) {
            error = error_;
            return false;
        }
        SkipSpaces();
        if (i_ != s_.size()) {
            error = "unexpected '" + s_.substr(i_) + "'";
            return false;
        }
        if (!ToDimension(result, out.dimension)) {
            error = "'" + s_ + "' is not a unit of length, mass, time, temperature or amount of substance";
            return false;
        }
        // A temperature with an offset (°C, °F) is absolute only when the unit
        // is that atom alone; inside a compound unit it is a difference, which
        // the multiplier alone describes.
        if (atoms_ == 1 && lone_offset_) {
            out.offset = llnl::convert(0.0, result, llnl::precise::K);
            out.factor = llnl::convert(1.0, result, llnl::precise::K) - out.offset;
        } else {
            out.offset = 0.0;
            out.factor = result.multiplier();
        }
        return true;
    }

private:
    const std::string& s_;
    size_t i_ = 0;
    std::string error_;
    int atoms_ = 0;
    bool lone_offset_ = false;

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
    bool Quotient(llnl::precise_unit& out) {
        if (!Product(out)) return false;
        for (;;) {
            SkipSpaces();
            if (i_ >= s_.size() || s_[i_] != '/') return true;
            ++i_;
            llnl::precise_unit denominator;
            if (!Product(denominator)) return false;
            out = out / denominator;
        }
    }

    // product := factor (separator factor)*, separator = '*' '.' '·' '-' or spaces.
    // Only the first factor of a product may be a number ("1/s"), so "s-1" is
    // rejected rather than read as s·1.
    bool Product(llnl::precise_unit& out) {
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
            llnl::precise_unit next;
            if (!Factor(next, /*allow_number=*/false)) return false;
            out = out * next;
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
    bool Factor(llnl::precise_unit& out, bool allow_number) {
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
            out = llnl::precise_unit(value, llnl::precise::one);
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
                out = llnl::precise_unit(1.0 / kStandardMolarVolume, llnl::precise::mol);
                ++atoms_;
                return true;
            }
            Atom atom;
            if (!LookupAtom(name, atom)) return Fail("unknown unit '" + name + "'");
            ++atoms_;
            lone_offset_ = atom.offset_unit && power == 1;
            out = power == 1 ? atom.unit : atom.unit.pow(power);
            return true;
        }
        int power = 1;
        if (!Exponent(power)) return false;
        if (power != 1) out = out.pow(power);
        return true;
    }
};


// The UnitCategory of a unit: by its dimension, with %, ppm and mol/mol as
// mole fractions. Dimensions no category describes (power, inverse
// pressure...) are Custom.
UnitCategory CategoryOf(const std::string& symbol, const Dimension& dimension) {
    if (dimension.IsDimensionless()) {
        const std::string t = Trim(symbol);
        return (t == "%" || t == "ppm" || t == "mol/mol") ? UnitCategory::MoleFraction : UnitCategory::Dimensionless;
    }
    static const std::pair<Dimension, UnitCategory> kCategories[] = {
        {kTime, UnitCategory::Time},
        {kLength, UnitCategory::Length},
        {kTemperature, UnitCategory::Temperature},
        {kPressure, UnitCategory::Pressure},
        {kEnergy, UnitCategory::Energy},
        {{2, 1, -2, 0, -1}, UnitCategory::SpecificMolarEnergy},
        {kInverseTime, UnitCategory::InverseTime},
        {kVelocity, UnitCategory::Velocity},
        {kVolumetricFlow, UnitCategory::VolumetricFlowRate},
        {kMolarFlow, UnitCategory::MolarFlowRate},
        {{0, 1, -1, 0, 0}, UnitCategory::MassFlowRate},
        {{-1, 1, -1, 0, 0}, UnitCategory::Viscosity},
        {{-3, 1, 0, 0, 0}, UnitCategory::Density},
        {{-3, 0, 0, 0, 1}, UnitCategory::MolarConcentration},
        {{2, 0, -2, -1, 0}, UnitCategory::HeatCapacity},
        {{0, 1, 0, 0, -1}, UnitCategory::MolarMass},
        {{1, 1, -3, -1, 0}, UnitCategory::ThermalConductivity},
        {{0, 1, -3, -1, 0}, UnitCategory::HeatTransferCoefficient},
        {{-2, 1, -2, 0, 0}, UnitCategory::PressureGradient},
        {{0, -1, 0, 0, 1}, UnitCategory::AdsorbedConcentration},
        {kMass, UnitCategory::Mass},
        {kVolume, UnitCategory::Volume},
        {kAmount, UnitCategory::MolarAmount},
        {{2, 0, 0, 0, 0}, UnitCategory::Area},
        {{2, 0, -1, 0, 0}, UnitCategory::Diffusivity},
        {{2, 0, -2, 0, 0}, UnitCategory::SpecificEnergy},
    };
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
    // Expressions re-parse their unit suffix on every evaluation; cache the
    // result per thread (the set of unit strings in a run is small).
    thread_local std::unordered_map<std::string, std::pair<bool, UnitInfo>> cache;
    thread_local std::unordered_map<std::string, std::string> errors;
    if (auto it = cache.find(t); it != cache.end()) {
        if (!it->second.first && error) *error = errors[t];
        if (it->second.first) out = it->second.second;
        return it->second.first;
    }
    UnitExpressionParser parser(t);
    std::string why;
    UnitInfo parsed;
    const bool ok = parser.Parse(parsed, why);
    if (cache.size() > 4096) {
        cache.clear();
        errors.clear();
    }
    cache.emplace(t, std::make_pair(ok, parsed));
    if (!ok) {
        errors[t] = why;
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
    UnitInfo info;
    if (!ParseUnitExpression(symbol, info)) {
        return false;  // unknown symbol — do not apply conversion
    }
    const UnitXform xform{info.factor, info.offset};
    const UnitCategory category = CategoryOf(symbol, info.dimension);

    // math_out is everything before the [symbol] bracket
    math_out  = expr.substr(0, static_cast<size_t>(m.position(0)));
    xform_out = xform;
    if (category_out) *category_out = category;
    return true;
}

} // namespace adsim
