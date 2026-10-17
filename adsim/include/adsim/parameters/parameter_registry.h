/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// The parameter registry: the single engine-side description of every value a
// user can read or set on a deck by name.
//
// Two kinds of entry share one ParameterInfo description:
//
//   * Static fields (StaticFields()) sit at a fixed place in every deck:
//     solver settings, cyclic-steady-state settings and the deck's info block.
//     Their accessors take a UserData, so they can be listed without a deck.
//     Front ends can build their flat settings forms from this table.
//
//   * Model parameters depend on the deck's structure (steps, streams, layers,
//     adsorbents, isotherm and rate models, user variables). ParameterRegistry
//     enumerates them for one UserData and binds each to the live object. The
//     optimiser (AvailableParameters in solver/optimisation.h), the Simulation
//     facade, the Python package and the CLI all read this list.
//
// Paths. Every entry has a stable `path` built from object IDs, so renaming a
// step, stream, layer or adsorbent does not change it:
//
//   cycle.total_time                       step[2].time
//   endpoint[5].temperature                endpoint[5].y[Nitrogen]
//   endpoint[3].pressure                   endpoint[5].pump.efficiency
//   bed.diameter      wall.thickness       layer[1].bed_void_fraction
//   adsorbent[1].solid_density             adsorbent[1].isotherm[CarbonDioxide].b1
//   adsorbent[1].kinetics[CarbonDioxide].k variable[feed_flow]
//   settings.reltol   css.min_cycles       info.title
//   energy_balance.isothermal_temperature  layer[1].momentum_model
//
// An endpoint is a flow stream or a pressure boundary (the same thing the KPI
// configuration calls an endpoint). Each entry also has an `alias` that uses
// display names instead of IDs (step[Feed].bottom[Stream 1].temperature,
// layer[Layer 1].length, ...). ParameterRegistry::Find accepts either; an alias
// that matches more than one entry (two streams with the same name) is
// rejected as ambiguous. Component names are CoolProp fluid names, which never
// change, so they appear as-is.
//
// Lifetime: a bound Parameter points into the UserData as it was when the
// registry was (re)built. Structural edits (changing a model type, adding or
// removing steps, streams, layers or user variables, loading a deck) destroy
// the referenced objects; call Refresh() after them and look entries up again
// by path.

#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

namespace adsim {

  class UserData;
  class ExpressionParam;
  enum class UnitCategory;

  enum class ParameterKind {
    Real,        // a plain double
    Expression,  // an ExpressionParam: a number or a formula of user variables
    Integer,     // an int / enum value
    Boolean,     // a switch (read and written as 0 / 1)
    String       // free text (deck info)
  };

  const char* ParameterKindName(ParameterKind kind);

  // Description of one parameter (static or bound).
  struct ParameterInfo {
    std::string path;         // stable, ID-based address (see above)
    std::string alias;        // the same address with display names (may equal path)
    std::string category;     // coarse group: cycle, bed, layer, adsorbent, isotherm,
                              // User Parameters, settings, css, info
    std::string subcategory;  // display names of the owning objects, e.g. "Feed, Stream 1"
    std::string name;         // leaf display name, e.g. "inlet temperature"
    std::string unit;         // SI unit symbol; "-" when dimensionless, "" when not numeric
    ParameterKind kind = ParameterKind::Real;
    // A numeric model value the optimiser and parametric studies may vary.
    // False for solver settings, integers, text, and values written as a
    // formula of user variables (overwriting those with a number would drop
    // the formula).
    bool tunable = false;
    // Suggested search bounds for tunable entries (NaN when there are none).
    double lower = std::numeric_limits<double>::quiet_NaN();
    double upper = std::numeric_limits<double>::quiet_NaN();
  };

  // A fixed-position deck field, accessed through a UserData. `section` and
  // `key` are the web contract's names (document section and JSON key); the
  // registry path is "<section>.<key>". Exactly one accessor group is set,
  // matching `kind` (Expression -> expr; Integer -> get_int/set_int; Boolean
  // -> get_bool/set_bool; String -> get_str/set_str). `web_unit` is the
  // optional unit-combo category the web document attaches to an Expression
  // field (e.g. "temperature"). Static fields are never tunable.
  struct StaticField {
    const char* section;
    const char* key;
    ParameterKind kind;
    const char* unit;
    const char* name;  // display name
    ExpressionParam& (*expr)(UserData&) = nullptr;
    long long (*get_int)(const UserData&) = nullptr;
    void (*set_int)(UserData&, long long) = nullptr;
    std::string (*get_str)(const UserData&) = nullptr;
    void (*set_str)(UserData&, const std::string&) = nullptr;
    bool (*get_bool)(const UserData&) = nullptr;
    void (*set_bool)(UserData&, bool) = nullptr;
    const char* web_unit = nullptr;
  };

  // Every static field, grouped by section in a fixed order (info, settings,
  // css, energy_balance). The order within a section is the web contract's
  // field order.
  const std::vector<StaticField>& StaticFields();
  // The static fields of one section, in order.
  std::vector<const StaticField*> StaticFieldsInSection(std::string_view section);

  // A parameter bound to a live UserData.
  class Parameter {
  public:
    const ParameterInfo& Info() const { return info_; }
    const std::string& Path() const { return info_.path; }

    // Numeric value (Real, Expression, Integer, Boolean as 0/1). Throws
    // std::logic_error for a String parameter.
    double GetValue() const;
    // Sets a numeric value. For an Expression parameter the expression text
    // becomes the number too, so a later re-evaluation keeps it. Integers are
    // rounded to the nearest whole number; a Boolean is true for any nonzero
    // value. Throws for a String parameter.
    void SetValue(double value);

    // The text form: the expression for an Expression parameter, the string
    // for a String parameter, the formatted number otherwise.
    std::string GetText() const;
    // Sets from text: an Expression parameter takes the expression (a number,
    // a formula of user variables, or a value with an inline [unit]); call
    // UserData::EvaluateAllExpressions afterwards to refresh its value. A
    // String parameter takes the text. A Boolean takes true/false/1/0. Real
    // and Integer parameters parse a number (std::invalid_argument when it is
    // not one).
    void SetText(const std::string& text);

    // Raw access for the optimiser adapter: the backing double (Real and
    // Expression parameters; nullptr otherwise) and the ExpressionParam
    // (Expression parameters only).
    double* ValuePointer() const { return value_; }
    ExpressionParam* ExpressionPointer() const { return expr_; }

  private:
    friend class ParameterRegistry;
    friend struct ParameterBuilder;
    ParameterInfo info_;
    double* value_ = nullptr;
    ExpressionParam* expr_ = nullptr;
    int* int_value_ = nullptr;
    std::function<long long()> get_int_;
    std::function<void(long long)> set_int_;
    std::function<std::string()> get_str_;
    std::function<void(const std::string&)> set_str_;
    std::function<bool()> get_bool_;
    std::function<void(bool)> set_bool_;
  };

  class ParameterRegistry {
  public:
    explicit ParameterRegistry(UserData& data);

    // Re-enumerates against the current state of the deck.
    void Refresh();

    // Every parameter: static fields first (in StaticFields() order), then the
    // model parameters (cycle, bed, layers, adsorbents, user variables).
    const std::vector<Parameter>& All() const { return parameters_; }
    std::vector<Parameter>& All() { return parameters_; }
    size_t Size() const { return parameters_.size(); }

    // Looks a parameter up by path or alias; nullptr if absent. Throws
    // std::invalid_argument when `path` is an alias shared by several entries.
    Parameter* Find(std::string_view path);
    const Parameter* Find(std::string_view path) const;
    // As Find, but throws std::out_of_range naming the path when it is absent.
    Parameter& Get(std::string_view path);

    // The tunable model parameters, in enumeration order.
    std::vector<Parameter*> Tunable();

    UserData& Data() { return data_; }

  private:
    UserData& data_;
    std::vector<Parameter> parameters_;
  };

  // SI symbol used for a user variable's unit category.
  const char* UnitCategorySymbol(UnitCategory category);

}  // namespace adsim
