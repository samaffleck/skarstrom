/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "adsim/parameters/parameter_registry.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include "adsim/math/unit_parser.h"
#include "adsim/objects/user_data.h"
#include "adsim/solver/cyclic_steady_state.h"

namespace adsim {

  const char* ParameterKindName(ParameterKind kind) {
    switch (kind) {
      case ParameterKind::Real:       return "real";
      case ParameterKind::Expression: return "expression";
      case ParameterKind::Integer:    return "integer";
      case ParameterKind::Boolean:    return "boolean";
      case ParameterKind::String:     return "string";
    }
    return "real";
  }

  const char* UnitCategorySymbol(UnitCategory cat) {
    switch (cat) {
      case UnitCategory::Time:                    return "s";
      case UnitCategory::Length:                  return "m";
      case UnitCategory::Temperature:             return "K";
      case UnitCategory::Pressure:                return "Pa";
      case UnitCategory::Energy:                  return "J";
      case UnitCategory::SpecificMolarEnergy:     return "J/mol";
      case UnitCategory::InverseTime:             return "1/s";
      case UnitCategory::Velocity:                return "m/s";
      case UnitCategory::VolumetricFlowRate:      return "m3/s";
      case UnitCategory::MolarFlowRate:           return "mol/s";
      case UnitCategory::MassFlowRate:            return "kg/s";
      case UnitCategory::MoleFraction:            return "mol/mol";
      case UnitCategory::Viscosity:               return "Pa.s";
      case UnitCategory::Density:                 return "kg/m3";
      case UnitCategory::MolarConcentration:      return "mol/m3";
      case UnitCategory::HeatCapacity:            return "J/(kg.K)";
      case UnitCategory::MolarMass:               return "kg/mol";
      case UnitCategory::ThermalConductivity:     return "W/(m.K)";
      case UnitCategory::HeatTransferCoefficient: return "W/(m2.K)";
      case UnitCategory::PressureGradient:        return "Pa/m";
      case UnitCategory::AdsorbedConcentration:   return "mol/kg";
      case UnitCategory::Mass:                    return "kg";
      case UnitCategory::Volume:                  return "m3";
      case UnitCategory::MolarAmount:             return "mol";
      case UnitCategory::Area:                    return "m2";
      case UnitCategory::Diffusivity:             return "m2/s";
      case UnitCategory::SpecificEnergy:          return "J/kg";
      default:                                    return "-";
    }
  }

  // -------------------------------------------------------------------------
  // Static fields
  // -------------------------------------------------------------------------

  namespace {

    // "none,aitken" -> {"none", "aitken"}
    std::vector<std::string> SplitChoices(const char* list) {
      std::vector<std::string> out;
      std::string item;
      for (const char* c = list; ; ++c) {
        if (*c == ',' || *c == '\0') {
          out.push_back(item);
          item.clear();
          if (*c == '\0') break;
        } else {
          item += *c;
        }
      }
      return out;
    }

    StaticField ExprField(const char* section, const char* key, const char* unit, const char* name,
                          ExpressionParam& (*expr)(UserData&), const char* web_unit = nullptr) {
      StaticField f{section, key, ParameterKind::Expression, unit, name};
      f.expr = expr;
      f.web_unit = web_unit;
      return f;
    }
    StaticField BoolField(const char* section, const char* key, const char* name,
                          bool (*get)(const UserData&), void (*set)(UserData&, bool)) {
      StaticField f{section, key, ParameterKind::Boolean, "-", name};
      f.get_bool = get;
      f.set_bool = set;
      return f;
    }
    StaticField IntField(const char* section, const char* key, const char* name,
                         long long (*get)(const UserData&), void (*set)(UserData&, long long),
                         const char* choices = nullptr) {
      StaticField f{section, key, ParameterKind::Integer, "-", name};
      f.get_int = get;
      f.set_int = set;
      f.choices = choices;
      return f;
    }
    StaticField StrField(const char* section, const char* key, const char* name,
                         std::string (*get)(const UserData&), void (*set)(UserData&, const std::string&)) {
      StaticField f{section, key, ParameterKind::String, "", name};
      f.get_str = get;
      f.set_str = set;
      return f;
    }

  }  // namespace

  const std::vector<StaticField>& StaticFields() {
    static const std::vector<StaticField> fields = {
      // info. sim_version is deliberately absent: it is a save counter, not an
      // editable field (the web document emits it directly).
      StrField("info", "title", "title",
               [](const UserData& u) { return u.info.title; },
               [](UserData& u, const std::string& s) { u.info.title = s; }),
      StrField("info", "description", "description",
               [](const UserData& u) { return u.info.description; },
               [](UserData& u, const std::string& s) { u.info.description = s; }),
      StrField("info", "author", "author",
               [](const UserData& u) { return u.info.author; },
               [](UserData& u, const std::string& s) { u.info.author = s; }),

      // settings (IDA integrator and logging)
      ExprField("settings", "init_step", "s", "initial step",
                [](UserData& u) -> ExpressionParam& { return u.settings.init_step; }),
      ExprField("settings", "min_step", "s", "minimum step",
                [](UserData& u) -> ExpressionParam& { return u.settings.min_step; }),
      ExprField("settings", "max_step", "s", "maximum step",
                [](UserData& u) -> ExpressionParam& { return u.settings.max_step; }),
      ExprField("settings", "reltol", "-", "relative tolerance",
                [](UserData& u) -> ExpressionParam& { return u.settings.reltol; }),
      ExprField("settings", "nonlin_conv_coef", "-", "nonlinear convergence coefficient",
                [](UserData& u) -> ExpressionParam& { return u.settings.nonlin_conv_coef; }),
      IntField("settings", "max_order", "maximum order",
               [](const UserData& u) { return (long long)u.settings.max_order; },
               [](UserData& u, long long x) { u.settings.max_order = (int)x; }),
      IntField("settings", "max_nonlin_iters", "maximum nonlinear iterations",
               [](const UserData& u) { return (long long)u.settings.max_nonlin_iters; },
               [](UserData& u, long long x) { u.settings.max_nonlin_iters = (int)x; }),
      IntField("settings", "max_conv_fails", "maximum convergence failures",
               [](const UserData& u) { return (long long)u.settings.max_conv_fails; },
               [](UserData& u, long long x) { u.settings.max_conv_fails = (int)x; }),
      IntField("settings", "max_err_test_fails", "maximum error-test failures",
               [](const UserData& u) { return (long long)u.settings.max_err_test_fails; },
               [](UserData& u, long long x) { u.settings.max_err_test_fails = (int)x; }),
      IntField("settings", "max_num_steps", "maximum internal steps",
               [](const UserData& u) { return (long long)u.settings.max_num_steps; },
               [](UserData& u, long long x) { u.settings.max_num_steps = (long int)x; }),
      IntField("settings", "max_iters_per_step", "maximum iterations per step",
               [](const UserData& u) { return (long long)u.settings.max_iters_per_step; },
               [](UserData& u, long long x) { u.settings.max_iters_per_step = (long int)x; }),
      IntField("settings", "advection_scheme", "advection scheme",
               [](const UserData& u) { return (long long)static_cast<int>(u.settings.advection_scheme_); },
               [](UserData& u, long long x) { u.settings.advection_scheme_ = static_cast<AdvectionScheme>((int)x); },
               "upwind,weno,van_leer"),
      IntField("settings", "max_cycles_to_keep", "cycles kept in the log",
               [](const UserData& u) { return (long long)u.settings.max_cycles_to_keep_; },
               [](UserData& u, long long x) { u.settings.max_cycles_to_keep_ = (int)x; }),
      // <= 0 disables the cap (log every internal step); see IDASettings.
      IntField("settings", "max_log_points_per_step", "maximum logged points per step",
               [](const UserData& u) { return (long long)u.settings.max_log_points_per_step_; },
               [](UserData& u, long long x) { u.settings.max_log_points_per_step_ = (int)x; }),
      IntField("settings", "acceleration_method", "CSS acceleration method",
               [](const UserData& u) { return (long long)static_cast<int>(u.settings.acc_css.method_); },
               [](UserData& u, long long x) { u.settings.acc_css.method_ = static_cast<AccelerationMethod>((int)x); },
               "none,aitken"),
      // Advisory warnings a run reports (warning_settings.h).
      BoolField("settings", "warn_pressure_step_change", "warn on a constant pressure boundary far from the bed pressure",
                [](const UserData& u) { return u.settings.warnings_.pressure_step_change; },
                [](UserData& u, bool b) { u.settings.warnings_.pressure_step_change = b; }),
      BoolField("settings", "warn_wall_thickness", "warn when the wall is thinner than its pressure design",
                [](const UserData& u) { return u.settings.warnings_.wall_thickness; },
                [](UserData& u, bool b) { u.settings.warnings_.wall_thickness = b; }),
      BoolField("settings", "warn_design_pressure", "warn when the pressure passes the wall's design pressure",
                [](const UserData& u) { return u.settings.warnings_.design_pressure; },
                [](UserData& u, bool b) { u.settings.warnings_.design_pressure = b; }),
      BoolField("settings", "warn_vacuum", "warn below atmospheric without vacuum service",
                [](const UserData& u) { return u.settings.warnings_.vacuum; },
                [](UserData& u, bool b) { u.settings.warnings_.vacuum = b; }),
      BoolField("settings", "warn_design_temperature", "warn when the wall passes its design temperature",
                [](const UserData& u) { return u.settings.warnings_.design_temperature; },
                [](UserData& u, bool b) { u.settings.warnings_.design_temperature = b; }),
      BoolField("settings", "warn_fluidisation", "warn when the approach to fluidisation passes 1",
                [](const UserData& u) { return u.settings.warnings_.fluidisation; },
                [](UserData& u, bool b) { u.settings.warnings_.fluidisation = b; }),

      // css (cyclic steady state)
      IntField("css", "method", "CSS method",
               [](const UserData& u) { return (long long)static_cast<int>(u.settings.css.method_); },
               [](UserData& u, long long x) { u.settings.css.method_ = static_cast<CSSMethod>((int)x); },
               "none,molar_profile,mass_balance_total,mass_balance_component"),
      IntField("css", "min_cycles", "minimum cycles",
               [](const UserData& u) { return (long long)u.settings.css.min_cycles_; },
               [](UserData& u, long long x) { u.settings.css.min_cycles_ = (int)x; }),
      ExprField("css", "molar_absolute_tolerance", "-", "molar absolute tolerance",
                [](UserData& u) -> ExpressionParam& { return u.settings.css.molar_absolute_tolerance; }),
      ExprField("css", "molar_relative_tolerance", "-", "molar relative tolerance",
                [](UserData& u) -> ExpressionParam& { return u.settings.css.molar_relative_tolerance; }),
      ExprField("css", "mass_balance_tolerance", "-", "mass balance tolerance",
                [](UserData& u) -> ExpressionParam& { return u.settings.css.mass_balance_tolerance_; }),

      // energy_balance: the isothermal switch (the web document emits these
      // inside its "energy_balance" object, next to hand-written fields).
      // 0 = non-isothermal (solve the energy balance), 1 = isothermal.
      IntField("energy_balance", "energy_balance_mode", "energy balance mode",
               [](const UserData& u) { return (long long)static_cast<int>(u.bed.equations_.energy_balance_mode_); },
               [](UserData& u, long long x) { u.bed.equations_.energy_balance_mode_ = static_cast<EnergyBalanceMode>((int)x); },
               "non_isothermal,isothermal"),
      // true: T_iso is the initial-condition temperature; false: isothermal_temperature.
      BoolField("energy_balance", "isothermal_use_initial_temperature",
                "isothermal temperature from initial conditions",
                [](const UserData& u) { return u.bed.equations_.isothermal_use_initial_temperature_; },
                [](UserData& u, bool b) { u.bed.equations_.isothermal_use_initial_temperature_ = b; }),
      ExprField("energy_balance", "isothermal_temperature", "K", "isothermal temperature",
                [](UserData& u) -> ExpressionParam& { return u.bed.equations_.isothermal_temperature_; },
                "temperature"),
    };
    return fields;
  }

  std::vector<const StaticField*> StaticFieldsInSection(std::string_view section) {
    std::vector<const StaticField*> out;
    for (const auto& f : StaticFields())
      if (section == f.section) out.push_back(&f);
    return out;
  }

  // -------------------------------------------------------------------------
  // Parameter
  // -------------------------------------------------------------------------

  double Parameter::GetValue() const {
    switch (info_.kind) {
      case ParameterKind::Real:       return *value_;
      case ParameterKind::Expression: return expr_->Get();
      case ParameterKind::Integer:    return static_cast<double>(get_int_());
      case ParameterKind::Boolean:    return get_bool_() ? 1.0 : 0.0;
      case ParameterKind::String:     break;
    }
    throw std::logic_error("parameter " + info_.path + " is text, not a number");
  }

  void Parameter::SetValue(double v) {
    switch (info_.kind) {
      case ParameterKind::Real:       *value_ = v; return;
      case ParameterKind::Expression: *expr_ = v; return;  // keeps the expression text in sync
      case ParameterKind::Integer: {
        const long long n = static_cast<long long>(std::llround(v));
        if (!info_.choices.empty() && (n < 0 || n >= static_cast<long long>(info_.choices.size())))
          throw std::invalid_argument("parameter " + info_.path + " takes 0.." +
                                      std::to_string(info_.choices.size() - 1) + " (or a name); got " +
                                      std::to_string(n));
        set_int_(n);
        return;
      }
      case ParameterKind::Boolean:    set_bool_(v != 0.0); return;
      case ParameterKind::String:     break;
    }
    throw std::logic_error("parameter " + info_.path + " is text, not a number");
  }

  std::string Parameter::GetText() const {
    switch (info_.kind) {
      case ParameterKind::Real:       return ExpressionParam::FormatScalar(*value_);
      case ParameterKind::Expression: return expr_->GetExpr();
      case ParameterKind::Integer:    return std::to_string(get_int_());
      case ParameterKind::Boolean:    return get_bool_() ? "true" : "false";
      case ParameterKind::String:     return get_str_();
    }
    return {};
  }

  namespace {
    // The unit text of a quantity ("2 min", "2 [min]") or of a formula with a
    // trailing "[unit]" ("2*t_feed [min]"); false when `text` has no unit.
    bool UnitOfText(const std::string& text, std::string& unit, bool& is_quantity, double& number) {
      if (SplitQuantity(text, number, unit)) {
        is_quantity = true;
        return true;
      }
      is_quantity = false;
      const size_t close = text.find_last_not_of(" \t");
      if (close == std::string::npos || text[close] != ']') return false;
      const size_t open = text.rfind('[', close);
      if (open == std::string::npos) return false;
      unit = text.substr(open + 1, close - open - 1);
      return true;
    }
  }  // namespace

  void Parameter::SetText(const std::string& text) {
    const bool numeric = info_.kind == ParameterKind::Real || info_.kind == ParameterKind::Expression;
    std::string unit;
    bool is_quantity = false;
    double number = 0.0;
    UnitInfo target, given;
    if (numeric && info_.unit_checked && UnitOfText(text, unit, is_quantity, number) &&
        ParseUnitExpression(info_.unit, target) && ParseUnitExpression(unit, given)) {
      if (!(given.dimension == target.dimension)) {
        const std::string want = DimensionName(target.dimension);
        const std::string got = DimensionName(given.dimension);
        throw std::invalid_argument("parameter " + info_.path + " is " +
                                    (target.dimension.IsDimensionless() ? std::string("dimensionless")
                                                                        : "a " + want + " [" + info_.unit + "]") +
                                    "; '" + text + "' is " +
                                    (given.dimension.IsDimensionless() ? std::string("dimensionless") : "a " + got));
      }
      const bool coherent = target.factor == 1.0 && target.offset == 0.0;
      if (info_.kind == ParameterKind::Expression && coherent) {
        // The expression evaluates "[unit]" to SI, which is this parameter's unit.
        expr_->SetExpr(is_quantity ? ExpressionParam::FormatScalar(number) + " [" + unit + "]" : text);
        return;
      }
      if (!is_quantity)
        throw std::invalid_argument("parameter " + info_.path + " is stored in " + info_.unit +
                                    ", not SI; give a number in " + info_.unit + " or a single quantity such as '1 [" +
                                    info_.unit + "]'");
      SetValue(ConvertUnit(number, unit, info_.unit));
      return;
    }
    switch (info_.kind) {
      case ParameterKind::String:     set_str_(text); return;
      case ParameterKind::Expression: expr_->SetExpr(text); return;
      case ParameterKind::Boolean: {
        std::string t;
        for (char c : text)
          if (!std::isspace(static_cast<unsigned char>(c)))
            t += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (t == "true" || t == "1" || t == "yes" || t == "on") { set_bool_(true); return; }
        if (t == "false" || t == "0" || t == "no" || t == "off") { set_bool_(false); return; }
        throw std::invalid_argument("parameter " + info_.path + " needs true or false, got '" + text + "'");
      }
      case ParameterKind::Real:
      case ParameterKind::Integer: {
        if (info_.kind == ParameterKind::Integer && !info_.choices.empty()) {
          // A choice by name ("mass_balance_component", "Ergun", "van-leer").
          std::string key;
          for (char c : text) {
            if (std::isspace(static_cast<unsigned char>(c))) continue;
            key += (c == '-' || c == ' ') ? '_' : static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
          }
          for (size_t i = 0; i < info_.choices.size(); ++i) {
            if (key == info_.choices[i]) {
              set_int_(static_cast<long long>(i));
              return;
            }
          }
          bool numeric = !key.empty();
          for (char c : key) numeric = numeric && (std::isdigit(static_cast<unsigned char>(c)) || c == '.');
          if (!numeric) {
            std::string list;
            for (const auto& c : info_.choices) list += (list.empty() ? "" : ", ") + c;
            throw std::invalid_argument("parameter " + info_.path + " is one of: " + list + "; got '" + text + "'");
          }
        }
        size_t used = 0;
        double v = 0.0;
        try {
          v = std::stod(text, &used);
        } catch (const std::exception&) {
          used = 0;
        }
        while (used < text.size() && std::isspace(static_cast<unsigned char>(text[used]))) ++used;
        if (used == 0 || used != text.size())
          throw std::invalid_argument("parameter " + info_.path + " needs a number, got '" + text + "'");
        SetValue(v);
        return;
      }
    }
  }

  // -------------------------------------------------------------------------
  // Enumeration
  // -------------------------------------------------------------------------

  // Model parameters are addressed by the stored IDs of the objects that own
  // them: Step::id_, the endpoint_id_ of streams and pressure boundaries
  // (Cycle::EnsureIds assigns any missing one first), MixedAdsorbentLayer::id_
  // and the AdsorbentID.
  struct ParameterBuilder {
    std::vector<Parameter>& out;

    static std::string Bracket(const std::string& head, const std::string& id) {
      return head + "[" + id + "]";
    }

    Parameter& Add(ParameterInfo info) {
      Parameter p;
      p.info_ = std::move(info);
      out.push_back(std::move(p));
      return out.back();
    }

    ParameterInfo Info(std::string path, std::string alias, const std::string& category,
                       const std::string& subcategory, const std::string& name,
                       const std::string& unit, ParameterKind kind, double lo, double hi) {
      ParameterInfo i;
      i.path = std::move(path);
      i.alias = alias.empty() ? i.path : std::move(alias);
      i.category = category;
      i.subcategory = subcategory;
      i.name = name;
      i.unit = unit;
      i.kind = kind;
      i.lower = lo;
      i.upper = hi;
      return i;
    }

    void Emplace(const std::string& path, const std::string& alias, const std::string& category,
                 const std::string& subcategory, const std::string& name, const std::string& unit,
                 ExpressionParam& param, double lo, double hi) {
      auto& p = Add(Info(path, alias, category, subcategory, name, unit, ParameterKind::Expression, lo, hi));
      p.expr_ = &param;
      p.value_ = &param.Get();
      p.info_.tunable = !param.HasExpression();
      WidenBoundsTo(p.info_, param.Get());
    }

    // Suggested bounds always contain the current value (the defaults are
    // typical ranges; a deck may sit outside them).
    static void WidenBoundsTo(ParameterInfo& i, double v) {
      if (!std::isfinite(v)) return;
      if (v < i.lower) i.lower = v;
      if (v > i.upper) i.upper = v;
    }

    void Emplace(const std::string& path, const std::string& alias, const std::string& category,
                 const std::string& subcategory, const std::string& name, const std::string& unit,
                 double& value, double lo, double hi) {
      auto& p = Add(Info(path, alias, category, subcategory, name, unit, ParameterKind::Real, lo, hi));
      p.value_ = &value;
      p.info_.tunable = true;
      WidenBoundsTo(p.info_, value);
    }

    void EmplaceInt(const std::string& path, const std::string& category, const std::string& subcategory,
                    const std::string& name, int& value, const std::string& alias = "") {
      auto& p = Add(Info(path, alias, category, subcategory, name, "-", ParameterKind::Integer,
                         std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));
      int* ptr = &value;
      p.int_value_ = ptr;
      p.get_int_ = [ptr] { return static_cast<long long>(*ptr); };
      p.set_int_ = [ptr](long long x) { *ptr = static_cast<int>(x); };
    }

    // An integer switch stored as an enum (or any int-like field).
    template <class E>
    void EmplaceEnum(const std::string& path, const std::string& alias, const std::string& category,
                     const std::string& subcategory, const std::string& name, E& value, const char* choices) {
      auto& p = Add(Info(path, alias, category, subcategory, name, "-", ParameterKind::Integer,
                         std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));
      p.info_.choices = SplitChoices(choices);
      E* ptr = &value;
      p.get_int_ = [ptr] { return static_cast<long long>(static_cast<int>(*ptr)); };
      p.set_int_ = [ptr](long long x) { *ptr = static_cast<E>(static_cast<int>(x)); };
    }

    void EmplaceBool(const std::string& path, const std::string& alias, const std::string& category,
                     const std::string& subcategory, const std::string& name, bool& value) {
      auto& p = Add(Info(path, alias, category, subcategory, name, "-", ParameterKind::Boolean,
                         std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));
      bool* ptr = &value;
      p.get_bool_ = [ptr] { return *ptr; };
      p.set_bool_ = [ptr](bool b) { *ptr = b; };
    }

    // ---- static fields ---------------------------------------------------
    void AddStatic(UserData& data) {
      for (const auto& f : StaticFields()) {
        const std::string path = std::string(f.section) + "." + f.key;
        auto& p = Add(Info(path, path, f.section, f.section, f.name, f.unit, f.kind,
                           std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));
        if (f.choices) p.info_.choices = SplitChoices(f.choices);
        UserData* u = &data;
        switch (f.kind) {
          case ParameterKind::Expression:
            p.expr_ = &f.expr(data);
            p.value_ = &p.expr_->Get();
            break;
          case ParameterKind::Integer: {
            auto get = f.get_int;
            auto set = f.set_int;
            p.get_int_ = [u, get] { return get(*u); };
            p.set_int_ = [u, set](long long x) { set(*u, x); };
            break;
          }
          case ParameterKind::String: {
            auto get = f.get_str;
            auto set = f.set_str;
            p.get_str_ = [u, get] { return get(*u); };
            p.set_str_ = [u, set](const std::string& s) { set(*u, s); };
            break;
          }
          case ParameterKind::Boolean: {
            auto get = f.get_bool;
            auto set = f.set_bool;
            p.get_bool_ = [u, get] { return get(*u); };
            p.set_bool_ = [u, set](bool b) { set(*u, b); };
            break;
          }
          case ParameterKind::Real:
            break;
        }
      }
    }

    // ---- cycle -----------------------------------------------------------
    const Bed* bed_ = nullptr;  // layer names for energy-source aliases

    void AddCycle(Cycle& cycle) {
      cycle.EnsureIds();
      EmplaceInt("cycle.number_of_cycles", "cycle", "cycle", "number of cycles", cycle.number_of_cycles_);
      // 0 = step times in seconds, 1 = fractions of cycle.total_time.
      EmplaceEnum("cycle.step_time_mode", "", "cycle", "cycle", "step time mode", cycle.step_time_mode_,
                  "time,fraction");
      if (cycle.step_time_mode_ == StepTimeMode::FRACTION) {
        Emplace("cycle.total_time", "", "cycle", "cycle", "total cycle time", "s",
                cycle.total_cycle_time_, 1.0, 10000.0);
      }
      for (auto& step : cycle.step_list_) AddStep(step, cycle.step_time_mode_);
    }

    struct EndpointScope {
      std::string path;     // "endpoint[5]"
      std::string alias;    // "step[Feed].bottom[Stream 1]" or "step[Feed].bottom"
      std::string subcategory;
    };

    void AddStep(Step& step, StepTimeMode mode) {
      const std::string sp = Bracket("step", std::to_string(step.id_));
      const std::string sa = Bracket("step", step.step_name_);
      if (mode == StepTimeMode::FRACTION) {
        Emplace(sp + ".fraction", sa + ".fraction", "cycle", step.step_name_, "step fraction", "-",
                step.step_fraction_, 0.0, 1.0);
      } else {
        Emplace(sp + ".time", sa + ".time", "cycle", step.step_name_, "step time", "s",
                step.step_time_, 0.0, 1000.0);
      }
      AddBoundary(step.bottom_boundary_condition_.get(), step.step_name_, sa + ".bottom");
      AddBoundary(step.top_boundary_condition_.get(), step.step_name_, sa + ".top");
      // Energy sources (heaters and coolers) acting on a layer during this step.
      for (auto& [layer_id, source] : step.energy_sources_) {
        std::string layer_name = std::to_string(layer_id);
        if (bed_)
          for (const auto& l : bed_->layers_.layers_)
            if (l.id_ == layer_id) layer_name = l.name_;
        const std::string ep = sp + ".energy[" + std::to_string(layer_id) + "]";
        const std::string ea = sa + ".energy[" + layer_name + "]";
        const std::string sub = step.step_name_ + ", " + layer_name;
        EmplaceBool(ep + ".enabled", ea + ".enabled", "cycle", sub, "energy source enabled", source.enabled_);
        Emplace(ep + ".power", ea + ".power", "cycle", sub, "energy source power", "W", source.power_, -1e4, 1e4);
      }
    }

    void AddBoundary(IBoundaryCondition* boundary, const std::string& step_name, const std::string& side_alias) {
      if (auto* bc = dynamic_cast<PressureBoundaryCondition*>(boundary)) {
        EndpointScope ep{Bracket("endpoint", std::to_string(bc->endpoint_id_)), side_alias, step_name};
        AddPressureProfile(bc->pressure_.get(), ep);
        // An outlet's temperature and composition apply when the flow reverses.
        const bool inlet = bc->direction_ == BoundaryDirection::INLET;
        AddInletTemperature(bc->temperature_, bc->temperature_profile_.get(), ep, inlet);
        AddMoleFractions(bc->mole_fraction_, ep, inlet);
        AddPump(bc->pump_.get(), ep);
        AddBufferTank(bc->buffer_tank_, ep);
      } else if (auto* bc = dynamic_cast<FlowBoundaryCondition*>(boundary)) {
        for (auto& stream : bc->streams_) {
          if (!stream) continue;
          EndpointScope ep{Bracket("endpoint", std::to_string(stream->endpoint_id_)),
                           side_alias + "[" + stream->name_ + "]", step_name + ", " + stream->name_};
          AddStream(stream.get(), bc->direction_, ep);
        }
      }
    }

    void AddStream(IBoundaryStream* stream, BoundaryDirection direction, const EndpointScope& ep) {
      const bool inlet = direction == BoundaryDirection::INLET;
      AddInletTemperature(stream->temperature_, stream->temperature_profile_.get(), ep, inlet);
      AddMoleFractions(stream->mole_fraction, ep, inlet);
      if (auto* s = dynamic_cast<BoundaryStream*>(stream)) {
        AddFlow(s->boundary_flow_.get(), ep);
        AddPump(s->pump_.get(), ep);
      } else if (auto* v = dynamic_cast<ValveBoundaryStream*>(stream)) {
        AddPressureProfile(v->boundary_value_.get(), ep);
        Emplace(ep.path + ".valve.cv", ep.alias + ".valve.cv", "cycle", ep.subcategory,
                "valve flow coefficient (Cv, US gpm/psi^0.5)", "-", v->valve_.cv_, 1e-3, 1e3);
        Emplace(ep.path + ".valve.xt", ep.alias + ".valve.xt", "cycle", ep.subcategory,
                "valve pressure-drop ratio factor (xT)", "-", v->valve_.xt_, 0.1, 1.0);
      } else if (auto* p = dynamic_cast<PumpBoundaryStream*>(stream)) {
        AddPumpPower(p->power_.get(), ep);
        AddPump(p->pump_.get(), ep);
      }
      AddBufferTank(stream->buffer_tank_, ep);
    }

    // The shaft power that drives a pump stream: a constant, or a polynomial
    // of the boundary pressure P [Pa], a + b P + ... + f P^5.
    void AddPumpPower(IPumpPowerModel* power, const EndpointScope& ep) {
      auto E = [&](const std::string& leaf, const std::string& name, const std::string& unit, ExpressionParam& v,
                   double lo, double hi) {
        Emplace(ep.path + "." + leaf, ep.alias + "." + leaf, "cycle", ep.subcategory, name, unit, v, lo, hi);
      };
      if (auto* c = dynamic_cast<ConstantPumpPowerModel*>(power)) {
        E("power", "pump shaft power", "W", c->power_, 0.0, 1e6);
      } else if (auto* poly = dynamic_cast<PolyPumpPowerModel*>(power)) {
        E("power_poly_a", "pump power polynomial a", "W", poly->poly_.a, -1e6, 1e6);
        E("power_poly_b", "pump power polynomial b (x P)", "W/Pa", poly->poly_.b, -1e3, 1e3);
        E("power_poly_c", "pump power polynomial c (x P^2)", "W/Pa^2", poly->poly_.c, -1.0, 1.0);
        E("power_poly_d", "pump power polynomial d (x P^3)", "W/Pa^3", poly->poly_.d, -1e-3, 1e-3);
        E("power_poly_e", "pump power polynomial e (x P^4)", "W/Pa^4", poly->poly_.e, -1e-6, 1e-6);
        E("power_poly_f", "pump power polynomial f (x P^5)", "W/Pa^5", poly->poly_.f, -1e-9, 1e-9);
      }
    }

    // The share of its source's moles an inlet draws from a buffer tank.
    void AddBufferTank(BufferTankSource& tank, const EndpointScope& ep) {
      if (!tank.IsActive()) return;
      Emplace(ep.path + ".buffer_tank.percentage", ep.alias + ".buffer_tank.percentage", "cycle", ep.subcategory,
              "buffer tank: share of the source endpoint's moles drawn", "%", tank.percentage_, 0.0, 100.0);
    }

    void AddMoleFractions(ExpressionMap& y, const EndpointScope& ep, bool inlet) {
      for (auto& [component, yi] : y) {
        Emplace(ep.path + ".y[" + component + "]", ep.alias + ".y[" + component + "]", "cycle", ep.subcategory,
                component + (inlet ? " inlet mole fraction" : " backflow mole fraction"), "-", yi, 0.0, 1.0);
      }
    }

    // The start, end and duration of a linear, exponential or log-linear profile,
    // as <head>_linear_start/_end/_time (or _exponential_..., _log_linear_...); a
    // constant one as <head>.
    void AddProfile(IBoundaryValue* value, const EndpointScope& ep, const std::string& head,
                    const std::string& name, const std::string& unit, double lo, double hi) {
      auto E = [&](const std::string& leaf, const std::string& label, const std::string& u, ExpressionParam& v,
                   double l, double h) {
        Emplace(ep.path + "." + leaf, ep.alias + "." + leaf, "cycle", ep.subcategory, label, u, v, l, h);
      };
      if (auto* c = dynamic_cast<ConstantBoundaryValue*>(value)) {
        E(head, name, unit, c->constant_value_, lo, hi);
      } else if (auto* lin = dynamic_cast<LinearBoundaryValue*>(value)) {
        E(head + "_linear_start", name + " (linear start)", unit, lin->start_value_, lo, hi);
        E(head + "_linear_end", name + " (linear end)", unit, lin->end_value_, lo, hi);
        E(head + "_linear_time", name + " (linear ramp time)", "s", lin->ramp_time_, 0.0, 1e4);
      } else if (auto* ex = dynamic_cast<ExponentialBoundaryValue*>(value)) {
        E(head + "_exponential_start", name + " (exponential start)", unit, ex->start_value_, lo, hi);
        E(head + "_exponential_end", name + " (exponential end)", unit, ex->end_value_, lo, hi);
        E(head + "_exponential_time", name + " (exponential time)", "s", ex->step_time_, 0.0, 1e4);
      } else if (auto* lg = dynamic_cast<LogLinearBoundaryValue*>(value)) {
        E(head + "_log_linear_start", name + " (log-linear start)", unit, lg->start_value_, lo, hi);
        E(head + "_log_linear_end", name + " (log-linear end)", unit, lg->end_value_, lo, hi);
        E(head + "_log_linear_time", name + " (log-linear ramp time)", "s", lg->ramp_time_, 0.0, 1e4);
      } else if (auto* de = dynamic_cast<DoubleExponentialBoundaryValue*>(value)) {
        E(head + "_double_exponential_start", name + " (double-exponential start)", unit, de->start_value_, lo, hi);
        E(head + "_double_exponential_end", name + " (double-exponential end)", unit, de->end_value_, lo, hi);
        E(head + "_double_exponential_rate_1", name + " (double-exponential rate 1)", "1/s", de->rate_1_, 0.0, 10.0);
        E(head + "_double_exponential_rate_2", name + " (double-exponential rate 2)", "1/s", de->rate_2_, 0.0, 10.0);
        E(head + "_double_exponential_weight", name + " (double-exponential weight of term 1)", "-", de->weight_, 0.0, 1.0);
      }
    }

    // A constant inlet temperature, or the start and end of a ramp (the ramp
    // replaces the constant value, so the constant is not listed then).
    void AddInletTemperature(ExpressionParam& constant, IBoundaryValue* profile, const EndpointScope& ep,
                             bool inlet = true) {
      const std::string name = inlet ? "inlet temperature" : "backflow temperature";
      if (profile && !dynamic_cast<ConstantBoundaryValue*>(profile)) {
        AddProfile(profile, ep, "temperature", name, "K", 273.0, 473.0);
      } else {
        Emplace(ep.path + ".temperature", ep.alias + ".temperature", "cycle", ep.subcategory, name, "K", constant,
                273.0, 473.0);
      }
    }

    void AddFlow(IFlowBoundaryValue* flow, const EndpointScope& ep) {
      if (!flow) return;
      IBoundaryValue* profile = flow->profile_.get();
      switch (flow->GetType()) {
        case FlowType::VELOCITY: AddProfile(profile, ep, "velocity", "velocity", "m/s", 0.0, 10.0); break;
        case FlowType::MOLAR:    AddProfile(profile, ep, "molar_flow", "molar flow", "mol/s", 0.0, 100.0); break;
        case FlowType::VOLUME:   AddProfile(profile, ep, "volumetric_flow", "volumetric flow", "m3/s", 0.0, 1.0); break;
        case FlowType::MASS:     AddProfile(profile, ep, "mass_flow", "mass flow", "kg/s", 0.0, 10.0); break;
      }
    }

    void AddPump(IPump* pump, const EndpointScope& ep) {
      if (!pump || !pump->use_pump_) return;
      auto E = [&](const std::string& leaf, const std::string& name, const std::string& unit, ExpressionParam& v,
                   double lo, double hi) {
        Emplace(ep.path + ".pump." + leaf, ep.alias + ".pump." + leaf, "cycle", ep.subcategory, name, unit, v, lo, hi);
      };
      E("polytropic_index", "pump polytropic index", "-", pump->poly_index_, 1.0, 1.7);
      if (auto* eff = dynamic_cast<ConstantPumpEfficiencyModel*>(pump->efficiency_model_.get())) {
        E("efficiency", "pump efficiency", "-", eff->efficiency_, 0.01, 1.0);
      } else if (auto* eff = dynamic_cast<VacuumPumpEfficiencyModel*>(pump->efficiency_model_.get())) {
        E("eta_max", "vacuum pump peak efficiency", "-", eff->eta_max_, 0.01, 1.0);
        E("alpha", "vacuum pump efficiency decay", "-", eff->alpha_, 0.0, 100.0);
      } else if (auto* eff = dynamic_cast<PolyPumpEfficiencyModel*>(pump->efficiency_model_.get())) {
        E("poly_a", "pump efficiency polynomial a", "-", eff->a_, -1e3, 1e3);
        E("poly_b", "pump efficiency polynomial b (x P)", "1/Pa", eff->b_, -1e3, 1e3);
        E("poly_c", "pump efficiency polynomial c (x P^2)", "1/Pa^2", eff->c_, -1e3, 1e3);
        E("poly_d", "pump efficiency polynomial d (x P^3)", "1/Pa^3", eff->d_, -1e3, 1e3);
        E("poly_e", "pump efficiency polynomial e (x P^4)", "1/Pa^4", eff->e_, -1e3, 1e3);
        E("poly_f", "pump efficiency polynomial f (x P^5)", "1/Pa^5", eff->f_, -1e3, 1e3);
      }
      if (auto* c = dynamic_cast<Compressor*>(pump)) {
        Emplace(ep.path + ".pump.suction_pressure", ep.alias + ".pump.suction_pressure", "cycle",
                ep.subcategory, "compressor suction pressure", "Pa", c->suction_pressure_, 1e3, 1e7);
      } else if (auto* v = dynamic_cast<VacuumPump*>(pump)) {
        Emplace(ep.path + ".pump.discharge_pressure", ep.alias + ".pump.discharge_pressure", "cycle",
                ep.subcategory, "vacuum pump discharge pressure", "Pa", v->discharge_pressure_, 1e3, 1e7);
      }
    }

    void AddPressureProfile(IBoundaryValue* value, const EndpointScope& ep) {
      if (auto* bv = dynamic_cast<ConstantBoundaryValue*>(value)) {
        Emplace(ep.path + ".pressure", ep.alias + ".pressure", "cycle", ep.subcategory,
                "pressure (constant)", "Pa", bv->constant_value_, 1e4, 1e6);
      } else if (auto* bv = dynamic_cast<LinearBoundaryValue*>(value)) {
        Emplace(ep.path + ".pressure_linear_end", ep.alias + ".pressure_linear_end", "cycle",
                ep.subcategory, "pressure (linear end)", "Pa", bv->end_value_, 1e4, 1e6);
      } else if (auto* bv = dynamic_cast<ExponentialBoundaryValue*>(value)) {
        Emplace(ep.path + ".pressure_exponential_end", ep.alias + ".pressure_exponential_end", "cycle",
                ep.subcategory, "pressure (exponential end)", "Pa", bv->end_value_, 1e4, 1e6);
      } else if (auto* bv = dynamic_cast<LogLinearBoundaryValue*>(value)) {
        Emplace(ep.path + ".pressure_log_linear_end", ep.alias + ".pressure_log_linear_end", "cycle",
                ep.subcategory, "pressure (log-linear end)", "Pa", bv->end_value_, 1e4, 1e6);
      } else if (auto* bv = dynamic_cast<DoubleExponentialBoundaryValue*>(value)) {
        Emplace(ep.path + ".pressure_double_exponential_end", ep.alias + ".pressure_double_exponential_end", "cycle",
                ep.subcategory, "pressure (double-exponential end)", "Pa", bv->end_value_, 1e4, 1e6);
        Emplace(ep.path + ".pressure_double_exponential_rate_1", ep.alias + ".pressure_double_exponential_rate_1", "cycle",
                ep.subcategory, "pressure (double-exponential rate 1)", "1/s", bv->rate_1_, 0.0, 10.0);
        Emplace(ep.path + ".pressure_double_exponential_rate_2", ep.alias + ".pressure_double_exponential_rate_2", "cycle",
                ep.subcategory, "pressure (double-exponential rate 2)", "1/s", bv->rate_2_, 0.0, 10.0);
        Emplace(ep.path + ".pressure_double_exponential_weight", ep.alias + ".pressure_double_exponential_weight", "cycle",
                ep.subcategory, "pressure (double-exponential weight of term 1)", "-", bv->weight_, 0.0, 1.0);
      }
    }

    // ---- bed, layers, adsorbents ----------------------------------------
    std::set<AdsorbentID> emitted_adsorbents;

    void AddBed(Bed& bed, AdsorbentDirectory& dir) {
      if (auto* cyl = std::get_if<CylinderBed>(&bed.shape_)) {
        Emplace("bed.diameter", "", "bed", "geometry", "diameter", "m", cyl->diameter, 0.01, 1.0);
      } else if (auto* sq = std::get_if<SquareBed>(&bed.shape_)) {
        Emplace("bed.width", "", "bed", "geometry", "width", "m", sq->width, 0.01, 1.0);
      } else if (auto* rect = std::get_if<RectangularBed>(&bed.shape_)) {
        Emplace("bed.width", "", "bed", "geometry", "width", "m", rect->width, 0.01, 1.0);
        Emplace("bed.length", "", "bed", "geometry", "length", "m", rect->length, 0.01, 1.0);
      }

      emitted_adsorbents.clear();
      for (auto& layer : bed.layers_.layers_) AddLayer(layer, dir);
      // Adsorbents no layer uses yet (a deck being built, or spare materials).
      std::vector<AdsorbentID> spare;
      for (auto& [id, ads] : dir.GetAll()) {
        (void)ads;
        if (!emitted_adsorbents.count(id)) spare.push_back(id);
      }
      std::sort(spare.begin(), spare.end());
      for (AdsorbentID id : spare) {
        emitted_adsorbents.insert(id);
        AddAdsorbent(id, dir.Get(id), "");
      }

      Emplace("bed.ambient_temperature", "", "bed", "environment", "ambient temperature", "K",
              bed.ambient_temperature_, 273.0, 373.0);
      Emplace("wall.density", "", "bed", "wall", "wall density", "kg/m3", bed.wall_.density_, 1000.0, 10000.0);
      Emplace("wall.heat_capacity", "", "bed", "wall", "wall heat capacity", "J/(kg-K)",
              bed.wall_.heat_capacity_, 100.0, 2000.0);
      Emplace("wall.thermal_conductivity", "", "bed", "wall", "wall thermal conductivity", "W/(m-K)",
              bed.wall_.thermal_conductivity_, 0.1, 100.0);
      Emplace("wall.thickness", "", "bed", "wall", "wall thickness", "m", bed.wall_.thickness_, 0.001, 0.1);
      // Pressure design (wall_design.h): sizes the thickness, sets the
      // operating-limit warnings; the solver does not read it.
      EmplaceEnum("wall.design_code", "", "bed", "wall design", "wall design code", bed.wall_.design_code_,
                  "asme_viii_1,asme_b31_3");
      Emplace("wall.design_pressure", "", "bed", "wall design", "wall design pressure (gauge)", "Pa",
              bed.wall_.design_pressure_, 0.0, 1e7);
      Emplace("wall.design_temperature", "", "bed", "wall design", "wall design temperature", "K",
              bed.wall_.design_temperature_, 273.0, 873.0);
      Emplace("wall.weld_joint_efficiency", "", "bed", "wall design", "weld joint efficiency", "-",
              bed.wall_.weld_joint_efficiency_, 0.5, 1.0);
      Emplace("wall.corrosion_allowance", "", "bed", "wall design", "corrosion allowance", "m",
              bed.wall_.corrosion_allowance_, 0.0, 0.006);
      EmplaceBool("wall.vacuum_service", "", "bed", "wall design", "vacuum service", bed.wall_.vacuum_service_);
      auto& wall_bc = bed.equations_.wall_temperature_boundary_condition_;
      if (auto* fixed = dynamic_cast<WallTemperatureBoundaryConditionFixed*>(wall_bc.get())) {
        Emplace("wall.fixed_temperature", "", "bed", "wall", "fixed wall temperature", "K",
                fixed->temperature_, 273.0, 473.0);
      }

      // 0 = local thermal equilibrium (LTE), 1 = non-equilibrium (LTNE).
      EmplaceEnum("bed.thermal_equilibrium_mode", "", "bed", "thermal", "thermal equilibrium mode",
                  bed.equations_.thermal_equilibrium_mode_, "lte,ltne");
      if (auto* h = dynamic_cast<ConstantParticleHeatTransferModel*>(bed.equations_.particle_heat_transfer_model_.get()))
        Emplace("bed.particle_htc", "", "bed", "thermal", "gas-to-particle heat transfer coefficient", "W/m2-K",
                h->heat_transfer_coefficient_, 1.0, 1e5);
      auto void_region = [&](const std::string& head, VoidRegion& v, const std::string& where) {
        EmplaceBool(head + ".enabled", "", "bed", where, where + " enabled", v.enabled_);
        Emplace(head + ".volume", "", "bed", where, where + " volume", "m3", v.volume_, 1e-7, 1.0);
        Emplace(head + ".resistance_coefficient", "", "bed", where, where + " flow resistance coefficient",
                "Pa.s/m", v.resistance_coefficient_, 0.0, 1e6);
      };
      void_region("bed.bottom_void", bed.bottom_void_, "bottom void");
      void_region("bed.top_void", bed.top_void_, "top void");
    }

    // ---- initial conditions and fluid -------------------------------------
    void AddInitialConditions(InitialConditions& ic) {
      Emplace("initial.pressure", "", "initial", "initial", "initial pressure", "Pa", ic.pressure, 1e3, 1e7);
      Emplace("initial.temperature", "", "initial", "initial", "initial temperature", "K", ic.temperature, 200.0, 800.0);
      Emplace("initial.wall_temperature", "", "initial", "initial", "initial wall temperature", "K",
              ic.wall_temperature, 200.0, 800.0);
      Emplace("initial.velocity", "", "initial", "initial", "initial velocity", "m/s", ic.velocity, -10.0, 10.0);
      for (auto& [component, yi] : ic.mole_fraction)
        Emplace("initial.y[" + component + "]", "", "initial", "initial", component + " initial mole fraction", "-",
                yi, 0.0, 1.0);
    }

    void AddFluid(Fluid& fluid) {
      // 0 = constant (fluid.adsorbed_heat_capacity), 1 = equal to the gas-phase value.
      EmplaceEnum("fluid.adsorbed_heat_capacity_mode", "", "fluid", "fluid", "adsorbed-phase heat capacity mode",
                  fluid.adsorbed_heat_capacity_mode_, "constant,gas_phase");
      // 0 = isosteric heat of adsorption, 1 = internal-energy convention.
      EmplaceEnum("fluid.heat_of_adsorption_convention", "", "fluid", "fluid", "heat of adsorption convention",
                  fluid.heat_of_adsorption_convention_, "isosteric,internal_energy");
      for (const auto& component : fluid.components_) {
        auto it = fluid.adsorbed_heat_capacity_.find(component);
        if (it == fluid.adsorbed_heat_capacity_.end()) continue;
        Emplace("fluid.adsorbed_heat_capacity[" + component + "]", "", "fluid", "fluid",
                component + " adsorbed-phase heat capacity", "J/(mol-K)", it->second, 0.0, 500.0);
      }
    }

    void AddLayer(MixedAdsorbentLayer& layer, AdsorbentDirectory& dir) {
      const std::string& ln = layer.name_;
      const std::string lp = Bracket("layer", std::to_string(layer.id_));
      const std::string la = Bracket("layer", ln);
      auto E = [&](const std::string& leaf, const std::string& name, const std::string& unit,
                   auto& value, double lo, double hi) {
        Emplace(lp + "." + leaf, la + "." + leaf, "layer", ln, name, unit, value, lo, hi);
      };
      E("length", "length", "m", layer.domain_.length, 0.01, 10.0);
      EmplaceInt(lp + ".cells", "layer", ln, "number of cells", layer.domain_.n_cells, la + ".cells");
      for (auto& [ads_id, fraction] : layer.composition_) {
        const std::string an = dir.Contains(ads_id) ? dir.Get(ads_id).name_ : std::to_string(ads_id);
        Emplace(lp + ".mass_fraction[" + std::to_string(ads_id) + "]", la + ".mass_fraction[" + an + "]", "layer", ln,
                an + " mass fraction", "-", fraction, 0.0, 1.0);
      }
      E("bed_void_fraction", "bed void fraction", "-", layer.eb_, 0.2, 0.8);
      E("particle_diameter", "particle diameter", "m", layer.particle_diameter_, 1e-5, 0.01);
      {
        // 0 = Darcy (Kozeny-Carman), 1 = Ergun, 2 = Darcy (Blake-Kozeny).
        auto& p = Add(Info(lp + ".momentum_model", la + ".momentum_model", "layer", ln, "momentum model", "-",
                           ParameterKind::Integer, std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::quiet_NaN()));
        p.info_.choices = {"kozeny_carman", "ergun", "blake_kozeny"};
        MixedAdsorbentLayer* l = &layer;
        p.get_int_ = [l] { return static_cast<long long>(static_cast<int>(l->momentum_model_)); };
        p.set_int_ = [l](long long x) { l->momentum_model_ = static_cast<MomentumModel>(static_cast<int>(x)); };
      }

      auto& eqs = layer.equations_;
      if (auto* d = dynamic_cast<DispersionModelConstant*>(eqs.dispersion_model_.get())) {
        E("dispersion_coefficient", "dispersion coefficient", "m2/s", d->dispersion_coefficient_, 1e-8, 1e-2);
      } else if (auto* d = dynamic_cast<DispersionModelRuthven*>(eqs.dispersion_model_.get())) {
        E("molecular_diffusion_coefficient", "molecular diffusion coefficient", "m2/s", d->Dm, 1e-8, 1e-4);
      }
      if (auto* d = dynamic_cast<HeatDispersionModelConstant*>(eqs.heat_dispersion_model_.get())) {
        // The effective axial thermal conductivity lambda of the energy balance
        // (HeatDispersionModelConstant::Update), so W/(m K), not a diffusivity.
        E("heat_dispersion_coefficient", "axial thermal conductivity", "W/(m-K)", d->dispersion_coefficient_,
          1e-3, 100.0);
      }
      if (auto* h = dynamic_cast<ConstantInternalHeatTransferModel*>(eqs.internal_heat_transfer_model_.get())) {
        E("internal_htc", "internal heat transfer coefficient", "W/m2-K", h->heat_transfer_coefficient_, 1.0, 1e4);
      } else if (auto* h = dynamic_cast<InternalHeatTransferModel_Leva*>(eqs.internal_heat_transfer_model_.get())) {
        // Nu = eta1 * Re^eta2 * exp(-exp_coeff * dp / D_t). Bounds span Leva's
        // canonical (0.813, 0.90) and Casas 2012's fitted (41.13, 0.32) values.
        E("leva_eta1", "Leva Re coefficient (eta1)", "-", h->eta1, 0.01, 100.0);
        E("leva_eta2", "Leva Re exponent (eta2)", "-", h->eta2, 0.1, 1.5);
        E("leva_exp_coeff", "Leva wall-decay coefficient (exp_coeff)", "-", h->exp_coeff, 0.0, 20.0);
      }
      if (auto* h = dynamic_cast<ConstantExternalHeatTransferModel*>(eqs.external_heat_transfer_model_.get())) {
        E("external_htc", "external heat transfer coefficient", "W/m2-K", h->heat_transfer_coefficient_, 1.0, 1e4);
      } else if (auto* h = dynamic_cast<ExternalHeatTransferModel_Holman*>(eqs.external_heat_transfer_model_.get())) {
        E("holman_c", "Holman coefficient (c)", "-", h->c, 0.001, 1.0);
        E("holman_m", "Holman Re exponent (m)", "-", h->m, 0.1, 1.0);
        E("holman_n", "Holman Pr exponent (n)", "-", h->n, 0.1, 1.0);
        E("external_velocity", "external velocity", "m/s", h->external_velocity, 0.0, 10.0);
      }

      // Adsorbents are listed once, after the first layer that uses them; the
      // display subcategory names that layer too ("<layer>, <adsorbent>").
      for (auto& [ads_id, mass_frac] : layer.composition_) {
        (void)mass_frac;
        if (!dir.Contains(ads_id) || !emitted_adsorbents.insert(ads_id).second) continue;
        AddAdsorbent(ads_id, dir.Get(ads_id), ln);
      }
    }

    void AddAdsorbent(AdsorbentID id, Adsorbent& ads, const std::string& layer_name) {
      const std::string ap = Bracket("adsorbent", std::to_string(id));
      const std::string aa = Bracket("adsorbent", ads.name_);
      const std::string an = layer_name.empty() ? ads.name_ : layer_name + ", " + ads.name_;
      auto E = [&](const std::string& leaf, const std::string& category, const std::string& sub,
                   const std::string& name, const std::string& unit, auto& value, double lo, double hi) {
        Emplace(ap + "." + leaf, aa + "." + leaf, category, sub, name, unit, value, lo, hi);
      };
      E("solid_density", "adsorbent", an, "solid density", "kg/m3", ads.solid_density_, 100.0, 5000.0);
      E("heat_capacity", "adsorbent", an, "heat capacity", "J/(kg-K)", ads.heat_capacity_, 50.0, 2000.0);
      E("thermal_conductivity", "adsorbent", an, "thermal conductivity", "W/(m-K)",
        ads.thermal_conductivity_, 0.01, 10.0);
      E("particle_void_fraction", "adsorbent", an, "particle void fraction", "-", ads.ep_, 0.0, 0.8);

      auto& source = ads.adsorption_source_;
      for (auto& [component, dH] : source.heat_of_adsorption_) {
        E("heat_of_adsorption[" + component + "]", "adsorbent", an, component + " heat of adsorption",
          "J/mol", dH, 0.0, 1e6);
      }

      for (auto& [component, mtm] : source.mass_transfer_models_) {
        if (!mtm || !mtm->mtc_model_) continue;
        auto* mtc = mtm->mtc_model_.get();
        const std::string sub = an + ", " + component;
        const std::string kp = "kinetics[" + component + "].";
        // The rate coefficient's units depend on the rate model: LDF
        // dq/dt = k (q* - q) has k in 1/s, LDF2 dq/dt = k (q* - q)|q* - q|
        // has k in 1/(s mol/kg) = kg/(mol s). An LDF2 k ~ k_LDF / |q* - q|
        // with |q* - q| typically 0.1-10 mol/kg, so its bounds are a decade
        // wider than the LDF ones on each side.
        const bool ldf2 = mtm->GetModelType() == MassTransferModelType::LDF2;
        const char* k_units = ldf2 ? "kg/(mol-s)" : "1/s";
        if (auto* m = dynamic_cast<MassTransferCoefficientConstant*>(mtc)) {
          E(kp + "k", "adsorbent", sub, "mass transfer coefficient", k_units, m->constant_mtc_,
            ldf2 ? 1e-5 : 1e-4, ldf2 ? 1e3 : 100.0);
        } else if (auto* m = dynamic_cast<MassTransferCoefficientArrhenius*>(mtc)) {
          E(kp + "pre_exp_factor", "adsorbent", sub, "MTC pre-exponential factor", k_units, m->pre_exp_factor_,
            ldf2 ? 1e-3 : 0.01, ldf2 ? 1e5 : 1e4);
          E(kp + "activation_energy", "adsorbent", sub, "MTC activation energy", "J/mol",
            m->activation_energy_, 0.0, 1e5);
        } else if (auto* m = dynamic_cast<MassTransferCoefficientMacropore*>(mtc)) {
          E(kp + "tortuosity", "adsorbent", sub, "tortuosity", "-", m->tortuosity_, 1.0, 10.0);
          E(kp + "macropore_void_fraction", "adsorbent", sub, "macropore void fraction", "-",
            m->macropore_void_fraction_, 0.01, 0.8);
        } else if (auto* m = dynamic_cast<MassTransferCoefficientMacroporeKnudsen*>(mtc)) {
          E(kp + "pore_radius", "adsorbent", sub, "pore radius", "m", m->pore_radius_, 1e-10, 1e-6);
          E(kp + "tortuosity", "adsorbent", sub, "tortuosity", "-", m->tortuosity_, 1.0, 10.0);
          E(kp + "macropore_void_fraction", "adsorbent", sub, "macropore void fraction", "-",
            m->macropore_void_fraction_, 0.01, 0.8);
        } else if (auto* m = dynamic_cast<MassTransferCoefficient_Haynes*>(mtc)) {
          E(kp + "tortuosity", "adsorbent", sub, "tortuosity", "-", m->tortuosity_, 1.0, 10.0);
          E(kp + "macropore_void_fraction", "adsorbent", sub, "macropore void fraction", "-",
            m->macropore_void_fraction_, 0.01, 0.8);
        } else if (auto* m = dynamic_cast<MassTransferCoefficientMacroporeCapacity*>(mtc)) {
          E(kp + "tortuosity", "adsorbent", sub, "tortuosity", "-", m->tortuosity_, 1.0, 10.0);
          E(kp + "macropore_void_fraction", "adsorbent", sub, "macropore void fraction", "-",
            m->macropore_void_fraction_, 0.01, 0.8);
          E(kp + "concentration_floor", "adsorbent", sub, "capacity-factor concentration floor", "mol/m3",
            m->concentration_floor_, 0.0, 10.0);
        }
      }

      AddIsotherms(source, ap, aa, an);
    }

    void AddIsotherms(AdsorptionSource& source, const std::string& ap, const std::string& aa, const std::string& an) {
      for (const auto& entry : source.mass_transfer_models_) {
        // A named reference, not a structured binding: the lambdas below capture
        // it, and Apple clang before 16 rejects lambdas that capture a binding.
        const std::string& component = entry.first;
        auto& isotherm = source.GetIsotherm(component);
        if (!isotherm) continue;
        const std::string sub = an + ", " + component;
        const std::string leaf_head = ".isotherm[" + component + "].";
        auto E = [&](const char* leaf, const char* unit, ExpressionParam& value, double lo, double hi) {
          Emplace(ap + leaf_head + leaf, aa + leaf_head + leaf, "isotherm", sub, leaf, unit, value, lo, hi);
        };
        // Finds this component's parameter block in an isotherm of type T.
        auto params = [&](auto* typed) -> decltype(&typed->params_.begin()->second) {
          if (!typed) return nullptr;
          auto it = typed->params_.find(component);
          return it == typed->params_.end() ? nullptr : &it->second;
        };
        auto langmuir_p = [&](auto* p) {
          if (!p) return;
          E("m1", "mol/kg", p->m1, 0.01, 50.0);
          E("b1", "1/Pa", p->b1, 1e-15, 1e-3);
          E("q1", "J/mol", p->q1, 0.0, 1e5);
        };
        auto langmuir_c = [&](auto* p) {
          if (!p) return;
          E("m1", "mol/kg", p->m1, 0.01, 50.0);
          E("b1", "m3/mol", p->b1, 1e-10, 1e-3);
          E("q1", "J/mol", p->q1, 0.0, 1e5);
        };
        auto dual_p = [&](auto* p) {
          if (!p) return;
          langmuir_p(p);
          E("m2", "mol/kg", p->m2, 0.01, 50.0);
          E("b2", "1/Pa", p->b2, 1e-15, 1e-3);
          E("q2", "J/mol", p->q2, 0.0, 1e5);
        };
        auto dual_c = [&](auto* p) {
          if (!p) return;
          langmuir_c(p);
          E("m2", "mol/kg", p->m2, 0.01, 50.0);
          E("b2", "m3/mol", p->b2, 1e-10, 1e-3);
          E("q2", "J/mol", p->q2, 0.0, 1e5);
        };
        auto sips2 = [&](auto* p) {
          if (!p) return;
          E("A1", "mol/kg", p->A1, 0.01, 50.0);
          E("A2", "-", p->A2, -5.0, 5.0);
          // b_i = B1 exp(B2/RT) multiplies the partial pressure in Pa (IsothermFunction::Sips2).
          E("B1", "1/Pa", p->B1, 1e-13, 1e-5);
          E("B2", "J/mol", p->B2, 0.0, 1e5);
          E("C1", "-", p->C1, 0.1, 5.0);
          E("C2", "-", p->C2, -5.0, 5.0);
          E("T_ref", "K", p->T_ref, 273.0, 473.0);
        };
        auto* iso = isotherm.get();
        switch (iso->GetIsothermType()) {
          case IsothermType::HENRY:
            if (auto* p = params(dynamic_cast<HenryIsotherm*>(iso))) {
              E("k", "-", p->k, 1e-15, 1e-3);
              E("Q", "J/mol", p->Q, 0.0, 1e5);
            }
            break;
          case IsothermType::SINGLE_SITE_LANGMUIR_P:
            langmuir_p(params(dynamic_cast<LangmuirIsotherm*>(iso)));
            break;
          case IsothermType::EXTENDED_SINGLE_SITE_LANGMUIR_P:
            langmuir_p(params(dynamic_cast<ExtendedLangmuirIsotherm*>(iso)));
            break;
          case IsothermType::SINGLE_SITE_LANGMUIR_C:
            langmuir_c(params(dynamic_cast<LangmuirIsotherm_Concentration*>(iso)));
            break;
          case IsothermType::EXTENDED_SINGLE_SITE_LANGMUIR_C:
            langmuir_c(params(dynamic_cast<ExtendedLangmuirIsotherm_Concentration*>(iso)));
            break;
          case IsothermType::DUAL_SITE_LANGMUIR_P:
            dual_p(params(dynamic_cast<DualSiteLangmuirIsotherm*>(iso)));
            break;
          case IsothermType::EXTENDED_DUAL_SITE_LANGMUIR_P:
            dual_p(params(dynamic_cast<ExtendedDualSiteLangmuirIsotherm*>(iso)));
            break;
          case IsothermType::DUAL_SITE_LANGMUIR_C:
            dual_c(params(dynamic_cast<DualSiteLangmuirIsotherm_Concentration*>(iso)));
            break;
          case IsothermType::EXTENDED_DUAL_SITE_LANGMUIR_C:
            dual_c(params(dynamic_cast<ExtendedDualSiteLangmuirIsotherm_Concentration*>(iso)));
            break;
          case IsothermType::ANTI_LANGMUIR:
            if (auto* p = params(dynamic_cast<AntiLangmuirIsotherm*>(iso))) {
              E("a", "-", p->a, 1e-10, 1.0);
              E("b", "-", p->b, 1e-10, 1.0);
            }
            break;
          case IsothermType::FREUNDLICH:
            if (auto* p = params(dynamic_cast<FreundlichIsotherm*>(iso))) {
              E("a", "-", p->a, 0.01, 100.0);
              E("v", "-", p->v, 0.1, 10.0);
            }
            break;
          case IsothermType::SIPS:
            if (auto* p = params(dynamic_cast<SipsIsotherm*>(iso))) {
              E("q_sat", "mol/kg", p->q_sat, 0.01, 50.0);
              E("b", "1/Pa", p->b, 1e-10, 1.0);
              E("v", "-", p->v, 0.1, 5.0);
            }
            break;
          case IsothermType::SIPS_2:
            sips2(params(dynamic_cast<Sips2Isotherm*>(iso)));
            break;
          case IsothermType::EXTENDED_SIPS_2:
            sips2(params(dynamic_cast<ExtendedSips2Isotherm*>(iso)));
            break;
          case IsothermType::DUAL_SITE_SIPS:
            if (auto* p = params(dynamic_cast<DualSiteSipsIsotherm*>(iso))) {
              E("q_sat1", "mol/kg", p->q_sat1, 0.01, 50.0);
              E("b1", "1/Pa", p->b1, 1e-10, 1.0);
              E("v1", "-", p->v1, 0.1, 5.0);
              E("q_sat2", "mol/kg", p->q_sat2, 0.01, 50.0);
              E("b2", "1/Pa", p->b2, 1e-10, 1.0);
              E("v2", "-", p->v2, 0.1, 5.0);
            }
            break;
          case IsothermType::REDLICH_PETERSON:
            if (auto* p = params(dynamic_cast<RedlichPetersonIsotherm*>(iso))) {
              E("a", "-", p->a, 0.01, 100.0);
              E("b", "1/Pa", p->b, 1e-10, 1.0);
              E("v", "-", p->v, 0.1, 5.0);
            }
            break;
          case IsothermType::TOTH:
            if (auto* p = params(dynamic_cast<TothIsotherm*>(iso))) {
              E("q_sat", "mol/kg", p->q_sat, 0.01, 50.0);
              E("b", "1/Pa", p->b, 1e-10, 1.0);
              E("v", "-", p->v, 0.1, 5.0);
            }
            break;
          case IsothermType::TOTH_2:
            if (auto* p = params(dynamic_cast<Toth2Isotherm*>(iso))) {
              E("q_sat", "mol/kg", p->q_sat, 0.01, 50.0);
              E("x", "-", p->x, -1.0, 1.0);
              E("T_0", "K", p->T_0, 273.0, 473.0);
              E("b0", "1/Pa", p->b0, 1e-6, 1.0);
              E("dH0", "J/mol", p->dH0, 0.0, 1e6);
              E("w0", "-", p->w0, 0.01, 5.0);
              E("alpha", "-", p->alpha, 0.0, 2.0);
            }
            break;
          case IsothermType::ARRHENIUS_SIPS:
            if (auto* p = params(dynamic_cast<ArrheniusSipsIsotherm*>(iso))) {
              E("n_inf_a", "mol/kg", p->n_inf_a, 0.01, 50.0);
              E("n_inf_b", "J/mol", p->n_inf_b, 0.0, 1e5);
              E("k_A", "1/Pa", p->k_A, 1e-15, 1e-3);
              E("k_B", "J/mol", p->k_B, 0.0, 1e5);
              E("c_ref", "-", p->c_ref, 0.1, 5.0);
              E("c_amp", "-", p->c_amp, -1.0, 1.0);
              E("c_beta", "1/K", p->c_beta, 0.0, 1.0);
              E("c_Tref", "K", p->c_Tref, 273.0, 473.0);
            }
            break;
          case IsothermType::INERT:
          default:
            break;
        }
      }
    }

    // ---- user variables --------------------------------------------------
    void AddUserVariables(UserVariableRegistry& registry) {
      for (auto& var : registry.GetAll()) {
        double lo, hi;
        if (var.value > 0.0)      { lo = var.value * 0.1;  hi = var.value * 10.0; }
        else if (var.value < 0.0) { lo = var.value * 10.0; hi = var.value * 0.1; }
        else                      { lo = -1.0;             hi = 1.0; }
        const std::string path = "variable[" + var.name + "]";
        Emplace(path, path, "User Parameters", "Parameters", var.name, UnitCategorySymbol(var.unit_category),
                var.value, lo, hi);
        if (var.unit_category == UnitCategory::Custom) out.back().info_.unit_checked = false;
        // A variable defined by an expression is recomputed from it, so a
        // number written into var.value would not stick.
        if (!var.value_expr.empty()) out.back().info_.tunable = false;
      }
    }
  };

  // -------------------------------------------------------------------------
  // ParameterRegistry
  // -------------------------------------------------------------------------

  ParameterRegistry::ParameterRegistry(UserData& data) : data_(data) { Refresh(); }

  void ParameterRegistry::Refresh() {
    parameters_.clear();
    ParameterBuilder b{parameters_};
    b.AddStatic(data_);
    b.bed_ = &data_.bed;
    b.AddCycle(data_.cycle);
    b.AddBed(data_.bed, data_.adsorbent_directory);
    b.AddInitialConditions(data_.ic);
    b.AddFluid(data_.fluid);
    b.AddUserVariables(data_.variable_registry);

    // Paths are unique by construction once every owner has an ID; this guard
    // keeps Find well defined regardless (a repeat gets a " #n" suffix).
    std::unordered_set<std::string> used;
    for (auto& p : parameters_) {
      if (used.insert(p.info_.path).second) continue;
      const std::string base = p.info_.path;
      for (int n = 2; !used.insert(p.info_.path = base + " #" + std::to_string(n)).second; ++n) {}
    }
  }

  Parameter* ParameterRegistry::Find(std::string_view path) {
    for (auto& p : parameters_)
      if (p.info_.path == path) return &p;
    Parameter* match = nullptr;
    for (auto& p : parameters_) {
      if (p.info_.alias != path) continue;
      if (match)
        throw std::invalid_argument("parameter name '" + std::string(path) +
                                    "' is ambiguous (two objects share a name); use its ID path, e.g. " +
                                    match->info_.path + " or " + p.info_.path);
      match = &p;
    }
    return match;
  }

  const Parameter* ParameterRegistry::Find(std::string_view path) const {
    return const_cast<ParameterRegistry*>(this)->Find(path);
  }

  Parameter& ParameterRegistry::Get(std::string_view path) {
    if (Parameter* p = Find(path)) return *p;
    throw std::out_of_range("no parameter '" + std::string(path) + "' in this deck");
  }

  std::vector<Parameter*> ParameterRegistry::Tunable() {
    std::vector<Parameter*> out;
    for (auto& p : parameters_)
      if (p.info_.tunable) out.push_back(&p);
    return out;
  }

}  // namespace adsim
