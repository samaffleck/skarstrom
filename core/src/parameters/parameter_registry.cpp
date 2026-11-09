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

#include "adsim/equations/reaction_kinds.h"
#include "adsim/math/unit_parser.h"
#include "adsim/equations/isotherm_kinds.h"
#include "adsim/objects/deck_edit.h"
#include "adsim/objects/model_kinds.h"
#include "adsim/objects/reaction.h"
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
      // Bed side of every valve stream's linear valve law (ida_settings.h).
      IntField("settings", "valve_bed_pressure", "valve bed-side pressure",
               [](const UserData& u) { return (long long)static_cast<int>(u.settings.valve_bed_pressure_); },
               [](UserData& u, long long x) { u.settings.valve_bed_pressure_ = static_cast<ValveBedPressure>((int)x); },
               "face,end_cell"),
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
               "none,aitken,anderson,newton_krylov"),
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
      BoolField("settings", "warn_channel_reynolds", "warn when the channel Reynolds number of a structured layer passes 2000",
                [](const UserData& u) { return u.settings.warnings_.channel_reynolds; },
                [](UserData& u, bool b) { u.settings.warnings_.channel_reynolds = b; }),
      BoolField("settings", "warn_outlet_backflow",
                "warn when gas flows back into the bed through an unconnected outlet",
                [](const UserData& u) { return u.settings.warnings_.outlet_backflow; },
                [](UserData& u, bool b) { u.settings.warnings_.outlet_backflow = b; }),
      // Run health, after every cycle (run_health.h).
      BoolField("settings", "warn_mole_fraction_sum", "warn when the gas mole fractions do not sum to 1 (MOLE_FRACTION_SUM)",
                [](const UserData& u) { return u.settings.warnings_.mole_fraction_sum; },
                [](UserData& u, bool b) { u.settings.warnings_.mole_fraction_sum = b; }),
      BoolField("settings", "warn_mass_balance_stuck", "warn when a cycle mole balance stays open as the run settles (MASS_BALANCE_STUCK)",
                [](const UserData& u) { return u.settings.warnings_.mass_balance_stuck; },
                [](UserData& u, bool b) { u.settings.warnings_.mass_balance_stuck = b; }),
      BoolField("settings", "warn_step_retries", "warn on many rejected integrator attempts in a cycle (STEP_RETRIES)",
                [](const UserData& u) { return u.settings.warnings_.step_retries; },
                [](UserData& u, bool b) { u.settings.warnings_.step_retries = b; }),

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
      ExprField("css", "state_tolerance", "-", "state change tolerance",
                [](UserData& u) -> ExpressionParam& { return u.settings.css.state_tolerance_; }),
      IntField("css", "consecutive_cycles", "consecutive cycles",
               [](const UserData& u) { return (long long)u.settings.css.consecutive_cycles_; },
               [](UserData& u, long long x) { u.settings.css.consecutive_cycles_ = (int)x; }),

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
          // Choices are compared normalised the same way (a component name
          // such as "CarbonDioxide" is a choice of kpi.reaction.key_reactant).
          auto normalised = [](const std::string& name) {
            std::string n;
            for (char c : name) {
              if (std::isspace(static_cast<unsigned char>(c))) continue;
              n += (c == '-') ? '_' : static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            }
            return n;
          };
          for (size_t i = 0; i < info_.choices.size(); ++i) {
            if (key == normalised(info_.choices[i])) {
              set_int_(static_cast<long long>(i));
              return;
            }
          }
          for (const auto& [old_name, value] : info_.choice_aliases) {
            if (key == old_name) {
              set_int_(static_cast<long long>(value));
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

    // One entry per parameter of a live model, as its kind's descriptor
    // (objects/model_kinds.h) gives it: <path_head><name> and
    // <alias_head><name>, with the descriptor's label, unit, bounds and
    // tunable flag (`adjust` may rewrite the descriptor parameter first).
    void EmplaceKind(const KindDescriptor* kind, const ModelFieldList& fields, const std::string& path_head,
                     const std::string& alias_head, const std::string& category, const std::string& sub,
                     KindParameter (*adjust)(const KindParameter&) = nullptr) {
      if (!kind) return;
      for (const auto& [name, param] : fields) {
        const KindParameter* d = kind->Find(name);
        if (!d || !param) continue;
        const KindParameter p = adjust ? adjust(*d) : *d;
        Emplace(path_head + name, alias_head + name, category, sub, p.Label(), p.unit, *param, p.lower, p.upper);
        if (!p.tunable) out.back().info_.tunable = false;
      }
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
    UserData* data_ = nullptr;  // the deck, for the leaves that set through a deck-edit cascade

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
      AddHeatingMedium(step.heating_medium_, sp, sa, step.step_name_);
    }

    // The heating medium (jacket) on the wall's outer face during a step:
    // step[ID].medium.enabled, .heat_transfer_coefficient, the medium
    // temperature, a constant (.temperature) or a ramp (.temperature_linear_*,
    // ...) like an inlet temperature, and the internal heat exchanger's
    // .exchanger_area_per_volume and .exchanger_heat_transfer_coefficient. Listed for every step, enabled or not,
    // so a sweep or the CLI can switch a medium on.
    void AddHeatingMedium(HeatingMedium& medium, const std::string& sp, const std::string& sa,
                          const std::string& step_name) {
      const EndpointScope ms{sp + ".medium", sa + ".medium", step_name + ", heating medium"};
      EmplaceBool(ms.path + ".enabled", ms.alias + ".enabled", "cycle", ms.subcategory, "heating medium enabled",
                  medium.enabled_);
      const std::string name = "heating medium temperature";
      if (medium.temperature_profile_ && !dynamic_cast<ConstantBoundaryValue*>(medium.temperature_profile_.get())) {
        AddProfile(medium.temperature_profile_.get(), ms, "temperature", name, "K", 200.0, 773.0);
      } else {
        Emplace(ms.path + ".temperature", ms.alias + ".temperature", "cycle", ms.subcategory, name, "K",
                medium.temperature_, 200.0, 773.0);
      }
      Emplace(ms.path + ".heat_transfer_coefficient", ms.alias + ".heat_transfer_coefficient", "cycle",
              ms.subcategory, "heating medium heat transfer coefficient (wall outer face)", "W/m2-K",
              medium.heat_transfer_coefficient_, 0.0, 5000.0);
      // Internal heat exchanger (medium -> bed), off while the area is 0.
      Emplace(ms.path + ".exchanger_area_per_volume", ms.alias + ".exchanger_area_per_volume", "cycle",
              ms.subcategory, "internal heat exchanger area per bed volume", "1/m",
              medium.exchanger_area_per_volume_, 0.0, 500.0);
      Emplace(ms.path + ".exchanger_heat_transfer_coefficient", ms.alias + ".exchanger_heat_transfer_coefficient",
              "cycle", ms.subcategory, "internal heat exchanger heat transfer coefficient", "W/m2-K",
              medium.exchanger_heat_transfer_coefficient_, 0.0, 5000.0);
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
        EmplaceKind(FindKind(ModelFamily::Stream, "valve"), {{"cv", &v->valve_.cv_}, {"xt", &v->valve_.xt_}},
                    ep.path + ".valve.", ep.alias + ".valve.", "cycle", ep.subcategory);
        EmplaceBool(ep.path + ".valve.check_valve", ep.alias + ".valve.check_valve", "cycle", ep.subcategory,
                    "check valve (no flow against the stream's direction; unconnected streams only)",
                    v->check_valve_);
      } else if (auto* p = dynamic_cast<PumpBoundaryStream*>(stream)) {
        AddPumpPower(p->power_.get(), ep);
        AddPump(p->pump_.get(), ep);
      }
      AddBufferTank(stream->buffer_tank_, ep);
    }

    // The shaft power that drives a pump stream: a constant, or a polynomial
    // of the boundary pressure P [Pa], a + b P + ... + f P^5.
    void AddPumpPower(IPumpPowerModel* power, const EndpointScope& ep) {
      if (!power) return;
      EmplaceKind(KindAt(ModelFamily::PumpPower, static_cast<size_t>(power->GetType())), ModelFields(power),
                  ep.path + ".", ep.alias + ".", "cycle", ep.subcategory);
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
      if (!value) return;
      const KindDescriptor* kind = KindAt(ModelFamily::Profile, static_cast<size_t>(value->GetType()));
      if (!kind) return;
      const bool constant = value->GetType() == BoundaryValueType::CONSTANT;
      for (const auto& [pname, param] : ModelFields(value)) {
        const KindParameter* d = kind->Find(pname);
        if (!d) continue;
        const bool quantity = d->unit == kProfileQuantityUnit;  // the profiled value: the caller's unit and bounds
        const std::string leaf = constant ? head : head + "_" + kind->tag + "_" + pname;
        Emplace(ep.path + "." + leaf, ep.alias + "." + leaf, "cycle", ep.subcategory,
                constant ? name : name + " (" + d->Label() + ")", quantity ? unit : d->unit, *param,
                quantity ? lo : d->lower, quantity ? hi : d->upper);
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
      const KindDescriptor* kind = KindAt(ModelFamily::Flow, static_cast<size_t>(flow->GetType()));
      if (!kind || kind->parameters.empty()) return;
      const KindParameter& q = kind->parameters.front();  // the flow quantity
      AddProfile(flow->profile_.get(), ep, q.name, q.Label(), q.unit, q.lower, q.upper);
    }

    void AddPump(IPump* pump, const EndpointScope& ep) {
      if (!pump || !pump->use_pump_) return;
      auto E = [&](const std::string& leaf, const std::string& name, const std::string& unit, ExpressionParam& v,
                   double lo, double hi) {
        Emplace(ep.path + ".pump." + leaf, ep.alias + ".pump." + leaf, "cycle", ep.subcategory, name, unit, v, lo, hi);
      };
      E("polytropic_index", "pump polytropic index", "-", pump->poly_index_, 1.0, 1.7);
      if (IPumpEfficiencyModel* eff = pump->efficiency_model_.get())
        EmplaceKind(KindAt(ModelFamily::PumpEfficiency, static_cast<size_t>(eff->GetType())), ModelFields(eff),
                    ep.path + ".pump.", ep.alias + ".pump.", "cycle", ep.subcategory);
      // Pump kinds: index 0 is "none", PumpType k is k + 1.
      EmplaceKind(KindAt(ModelFamily::Pump, static_cast<size_t>(pump->GetType()) + 1), ModelFields(pump),
                  ep.path + ".pump.", ep.alias + ".pump.", "cycle", ep.subcategory);
    }

    // A pressure (pressure boundary or valve) profile starts from the bed
    // pressure, so its start is not listed; a ramp lasts the step unless
    // <kind>_use_step_time is off, when <kind>_time is listed too (it must not
    // exceed the step time; Cycle::PerformChecks).
    void AddPressureProfile(IBoundaryValue* value, const EndpointScope& ep) {
      if (!value) return;
      if (auto* bv = dynamic_cast<ConstantBoundaryValue*>(value)) {
        Emplace(ep.path + ".pressure", ep.alias + ".pressure", "cycle", ep.subcategory,
                "pressure (constant)", "Pa", bv->constant_value_, 1e4, 1e6);
        return;
      }
      const KindDescriptor* kind = KindAt(ModelFamily::Profile, static_cast<size_t>(value->GetType()));
      if (!kind) return;
      bool* use_step_time = nullptr;
      if (auto* v = dynamic_cast<LinearBoundaryValue*>(value)) use_step_time = &v->use_step_time_;
      if (auto* v = dynamic_cast<ExponentialBoundaryValue*>(value)) use_step_time = &v->use_step_time_;
      if (auto* v = dynamic_cast<LogLinearBoundaryValue*>(value)) use_step_time = &v->use_step_time_;
      std::string short_label = kind->label;
      for (auto& c : short_label) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
      for (const auto& [pname, param] : ModelFields(value)) {
        if (pname == "start") continue;  // a pressure ramp starts from the bed pressure
        const KindParameter* d = kind->Find(pname);
        if (!d) continue;
        const std::string head = ep.path + ".pressure_" + kind->tag + "_";
        const std::string ahead = ep.alias + ".pressure_" + kind->tag + "_";
        if (pname == "time") {
          if (!use_step_time) continue;
          EmplaceBool(head + "use_step_time", ahead + "use_step_time", "cycle", ep.subcategory,
                      "pressure (" + short_label + " ramp over the step time)", *use_step_time);
          if (*use_step_time) continue;
        }
        const bool quantity = d->unit == kProfileQuantityUnit;
        Emplace(head + pname, ahead + pname, "cycle", ep.subcategory, "pressure (" + d->Label() + ")",
                quantity ? "Pa" : d->unit, *param, quantity ? 1e4 : d->lower, quantity ? 1e6 : d->upper);
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
      if (IParticleHeatTransferModel* h = bed.equations_.particle_heat_transfer_model_.get())
        EmplaceKind(KindAt(ModelFamily::ParticleHeatTransfer, static_cast<size_t>(h->GetModelType())),
                    ModelFields(h), "bed.", "bed.", "bed", "thermal");
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
      // 0 = isosteric heat of adsorption, 1 = internal-energy convention,
      // 2 = isosteric heat from -dU + R*T at the local temperature.
      EmplaceEnum("fluid.heat_of_adsorption_convention", "", "fluid", "fluid", "heat of adsorption convention",
                  fluid.heat_of_adsorption_convention_, "isosteric,internal_energy,isosteric_rt");
      // 0 = ideal-gas Cp(T) polynomials (held CoolProp value for a component
      // without one), 1 = constant (fluid.gas_heat_capacity). "coolprop" is
      // the old name of 0.
      EmplaceEnum("fluid.gas_heat_capacity_mode", "", "fluid", "fluid", "gas heat capacity mode",
                  fluid.gas_heat_capacity_mode_, "ideal_gas,constant");
      out.back().info_.choice_aliases = {{"coolprop", 0}};
      // 0 = pure-component polynomials, 1 = constant mixture value (fluid.viscosity).
      EmplaceEnum("fluid.viscosity_mode", "", "fluid", "fluid", "gas viscosity mode",
                  fluid.viscosity_mode_, "polynomial,constant");
      Emplace("fluid.viscosity", "", "fluid", "fluid", "constant gas viscosity", "Pa.s", fluid.viscosity_, 1e-6, 1e-4);
      for (const auto& component : fluid.components_) {
        auto it = fluid.adsorbed_heat_capacity_.find(component);
        if (it == fluid.adsorbed_heat_capacity_.end()) continue;
        Emplace("fluid.adsorbed_heat_capacity[" + component + "]", "", "fluid", "fluid",
                component + " adsorbed-phase heat capacity", "J/(mol-K)", it->second, 0.0, 500.0);
      }
      for (const auto& component : fluid.components_) {
        auto it = fluid.gas_heat_capacity_.find(component);
        if (it == fluid.gas_heat_capacity_.end()) continue;
        Emplace("fluid.gas_heat_capacity[" + component + "]", "", "fluid", "fluid",
                component + " constant gas heat capacity", "J/(mol-K)", it->second, 10.0, 500.0);
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
      if (!layer.IsStructured()) {
        // A packed bed.
        E("bed_void_fraction", "bed void fraction", "-", layer.eb_, 0.2, 0.8);
        E("particle_diameter", "particle diameter", "m", layer.particle_diameter_, 1e-5, 0.01);
        // 0 = Darcy (Kozeny-Carman), 1 = Ergun, 2 = Darcy (Blake-Kozeny).
        auto& p = Add(Info(lp + ".momentum_model", la + ".momentum_model", "layer", ln, "momentum model", "-",
                           ParameterKind::Integer, std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::quiet_NaN()));
        p.info_.choices = KindTags(ModelFamily::Momentum);
        MixedAdsorbentLayer* l = &layer;
        p.get_int_ = [l] { return static_cast<long long>(static_cast<int>(l->momentum_model_)); };
        p.set_int_ = [l](long long x) { l->momentum_model_ = static_cast<MomentumModel>(static_cast<int>(x)); };
        // Gas-particle film correlation (film_correlation.h): the Nu of the
        // film particle heat transfer and the Sh of the Haynes k_f.
        auto& f = Add(Info(lp + ".film_correlation", la + ".film_correlation", "layer", ln,
                           "film correlation (particle Nu and Sh)", "-", ParameterKind::Integer,
                           std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));
        f.info_.choices = FilmCorrelationKinds();
        f.get_int_ = [l] { return static_cast<long long>(static_cast<int>(l->film_correlation_)); };
        f.set_int_ = [l](long long x) { l->film_correlation_ = static_cast<FilmCorrelation>(static_cast<int>(x)); };
        // Particle shape (ParticleShape): the LDF diffusion geometry; for a
        // cylinder the particle diameter is the extrudate diameter 2 r_f.
        auto& s = Add(Info(lp + ".particle_shape", la + ".particle_shape", "layer", ln,
                           "particle shape (diffusion geometry)", "-", ParameterKind::Integer,
                           std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));
        s.info_.choices = ParticleShapeKinds();
        s.get_int_ = [l] { return static_cast<long long>(static_cast<int>(l->particle_shape_)); };
        s.set_int_ = [l](long long x) { l->particle_shape_ = static_cast<ParticleShape>(static_cast<int>(x)); };
      } else {
        // A structured layer: only the fields of its kind. Bed void fraction
        // and momentum model are not used (the open frontal area and laminar
        // channel flow replace them). The particle diameter stays listed (its
        // path is kept for existing sweeps and front ends), although no
        // closure a structured layer accepts reads it: the deck checks refuse
        // the packed-bed correlations that did.
        E("particle_diameter", "particle diameter (not used by a structured layer; kept for switching back to packed)",
          "m", layer.particle_diameter_, 1e-5, 0.01);
        auto& g = layer.structured_;
        if (layer.kind_ == LayerKind::Monolith) {
          // 0 = square, 1 = circular.
          auto& p = Add(Info(lp + ".channel_shape", la + ".channel_shape", "layer", ln, "channel shape", "-",
                             ParameterKind::Integer, std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::quiet_NaN()));
          p.info_.choices = {"square", "circular"};
          StructuredGeometry* sg = &g;
          p.get_int_ = [sg] { return static_cast<long long>(static_cast<int>(sg->channel_shape_)); };
          p.set_int_ = [sg](long long x) { sg->channel_shape_ = static_cast<ChannelShape>(static_cast<int>(x)); };
          // 100 to 1300 cells per square inch.
          E("cell_density", "cell density", "1/m2", g.cell_density_, 1.55e5, 2.0e6);
          E("wall_thickness", "wall thickness", "m", g.wall_thickness_, 5e-5, 1e-3);
        } else {
          E("sheet_thickness", "sheet thickness", "m", g.sheet_thickness_, 5e-5, 2e-3);
          E("gap", "gap between sheets", "m", g.gap_, 5e-5, 2e-3);
        }
        E("washcoat_thickness", "washcoat thickness", "m", g.washcoat_thickness_, 0.0, 2e-4);
        // Hawthorn's developing-flow film coefficients over the layer length.
        EmplaceBool(lp + ".entrance_effects", la + ".entrance_effects", "layer", ln,
                    "entrance effects on the film coefficients", g.entrance_effects_);
      }

      auto& eqs = layer.equations_;
      auto models = [&](ModelFamily family, auto* model) {
        if (model)
          EmplaceKind(KindAt(family, static_cast<size_t>(model->GetModelType())), ModelFields(model), lp + ".",
                      la + ".", "layer", ln);
      };
      models(ModelFamily::Dispersion, eqs.dispersion_model_.get());
      models(ModelFamily::HeatDispersion, eqs.heat_dispersion_model_.get());
      models(ModelFamily::InternalHeatTransfer, eqs.internal_heat_transfer_model_.get());
      models(ModelFamily::ExternalHeatTransfer, eqs.external_heat_transfer_model_.get());

      // Adsorbents are listed once, after the first layer that uses them; the
      // display subcategory names that layer too ("<layer>, <adsorbent>").
      for (auto& [ads_id, mass_frac] : layer.composition_) {
        (void)mass_frac;
        if (!dir.Contains(ads_id) || !emitted_adsorbents.insert(ads_id).second) continue;
        AddAdsorbent(ads_id, dir.Get(ads_id), ln);
      }
    }

    // One coefficient model of a component (fast site, or the slow site of a
    // dual-rate model, whose leaves take the prefix slow_). The generic model
    // (series_resistance) lists, from its kind descriptor: a Boolean switch
    // per resistance (resistances.<tag>; knudsen only with macropore on), the
    // capacity-factor mode while a gas-side term is on, then the parameters of
    // each enabled term and of the mode under their flat names. The switches
    // and the mode set through the deck-edit cascades (SetMassTransferResistance,
    // SetCapacityFactorMode), never by a raw field write. An older coefficient
    // class (built by hand) lists its own kind's parameters.
    void AddCoefficient(IMassTransferCoefficientModel* mtc, AdsorbentID id, const std::string& component,
                        MassTransferSite site, bool ldf2, const std::string& head, const std::string& ahead,
                        const std::string& sub) {
      const bool slow = site == MassTransferSite::Slow;
      const std::string pre = slow ? kSlowSitePrefix : "";
      auto* g = dynamic_cast<MassTransferCoefficientSeriesResistance*>(mtc);
      if (!g) {
        const KindDescriptor* kind = KindAt(ModelFamily::RateCoefficient, static_cast<size_t>(mtc->GetModelType()));
        ModelFieldList fields = ModelFields(mtc);
        if (slow) {
          for (auto& f : fields) f.first = pre + f.first;
        }
        static const std::vector<KindDescriptor> renamed = [] {
          std::vector<KindDescriptor> out = KindDescriptors(ModelFamily::RateCoefficient);
          for (auto& k : out)
            for (auto& p : k.parameters) p = SlowSiteParameter(p);
          return out;
        }();
        if (slow && kind) kind = &renamed[static_cast<size_t>(mtc->GetModelType())];
        EmplaceKind(kind, fields, head, ahead, "adsorbent", sub, ldf2 ? &Ldf2RateConstant : nullptr);
        return;
      }
      const KindDescriptor& kind = ResistanceKind();
      const ModelFieldList all = ResistanceFields(g);
      auto emit_params = [&](const std::vector<std::string>& names) {
        for (const std::string& name : names) {
          const KindParameter* d = kind.Find(name);
          if (!d) continue;
          for (const auto& [field, param] : all) {
            if (field != name || !param) continue;
            KindParameter p = ldf2 ? Ldf2RateConstant(*d) : *d;
            if (slow) p = SlowSiteParameter(p);
            Emplace(head + pre + name, ahead + pre + name, "adsorbent", sub, p.Label(), p.unit, *param, p.lower,
                    p.upper);
            if (!p.tunable) out.back().info_.tunable = false;
          }
        }
      };
      UserData* data = data_;
      const std::string label_pre = slow ? "slow-site " : "";
      for (const ResistanceDescriptor& r : kind.resistances) {
        if (!r.parent.empty() && !g->Enabled(r.parent)) continue;
        auto& p = Add(Info(head + pre + "resistances." + r.tag, ahead + pre + "resistances." + r.tag, "adsorbent", sub,
                           label_pre + r.label + " resistance", "-", ParameterKind::Boolean,
                           std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));
        const std::string tag = r.tag;
        p.get_bool_ = [g, tag] { return g->Enabled(tag); };
        p.set_bool_ = [data, id, component, site, tag](bool on) {
          std::string err;
          if (!data || !SetMassTransferResistance(*data, id, component, site, tag, on, err))
            throw std::invalid_argument(err.empty() ? "cannot set the " + tag + " resistance" : err);
        };
        if (g->Enabled(r.tag)) emit_params(r.parameters);
      }
      if (g->UsesCapacityFactor()) {
        auto& p = Add(Info(head + pre + "capacity_factor_mode", ahead + pre + "capacity_factor_mode", "adsorbent", sub,
                           label_pre + "capacity factor mode", "-", ParameterKind::Integer,
                           std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));
        for (const auto& m : kind.capacity_modes) p.info_.choices.push_back(m.tag);
        p.get_int_ = [g] { return static_cast<long long>(g->CapacityMode()); };
        p.set_int_ = [data, id, component, site](long long x) {
          std::string err;
          if (!data || !SetCapacityFactorMode(*data, id, component, site, static_cast<CapacityFactorMode>(x), err))
            throw std::invalid_argument(err.empty() ? "cannot set the capacity-factor mode" : err);
        };
      }
      emit_params(g->CapacityParameters());
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
      // 0 = none (each isotherm on its own), 1 = IAST over the pure isotherms.
      EmplaceEnum(ap + ".mixing_rule", aa + ".mixing_rule", "adsorbent", an, "mixing rule",
                  source.mixing_rule_, "none,iast");
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
        AddCoefficient(mtc, id, component, MassTransferSite::Fast, ldf2, ap + "." + kp, aa + "." + kp, sub);
        // The rate model's own parameters (dual_rate_ldf: split) and the slow
        // site's coefficient model (dual_rate_ldf: slow_<name>).
        EmplaceKind(KindAt(ModelFamily::RateModel, static_cast<size_t>(mtm->GetModelType())),
                    ModelFields(mtm.get()), ap + "." + kp, aa + "." + kp, "adsorbent", sub);
        if (IMassTransferCoefficientModel* slow = mtm->SlowCoefficientModel())
          AddCoefficient(slow, id, component, MassTransferSite::Slow, ldf2, ap + "." + kp, aa + "." + kp, sub);
      }

      AddIsotherms(source, ap, aa, an);
      AddReactions(ads.reaction_source_, ap, aa, an);
    }

    // Reactions on a solid: adsorbent[A].reaction[R].*, category "reaction",
    // subcategory "<layer>, <solid>, <reaction>". The rate constants' units
    // depend on the orders and the basis, so they carry their literature unit
    // text unchecked (unit_checked = false). Pre-exponential factors,
    // activation energies, orders and the effectiveness factor are tunable.
    void AddReactions(ReactionSource& source, const std::string& ap, const std::string& aa, const std::string& an) {
      const double nan = std::numeric_limits<double>::quiet_NaN();
      for (auto& r : source.reactions_) {
        const std::string rp = ap + ".reaction[" + std::to_string(r.id_) + "].";
        const std::string ra = aa + ".reaction[" + r.name_ + "].";
        const std::string sub = an + ", " + r.name_;
        const bool pp = r.basis_ == RateBasis::PartialPressure;
        const std::string act = pp ? "bar" : "(mol/m3)";
        auto power = [&](double n) { return act + "^" + ExpressionParam::FormatScalar(n); };
        // Emplace an expression and set its tunable flag and unit check.
        auto E = [&](const std::string& leaf, const std::string& name, const std::string& unit,
                     ExpressionParam& value, double lo, double hi, bool tunable, bool checked = true) {
          Emplace(rp + leaf, ra + leaf, "reaction", sub, name, unit, value, lo, hi);
          out.back().info_.tunable = tunable && !value.HasExpression();
          out.back().info_.unit_checked = checked;
        };

        for (auto& entry : r.nu_) {
          const std::string leaf = "nu[" + entry.first + "]";
          Emplace(rp + leaf, ra + leaf, "reaction", sub, entry.first + " stoichiometric coefficient", "-",
                  entry.second, nan, nan);
          out.back().info_.tunable = false;
        }
        EmplaceEnum(rp + "basis", ra + "basis", "reaction", sub, "rate basis", r.basis_,
                    "partial_pressure,concentration");

        if (r.rate_) {
          IReactionRate* rate = r.rate_.get();
          // Sum of the orders (the rate constant's activity power) and of the
          // stoichiometric coefficients (the equilibrium constant's).
          double order_sum = 0.0, nu_sum = 0.0;
          for (auto& [leaf, param] : ReactionRateFields(rate))
            if (leaf.rfind("order[", 0) == 0) order_sum += param->Get();
          for (const auto& [c, v] : r.nu_) nu_sum += v;
          const std::string k_unit = order_sum == 0.0 ? "mol/(kg-s)" : "mol/(kg-s-" + power(order_sum) + ")";
          if (auto* l = dynamic_cast<LhhwRate*>(rate)) {
            EmplaceBool(rp + "rate.reversible", ra + "rate.reversible", "reaction", sub, "reversible",
                        l->reversible_);
          }
          for (auto& [leaf, param] : ReactionRateFields(rate)) {
            const std::string path_leaf = "rate." + leaf;
            if (leaf == "pre_exp_factor") {
              E(path_leaf, "pre-exponential factor", k_unit, *param, 0.0, 1e16, true, false);
            } else if (leaf == "activation_energy") {
              E(path_leaf, "activation energy", "J/mol", *param, 0.0, 3e5, true);
            } else if (leaf.rfind("order[", 0) == 0) {
              E(path_leaf, leaf.substr(6, leaf.size() - 7) + " order", "-", *param, -3.0, 3.0, true);
            } else if (leaf == "equilibrium_pre_exp_factor") {
              E(path_leaf, "equilibrium pre-exponential factor", nu_sum == 0.0 ? "-" : power(nu_sum), *param,
                1e-30, 1e30, true, false);
            } else if (leaf == "equilibrium_energy") {
              E(path_leaf, "equilibrium energy", "J/mol", *param, -5e5, 5e5, false);
            } else if (leaf == "denominator_exponent") {
              E(path_leaf, "denominator exponent", "-", *param, 0.0, 4.0, false);
            } else if (leaf == "expression") {
              E(path_leaf, "rate expression", "mol/(kg-s)", *param, nan, nan, false, false);
            } else if (leaf.rfind("adsorption[", 0) == 0) {
              const size_t dot = leaf.find("].") + 2;
              const std::string term = leaf.substr(0, dot - 1);   // adsorption[j]
              const std::string field = leaf.substr(dot);
              if (field == "pre_exp_factor") {
                // The term's activity power is the sum of its exponents.
                double e_sum = 0.0;
                const size_t j = std::stoul(term.substr(11, term.size() - 12));
                if (auto* l = dynamic_cast<LhhwRate*>(rate); l && j >= 1 && j <= l->adsorption_.size())
                  for (const auto& [c, e] : l->adsorption_[j - 1].exponent_) e_sum += e.Get();
                E(path_leaf, term + " pre-exponential factor", e_sum == 0.0 ? "-" : "1/" + power(e_sum), *param,
                  0.0, 1e10, true, false);
              } else if (field == "energy") {
                E(path_leaf, term + " energy", "J/mol", *param, -3e5, 3e5, false);
              } else {
                const std::string c = field.substr(9, field.size() - 10);  // exponent[C]
                E(path_leaf, term + " " + c + " exponent", "-", *param, -3.0, 3.0, false);
              }
            } else {
              E(path_leaf, leaf, "-", *param, nan, nan, false, false);
            }
          }
        }

        E("heat_of_reaction", "heat of reaction", "J/mol", r.heat_of_reaction_, -1e6, 1e6, false);
        EmplaceEnum(rp + "heat_mode", ra + "heat_mode", "reaction", sub, "heat of reaction mode", r.heat_mode_,
                    "constant,kirchhoff");
        E("heat_reference_temperature", "heat of reaction reference temperature", "K",
          r.heat_reference_temperature_, 200.0, 1500.0, false);
        if (r.effectiveness_.mode_ == EffectivenessMode::Constant) {
          E("effectiveness", "effectiveness factor", "-", r.effectiveness_.effectiveness_, 0.0, 1.0, true);
        } else {
          E("tortuosity", "tortuosity", "-", r.effectiveness_.tortuosity_, 1.0, 10.0, false);
        }
      }
    }

    void AddIsotherms(AdsorptionSource& source, const std::string& ap, const std::string& aa, const std::string& an) {
      for (const auto& entry : source.mass_transfer_models_) {
        const std::string& component = entry.first;
        auto& isotherm = source.GetIsotherm(component);
        if (!isotherm) continue;
        const std::string head = ".isotherm[" + component + "].";
        EmplaceKind(&IsothermKind(isotherm->GetIsothermType()), IsothermFields(isotherm.get(), component),
                    ap + head, aa + head, "isotherm", an + ", " + component);
      }
    }

    // ---- KPI configuration, cost data, the sensitivity study ----------------
    // Settings, not model values: none is tunable.
    void EmplaceSetting(const std::string& path, const std::string& category, const std::string& sub,
                        const std::string& name, const std::string& unit, double& value, double lo, double hi) {
      Emplace(path, "", category, sub, name, unit, value, lo, hi);
      out.back().info_.tunable = false;
    }

    // Which components' KPIs are computed: kpi[Component].tracked. The
    // feed / product / recycle endpoints are structure (build::SetKPI).
    void AddKpi(KPIConfig& kpi) {
      for (auto& c : kpi.per_component) {
        if (c.product_component.empty()) continue;
        EmplaceBool("kpi[" + c.product_component + "].tracked", "", "kpi", c.product_component,
                    c.product_component + " KPIs tracked", c.tracked);
      }

      // Reaction KPIs: kpi.reaction.key_reactant (a choice: "none" or a component),
      // kpi.reaction.product[C] and kpi.reaction.basis.
      ReactionKPIConfig* rx = &kpi.reaction;
      std::vector<std::string> names;
      for (const auto& c : kpi.per_component)
        if (!c.product_component.empty()) names.push_back(c.product_component);
      {
        constexpr double nan = std::numeric_limits<double>::quiet_NaN();
        auto& p = Add(Info("kpi.reaction.key_reactant", "", "kpi", "reaction", "reaction KPIs key reactant", "-",
                           ParameterKind::Integer, nan, nan));
        p.info_.choices = {"none"};
        for (const auto& n : names) p.info_.choices.push_back(n);
        p.info_.tunable = false;
        p.get_int_ = [rx, names] {
          for (size_t i = 0; i < names.size(); ++i)
            if (names[i] == rx->key_reactant) return static_cast<long long>(i + 1);
          return 0LL;
        };
        p.set_int_ = [rx, names](long long x) {
          rx->key_reactant = (x >= 1 && x <= static_cast<long long>(names.size()))
                                 ? names[static_cast<size_t>(x - 1)] : std::string{};
        };
      }
      for (const auto& n : names) {
        constexpr double nan = std::numeric_limits<double>::quiet_NaN();
        auto& p = Add(Info("kpi.reaction.product[" + n + "]", "", "kpi", "reaction", n + " is a reaction product",
                           "-", ParameterKind::Boolean, nan, nan));
        p.info_.tunable = false;
        p.get_bool_ = [rx, n] { return rx->IsProduct(n); };
        p.set_bool_ = [rx, n](bool b) { rx->SetProduct(n, b); };
      }
      EmplaceEnum("kpi.reaction.basis", "", "kpi", "reaction", "space-time yield basis", rx->basis,
                  "catalyst_mass,bed_volume");
      out.back().info_.tunable = false;
    }

    // The cost estimate (solver/cost_estimator.h): cost.*. Money in the
    // deck's currency (unit "GBP", not checked).
    void AddCost(CostEstimator& cost, AdsorbentDirectory& dir) {
      auto money = [&](const std::string& path, const std::string& sub, const std::string& name, double& v,
                       double lo, double hi) {
        EmplaceSetting(path, "cost", sub, name, "GBP", v, lo, hi);
        out.back().info_.unit_checked = false;
      };
      auto count = [&](const std::string& path, const std::string& sub, const std::string& name, int& v) {
        EmplaceInt(path, "cost", sub, name, v);
      };
      count("cost.vessel_count", "plant", "vessels per train", cost.vessel_count);
      count("cost.train_count", "plant", "trains", cost.train_count);
      count("cost.operating_days_per_year", "plant", "operating days per year", cost.operating_days_per_year);
      EmplaceSetting("cost.plant_life_years", "cost", "finance", "plant life (years)", "-", cost.plant_life_years, 1.0,
                     50.0);
      EmplaceSetting("cost.discount_rate", "cost", "finance", "discount rate", "-", cost.discount_rate, 0.0, 0.3);
      money("cost.electricity_cost_per_kWh", "operating", "electricity cost per kWh", cost.electricity_cost_per_kWh,
            0.0, 1.0);
      EmplaceSetting("cost.maintenance_fraction", "cost", "operating", "maintenance (fraction of CAPEX per year)", "-",
                     cost.maintenance_fraction, 0.0, 0.2);
      EmplaceSetting("cost.overhead_fraction", "cost", "operating", "overhead (fraction of CAPEX per year)", "-",
                     cost.overhead_fraction, 0.0, 0.2);
      money("cost.labour_cost_per_year", "operating", "labour cost per year", cost.labour_cost_per_year, 0.0, 1e7);
      EmplaceSetting("cost.adsorbent_replacement_years", "cost", "adsorbent", "adsorbent replacement interval (years)",
                     "-", cost.adsorbent_replacement_years, 0.5, 50.0);
      std::vector<AdsorbentID> ids;
      for (const auto& [id, ads] : dir.GetAll()) ids.push_back(id);
      std::sort(ids.begin(), ids.end());
      for (AdsorbentID id : ids) {
        auto it = cost.adsorbent_cost_per_kg.find(id);
        if (it == cost.adsorbent_cost_per_kg.end()) continue;
        money("cost.adsorbent[" + std::to_string(id) + "].cost_per_kg", "adsorbent",
              dir.Get(id).name_ + " cost per kg", it->second, 0.0, 1e4);
      }
      EmplaceSetting("cost.vessel.installation_factor", "cost", "vessel", "vessel installation factor", "-",
                     cost.vessel_cost_.installation_factor, 1.0, 10.0);
      EmplaceBool("cost.heat_exchanger.include", "", "cost", "heat exchanger", "include a heat exchanger",
                  cost.heat_exchanger_cost_.include);
      money("cost.heat_exchanger.cost", "heat exchanger", "heat exchanger cost", cost.heat_exchanger_cost_.cost, 0.0,
            1e7);
      auto units = [&](auto& bank, const std::string& head, const std::string& what) {
        for (size_t i = 0; i < bank.units.size(); ++i) {
          auto& u = bank.units[i];
          const std::string up = "cost." + head + "[" + std::to_string(i + 1) + "]";
          const std::string sub = what + " " + std::to_string(i + 1);
          auto& lbl = Add(Info(up + ".label", "", "cost", sub, what + " label", "", ParameterKind::String,
                               std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));
          std::string* text = &u.label;
          lbl.get_str_ = [text] { return *text; };
          lbl.set_str_ = [text](const std::string& t) { *text = t; };
          EmplaceSetting(up + ".rated_power", "cost", sub, what + " rated shaft power", "kW", u.rated_power_kW, 0.0,
                         1e5);
          money(up + ".reference_cost", sub, what + " reference purchase cost", u.reference_cost, 0.0, 1e8);
          EmplaceSetting(up + ".reference_power", "cost", sub, what + " reference power", "kW", u.reference_power_kW,
                         0.0, 1e5);
          EmplaceSetting(up + ".cost_exponent", "cost", sub, what + " cost exponent", "-", u.cost_exponent, 0.0, 1.5);
          count(up + ".count", sub, what + " count", u.count);
        }
      };
      units(cost.compressors_, "compressor", "compressor");
      units(cost.vacuum_pumps_, "vacuum_pump", "vacuum pump");
    }

    // The settings of the deck's sensitivity study (solver/sensitivity_study.h):
    // sensitivity.*. Its inputs and outputs are lists edited as a whole.
    void AddSensitivity(sensitivity::Study& study) {
      const std::string c = "sensitivity";
      auto& seed = Add(Info("sensitivity.seed", "", c, c, "random seed", "-", ParameterKind::Integer,
                            std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));
      std::uint64_t* sp = &study.seed;
      seed.get_int_ = [sp] { return static_cast<long long>(*sp); };
      seed.set_int_ = [sp](long long x) { *sp = static_cast<std::uint64_t>(x); };
      EmplaceInt("sensitivity.trajectories", c, c, "Morris trajectories", study.trajectories);
      EmplaceInt("sensitivity.levels", c, c, "Morris grid levels", study.levels);
      EmplaceSetting("sensitivity.threshold", c, c, "screening threshold (fraction of the largest mu*)", "-",
                     study.threshold, 0.0, 1.0);
      EmplaceInt("sensitivity.runs", c, c, "Monte Carlo first batch", study.runs);
      EmplaceInt("sensitivity.batch", c, c, "Monte Carlo runs per added batch", study.batch);
      EmplaceInt("sensitivity.max_runs", c, c, "Monte Carlo run limit", study.max_runs);
      EmplaceSetting("sensitivity.rel_tol", c, c, "percentile tolerance (fraction of the 5-95 % band)", "-",
                     study.rel_tol, 0.0, 1.0);
      EmplaceInt("sensitivity.bootstrap", c, c, "bootstrap resamples", study.bootstrap);
      EmplaceBool("sensitivity.require_css", "", c, c, "runs must reach cyclic steady state", study.policy.require_css);
      EmplaceInt("sensitivity.cycles", c, c, "cycles per run (0: the deck's)", study.policy.cycles);
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
    b.data_ = &data_;
    b.AddCycle(data_.cycle);
    b.AddBed(data_.bed, data_.adsorbent_directory);
    b.AddInitialConditions(data_.ic);
    b.AddFluid(data_.fluid);
    b.AddUserVariables(data_.variable_registry);
    b.AddKpi(data_.kpi_config);
    b.AddCost(data_.cost_estimator, data_.adsorbent_directory);
    b.AddSensitivity(data_.sensitivity_study);

    // Paths are unique by construction once every owner has an ID; this guard
    // keeps Find well defined regardless (a repeat gets a " #n" suffix).
    std::unordered_set<std::string> used;
    for (auto& p : parameters_) {
      if (used.insert(p.info_.path).second) continue;
      const std::string base = p.info_.path;
      for (int n = 2; !used.insert(p.info_.path = base + " #" + std::to_string(n)).second; ++n) {}
    }
  }

  namespace {
    // The former arrhenius coefficient kind named the lumped coefficient
    // "pre_exp_factor": kinetics[C].pre_exp_factor and slow_pre_exp_factor
    // resolve to kinetics[C].k and slow_k (shown while the lumped term is on).
    // Empty when the path has no such leaf.
    std::string LumpedAlias(std::string_view path) {
      for (const char* old_leaf : {"].slow_pre_exp_factor", "].pre_exp_factor"}) {
        const std::string_view leaf(old_leaf);
        if (path.size() <= leaf.size() || path.substr(path.size() - leaf.size()) != leaf) continue;
        const std::string head(path.substr(0, path.size() - leaf.size()));
        if (head.find(".kinetics[") == std::string::npos) return {};
        return head + (leaf.find("slow_") != std::string_view::npos ? "].slow_k" : "].k");
      }
      return {};
    }
  }  // namespace

  Parameter* ParameterRegistry::Find(std::string_view path) {
    Parameter* p = FindExact(path);
    if (p) return p;
    const std::string alias = LumpedAlias(path);
    return alias.empty() ? nullptr : FindExact(alias);
  }

  Parameter* ParameterRegistry::FindExact(std::string_view path) {
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
