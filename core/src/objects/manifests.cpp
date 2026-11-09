/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "adsim/objects/manifests.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

#include "CoolProp.h"
#include "adsim/math/unit_parser.h"
#include "adsim/math/user_variables.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/fluid.h"
#include "adsim/objects/model_kinds.h"
#include "adsim/objects/porous_media.h"
#include "adsim/parameters/parameter_registry.h"
#include "adsim/simulation.h"
#include "adsim/solver/log.h"
#include "adsim/solver/operating_limits.h"
#include "adsim/solver/warning_settings.h"
#include "view_json.h"

namespace adsim {

  using namespace view_json;

  // ---- KPIs --------------------------------------------------------------------

  const std::vector<KpiDefinition>& KpiDefinitions() {
    static const std::vector<KpiDefinition> kpis = {
      {"purity", "Purity", "purity", "mol/mol", "purity",
       "Mole fraction of the component in the net product, averaged over the cycle."},
      {"recovery", "Recovery", "recovery", "mol/mol", "recovery",
       "Moles of the component leaving in the net product per mole fed, over the cycle."},
      {"productivity_mass", "Productivity (per adsorbent mass)", "productivity", "kg/(kg s)", "productivity",
       "Mass of the component in the net product per second per kilogram of adsorbent."},
      {"productivity_volume", "Productivity (per bed volume)", "productivity", "kg/(m3 s)", "productivity",
       "Mass of the component in the net product per second per cubic metre of packed bed (no ghost or void end "
       "cells)."},
      {"energy_consumption", "Specific energy", "energy", "J/kg", "energy",
       "Compressor and vacuum-pump work per kilogram of net product (the original name of the electrical energy, "
       "which holds the same value).",
       false},
      {"production_rate", "Production rate", "production", "mol/s", "productivity",
       "Moles of the component in the net product per second of cycle."},
      {"electrical_energy", "Electrical energy", "electrical energy", "J/kg", "energy",
       "Compressor and vacuum-pump work per kilogram of net product."},
      {"thermal_energy", "Thermal energy", "thermal energy", "J/kg", "thermal",
       "Heat supplied to the bed over the cycle (energy sources, heating medium) per kilogram of net product."},
      {"cooling_duty", "Cooling duty", "cooling duty", "J/kg", "thermal",
       "Heat removed from the bed over the cycle per kilogram of net product, as a positive magnitude; for "
       "information, not part of the thermal energy."},
      {"step_thermal_energy", "Thermal energy by step", "step thermal energy", "J/kg", "thermal",
       "Each step's share of the thermal energy, in the deck's step order.", true, true},
      {"step_cooling_duty", "Cooling duty by step", "step cooling duty", "J/kg", "thermal",
       "Each step's share of the cooling duty, in the deck's step order.", true, true},
      // Reaction KPIs (KPIConfig::reaction): only on the key reactant's and the products' records.
      {"conversion", "Conversion", "conversion", "mol/mol", "reaction",
       "Moles of the key reactant the reactions consumed per mole of it fed, over the cycle (key reactant's "
       "record).", true, false, false, true},
      {"selectivity", "Selectivity", "selectivity", "mol/mol", "reaction",
       "Moles of the product the reactions made per mole of key reactant they consumed, over the cycle; no "
       "stoichiometric normalisation.", true, false, false, true},
      {"reaction_yield", "Yield", "yield", "mol/mol", "reaction",
       "Moles of the product the reactions made per mole of key reactant fed, over the cycle; no stoichiometric "
       "normalisation.", true, false, false, true},
      {"space_time_yield", "Space-time yield", "space-time yield", "kg/(kg s)", "reaction",
       "Mass of the product the reactions made per second of cycle per kilogram of catalyst (the solids that "
       "carry the reactions), or per cubic metre of bed (kg/(m3 s), space_time_yield_unit) with the bed-volume "
       "basis.", true, false, false, true},
      {"peak_temperature", "Peak temperature", "peak T", "K", "reaction",
       "Highest bed temperature over the cycle (the solid temperature in an LTNE run).", true, false, false, true},
      {"peak_temperature_position", "Peak temperature position", "peak T position", "m", "reaction",
       "Centre of the cell of the peak temperature, from the bottom of the column.", true, false, false, true},
      {"peak_temperature_time", "Peak temperature time", "peak T time", "s", "reaction",
       "Time of the peak temperature from the start of the cycle.", true, false, false, true},
      {"reaction_generation", "Reaction generation", "generated", "mol", "reaction",
       "Moles of the component the reactions made in the bed over the cycle (negative when consumed).", true,
       false, false, true},
      // Diagnostics: the same KPIs from the trapezoid over the logged samples, and the cycle mole balance.
      {"trapezoid_purity", "Purity (trapezoid)", "purity (trap.)", "mol/mol", "purity",
       "Purity from the trapezoid over the logged endpoint samples instead of the KPI quadratures (the reported "
       "basis); equal to purity when the quadratures were off.", false, false, true},
      {"trapezoid_recovery", "Recovery (trapezoid)", "recovery (trap.)", "mol/mol", "recovery",
       "Recovery from the trapezoid over the logged endpoint samples.", false, false, true},
      {"trapezoid_productivity_mass", "Productivity per adsorbent mass (trapezoid)", "productivity mass (trap.)",
       "kg/(kg s)", "productivity", "Productivity per kilogram of adsorbent from the trapezoid over the logged "
       "samples.", false, false, true},
      {"trapezoid_productivity_volume", "Productivity per bed volume (trapezoid)", "productivity volume (trap.)",
       "kg/(m3 s)", "productivity", "Productivity per cubic metre of packed bed from the trapezoid over the logged "
       "samples.", false, false, true},
      {"trapezoid_energy_consumption", "Specific energy (trapezoid)", "energy (trap.)", "J/kg", "energy",
       "Compressor and vacuum-pump work per kilogram of net product from the trapezoid over the logged samples.",
       false, false, true},
      {"trapezoid_production_rate", "Production rate (trapezoid)", "production rate (trap.)", "mol/s",
       "productivity", "Production rate from the trapezoid over the logged samples.", false, false, true},
      {"mass_balance_error", "Mole balance error", "mass balance error", "mol/mol", "balance",
       "Cycle mole balance of the component as a fraction of the moles that entered (zero for exact conservation), "
       "with the endpoint flows from the KPI quadratures; NaN without the quadratures.", false, false, true},
      {"trapezoid_mass_balance_error", "Mole balance error (trapezoid)", "mass balance error (trap.)", "mol/mol",
       "balance", "Cycle mole balance error with the endpoint flows from the trapezoid over the logged samples.",
       false, false, true},
      {"uncounted_inflow", "Uncounted inflow", "uncounted inflow", "mol", "balance",
       "Moles of the component that entered the bed in the cycle through unconnected endpoints that are neither "
       "feed nor recycle endpoints of it (gross inflow, e.g. a vent that backflows).", false, false, true},
      {"uncounted_inflow_fraction", "Uncounted inflow fraction", "uncounted inflow fraction", "mol/mol", "balance",
       "Uncounted inflow over the counted feed of the component.", false, false, true},
    };
    return kpis;
  }

  std::string KpiManifest(const ManifestOptions& options) {
    Document doc;
    auto& a = doc.GetAllocator();
    Value list(rapidjson::kArrayType);
    for (const auto& k : KpiDefinitions()) {
      Value v(rapidjson::kObjectType);
      Set(v, "key", k.key, a);
      Set(v, "label", k.label, a);
      Set(v, "short_label", k.short_label, a);
      Set(v, "unit", k.unit, a);
      Set(v, "family", k.family, a);
      Set(v, "description", k.description, a);
      Set(v, "display", k.display, a);
      Set(v, "per_step", k.per_step, a);
      Set(v, "diagnostic", k.diagnostic, a);
      Set(v, "reaction", k.reaction, a);
      list.PushBack(v, a);
    }
    return Write(list, options.pretty);
  }

  // ---- units -------------------------------------------------------------------

  namespace {

    // UnitCategory values by their snake_case name, in the builder's order.
    const std::vector<std::pair<std::string, UnitCategory>>& Categories() {
      static const std::vector<std::pair<std::string, UnitCategory>> out = [] {
        std::map<std::string, UnitCategory> by_name;
        for (int i = 0; i <= static_cast<int>(UnitCategory::SpecificEnergy); ++i) {
          const auto c = static_cast<UnitCategory>(i);
          by_name.emplace(build::VariableCategoryKind(c), c);
        }
        std::vector<std::pair<std::string, UnitCategory>> v;
        for (const auto& name : build::VariableCategoryKinds())
          if (auto it = by_name.find(name); it != by_name.end()) v.emplace_back(name, it->second);
        return v;
      }();
      return out;
    }

    // Display units a front end offers for each category, most common first.
    // Only the symbols are listed: factor, offset and category come from the
    // unit parser (a test checks each parses into its category). Any other
    // unit the parser accepts can be typed; check it with CheckUnit.
    const std::map<std::string, std::vector<std::string>>& SuggestedUnits() {
      static const std::map<std::string, std::vector<std::string>> units = {
        {"dimensionless", {"-"}},
        {"time", {"s", "ms", "min", "h", "day"}},
        {"length", {"m", "cm", "mm", "um", "nm", "in", "ft"}},
        {"temperature", {"K", "degC", "degF"}},
        {"pressure", {"Pa", "kPa", "MPa", "bar", "mbar", "atm", "psi", "Torr", "mmHg"}},
        {"energy", {"J", "kJ", "MJ", "kWh"}},
        {"specific_molar_energy", {"J/mol", "kJ/mol", "kcal/mol"}},
        {"inverse_time", {"1/s", "1/min", "1/h"}},
        {"velocity", {"m/s", "cm/s", "m/min", "m/h"}},
        {"volumetric_flow_rate", {"m3/s", "m3/h", "L/s", "L/min", "cm3/s", "ccm", "mL/min"}},
        {"molar_flow_rate", {"mol/s", "mmol/s", "mol/min", "kmol/h", "SCCM", "SCFM", "Nm3/h"}},
        {"mass_flow_rate", {"kg/s", "kg/h", "g/s", "g/min"}},
        {"mole_fraction", {"mol/mol", "%", "ppm"}},
        {"viscosity", {"Pa.s", "mPa.s", "cP"}},
        {"density", {"kg/m3", "g/cm3", "g/L"}},
        {"molar_concentration", {"mol/m3", "mol/L", "mmol/L"}},
        {"heat_capacity", {"J/(kg-K)", "kJ/(kg-K)"}},
        {"molar_mass", {"kg/mol", "g/mol"}},
        {"thermal_conductivity", {"W/(m-K)"}},
        {"heat_transfer_coefficient", {"W/(m2-K)", "kW/(m2-K)"}},
        {"pressure_gradient", {"Pa/m", "kPa/m", "bar/m"}},
        {"adsorbed_concentration", {"mol/kg", "mmol/g"}},
        {"mass", {"kg", "g", "mg"}},
        {"volume", {"m3", "L", "cm3", "mL"}},
        {"molar_amount", {"mol", "mmol", "kmol"}},
        {"area", {"m2", "cm2", "mm2"}},
        {"diffusivity", {"m2/s", "cm2/s"}},
        {"specific_energy", {"J/kg", "kJ/kg", "MJ/kg", "kWh/kg"}},
        {"custom", {}},
      };
      return units;
    }

    bool Interchangeable(const std::string& a, const std::string& b) {
      auto dimensionless = [](const std::string& c) { return c == "dimensionless" || c == "mole_fraction"; };
      return a == b || (dimensionless(a) && dimensionless(b));
    }

  }  // namespace

  UnitCheck CheckUnit(const std::string& text, const std::string& expected_category) {
    UnitCheck out;
    UnitInfo info;
    std::string why;
    if (!ParseUnitExpression(text, info, &why)) {
      out.message = "unknown unit '" + text + "': " + why;
      return out;
    }
    out.factor = info.factor;
    out.offset = info.offset;
    out.dimension = DimensionName(info.dimension);
    out.category = build::VariableCategoryKind(UnitCategoryOf(text, info.dimension));
    if (!expected_category.empty() && expected_category != "custom" &&
        !Interchangeable(out.category, expected_category)) {
      out.message = "'" + text + "' is a unit of " + out.dimension + " (" + out.category + "), not " + expected_category;
      return out;
    }
    out.ok = true;
    return out;
  }

  std::string UnitManifest(const ManifestOptions& options) {
    Document doc;
    auto& a = doc.GetAllocator();
    Value root(rapidjson::kObjectType);

    auto unit_entry = [&](const std::string& symbol) {
      const UnitCheck c = CheckUnit(symbol);
      Value u(rapidjson::kObjectType);
      Set(u, "symbol", symbol, a);
      Set(u, "category", c.category, a);
      Set(u, "dimension", c.dimension, a);
      Set(u, "factor", Num(c.factor), a);
      Set(u, "offset", Num(c.offset), a);
      return u;
    };

    Value categories(rapidjson::kArrayType);
    for (const auto& [name, category] : Categories()) {
      Value c(rapidjson::kObjectType);
      Set(c, "key", name, a);
      Set(c, "name", UnitCategoryName(category), a);
      Set(c, "si_unit", category == UnitCategory::Custom ? "" : UnitCategorySymbol(category), a);
      Value units(rapidjson::kArrayType);
      if (auto it = SuggestedUnits().find(name); it != SuggestedUnits().end())
        for (const auto& symbol : it->second) units.PushBack(unit_entry(symbol), a);
      Set(c, "units", std::move(units), a);
      categories.PushBack(c, a);
    }
    Set(root, "categories", std::move(categories), a);

    Value overrides(rapidjson::kArrayType);
    for (const auto& symbol : UnitOverrideSymbols()) overrides.PushBack(unit_entry(symbol), a);
    Set(root, "overrides", std::move(overrides), a);

    Set(root, "syntax",
        StrList({
            "Any unit symbol of the LLNL units library with an SI prefix (kPa, mbar, um, min, h, atm, psi, degC, "
            "kWh), read strictly: a misspelt symbol is an error.",
            "Only units of length, mass, time, temperature and amount of substance (and their products).",
            "C, F and A are degC, degF and the angstrom; ccm, SCCM, SCFM and Nm3 are gas-flow units at 0 degC and "
            "1 atm.",
            "Exponents: m2, m^2, m**2, m\xC2\xB2, s^-1.",
            "Multiplication: '*', '.', '\xC2\xB7', '-' or a space. Division: everything after a '/' up to the next "
            "'/' is in the denominator (W/m2-K = W/(m2 K)); parentheses group.",
            "'-', '1' and '' are dimensionless; '%' is 0.01 and 'ppm' 1e-6.",
            "A temperature unit with an offset (degC, degF) is absolute only on its own; inside a compound unit it "
            "is a temperature difference.",
            "In a value, a unit follows the number ('2 min') or ends an expression in brackets ('2*t_feed [min]').",
        }, a), a);
    return Write(root, options.pretty);
  }

  // ---- warnings ------------------------------------------------------------------

  const std::vector<WarningDefinition>& WarningDefinitions() {
    static const std::vector<WarningDefinition> warnings = [] {
      const WarningSettings defaults{};
      std::vector<WarningDefinition> w = {
        // Run start (deck_warnings.h), switchable.
        {"pressure_step_change", "Pressure step change at a step start", "start",
         "settings.warn_pressure_step_change", defaults.pressure_step_change,
         "step \"{step}\" ({end}) starts with a step change in boundary pressure: its constant {pressure} against "
         "{expected}. The solver has to resolve an instantaneous pressure change there; a linear or exponential ramp "
         "to {pressure} avoids it."},
        {"wall_thickness", "Wall thinner than its pressure design", "start", "settings.warn_wall_thickness",
         defaults.wall_thickness,
         "the wall is {thickness} mm thick, thinner than the {required} mm {code} requires for its design. Use "
         "Calculate in the wall design, or a thicker wall."},
        // Run start, always reported.
        {"mixing_rule", "IAST with a competitive isotherm", "start", "", true,
         "adsorbent '{adsorbent}' uses the IAST mixing rule with the competitive isotherm '{isotherm}' for "
         "{components}. IAST takes every isotherm in its pure-component limit, so that kind's own competition term "
         "is not applied."},
        {"gas_heat_capacity", "Gas heat capacity held at the initial conditions", "start", "", true,
         "{component} has no ideal-gas heat capacity polynomial, so its gas heat capacity is CoolProp's value at the "
         "initial pressure and temperature, held for the whole run."},
        {"layer_composition", "Solid at mass fraction 0", "start", "", true,
         "adsorbent {adsorbent} has mass fraction 0 in layer {layer}, so it takes no part in the run."},
        {"film_correlation", "Single-sphere film correlation in a bed", "start", "", true,
         "layer '{layer}' uses the Ranz-Marshall film correlation, which is for a single sphere, not a packed bed: "
         "in a bed it under-predicts the particle heat and mass transfer coefficients at moderate Reynolds number."},
        {"structured_layer", "Structured layer heat transfer", "start", "", true,
         "the bed has a structured layer and uses local thermal equilibrium (LTE); LTNE is recommended. Or: layer "
         "'{layer}': the bed-to-wall heat transfer coefficient is 0, so the column wall exchanges no heat with it."},
        {"reaction", "Reaction not applied as written", "start", "", true,
         "solid '{solid}' has reactions but no layer uses it; or reaction '{reaction}' has a nonzero heat of "
         "reaction in an isothermal run; or its element balance cannot be checked."},
        // End of a run, over the last logged cycle (operating_limits.h), switchable.
        {"design_pressure", "Pressure above the design pressure", "end", "settings.warn_design_pressure",
         defaults.design_pressure,
         "the pressure reaches {value} in step {step}, above the wall's design pressure of {limit}."},
        {"vacuum", "Vacuum without vacuum service", "end", "settings.warn_vacuum", defaults.vacuum,
         "the pressure falls to {value} in step {step}, below atmospheric, but the wall is not designed for vacuum "
         "service (external pressure)."},
        {"design_temperature", "Wall above its design temperature", "end", "settings.warn_design_temperature",
         defaults.design_temperature,
         "the wall temperature reaches {value} in step {step}, above the wall's design temperature of {limit}."},
        {"fluidisation", "Approach to fluidisation above 1", "end", "settings.warn_fluidisation",
         defaults.fluidisation,
         "the approach to fluidisation reaches {value} in step {step}: the upward pressure drop exceeds the weight "
         "of the bed, which would lift unless it is held down."},
        {"channel_reynolds", "Channel flow no longer laminar", "end", "settings.warn_channel_reynolds",
         defaults.channel_reynolds,
         "the channel Reynolds number of a structured layer reaches {value} in step {step}, above {limit}: the "
         "flow is no longer laminar."},
        {"outlet_backflow", "Backflow into the bed through an outlet", "end", "settings.warn_outlet_backflow",
         defaults.outlet_backflow,
         "gas flows back into the bed through the outlet {endpoint} in step {step} of cycle {cycle}: {moles} mol "
         "in, {value} of the cycle feed (more than {limit}). It carries the bed-end gas, not the stream's, and is "
         "not counted as feed. Turn on the valve's check valve (endpoint[ID].valve.check_valve) to stop it."},
        // Run health, checked after every cycle, reported at the end (run_health.h), switchable.
        {"mole_fraction_sum", "Mole fractions do not sum to 1", "end", "settings.warn_mole_fraction_sum",
         defaults.mole_fraction_sum,
         "the gas mole fractions sum to 1 +/- {value} in step {step} of cycle {cycle} (more than {limit}): the "
         "composition is not conserved.",
         "MOLE_FRACTION_SUM"},
        {"mass_balance_stuck", "Cycle mole balance stays open", "end", "settings.warn_mass_balance_stuck",
         defaults.mass_balance_stuck,
         "the cycle mole balance of {component} stays at {value} in cycle {cycle} (above {limit}, no longer "
         "falling) while the end-of-cycle state has settled: moles are gained or lost.",
         "MASS_BALANCE_STUCK"},
        {"step_retries", "Many rejected integrator steps", "end", "settings.warn_step_retries",
         defaults.step_retries,
         "the integrator rejected {value} attempts per accepted step in cycle {cycle} (more than {limit}).",
         "STEP_RETRIES"},
        {"history_truncated", "Log shorter than the diagnostics need", "end", "settings.warn_history_truncated",
         defaults.history_truncated,
         "the log kept {value} cycles, fewer than {limit}: series and profiles cover only the logged cycles "
         "(settings.max_cycles_to_keep).",
         "HISTORY_TRUNCATED"},
        {"reaction_stiff", "Step collapse at a reacting cell", "end", "settings.warn_reaction_stiff",
         defaults.reaction_stiff,
         "the integrator's step collapsed at a reacting cell in step {step} of cycle {cycle}: {value} steps shorter "
         "than 1e-4 s (at least {limit}), with near-zero {component}. Start the bed nearer the reaction's state "
         "(docs/reactions.md, \"Start-up of fast reactions\").",
         "REACTION_STIFF"},
      };
      return w;
    }();
    return warnings;
  }

  std::string ModelKindManifest(const ManifestOptions& options) {
    Document doc;
    auto& a = doc.GetAllocator();
    Value out(rapidjson::kObjectType);
    auto fragments = [&](const ResistanceDescriptor::Latex& l) {
      Value f(rapidjson::kObjectType);
      Set(f, "term", l.term, a);
      Set(f, "where", l.where, a);
      Set(f, "where_isothermal", l.where_isothermal, a);
      return f;
    };
    for (ModelFamily family : AllModelFamilies()) {
      Value kinds(rapidjson::kArrayType);
      for (const KindDescriptor& k : KindDescriptors(family)) {
        Value v(rapidjson::kObjectType);
        Set(v, "tag", k.tag, a);
        Set(v, "label", k.label, a);
        Set(v, "legacy", k.legacy, a);
        Value ps(rapidjson::kArrayType);
        for (const KindParameter& p : k.parameters) {
          Value pv(rapidjson::kObjectType);
          Set(pv, "name", p.name, a);
          Set(pv, "unit", p.unit, a);
          Set(pv, "default", Num(p.default_value), a);
          Set(pv, "lower", Num(p.lower), a);
          Set(pv, "upper", Num(p.upper), a);
          Set(pv, "tunable", p.tunable, a);
          Set(pv, "label", p.Label(), a);
          ps.PushBack(pv, a);
        }
        Set(v, "parameters", std::move(ps), a);
        if (!k.resistances.empty()) {
          Value rs(rapidjson::kArrayType);
          for (const ResistanceDescriptor& r : k.resistances) {
            Value rv(rapidjson::kObjectType);
            Set(rv, "tag", r.tag, a);
            Set(rv, "label", r.label, a);
            if (r.parent.empty()) SetNull(rv, "parent", a);
            else Set(rv, "parent", r.parent, a);
            Set(rv, "scaled_by_capacity", r.scaled_by_capacity, a);
            Set(rv, "parameters", StrList(r.parameters, a), a);
            Value fr(rapidjson::kObjectType);
            Set(fr, "sphere", fragments(r.sphere), a);
            Set(fr, "slab", fragments(r.slab), a);
            Set(fr, "cylinder", fragments(r.cylinder), a);
            Set(rv, "fragments", std::move(fr), a);
            Set(rv, "source", r.source, a);
            rs.PushBack(rv, a);
          }
          Set(v, "resistances", std::move(rs), a);
          Value ms(rapidjson::kArrayType);
          for (const CapacityModeDescriptor& m : k.capacity_modes) {
            Value mv(rapidjson::kObjectType);
            Set(mv, "tag", m.tag, a);
            Set(mv, "label", m.label, a);
            Set(mv, "latex", m.latex, a);
            ms.PushBack(mv, a);
          }
          Set(v, "capacity_modes", std::move(ms), a);
          Value pr(rapidjson::kArrayType);
          for (const PresetDescriptor& p : k.presets) {
            Value pv(rapidjson::kObjectType);
            Set(pv, "tag", p.tag, a);
            Set(pv, "on", StrList(p.on, a), a);
            const auto m = static_cast<std::size_t>(p.capacity_mode);
            Set(pv, "capacity_mode", m < k.capacity_modes.size() ? k.capacity_modes[m].tag : std::string("none"), a);
            Value vals(rapidjson::kObjectType);
            for (const auto& [name, x] : p.values) Set(vals, name, Num(x), a);
            Set(pv, "values", std::move(vals), a);
            pr.PushBack(pv, a);
          }
          Set(v, "presets", std::move(pr), a);
        }
        kinds.PushBack(v, a);
      }
      Set(out, ModelFamilyName(family), std::move(kinds), a);
    }
    // The packed particle shapes (ParticleShape), in the families' form.
    Value shapes(rapidjson::kArrayType);
    for (std::size_t i = 0; i < ParticleShapeKinds().size(); ++i) {
      Value v(rapidjson::kObjectType);
      Set(v, "tag", ParticleShapeKinds()[i], a);
      Set(v, "label", std::string(ParticleShapeList[i]), a);
      Set(v, "legacy", false, a);
      Set(v, "parameters", Value(rapidjson::kArrayType), a);
      shapes.PushBack(v, a);
    }
    Set(out, "particle_shape", std::move(shapes), a);
    return Write(out, options.pretty);
  }

  std::string WarningManifest(const ManifestOptions& options) {
    Document doc;
    auto& a = doc.GetAllocator();
    Value list(rapidjson::kArrayType);
    for (const auto& w : WarningDefinitions()) {
      Value v(rapidjson::kObjectType);
      Set(v, "key", w.key, a);
      Set(v, "label", w.label, a);
      Set(v, "phase", w.phase, a);
      if (w.setting.empty()) SetNull(v, "setting", a);
      else Set(v, "setting", w.setting, a);
      Set(v, "default_on", w.default_on, a);
      Set(v, "message", w.message, a);
      if (w.code.empty()) SetNull(v, "code", a);
      else Set(v, "code", w.code, a);
      list.PushBack(v, a);
    }
    return Write(list, options.pretty);
  }

  // ---- example decks -------------------------------------------------------------

  namespace {

    std::string ReadFile(const std::filesystem::path& path) {
      std::ifstream in(path, std::ios::binary);
      if (!in) throw std::runtime_error("cannot read " + path.string());
      std::ostringstream ss;
      ss << in.rdbuf();
      return ss.str();
    }

    std::string MemberString(const Value& o, const char* key) {
      if (!o.IsObject()) return {};
      auto it = o.FindMember(key);
      return it != o.MemberEnd() && it->value.IsString() ? std::string(it->value.GetString()) : std::string{};
    }

  }  // namespace

  std::string ExampleDeckManifest(const std::string& manifest_path, const std::string& decks_dir,
                                  const ManifestOptions& options) {
    namespace fs = std::filesystem;
    Document manifest;
    const std::string text = ReadFile(manifest_path);
    manifest.Parse(text.c_str());
    if (manifest.HasParseError() || !manifest.IsObject() || !manifest.HasMember("decks") ||
        !manifest["decks"].IsArray())
      throw std::runtime_error("example deck manifest " + manifest_path + ": expected an object with a \"decks\" array");

    std::set<std::string> on_disk;
    for (const auto& entry : fs::directory_iterator(decks_dir))
      if (entry.is_regular_file() && entry.path().extension() == ".json") on_disk.insert(entry.path().filename().string());

    Document out;
    auto& a = out.GetAllocator();
    Value list(rapidjson::kArrayType);
    std::set<std::string> listed;
    for (const auto& d : manifest["decks"].GetArray()) {
      const std::string file = MemberString(d, "file");
      if (file.empty()) throw std::runtime_error("example deck manifest: an entry has no \"file\"");
      if (!listed.insert(file).second) throw std::runtime_error("example deck manifest: " + file + " is listed twice");
      if (!on_disk.count(file)) throw std::runtime_error("example deck manifest: " + file + " is not in " + decks_dir);
      const std::string family = MemberString(d, "family");
      if (family.empty()) throw std::runtime_error("example deck manifest: " + file + " has no \"family\"");

      // Title and description from the deck itself.
      Document deck;
      const std::string deck_text = ReadFile(fs::path(decks_dir) / file);
      deck.Parse(deck_text.c_str());
      std::string title, description;
      if (!deck.HasParseError() && deck.IsObject()) {
        for (const char* root : {"value0", "simulation"}) {
          auto r = deck.FindMember(root);
          if (r == deck.MemberEnd() || !r->value.IsObject()) continue;
          auto info = r->value.FindMember("System_Information");
          if (info == r->value.MemberEnd()) continue;
          title = MemberString(info->value, "Title");
          description = MemberString(info->value, "Description");
        }
      }

      Value v(rapidjson::kObjectType);
      Set(v, "file", file, a);
      Set(v, "name", fs::path(file).stem().string(), a);
      Set(v, "title", title.empty() ? fs::path(file).stem().string() : title, a);
      Set(v, "family", family, a);
      const std::string family_key = MemberString(d, "family_key");
      Set(v, "family_key", family_key.empty() ? family : family_key, a);
      const auto hidden = d.FindMember("hidden");
      Set(v, "hidden", hidden != d.MemberEnd() && hidden->value.IsBool() && hidden->value.GetBool(), a);
      Set(v, "summary", MemberString(d, "summary"), a);
      Set(v, "description", description, a);
      list.PushBack(v, a);
    }
    for (const auto& file : on_disk)
      if (!listed.count(file))
        throw std::runtime_error("example deck manifest: " + file + " in " + decks_dir + " is not listed in " +
                                 manifest_path);
    return Write(list, options.pretty);
  }

  // ---- components ----------------------------------------------------------------

  namespace {

    std::string CoolPropString(const std::string& fluid, const std::string& param) {
      try {
        return CoolProp::get_fluid_param_string(fluid, param);
      } catch (...) {
        return {};
      }
    }

    // CoolProp writes formulas with subscripts in braces, "C_{1}O_{2}": drop
    // the ones and the markup ("CO2").
    std::string PlainFormula(std::string f) {
      std::string out;
      for (size_t i = 0; i < f.size(); ++i) {
        if (f[i] == '_' && i + 1 < f.size() && f[i + 1] == '{') {
          const size_t close = f.find('}', i);
          if (close == std::string::npos) break;
          const std::string n = f.substr(i + 2, close - i - 2);
          if (n != "1") out += n;
          i = close;
        } else if (f[i] != '{' && f[i] != '}') {
          out += f[i];
        }
      }
      return out;
    }

  }  // namespace

  ComponentInfo DescribeComponent(const std::string& name) {
    ComponentInfo info;
    info.name = name;
    info.formula = PlainFormula(CoolPropString(name, "formula"));
    info.cas = CoolPropString(name, "CAS");
    const std::string aliases = CoolPropString(name, "aliases");
    std::stringstream ss(aliases);
    for (std::string item; std::getline(ss, item, ',');) {
      item.erase(0, item.find_first_not_of(' '));
      item.erase(item.find_last_not_of(' ') + 1);
      if (!item.empty() && item != name) info.aliases.push_back(item);
    }
    try {
      info.molar_mass = CoolProp::Props1SI(name, "molar_mass");
    } catch (...) {
      info.molar_mass = 0.0;
    }
    if (!std::isfinite(info.molar_mass) || info.molar_mass <= 0.0 || info.molar_mass > 1e3) info.molar_mass = 0.0;
    if (const NasaPolynomial* p = Fluid::IdealGasCoeffs(name)) {
      info.ideal_gas_heat_capacity = true;
      info.cp_species = std::string(p->species);
      info.cp_source = std::string(p->source);
    }
    return info;
  }

  std::string ComponentManifest(const ManifestOptions& options) {
    Document doc;
    auto& a = doc.GetAllocator();
    Value list(rapidjson::kArrayType);
    for (const auto& name : Components::componentNames) {
      const ComponentInfo c = DescribeComponent(name);
      Value v(rapidjson::kObjectType);
      Set(v, "name", c.name, a);
      Set(v, "formula", c.formula, a);
      Set(v, "molar_mass", Num(c.molar_mass > 0 ? c.molar_mass : std::numeric_limits<double>::quiet_NaN()), a);
      Set(v, "cas", c.cas, a);
      Set(v, "properties_source", "CoolProp", a);
      Set(v, "coolprop_name", c.name, a);
      Set(v, "aliases", StrList(c.aliases, a), a);
      Set(v, "ideal_gas_heat_capacity", c.ideal_gas_heat_capacity, a);
      Set(v, "heat_capacity_held", !c.ideal_gas_heat_capacity, a);
      Set(v, "cp_species", c.cp_species, a);
      Set(v, "cp_source", c.cp_source, a);
      list.PushBack(v, a);
    }
    Value out(rapidjson::kObjectType);
    Set(out, "cp_reference", std::string(kIdealGasHeatCapacityReference), a);
    out.AddMember("components", list, a);
    return Write(out, options.pretty);
  }

  std::string SeriesManifest(const ManifestOptions& options) {
    Document doc;
    auto& a = doc.GetAllocator();
    Value list(rapidjson::kArrayType);
    // BulkVariableNames() / ComponentVariableNames() follow the enum order of
    // VariableTypeList / ComponentVariableTypeList (static_asserts in
    // simulation.cpp), so entry i names enum value i.
    const auto bulk = BulkVariableNames();
    for (size_t i = 0; i < bulk.size(); ++i) {
      const auto t = static_cast<VariableType>(i);
      const std::string label = VariableTypeList[i];
      Value v(rapidjson::kObjectType);
      Set(v, "key", bulk[i], a);
      Set(v, "label", label, a);
      Set(v, "name", label, a);
      std::string unit = GetVariableUnits(t);
      if (unit.size() >= 2 && unit.front() == '[' && unit.back() == ']') unit = unit.substr(1, unit.size() - 2);
      Set(v, "unit", unit, a);
      Set(v, "kind", "bulk", a);
      Set(v, "porous_only", IsPorousOnlyVariable(t), a);
      Set(v, "per_adsorbent", false, a);
      Set(v, "adsorbent_name", "", a);
      list.PushBack(v, a);
    }
    const auto comp = ComponentVariableNames();
    for (size_t i = 0; i < comp.size(); ++i) {
      const auto t = static_cast<ComponentVariableType>(i);
      const std::string label = ComponentVariableTypeList[i];
      Value v(rapidjson::kObjectType);
      Set(v, "key", comp[i], a);
      Set(v, "label", label, a);
      Set(v, "name", "{component} " + label, a);
      std::string unit = GetComponentVariableUnits(t);
      if (unit.size() >= 2 && unit.front() == '[' && unit.back() == ']') unit = unit.substr(1, unit.size() - 2);
      Set(v, "unit", unit, a);
      Set(v, "kind", "component", a);
      Set(v, "porous_only", IsPorousOnlyComponentVariable(t), a);
      Set(v, "per_adsorbent", IsPerAdsorbentVariable(t), a);
      Set(v, "adsorbent_name",
          IsPerAdsorbentVariable(t) ? PerAdsorbentVariableName("{component}", t, "{adsorbent}") : std::string(), a);
      list.PushBack(v, a);
    }
    // The rate of each reaction on each solid (logged only for a deck with
    // reactions; EquationLog::GetReactionRate).
    Value reactions(rapidjson::kArrayType);
    {
      Value v(rapidjson::kObjectType);
      Set(v, "key", "reaction_rate", a);
      Set(v, "label", "Rate", a);
      Set(v, "name", ReactionRateVariableName("{reaction}", "{adsorbent}"), a);
      std::string unit = kReactionRateUnits;
      unit = unit.substr(1, unit.size() - 2);
      Set(v, "unit", unit, a);
      Set(v, "kind", "reaction", a);
      Set(v, "argument", "reaction", a);  // Series / Profile take reaction= (and adsorbent= to pick the solid)
      Set(v, "porous_only", true, a);
      Set(v, "description",
          "effectiveness factor times the intrinsic rate of one reaction, per kg of the solid that carries it, "
          "in moles of reaction extent; 0 in the cells without that solid",
          a);
      reactions.PushBack(v, a);
    }
    Value out(rapidjson::kObjectType);
    out.AddMember("variables", list, a);
    out.AddMember("reaction_variables", reactions, a);
    Set(out, "reactions_logged_when", "some solid of some layer has a reaction (reaction_heat, reaction_source and "
        "the reaction_variables)", a);
    Set(out, "per_adsorbent_logged_when",
        "some layer of the bed holds more than one solid; otherwise each cell's series is that of its one "
        "solid and the per-adsorbent series are read off it",
        a);
    return Write(out, options.pretty);
  }

}  // namespace adsim
