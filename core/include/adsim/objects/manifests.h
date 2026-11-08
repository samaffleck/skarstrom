/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Display manifests: the lists the engine defines and front ends show, each
// emitted as JSON so a front end reads them instead of keeping a copy.
//
//   KpiManifest()        the KPIs the engine reports (KPIRecord fields, in the
//                        order of sensitivity::KPINames()): key, label, unit,
//                        description, family, and the reaction KPIs
//   UnitManifest()       the unit categories with their SI unit and suggested
//                        display units, each resolved by the unit parser
//                        (factor, offset), the parser's own symbol overrides
//                        and syntax; CheckUnit() checks any unit text
//   WarningManifest()    the advisory warnings: run-start deck warnings and
//                        end-of-run operating limits, their switch in the
//                        deck's settings, default and message template
//   ExampleDeckManifest  the example decks: file, title, family (display label)
//                        and family_key (AuthorYear), hidden flag, summary (declared in examples/example_decks.json)
//   ComponentManifest()  the components a deck can use: formula, molar mass,
//                        CAS number, property source and whether the gas heat
//                        capacity follows the temperature
//   SeriesManifest()     the bed series variables (Simulation::Series and
//                        Profile, the solver log): key, log name pattern,
//                        unit, bulk or per component, and which are also
//                        logged per adsorbent in a bed with a mixed layer;
//                        the per-reaction series (Series / Profile with a
//                        reaction) in their own reaction_variables list
//
// The results variables (series and profiles) are listed by
// BulkVariableNames() and its siblings in simulation.h (SeriesManifest
// describes the same list); the model kinds by
// the builder (deck_builder.h); the editable values by the parameter registry
// and the deck view (deck_view.h).

#include <string>
#include <vector>

namespace adsim {

  struct ManifestOptions {
    bool pretty = false;  // indented JSON
  };

  // ---- KPIs ------------------------------------------------------------------

  struct KpiDefinition {
    std::string key;          // KPIRecord field name
    std::string label;        // display label, e.g. "Purity"
    std::string short_label;  // table column head, e.g. "purity"
    std::string unit;         // SI unit symbol (parses with the unit parser)
    std::string family;       // purity, recovery, productivity, energy, thermal
    std::string description;
    bool display = true;      // shown in summary tables (false for a duplicate kept for compatibility)
    bool per_step = false;    // one value per step (a list), not one per cycle
    bool diagnostic = false;  // not a reported KPI: a cross-check (quadrature value, mole balance error)
    // A reaction KPI (KPIConfig::reaction): on the records with
    // KPIRecord::reaction_kpis only, absent from the others' as_dict and from
    // a kpis.csv without reactions.
    bool reaction = false;
  };

  // Every KPI, the per-cycle scalars first, in sensitivity::KPINames() order,
  // then the per-step lists, then the reaction KPIs (reaction = true), then
  // the diagnostics (diagnostic = true; the quadrature values and mole
  // balance errors of KPIRecord).
  const std::vector<KpiDefinition>& KpiDefinitions();
  std::string KpiManifest(const ManifestOptions& options = {});

  // ---- units -----------------------------------------------------------------

  struct UnitCheck {
    bool ok = false;
    std::string category;   // UnitCategory as a snake_case name ("pressure", ... "custom")
    std::string dimension;  // DimensionName, e.g. "pressure"
    double factor = 1.0;    // SI = value * factor + offset
    double offset = 0.0;
    std::string message;    // why it failed; empty when ok
  };

  // Parses `text` with the engine's unit parser (ParseUnitExpression). With
  // `expected_category` (a snake_case UnitCategory name such as "pressure"),
  // a unit of another category fails too. "" and "-" are dimensionless.
  UnitCheck CheckUnit(const std::string& text, const std::string& expected_category = "");
  std::string UnitManifest(const ManifestOptions& options = {});

  // ---- warnings --------------------------------------------------------------

  struct WarningDefinition {
    std::string key;       // stable key: the WarningSettings field, or the deck-warning family
    std::string label;
    std::string phase;     // "start" (deck checks before the run) or "end" (operating limits of the last cycle)
    std::string setting;   // registry path of its switch (settings.warn_<key>), "" when always on
    bool default_on = true;
    std::string message;   // message template; {name} placeholders
    // Run-health checks (run_health.h; phase "end", checked after every cycle
    // and reported at the end of the run): the structured warning's code, e.g.
    // "MOLE_FRACTION_SUM". Empty for the other warnings.
    std::string code = {};
  };

  const std::vector<WarningDefinition>& WarningDefinitions();
  std::string WarningManifest(const ManifestOptions& options = {});

  // ---- example decks ---------------------------------------------------------

  // Reads the example-deck manifest (examples/example_decks.json: file,
  // family label, family_key, hidden, summary) and adds each deck's title and description from
  // the deck in `decks_dir`. Throws std::runtime_error when the manifest is
  // malformed, a listed deck is missing, or a deck in `decks_dir` is not
  // listed.
  std::string ExampleDeckManifest(const std::string& manifest_path, const std::string& decks_dir,
                                  const ManifestOptions& options = {});

  // ---- components ------------------------------------------------------------

  struct ComponentInfo {
    std::string name;                  // CoolProp fluid name, as decks write it
    std::string formula;               // e.g. "CO2" ("" when CoolProp has none)
    std::string cas;                   // CAS registry number ("" when none)
    std::vector<std::string> aliases;  // CoolProp's other names for the fluid
    double molar_mass = 0.0;           // [kg/mol]
    // An ideal-gas Cp(T) polynomial exists, so in the ideal-gas mode the gas
    // heat capacity follows the temperature; without one it is held at the
    // initial conditions (Fluid::ComponentsWithHeldHeatCapacity).
    bool ideal_gas_heat_capacity = false;
    // Provenance of that polynomial (NasaPolynomial::species / source); both
    // empty when the heat capacity is held.
    std::string cp_species;
    std::string cp_source;
  };

  // Formula, CAS and molar mass come from CoolProp; empty / 0 for a name
  // CoolProp does not know.
  ComponentInfo DescribeComponent(const std::string& name);
  // Every component a deck can add (Components::componentNames): an object
  // with `cp_reference` (the source table of the polynomials) and `components`.
  std::string ComponentManifest(const ManifestOptions& options = {});

  // ---- series ----------------------------------------------------------------

  // The bed series variables, BulkVariableNames() then ComponentVariableNames()
  // in that order: an object with `variables`, each
  //   key            the name Series / Profile take ("adsorbed_concentration")
  //   label          the log's display name ("Adsorbed Concentration")
  //   name           the log name pattern: "Pressure", "{component} Adsorbed
  //                  Concentration"
  //   unit           SI unit symbol
  //   kind           "bulk" or "component"
  //   porous_only    meaningful in packed cells only
  //   per_adsorbent  also logged per adsorbent (the `adsorbent` argument)
  //   adsorbent_name the per-adsorbent log name pattern, "{component} Adsorbed
  //                  Concentration ({adsorbent})" ({adsorbent} is a label of
  //                  Simulation::SeriesAdsorbents()), "" when not per adsorbent
  // and `per_adsorbent_logged_when`, the rule for when a run logs them;
  // `reaction_variables`, the series logged per reaction and solid (key,
  // label, name pattern "{reaction} Rate ({adsorbent})", unit, kind
  // "reaction", porous_only, description; read with
  // EquationLog::GetReactionRate, not yet with Series / Profile), and
  // `reactions_logged_when`. reaction_heat and reaction_source (in
  // `variables`) and the reaction_variables are logged only for a deck with
  // reactions.
  std::string SeriesManifest(const ManifestOptions& options = {});

}  // namespace adsim
