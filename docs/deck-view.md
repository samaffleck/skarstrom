# Deck View and Manifests

A front end (a web workspace, a notebook, a script that renders a form) needs two things from the
engine: the deck's structure with its editable values, and the fixed lists the engine defines
(KPIs, units, warnings, example decks, components). The engine emits both as JSON, so a front end
does not need its own copy of either.

## The deck view

```python
import skarstrom as sk

sim = sk.examples.Perez2019_A.load()
view = sim.deck_view()          # dict
schema = sk.deck_view_schema()  # JSON Schema (draft 4) of the view
```

In C++ the same calls are `adsim::DeckView(UserData&)` and `adsim::DeckViewSchema()`
(`adsim/objects/deck_view.h`). Both return JSON text.

The view has the deck's structure: components (name, formula, molar mass), adsorbents (each
component's isotherm kind and rate-model kinds, reactions), layers (kind, transport models,
composition), the bed and its wall, steps (time mode, the bottom and top ends with their streams,
profiles, pumps and buffer tanks, energy sources, the heating medium), connections, the KPI
configuration (each component's entry also lists its `kpi[C].tracked` value), user variables, the
`cost` (cost estimate) and `sensitivity` (the stored study's settings) sections and the `info`,
`settings`, `css` and `energy_balance` sections.
Every node has its stable id (step, endpoint, layer, adsorbent and reaction id). Model kinds are
given by the builder's kind names (`skarstrom.build_kinds()`).

Each node's `values` list holds the editable values that belong to it. Every value is an entry of
the [parameter registry](parameterisation.md), with its fields copied from the registry:

| Field | Meaning |
|---|---|
| `path` | Stable, id-based registry path (`endpoint[5].temperature`) |
| `alias` | The same address with display names (`step[Feed].bottom[Feed inlet].temperature`) |
| `key` | The path relative to its node (`temperature`, `pump.efficiency`) |
| `label`, `category`, `subcategory` | Display names |
| `kind` | `real`, `expression`, `integer`, `boolean` or `string` |
| `value` | Number in SI units (boolean or string for those kinds; `null` when not finite) |
| `text` | The text form: the expression for an expression entry |
| `formula` | True when an expression entry is a formula of user variables |
| `unit`, `unit_checked` | SI unit symbol (`-` dimensionless); whether values given with a unit are converted |
| `tunable`, `lower`, `upper` | For optimisers and studies |
| `choices`, `choice` | For an enumerated value: every `{name, index}` and the current name |

Each node claims the registry entries under its path prefix (`endpoint[5].pump.` for a pump,
`endpoint[5].` for the rest of the stream). So every registry entry appears exactly once. An entry
no node claims is listed under `other`, so a newer registry family is never hidden from the view.

The view is read-only. Change a value by its path (`sim.set(path, value)`, `Parameter::SetValue` or
`SetText` in C++), change structure with the [builder](python.md), then take a new view. The test
suite checks, for every example deck, that the view validates against the schema and that setting
every value back by path leaves the deck unchanged.

`view_version` changes when the view's shape changes incompatibly.

## Manifests

| Python | C++ (`adsim/objects/manifests.h`) | Contents |
|---|---|---|
| `sk.kpi_manifest()` | `KpiManifest()`, `KpiDefinitions()` | Each KPI the engine reports, in report order: `key` (the `KPIRecord` field), `label`, `short_label`, SI `unit`, `family` (purity, recovery, productivity, energy, thermal), `description`, `display`, `per_step`, `diagnostic`. The reported KPIs come first, then the diagnostics (`diagnostic: true`: the `trapezoid_*` values and the `mass_balance_error` fractions, [mol/mol]), which a front end can show in a separate section |
| `sk.unit_manifest()` | `UnitManifest()` | Unit categories (`key`, `name`, `si_unit`) with suggested display units, the parser's symbol overrides and syntax notes |
| `sk.check_unit(text, category)` | `CheckUnit(text, category)` | `{ok, category, dimension, factor, offset, message}` from the engine's unit parser |
| `sk.warning_manifest()` | `WarningManifest()`, `WarningDefinitions()` | Advisory warnings: `key`, `label`, `phase` (`start` of a run or `end`), `setting` (registry path of its switch, or null), `default_on`, `message` template |
| `sk.example_manifest()` | `ExampleDeckManifest(manifest, decks_dir)` | Example decks: `file`, `name`, `title`, `family` (display label, "Perez (2019)"), `family_key` ("Perez2019"), `hidden`, `summary`, `description` |
| `sk.component_manifest()` | `ComponentManifest()`, `DescribeComponent(name)` | Components: `name` (CoolProp), `formula`, `molar_mass` [kg/mol], `cas`, `properties_source`, `aliases`, `ideal_gas_heat_capacity`, `heat_capacity_held`, `cp_species`, `cp_source` (provenance of the Cp polynomial, empty when held). The manifest is `{cp_reference, components}`: `cp_reference` names the table the polynomials come from |

Units: the parser accepts any LLNL-units symbol with an SI prefix, combined as described in
[Global Variables and Units](global-variables.md). So the manifest lists suggested units for each
category, not every accepted unit. Their factor and offset come from the parser itself. To
validate a unit someone typed, use `check_unit` rather than a list.

Example decks: `examples/example_decks.json` declares each deck's family label and key, hidden flag and summary.
Titles and descriptions come from the decks themselves. The manifest must list every deck in
`examples/decks/`; the test suite checks this.

Results series and profiles are listed by `sk.bulk_variables()` and its siblings. Model kinds are
listed by `sk.build_kinds()`.
