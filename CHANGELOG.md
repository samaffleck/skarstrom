# Changelog

All notable changes to Skarstrom are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/). The single source of the version is
`ADSIM_PROJECT_VERSION` (plus `ADSIM_VERSION_SUFFIX` for pre-releases) in the
root `CMakeLists.txt`.

## [Unreleased]

### Changed
- Contact is through LinkedIn instead of email (README "Get involved",
  SECURITY.md, CODE_OF_CONDUCT.md); the README and docs gain the logo, a
  "Why Skarstrom" introduction and a call for steering-committee members.
- README: a "Browser interface" section introducing the Skarstrom website
  (https://www.skarstrom.co.uk/), the browser front end built on this engine.

### Removed
- `adsim/tools/eigenlldb.py`, Eigen's LLDB pretty-printer (unused).
- Example decks renamed so every name is a Python identifier:
  `Ward_2022_Q=50[ccm]` → `Ward2022_Q50ccm` (and `Q100ccm`, `Q200ccm`),
  `Nguyen_2022_Fig4`/`Fig5` → `Nguyen2022_Fig4`/`Fig5`, `streb_2021_A1`/`A3`
  → `Streb2021_A1`/`A3`. The shell-safe-name allowlist is now empty.
- Python: `Simulation.steps` is a collection of step objects (names:
  `sim.steps.names`); `skarstrom.examples()`, `example_path()` and
  `copy_example()` are replaced by the `skarstrom.examples` module; the CLI's
  `examples` command lists names without `.json`.

### Fixed
- A value with an inline unit of the wrong dimension (`step time = "2 [bar]"`)
  was converted silently (to 200 000 s); it is now rejected.
- Competitive isotherms summed their per-component terms in `std::unordered_map`
  order, which depends on the order of the entries in the deck file and on the
  standard library. A deck and its re-saved copy could therefore run to
  slightly different numbers (seen on macOS arm64, where FMA contraction makes
  even a two-component sum order-dependent). Isotherm parameters are now kept
  in key order (`std::map`); the deck format is unchanged.
- CI: the "example decks unchanged" check no longer uses `A && B || C`
  (shellcheck SC2015) and also fails when the test suite deletes a file.

### Added
- Python: the deck as typed objects (`skarstrom.model`). `sim.steps["Feed"].time`,
  `sim.steps[0].bottom.temperature`, `sim.steps["Feed"].bottom.y["CO2"]`,
  `sim.layers[0].length`, `sim.adsorbents["13X"].isotherm["CO2"].b1`,
  `sim.variables[...]`, `sim.cycle`, `sim.bed`, `sim.wall`, `sim.settings`,
  `sim.css`, `sim.energy_balance`, `sim.info`. They are views onto the
  parameter registry (string paths keep working); misspelt or absent fields
  raise `AttributeError` with suggestions; editors complete the fields and
  show their SI units. A test checks every registry parameter of every
  example deck is reachable this way.
- Python: `skarstrom.examples`, one attribute per bundled deck
  (`skarstrom.examples.Perez2019_A.load()`, `.path`, `.title`, `.copy_to()`),
  generated from `examples/decks/` by `scripts/gen_examples_module.py`
  (checked in CI).
- Units. A bare number is SI, as before; a value can now carry its unit as a
  string (`"2 min"`, `"25 degC"`, `"35 kJ/mol"`, also `"2 [min]"`), as a pint
  Quantity, or with `unit=`; `get(..., unit=...)` reads in any unit;
  `Simulation(deck, require_units=True)` rejects bare numbers for dimensional
  values; `skarstrom.units.convert()`. In the engine: `ParseUnitExpression`
  reads compound units (`W/(m^2*K)`, `W/m2-K`, `kJ/(kg-K)`, `1/bar`,
  `Nm3/h`...), `ConvertUnit` and `SplitQuantity`; `Parameter::SetText` (and
  so the CLI's `--set`) takes `"<number> <unit>"` and checks the unit's
  dimension against the parameter's. `ParseUnitSuffix` falls back to unit
  expressions for symbols missing from its dictionary.
- `DeckRoundTrip` tests: every example deck's numbers parse exactly and
  survive save/reload bit for bit, and a run is unchanged when every map in
  the deck is listed in reverse order.

## [0.1.0] — initial open-source release (unreleased)

First public release of the Skarstrom engine under the Apache License 2.0.

### Added
- The `adsim` C++20 simulation engine: 1D finite-volume packed-bed model
  (mass, momentum, energy and wall energy balances; 20 isotherm models;
  linear driving force kinetics with several coefficient correlations;
  dispersion and heat-transfer correlations; CoolProp fluid properties),
  multi-step and multi-bed cycles with valves, pumps, buffer tanks and step
  connections, SUNDIALS IDA integration, cyclic-steady-state detection and
  Aitken acceleration, KPIs, a parameter registry and a dlib optimiser.
- Public headers for front ends: `adsim/objects/deck_io.h` (load decks with
  either cereal root), `adsim/objects/deck_edit.h` (component, adsorbent and
  connection edits that keep a deck consistent),
  `adsim/equations/isotherm_kinds.h` (isotherm kind tags, parameter access,
  evaluation) and `adsim/solver/isotherm_fit.h` (the isotherm fitter).
- CMake build with presets for Linux, macOS and Windows; `find_package(Adsim)`
  install/export; use as a subproject (`add_subdirectory`) with front-end
  friendly defaults; builds with Emscripten.
- GoogleTest suite (~670 tests): unit, integration, analytic verification and
  quantitative validation against ten published studies.
- The `skarstrom` Python package (scikit-build-core + pybind11) and CLI
  (`run`, `params`, `examples list|run|path|copy`, JSON output), bundling all
  example decks and the adsorbent library; structural edits from Python
  (`add_component`, `remove_component`, `add_library_adsorbent`).
- 33 example decks, a runnable Python script per deck family in `examples/`,
  and three executed notebooks.
- A curated, cited adsorbent library (`data/adsorbents/`) and cycle templates
  (`data/cycles/`).
- MkDocs documentation with generated validation and adsorbent-library pages.
- GitHub Actions: CI on Linux, macOS and Windows, the Python wheel, the
  Emscripten build and the docs; release workflow with cibuildwheel wheels,
  sdist, source tarball and gated (Test)PyPI trusted publishing; docs
  deployment to GitHub Pages.
