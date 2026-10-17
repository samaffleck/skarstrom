# AGENTS.md

Guidance for coding agents working **on** this repository and for agents **using** Skarstrom as a tool. Humans: the same facts are in the [documentation](https://samaffleck.github.io/skarstrom/).

Skarstrom simulates pressure, vacuum and temperature swing adsorption (PSA/VSA/TSA) cycles. The engine is the C++20 library `adsim`; the `skarstrom` Python package and CLI wrap it.

## Using Skarstrom

### Install

```bash
pip install skarstrom                    # once published on PyPI
# or, from a checkout (needs CMake >= 3.15 and a C++20 compiler; several minutes):
git clone --recursive https://github.com/samaffleck/skarstrom.git && cd skarstrom && pip install .
```

### Run a deck

A **deck** is a JSON simulation file. 33 example decks ship with the package (`examples/decks/` in the repository).

CLI with machine-readable output (exit code 0 = completed or cyclic steady state reached, 1 = failed, 130 = interrupted):

```bash
skarstrom examples --json                                   # [{name, family, path}, ...]
skarstrom examples run "Ward_2022_Q=100[ccm]" --json -q     # run a bundled deck
skarstrom run path/to/deck.json --cycles 20 --set "step[Feed].time=120" --json -q
skarstrom run deck.json --csv out/                          # kpis.csv, css.csv, run.json, logs/
```

The `--json` summary has `status` (`completed`, `reached_css`, `failed`, `cancelled`), `message`, `css_reached`, `css_cycle` (0-based, -1 if not reached), `cycles_run`, `wall_time_s`, `stats` (integrator counters), `kpis` (one object per tracked component per logged cycle; undefined values are `null`) and `css` (the CSS check per cycle). Progress goes to stderr; use `-q` to silence it.

Python:

```python
import skarstrom

sim = skarstrom.Simulation(skarstrom.example_path("Perez2019_A.json"))   # or a path / JSON text
result = sim.run()                       # blocking; run(cycles=N, progress=f, log_level=..., on_message=f)
result.status, result.ok, result.css_cycle, result.message
```

### Parameters (the registry)

Every settable value has a stable **path** (ID-based, survives renames) and an **alias** (display names):

```python
for p in sim.parameters(tunable_only=True, search="time"):
    print(p.path, p.alias, p.value, p.unit, p.kind, p.bounds)
sim["step[1].time"] = 45.0               # number -> value (SI units)
sim["step[Feed].time"] = "2 * t_feed"    # string -> expression (may use variables and units, "1.5 [bar]")
sim.get("step[1].time"); sim.get_text("step[1].time")
sim.update({"css.min_cycles": 3, "settings.reltol": 1e-6})
sim.save("modified.json")
```

CLI: `skarstrom params deck.json [--all] [--search TEXT] [--json]`; `--set PATH=VALUE` on `run` (repeatable). Path families: `cycle.*`, `step[ID].time|fraction`, `endpoint[ID].<field>` (flows, pressures, temperatures, `y[Component]`, pump settings), `bed.*`, `wall.*`, `layer[ID].*`, `adsorbent[ID].<property>|isotherm[C].<param>|kinetics[C].<param>|heat_of_adsorption[C]`, `variable[name]`, `settings.*`, `css.*`, `energy_balance.*`, `info.*`. Unknown paths raise `KeyError` (CLI: exit with "unknown parameter").

Structural edits: `sim.add_component(name)`, `sim.remove_component(name)`, `sim.add_library_adsorbent(key)` (returns the adsorbent ID).

### Results

```python
sim.kpis(cycle=-1)        # KPIRecord: component, cycle, valid, purity, recovery, productivity_mass,
                          #            productivity_volume, energy_consumption, production_rate (nan if undefined)
sim.kpi_table()           # pandas DataFrame
sim.css_history()         # CSSRecord per cycle: converged, error, tolerance, method
sim.series("mole_fraction", "CarbonDioxide", location="top", cycle=-1)   # .time, .values (NumPy), .unit
sim.series("temperature", z=0.15)                                        # at a position [m]
sim.profile("pressure", cycle=-1, step=0)  # .z, .time, .values with shape (len(time), len(z))
sim.streams(); sim.stream_series(step, name, "molar_flow", "CarbonDioxide")
skarstrom.bulk_variables(); skarstrom.component_variables(); skarstrom.stream_variables()
```

KPIs exist only for components tracked in the deck's KPI configuration; some decks have none (Casas, Ferré, Nguyen 2022, Streb, Vysyaraju, Ward 2022, FiveStepProductPurge) — integrate `stream_series` instead (see `examples/_common.py:endpoint_moles`). Only the last `settings.max_cycles_to_keep` (default 10) cycles are logged.

### Isotherms and the adsorbent library

```python
skarstrom.isotherm_kinds(); skarstrom.isotherm_parameters("sips")
skarstrom.evaluate_isotherm("sips", {"q_sat": 5, "b": 1e-5, "v": 0.8}, pressure=p_pa, temperature=298.15)
skarstrom.fit_isotherm("dual_site_langmuir_p", p_pa, q_molkg, temperature=T_K)   # params, std_errors, r2, rmse, ...
skarstrom.adsorbent_library()["entries"]; skarstrom.adsorbent_path("zeolite_13x_perez2019")
```

All values are SI base units (Pa, K, mol/kg, J/mol, m, s).

### Deck schema and versioning

Decks are cereal JSON archives of `adsim::UserData` (root `value0` for files, `simulation` for string saves; both load). Every serialised class carries a `cereal_class_version`; loading runs the migration for older versions, so old decks keep loading, but a deck written by a newer engine may not load in an older one. Do not hand-edit structure you do not understand — change values through the parameter registry. The package version is SemVer (`skarstrom.__version__`, `adsim::kVersionString`).

## Working on the repository

### Layout

```
adsim/include/adsim/   public C++ API (source of truth): simulation.h, objects/, equations/, solver/, parameters/
adsim/src/             implementations (mirrors include/)
adsim/test/            GoogleTest: unit/, integration/, validation/ (analytic verification), papers/ (literature)
adsim/data/            digitised experimental curves used by the paper tests
adsim/bench/           bench_solver timing tool
adsim/install_consumer/, adsim/subproject_consumer/   find_package / add_subdirectory smoke projects
python/                pybind11 bindings (src/bindings.cpp), the skarstrom package, pytest tests
examples/              decks/ (bundled decks), one script per deck family, notebooks/
data/adsorbents/       curated adsorbent library: library.json + <key>.json
data/cycles/           cycle templates (written by test_cycle_designs.cpp)
docs/                  MkDocs sources (validation.md and adsorbent-library.md are GENERATED)
scripts/               docs generators (scripts/docs/), release and consumer test scripts
external/              submodules (sundials, CoolProp, cereal, muparser, dlib, googletest) + vendored eigen-3.4.0
```

### Build and test

```bash
git submodule update --init --recursive
cmake --preset linux-release            # macos-release, windows; tests are on in the presets
cmake --build --preset linux-release
ctest --preset linux-release -j3        # full suite ~10 min; adsim/ is the working directory
ctest --preset linux-release -E "^(Arora2021|Casas2012|Ferre2026|Nguyen2022|Nguyen2023|Perez2019|Streb2021|Vysyaraju2024|Ward2022|WardPini2024)Test\.|^CSSAcceleration\."   # fast subset
rm -rf adsim/Results                    # the suite writes ~4 GB there
pip install . && pytest python/tests    # Python package + example smoke tests
pip install -r requirements-docs.txt && mkdocs build --strict
```

### Conventions

- **Never commit `adsim/Results/`** (test output, git-ignored) or build trees.
- **Keep generated docs generated.** `docs/validation.md` comes from `adsim/test/papers/validation_summary.json` (`python scripts/docs/gen_validation_page.py`, `--merge <build>/validation` after running the paper tests); `docs/adsorbent-library.md` from `data/adsorbents/library.json` (`python scripts/docs/gen_adsorbent_library.py`, which also validates the catalogue). Never edit the generated pages; CI runs both with `--check`.
- **Cite every adsorbent parameter.** A library entry needs a DOI and table/figure reference for each value, or an explicit `Not stated in source (checked: ...)`. Paper reference values (`reference_values.csv`) must be numbers printed in the paper; never eyeball values off a plot, never commit PDFs or full-text extractions.
- **Tolerances in paper tests are regression bounds** set from the measured error; never loosen one to hide a disagreement — document it in the paper's README.
- **Regenerated example decks are expected.** The test suite re-saves `examples/decks/*.json` and some `data/adsorbents/*.json`; embedded result numbers move in the last digits. Commit those changes with the code change that caused them.
- **Data file names must be shell-safe** (`[A-Za-z0-9._,=+@-]`; checked by the `DataFileNamesAreShellSafe` test).
- **Bump `cereal_class_version`** for any change to a serialised class and handle older versions on load.
- **The engine stays front-end free:** no GUI, web or JSON-schema code in `adsim/`; it must keep building as a subproject and with Emscripten (`scripts/test_subproject_consumer.sh [--emscripten]`).
- Licence headers point to the root `LICENSE`; do not add per-file licence text. Leave third-party code in `external/` untouched.
- Commits need a DCO sign-off (`git commit -s`); see CONTRIBUTING.md.
