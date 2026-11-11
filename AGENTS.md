# AGENTS.md

This file gives instructions to coding agents that work **on** this repository. It also gives instructions to agents that **use** Skarstrom as a tool. This file uses ASD-STE100 Simplified Technical English. Humans: the [documentation](https://samaffleck.github.io/skarstrom/) contains the same data.

Skarstrom simulates pressure, vacuum and temperature swing adsorption (PSA/VSA/TSA) cycles. The engine is the C++20 library `adsim`. The `skarstrom` Python package and CLI wrap the engine.


<!-- doc-setup
import numpy as np
import skarstrom
import skarstrom as sk
sim = skarstrom.examples.Perez2019_A.load()
sim.add_variable("t_feed", 60, category="time")
sim.save("my_deck.json")
node, ads, step = sim.steps[0], sim.adsorbents["13X"], sim.steps[0]
name = "time"
path, mean, sd, lo, mode, hi = "layer[1].bed_void_fraction", 0.41, 0.01, 0.39, 0.41, 0.43
lower, upper, x, cand = [60], [150], [90], sim.copy()
p_pa = np.array([0.1, 0.25, 0.5, 1.0, 2.0, 4.0]) * 1e5
q_molkg = np.array([1.06, 1.94, 2.76, 3.47, 4.01, 4.36])
T_K = 298.15
-->

## Using Skarstrom

### Before you start (for agents)

- **Check the version, read the docs for it.**
  - Run `skarstrom --version` first.
  - The site root <https://samaffleck.github.io/skarstrom/> documents the latest release. This is the release that `pip install skarstrom` installs.
  - <https://samaffleck.github.io/skarstrom/dev/> documents `main`. These pages can describe parameters that the installed wheel does not accept.
- **Prefer the wheel.** Ask the user before you build from source. A source build takes several minutes. It must have CMake >= 3.15 and a C++20 compiler.
- **Run a few cycles before you run hundreds.** Start with `cycles=3`. After each run, examine `result.status`, `result.message`, `result.healthy` and `result.warnings`.
- **A run-health warning is a wrong answer until you can explain it.**
  - The status stays `completed` intentionally. Thus, you must examine the warnings yourself.
  - The codes are `MOLE_FRACTION_SUM`, `MASS_BALANCE_STUCK`, `STEP_RETRIES` and `REACTION_STIFF`.
  - Each warning has its value, threshold, suspects and a hint. Refer to the docs page `run-health.md`. The names and hints are in `core/include/adsim/solver/run_health.h`.
- **Do you compare with a different model or a paper?** First, read "Matching Another Model" (`docs/model-conventions.md`). It describes these items: bulk density and solid density, total voidage, the sign of the heat of adsorption, the composition that the rate uses, and the boundary profiles.
- **Machine-readable pointers:**
  - `skarstrom examples --json`.
  - `skarstrom run ... --json -q`.
  - The parameter registry (`sim.parameters()`, `skarstrom params deck.json --json`).
  - The deck view and the manifests (`docs/deck-view.md`).
  - `llms.txt` and `llms-full.txt` at the root of the docs site. `llms.txt` is an index of the pages. `llms-full.txt` is the full documentation as one Markdown file.

```python
result = sim.run(cycles=3)
print(result.status, result.message, result.healthy)
for w in result.warnings:
    print(w.code, w.value, w.threshold, w.hint)
```

### Install

<!-- no-run: installs the package -->
```bash
pip install skarstrom                    # wheels for CPython 3.10–3.15
# or, from a checkout (needs CMake >= 3.15 and a C++20 compiler; several minutes):
git clone --recursive https://github.com/samaffleck/skarstrom.git && cd skarstrom && pip install .
```

### Run a deck

A **deck** is a JSON simulation file. The package contains 47 example decks (`examples/decks/` in the repository).

The CLI can give machine-readable output. The exit codes are:

- 0: the run completed, or it reached cyclic steady state (CSS).
- 1: the run failed.
- 130: the run was interrupted.

<!-- no-run: builds from source -->
```bash
skarstrom examples --json                                   # [{name, family, path}, ...]
skarstrom examples run "Ward2022_Q100ccm" --json -q     # run a bundled deck
skarstrom run my_deck.json --cycles 20 --set "step[Feed].time=120" --json -q
skarstrom run deck.json --csv out/                          # kpis.csv, css.csv, run.json, logs/
```

The `--json` summary has these fields:

- `status`: `completed`, `reached_css`, `failed` or `cancelled`.
- `message`.
- `css_reached`.
- `css_cycle`: 0-based. The value is -1 if the run did not reach CSS.
- `cycles_run`.
- `wall_time_s`.
- `stats`: the integrator counters.
- `kpis`: one object for each tracked component in each cycle of the run. Undefined values are `null`.
- `css`: the CSS check for each cycle.
- `warnings`: the operating limits of the last cycle (wall design pressure and temperature, vacuum, approach to fluidisation; refer to `sim.operating_limits()`). The run-health messages are added at the end.
- `healthy`.
- `run_health`: the structured run-health warnings `{code, value, threshold, cycle, step, suspects, hint}`. The codes are MOLE_FRACTION_SUM, MASS_BALANCE_STUCK, STEP_RETRIES and REACTION_STIFF (refer to docs/run-health.md).

The warnings do not change the status or the exit code. The progress output goes to stderr. Use `-q` to stop the progress output. The run-health warnings continue to go to stderr.

Python:

```python
import skarstrom

sim = skarstrom.examples.Perez2019_A.load()   # or skarstrom.Simulation(path | JSON text | Example)
result = sim.run()                       # blocking; run(cycles=N, progress=f, log_level=..., on_message=f)
result.status, result.ok, result.css_cycle, result.message
result.healthy, result.warnings          # run-health checks (status unchanged); sim.run_health() per cycle
```

### Parameters: objects

The parts of a deck are typed objects over the parameter registry (`python/src/skarstrom/model.py`):

```python
from skarstrom import ureg as u
sim.add_variable("t_feed", 60 * u.s, category="time")   # a user variable, usable in expressions
sim.steps["Feed"].time = 2 * u.min             # steps by name or position; .by_id(n) by deck ID
sim.steps[0].bottom.temperature = 298.15     # step end -> its stream/boundary (end[name|i] if several)
sim.steps["Feed"].bottom.y["CO2"] = 0.15     # components: CoolProp names or formulas
sim.layers[0].length; sim.adsorbents["13X"].isotherm["CO2"].b1; sim.adsorbents[0].kinetics["CO2"]
sim.variables["t_feed"]; sim.cycle.number_of_cycles; sim.bed; sim.wall; sim.settings; sim.css
sim.energy_balance; sim.info
node.fields; node.get(name, unit="min"); node.set(name, 2, unit="min"); node.path(name); node.parameters()
```

Units:

- A bare number is in SI units. The docstring of each field and `Parameter.unit` give the unit.
- To use a different unit, give a pint Quantity (`2 * u.min`, `u.Quantity(25, "degC")`; `from skarstrom import ureg as u`). You can also use `node.set(name, 2, "min")`.
- A text value is an expression or a quantity in deck syntax (`"2 [min]"`). Decks and the CLI use text values. The engine parses the unit text with LLNL units.
- The engine checks the dimension. For example, `2 * u.bar` for a time causes a ValueError.
- `Simulation(deck, require_units=True)` does not accept bare numbers for dimensional fields.
- A misspelt or absent field causes an AttributeError. The package does not make new attributes.

`skarstrom.examples.<Name>` gives the bundled decks (`.load()`, `.path`, `.title`, `.copy_to()`).

### Parameters: paths (the registry)

Each settable value also has a stable **path** and an **alias**. The path uses IDs, and it does not change when you change a name. The alias uses the display names. Use paths for sweeps, optimisers and the CLI:

```python
from skarstrom import ureg as u
for p in sim.parameters(tunable_only=True, search="time"):
    print(p.path, p.alias, p.value, p.unit, p.kind, p.bounds)
sim["step[1].time"] = 45.0               # number -> value (SI units)
sim["step[1].time"] = 2 * u.min            # quantity -> converted, dimension-checked
sim["step[Feed].time"] = "2 * t_feed"    # string -> expression (may use variables and a trailing [unit])
sim.set("step[1].time", 2, unit="min"); sim.get("step[1].time", unit="min")
sim.get("step[1].time"); sim.get_text("step[1].time")
sim.update({"css.min_cycles": 3, "settings.reltol": 1e-6})
sim.save("modified.json")
```

CLI: `skarstrom params deck.json [--all] [--search TEXT] [--json]`. On `run`, use `--set PATH=VALUE` (you can use it more than one time).

The path families are:

- `cycle.*`
- `step[ID].time|fraction`
- `endpoint[ID].<field>` (flows, pressures, temperatures, `y[Component]`, pump settings)
- `bed.*`, `wall.*`, `layer[ID].*`
- `adsorbent[ID].<property>|isotherm[C].<param>|kinetics[C].<param>|heat_of_adsorption[C]`
- `variable[name]`
- `settings.*`, `css.*`, `energy_balance.*`, `info.*`
- `kpi[C].tracked`
- `cost.*` (the cost estimate)
- `sensitivity.*` (the settings of the stored study)

An unknown path causes a `KeyError`. The CLI stops with the message "unknown parameter".

Structural edits:

- `sim.edit.*`: the `deck_edit.h` cascades of the engine. They add, remove and move layers, steps and streams. They also do kind switches that keep the shared parameters, `ensure_ids`, `repair_connections` and `fill_gaps`. `skarstrom.model_kinds()` gives the kind descriptors.
- `sim.add_component(name)` and `sim.remove_component(name)`.
- `sim.add_library_adsorbent(key)`. The key is an entry key, or a material ID (for example `"zeolite-13x"`) for the recommended entry of that material. The function returns the adsorbent ID.
- The builder in the section below.

### Building a simulation in Python

The C++ code makes the structure from kind names (`core/include/adsim/objects/deck_builder.h`, Python `sim._sim.build(op, args)`). Then you set the values through the registry. `skarstrom.build_kinds()` gives the kinds of each family. An unknown kind causes a ValueError that gives the valid kinds.

<!-- no-run: abbreviated sketch with placeholders; docs/python.md has the complete, executed version -->
```python
sim = sk.Simulation.new(title="...")
sim.add_component("CarbonDioxide")
ads = sim.add_adsorbent("13X", solid_density=1138, heat_capacity=900, thermal_conductivity=0.63, particle_void_fraction=0.35)
ads.set_isotherm("extended_dual_site_langmuir_c", {"CO2": {"m1": ..., ...}, "N2": {...}})   # components given together are coupled
ads.set_rate_model("CO2", "ldf", "macropore", tortuosity=3, macropore_void_fraction=0.35); ads.heat_of_adsorption["CO2"] = 30850
sim.bed.set_shape("cylinder", diameter=0.032); sim.wall.set_material("ss316l"); sim.wall.set_boundary("zero_flux")
layer = sim.add_layer("13X", composition={"13X": 1.0}, momentum="kozeny_carman", length=0.3, cells=30, bed_void_fraction=0.41, particle_diameter=1.5e-3)
layer.set_dispersion("ruthven", molecular_diffusion_coefficient=1.6e-5)   # + set_heat_dispersion / set_internal_heat_transfer / set_external_heat_transfer
sim.initial_conditions.pressure = 1e5; sim.initial_conditions.y["N2"] = 1.0
feed = sim.add_step("Feed", time=99)                                     # both ends closed until set
feed.bottom.set_streams("inlet", [sk.Stream("Feed inlet", flow=sk.velocity(0.085), y={"CO2": 0.15, "N2": 0.85}, temperature=298)])
feed.top.set_pressure("outlet", 1e5)                                     # or sk.linear / sk.exponential / sk.log_linear(end); pump=sk.Pump(...)
step.top.set_streams("outlet", [sk.ValveStream("Vent", pressure=1e5, cv=10, xt=0.53)]); step.bottom.set_closed()
sim.connect(feed.top["LPP out"], other.top["LPP in"]); step.add_energy_source(layer, power=500)
sim.kpi.track("CO2", feed=[feed.bottom], product=[evac.bottom])
```

- **Flows:** `sk.velocity|volumetric_flow|molar_flow|mass_flow(value or profile)`.
- **Profiles:** `sk.constant|linear|exponential|log_linear(end, start=, time=)`, `sk.exponential(..., factor=100)`.
  - Pressure ramps start from the bed pressure. Do not give `start=` for them.
  - `time=` gives a ramp that is shorter than the step. The default is the step time.
- **Structure and scripts:** `sim.describe()` gives the structure. `sk.to_script(sim)` writes builder code for any deck.
- **Example scripts:** each example script `examples/<family>.py` builds its cases from the start (`build(case)`). `python/tests/test_examples_build.py` makes sure that each script makes the same deck as the deck in `examples/decks/`.
- **User variables:** `sim.add_variable("t_feed", 60 * u.s, category="time")`, then `time="t_feed"`. `formula="t_feed / 3"` makes a variable from other variables.
- **Pump-driven streams:** `sk.PumpStream(name, power=, pump=)`. The power is a constant or `sk.polynomial_power(a, b, ...)`.
- **Buffer tanks:** `inlet.set_buffer_tank(source, percentage=)`.
- **Measured data in the deck:** `sim.experimental_data` (`add`, `load_csv`, `[name].compare(x, y)`).
- The builder does not make cost data.

### Results

```python
sim.kpis(cycle=-1)        # KPIRecord: component, cycle, valid, purity, recovery, productivity_mass,
                          #            productivity_volume, energy_consumption, production_rate,
                          #            electrical_energy (= energy_consumption: compressor/pump work),
                          #            thermal_energy (heat supplied), cooling_duty (heat removed, >= 0),
                          #            step_thermal_energy, step_cooling_duty (per step) (nan if undefined)
sim.kpi_table()           # pandas DataFrame
sim.css_history()         # CSSRecord per cycle: converged, error, tolerance, method
sim.series("mole_fraction", "CarbonDioxide", location="top", cycle=-1)   # .time, .values (NumPy), .unit
sim.series("temperature", z=0.15)                                        # at a position [m]
sim.profile("pressure", cycle=-1, step=0)  # .z, .time, .values with shape (len(time), len(z))
for s in sim.streams(): f = sim.stream_series(s["step"], s["name"], "molar_flow_rate")   # per stream: step, name, endpoint_id, location
skarstrom.bulk_variables(); skarstrom.component_variables(); skarstrom.stream_variables()
```

- KPIs are available only for the components that the KPI configuration of the deck tracks.
- Some decks have no tracked components: Casas, Ferré, Nguyen 2022, Streb, Vysyaraju, Ward 2022 and FiveStepProductPurge. For these decks, integrate `stream_series` (refer to `examples/_common.py:endpoint_moles`).
- The engine logs only the last `settings.max_cycles_to_keep` cycles (default 1). This applies to profiles and series. The engine keeps the KPIs and the CSS history for all cycles.

### Sensitivity and error bars

<!-- long: 50 + 100 or more cycle simulations -->
```python
from skarstrom import sensitivity as sa
study = sa.Study(sim, [sa.relative("step[1].time", 0.1), sa.uniform("layer[1].bed_void_fraction", 0.39, 0.43),
                       sa.normal(path, mean, sd, lower=0), sa.triangular(path, lo, mode, hi)],
                 outputs=None, settings={"css.state_tolerance": 5e-6}, require_css=True, workers=4)
screen = study.screen(trajectories=10)            # Morris: .table() (mu*, mu, sigma), .ranking(), .important()
spread = study.propagate(inputs=screen.important())   # LHS Monte Carlo: .table() (p05/p50/p95, converged), .shares()
```

- **CLI:** `skarstrom sensitivity deck.json --input "PATH=10%" --input "PATH=LO:HI" --input "PATH=normal:MEAN:SD" [--output Component.kpi] [--method screen|propagate|both] [--workers N] [--json] [--save FILE]`.
- **Outputs:** `Component.kpi` at the last cycle, for any KPI of `sa.kpi_names()`. The default is purity, recovery, productivity_mass and energy_consumption of each tracked component.
- **CSS:** each run must reach CSS, unless you set `require_css=False` or `--no-require-css`. Some decks have `css.method` set to none (for example the Arora decks). For these decks, use `settings={"css.method": 3}` or the flag.
- **Rejected runs:** `.failed` and `.warnings` count the runs that the study does not accept. The study does not remove them silently.
- **C++:** `adsim/solver/sensitivity.h` (plans, `EvaluateSample`, analyses).
- **Guide:** `docs/sensitivity.md`.

### Optimisation and warm starts

<!-- long: 20 cycle simulations to the optimum -->
```python
from skarstrom import optimisation as opt
r = opt.optimise(sim, {"step[Feed].time": (60, 150)}, objectives=[opt.maximise("CarbonDioxide", "recovery")],
                 constraints=[opt.at_least("CarbonDioxide", "purity", 0.9, penalty=10)], evaluations=20)
r.x, r.cost, r.search_cost, r.kpis, r.table()   # cost = cold re-run of the optimum; search_cost = warm search value
cache = opt.StateCache(lower, upper); run, info = opt.run_warm(cand, cache, x)   # your own loop
```

- Candidates and study runs do a warm start from the CSS of the nearest earlier point.
- To stop warm starts, use `warm_start=False` or `--no-warm-start`.
- The cache keeps only CSS states.
- C++: `adsim/solver/state_cache.h`, and `OptimizeSimulation` in `adsim/solver/optimisation.h`.

### Isotherms and the adsorbent library

```python
skarstrom.isotherm_kinds(); skarstrom.isotherm_parameters("sips")
skarstrom.evaluate_isotherm("sips", {"q_sat": 5, "b": 1e-5, "v": 0.8}, pressure=p_pa, temperature=298.15)
skarstrom.evaluate_iast(["sips", "toth"], [{"q_sat": 5, "b": 1e-5, "v": 0.8}, {"q_sat": 4, "b": 2e-5, "v": 0.6}], pressure=2e5, temperature=298.15, mole_fraction=[0.3, 0.7])
ads.set_mixing_rule("iast")   # per adsorbent: "none" (default) or "iast" over the components' pure isotherms
skarstrom.fit_isotherm("dual_site_langmuir_p", p_pa, q_molkg, temperature=T_K)   # params, std_errors, r2, rmse, ...
skarstrom.adsorbent_library()["materials"]; skarstrom.adsorbent_library()["entries"]; skarstrom.adsorbent_path("zeolite-13x")
```

All values are in SI base units (Pa, K, mol/kg, J/mol, m, s).

### Deck schema and versioning

- Decks are cereal JSON archives of `adsim::UserData`. The root is `value0` for files and `simulation` for string saves. The engine can load the two types.
- Each serialised class has a `cereal_class_version`. For older versions, the load function does a migration. Thus old decks continue to load. But an older engine possibly cannot load a deck from a newer engine.
- Do not edit a structure manually if you do not understand it. Change values through the parameter registry.
- The package version is SemVer (`skarstrom.__version__`, `adsim::kVersionString`).

## Working on the repository

### Layout

```
core/include/adsim/   public C++ API (source of truth): simulation.h, objects/, equations/, solver/, parameters/
core/src/             implementations (mirrors include/)
core/test/            GoogleTest: unit/, integration/, validation/ (analytic verification), papers/ (literature)
core/data/            digitised experimental curves used by the paper tests
core/bench/           bench_solver timing tool
core/install_consumer/, core/subproject_consumer/   find_package / add_subdirectory smoke projects
python/                pybind11 bindings (src/bindings.cpp), the skarstrom package, pytest tests
examples/              decks/ (bundled decks), one script per deck family, notebooks/
data/adsorbents/       curated adsorbent library: library.json + <key>.json
data/cycles/           cycle templates (written by test_cycle_designs.cpp)
docs/                  MkDocs sources (validation.md and adsorbent-library.md are GENERATED)
scripts/               docs generators (scripts/docs/), release and consumer test scripts
external/              submodules (sundials, cereal, muparser, dlib, units, googletest) + vendored eigen-3.4.0
```

### Build and test

<!-- no-run: builds and tests the engine -->
```bash
git submodule update --init --recursive
cmake --preset linux-release            # macos-release, windows; tests are on in the presets
cmake --build --preset linux-release
ctest --preset linux-release -j3        # full suite ~10 min; core/ is the working directory
ctest --preset linux-release -j3 -LE "paper|slow"   # fast tier (what pull-request CI runs), ~1 min
ctest --preset linux-release -L slow    # slow verification tier (main and nightly CI)
pip install . && pytest python/tests    # Python package + example smoke tests
pip install -r requirements-docs.txt && mkdocs build --strict
```

### Conventions

- **Mark new names "New in X.Y".**
  - Each new kind, parameter, Python argument, CLI option and series variable in `docs/` must have a `!!! note "New in X.Y"` admonition. If a table row is the only location, use an inline `*New in X.Y.*`.
  - X.Y is the first release that has the name. This is the NEXT release: the version that the daily release workflow (`daily-release.yml`) selects. The workflow selects a minor version if `changelog.d/` contains an Added or Changed fragment.
  - Remove the note when that release is two minor versions old.
  - The site root shows the latest release, and `/dev/` shows `main`. The note tells the reader of each site which wheel has the name.
- **Never commit `core/Results/`** (test output, git-ignored) or build trees.
- **Keep generated docs generated.**
  - `python scripts/docs/gen_validation_page.py` makes `docs/validation.md` from `core/test/papers/validation_summary.json`. After you run the paper tests, add `--merge <build>/validation`.
  - `python scripts/docs/gen_adsorbent_library.py` makes `docs/adsorbent-library.md` from `data/adsorbents/library.json`. This script also validates the catalogue.
  - Do not edit the generated pages. CI runs the two scripts with `--check`.
- **Cite every adsorbent parameter.** A library entry must have a DOI and a table or figure reference for each value. If the source does not give a value, write `Not stated in source (checked: ...)`.
- **One recommended entry per material.**
  - Each library entry belongs to a material (`materials` in `library.json`).
  - Each material has exactly one `recommended` entry. Select it on data quality. The `data_quality` of each isotherm has this order: measured with tabulated data > measured, fit only > refit of the data of others > copied from an earlier paper.
  - Add a second entry for a material only if it is a different product, form or data range that a user can name (its `variant`). If not, replace the weaker entry and put its key in the `aliases` of the material.
  - A label is the material name and the variant in brackets. Do not use author names.
  - Paper reference values (`reference_values.csv`) must be numbers that the paper prints. Do not read values from a plot by eye. Do not commit PDFs or full-text extractions.
- **Documentation code blocks run in CI.**
  - `python scripts/docs/run_doc_blocks.py` runs each ```` ```python ````, ```` ```bash ```` and ```` ```console ```` block in `docs/*.md` and `AGENTS.md` against the built package. CI runs it in the `Python wheel` job. Locally, run it after `pip install .`. Thus a doc cannot show a call that the engine does not accept.
  - The Python blocks of a page share one interpreter and one new temporary directory. They run in page order.
  - Each shell block is a separate bash process in that directory. The shell blocks run after the Python blocks.
  - The script limits runs to 3 cycles (`--cycles N`).
  - Some blocks are only fragments of a snippet. A hidden `<!-- doc-setup` ... `-->` comment of Python near the top of the page gives them what they use (`sim`, `ads`, a deck file, data arrays). This comment runs first.
  - Put markers on the line before the fence:
    - `<!-- no-run: reason -->`: for a block that must not run (installs, builds, placeholders). For a known engine bug, write `<!-- no-run: engine bug #N -->`.
    - `<!-- long: reason -->`: for a block that is too slow for pull requests. `--full` still runs it. `--full` is the run without the cycle limit, as at night, and it must have `pymoo`.
  - Do not edit a generated page (`validation.md`, `adsorbent-library.md`) to add a marker. Change its generator.
  - `--list` shows the blocks that run.
- **Keep the fast tier fast.**
  - Pull-request CI runs `ctest -LE "paper|slow"` on three platforms. A test in this tier must answer only two questions: does it run, and does the feature work?
  - Use approximately 20-30 cells, a maximum of approximately 30 cycles and short steps. Set CSS off, unless the test examines CSS.
  - Some tests compare with a closed-form solution, verify a convergence order or examine CSS behaviour. Only these tests can use the grid or the number of cycles that their assertion must have.
  - If such a test takes more than a few seconds, add it to `ADSIM_SLOW_TESTS` in `core/test/CMakeLists.txt` (ctest label `slow`, runs on main and at night). Write the reason in the comment of the test.
  - Put a new paper suite in `ADSIM_PAPER_TESTS` in the same file (label `paper`).
- **Tolerances in paper tests are regression bounds.** The measured error sets them. Do not make a tolerance larger to hide a disagreement. Write the disagreement in the README of the paper.
- **The test suite re-saves the example decks. The files stay equal to its output.**
  - The tests write `examples/decks/*.json`, some `data/adsorbents/*.json` and `data/cycles/*.json`. We commit these files, so that the decks that ship are the decks that the tests run.
  - The save is deterministic: the maps are written in key order and the decks are written in the form that a load gives (`test::SaveDeckAtomically`). Thus a second run of the suite makes no change.
  - After a change that moves the numbers in a deck, run the tests and commit the changed files with the change. CI fails if the tests change a file in `examples/` or `data/`.
  - A change of the order of entries in a file is never a change of a result. Do not change a tolerance or a value by hand to remove a diff.
- **Data file names must be shell-safe** (`[A-Za-z0-9._,=+@-]`). The `DataFileNamesAreShellSafe` test examines them.
- **Bump `cereal_class_version`** for each change to a serialised class. Make the load function accept the older versions.
- **The engine stays front-end free.**
  - Do not put GUI, web or JSON-schema code in `core/`.
  - The engine must continue to build as a subproject and with Emscripten (`scripts/test_subproject_consumer.sh [--emscripten]`).
- Licence headers refer to the root `LICENSE`. Do not add licence text to each file. Do not change third-party code in `external/`.
- Each commit must have a DCO sign-off (`git commit -s`). Refer to CONTRIBUTING.md.
- Write this file in ASD-STE100 Simplified Technical English. Use short sentences, the active voice and one instruction in each sentence. Do not change code, commands, paths or names to obey these rules.
