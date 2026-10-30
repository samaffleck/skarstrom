# AGENTS.md

Guidance for coding agents working **on** this repository and for agents **using** Skarstrom as a tool. Humans: the same facts are in the [documentation](https://samaffleck.github.io/skarstrom/).

Skarstrom simulates pressure, vacuum and temperature swing adsorption (PSA/VSA/TSA) cycles. The engine is the C++20 library `adsim`; the `skarstrom` Python package and CLI wrap it.

## Using Skarstrom

### Install

```bash
pip install skarstrom                    # wheels for CPython 3.10–3.15
# or, from a checkout (needs CMake >= 3.15 and a C++20 compiler; several minutes):
git clone --recursive https://github.com/samaffleck/skarstrom.git && cd skarstrom && pip install .
```

### Run a deck

A **deck** is a JSON simulation file. 33 example decks ship with the package (`examples/decks/` in the repository).

CLI with machine-readable output (exit code 0 = completed or cyclic steady state reached, 1 = failed, 130 = interrupted):

```bash
skarstrom examples --json                                   # [{name, family, path}, ...]
skarstrom examples run "Ward2022_Q100ccm" --json -q     # run a bundled deck
skarstrom run path/to/deck.json --cycles 20 --set "step[Feed].time=120" --json -q
skarstrom run deck.json --csv out/                          # kpis.csv, css.csv, run.json, logs/
```

The `--json` summary has `status` (`completed`, `reached_css`, `failed`, `cancelled`), `message`, `css_reached`, `css_cycle` (0-based, -1 if not reached), `cycles_run`, `wall_time_s`, `stats` (integrator counters), `kpis` (one object per tracked component per logged cycle; undefined values are `null`) `css` (the CSS check per cycle) and `warnings` (operating limits of the last cycle: wall design pressure and temperature, vacuum, approach to fluidisation; see `sim.operating_limits()`). Progress goes to stderr; use `-q` to silence it.

Python:

```python
import skarstrom

sim = skarstrom.examples.Perez2019_A.load()   # or skarstrom.Simulation(path | JSON text | Example)
result = sim.run()                       # blocking; run(cycles=N, progress=f, log_level=..., on_message=f)
result.status, result.ok, result.css_cycle, result.message
```

### Parameters: objects

The deck's parts are typed objects over the parameter registry (`python/src/skarstrom/model.py`):

```python
from skarstrom import ureg as u
sim.steps["Feed"].time = 2 * u.min             # steps by name or position; .by_id(n) by deck ID
sim.steps[0].bottom.temperature = 298.15     # step end -> its stream/boundary (end[name|i] if several)
sim.steps["Feed"].bottom.y["CO2"] = 0.15     # components: CoolProp names or formulas
sim.layers[0].length; sim.adsorbents["13X"].isotherm["CO2"].b1; sim.adsorbents[0].kinetics["CO2"]
sim.variables["t_feed"]; sim.cycle.number_of_cycles; sim.bed; sim.wall; sim.settings; sim.css
sim.energy_balance; sim.info
node.fields; node.get(name, unit="min"); node.set(name, 2, unit="min"); node.path(name); node.parameters()
```

Units: a bare number is SI (the unit is in each field's docstring and `Parameter.unit`). Give another
unit with a pint Quantity (`2 * u.min`, `u.Quantity(25, "degC")`; `from skarstrom import ureg as u`) or
`node.set(name, 2, "min")`. Text is an expression or a deck-syntax quantity (`"2 [min]"`), for decks and
the CLI; unit text is parsed by the engine with LLNL units. The dimension is checked (`2 * u.bar` for a
time raises ValueError). `Simulation(deck, require_units=True)` rejects bare numbers for dimensional fields. Misspelt or absent fields raise AttributeError (never create attributes).
`skarstrom.examples.<Name>` lists the bundled decks (`.load()`, `.path`, `.title`, `.copy_to()`).

### Parameters: paths (the registry)

Every settable value also has a stable **path** (ID-based, survives renames) and an **alias** (display names);
use paths for sweeps, optimisers and the CLI:

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

CLI: `skarstrom params deck.json [--all] [--search TEXT] [--json]`; `--set PATH=VALUE` on `run` (repeatable). Path families: `cycle.*`, `step[ID].time|fraction`, `endpoint[ID].<field>` (flows, pressures, temperatures, `y[Component]`, pump settings), `bed.*`, `wall.*`, `layer[ID].*`, `adsorbent[ID].<property>|isotherm[C].<param>|kinetics[C].<param>|heat_of_adsorption[C]`, `variable[name]`, `settings.*`, `css.*`, `energy_balance.*`, `info.*`. Unknown paths raise `KeyError` (CLI: exit with "unknown parameter").

Structural edits: `sim.add_component(name)`, `sim.remove_component(name)`, `sim.add_library_adsorbent(key)` (returns the adsorbent ID), and the builder below.

### Building a simulation in Python

Structure is created by kind name in C++ (`core/include/adsim/objects/deck_builder.h`, Python `sim._sim.build(op, args)`); values then go through the registry. `skarstrom.build_kinds()` lists every family's kinds; unknown kinds raise ValueError listing the valid ones.

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

Flows: `sk.velocity|volumetric_flow|molar_flow|mass_flow(value or profile)`; profiles `sk.constant|linear|exponential|log_linear(end, start=, time=)` (pressure ramps take only the end). `sim.describe()` gives the structure; `sk.to_script(sim)` writes builder code for any deck. Every example script `examples/<family>.py` builds its cases from scratch (`build(case)`); `python/tests/test_examples_build.py` checks each reproduces its deck in `examples/decks/`. User variables: `sim.add_variable("t_feed", 60 * u.s, category="time")`, then `time="t_feed"`; `formula="t_feed / 3"` defines one from others. Pump-driven streams: `sk.PumpStream(name, power=, pump=)` (power constant or `sk.polynomial_power(a, b, ...)`). Buffer tanks: `inlet.set_buffer_tank(source, percentage=)`. Measured data in the deck: `sim.experimental_data` (`add`, `load_csv`, `[name].compare(x, y)`). Not created by the builder: cost data.

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

### Sensitivity and error bars

```python
from skarstrom import sensitivity as sa
study = sa.Study(sim, [sa.relative("step[1].time", 0.1), sa.uniform("layer[1].bed_void_fraction", 0.39, 0.43),
                       sa.normal(path, mean, sd, lower=0), sa.triangular(path, lo, mode, hi)],
                 outputs=None, settings={"css.mass_balance_tolerance": 2e-3}, require_css=True, workers=4)
screen = study.screen(trajectories=10)            # Morris: .table() (mu*, mu, sigma), .ranking(), .important()
spread = study.propagate(inputs=screen.important())   # LHS Monte Carlo: .table() (p05/p50/p95, converged), .shares()
```

CLI: `skarstrom sensitivity deck.json --input "PATH=10%" --input "PATH=LO:HI" --input "PATH=normal:MEAN:SD" [--output Component.kpi] [--method screen|propagate|both] [--workers N] [--json] [--save FILE]`. Outputs are `Component.kpi` at the last cycle (default: purity, recovery, productivity_mass, energy_consumption of every tracked component). Runs must reach CSS unless `require_css=False` / `--no-require-css`; a deck whose `css.method` is none (e.g. the Arora decks) needs `settings={"css.method": 3}` or the flag. Not-accepted runs are counted in `.failed` / `.warnings`, never dropped silently. C++: `adsim/solver/sensitivity.h` (plans, `EvaluateSample`, analyses). Guide: `docs/sensitivity.md`.

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
external/              submodules (sundials, CoolProp, cereal, muparser, dlib, googletest) + vendored eigen-3.4.0
```

### Build and test

```bash
git submodule update --init --recursive
cmake --preset linux-release            # macos-release, windows; tests are on in the presets
cmake --build --preset linux-release
ctest --preset linux-release -j3        # full suite ~10 min; core/ is the working directory
ctest --preset linux-release -E "^(Arora2021|Casas2012|Ferre2026|Nguyen2022|Nguyen2023|Perez2019|Streb2021|Vysyaraju2024|Ward2022|WardPini2024)Test\.|^CSSAcceleration\."   # fast subset
rm -rf core/Results                    # the suite writes ~4 GB there
pip install . && pytest python/tests    # Python package + example smoke tests
pip install -r requirements-docs.txt && mkdocs build --strict
```

### Conventions

- **Never commit `core/Results/`** (test output, git-ignored) or build trees.
- **Keep generated docs generated.** `docs/validation.md` comes from `core/test/papers/validation_summary.json` (`python scripts/docs/gen_validation_page.py`, `--merge <build>/validation` after running the paper tests); `docs/adsorbent-library.md` from `data/adsorbents/library.json` (`python scripts/docs/gen_adsorbent_library.py`, which also validates the catalogue). Never edit the generated pages; CI runs both with `--check`.
- **Cite every adsorbent parameter.** A library entry needs a DOI and table/figure reference for each value, or an explicit `Not stated in source (checked: ...)`. Paper reference values (`reference_values.csv`) must be numbers printed in the paper; never eyeball values off a plot, never commit PDFs or full-text extractions.
- **Tolerances in paper tests are regression bounds** set from the measured error; never loosen one to hide a disagreement — document it in the paper's README.
- **Regenerated example decks are expected.** The test suite re-saves `examples/decks/*.json` and some `data/adsorbents/*.json`; embedded result numbers move in the last digits. Commit those changes with the code change that caused them.
- **Data file names must be shell-safe** (`[A-Za-z0-9._,=+@-]`; checked by the `DataFileNamesAreShellSafe` test).
- **Bump `cereal_class_version`** for any change to a serialised class and handle older versions on load.
- **The engine stays front-end free:** no GUI, web or JSON-schema code in `core/`; it must keep building as a subproject and with Emscripten (`scripts/test_subproject_consumer.sh [--emscripten]`).
- Licence headers point to the root `LICENSE`; do not add per-file licence text. Leave third-party code in `external/` untouched.
- Commits need a DCO sign-off (`git commit -s`); see CONTRIBUTING.md.
