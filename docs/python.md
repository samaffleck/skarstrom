# Python and the command line

The `skarstrom` Python package wraps the C++ engine (`adsim`) and runs deck files. Use it for parametric studies, sensitivity analysis, isotherm fitting, your own cost functions and your own plots, or run decks from a terminal with the `skarstrom` command.

## Install

See [Getting Started](getting-started.md) for the details. In short, once the package is on PyPI:

```bash
pip install skarstrom
pip install "skarstrom[examples]"   # also pandas, matplotlib, scipy and Jupyter
```

Until then, or for the development version, build it from a checkout (Python 3.9+, CMake 3.15+, a C++20 compiler):

```bash
git clone --recursive https://github.com/samaffleck/skarstrom.git
cd skarstrom
pip install .                       # builds the engine; several minutes the first time
```

Check the install:

```bash
skarstrom --version
skarstrom examples                  # the bundled decks
skarstrom run Perez2019_A --cycles 3
```

The package bundles every example deck of the repository and the curated adsorbent library (see [Bundled data](#bundled-data)). Any deck file written by `sim.save()`, by the C++ API or by a front end built on the engine also works.

## A first run

```python
import skarstrom

sim = skarstrom.examples.Perez2019_A.load()     # a bundled deck; or skarstrom.Simulation("my_deck.json")
sim.cycle.number_of_cycles = 20
sim.steps["Feed"].time = "2 min"                # units are explicit; a bare number is SI (seconds)
result = sim.run()

print(result.status, result.message)            # e.g. "reached_css", "cyclic steady state reached at cycle 14"
for k in sim.kpis(cycle=-1):                    # KPIs at the last cycle
    print(k.component, k.purity, k.recovery, k.energy_consumption)
```

`run()` blocks until the run finishes. It returns a `RunResult` with:

| Field | Meaning |
|---|---|
| `status` | `completed`, `reached_css`, `failed` or `cancelled` |
| `ok` | `True` for `completed` and `reached_css` |
| `message` | why the run stopped (for a failure, the solver's reason) |
| `css_reached`, `css_cycle` | whether [cyclic steady state](cyclic-steady-state.md) was reached, and at which cycle (0-based, `-1` if not) |
| `cycles_run`, `wall_time_s` | cycles completed and wall-clock time |
| `stats` | integrator counters: steps, residual evaluations, error-test and convergence failures, Jacobian evaluations |

Options: `run(cycles=N)` overrides the number of cycles. `progress=f` calls `f(info)` after every step with a dict holding `cycle`, `step`, `total_cycles`, `total_steps` and `fraction`. Return `False` from it to cancel. `log_level="medium"` and `on_message=f` pass the solver's log messages to `f(level, text)`. Ctrl-C cancels a run and raises `KeyboardInterrupt`.

## Example decks

`skarstrom.examples` has one attribute per bundled deck, so your editor completes the names and shows each deck's title and description:

```python
ex = skarstrom.examples.Ward2022_Q100ccm
ex.title, ex.description, ex.path
sim = ex.load()                                 # same as skarstrom.Simulation(ex)
ex.copy_to("my_deck.json")                      # an editable copy
skarstrom.examples.names()                      # every name; skarstrom.examples.ALL holds the objects
skarstrom.examples.get("Perez2019_A")           # by name, when the name is in a variable
```

The decks are those of `examples/decks/` in the repository ([Examples](examples.md)).

## The deck as objects

A `Simulation` exposes the parts of its deck as objects. Reading an attribute reads the deck and assigning one sets it:

```python
sim = skarstrom.examples.Perez2019_A.load()
sim.steps                                       # the cycle's steps, in order
sim.steps["Feed"].time = 120                    # a step by name...
sim.steps[0].time                               # ...or by position (0-based)
sim.steps["Feed"].bottom.temperature = "25 degC"  # the feed-end stream of that step
sim.steps["Feed"].bottom.y["CO2"] = 0.15        # inlet mole fractions by component
sim.steps["Blowdown"].top.pressure_exponential_end = "0.1 bar"
sim.layers[0].length = "1.2 m"                  # adsorbent layers, feed end first
sim.layers[0].bed_void_fraction = 0.37
ads = sim.adsorbents["13X"]
ads.heat_of_adsorption["CO2"] = "35 kJ/mol"
ads.isotherm["CO2"].b1 = 3.9e-7                 # isotherm parameters (fields depend on the isotherm kind)
ads.kinetics["CO2"].tortuosity = 3.0           # rate-model parameters (fields depend on the model)
sim.variables                                   # user variables by name: sim.variables["feed_flow"] = "2 L/min"
sim.bed.diameter, sim.wall.thickness
sim.settings.reltol = 1e-6                      # solver, CSS and energy-balance settings
sim.css.min_cycles = 3
sim.energy_balance.energy_balance_mode = 1      # isothermal
sim.info.title = "my variant"
```

The rules:

- **Attributes are the fields every object of that kind has** (`time`, `length`, `temperature`). **Brackets select things that are named in the deck**: steps, streams, layers and adsorbents by name, and components. An integer in brackets is a position; `sim.steps.by_id(3)` selects by the ID stored in the deck.
- **Each step has a `bottom` (feed end, z = 0) and a `top`.** When that end has one stream or pressure boundary, its fields are available directly (`step.bottom.temperature`). When it has several, choose one: `step.top["Purge stream"]` or `step.top[0]`.
- **Typos fail loudly.** A misspelt field, or one this deck does not have (`fraction` when the steps are timed, `velocity` on a pressure boundary), raises `AttributeError` with the fields it does have. Assigning never creates a new attribute.
- **Components** are CoolProp names (`"CarbonDioxide"`); the usual formulas (`"CO2"`, `"N2"`, `"CH4"`, `"He"`...) work too.
- **Printing an object** lists its values with units and string paths:

```text
>>> sim.steps["Feed"].bottom
<steps['Feed'].bottom['Feed inlet'] (endpoint id 1)>
  temperature         298 K        endpoint[1].temperature
  velocity            0.08487 m/s  endpoint[1].velocity
  y['Nitrogen']       0.85         endpoint[1].y[Nitrogen]
  y['CarbonDioxide']  0.15         endpoint[1].y[CarbonDioxide]
```

Each object also has `fields` (its field names), `get(name, unit=None)`, `set(name, value, unit=None)`, `text(name)` (the expression behind a value), `path(name)` (the string path, below) and `parameters()` (its entries of `sim.parameters()`, with units and bounds). The objects are typed, so editors complete the fields and show each one's unit.

## Units

**A bare number is in the field's SI unit**: seconds, kelvin, pascal, metre, m/s, mol/s, kg/m³, J/mol and so on. The unit is in each field's documentation, in the printed objects, and in `Parameter.unit`. To use another unit, say so:

```python
step = sim.steps["Feed"]
step.time = 120                       # 120 s
step.time = "2 min"                   # a quantity: number and unit
step.time = "2 [min]"                 # the same, in the deck's own syntax
step.set("time", 2, unit="min")       # an explicit unit argument
step.get("time", unit="min")          # 2.0: read back in any unit

import pint
ureg = pint.UnitRegistry()
step.time = 2 * ureg.minute           # pint Quantities work too
step.bottom.temperature = ureg.Quantity(25, ureg.degC)
```

The unit is checked against the field's dimension, so a mistake is an error rather than a silently wrong number:

```text
>>> step.time = "2 bar"
ValueError: parameter step[1].time is a time [s]; '2 bar' is a pressure
```

String units are read by the engine, the same parser decks and `skarstrom run --set` use: SI symbols with prefixes (`kPa`, `mmol`, `cm`, `kJ`), compound units (`kJ/mol`, `W/(m^2*K)`, `mmol/g`, `kg/m3`, `1/bar`) and common named units (`bar`, `atm`, `psi`, `degC`, `degF`, `min`, `h`, `L/min`, `SCCM`, `Nm3/h`, `%`, `ppm`). Everything after a `/` is in the denominator, so `W/m2-K` means W/(m²·K). A temperature in `degC` is absolute on its own and a temperature difference inside a compound unit (`J/(kg*degC)`). pint Quantities are converted by pint itself, so any unit pint knows works. `skarstrom.units.convert(25, "degC", "K")` converts plain numbers.

To make units mandatory, create the simulation with `require_units=True`. A bare number for a dimensional field then raises `TypeError` (dimensionless values, counts and switches still take plain numbers):

```python
sim = skarstrom.Simulation(skarstrom.examples.Perez2019_A, require_units=True)
sim.steps["Feed"].time = 120          # TypeError: step[1].time is a time [s]; ... give the unit
sim.steps["Feed"].time = "120 s"      # fine
```

## Parameter paths

Every value also has a stable string **path**, which is what the object layer uses underneath. Paths are the right tool when the parameter is data: sweeps, sensitivity studies, the optimiser, configuration files and the command line. `sim.parameters()` lists them, and `skarstrom params deck.json` prints the same list:

```python
for p in sim.parameters(tunable_only=True):
    print(p.path, p.value, p.unit, p.bounds, p.alias)
```

Paths look like this:

| Path | Value |
|---|---|
| `cycle.number_of_cycles`, `cycle.total_time` | number of cycles; total cycle time (fraction mode) |
| `step[2].time`, `step[2].fraction` | a step's duration or cycle fraction |
| `endpoint[5].temperature`, `endpoint[5].y[CarbonDioxide]` | inlet temperature and mole fraction of a stream or pressure boundary |
| `endpoint[5].velocity`, `.molar_flow`, `.volumetric_flow`, `.mass_flow` | a stream's constant flow |
| `endpoint[3].pressure` | a pressure boundary's constant pressure (`pressure_linear_end`, `pressure_exponential_end` for ramps) |
| `endpoint[5].pump.efficiency` | pump, compressor and vacuum-pump settings |
| `bed.diameter`, `bed.ambient_temperature`, `wall.thickness` | column and wall |
| `layer[1].length`, `layer[1].bed_void_fraction`, `layer[1].particle_diameter` | layer properties and correlation coefficients |
| `adsorbent[1].solid_density`, `adsorbent[1].heat_of_adsorption[CarbonDioxide]` | adsorbent properties |
| `adsorbent[1].isotherm[CarbonDioxide].b1`, `adsorbent[1].kinetics[CarbonDioxide].k` | isotherm and rate-model parameters |
| `variable[feed_flow]` | a [user variable](cycles.md) |
| `layer[1].momentum_model`, `energy_balance.energy_balance_mode`, `energy_balance.isothermal_temperature` | model switches (integers or true/false) and the isothermal temperature |
| `settings.reltol`, `css.min_cycles`, `info.title` | solver settings, CSS settings, deck information |

The numbers in brackets are IDs stored in the deck, so a path keeps working when you rename a step, stream, layer or adsorbent. Each parameter also has an **alias** that uses the display names instead, for example `step[Feed].time`, `step[Feed].bottom[Feed gas].temperature` or `layer[Layer 1].length`. You can use either form. An alias that matches two objects with the same name is rejected, and you then need the ID path. `obj.path("time")` gives the path of any object's field.

```python
sim["step[Feed].time"] = 120.0                # a number (SI)
sim["endpoint[2].pressure"] = "1.5 bar"       # a quantity
sim["step[2].time"] = "0.2 * 95"              # an expression (it may use the deck's user variables)
sim.set("step[1].time", 2, unit="min")
sim.update({"css.min_cycles": 3, "settings.reltol": 1e-6})
print(sim["step[Feed].time"], sim.get("step[1].time", unit="min"), sim.get_text("endpoint[2].pressure"))
```

Values follow the same unit rules as the objects. A string sets the parameter's expression, which is evaluated straight away; an invalid expression raises `ValueError` and leaves the old value in place. `tunable_only=True` keeps the numeric model values that the optimiser can vary, each with suggested search bounds. A value written as a formula is not tunable, because setting a number would replace the formula. `sim.save("new_deck.json")` writes the deck with your changes.

### Structural edits

Parameters change values; a few methods change the deck's structure, with the cascades that keep it consistent (they call the engine's `adsim/objects/deck_edit.h`):

```python
sim.add_component("Argon")            # 0 mole fraction everywhere, default physics on every adsorbent
sim.remove_component("Argon")
ads = sim.add_library_adsorbent("zeolite_13x_perez2019")   # returns the new adsorbent's ID
```

A library adsorbent is added to the deck's adsorbent directory; it is used (and its values appear in `sim.parameters()`) once a layer contains it — automatically for a layer that had no adsorbent yet.

Other structural changes (steps, boundary types, layers, connections) are made in C++ or by a front end, and saved in the deck.

## Results

After a run:

```python
sim.cycles                 # logged cycles, e.g. [0, 1, ..., 13]
sim.kpis()                 # KPIRecord per tracked component per cycle
sim.kpi_table()            # the same as a pandas DataFrame
sim.css_history()          # the deck's CSS check at every cycle: cycle, converged, error, tolerance
y = sim.series("mole_fraction", "CarbonDioxide", location="top")   # NumPy arrays y.time, y.values
T = sim.series("temperature", z=0.25, cycle=0)                      # at z = 0.25 m in the first cycle
prof = sim.profile("adsorbed_concentration", "CarbonDioxide", cycle=-1, step=0)
prof.z, prof.time, prof.values                                       # values has shape (len(time), len(z))
for s in sim.streams():                                              # every endpoint: step, name, endpoint_id, location
    f = sim.stream_series(s["step"], s["name"], "molar_flow_rate")
    f = sim.endpoint_series(s["endpoint_id"], "molar_flow_rate")    # the same, by ID
```

- **KPIs** are computed for the components ticked in the deck's [KPI configuration](performance-metrics.md): purity, recovery, productivity (mass and volume), specific energy consumption and production rate. A KPI that is undefined (no product drawn) is `nan`.
- **Series** are one cycle long and start at `t = 0` at the start of that cycle. `location` is `"top"`, `"bottom"` or `"position"` (implied by `z=`). `cycle=-1` is the last logged cycle.
- **Variable names**: `skarstrom.bulk_variables()` (pressure, temperature, velocity, ...) and `skarstrom.component_variables()` (mole_fraction, adsorbed_concentration, ...) list what `series` and `profile` accept. `stream_series` accepts `temperature`, `molar_flow_rate` and `compressor_power`, or `mole_fraction` and `molar_flow` with a component.
- `sim.save_logs("out")` writes the full solver log as CSV files.

The solver keeps the last `settings.max_cycles_to_keep` cycles (10 by default) in memory.

## Isotherms

The isotherm functions use the engine's own isotherm code and the fitter behind the [Isotherm Fitting](isotherm-fitting.md) tool:

```python
skarstrom.isotherm_kinds()                        # e.g. "single_site_langmuir_p", "sips", "toth", ...
skarstrom.isotherm_parameters("single_site_langmuir_p")   # ["m1", "b1", "q1"]
q = skarstrom.evaluate_isotherm("single_site_langmuir_p", {"m1": 4.0, "b1": 2e-9, "q1": 2.5e4},
                                pressure=p_pa, temperature=303.15)
fit = skarstrom.fit_isotherm("sips", p_pa, q_molkg, temperature=T_K)
fit["params"], fit["std_errors"], fit["r2"], fit["rmse"]
```

`fit_isotherm` accepts `fixed={...}` to hold parameters, and `initial={...}` and `bounds={...}` to override the starting guesses and search bounds derived from the data.

## Command line

```bash
skarstrom run deck.json [--cycles N] [--set PATH=VALUE ...] [--csv DIR] [--json] [--log low|medium|high] [-q]
skarstrom params deck.json [--all] [--search TEXT] [--json]
skarstrom examples [list] [--json]          # the bundled decks
skarstrom examples run NAME [run options]   # run a bundled deck
skarstrom examples path NAME                # where it is installed
skarstrom examples copy NAME [DEST]         # copy it to edit
```

`skarstrom run` prints the run status, the CSS outcome and the KPIs at the last cycle. The deck may be a file or the name of a bundled example (with or without `.json`). `--set` takes the same paths and values as Python, units included, and can be repeated, for example `--set "step[Feed].time=2 min" --set "endpoint[3].pressure=1.5 bar"`. `--csv DIR` writes `kpis.csv` (every cycle), `css.csv`, `run.json` and the solver's CSV logs. `--json` prints the summary as JSON (status, message, CSS, cycles, wall time, integrator statistics, KPIs per cycle with `null` for undefined values, and the CSS history), which is the easiest way to drive the engine from another program:

```bash
skarstrom examples run "Ward2022_Q100ccm" --json -q | jq '.status, .css_reached'
```

The exit code is 0 when the run completes or reaches CSS, 1 when it fails and 130 when it is interrupted.

`skarstrom params` lists the tunable parameters of a deck with their values, units, suggested bounds and aliases. `--all` adds solver settings, integers and text fields; `--json` prints them as a list of objects.

## Bundled data

```python
skarstrom.examples.Perez2019_A.path        # the example decks (see Example decks above)
lib = skarstrom.adsorbent_library()        # data/adsorbents/library.json as a dict
skarstrom.adsorbent_path("zeolite_13x_perez2019")
```

The decks are those of `examples/decks/` in the repository and the library is `data/adsorbents/` ([Adsorbent Library](adsorbent-library.md)).

## Examples and notebooks

[Examples](examples.md) lists a runnable script for every deck family in `examples/` (each runs its decks, prints KPIs or the CSS status, and compares with the published or digitised data). The repository also has three executed notebooks in `examples/notebooks/`:

- `parametric_study.ipynb` sweeps the feed-step time of a CO₂ VSA cycle and plots purity and recovery.
- `sensitivity_analysis.ipynb` computes local sensitivities of the KPIs to a set of tunable parameters by finite differences.
- `isotherm_fit.ipynb` fits single-temperature and temperature-dependent isotherm models with the engine's fitter.

## Licence

Skarstrom is open source under the Apache License 2.0. The wheels contain the engine and its third-party dependencies; `LICENSE`, `NOTICE` and `THIRD_PARTY_NOTICES.md` are installed in `skarstrom/licenses/`.
