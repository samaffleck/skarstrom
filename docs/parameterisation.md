# Parametric Studies and Parameterisation

<!-- doc-setup
import os
import skarstrom
base = skarstrom.examples.Perez2019_A.load()
os.symlink(os.environ["SKARSTROM_DOC_REPO"] + "/core", "core")   # the page reads the repository's measured data
-->

Every settable deck value is a parameter with a stable path ([Python and the Command Line](python.md#parameter-paths)), so studies are ordinary Python loops around `sim.set(...)` and `sim.run()`. This page shows the three common patterns: a **parametric study** (sweep one input and watch a KPI), an **optimisation** (search inputs for the best KPI), and a **parameterisation** (fit deck parameters to measured data with `skarstrom.calibrate`). If you only need isotherm parameters from equilibrium uptake data, use [Isotherm Fitting](isotherm-fitting.md) instead: it fits the isotherm directly, with no column simulation, and is far faster.

Each trial is a full cycle simulation. Start from a copy of the deck (`sim.copy()`) so every trial begins from the same inputs, and keep the number of trials small.

## Choosing parameters

```python
import skarstrom

base = skarstrom.examples.Perez2019_A.load()
for p in base.parameters(tunable_only=True, search="time"):
    print(p.path, p.alias, p.value, p.unit, p.bounds)
```

`tunable_only=True` keeps the numeric model values that can be varied, each with suggested search bounds (`p.bounds`, by default ±50 % of the current value). Parameters are in SI base units (give another unit as a string, `"2 min"`, or with `unit=`; see [Units](python.md#units)). A parameter written as a formula is not tunable, because setting a number would replace the formula; vary the [global variable](global-variables.md) it refers to instead (`variable[<name>]`).

## Parametric study

Sweep the feed-step time of a CO₂ VSA cycle and record purity and recovery at cyclic steady state:

```python
import numpy as np

rows = []
for t_feed in np.linspace(30, 90, 5):
    sim = base.copy()
    sim.steps["Feed"].time = t_feed            # seconds
    r = sim.run()
    k = next(k for k in sim.kpis(cycle=-1) if k.component == "CarbonDioxide")
    rows.append((t_feed, k.purity, k.recovery, r.css_cycle))
```

`examples/notebooks/parametric_study.ipynb` is the executed version, with plots, and [Sensitivity Analysis and Error Bars](sensitivity.md) covers screening uncertain inputs and putting error bars on the KPIs. Independent trials parallelise well: run them in a `concurrent.futures.ProcessPoolExecutor`, one `Simulation` per process.

## Optimisation

Maximise CO₂ recovery subject to a purity floor, varying two inputs, with SciPy (any optimiser that takes a Python objective works):

<!-- long: 30 cycle simulations to optimum -->
```python
from scipy.optimize import minimize

def objective(x):
    sim = base.copy()
    sim.update({"step[Feed].time": x[0], "step[Evacuation].time": x[1]})
    if not sim.run().ok:
        return 1e3                              # a failed trial scores worst
    k = next(k for k in sim.kpis(cycle=-1) if k.component == "CarbonDioxide")
    penalty = 100.0 * max(0.0, 0.95 - k.purity)  # purity >= 95 %
    return -k.recovery + penalty

res = minimize(objective, x0=[60.0, 120.0], method="Nelder-Mead", options={"maxfev": 30})
```

### Multi-objective optimisation (NSGA-II)

Purity and recovery pull against each other, so a design study usually wants the whole trade-off (the Pareto front) rather than one optimum. `examples/optimise_nsga2.py` does this with [pymoo](https://pymoo.org)'s NSGA-II (`pip install pymoo`): two decision variables (feed time and evacuation pressure of `Perez2019_A`), two objectives (CO₂ purity and recovery at cyclic steady state), and every generation's candidates run in parallel processes. The core of it:

<!-- long: needs pymoo; runs 20 generations of 24 cycle simulations -->
```python
from concurrent.futures import ProcessPoolExecutor
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.problem import Problem
from pymoo.optimize import minimize

def evaluate(x):                                   # one candidate, in a worker process
    sim = base.copy()
    sim.set("step[Feed].time", x[0], unit="s")
    sim.set("step[Evacuation].bottom.pressure_double_exponential_end", x[1], unit="kPa")
    if not sim.run().ok:
        return [1.0, 1.0]                          # a failed run scores worst
    k = next(k for k in sim.kpis(cycle=-1) if k.component == "CarbonDioxide")
    return [-k.purity, -k.recovery]                # pymoo minimises

class CycleDesign(Problem):
    def __init__(self):
        super().__init__(n_var=2, n_obj=2, xl=[60, 2], xu=[150, 6])
    def _evaluate(self, X, out, *args, **kwargs):
        out["F"] = list(pool.map(evaluate, X))

with ProcessPoolExecutor() as pool:
    res = minimize(CycleDesign(), NSGA2(pop_size=24), ("n_gen", 20))
res.X, -res.F                                      # Pareto-optimal designs and their KPIs
```

Each run is single-threaded, so a population evaluates on as many cores as you give it. Constraints (a purity floor, a power limit) go in pymoo's `G` output; other optimisers (SciPy, Optuna, BoTorch) wrap the same `evaluate` function. Use the per-component CSS check for the KPIs a search compares ([Cyclic steady state](cyclic-steady-state.md)).

The C++ engine also has a derivative-free global optimiser (`adsim/solver/optimisation.h`: dlib's MaxLIPO with a trust region), which front ends use; from Python, use any Python optimiser as above.

## Fitting parameters to measured data

!!! note "New in 0.12"
    `skarstrom.calibrate`, `Simulation.calibrate`, `skarstrom.calibration` (`Case`, `Parameter`,
    `CalibrationResult`, `CancelToken`), `skarstrom calibrate` on the command line, and the C++
    object `adsim/solver/calibration.h`.

To fit deck parameters to measurements (a saturation capacity to a breakthrough curve, a
dispersion coefficient per flow rate, a heat-transfer coefficient to a thermocouple), call
`sim.calibrate` for one deck or `skarstrom.calibrate` for several. Do not write your own optimiser
loop: the calibration object handles bounds, log scales, several datasets and several experiments,
runs the simulations in parallel threads, and gives standard errors and correlations.

It fits the measured data stored in the deck ([Experimental data](python.md#experimental-data)).
Each dataset must be bound to the quantity it measures (`dataset.observable`, see
[What a dataset measures](python.md#what-a-dataset-measures)); the bundled decks with data are
bound already.

### One deck

```python
import skarstrom
from skarstrom import calibration as cal

sim = skarstrom.examples.Illustrative_Diluted_Bed.load()
sim["layer[1].cells"] = 12                     # a coarse grid keeps this example quick
sim.experimental_data                          # outlet CO2, outlet flow and a thermocouple, all bound

result = sim.calibrate(
    [cal.Parameter("adsorbent[1].isotherm[CarbonDioxide].m2", 4, 8)],   # path, lower, upper (SI)
    datasets=["outlet_CO2_mole_fraction[-]_Qin=100[ccm]"],               # default: every bound dataset
    max_simulations=20)
print(result.table())
```

<!-- no-run: printed output -->
```console
parameter                                estimate  std error  rel.     unit    start  bounds  scale
---------------------------------------  --------  ---------  -------  ------  -----  ------  ------
adsorbent[1].isotherm[CarbonDioxide].m2  6.86277   0.0147849  0.215 %  mol/kg  6.04   4 .. 8  linear

dataset                                   case    points  rmse      nrmse      max |error|  weight
----------------------------------------  ------  ------  --------  ---------  -----------  ------
outlet_CO2_mole_fraction[-]_Qin=100[ccm]  case 1  29      0.059729  0.0599598  0.139976     1

objective S = 0.00359518, sigma^2 = 1.2840e-04, 29 residuals
least_squares: budget after 11 iterations, 20 simulations, 0.9 s
reached the budget of 20 runs
```

The numbers depend on the grid: a 12-cell bed is too coarse for a real fit. Read the result:

```python
est = result["adsorbent[1].isotherm[CarbonDioxide].m2"]
est.value, est.std_error, est.at_bound          # SI units (mol/kg here)
fit = result.datasets[0]
fit.rmse, fit.nrmse, fit.points                 # in the dataset's own units
fit.x, fit.measured, fit.simulated, fit.residuals   # NumPy arrays at the measured points
result.objective_trace, result.simulations, result.stop_reason, result.notes
result.apply(sim)                               # write the estimate into the deck (a normal edit)
calibrated = result.calibrated()                # or a new Simulation with it
result.to_json("fit.json")                      # everything, including the calibrated deck
```

`result.plot()` (matplotlib) shows each dataset with the simulation at the start and at the
estimate, and the residuals below. The deck passed to `calibrate` is never changed by the fit.

A parameter is a registry path (`sim.parameters(tunable_only=True)` lists them). A bare path,
`sim.calibrate(["adsorbent[1].isotherm[CarbonDioxide].m2"])`, takes the registry's suggested bounds;
`cal.Parameter(path, lower, upper, scale="linear"|"log", shared=True, start=None)` sets them.
Bounds and the start are in SI units, or carry their own like everywhere else in the package:
`cal.Parameter("step[1].time", "10 min", "1 h")`, `cal.Parameter("step[1].time", 10, 60, unit="min")`
or a pint quantity. The start defaults to the deck's value.

### The objective

The objective is a weighted sum of squares over every measured point $k$ of every dataset $d$:

$$
S = \sum_k r_k^2, \qquad r_k = \sqrt{\frac{w_d}{n_d}}\;\frac{y_{sim}(x_k) - y_k}{s_k}
$$

with $y_{sim}$ interpolated linearly at each measured $x_k$ (a KPI is matched by cycle). $s_k$ is
the point's uncertainty when the dataset gives one, else the range of the dataset's $y$ values, so
datasets in different units are comparable. $w_d$ is the dataset's weight and $n_d$ its number of
points, so a densely sampled curve does not outweigh a sparse one. Points that the simulation does
not cover (measured after the end of the simulated cycle) are left out and counted
(`fit.skipped`, with a note).

For a cyclic process measured once per cycle,
[`examples/calibrate_cycle_kpis.py`](https://github.com/samaffleck/skarstrom/blob/main/examples/calibrate_cycle_kpis.py)
binds purity and recovery histories to KPI observables (`sk.Observable("recovery", "CarbonDioxide",
cycle_offset=-1)`, x the cycle number) and fits two parameters of `Perez2019_A` to them. *New in 0.12.*

### Several experiments: shared and per-case parameters

Experiments at different conditions (flow rates, feed compositions) are fitted together as
*cases*. A **shared** parameter takes one value for every case (an isotherm capacity); a
**per-case** parameter takes one value in each (a dispersion coefficient, which depends on the
flow). The three Ward 2022 breakthrough runs, with the outlet CO₂ and a thermocouple at half weight:

```python
cases = []
for q in (50, 100, 200):
    deck = skarstrom.examples.get(f"Ward2022_Q{q}ccm").load()
    deck["layer[1].cells"] = 20                                    # coarse, for this example
    cases.append(cal.Case(deck, [f"outlet_CO2_mole_fraction[-]_Qin={q}[ccm]",
                                 (f"temperature[K]_z=0.0715[m]_Qin={q}[ccm]", 0.5)],   # (name, weight)
                          name=f"{q} ccm"))

result = skarstrom.calibrate(
    cases,
    [cal.Parameter("adsorbent[1].isotherm[CarbonDioxide].m2", 3, 8),             # shared
     cal.Parameter("layer[1].dispersion_coefficient", 1e-6, 1e-3, scale="log", shared=False)],
    differences="central", max_simulations=60)
for p in result.parameters:
    print(f"{p.name}  {p.value:.4g} +- {p.std_error:.2g}")
result.correlation                              # NumPy matrix, in the order of result.parameters
```

<!-- no-run: printed output -->
```console
adsorbent[1].isotherm[CarbonDioxide].m2  5.9 +- 0.0036
layer[1].dispersion_coefficient [50 ccm]  3.411e-05 +- 2.6e-06
layer[1].dispersion_coefficient [100 ccm]  3.656e-05 +- 5.5e-06
layer[1].dispersion_coefficient [200 ccm]  9.82e-05 +- 1.2e-05
```

This took 60 simulations (about 3 s on four threads). A case is a `cal.Case(sim, datasets, weight=,
name=, cycles=, require_css=, settings=)`, a `(sim, [dataset names], weight)` tuple or just a
Simulation (every bound dataset). A per-case parameter can have bounds and a start per case
(`case_bounds=`, `case_start=`). `result.apply(i)` writes case `i`'s values into its own deck, and
`result.calibrated()` returns one calibrated Simulation per case. The deck is not changed by the
fit, and the shipped Ward decks are already calibrated against these curves (see
[Validation](validation.md#ward-2022)): the example is a demonstration, not a recommendation.
The full fit that calibrated those decks (shared q_s,2 and dispersion factor, one inside heat
transfer coefficient per flow, the digitising errors as uncertainties) is
[`examples/calibrate_ward_2022.py`](https://github.com/samaffleck/skarstrom/blob/main/examples/calibrate_ward_2022.py),
which reproduces the published values and standard errors (with a notebook,
`examples/notebooks/calibrate_ward_2022.ipynb`). *New in 0.12.*

### Reading the standard errors and correlations

The standard errors come from the Jacobian $J$ at the estimate, $\mathrm{cov} = \sigma^2 (J^\mathsf{T}
J)^{-1}$ with $\sigma^2 = S / (N - p)$ for $N$ residuals and $p$ estimated values (`result.sigma2`).
When the datasets give uncertainties, $\sigma^2$ near 1 means that the residuals match them. For a
log-scale parameter the standard error is the value times the standard error of its logarithm.

- A **large relative error** (`est.relative_error`) means the data barely fix the value. Fit fewer
  parameters, or add a dataset that responds to it.
- A **correlation** near ±1 (`result.correlation`) means two parameters change the fit in the same
  way, and the data fix only a combination of them (an isotherm capacity and affinity, typically).
  Their individual standard errors are then large; fix one of them, or add data at other conditions.
- A standard error of **NaN** means the data do not depend on a parameter at all, or $N \le p$:
  `result.notes` says which.
- At a bound (`est.at_bound`), the standard error is only a local statement: widen the bounds if
  the bound is not physical.

A fitted value absorbs whatever else moves the measured curve (dead volume, the detector, the bed
temperature), so check it against independent data: fit one flow rate and compare with another.

### Choosing the method

- `method="least_squares"` (the default) is a bounded Levenberg–Marquardt method with a
  finite-difference Jacobian. Each Jacobian column is one run per case that the parameter affects,
  and the columns run in parallel threads (`threads=`, default the cores minus one). It is the fast
  choice from a reasonable start, and the only one that gives standard errors along the way.
- `method="global"` is a derivative-free MaxLIPO search over the bounds (`global_evaluations=`
  candidates, default 30; the cases of a candidate run in parallel). Choose it when the start is
  poor or the objective may have several minima. `refine=True` then polishes its best point with
  least squares, which also gives the standard errors.
- `differences="central"` doubles the cost of each Jacobian but is needed for strongly correlated
  parameters: forward differences can stall the fit in the narrow valley of the objective.
- `max_simulations` caps the runs (one case run is one simulation); `ftol` and `xtol` stop least
  squares when the objective or the values stop changing.

```python
result = sim.calibrate([("adsorbent[1].isotherm[CarbonDioxide].m2", 4, 8)], method="global",
                       global_evaluations=4, threads=2)      # (path, lower, upper) also works
result.stop_reason                              # "global_done"
```

`progress=f` calls `f(info)` after each iteration (`info["iteration"]`, `["objective"]`,
`["values"]`, `["simulations"]`); return `False` to stop. Ctrl-C stops the fit and raises
`KeyboardInterrupt`; a `cal.CancelToken()` passed as `cancel=` stops it from another thread, and
the result is then the best point so far (`result.cancelled`). For a cyclic deck, runs start from
the cyclic steady state of the nearest point already run (`warm_start=True`); tighten the CSS
tolerance in the case's `settings` so that it does not add noise to the finite differences.

From a terminal, `skarstrom calibrate` does the same
([Command line](python.md#command-line)):

```bash
skarstrom examples copy Illustrative_Diluted_Bed .
skarstrom calibrate Illustrative_Diluted_Bed.json --set "layer[1].cells=12" \
    --param "adsorbent[1].isotherm[CarbonDioxide].m2:4:8" \
    --dataset "outlet_CO2_mole_fraction[-]_Qin=100[ccm]" \
    --max-simulations 20 --out calibrated.json --json fit.json -q
```

### Calibration from C++

The engine's object is `adsim/solver/calibration.h` (namespace `adsim::calibration`), which the
Python API wraps:

| Object | What it holds |
|---|---|
| `Case` | One experiment: the deck (`Case::FromFile(path)`, `Case::FromData(user_data)`, or the JSON text in `deck_json`), the datasets to fit (`Dataset("name", weight)` or `Dataset(index)`; empty means every bound dataset) and a run policy (`sensitivity::RunPolicy`: a cycle cap, whether a run must reach cyclic steady state, registry text sets applied before the parameters). CSS is not required by default. |
| `Parameter` | A tunable registry path, its bounds, `Scale::Linear` or `Scale::Log`, an optional start (default: the deck's value) and the sharing: `Shared(path, lower, upper)` gives one value for every case, `PerCase(path, lower, upper)` one value per case (with optional per-case bounds and starts). |
| `Problem` | The cases and parameters, checked at construction against every deck: an unknown path, or a dataset that is missing or unbound, is an error that names the case. |
| `Options` | The solver (`Solver::LeastSquares`, the default, or `Solver::Global` with `refine`), the run budget (`max_simulations`), stopping tolerances, the finite-difference scheme and step, the number of threads, warm starts, and the `progress` and `cancel` callbacks. |
| `Result` | Per parameter: the estimate, its standard error, whether it is at a bound. The correlation matrix. Per dataset: RMSE, NRMSE, the maximum error, the points used, the residuals and the simulated series at the estimate. The objective trace, the runs made and failed, the solver and why it stopped, and one calibrated deck per case. |

```cpp
#include "adsim/solver/calibration.h"
namespace cal = adsim::calibration;

std::vector<cal::Case> cases = {
    cal::Case::FromFile("Q50.json", {cal::Dataset("outlet_CO2"), cal::Dataset("T_z=0.0715", 0.5)}),
    cal::Case::FromFile("Q100.json", {cal::Dataset("outlet_CO2"), cal::Dataset("T_z=0.0715", 0.5)})};
std::vector<cal::Parameter> parameters = {
    cal::Shared("adsorbent[1].isotherm[CarbonDioxide].m2", 3.0, 8.0),
    cal::PerCase("layer[1].dispersion_coefficient", 1e-6, 1e-3, cal::Scale::Log)};
cal::Problem problem(cases, parameters);

cal::Options options;
options.max_simulations = 200;
options.progress = [](const cal::Progress& p) { std::printf("%d  S = %g\n", p.iteration, p.objective); };
cal::Result r = cal::Calibrate(problem, options);
for (const auto& p : r.parameters)
  std::printf("%s = %g +- %g%s\n", p.label.c_str(), p.estimate, p.standard_error, p.at_lower || p.at_upper ? " (at a bound)" : "");
```

The deck and dataset names in the C++ example are illustrative. The header documents every option,
the treatment of points the simulation does not cover, warm starts and determinism. The engine also
has the lower-level extraction and RMSE helpers (`adsim/solver/parameterisation.h`:
`SimulatedSeries`, `CompareWithSimulation`, `ExtractBulkTimeSeries`, `ComputeRMSE`), which front
ends use.

### Rolling your own

For an objective the calibration object does not cover, wrap a run in any optimiser. Read the
simulated counterpart of a bound dataset with `dataset.simulated(sim)` or compare with
`dataset.compare(sim)`:

```python
from scipy.optimize import minimize_scalar

base = skarstrom.examples.Ward2022_Q100ccm.load()
base["layer[1].cells"] = 20
path = "adsorbent[1].isotherm[CarbonDioxide].m2"

def rmse(m2):
    trial = base.copy()
    trial[path] = m2
    if not trial.run().ok:
        return 1.0
    return trial.experimental_data[0].compare(trial)["rmse"]

res = minimize_scalar(rmse, bounds=(5.0, 7.0), method="bounded", options={"maxiter": 8})
```

Search positive scale parameters in log space, check that the objective responds to a parameter
before fitting it, and keep the number of runs small: each one is a full simulation.
