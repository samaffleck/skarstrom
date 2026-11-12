# Parametric Studies and Parameterisation

<!-- doc-setup
import os
import skarstrom
base = skarstrom.examples.Perez2019_A.load()
os.symlink(os.environ["SKARSTROM_DOC_REPO"] + "/core", "core")   # the page reads the repository's measured data
-->

Every settable deck value is a parameter with a stable path ([Python and the Command Line](python.md#parameter-paths)), so studies are ordinary Python loops around `sim.set(...)` and `sim.run()`. This page shows the three common patterns: a **parametric study** (sweep one input and watch a KPI), an **optimisation** (search inputs for the best KPI), and a **parameterisation** (fit deck parameters to measured data). If you only need isotherm parameters from equilibrium uptake data, use [Isotherm Fitting](isotherm-fitting.md) instead: it fits the isotherm directly, with no column simulation, and is far faster.

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

To fit, for example, a saturation capacity of the isotherm to a measured breakthrough curve, minimise the root-mean-square error between the simulated series and the data:

$$
\mathrm{RMSE} = \sqrt{\frac{1}{N}\sum_{k=1}^{N}\left(y_{sim}(x_k) - y_k\right)^2}
$$

with $y_{sim}$ interpolated linearly at each data time $x_k$:

```python
import numpy as np
from scipy.optimize import minimize_scalar

base = skarstrom.examples.Ward2022_Q100ccm.load()
t_data, y_data = np.loadtxt(
    "core/data/ward_2022/Qin=100/outlet_CO2_mole_fraction[-]_Qin=100[ccm].csv",
    delimiter=",", skiprows=1, unpack=True)

path = "adsorbent[1].isotherm[CarbonDioxide].m2"     # see base.parameters(tunable_only=True)

def rmse(log10_k):
    sim = base.copy()
    sim[path] = 10.0 ** log10_k
    if not sim.run().ok:
        return 1.0
    y = sim.series("mole_fraction", "CarbonDioxide", location="top", cycle=0)
    return float(np.sqrt(np.mean((np.interp(t_data, y.time, y.values) - y_data) ** 2)))

base_value = base[path]
res = minimize_scalar(rmse, bounds=(np.log10(base_value) - 0.3, np.log10(base_value) + 0.3),
                      method="bounded", options={"maxiter": 20})
```

In our run this took 17 simulations (about 6 s) and moved `m2` from 5.89 to 5.90 mol/kg, lowering the RMSE from 0.027 to 0.021. The shipped deck is already calibrated against these curves (`m2` among the calibrated values; see [Validation](validation.md#ward-2022)), so the search only refines it; on an uncalibrated deck it moves further. That is a demonstration, not a recommendation: a fitted capacity also absorbs whatever else moves the measured front (dead volume, the detector, the bed temperature), so check a fitted value against independent data.

Points to keep in mind:

- **What to compare.** `sim.series(variable, component, location="top" | "bottom", z=...)` gives bulk and component variables at the bed ends or at a position $z$ [m] from the bottom of the column, with time from the start of the cycle. At the bed ends these are the boundary values (an outlet composition is the composition leaving the bed; see [Results](python.md#results)); a thermocouple at height $z$ reads the profile interpolated between the cell centres. Stream quantities (outlet molar flow, compressor power) come from `sim.stream_series(...)`, KPIs per cycle from `sim.kpis()`.
- **Which cycle.** A breakthrough deck has one cycle (`cycle=0`); for a cyclic process compare the last cycle (`cycle=-1`) at cyclic steady state, or KPIs against the cycle number.
- **Search in log space** for positive scale parameters that span decades (rate coefficients, affinities, capacities).
- **Check the sensitivity first.** A parameter the measured quantity barely responds to cannot be fitted; vary it by a factor of two and look at the change in the objective before starting a search.
- **Identifiability.** Fit as few parameters as the data can determine; check the result by fitting a second dataset (another flow rate or feed composition) with the same parameters.

The engine also exposes the same extraction and RMSE helpers in C++ (`adsim/solver/parameterisation.h`: `ExtractBulkTimeSeries`, `ExtractComponentTimeSeries`, `ComputeRMSE`, `ComputeKPIRMSE`), which front ends use to offer fitting with the same objective. Pass the run's boundary log and steps (`BedBoundarySources`) to the extraction helpers so their bed-end values are the endpoint logs, as `sim.series` does (`adsim/solver/bed_series.h`).

## Calibration from C++

!!! note "New in 0.11"
    `adsim/solver/calibration.h`. The Python interface follows in a later release.

The engine has a calibration object for problems larger than a one-parameter search: several datasets, several experiments (decks) at different conditions that share some parameters, and standard errors on the estimates. It works on the datasets stored in each deck, each bound to the quantity it measures ([Results](python.md#results), `XYData::observable`).

The objects (namespace `adsim::calibration`):

| Object | What it holds |
|---|---|
| `Case` | One experiment: the deck (`Case::FromFile(path)`, `Case::FromData(user_data)`, or the JSON text in `deck_json`), the datasets to fit (`Dataset("name", weight)` or `Dataset(index)`; empty means every bound dataset) and a run policy (`sensitivity::RunPolicy`: a cycle cap, whether a run must reach cyclic steady state, registry text sets applied before the parameters). CSS is not required by default. |
| `Parameter` | A tunable registry path, its bounds, `Scale::Linear` or `Scale::Log`, an optional start (default: the deck's value) and the sharing: `Shared(path, lower, upper)` gives one value for every case, `PerCase(path, lower, upper)` one value per case (with optional per-case bounds and starts). |
| `Problem` | The cases and parameters, checked at construction against every deck: an unknown path, or a dataset that is missing or unbound, is an error that names the case. |
| `Options` | The solver (`Solver::LeastSquares`, the default, or `Solver::Global` with `refine`), the run budget (`max_simulations`), stopping tolerances, the finite-difference scheme and step, the number of threads, warm starts, and the `progress` and `cancel` callbacks. |
| `Result` | Per parameter: the estimate, its standard error, whether it is at a bound. The correlation matrix. Per dataset: RMSE, NRMSE, the maximum error, the points used, the residuals and the simulated series at the estimate. The objective trace, the runs made and failed, the solver and why it stopped, and one calibrated deck per case. |

The objective is a weighted sum of squares over every measured point $k$ of every dataset $d$:

$$
S = \sum_k r_k^2, \qquad r_k = \sqrt{\frac{w_d}{n_d}}\;\frac{y_{sim}(x_k) - y_k}{s_k}
$$

where $s_k$ is the point's uncertainty when the dataset gives one, else the range of the dataset's $y$ values, $w_d$ is the dataset's weight and $n_d$ its number of points (so a densely sampled curve does not outweigh a sparse one). The default solver is a bounded Levenberg–Marquardt method with a finite-difference Jacobian; each Jacobian column is one run per case the parameter affects, and the columns run in parallel threads. The global solver is the optimiser's MaxLIPO search on $S$; `refine = true` polishes its best point with least squares.

Standard errors come from the Jacobian $J$ at the estimate, $\mathrm{cov} = \sigma^2 (J^\mathsf{T} J)^{-1}$ with $\sigma^2 = S / (N - p)$ for $N$ residuals and $p$ estimated values. For a log-scale parameter the covariance is in $\ln$ of the value, and the standard error in model units is the value times $\mathrm{se}(\ln x)$. When $N \le p$ or $J^\mathsf{T} J$ is singular, the standard errors are NaN and `Result::notes` says why. A run that fails never stops the fit: it is counted in `failed_simulations` and the step that needed it is rejected.

For strongly correlated parameters (two that the data mostly fix in combination, such as an isotherm's capacity and affinity), set `options.difference = Difference::Central`: forward differences are cheaper but their error can stall the fit in the narrow valley of the objective.

A fit of three breakthrough runs with a shared isotherm capacity and a heat-transfer coefficient per run:

```cpp
#include "adsim/solver/calibration.h"
namespace cal = adsim::calibration;

std::vector<cal::Case> cases = {
    cal::Case::FromFile("Q50.json", {cal::Dataset("outlet_CO2"), cal::Dataset("T_z=0.0715", 0.5)}),
    cal::Case::FromFile("Q100.json", {cal::Dataset("outlet_CO2"), cal::Dataset("T_z=0.0715", 0.5)}),
    cal::Case::FromFile("Q200.json", {cal::Dataset("outlet_CO2"), cal::Dataset("T_z=0.0715", 0.5)})};
std::vector<cal::Parameter> parameters = {
    cal::Shared("adsorbent[1].isotherm[CarbonDioxide].m2", 3.0, 8.0),
    cal::PerCase("step[1].medium.heat_transfer_coefficient", 1.0, 500.0, cal::Scale::Log)};
cal::Problem problem(cases, parameters);

cal::Options options;
options.max_simulations = 200;
options.progress = [](const cal::Progress& p) { std::printf("%d  S = %g\n", p.iteration, p.objective); };
cal::Result r = cal::Calibrate(problem, options);
for (const auto& p : r.parameters)
  std::printf("%s = %g +- %g%s\n", p.label.c_str(), p.estimate, p.standard_error, p.at_lower || p.at_upper ? " (at a bound)" : "");
```

The deck and dataset names in the example are illustrative: list a deck's tunable paths with `skarstrom params deck.json`. The cases run on fresh copies of their decks, so the decks you pass are not changed; `Result::calibrated_decks` holds a copy of each with the estimates applied. The header (`core/include/adsim/solver/calibration.h`) documents every option, the treatment of points the simulation does not cover, warm starts and determinism.
