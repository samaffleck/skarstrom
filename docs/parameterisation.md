# Parametric Studies and Parameterisation

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

`examples/notebooks/parametric_study.ipynb` is the executed version, with plots, and `sensitivity_analysis.ipynb` computes local sensitivities of the KPIs to several tunable parameters by finite differences. Independent trials parallelise well: run them in a `concurrent.futures.ProcessPoolExecutor`, one `Simulation` per process.

## Optimisation

Maximise CO₂ recovery subject to a purity floor, varying two inputs, with SciPy (any optimiser that takes a Python objective works):

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

The C++ engine also has a derivative-free global optimiser for this purpose (`adsim/solver/optimisation.h`: dlib's MaxLIPO with a trust region, over the registry's tunable parameters), which front ends use; it is not wrapped in Python.

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
    "adsim/data/ward_2022/Qin=100/outlet_CO2_mole_fraction[-]_Qin=100[ccm].csv",
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

In our run this took 18 simulations (about 15 s) and moved `m2` from 6.04 to 6.58 mol/kg, lowering the RMSE from 0.38 to 0.095. That is a demonstration, not a recommendation: raising the capacity also absorbs whatever else delays the measured front (dead volume, the detector), so check a fitted value against independent data.

Points to keep in mind:

- **What to compare.** `sim.series(variable, component, location="top" | "bottom", z=...)` gives bulk and component variables at the bed ends or at a position $z$ [m], with time from the start of the cycle. Stream quantities (outlet molar flow, compressor power) come from `sim.stream_series(...)`, KPIs per cycle from `sim.kpis()`.
- **Which cycle.** A breakthrough deck has one cycle (`cycle=0`); for a cyclic process compare the last cycle (`cycle=-1`) at cyclic steady state, or KPIs against the cycle number.
- **Search in log space** for positive scale parameters that span decades (rate coefficients, affinities, capacities).
- **Check the sensitivity first.** A parameter the measured quantity barely responds to cannot be fitted; vary it by a factor of two and look at the change in the objective before starting a search.
- **Identifiability.** Fit as few parameters as the data can determine; check the result by fitting a second dataset (another flow rate or feed composition) with the same parameters.

The engine also exposes the same extraction and RMSE helpers in C++ (`adsim/solver/parameterisation.h`: `ExtractBulkTimeSeries`, `ExtractComponentTimeSeries`, `ComputeRMSE`, `ComputeKPIRMSE`), which front ends use to offer fitting with the same objective.
