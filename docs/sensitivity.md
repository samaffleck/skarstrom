# Sensitivity Analysis and Error Bars

A sensitivity study answers two questions about a deck:

1. **Which inputs matter?** Of the uncertain inputs (step times, flows, bed voidage, isotherm constants...), which ones move purity, recovery, productivity and energy, and which can be left alone?
2. **How sure are the KPIs?** If each input is only known within some range, what range should each KPI be quoted with? These are the error bars.

Skarstrom answers them in two stages, the approach recommended for slow models (Saltelli et al., *Global Sensitivity Analysis: The Primer*, Wiley 2008, [doi:10.1002/9780470725184](https://doi.org/10.1002/9780470725184)):

| Stage | Method | Runs | Answers |
|---|---|---|---|
| Screening | Morris elementary effects | trajectories × (inputs + 1), e.g. 10 × 11 = 110 | which inputs matter, and whether they act alone |
| Propagation | Monte Carlo with Latin hypercube sampling | 100–400, independent of the number of inputs | the KPIs' spread: median, 90 % band, each input's share |

Both vary several inputs together. Changing one input at a time around the nominal design (a classic parametric study) only looks along a few straight lines through the design space and misses inputs whose effect depends on the others.

## Describing the inputs

Each input is a parameter path (see [Python and the Command Line](python.md#parameter-paths)) with a distribution, in the parameter's SI unit:

| Distribution | Python | CLI `--input PATH=` | Use for |
|---|---|---|---|
| uniform ± fraction | `sa.relative(path, 0.1)` | `10%` | an operating value known to ±10 % |
| uniform | `sa.uniform(path, lo, hi)` | `LO:HI` | a value somewhere in a range |
| triangular | `sa.triangular(path, lo, mode, hi)` | `triangular:LO:MODE:HI` | a range with a most likely value |
| normal | `sa.normal(path, mean, sd, lower=None, upper=None)` | `normal:MEAN:SD[:LO:HI]` | a measured or fitted value and its standard error |

The parameter must be tunable (`skarstrom params deck.json` lists those). A value written as a formula of user variables is not: vary the variable instead. Inputs are independent. Isotherm constants fitted together are correlated, so giving each its own standard error overstates the spread a little; a joint treatment is planned.

The ranges decide the answer. Record where each came from, and remember that the error bars cover uncertain inputs only, not model error (see [Validation](validation.md)).

## Running a study in Python

```python
import skarstrom
from skarstrom import sensitivity as sa

sim = skarstrom.examples.Perez2019_A.load()
study = sa.Study(sim, [
    sa.relative("step[1].time", 0.10),
    sa.relative("endpoint[1].velocity", 0.05),
    sa.uniform("layer[1].bed_void_fraction", 0.36, 0.40),
    sa.normal("adsorbent[1].isotherm[CarbonDioxide].b1", 1.0e-5, 0.1e-5, lower=0),
], settings={"css.molar_relative_tolerance": 1e-4}, workers=4)

screen = study.screen(trajectories=10)       # 10 x (4 + 1) = 50 runs
print(screen.table())                        # mu*, its 95 % interval, mu, sigma
keep = screen.important(threshold=0.1)       # inputs with mu* >= 10 % of the largest

spread = study.propagate(inputs=keep)        # 100 runs, then +50 until converged (max 400)
print(spread.table(nominal=True))            # nominal, mean, sd, p05, p50, p95, converged, valid
print(spread.shares())                       # each input's share of each output's spread
spread.save("perez_spread.json")             # study, every run and the results
```

`outputs` (default: purity, recovery, mass productivity and specific energy of every tracked component) are `"Component.kpi"` names; `sa.kpi_names()` lists the KPIs. `workers` runs that many processes. `seed` (default 1) makes a study repeatable on any platform. Ctrl-C stops a study.

## Running a study from the command line

```bash
skarstrom sensitivity Perez2019_A \
    --input "step[1].time=10%" --input "layer[1].bed_void_fraction=0.36:0.40" \
    --input "adsorbent[1].isotherm[CarbonDioxide].b1=normal:1e-5:1e-6:0:" \
    --output CarbonDioxide.purity --output CarbonDioxide.recovery \
    --setting css.molar_relative_tolerance=1e-4 --workers 4 --save study.json
```

`--method both` (the default) screens all inputs, then propagates those whose μ\* reaches `--threshold` (0.1) of the largest; `screen` or `propagate` runs one stage. `--json` prints everything as JSON.

## Reading the results

**Screening.** For each input and output:

- **μ\*** (mu-star) ranks the inputs. It is in the output's units: the typical change of the KPI when the input moves across its whole range. An input with a μ\* near zero can be fixed at its nominal value.
- **σ** (sigma) is how much that effect changes from one point of the ranges to another. σ small next to μ\*: the input acts alone and roughly linearly. σ comparable to or larger than μ\*: its effect depends on the other inputs, or is curved.
- **μ** is the signed mean: the direction of the effect, when it is consistent.

**Propagation.** For each output:

- **p05, p50, p95**: the 5th, 50th and 95th percentiles of the KPI over the runs. Quote the KPI as the median with the 90 % band p05–p95.
- **converged**: each percentile is known to within `rel_tol` (10 %) of the band's width (from a bootstrap). Runs are added in batches until every output converges or `max_runs` is reached.
- **shares**: each input's share of the output's spread, from a regression on ranks; they add up to `r2`. An `r2` well below 1 means the inputs act through interactions, and the shares are only part of the story.

## Rules that keep the numbers honest

- **Every run must reach cyclic steady state.** A run that stops at the cycle limit is not settled, and the difference would read as sensitivity. Such runs are counted as not accepted (`require_css=True`, the default). Tighten the CSS tolerance for a study with `settings`, so that solver noise stays well below the effects being measured.
- **Failed runs are reported, never dropped silently.** `failed` and `warnings` say how many runs did not count and why. Failures tend to cluster in one corner of the ranges, so the results then describe the rest of it.
- **Fixed-cycle screening** (`require_css=False`, `cycles=N`) is fine for a quick ranking when every run gets the same number of cycles. Do not use it for error bars.

## The C++ API

The engine plans the runs and analyses the results; the caller executes them, so runs can go to any number of processes or web workers. Everything is in `adsim/solver/sensitivity.h`:

```cpp
#include "adsim/solver/sensitivity.h"
namespace sa = adsim::sensitivity;

std::vector<sa::Input> inputs = {sa::Relative("step[1].time", 99.0, 0.1),
                                 sa::Uniform("layer[1].bed_void_fraction", 0.36, 0.40)};
std::vector<sa::Output> outputs = {sa::ParseOutput("CarbonDioxide.purity")};
sa::CheckStudy(deck_json, inputs, outputs);

sa::Plan plan = sa::MorrisPlan(inputs, 10);        // or sa::MonteCarloPlan(inputs, 100, seed, batch)
std::vector<std::vector<double>> y;
for (const auto& row : plan.values)
  y.push_back(sa::EvaluateSample(deck_json, inputs, row, outputs).outputs);
sa::MorrisResult m = sa::AnalyseMorris(plan, {"CarbonDioxide.purity"}, y);
```

`AnalyseMonteCarlo(inputs, values, names, y)` takes the runs of every batch stacked. Plans use `std::mt19937_64` with every draw done by the engine, so a seed gives the same runs natively, in the Python wheels and in WebAssembly.
