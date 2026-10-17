<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/assets/logo-dark.png">
    <img src="docs/assets/logo.png" alt="Skarstrom" width="480">
  </picture>
</p>

# Skarstrom

[![CI](https://github.com/samaffleck/skarstrom/actions/workflows/ci.yml/badge.svg)](https://github.com/samaffleck/skarstrom/actions/workflows/ci.yml)
[![Docs](https://github.com/samaffleck/skarstrom/actions/workflows/docs.yml/badge.svg)](https://samaffleck.github.io/skarstrom/)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)

**Skarstrom** is an open-source simulator for **pressure, vacuum and temperature swing adsorption** (PSA, VSA and TSA) processes. It solves a one-dimensional finite-volume model of a packed adsorption column — coupled mass, momentum and energy balances with non-linear multi-component isotherms and mass-transfer kinetics, a wall energy balance and real-gas properties from CoolProp — as a differential-algebraic system with SUNDIALS IDA. It assembles multi-step, multi-bed cycles (valves, pumps, pressure equalisation, light-product pressurisation, reflux), runs them to cyclic steady state with optional acceleration, and reports purity, recovery, productivity and energy.

- **Engine:** a C++20 library (`adsim`) with a CMake package (`find_package(Adsim)`), usable as a subproject and buildable with Emscripten.
- **Python:** the `skarstrom` package and command-line tool (pybind11), with NumPy results, a parameter registry, and the engine's isotherm fitter.
- **Data:** 33 example decks reproducing published studies, a curated and cited adsorbent library, and automated validation against the literature.

## Why Skarstrom

Adsorption processes such as PSA, VSA and TSA are central to gas separation, carbon capture and gas purification, but there has been no good open-source tool for simulating them. Researchers and engineers write their own column models from scratch, or work in closed commercial packages. That costs time, and it makes results hard to compare: two groups modelling the same cycle rarely share an implementation, a set of parameters or a way to check one another's numbers.

Skarstrom aims to be the shared, open, tested foundation. A standard engine, deck format and cited adsorbent library let the community reproduce published work, compare like with like, and build on each other's models instead of starting again.

## Performance

Skarstrom is built for speed. Optimisation, parameter fitting and sensitivity studies run a cycle to cyclic steady state (CSS) hundreds or thousands of times, so the time to CSS decides what is practical. A full four-step CO₂/N₂ vacuum swing cycle with light-product pressurisation reaches CSS in about **10–17 seconds** on a single core: one core evaluates roughly 200–350 cycle designs an hour, and independent runs parallelise across cores.

| Case (bundled example) | Grid | Cycles to CSS | Time to CSS |
|---|---|---|---|
| 4-step VSA with light-product pressurisation, 15 % CO₂ / N₂ on zeolite 13X (`Perez2019_E`) | 30 cells, 5 steps (feed split in two) | 91 (54 with Aitken acceleration) | 17.4 s (10.1 s) |
| Basic 4-step VSA, 15 % CO₂ / N₂ on zeolite 13X (`Perez2019_A`) | 30 cells, 4 steps | 34 (27) | 6.6 s (5.9 s) |
| 4-step VPSA with light-product pressurisation, CO₂ / N₂ on CALF-20 (`Nguyen2023_Exp_VI`) | 100 cells, 5 steps | 52 (49) | 28.9 s (27.5 s) |
| CO₂ breakthrough on activated carbon (`Ward2022_Q100ccm`) | 60 cells | single run | 0.67 s |
| CO₂ / CH₄ breakthrough on 13X with thermal effects (`Streb2021_A1`) | 60 cells | single run | 2.95 s |

*Non-isothermal model (mass, momentum, energy and wall energy balances) with real-gas properties. CSS when the total mass balance over a cycle closes within 0.5 % (at least 5 cycles). Median of 3 runs, one core of a 2.1 GHz Intel Xeon cloud VM (4 cores), Release build; reproduce with `python scripts/benchmark_css.py`.*

The speed comes from modern C++20 and SUNDIALS:

- **SUNDIALS IDA**, a variable-order, variable-step BDF integrator, solves the discretised column as one differential-algebraic system, taking large steps when the bed is quiet and small ones at step switches.
- **A narrow banded system.** Unknowns are ordered cell by cell, so each Newton iteration solves a banded linear system rather than a general sparse one.
- **A coloured Jacobian.** The Jacobian is built by difference quotients over a sparsity pattern the solver learns and audits, which needs 2.3× fewer residual evaluations than SUNDIALS' default banded Jacobian. That makes the light-product-pressurisation cycle above 1.6× faster.
- **Compiled residuals.** The balances, isotherms and transport correlations are evaluated in compiled C++, and user expressions are compiled once per run instead of being re-parsed at every evaluation.
- **Cyclic steady state acceleration.** Aitken (Graves–Morris) extrapolation of the end-of-cycle state cuts the cycle count further.

## Browser interface

[![The Skarstrom web interface: the simulation tree, solver settings and a plotted pressure profile of the Arora (2021) PSA example](docs/assets/skarstrom-web.png)](https://www.skarstrom.co.uk/)

The [Skarstrom website](https://www.skarstrom.co.uk/) is a browser-based user interface that runs this engine. You can create simulations or import them as JSON deck files, visualise the results, and run optimisation and parameterisation studies.

## Install

Once published on PyPI (Linux, macOS 13.3+ and Windows wheels for Python 3.10–3.13):

```bash
pip install skarstrom
```

From source now (Python 3.9+, CMake 3.15+, a C++20 compiler; Ninja recommended):

```bash
git clone --recursive https://github.com/samaffleck/skarstrom.git
cd skarstrom
pip install .              # builds the engine (several minutes the first time)
```

## Quick start

```python
import skarstrom

sim = skarstrom.examples.Perez2019_A.load()   # 33 bundled decks; or skarstrom.Simulation("deck.json")
sim.steps["Feed"].time = "2 min"          # units are explicit; a bare number is SI
sim.steps["Feed"].bottom.y["CO2"] = 0.15  # the feed stream's CO2 mole fraction
result = sim.run()                        # to cyclic steady state
print(result.status, result.css_cycle)

for k in sim.kpis(cycle=-1):
    print(k.component, f"purity {k.purity:.3f}", f"recovery {k.recovery:.3f}")

y = sim.series("mole_fraction", "CarbonDioxide", location="top")   # NumPy: y.time, y.values
fit = skarstrom.fit_isotherm("sips", p_pa, q_molkg, temperature=298.15)
```

## Command line

```bash
skarstrom examples                                     # list the bundled decks
skarstrom examples run "Ward2022_Q100ccm" --json   # run one; JSON summary on stdout
skarstrom run deck.json --cycles 20 --set "step[Feed].time=120" --csv out/
skarstrom params deck.json [--all] [--json]            # settable parameters and their paths
```

## Documentation

<https://samaffleck.github.io/skarstrom/> — getting started, the Python API and CLI, tutorials, the model equations, numerics, the adsorbent library and the validation results. The sources are in [`docs/`](docs/); `examples/` has a runnable script for every deck family.

Build the C++ library and tests with the CMake presets:

```bash
cmake --preset linux-release        # or macos-release / windows
cmake --build --preset linux-release
ctest --preset linux-release -j3    # ~670 tests; writes ~4 GB to core/Results (delete afterwards)
```

## Validation

Every bundled paper deck is an automated test that computes an explicit error metric against numbers transcribed from the paper or curves digitised from its figures. The [validation page](https://samaffleck.github.io/skarstrom/validation/) is generated from those tests (48 cases from 10 publications; 194 asserted comparisons):

| Publication | Against | Median error | Worst error |
|---|---|---|---|
| Arora & Hasan (2021), O₂ PSA | paper model | 4.09 pp | 45.6 pp |
| Casas et al. (2012), CO₂/H₂ breakthrough | paper model, isotherm fit | 3.68 % (rel.) | 18.8 % (rel.) |
| Ferré et al. (2026), biogas PVSA | paper model | 1.11 pp / 7.37 % | 10.7 pp / 58.6 % |
| Nguyen et al. (2022), CALF-20 breakthrough | experiment | 9.94 % NRMSE | 12.4 % NRMSE |
| Nguyen et al. (2023), CALF-20 VPSA | experiment, paper model | 2.07 pp | 9.76 pp |
| Estupiñan Perez et al. (2019), 13X VSA | experiment, paper model | 4.28 pp | 11.1 pp |
| Streb & Mazzotti (2021), 13X breakthrough | experiment | 36.1 % NRMSE, 16.5 K | 82.9 % NRMSE, 32.6 K |
| Vysyaraju et al. (2024), He/N₂ PSA | experiment, paper model | 1.91 pp | 14 pp |
| Ward & Pini (2022), activated carbon breakthrough | experiment | 37.9 % NRMSE, 3.54 K | 38.2 % NRMSE, 4.21 K |
| Ward & Pini (2024), Mg-MOF-74 VSA | paper model | 1.79 pp | 2.3 pp |

Caveats: the tolerances are **regression bounds** set from the current model's measured error, not acceptance criteria, so a passing test can still be a large, documented disagreement (Arora 2021 net recovery, Streb 2021 curves). Many comparisons are against the authors' **own models** rather than measurements. Every deck is a prediction with the paper's published parameters — nothing is refitted per experiment — and the sampling lines between column and detector are not modelled. Each paper folder in `core/test/papers/` documents its inputs and the known causes of disagreement.

## Citing

If you use Skarstrom in published work, please cite it using [`CITATION.cff`](CITATION.cff) (GitHub shows a "Cite this repository" button), and cite the original sources of any adsorbent data or validation case you rely on.

## Get involved

We are forming a **steering committee** to guide Skarstrom's direction: which models and features come next, how validation cases and adsorbent data are curated, and how the project is governed as it grows. People from **academia and industry** are both welcome. If you would like to join, contact [Sam Affleck on LinkedIn](https://www.linkedin.com/in/sam-affleck).

## Contributing

Contributions are welcome — see [CONTRIBUTING.md](CONTRIBUTING.md) (DCO sign-off, tests, data-integrity rules) and the [Code of Conduct](CODE_OF_CONDUCT.md). Security issues: [SECURITY.md](SECURITY.md). Coding agents: [AGENTS.md](AGENTS.md).

## Licence

Skarstrom is licensed under the [Apache License 2.0](LICENSE); see also [NOTICE](NOTICE). It builds on third-party components under their own permissive licences, listed in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
