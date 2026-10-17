# Skarstrom

![Skarstrom](assets/logo.png#only-light){ width="420" }
![Skarstrom](assets/logo-dark.png#only-dark){ width="420" }

**Skarstrom** is an open-source simulator for pressure, vacuum and temperature swing adsorption (PSA, VSA and TSA) processes. It solves a one-dimensional finite-volume model of a packed adsorption column — coupled mass, momentum and energy balances with non-linear isotherms and mass-transfer kinetics — integrates it with SUNDIALS IDA, assembles multi-step, multi-bed cycles, runs them to cyclic steady state, and reports purity, recovery, productivity and energy.

The engine (`adsim`) is a C++20 library. The `skarstrom` Python package and command-line tool wrap it, and front ends (desktop or browser) can build it as a subproject.

```python
import skarstrom

sim = skarstrom.examples.Perez2019_A.load()
sim.steps["Feed"].time = "2 min"
result = sim.run()                                   # to cyclic steady state
print(result.status, result.css_cycle)
for k in sim.kpis(cycle=-1):
    print(k.component, k.purity, k.recovery)
```

---

## Why Skarstrom

There has been no good open-source tool for simulating adsorption processes, so researchers and engineers write their own column models from scratch or work in closed commercial packages. Results are then hard to reproduce or compare. Skarstrom aims to be the shared, tested foundation: a standard engine, deck format and cited adsorbent library, so the community can reproduce published work and build on each other's models.

## Start here

- [Getting Started](getting-started.md) — install the Python package or build from source, run a first deck from Python and from the command line.
- [Python and the Command Line](python.md) — the deck as objects, units, parameter paths, results, isotherms, the CLI and JSON output.
- [Examples](examples.md) — a runnable script for every bundled deck family, each compared with its paper.
- Tutorials: the [O₂ PSA (Arora 2021)](arora-2021-o2-psa.md), the [CO₂ VPSA with light-product pressurisation (Perez 2019)](perez-2019-vpsa.md) and the illustrative [TSA regeneration](tsa-regeneration.md).
- [Using the C++ library](cpp-library.md) — `find_package(Adsim)`, `add_subdirectory`, and the public headers.

## The model

The **Model** section documents the physics behind the solver:

- [Governing Equations](equations.md) — an overview of the coupled system.
- [Mass Balance](mass-balance.md), [Momentum Balance](momentum-balance.md), [Energy Balance](energy-balance.md) and [Wall Energy Balance](wall-energy-balance.md).
- [Isotherms](isotherms.md) — the equilibrium models and their parameters.
- [Kinetic Models](kinetic-models.md) — linear driving force and the mass-transfer coefficient correlations.
- [Transport Correlations](transport-correlations.md) — axial dispersion, thermal conductivity and heat transfer coefficients.
- [Equation of State and Fluid Properties](equation-of-state.md) — the gas model and CoolProp.

## Simulating a process

- [Cycles and Initial Conditions](cycles.md), [Boundary Conditions](boundary-conditions.md) and [Solver Settings](solver-settings.md).
- [Connecting Beds](connecting-beds.md), [Pressure Equalisation](pressure-equalization.md), [Light-Product Pressurisation](light-product-pressurization.md) and [Light Reflux / Recycle](light-reflux-recycle.md).
- [Numerics](numerics.md) and [Cyclic Steady State](cyclic-steady-state.md).
- [Performance Metrics](performance-metrics.md) — purity, recovery, productivity and energy.
- [Global Variables, Expressions and Units](global-variables.md).

## Studies and data

- [Parameterisation](parameterisation.md) — parametric studies and fitting deck parameters to data from Python.
- [Isotherm Fitting](isotherm-fitting.md) — fitting isotherm models to equilibrium data with the engine's fitter.
- [Adsorbent Library](adsorbent-library.md) — the curated, cited adsorbent data (generated page).
- [Validation](validation.md) — every comparison with published studies, generated from the test suite, and [Verification](verification.md) — the analytic and numerical checks.

## Project

Source, issues and releases: <https://github.com/samaffleck/skarstrom>. Skarstrom is licensed under the Apache License 2.0; see the repository's `LICENSE`, `NOTICE` and `THIRD_PARTY_NOTICES.md`. If you use it in published work, please cite it (`CITATION.cff`).

**Get involved.** We are forming a steering committee to guide the project's direction; people from academia and industry are welcome. To join, contact [Sam Affleck on LinkedIn](https://www.linkedin.com/in/sam-affleck).
