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

**Skarstrom** is an open-source library for adsorption process simulations including **pressure, vacuum and temperature swing adsorption** (PSA, VSA and TSA). Written in C++, the engine solves a 1D finite-volume model of a packed-bed adsorption column with multiple layers. Skarstrom is designed for engineers and scientists to rapidly develop new adsorption-based gas separation systems, with modules for accelerating convergence to cyclic steady state and for optimising purity, recovery, productivity and energy consumption.

## Quick start: optimise a VSA cycle from the command line

This example takes one of the 43 bundled [examples](https://samaffleck.github.io/skarstrom/examples/), [`Perez2019_A`](https://github.com/samaffleck/skarstrom/blob/main/examples/decks/Perez2019_A.json). It is a 4-step vacuum swing adsorption (VSA) cycle that captures CO₂ from a 15 % CO₂ / 85 % N₂ flue gas on zeolite 13X. You can find the feed time and vacuum level that give the most CO₂ for a purity target.

**1. Install** (wheels for CPython 3.10–3.15 on Linux, macOS and Windows):

```bash
pip install skarstrom matplotlib
```

**2. Run the example to cyclic steady state:**

```console
$ skarstrom examples run Perez2019_A
Deck:    Perez2019_A
Status:  reached_css (cyclic steady state reached at cycle 49)
Cycles:  49 run in 40.00 s (160801 integrator steps)

KPIs at cycle 49:
component      purity [mol/mol]  recovery [mol/mol]  productivity [kg/(kg s)]  ...
CarbonDioxide  0.995067          0.596507            3.6742e-05                ...
```

**3. Find the design variables, then sweep them.** Every value in a deck has a stable path, and `skarstrom params` lists them:

```console
$ skarstrom params Perez2019_A --search Evacuation      # columns trimmed
path                                          value             unit  alias
step[3].time                                  57                s     step[Evacuation].time
endpoint[4].pressure_double_exponential_end   3.2e+03 (= 3200)  Pa    step[Evacuation].bottom.pressure_double_exponential_end
...
```

`--set` changes any value on the command line (units are checked), and `--json` writes a machine-readable summary. This loop runs 18 cases (6 feed times × 3 vacuum levels), three at a time:

```bash
for vac in 2 4 6; do
  for feed in 60 80 100 120 140 160; do
    skarstrom examples run Perez2019_A -q --json \
      --set "step[Feed].time=$feed" \
      --set "step[Evacuation].bottom.pressure_double_exponential_end=$vac [kPa]" \
      > run_${vac}kPa_${feed}s.json
  done &
done
wait
```

**4. Plot the results.** Each JSON file holds the KPIs of every cycle; the last entry is the converged cycle:

```python
import glob, json, re
import matplotlib.pyplot as plt

runs = {}
for f in glob.glob("run_*.json"):
    vac, feed = map(int, re.findall(r"\d+", f))
    co2 = json.load(open(f))["kpis"][-1]          # CO2 KPIs of the last cycle
    runs.setdefault(vac, []).append((feed, co2["recovery"], co2["purity"]))

fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4), sharex=True)
for vac, points in sorted(runs.items()):
    feed, recovery, purity = zip(*sorted(points))
    left.plot(feed, recovery, "o-", label=f"{vac} kPa")
    right.plot(feed, purity, "o-")
left.set(xlabel="feed time [s]", ylabel="CO2 recovery")
right.set(xlabel="feed time [s]", ylabel="CO2 purity")
left.legend(title="vacuum")
fig.tight_layout()
fig.savefig("sweep.png", dpi=150)
```

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/assets/readme-sweep-dark.png">
    <img src="docs/assets/readme-sweep-light.png" alt="CO2 recovery and purity against feed time for vacuum levels of 2, 4 and 6 kPa. Deeper vacuum raises recovery from about 20 % to about 70 %; a longer feed raises purity and lowers recovery. The best point with at least 99 % purity is 2 kPa and a 120 s feed, at 70 % recovery." width="100%">
  </picture>
</p>

A deeper vacuum is the strongest lever on recovery. The feed time trades purity against recovery. For a 99 % purity target, the best of these 18 runs is a 2 kPa vacuum with a 120 s feed: 99.5 % purity at 70 % recovery. That is 10 points more recovery than the published cycle. All 18 runs reached cyclic steady state; the sweep took about 5 minutes on 4 cores. (The figure is a styled version of the plot above, drawn by [`scripts/docs/readme_figures.py`](scripts/docs/readme_figures.py).)

## Use Skarstrom with AI agents

Skarstrom is built to be driven by AI coding agents (Claude Code, Codex, Cursor, Copilot and others) as well as by people:

- **[`AGENTS.md`](AGENTS.md)** tells an agent how to use the package: check the version, run a few cycles first, and treat a run-health warning as a wrong answer until it can explain it.
- **[`llms.txt`](https://samaffleck.github.io/skarstrom/llms.txt)** and **[`llms-full.txt`](https://samaffleck.github.io/skarstrom/llms-full.txt)** give the whole documentation as plain Markdown, for the latest release (and under [`/dev/`](https://samaffleck.github.io/skarstrom/dev/llms.txt) for `main`).
- **Machine-readable output.** `--json` on `run`, `examples`, `params` and `sensitivity`, with exit codes 0 (done), 1 (failed) and 130 (interrupted).
- **A parameter registry.** Every value has a stable path, a unit and bounds, so an agent changes `step[Feed].time` instead of editing deck JSON. Unknown paths and wrong units are errors, not silent typos.
- **Run-health checks.** Each run reports structured warnings (`MOLE_FRACTION_SUM`, `MASS_BALANCE_STUCK`, `STEP_RETRIES`, `HISTORY_TRUNCATED`) with a value, a threshold and a hint, so an agent can tell a converged answer from a numerical artefact.
- **Documentation that runs.** CI executes every code block in the docs and in `AGENTS.md` against the built package, so the calls an agent reads are calls the engine accepts.

### Starter prompt

Copy this into any agent that can run shell commands, and replace the task at the end with your own:

```text
I want to use Skarstrom, an open-source simulator for pressure, vacuum and temperature swing adsorption (PSA/VSA/TSA), to answer a process design question.

Set-up:
- Install it with `pip install skarstrom matplotlib` (skip if already installed) and run `skarstrom --version`.
- Read the agent guide at https://github.com/samaffleck/skarstrom/blob/main/AGENTS.md and the docs index at https://samaffleck.github.io/skarstrom/llms.txt. Use the docs that match the installed version.

How to work:
- Use the `skarstrom` CLI with `--json -q`, or the Python API. Find parameter paths with `skarstrom params <deck> --search <text>`. Do not edit deck JSON by hand.
- Check a set-up with a short run (`--cycles 3`) before you run to cyclic steady state (CSS).
- After every run, check `status` and `warnings`. Do not trust a number until you can explain any warning.
- Give units, and say which numbers come from the simulation and which are your assumptions.

Task: Start from the bundled example `Perez2019_A` (a 4-step VSA that captures CO₂ from 15 % CO₂ / 85 % N₂ on zeolite 13X). Find the feed step time (60–160 s) and evacuation pressure (2–6 kPa) that give the highest CO₂ recovery with at least 95 % CO₂ purity at CSS. Plot the trade-off, mark your recommended point, and summarise the result in a short table.
```


## Graphical user interface

We host a user interface on the [Skarstrom website](https://www.skarstrom.co.uk/). It is a browser-based interface that runs the C++ engine and is the best starting point for users who don't want to write code. From the interface you can run examples, visualise the results, and run optimisation and parameterisation studies.

[![The Skarstrom web interface: the simulation tree, solver settings and a plotted pressure profile of the Arora (2021) PSA example](docs/assets/skarstrom-web.png)](https://www.skarstrom.co.uk/)

## Why Skarstrom

Several commercial simulation packages are available, but to date there has been no good open-source tool for adsorption modelling. It has become common practice for people to write their own models from scratch, resulting in a lot of repeated work and a lack of transparency and reproducibility.

Skarstrom aims to bridge this gap by being the first comprehensive open project for adsorption researchers to collaborate on. We also showcase exceptional solver speed and performance, and hope this can be improved further by collaboration with the community.

## Performance

Cycle optimisation and parameter fitting are important design tasks that require running to cyclic steady state (CSS) hundreds or thousands of times, so Skarstrom has been optimised for performance to make these tasks practical on common hardware. There are no well-established performance benchmarks; however, the table below can be used as a guide to the expected time to reach CSS. You are encouraged to run the examples yourself and benchmark them against your own code.

| Example | Cells | Cycles to CSS | Time to CSS |
|---|---|---|---|
| [Basic 4-step VSA, 15 % CO₂ / N₂ on zeolite 13X (`Perez2019_A`)](https://github.com/samaffleck/skarstrom/blob/main/examples/decks/Perez2019_A.json) | 30 | 50 | 4.4 s |
| [4-step VPSA with light-product pressurisation, CO₂ / N₂ on CALF-20 (`Nguyen2023_Exp_VI`)](https://github.com/samaffleck/skarstrom/blob/main/examples/decks/Nguyen2023_Exp_VI.json) | 100 | 68 | 16.5 s |
| [4-step PVSA for biogas upgrading, 40 % CO₂ / CH₄ on zeolite NaY (`Ferre2026_NaYF`)](https://github.com/samaffleck/skarstrom/blob/main/examples/decks/Ferre2026_NaYF.json) | 30 | 30 | 3.1 s |
| [4-step VSA, 15 % CO₂ / N₂ on Mg-MOF-74 (`Ward2024`)](https://github.com/samaffleck/skarstrom/blob/main/examples/decks/Ward2024.json) | 30 | 55 | 4.4 s |
| [4-step VSA, 15 % CO₂ / N₂ on Mg-MOF-74 (`Ward2024`), finer grid](https://github.com/samaffleck/skarstrom/blob/main/examples/decks/Ward2024.json) | 90 | 52 | 14.6 s |

*CSS is reached when every component's mass balance over a cycle closes within 0.5 %. Median of 3 runs on one core of an Apple M4, Release build; reproduce with `python scripts/benchmark_css.py`.*

## Documentation

The full documentation is at <https://samaffleck.github.io/skarstrom/>: getting started, the Python API and command line, tutorials, the model equations, the adsorbent library and the validation results.

## Citing

If you use Skarstrom in published work, please cite it using [`CITATION.cff`](CITATION.cff).

## Get involved

If you would like to be a part of the **steering committee** to guide Skarstrom's direction, please contact [Sam Affleck on LinkedIn](https://www.linkedin.com/in/sam-affleck). People from **academia and industry** are both welcome.

## Contributing

Contributions are welcome — see [CONTRIBUTING.md](CONTRIBUTING.md) (DCO sign-off, tests, data-integrity rules) and the [Code of Conduct](CODE_OF_CONDUCT.md). Security issues: [SECURITY.md](SECURITY.md). Coding agents: [AGENTS.md](AGENTS.md).

## Licence

Skarstrom is licensed under the [Apache License 2.0](LICENSE); see also [NOTICE](NOTICE). It builds on third-party components under their own permissive licences, listed in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
