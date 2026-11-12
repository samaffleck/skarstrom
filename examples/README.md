# Examples

One runnable Python script per family of example decks. Each script **builds
every case from scratch in code**: the fluid, the adsorbents with their
isotherm and rate models, the column and its heat transfer, every cycle step
and boundary condition, the KPI definition and the solver settings, each value
written out with its unit. So each script shows everything you can change and
how, and is a template for a process of your own. It then runs the cases,
prints the KPIs or the cyclic-steady-state status, compares the result with
the paper's published values (`core/test/papers/*/reference_values.csv`) or
digitised curves (`core/data/`), and saves a plot to `examples/output/` when
matplotlib is installed.

Each script's `build(case)` returns the configured simulation without running
it (`import perez_2019; sim = perez_2019.build("A")`), and
`python/tests/test_examples_build.py` checks that it reproduces the matching
deck in `decks/` exactly.

```bash
pip install skarstrom matplotlib        # or, in a checkout: pip install . matplotlib
python examples/ward_2022.py --case 100 # about a second
python examples/perez_2019.py --case D1 # one VSA configuration to cyclic steady state
```

Options common to every script (`optimise_nsga2.py` and `perez_2019_case_study_optimise.py` have `--pop`, `--gen` and
`--workers` instead of `--case`; they need `pip install pymoo`): `--case NAME` (repeatable) to run a subset,
`--cycles N` to override the number of cycles for a quick look, `--no-plot`,
and `-q`. `_common.py` holds the shared helpers (deck lookup, reference data,
NRMSE, KPI tables, stream-log integration for decks without a KPI
configuration).

| Script | Decks | Process | Compared with | Runtime |
|---|---|---|---|---|
| [`ward_2022.py`](ward_2022.py) | `Ward2022_Q50ccm`, `_Q100ccm`, `_Q200ccm` | CO₂ breakthrough on activated carbon, three flow rates | digitised outlet y_CO₂ and thermocouple; paper-model breakthrough-time interval | 2 s |
| [`casas_2012.py`](casas_2012.py) | `Casas2012_Fig4_*`, `Casas2012_Fig6_*` (3) | CO₂/H₂ breakthrough on AP3-60 at 15 bar | paper-model breakthrough times (Table 5) | 8 s |
| [`nguyen_2022.py`](nguyen_2022.py) | `Nguyen2022_Fig4`, `_Fig4_50pct`, `_Fig4_100pct`, `Nguyen2022_Fig5`, `_Fig5_5pct`, `_Fig5_30pct` | CO₂/He and CO₂/N₂ breakthrough on CALF-20, six runs | measured outlet y_CO₂ (NRMSE, half-feed crossing time) and thermocouple (RMSE) | 2 s |
| [`streb_2021.py`](streb_2021.py) | `Streb2021_A1`–`_A7`, `_B1`, `_B2` | CO₂/CH₄ and H₂/CO₂/CH₄ breakthrough on 13X with thermocouples, nine runs | digitised outlet mole fractions and five thermocouples | 16 s (A1, A3; all nine about 3 min) |
| [`arora_2021.py`](arora_2021.py) | `Arora2021_LiLSX_PSA`, `_PVSA` | O₂ PSA and PVSA on LiLSX (50 cycles) | paper-model O₂ purity and net recovery (Table 4) | 2 min |
| [`ferre_2026.py`](ferre_2026.py) | `Ferre2026_NaYF` | Biogas CO₂/CH₄ PVSA on NaY, to CSS | paper-model purity and recovery of both products (Table 5) | 4 min |
| [`nguyen_2023.py`](nguyen_2023.py) | `Nguyen2023_Exp_I` … `_X` (10) | CO₂/N₂ VPSA on CALF-20, with and without LPP, to CSS | measured and paper-model CO₂ purity/recovery (Table 2) | 22 min |
| [`haghpanah_2013.py`](haghpanah_2013.py) | `Haghpanah2013_Sec61`, `_PL002` … `_PL010`, `_Fig5Benchmark` | CO₂/N₂ 4-step PVSA on 13X with compression and vacuum energy, to CSS | paper-model CO₂ purity/recovery (Haghpanah 2013 Sect. 6.1; Glover 2026 Table 4) | 4 min |
| [`perez_2019.py`](perez_2019.py) | `Perez2019_A`, `_B`, `_C1`, `_D1`, `_E`, `_F` | CO₂/N₂ VSA on 13X, with and without LPP, to CSS | measured and paper-model CO₂ purity/recovery (Table 3) | 6 min |
| [`perez_2019_case_study.py`](perez_2019_case_study.py) | `Perez2019_CaseStudy_Base`, `_Base_LPP` | Perez 2019 A and D1 run by devices (blower, valves with a ramped opening for pressurisation, vacuum pump) on a 3-bed schedule, to CSS | the validated decks' step-end pressures, purity, recovery and specific energy | 2 min |
| [`perez_2019_case_study_optimise.py`](perez_2019_case_study_optimise.py) | `Perez2019_CaseStudy_Base`, `_Base_LPP` | NSGA-II (pymoo) of the case-study device settings under plant constraints (purity/recovery targets, fluidisation, pressure drop, evacuation floor), warm-started, resumable from `results.csv`; `--summary DIR` writes the tables in `case_study_results/` | — (energy–productivity and purity–recovery fronts; the paper's Table 3 points) | about 15 min per deck and objective set (48 runs of about 45 s, 4 cores) |
| [`vysyaraju_2024.py`](vysyaraju_2024.py) | `Vysyaraju2024_PSA_exp1` | He/N₂ 5-step PSA on binderless 5A, to CSS | measured and paper-model He purity/recovery (Table 3) | 30 s |
| [`ward_pini_2024.py`](ward_pini_2024.py) | `Ward2024` | CO₂/N₂ 4-step VSA on Mg-MOF-74, to CSS | paper-model CO₂ purity/recovery (Table 4) | 80 s |
| [`tsa_illustrative.py`](tsa_illustrative.py) | `Illustrative_TSA_13X` | Illustrative 4-step TSA heated and cooled by a medium (jacket + internal exchanger), to CSS | — (KPIs, CSS status, bed temperatures) | 3 s |
| [`illustrative_dac_breakthrough.py`](illustrative_dac_breakthrough.py) | `Illustrative_DAC_Breakthrough` | Illustrative 400 ppm CO₂ breakthrough on Lewatit VP OC 1065 with dual-rate LDF and series-resistance kinetics | — (outlet times to 5, 50 and 95 % of the feed) | 2 s |
| [`five_step_product_purge.py`](five_step_product_purge.py) | `FiveStepProductPurge` | CO₂/N₂ VSA template with product purge | — (CSS status, integrator statistics) | 4 s |
| [`optimise_nsga2.py`](optimise_nsga2.py) | `Perez2019_A` | NSGA-II search (pymoo) of feed time and evacuation pressure, each candidate run to CSS in a process pool | — (the CO₂ purity / recovery Pareto front) | about 10 min (estimated: 60 runs of about 30 s, 4 cores) |
| [`calibrate_ward_2022.py`](calibrate_ward_2022.py) | `Ward2022_Q50ccm`, `_Q100ccm`, `_Q200ccm` | `skarstrom.calibrate`: the three breakthrough runs fitted together (shared q_s,2 and dispersion factor, one h_in per flow), the reproducible method of the decks' calibration | the hand calibration (`core/test/papers/ward_2022/README.md`) | 40 s (`--quick`: 3 s) |
| [`calibrate_cycle_kpis.py`](calibrate_cycle_kpis.py) | `Perez2019_A` | `skarstrom.calibrate` against per-cycle purity and recovery (KPI observables): a twin experiment that fits the CO₂ affinity and the bed voidage back | the values that made the data | 30 s (`--quick`: 13 s) |

Runtimes: wall-clock time of a full run of each script (all cases, to cyclic
steady state where the deck checks it) on one core of a shared 4-core Linux
machine, so treat them as approximate. They are the times of the decks'
grid-converged WENO settings and feed-loaded starts (2026-11-05, from the
paper tests, which run the same cases); a quick look can lower `cells` or
pass `--cycles`. Every script completed. `python/tests/test_examples.py`
imports all scripts and runs each with fast options (`pytest python/tests`).

How the comparisons are made: purity and recovery are those of the validation
tests in `core/test/papers/` (for the VSA decks, the evacuation product; for
the Arora decks, the collected product purity and the net recovery), so the
numbers printed here match the [validation page](../docs/validation.md). The
published values are mostly the authors' optimised **model** results or
**measurements**, as labelled.

## Notebooks

`notebooks/` holds five executed notebooks:

- `parametric_study.ipynb` — sweeps the feed-step time of a CO₂ VSA cycle and plots purity and recovery;
- `sensitivity_analysis.ipynb` — screening of uncertain inputs and KPI error bars with `skarstrom.sensitivity`;
- `isotherm_fit.ipynb` — single- and multi-temperature isotherm fits with the engine's fitter;
- `case_study_vsa.ipynb` — the CO₂ VSA case study: the device-driven base decks, their operating-limit
  and KPI reports, a small optimisation, and the paper, unconstrained and constrained purity–recovery
  fronts (from `case_study_results/`);
- `calibrate_ward_2022.ipynb` — `calibrate_ward_2022.py` step by step with the quick settings: the
  cases, the parameters, the fit, the comparison with the hand calibration, and the overlay and
  residual plots.

## Decks

To turn any deck (yours, or one written by a front end) into builder code
like these scripts, use `skarstrom.to_script(sim)` (see
[Building a simulation in Python](../docs/python.md#building-a-simulation-in-python)).

`decks/` holds the 47 example decks (also bundled in the Python package:
`skarstrom examples`). Their file names must stay shell-safe (checked by the
`DataFileNamesAreShellSafe` test). The test suite re-saves several of them; small
changes to the embedded numbers are expected and should be committed with the
change that caused them. See the
[examples page](https://samaffleck.github.io/skarstrom/examples/) for what each
deck reproduces.
