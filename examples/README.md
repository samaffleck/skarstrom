# Examples

One runnable Python script per family of example decks. Each script loads its
decks (from `decks/`, or from the installed package when run elsewhere), runs
them, prints the KPIs or the cyclic-steady-state status, compares the result
with the paper's published values (`core/test/papers/*/reference_values.csv`)
or digitised curves (`core/data/`), and saves a plot to `examples/output/`
when matplotlib is installed.

```bash
pip install skarstrom matplotlib        # or, in a checkout: pip install . matplotlib
python examples/ward_2022.py --case 100 # about a second
python examples/perez_2019.py --case D1 # one VSA configuration to cyclic steady state
```

Options common to every script: `--case NAME` (repeatable) to run a subset,
`--cycles N` to override the number of cycles for a quick look, `--no-plot`,
and `-q`. `_common.py` holds the shared helpers (deck lookup, reference data,
NRMSE, KPI tables, stream-log integration for decks without a KPI
configuration).

| Script | Decks | Process | Compared with | Runtime |
|---|---|---|---|---|
| [`ward_2022.py`](ward_2022.py) | `Ward2022_Q50ccm`, `_Q100ccm`, `_Q200ccm` | CO₂ breakthrough on activated carbon, three flow rates | digitised outlet y_CO₂ and thermocouple; paper-model breakthrough-time interval | 3 s |
| [`casas_2012.py`](casas_2012.py) | `Casas2012_Fig4_*`, `Casas2012_Fig6_*` (3) | CO₂/H₂ breakthrough on AP3-60 at 15 bar | paper-model breakthrough times (Table 5) | 2 s |
| [`nguyen_2022.py`](nguyen_2022.py) | `Nguyen2022_Fig4`, `_Fig5` | CO₂/He breakthrough on CALF-20 | digitised outlet y_CO₂ (NRMSE, half-feed crossing time) | 1 s |
| [`streb_2021.py`](streb_2021.py) | `Streb2021_A1`, `_A3` | CO₂/CH₄ breakthrough on 13X with thermocouples | digitised outlet mole fractions and five thermocouples | 5 s |
| [`arora_2021.py`](arora_2021.py) | `Arora2021_LiLSX_PSA`, `_PVSA` | O₂ PSA and PVSA on LiLSX (50 cycles) | paper-model O₂ purity and net recovery (Table 4) | 16 s |
| [`ferre_2026.py`](ferre_2026.py) | `Ferre2026_NaYF` | Biogas CO₂/CH₄ PVSA on NaY, to CSS | paper-model purity and recovery of both products (Table 5) | 22 s |
| [`nguyen_2023.py`](nguyen_2023.py) | `Nguyen2023_Exp_I` … `_X` (10) | CO₂/N₂ VPSA on CALF-20, with and without LPP, to CSS | measured and paper-model CO₂ purity/recovery (Table 2) | 13 min |
| [`perez_2019.py`](perez_2019.py) | `Perez2019_A`, `_B`, `_C1`, `_D1`, `_E`, `_F` | CO₂/N₂ VSA on 13X, with and without LPP, to CSS | measured and paper-model CO₂ purity/recovery (Table 3) | 3 min |
| [`vysyaraju_2024.py`](vysyaraju_2024.py) | `Vysyaraju2024_PSA_exp1` | He/N₂ 5-step PSA on binderless 5A, to CSS | measured and paper-model He purity/recovery (Table 3) | 22 s |
| [`ward_pini_2024.py`](ward_pini_2024.py) | `Ward2024` | CO₂/N₂ 4-step VSA on Mg-MOF-74, to CSS | paper-model CO₂ purity/recovery (Table 4) | 18 s |
| [`tsa_illustrative.py`](tsa_illustrative.py) | `Illustrative_TSA_13X` | Illustrative 4-step TSA with heater and cooler, to CSS | — (KPIs, CSS status, bed temperatures) | 3 s |
| [`five_step_product_purge.py`](five_step_product_purge.py) | `FiveStepProductPurge` | CO₂/N₂ VSA template with product purge | — (CSS status, integrator statistics) | 2 s |

Runtimes: wall-clock time of a full run of each script (all cases, to cyclic
steady state where the deck checks it), measured on one core of a shared
4-core Linux machine for the 0.1.0 release with other jobs running, so treat
them as approximate. Every script completed. `python/tests/test_examples.py`
imports all scripts and runs each with fast options (`pytest python/tests`).

How the comparisons are made: purity and recovery are those of the validation
tests in `core/test/papers/` (for the VSA decks, the evacuation product; for
the Arora decks, the collected product purity and the net recovery), so the
numbers printed here match the [validation page](../docs/validation.md). The
published values are mostly the authors' optimised **model** results or
**measurements**, as labelled.

## Notebooks

`notebooks/` holds three executed notebooks:

- `parametric_study.ipynb` — sweeps the feed-step time of a CO₂ VSA cycle and plots purity and recovery;
- `sensitivity_analysis.ipynb` — local finite-difference sensitivities of the KPIs to tunable parameters;
- `isotherm_fit.ipynb` — single- and multi-temperature isotherm fits with the engine's fitter.

## Decks

`decks/` holds the 33 example decks (also bundled in the Python package:
`skarstrom examples`). Their file names must stay shell-safe (checked by the
`DataFileNamesAreShellSafe` test). The test suite re-saves several of them; small
changes to the embedded numbers are expected and should be committed with the
change that caused them. See the
[examples page](https://samaffleck.github.io/skarstrom/examples/) for what each
deck reproduces.
