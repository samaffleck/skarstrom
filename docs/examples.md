# Examples

Skarstrom ships 33 example decks in `examples/decks/`, most of them reproducing published adsorption studies. The Python package bundles all of them (`skarstrom examples`, `skarstrom.examples.<name>`), so they are the quickest way to see a complete process: run one, inspect how it is set up (`skarstrom params deck.json --all`), and adapt it.

Every paper-based deck is also an automated validation test; the [Validation](validation.md) page lists the simulated and published values and the error metric for each case, including where Skarstrom and the paper disagree, and why.

## Example scripts

`examples/` has one runnable script per deck family. Each **builds its cases from scratch in code** (fluid, adsorbents, column, every step and boundary condition, KPIs and solver settings, each value with its unit; `build(case)` returns the simulation), runs them, prints the KPIs or the cyclic-steady-state status, compares the result with the paper's published values (`core/test/papers/*/reference_values.csv`) or digitised curves (`core/data/`) when run from a repository checkout, and saves a plot to `examples/output/` when matplotlib is installed (`--no-plot` to skip). Common options: `--case NAME` (repeatable) to run a subset, `--cycles N` for a quick look.

```bash
pip install "skarstrom[examples]"        # or pip install . && pip install matplotlib
python examples/ward_2022.py --case 100  # the fastest: about a second
```

| Script | Decks | Process | Compared with | Runtime |
|---|---|---|---|---|
| [`ward_2022.py`](https://github.com/samaffleck/skarstrom/blob/main/examples/ward_2022.py) | `Ward2022_Q50ccm`, `_Q100ccm`, `_Q200ccm` | CO₂ breakthrough on activated carbon, three flow rates | digitised outlet y_CO₂ and thermocouple; paper-model breakthrough-time interval | 3 s |
| [`casas_2012.py`](https://github.com/samaffleck/skarstrom/blob/main/examples/casas_2012.py) | `Casas2012_Fig4_*`, `Casas2012_Fig6_*` (3) | CO₂/H₂ breakthrough on AP3-60 at 15 bar | paper-model breakthrough times (Table 5) | 2 s |
| [`nguyen_2022.py`](https://github.com/samaffleck/skarstrom/blob/main/examples/nguyen_2022.py) | `Nguyen2022_Fig4`, `_Fig5` | CO₂/He breakthrough on CALF-20 | digitised outlet y_CO₂ (NRMSE, half-feed crossing time) | 1 s |
| [`streb_2021.py`](https://github.com/samaffleck/skarstrom/blob/main/examples/streb_2021.py) | `Streb2021_A1`, `_A3` | CO₂/CH₄ breakthrough on 13X with thermocouples | digitised outlet mole fractions and five thermocouples | 5 s |
| [`arora_2021.py`](https://github.com/samaffleck/skarstrom/blob/main/examples/arora_2021.py) | `Arora2021_LiLSX_PSA`, `_PVSA` | O₂ PSA and PVSA on LiLSX (50 cycles) | paper-model O₂ purity and net recovery (Table 4) | 16 s |
| [`ferre_2026.py`](https://github.com/samaffleck/skarstrom/blob/main/examples/ferre_2026.py) | `Ferre2026_NaYF` | Biogas CO₂/CH₄ PVSA on NaY, to CSS | paper-model purity and recovery of both products (Table 5) | 22 s |
| [`nguyen_2023.py`](https://github.com/samaffleck/skarstrom/blob/main/examples/nguyen_2023.py) | `Nguyen2023_Exp_I` … `_X` (10) | CO₂/N₂ VPSA on CALF-20, with and without LPP, to CSS | measured and paper-model CO₂ purity/recovery (Table 2) | 13 min |
| [`perez_2019.py`](https://github.com/samaffleck/skarstrom/blob/main/examples/perez_2019.py) | `Perez2019_A`, `_B`, `_C1`, `_D1`, `_E`, `_F` | CO₂/N₂ VSA on 13X, with and without LPP, to CSS | measured and paper-model CO₂ purity/recovery (Table 3) | 3 min |
| [`vysyaraju_2024.py`](https://github.com/samaffleck/skarstrom/blob/main/examples/vysyaraju_2024.py) | `Vysyaraju2024_PSA_exp1` | He/N₂ 5-step PSA on binderless 5A, to CSS | measured and paper-model He purity/recovery (Table 3) | 22 s |
| [`ward_pini_2024.py`](https://github.com/samaffleck/skarstrom/blob/main/examples/ward_pini_2024.py) | `Ward2024` | CO₂/N₂ 4-step VSA on Mg-MOF-74, to CSS | paper-model CO₂ purity/recovery (Table 4) | 18 s |
| [`tsa_illustrative.py`](https://github.com/samaffleck/skarstrom/blob/main/examples/tsa_illustrative.py) | `Illustrative_TSA_13X` | Illustrative 4-step TSA with heater and cooler, to CSS | — (KPIs, CSS status, bed temperatures) | 3 s |
| [`five_step_product_purge.py`](https://github.com/samaffleck/skarstrom/blob/main/examples/five_step_product_purge.py) | `FiveStepProductPurge` | CO₂/N₂ VSA template with product purge | — (CSS status, integrator statistics) | 2 s |

Runtimes are wall-clock times of a full run of the script (all its cases, to cyclic steady state where the deck checks it) on one core of a shared 4-core Linux machine, measured for the 0.1.0 release. `python/tests/test_examples.py` imports every script and runs each once with fast options (one breakthrough case, or two cycles of one cyclic case).

Three executed notebooks are in `examples/notebooks/`: a parametric study of the feed time of a CO₂ VSA cycle, a finite-difference sensitivity analysis of its KPIs, and isotherm fitting with the engine's fitter.

## Worked tutorials

- [O₂ PSA (Arora 2021)](arora-2021-o2-psa.md): a four-step Skarstrom cycle for medical oxygen on LiLSX, and how to read purity and recovery when product is recycled.
- [CO₂ VPSA with LPP (Perez 2019)](perez-2019-vpsa.md): a five-step vacuum swing cycle on zeolite 13X with light-product pressurisation, run to cyclic steady state and compared with experiment.
- [TSA regeneration](tsa-regeneration.md): an illustrative temperature swing cycle with heater and cooler energy sources.

## Bundled example decks

| Family | Decks | What they are |
|---|---|---|
| **Arora (2021)** | `Arora2021_LiLSX_PSA`, `Arora2021_LiLSX_PVSA` | O₂ from air on LiLSX zeolite, 4-step PSA and vacuum PSA (Table 4) |
| **Casas (2012)** | `Casas2012_Fig4_50_50_15bar_25C`, `Casas2012_Fig6_25_75_15bar_25C`, `Casas2012_Fig6_75_25_15bar_25C` | CO₂/H₂ breakthrough on activated carbon AP3-60 at 15 bar and 25 °C, three feed compositions (Figs. 4 and 6) |
| **Ferré (2026)** | `Ferre2026_NaYF` | CO₂/CH₄ biogas upgrading, 4-step PVSA on zeolite NaY (Table 5) |
| **Nguyen (2022)** | `Nguyen2022_Fig4`, `Nguyen2022_Fig5` | CO₂/He breakthrough (70 % and 15 % CO₂) into a helium-filled bed of CALF-20 (Figs. 4 and 5) |
| **Nguyen (2023)** | `Nguyen2023_Exp_I` … `Nguyen2023_Exp_X` | CO₂/N₂ VPSA on structured CALF-20: experiments I–V use the basic 4-step cycle, VI–X the 4-step cycle with light-product pressurisation (Table 2) |
| **Perez (2019)** | `Perez2019_A`, `_B`, `_C1`, `_D1`, `_E`, `_F` | CO₂/N₂ VSA on zeolite 13X: A, B and C1 are the basic 4-step cycle, D1, E and F add light-product pressurisation (Table 3) |
| **Streb (2021)** | `Streb2021_A1`, `Streb2021_A3` | CO₂/CH₄ breakthrough on zeolite 13X with thermocouple data |
| **Vysyaraju (2024)** | `Vysyaraju2024_PSA_exp1` | He/N₂ separation on binderless zeolite 5A, 5-step PSA |
| **Ward (2022)** | `Ward2022_Q50ccm`, `Ward2022_Q100ccm`, `Ward2022_Q200ccm` | CO₂ breakthrough on activated carbon at three feed flow rates |
| **Ward (2024)** | `Ward2024` | CO₂/N₂ 4-step VSA on Mg-MOF-74 (Table 4) |
| **Illustrative** | `Illustrative_TSA_13X` | **Not a reproduction.** A four-step indirectly heated TSA for CO₂ on zeolite 13X, built for the [TSA tutorial](tsa-regeneration.md) |
| **Other** | `FiveStepProductPurge` | A CO₂/N₂ VSA cycle with feed, blowdown, evacuation, pressurisation and product purge on a 13X-type isotherm. It comes from an integration test, has no KPI configuration and no literature source; use it as a template, not as data |

Several decks also carry the experimental curves digitised from their paper (the deck's experimental-data section), the same data as `core/data/`.

## Working with an example

```python
import skarstrom

path = skarstrom.examples.Perez2019_A.copy_to(".")   # an editable copy
sim = skarstrom.Simulation(path)
print(sim.steps["Feed"])                  # the step's values, with units and string paths
sim.steps["Feed"].time = "80 s"
sim.run()
sim.save("Perez2019_A_80s.json")
```

The cycle templates in `data/cycles/` (breakthrough, 4-step VPSA with and without LPP, Skarstrom cycle, 6-step VPSA with rinse) are starting points for building a new process in C++ (see [Cycles](cycles.md)).
