# Paper validation tests

Each subfolder reproduces the experiments / cycle simulations of **one academic
publication** and compares the result with the paper **quantitatively**: every
case computes an explicit error metric against numbers transcribed from the
paper or digitised from its figures, asserts it against a tolerance, and writes
it to a machine-readable summary from which the docs page
[`docs/validation.md`](../../../docs/validation.md) is generated.

Every paper goes through the seven-step [validation protocol](#validation-protocol)
below; a paper counts as validated only when all seven steps are done and
written down.

## Folder convention

```
papers/
├── papers.json               # citation / DOI / system for every folder (Crossref-checked)
├── validation_summary.json   # GENERATED: merged metrics of the last validation run
└── <author_year>/
    ├── <author_year>.cpp     # GoogleTest cases — one TEST_F per reproduced case
    ├── <author_year>.h       # (optional) paper fixture + shared setup helper
    ├── reference_values.csv  # published numbers used by the tests (with table/figure/page)
    ├── README.md             # citation, DOI, reproduced cases, provenance table,
    │                         # model-form audit, attribution runs, measured errors
    └── assets/               # git-ignored: the paper PDF, SI, local text extractions
```

Experimental curves live outside this tree, one folder per paper in
`core/data/<paper>/`, each with a `README.md` recording the provenance of
every file (see [step 4](#4-data)).

- **Folder name**: `surname_year` (e.g. `perez_2019`, `arora_hasan_2021`). For
  two papers by the same first author in one year, add a discriminator.
- **Fixture**: one GoogleTest suite per paper, named `<Paper>Test`
  (e.g. `Perez2019Test`), deriving from `UserDataTest`.
- **Case name**: the figure / table / configuration reproduced
  (e.g. `ConfigD1`, `Fig5_TernaryBreakthrough`, `Q100ccm`).
- **`SaveResultsAsExample(...)` labels are stable** so the bundled example
  decks (`examples/decks/`) keep their names.

## Validation protocol

Each step's result is written into the paper's `README.md` (and, through
`validation_summary.json`, the generated validation page). A step that has not
been done is marked as open in the README, not skipped silently.

### 1. Source

The paper PDF and its supplementary information are read **in full**. They are
kept locally in the git-ignored `core/test/papers/<paper>/assets/` and never
committed (see [below](#copyrighted-paper-material-is-not-kept-in-the-repo)).
A paper whose text has never been read does not count as validated, however
good the agreement: the README says plainly when the text was not available
(as several READMEs do today).

### 2. Provenance audit

Every deck value is traced to a table, an equation or a sentence of the paper,
in a provenance table in the README (template [below](#provenance-table-template)).

- **Assumptions** (values the paper does not state) are listed. Each one is
  either closed from the text (and becomes an ordinary row) or kept and
  labelled as an assumption, with the reason for the value chosen.
- **Derived values** (a density from a mass and a volume, a voidage from two
  densities, a unit conversion) show their arithmetic in the note.
- **Deviations** from the paper (a deck value or form that differs from what the
  paper used) are in bold, so they are the first thing a reader sees.

### 3. Model-form audit

The model is compared line by line with the paper's, not only its parameters:

- isotherm form and the competition (mixing) rule;
- rate model, and how its coefficient is obtained (fitted, correlation, which one);
- axial dispersion, axial thermal conduction, film and wall heat-transfer closures;
- pressure-drop law;
- wall boundary conditions (ends and outside);
- gas-property treatment (constant values vs CoolProp mixture properties) and
  the carrier / initial gas;
- anything else the paper's balances contain (sampling lines, dead volumes,
  heat of adsorption form).

Where the engine lacks the paper's form, the form is added to the engine as a
selectable option; it is never approximated silently. Where the deck picked a
different option than the paper, the deck is changed. A form that cannot be
matched yet stays a bold deviation in the provenance table, with an open item.

### 4. Data

Experimental curves are re-checked against the figure: which figure and panel,
axis ranges and units, number of points. Their provenance is recorded per file
in `core/data/<paper>/README.md`. **Measured data supplied with the paper (SI
tables, data repositories) replace digitised points wherever they exist.**

This differs from the rule for the [adsorbent library](../../../data/adsorbents/)
(`AGENTS.md`): library parameters and `reference_values.csv` hold only numbers
printed in the paper, never values read off a plot. Paper-validation curves may
be digitised from figures when no measured data exist, and their README must
say that they are digitised, from which figure, and how.

### 5. Metrics

The paper's own indicator definitions are used where it gives them (e.g.
Ward & Pini 2022, Table 9: t_br at y_out = 0.5 and ΔT_ads against the feed
temperature). The standard set, each against the experiment and, where the
paper reports its own simulations, against the paper's model (code-to-code):

| Experiment type | Metrics |
|---|---|
| Breakthrough | front time (y = y_feed/2 unless the paper defines it otherwise), outlet curve NRMSE, thermocouple RMSE, peak temperature (rise) |
| Cycles | purity, recovery, productivity, specific energy |

How each metric is computed and asserted is under
[Quantitative validation](#quantitative-validation).

### 6. Attribution

Each remaining disagreement is attributed by one-factor sensitivity runs:
change one thing at a time (the carrier gas, the dispersion closure, the heat
of adsorption, the wall end condition, …), record the size of its effect on
each metric, and write the table into the README. The README then says *why*
the model differs from the paper, not only *that* it does. The Ward 2022 and
Streb 2021 READMEs ("Attribution runs") are worked examples.

### 7. Report

Tolerances are reset from the final errors (regression bounds, see below), the
validation page is regenerated, and the docs site gets a per-paper validation
page summarising steps 1–6.

### Order

1. Tier 1, breakthrough papers: Ward & Pini 2022, Streb & Mazzotti 2021,
   Nguyen et al. 2022, Casas et al. 2012.
2. Tier 2, CO₂/N₂ cycle papers: Estupiñan Perez et al. 2019, Nguyen et al.
   2023, Ward & Pini 2024.

## Provenance table template

One row per deck value. The deck parameter is the stable registry path
(`skarstrom params deck.json --all` lists them; IDs are the deck's IDs), so a
reader can check the value in the deck and a script can sweep it.

```markdown
| Value | Deck parameter (registry path) | Source (table / equation / page / sentence) | Note (derived arithmetic, assumption status) |
|---|---|---|---|
```

Filled example (two rows from Streb & Mazzotti 2021):

| Value | Deck parameter (registry path) | Source (table / equation / page / sentence) | Note (derived arithmetic, assumption status) |
|---|---|---|---|
| 1.2 m | `layer[1].length` | Table 1, column length | |
| 0.452 | `layer[1].bed_void_fraction` | Table 1, particle density 1085 kg/m³ and bed density 595 kg/m³ | derived: ε_b = 1 − 595/1085 = 0.452 |

An assumption row says so in the note, e.g. "**assumption** (not stated in the
paper; checked: Sect. 2, Table 1, SI): 1.5 mm, typical pellet size", and a row
that differs from the paper's model is in bold.

## Quantitative validation

Helpers live in [`../fixtures/validation_metrics.h`](../fixtures/validation_metrics.h):

```cpp
validation::CaseReport report("perez_2019", "ConfigA");
report.CheckAbs("CO2 purity vs experiment", purity,
                validation::Reference("perez_2019", "ConfigA", "CO2_purity_exp"), /*tol pp*/ 1.5);
report.CheckCurve("outlet y_CO2", validation::FaceMoleFraction(solver, "CarbonDioxide", 0, 0),
                  validation::LoadSeriesCsv("data/streb_2021/A3/A3_mole_fraction_CO2.csv"),
                  /*tol*/ 0.10, "-", "Fig. 6 (digitised)", /*normalise by y_feed*/ 0.5);
```

- **Metrics.** Breakthrough / outlet curves: RMSE at the digitised points,
  normalised by the feed mole fraction (NRMSE); thermocouples: RMSE in K.
  Purity / recovery: absolute error in percentage points. Productivity,
  breakthrough time: relative error. `Info(...)` records a comparison without
  asserting (e.g. approximate numbers read from the text).
- **Reference values** are read at run time from `reference_values.csv`
  (`case,quantity,value,unit,source,note`). Only numbers actually printed in the
  paper (tables, values stated in the text) go there; anything derived is marked
  as such in the `note`. Never eyeball values off a plot.
- **Tolerances are regression bounds, not targets.** Each is set from the
  measured error of the current model, rounded up modestly, with a comment
  giving the measured value. Never loosen one to hide a disagreement: if the
  model and the paper disagree, keep the honest bound and document the reason
  in the paper README.
- **Summary.** Each case writes `<paper>__<case>.json` into
  `$ADSIM_VALIDATION_DIR` (CMake points ctest at `<build>/validation`; without
  the variable, `./Results/validation`).

Regenerate the summary and the docs page after a run:

```bash
ctest --preset linux-release -R "(Arora2021|Casas2012|Ferre2026|Haghpanah2013|Nguyen2022|Nguyen2023|Perez2019|Streb2021|Vysyaraju2024|Ward2022|WardPini2024)" -j2
python scripts/docs/gen_validation_page.py --merge build/linux-release/validation
```

`python scripts/docs/gen_validation_page.py` (no arguments) re-renders the page
from the committed `validation_summary.json` without re-running the tests, and
`--check` (run by CI) fails if the page is out of date.

## Copyrighted paper material is not kept in the repo

Paper PDFs, supplementary files (DOCX, …), figure images and full-text
extractions of papers **must not be committed**. Keep them locally. What
belongs in the repo is what the tests need: the published *numbers* (in
`reference_values.csv`, each with its table / figure / page), digitised curves
(`core/data/`), and our own description of the setup in the README.

`.gitignore` blocks the usual culprits under `core/test/papers/` (any file in an
`assets/` folder, `*.pdf`, `*.docx`, `*.doc`, `*.xlsx`, extensionless files,
figure folders). If a paper's text is needed to check a parameter, read it
locally — `extract_pdfs.py` converts a local `assets/*.pdf` to markdown with
Docling for that purpose; its output stays in the git-ignored `assets/` folder.

> This repository has never contained paper PDFs, supplementary files or full-text
> extractions; keep it that way.

## Adding a new paper

1. `mkdir papers/<author_year>` and add an entry to `papers.json` (Crossref-check
   the DOI). Put the PDF and SI in its local `assets/` and read them
   ([step 1](#1-source)).
2. Add `<author_year>.cpp` (+ `.h` if a shared setup helper is needed) following
   an existing folder as a template; build the deck from the paper's tables,
   filling the provenance table as you go ([steps 2–3](#2-provenance-audit)).
3. Transcribe the published comparison values into `reference_values.csv`;
   put experimental curves under `core/data/<author_year>/` with a `README.md`
   recording each file's provenance ([step 4](#4-data)).
4. Record every comparison with a `validation::CaseReport`; run once, then set
   each tolerance from the measured error ([step 5](#5-metrics)).
5. Re-run CMake (the test glob picks up new `.cpp` files), build, `ctest`, and
   regenerate the docs page.
6. Fill in `README.md`: citation, reproduced cases, provenance table,
   model-form audit, attribution runs and measured errors
   ([steps 6–7](#6-attribution)).

## Note on material / isotherm tests

Pure **serialization / isotherm-function** tests that happen to use a paper's
published parameters live in `core/test/unit/test_adsorbents.cpp` — they test
the serialization machinery, not a simulation, so they stay as unit tests. The
relevant paper README cross-references them.

## Gas heat capacity per paper

The engine's default gas heat capacity is the ideal-gas C_p0(T) of each
component at the local temperature (NASA TM-4513 polynomials,
`docs/equation-of-state.md`). A paper deck follows its paper instead:

- a paper that uses a constant gas heat capacity and prints it gets the
  `constant` mode (`fluid.gas_heat_capacity_mode`) with the printed value per
  component (`fluid.gas_heat_capacity[<component>]`, J/(mol K); per-kilogram
  values are converted with the component's molar mass);
- a paper that says the value is constant but prints none gets the constant
  mode with C_p0 evaluated at the paper's stated temperature or averaged over
  its stated window;
- a paper that says nothing keeps the default.

| Paper | What the paper assumes | Deck setting | Effect on the metrics (2026-11) |
|---|---|---|---|
| Arora & Hasan 2021 | gas C_p appears only as a symbol, no value | default `ideal_gas` | none (unchanged) |
| Casas et al. 2012 | constant per species, averaged over 298–373 K at 15 bar, mixed by concentration; no value printed | `constant`: mean of C_p0(T) over 298–373 K, (H0(373 K) − H0(298 K)) / 75 K, to 4 significant figures: CO₂ 38.73, H₂ 29.06, He 20.79 J/(mol K) | breakthrough-time errors move by ≤ 0.02 pp (Fig. 4 CO₂ 2.137 → 2.146 %, bound 2.4 %) |
| Ferré et al. 2026 | constant 30.70 J/(mol K) for the gas and the adsorbed phase | `constant` 30.7 J/(mol K), CO₂ and CH₄ (adsorbed phase was already 30.7) | largest: MIL-160(Al) CO₂ purity 2.152 → 2.186 pp (bound 2.2), NaY(F) CO₂ purity 0.392 → 0.335 pp, NaY(W) CO₂ purity 1.265 → 1.214 pp, 13X CH₄ purity 0.005 → 0.014 pp (bound 0.076); others ≤ 0.1 pp or 0.4 % |
| Nguyen et al. 2022 | paper text not available here | default `ideal_gas` | none (unchanged) |
| Nguyen et al. 2023 | constant per species from NIST: N₂ 1040, CO₂ 841 J/(kg K) (H₂O 4165, not in these decks); adsorbed phase = gas | `constant`: CO₂ 841 × 0.0440098 = 37.01, N₂ 1040 × 0.02801348 = 29.13 J/(mol K); adsorbed mode `gas_phase` (was constant 37.435 / 29.171) | every metric moves by < 0.01 pp or 0.03 % (largest: Exp. I recovery vs paper model 0.688 → 0.695 pp, bound 1.1) |
| Estupiñan Perez et al. 2019 | constant 30.7 J/(mol K) for the gas ("standard"), adsorbed phase assumed equal | `constant` 30.7 J/(mol K), CO₂ and N₂ | ≤ 0.005 pp (largest: Config A recovery vs paper model 0.351 → 0.356 pp, bound 0.43; Config B purity vs experiment 0.691 → 0.693 pp, bound 0.7) |
| Streb & Mazzotti 2021 | constant, average gas heat capacity; no value printed | `constant`: C_p0 at the feed temperature (set earlier; verified, not changed) | none in this change |
| Vysyaraju et al. 2024 | constant per species from NIST: N₂ 1040, He 5193 J/(kg K); adsorbed phase = gas | `constant`: N₂ 1040 × 0.02801348 = 29.13, He 5193 × 0.004002602 = 20.79 J/(mol K) (adsorbed N₂ stays 29.1) | every metric moves by < 0.001 in its unit |
| Ward & Pini 2022 | one constant 840 J/(kg K) for the CO₂/He gas and the adsorbed phase | `constant`: 840 J/(kg K) per component with its own molar mass, CO₂ 36.97, He 3.362 J/(mol K); the engine's mixture C_p (Σ y_i C_p,i / Σ y_i M_i) is then exactly 840 J/(kg K) for every composition | ≤ 0.1 pp: Q200 front-time error 1.724 → 1.735 % (bound 1.8 %), Q200 outlet NRMSE 0.1402 → 0.1410 (0.17), Q200 thermocouple RMSE 2.58 → 2.56 K |
| Ward & Pini 2024 | constant 30.7 J/(mol K) for the gas and the adsorbed phase (Table 2) | `constant` 30.7 J/(mol K), CO₂ and N₂ | CO₂ purity error 2.406 → 2.391 pp (bound 2.7), productivity 2.13 → 2.12 % |

No bound was changed. Molar masses are CoolProp's (the engine's). The
Ward 2022 `MTC_Invariance` check (CO₂/air, not a paper case) keeps the default.
