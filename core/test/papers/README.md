# Paper validation tests

Each subfolder reproduces the experiments / cycle simulations of **one academic
publication** and compares the result with the paper **quantitatively**: every
case computes an explicit error metric against numbers transcribed from the
paper or digitised from its figures, asserts it against a tolerance, and writes
it to a machine-readable summary from which the docs page
[`docs/validation.md`](../../../docs/validation.md) is generated.

## Folder convention

```
papers/
├── papers.json               # citation / DOI / system for every folder (Crossref-checked)
├── validation_summary.json   # GENERATED: merged metrics of the last validation run
└── <author_year>/
    ├── <author_year>.cpp     # GoogleTest cases — one TEST_F per reproduced case
    ├── <author_year>.h       # (optional) paper fixture + shared setup helper
    ├── reference_values.csv  # published numbers used by the tests (with table/figure/page)
    └── README.md             # citation, DOI, reproduced cases, parameter provenance,
                              # known discrepancies and the measured errors
```

Digitised experimental curves live outside this tree, in `core/data/<paper>/`.

- **Folder name**: `surname_year` (e.g. `perez_2019`, `arora_hasan_2021`). For
  two papers by the same first author in one year, add a discriminator.
- **Fixture**: one GoogleTest suite per paper, named `<Paper>Test`
  (e.g. `Perez2019Test`), deriving from `UserDataTest`.
- **Case name**: the figure / table / configuration reproduced
  (e.g. `ConfigD1`, `Fig5_TernaryBreakthrough`, `Q100ccm`).
- **`SaveTest(...)` / `SaveResultsAsExample(...)` labels are stable** so the
  on-disk result artifacts and the bundled example decks (`examples/decks/`) keep their names.

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
ctest --preset linux-release -R "(Arora2021|Casas2012|Ferre2026|Nguyen2022|Nguyen2023|Perez2019|Streb2021|Vysyaraju2024|Ward2022|WardPini2024)" -j2
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
   the DOI).
2. Add `<author_year>.cpp` (+ `.h` if a shared setup helper is needed) following
   an existing folder as a template; build the deck from the paper's tables.
3. Transcribe the published comparison values into `reference_values.csv`;
   put digitised curves under `core/data/<author_year>/`.
4. Record every comparison with a `validation::CaseReport`; run once, then set
   each tolerance from the measured error.
5. Re-run CMake (the test glob picks up new `.cpp` files), build, `ctest`, and
   regenerate the docs page.
6. Fill in `README.md`: citation, reproduced cases, parameter provenance,
   assumptions, measured errors.

## Note on material / isotherm tests

Pure **serialization / isotherm-function** tests that happen to use a paper's
published parameters live in `core/test/unit/test_adsorbents.cpp` — they test
the serialization machinery, not a simulation, so they stay as unit tests. The
relevant paper README cross-references them.
