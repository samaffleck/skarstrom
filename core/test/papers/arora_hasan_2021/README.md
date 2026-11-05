# Arora & Hasan (2021) — Medical O₂ concentrator, LiLSX (PSA / PVSA)

## Citation

> Arora, A., Hasan, M.M.F. *Flexible oxygen concentrators for medical
> applications.* Sci Rep **11**, 14317 (2021).
> DOI: [10.1038/s41598-021-93796-3](https://doi.org/10.1038/s41598-021-93796-3)


## Files

| File | Purpose |
|------|---------|
| `reference_values.csv` | Table 4 purity / recovery / production / BSF |
| `arora_hasan_2021.cpp` | `Arora2021Test` fixture, `AddArora2021Adsorbents` / `SetUpArora2021BedLayer` helpers (anonymous namespace), and the two cycle cases |

## Reproduced cases

| Test (`Arora2021Test.*`) | Paper reference | Expected at CSS |
|--------------------------|-----------------|-----------------|
| `LiLSX_PSA`  | Table 4, PSA  | purity ≈ 84.72 %, recovery ≈ 30.14 % (P_O₂ ≥ 85 %) |
| `LiLSX_PVSA` | Table 4, PVSA | purity ≈ 89.97 %, recovery ≈ 65.51 % (P_O₂ ≥ 90 %) |

Both are 4-step Skarstrom cycles for air separation to O₂ on LiLSX zeolite.

## Key parameters

- Adsorbents LiX / LiLSX / 5A built from Table 2 by `AddArora2021Adsorbents`.
- Bed layer built by `SetUpArora2021BedLayer` with the Table 4 bulk density
  672.75 kg/m³ for both cases.

## Related unit tests

The Table 2 isotherm / mass-transfer parameters are also exercised as
serialization roundtrip unit tests in `core/test/unit/test_adsorbents.cpp`:
`AdsorbentTest.LiXZeolite`, `AdsorbentTest.LiLSXZeolite`, `AdsorbentTest.Zeolite5A`.

## Validation

The tests compare the simulation with the paper using the helpers in
`../../fixtures/validation_metrics.h`; the published numbers they use are in
[`reference_values.csv`](reference_values.csv) (each with its table / figure /
page). Every tolerance is the measured error rounded up, with the measured
value in a comment next to it. Results of the last run:
[validation page](../../../../docs/validation.md#arora-hasan-2021)
(generated from `../validation_summary.json`).

- O2 purity is within 3.9 (PVSA) and 5.4 (PSA) percentage points of Table 4, but recovery is far below the paper (7 % vs 30 % for PSA, 19 % vs 66 % for PVSA). This is an open disagreement, not a tuned result.
- The paper's recovery is a model value from an optimiser; the t1/t4 rows of the PSA column of Table 4 were lost in the text extraction (t1 = 3.0 s is an estimate). The isotherm form is the non-competitive dual-site Langmuir; the paper does not state whether its DSL is competitive, so it was left as is.
- The PSA deck now vents the depressurisation step counter-currently (feed end), as the paper states for every Table 4 case; it previously vented co-currently.

### Grid and advection scheme (2026-11-04)

The decks ran first-order upwind advection on 30 cells. A grid study (10-160
cells, Richardson extrapolation; [Numerics](../../../../docs/numerics.md#advection-schemes))
showed that the PSA deck's small net O2 recovery is very grid-sensitive: on
30 upwind cells it was 4.5 % against a grid limit of about 7.5 %, and 0.1 pp
needs about 550 upwind or 400 WENO cells (purity: 115 WENO cells). The PVSA
recovery needs 86 WENO cells. The decks now run WENO, the engine's default,
on 160 (PSA) and 100 (PVSA) cells; the PSA recovery still carries about 0.3 pp
of grid error there, accepted against the cost of 400 cells.

Before: upwind, 30 cells; after: WENO, 160 / 100 cells; both at the fixed 50
cycles (bold: changed bound). The PSA purity rose by 1.8 pp, further above the
paper's optimised model value, past its bound (4.5 → 6.5 pp), as the sharper
N2 front leaves less N2 in the product. The PSA recovery rose by 2.7 pp
and the PVSA recovery fell by 1.8 pp; the gap to the paper stays the open
disagreement noted above.

| Case | Quantity | Reference | Before | After | Error before → after | Bound before → after |
|---|---|---|---|---|---|---|
| LiLSX_PSA | O2 purity at cycle 50 | 84.72 | 88.25 | 90.08 | 3.53 → 5.36 | **4.5 → 6.5** |
| LiLSX_PSA | O2 recovery at cycle 50 | 30.14 | 4.54 | 7.25 | 25.60 → 22.89 | 30 |
| LiLSX_PVSA | O2 purity at cycle 50 | 89.97 | 93.96 | 93.84 | 3.99 → 3.87 | 4.7 |
| LiLSX_PVSA | O2 recovery at cycle 50 | 65.51 | 20.73 | 18.95 | 44.78 → 46.56 | 53 |
