# Ward & Pini (2024) — CO₂/N₂ PSA on Mg-MOF-74 / activated carbon

## Citation

> Ward & Pini (2024), *Ind. Eng. Chem. Res.*
> DOI: [10.1021/acs.iecr.4c01345](https://doi.org/10.1021/acs.iecr.4c01345)


## Files

| File | Purpose |
|------|---------|
| `ward_pini_2024.cpp` | `WardPini2024Test` fixture + reproduced cycle |

## Reproduced cases

| Test (`WardPini2024Test.*`) | Notes |
|-----------------------------|-------|
| `FourStepPSA` | Four-step PSA for CO₂/N₂ separation, Mg-MOF-74 + activated carbon layers |

## Key parameters

- Adsorbent: Mg-MOF-74 (single-adsorbent case of Table 4), extended
  (competitive) dual-site Langmuir (Table 1; non-competitive until the
  validation change). The activated-carbon material is defined but not used.
- Setup is inline in the test body.

## Related unit tests

The Table 1 isotherm parameters are also exercised as serialization roundtrip
unit tests in `core/test/unit/test_adsorbents.cpp`:
`AdsorbentTest.MgMOF74`, `AdsorbentTest.ActivatedCarbon`.

## Validation

The tests compare the simulation with the paper using the helpers in
`../../fixtures/validation_metrics.h`; the published numbers they use are in
[`reference_values.csv`](reference_values.csv) (each with its table / figure /
page). Every tolerance is the measured error rounded up, with the measured
value in a comment next to it. Results of the last run:
[validation page](../../../../docs/validation.md#ward-pini-2024)
(generated from `../validation_summary.json`).

- CO2 purity within 2.5 pp, recovery within 2.9 pp and productivity within 2.3 % of Table 4 (WENO, 100 cells; see below).
- The deck now uses the extended (competitive) DSL the paper states (Sect. 2.1). Energy is not compared (no pump model in the deck).

### Gas heat capacity as in the paper (2026-11)

Table 2 gives a constant 30.7 J/(mol K) for the gas and the adsorbed phase. The deck now uses the constant gas heat capacity mode with 30.7 J/(mol K) for CO₂ and N₂ (before: the default ideal-gas C_p0 at the local temperature). CO₂ purity error 2.406 → 2.391 pp (bound 2.7), productivity error 2.13 → 2.12 %. No bound changed. All papers' settings: [../README.md](../README.md#gas-heat-capacity-per-paper).

### Grid and advection scheme (2026-11-04)

The deck ran first-order upwind advection on 30 cells. A grid study (10-160
cells, Richardson extrapolation; [Numerics](../../../../docs/numerics.md#advection-schemes))
put its CO2 recovery 0.9 pp below the grid limit there; 0.1 pp needs about
360 upwind or 94 WENO cells. The deck now runs WENO, the engine's default, on
100 cells (CSS after 155 cycles).

Before: upwind, 30 cells; after: WENO, 100 cells, both at CSS (bold: changed
bound). Recovery rose by 0.86 pp and productivity by 1.3 %, both further above
the paper's model values, past their bounds (recovery 2.1 → 3.4 pp,
productivity 0.018 → 0.027): part of the earlier agreement was the upwind
scheme's numerical dispersion.

| Case | Quantity | Reference | Before | After | Error before → after | Bound before → after |
|---|---|---|---|---|---|---|
| FourStepPSA | CO2 purity | 98.50 | 96.21 | 96.08 | 2.29 → 2.42 | 2.7 |
| FourStepPSA | CO2 recovery | 89.00 | 90.97 | 91.83 | 1.97 → 2.83 | **2.1 → 3.4** |
| FourStepPSA | CO2 productivity | 0.828 | 0.8358 | 0.8465 | 0.009377 → 0.02229 | **0.018 → 0.027** |
