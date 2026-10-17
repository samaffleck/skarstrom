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
unit tests in `adsim/test/unit/test_adsorbents.cpp`:
`AdsorbentTest.MgMOF74`, `AdsorbentTest.ActivatedCarbon`.

## Validation

The tests compare the simulation with the paper using the helpers in
`../../fixtures/validation_metrics.h`; the published numbers they use are in
[`reference_values.csv`](reference_values.csv) (each with its table / figure /
page). Every tolerance is the measured error rounded up, with the measured
value in a comment next to it. Results of the last run:
[validation page](../../../../docs/validation.md#ward-pini-2024)
(generated from `../validation_summary.json`).

- CO2 purity within 2.3 pp, recovery within 1.8 pp and productivity within 0.7 % of Table 4.
- The deck now uses the extended (competitive) DSL the paper states (Sect. 2.1). Energy is not compared (no pump model in the deck).
