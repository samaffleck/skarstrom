# Streb et al. (2021) — CO₂/CH₄ separation breakthrough

## Citation

> Streb, A., Mazzotti, M. *Adsorption for efficient low carbon hydrogen
> production: part 1 — adsorption equilibrium and breakthrough studies for
> H₂/CO₂/CH₄ on zeolite 13X.* Adsorption 27:541–558 (2021).
> DOI: [10.1007/s10450-021-00306-y](https://doi.org/10.1007/s10450-021-00306-y)


## Files

| File | Purpose |
|------|---------|
| `streb_2021.cpp` | `Streb2021Test` fixture + reproduced breakthrough cases |

## Reproduced cases

| Test (`Streb2021Test.*`) | Paper reference | Notes |
|--------------------------|-----------------|-------|
| `ExperimentA1` | Table 2, A1 (25 °C, 1 bar, 20 cm³/s) | Binary CO₂/CH₄ 50:50, bundled example |
| `ExperimentA2` | A2 (25 °C, 10 bar, 10 cm³/s) | Binary |
| `ExperimentA3` | A3 (25 °C, 10 bar, 20 cm³/s) | Binary, reference experiment (Fig. 6), bundled example |
| `ExperimentA4` | A4 (25 °C, 10 bar, 40 cm³/s) | Binary |
| `ExperimentA5` | A5 (25 °C, 25 bar, 20 cm³/s) | Binary |
| `ExperimentA6` | A6 (45 °C, 10 bar, 20 cm³/s) | Binary |
| `ExperimentA7` | A7 (65 °C, 10 bar, 20 cm³/s) | Binary |
| `ExperimentB1` | B1 (25 °C, 1 bar, 20 cm³/s) | Ternary H₂:CO₂:CH₄ 75:20:5 |
| `ExperimentB2` | B2 (25 °C, 10 bar, 10 cm³/s) | Ternary |

## Key parameters

- Zeolite 13X, Table 1 column; extended Sips isotherm (Table 3, H₂ column added
  for B1/B2), LDF k = 0.06 / 0.2 / 1.0 s⁻¹ (CO₂/CH₄/H₂, Table 4), constant heats
  36.3 / 18.8 / 10.7 kJ/mol (Sect. 4.1), Dixon (1996) internal HTC,
  h_w = 220 W/(m²·K). Shared set-up in `Streb2021Test::SetUpStreb`.
- Digitised data in `core/data/streb_2021/<case>/` (outlet mole fractions and
  thermocouples at z = 10, 35, 60, 85, 110 cm).

## Validation

The tests compare the simulation with the paper using the helpers in
`../../fixtures/validation_metrics.h`; the published numbers they use are in
[`reference_values.csv`](reference_values.csv) (each with its table / figure /
page). Every tolerance is the measured error rounded up, with the measured
value in a comment next to it. Results of the last run:
[validation page](../../../../docs/validation.md#streb-2021)
(generated from `../validation_summary.json`).

- CO2 front positions agree with the digitised experiments to 0.5-9 %, CH4 fronts to 0.3-15 %. Curve NRMSE is large (13-83 % for CH4) because the fronts are nearly vertical: a 10-20 s shift on a 20 s front dominates the RMSE. Thermocouple RMSE is 11-33 K, driven by the timing and height of the 30-100 K temperature fronts at the centre line.
- adsim uses extended Sips competition and the MFC set-point velocity; the paper shows IAST/RAST predict the binary cases better than extended Sips and fits the feed velocity per experiment within the MFC uncertainty. The ~10 m of sampling line to the MS is not modelled.
- A1/A3 previously compared simulated CO2 against the digitised CH4 curve (alphabetical file order) with a bound of 0.8 in mole fraction; fixed.
