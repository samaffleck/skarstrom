# Nguyen et al. (2022) — CO₂ breakthrough on CALF-20

## Citation

> Nguyen, T.T.T., Lin, J.-B., Shimizu, G.K.H., Rajendran, A. *Separation of CO₂
> and N₂ on a hydrophobic metal organic framework CALF-20.* Chemical Engineering
> Journal 442:136263 (2022).
> DOI: [10.1016/j.cej.2022.136263](https://doi.org/10.1016/j.cej.2022.136263)

This is a different publication from `../nguyen_2023` (the VSA-cycle paper,
DOI 10.1016/j.cej.2022.139550), which cites it as its ref. [10] for the CALF-20
equilibrium and breakthrough characterisation. The citation was confirmed from
that reference list and Crossref; the paper text itself is not in the repo.

## Files

| File | Purpose |
|------|---------|
| `nguyen_2022.cpp` | GoogleTest cases (the two bundled examples + four further feed fractions) |
| `nguyen_2022.h`   | `Nguyen2022Test` fixture + `SetUpNguyen2022()` |

## Reproduced cases

| Test (`Nguyen2022Test.*`) | Paper reference | Notes |
|---------------------------|-----------------|-------|
| `Fig4_BinaryBreakthrough`  | Fig. 4 | Binary CO₂/He breakthrough, CALF-20, inlet y_CO₂ = 0.7 |
| `Fig5_TernaryBreakthrough` | Fig. 5 | Ternary CO₂/N₂/He breakthrough, CALF-20, inlet y_CO₂ = 0.15 |

Experimental data from the paper is bundled with the saved examples (Data
Plotter window). Varying the inlet CO₂ fraction reproduces the other cases in
the figures.

## Key parameters

- Adsorbent: CALF-20; setup in `SetUpNguyen2022()` (see `nguyen_2022.h`).
- Not stated in the material available here (the paper text was never added to
  the repo): particle diameter (1.5 mm assumed) and ambient temperature (295 K,
  the initial thermocouple reading of the digitised data). See the comments in
  `nguyen_2022.h`.
- Extra cases `Fig4_Binary_50pct`, `Fig4_Binary_100pct`, `Fig5_Ternary_5pct`,
  `Fig5_Ternary_30pct` use the digitised curves for the other feed fractions,
  assuming the same flow as the bundled case of that figure.

## Validation

The tests compare the simulation with the paper using the helpers in
`../../fixtures/validation_metrics.h`; the published numbers they use are in
[`reference_values.csv`](reference_values.csv) (each with its table / figure /
page). Every tolerance is the measured error rounded up, with the measured
value in a comment next to it. Results of the last run:
[validation page](../../../../docs/validation.md#nguyen-2022)
(generated from `../validation_summary.json`).

- Outlet CO2 curves agree with the digitised Figs. 4 and 5 to an NRMSE of 8-12 % and the CO2 front position (half the feed fraction) to 0.7-5 %.
- The paper text is not in the repo. The four additional feed fractions assume the same flow as the bundled case of the same figure (50 ccm for Fig. 4, 100 ccm for Fig. 5); Fig. 5 is described as a CO2/N2/He experiment but the deck feeds CO2 in He. Particle diameter and ambient temperature are assumptions (see the README).
