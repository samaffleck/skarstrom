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
`../../fixtures/validation_metrics.h`, against the digitised outlet curves in
[`core/data/nguyen_2022/`](../../../data/nguyen_2022/) (Figs. 4 and 5). No
number printed in the paper is compared, so this suite has no
`reference_values.csv`. Every tolerance is the measured error rounded up, with the measured
value in a comment next to it. Results of the last run:
[validation page](../../../../docs/validation.md#nguyen-2022)
(generated from `../validation_summary.json`).

- Outlet CO2 curves agree with the digitised Figs. 4 and 5 to an NRMSE of 8-12 % and the CO2 front position (half the feed fraction) to 0.7-5 %.
- The paper text is not in the repo. The four additional feed fractions assume the same flow as the bundled case of the same figure (50 ccm for Fig. 4, 100 ccm for Fig. 5); Fig. 5 is described as a CO2/N2/He experiment but the deck feeds CO2 in He. Particle diameter and ambient temperature are assumptions (see the README).

### Grid and advection scheme (2026-11-04)

The decks keep their 30 cells and now run WENO, the engine's default, instead
of first-order upwind. A grid study (10-160 cells, Richardson extrapolation;
[Numerics](../../../../docs/numerics.md#advection-schemes)) found the CO2
half-feed time within 0.5 % of its grid limit on 10 (Fig. 4) and 21 (Fig. 5)
WENO cells (11 and 26 upwind), so 30 cells are converged. The outlet curves
are sharper: every NRMSE rose by 0.001-0.017, the front times moved by less
than 0.3 % (19 s of 7400 for the 5 % feed, whose error fell from 0.31 % to
0.06 %; that bound was tightened). Before: upwind; after: WENO; both on 30
cells (bold: changed bound).

| Case | Quantity | Reference | Before | After | Error before → after | Bound before → after |
|---|---|---|---|---|---|---|
| Fig4_BinaryBreakthrough | outlet y_CO2 | 0.7 | 0.05689 | 0.06231 | 0.08128 → 0.08901 | 0.1 |
| Fig4_BinaryBreakthrough | CO2 front time (y_CO2 = y_feed/2) | 1747.30 | 1730.66 | 1733.98 | 0.009525 → 0.007626 | 0.027 |
| Fig4_Binary_100pct | outlet y_CO2 | 1.00 | 0.1021 | 0.1053 | 0.1021 → 0.1053 | 0.12 |
| Fig4_Binary_100pct | CO2 front time (y_CO2 = y_feed/2) | 1258.45 | 1248.50 | 1250.91 | 0.007906 → 0.00599 | 0.017 |
| Fig4_Binary_50pct | outlet y_CO2 | 0.5 | 0.04642 | 0.0486 | 0.09283 → 0.0972 | 0.14 |
| Fig4_Binary_50pct | CO2 front time (y_CO2 = y_feed/2) | 2441.30 | 2360.61 | 2364.22 | 0.03305 → 0.03158 | 0.06 |
| Fig5_TernaryBreakthrough | outlet y_CO2 | 0.15 | 0.01132 | 0.01275 | 0.0755 → 0.08497 | 0.097 |
| Fig5_TernaryBreakthrough | CO2 front time (y_CO2 = y_feed/2) | 3364.00 | 3206.85 | 3202.04 | 0.04672 → 0.04814 | 0.059 |
| Fig5_Ternary_30pct | outlet y_CO2 | 0.3 | 0.02911 | 0.03033 | 0.09704 → 0.1011 | 0.13 |
| Fig5_Ternary_30pct | CO2 front time (y_CO2 = y_feed/2) | 1708.31 | 1727.23 | 1727.87 | 0.01108 → 0.01145 | 0.039 |
| Fig5_Ternary_5pct | outlet y_CO2 | 0.05 | 0.006039 | 0.006875 | 0.1208 → 0.1375 | 0.15 |
| Fig5_Ternary_5pct | CO2 front time (y_CO2 = y_feed/2) | 7397.26 | 7420.14 | 7401.34 | 0.003093 → 0.0005517 | **0.026 → 0.0026** |
