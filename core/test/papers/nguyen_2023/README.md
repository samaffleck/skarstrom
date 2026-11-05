# Nguyen et al. (2023) — CO₂/N₂ VPSA cycles on CALF-20

## Citation

> Nguyen et al., *Chemical Engineering Journal* (2022/2023).
> DOI: [10.1016/j.cej.2022.139550](https://doi.org/10.1016/j.cej.2022.139550)


## Files

| File | Purpose |
|------|---------|
| `nguyen_2023.cpp` | GoogleTest cases + `RunNguyen2023Experiment` / `RunNguyen2023LPP` helpers |
| `nguyen_2023.h`   | `Nguyen2023Test` fixture + `SetUpNguyen2023()` |

## Reproduced cases

Ten experiments from the paper. Experiments I–V are the basic four-step cycle
(`RunNguyen2023Experiment`); VI–X add light-product pressurisation
(`RunNguyen2023LPP`). Each maps to a saved example `Nguyen2023_Exp_<n>`.

| Test (`Nguyen2023Test.*`) | Cycle | feed_time [s] | feed vel [m/s] | bd P [bar] | evac P [bar] |
|---------------------------|-------|---------------|----------------|------------|--------------|
| `ExperimentI`    | 4-step     | 176 | 0.093 | 0.119 | — |
| `ExperimentII`   | 4-step     | 174 | 0.065 | 0.148 | — |
| `ExperimentIII`  | 4-step     | 174 | 0.058 | 0.189 | — |
| `ExperimentIV`   | 4-step     | 174 | 0.055 | 0.233 | — |
| `ExperimentV`    | 4-step     | 174 | 0.046 | 0.263 | — |
| `ExperimentVI`   | 4-step LPP | 193 | 0.041 | 0.063 | 0.031 |
| `ExperimentVII`  | 4-step LPP | 157 | 0.024 | 0.046 | 0.031 |
| `ExperimentVIII` | 4-step LPP | 162 | 0.023 | 0.060 | 0.034 |
| `ExperimentIX`   | 4-step LPP | 156 | 0.021 | 0.066 | 0.031 |
| `ExperimentX`    | 4-step LPP | 156 | 0.022 | 0.080 | 0.031 |

> Values above are a human-readable summary; the authoritative parameters live
> in the `TEST_F` bodies in `nguyen_2023.cpp`.

## Key parameters

- Adsorbent: CALF-20; setup in `SetUpNguyen2023()` (see `nguyen_2023.h`).
- Competitive (extended) DSL with the SI Table S1 parameters (N₂ equal-energy
  scheme); the deck used the non-competitive form until the validation change.
- Table 1 is internally inconsistent (36 g at ρ_p = 570 kg/m³ and ε = 0.4 does
  not fit a 0.267 m column); the deck keeps the mass and derives L ≈ 0.304 m.
- Ambient and feed temperature are not stated in the paper or SI text; 293 K is
  assumed (documented in the code).
- Validation KPIs: purity of the evacuation product; recovery = CO₂ in that
  product / CO₂ fed in the feed (and, for the 4-step cycle, feed-pressurisation)
  steps; productivity per m³ of adsorbent (m_ads/ρ_p) per cycle time.

## Validation

The tests compare the simulation with the paper using the helpers in
`../../fixtures/validation_metrics.h`; the published numbers they use are in
[`reference_values.csv`](reference_values.csv) (each with its table / figure /
page). Every tolerance is the measured error rounded up, with the measured
value in a comment next to it. Results of the last run:
[validation page](../../../../docs/validation.md#nguyen-2023)
(generated from `../validation_summary.json`).

- Four-step cycles (I-V): CO2 purity within 1.7-4.9 pp and recovery within 1.2-2.7 pp of the experiments; productivity within 5-9 %.
- LPP cycles (VI-X): purity within 0.9-2.8 pp of the experiments, but recovery is 5-8 pp high, at the cyclic steady state of the per-component mass balance (VIII-X need about 140-440 cycles to reach it). Before the Danckwerts inlet and Darcy flow-face fixes (2026-11) recovery exceeded 100 % for IX and X (100.35 % and 100.78 %): the inlet faces let dispersion carry CO₂ into the bed that no stream log recorded, so more CO₂ left than was counted in. With the fixes they are 98.3 % and 97.6 %; purity moved by up to 1.2 pp (X: 80.84 → 79.67 %, now 0.93 pp from the experiment) and VIII's productivity by 1.4 %, and those tolerances were re-measured. The single-bed emulation couples the light-product pressurisation to the feed step by composition, not by moles, so CO2 enters through the LPP stream that is not drawn from the raffinate; this is a limitation of the deck/connection model, recorded rather than hidden.
- The decks now use the competitive DSL form (as the paper's ref. [10]); the paper's column length (0.267 m) is inconsistent with its adsorbent mass and particle density, and the deck keeps the mass (L = 0.304 m).

### Grid and advection scheme (2026-11-04)

The decks ran first-order upwind advection on 100 cells. A grid study (10-160
cells, Richardson extrapolation; [Numerics](../../../../docs/numerics.md#advection-schemes))
put the CO2 recovery 0.2 (I), 0.3 (V) and 0.35 pp (VI) below its grid limit
there; 0.1 pp needs 170-340 upwind cells but 40-76 WENO cells (I 68, V 76,
VI 40; II-IV lie between I and V, VII-X share VI's column). The decks now run
WENO, the engine's default, on 80 cells.

The decks also start from a feed-loaded bed (owner's decision, 2026-11-04):
the feed gas (15 % CO2) with the loadings in equilibrium with it, instead of a
clean N2 bed. The cyclic steady state is the same, reached in fewer cycles.

Before: upwind, 100 cells, clean start; after: WENO, 80 cells, feed-loaded
start (bold: changed bound). Every recovery rose by 0.15-0.40 pp and the
purities moved by less than 0.2 pp; no bound was crossed. Cycles to CSS from
the feed-loaded start: I-V 37, 42, 44, 46, 50; VI-X 80, 326, 447, 1166 and
1558 (from a clean bed with WENO: VI-X 104, 391, 532, 1283, 1669; under
upwind from a clean bed X had not reached CSS by the old cap of 1000, so its
earlier numbers were not at CSS: its purity error fell from 0.12 to 0.04 pp
and that bound was tightened, as was VI's productivity bound). The caps are
about 30 % above the feed-start stops (VI 110, VII 430, VIII 590, IX 1520,
X 2030; I-V keep 250). IX and X stay slow from any start: X's slowest mode
contracts by only about 0.3 % per cycle (0.997), because the deck couples the
light-product pressurisation to the feed step by composition, not by moles
(the limitation noted above), which leaves the CO2 inventory weakly
constrained; its cap stays high for that reason.

| Case | Quantity | Reference | Before | After | Error before → after | Bound before → after |
|---|---|---|---|---|---|---|
| ExperimentI | CO2 purity vs experiment | 97.70 | 95.70 | 95.71 | 2.00 → 1.99 | 2.4 |
| ExperimentI | CO2 recovery vs experiment | 66.90 | 69.27 | 69.41 | 2.37 → 2.51 | 2.8 |
| ExperimentI | CO2 productivity vs experiment | 0.65 | 0.6068 | 0.6081 | 0.06653 → 0.06451 | 0.076 |
| ExperimentI | CO2 purity vs paper model | 96.70 | 95.70 | 95.71 | 1.00 → 0.99 | 1.4 |
| ExperimentI | CO2 recovery vs paper model | 70.10 | 69.27 | 69.41 | 0.83 → 0.69 | 1.1 |
| ExperimentI | CO2 productivity vs paper model | 0.63 | 0.6068 | 0.6081 | 0.03689 → 0.03481 | 0.046 |
| ExperimentII | CO2 purity vs experiment | 91.20 | 89.33 | 89.36 | 1.87 → 1.84 | 2.2 |
| ExperimentII | CO2 recovery vs experiment | 71.70 | 72.81 | 73.01 | 1.11 → 1.31 | 1.6 |
| ExperimentII | CO2 productivity vs experiment | 0.52 | 0.4702 | 0.4716 | 0.09581 → 0.09313 | 0.11 |
| ExperimentII | CO2 purity vs paper model | 90.30 | 89.33 | 89.36 | 0.97 → 0.94 | 1.3 |
| ExperimentII | CO2 recovery vs paper model | 72.30 | 72.81 | 73.01 | 0.51 → 0.71 | 0.94 |
| ExperimentII | CO2 productivity vs paper model | 0.54 | 0.4702 | 0.4716 | 0.1293 → 0.1267 | 0.15 |
| ExperimentIII | CO2 purity vs experiment | 85.60 | 83.88 | 83.92 | 1.72 → 1.68 | 2.1 |
| ExperimentIII | CO2 recovery vs experiment | 71.80 | 74.14 | 74.37 | 2.34 → 2.57 | 2.9 |
| ExperimentIII | CO2 productivity vs experiment | 0.48 | 0.4383 | 0.4397 | 0.0869 → 0.084 | 0.098 |
| ExperimentIII | CO2 purity vs paper model | 85.70 | 83.88 | 83.92 | 1.82 → 1.78 | 2.2 |
| ExperimentIII | CO2 recovery vs paper model | 73.40 | 74.14 | 74.37 | 0.74 → 0.97 | 1.2 |
| ExperimentIII | CO2 productivity vs paper model | 0.51 | 0.4383 | 0.4397 | 0.1406 → 0.1379 | 0.16 |
| ExperimentIV | CO2 purity vs experiment | 82.30 | 79.10 | 79.15 | 3.20 → 3.15 | 3.7 |
| ExperimentIV | CO2 recovery vs experiment | 72.50 | 75.00 | 75.23 | 2.50 → 2.73 | 3.1 |
| ExperimentIV | CO2 productivity vs experiment | 0.46 | 0.4266 | 0.428 | 0.07254 → 0.06948 | 0.081 |
| ExperimentIV | CO2 purity vs paper model | 81.30 | 79.10 | 79.15 | 2.20 → 2.15 | 2.6 |
| ExperimentIV | CO2 recovery vs paper model | 74.30 | 75.00 | 75.23 | 0.70 → 0.93 | 1.2 |
| ExperimentIV | CO2 productivity vs paper model | 0.5 | 0.4266 | 0.428 | 0.1467 → 0.1439 | 0.17 |
| ExperimentIX | CO2 purity vs experiment | 85.67 | 84.50 | 84.40 | 1.17 → 1.27 | 1.6 |
| ExperimentIX | CO2 recovery vs experiment | 90.50 | 98.58 | 98.89 | 8.08 → 8.39 | 12 |
| ExperimentIX | CO2 productivity vs experiment | 0.13 | 0.1426 | 0.1431 | 0.09697 → 0.1005 | 0.14 |
| ExperimentIX | CO2 purity vs paper model | 88.30 | 84.50 | 84.40 | 3.80 → 3.90 | 4.2 |
| ExperimentIX | CO2 recovery vs paper model | 91.60 | 98.58 | 98.89 | 6.98 → 7.29 | 11 |
| ExperimentIX | CO2 productivity vs paper model | 0.14 | 0.1426 | 0.1431 | 0.01862 → 0.02187 | 0.047 |
| ExperimentV | CO2 purity vs experiment | 78.40 | 73.49 | 73.55 | 4.91 → 4.85 | 5.7 |
| ExperimentV | CO2 recovery vs experiment | 74.20 | 76.58 | 76.84 | 2.38 → 2.64 | 3 |
| ExperimentV | CO2 productivity vs experiment | 0.4 | 0.3791 | 0.3805 | 0.05215 → 0.04864 | 0.06 |
| ExperimentV | CO2 purity vs paper model | 75.60 | 73.49 | 73.55 | 2.11 → 2.05 | 2.4 |
| ExperimentV | CO2 recovery vs paper model | 75.40 | 76.58 | 76.84 | 1.18 → 1.44 | 1.7 |
| ExperimentV | CO2 productivity vs paper model | 0.44 | 0.3791 | 0.3805 | 0.1383 → 0.1351 | 0.16 |
| ExperimentVI | CO2 purity vs experiment | 97.20 | 96.20 | 96.21 | 1.00 → 0.99 | 1.3 |
| ExperimentVI | CO2 recovery vs experiment | 83.10 | 89.48 | 89.83 | 6.38 → 6.73 | 7.8 |
| ExperimentVI | CO2 productivity vs experiment | 0.3 | 0.279 | 0.2802 | 0.06987 → 0.06614 | 0.076 |
| ExperimentVI | CO2 purity vs paper model | 96.70 | 96.20 | 96.21 | 0.50 → 0.49 | 0.74 |
| ExperimentVI | CO2 recovery vs paper model | 84.80 | 89.48 | 89.83 | 4.68 → 5.03 | 5.9 |
| ExperimentVI | CO2 productivity vs paper model | 0.28 | 0.279 | 0.2802 | 0.003429 → 0.0005595 | **0.011 → 0.0026** |
| ExperimentVII | CO2 purity vs experiment | 94.30 | 93.45 | 93.46 | 0.85 → 0.84 | 1.1 |
| ExperimentVII | CO2 recovery vs experiment | 88.40 | 94.84 | 95.25 | 6.44 → 6.85 | 9 |
| ExperimentVII | CO2 productivity vs experiment | 0.14 | 0.1573 | 0.158 | 0.1236 → 0.1284 | 0.16 |
| ExperimentVII | CO2 purity vs paper model | 95.50 | 93.45 | 93.46 | 2.05 → 2.04 | 2.3 |
| ExperimentVII | CO2 recovery vs paper model | 88.10 | 94.84 | 95.25 | 6.74 → 7.15 | 9.3 |
| ExperimentVII | CO2 productivity vs paper model | 0.14 | 0.1573 | 0.158 | 0.1236 → 0.1284 | 0.16 |
| ExperimentVII | CO2 purity vs repeat experiment VII.b | 93.80 | 93.45 | 93.46 | 0.35 → 0.34 | – |
| ExperimentVII | CO2 recovery vs repeat experiment VII.b | 87.50 | 94.84 | 95.25 | 7.34 → 7.75 | – |
| ExperimentVIII | CO2 purity vs experiment | 92.70 | 89.91 | 89.90 | 2.79 → 2.80 | 3.1 |
| ExperimentVIII | CO2 recovery vs experiment | 87.80 | 95.75 | 96.14 | 7.95 → 8.34 | 11 |
| ExperimentVIII | CO2 productivity vs experiment | 0.17 | 0.1545 | 0.1552 | 0.09109 → 0.08733 | 0.11 |
| ExperimentVIII | CO2 purity vs paper model | 93.60 | 89.91 | 89.90 | 3.69 → 3.70 | 4.1 |
| ExperimentVIII | CO2 recovery vs paper model | 89.90 | 95.75 | 96.14 | 5.85 → 6.24 | 8.3 |
| ExperimentVIII | CO2 productivity vs paper model | 0.16 | 0.1545 | 0.1552 | 0.03428 → 0.03029 | 0.04 |
| ExperimentX | CO2 purity vs experiment | 80.60 | 80.72 | 80.56 | 0.12 → 0.04 | **1.1 → 0.13** |
| ExperimentX | CO2 recovery vs experiment | 92.40 | 99.15 | 99.43 | 6.75 → 7.03 | 9.7 |
| ExperimentX | CO2 productivity vs experiment | 0.15 | 0.1503 | 0.1507 | 0.001794 → 0.004638 | 0.029 |
| ExperimentX | CO2 purity vs paper model | 80.60 | 80.72 | 80.56 | 0.12 → 0.04 | **1.1 → 0.13** |
| ExperimentX | CO2 recovery vs paper model | 93.60 | 99.15 | 99.43 | 5.55 → 5.83 | 8.3 |
| ExperimentX | CO2 productivity vs paper model | 0.13 | 0.1503 | 0.1507 | 0.1559 → 0.1592 | 0.21 |
