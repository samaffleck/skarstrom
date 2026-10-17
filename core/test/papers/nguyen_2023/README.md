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
- LPP cycles (VI-X): purity within 0.1-2.7 pp of the experiments, but recovery is 6-10 pp high and exceeds 100 % for IX and X. The single-bed emulation couples the light-product pressurisation to the feed step by composition, not by moles, so CO2 enters through the LPP stream that is not drawn from the raffinate; this is a limitation of the deck/connection model, recorded rather than hidden.
- The decks now use the competitive DSL form (as the paper's ref. [10]); the paper's column length (0.267 m) is inconsistent with its adsorbent mass and particle density, and the deck keeps the mass (L = 0.304 m).
