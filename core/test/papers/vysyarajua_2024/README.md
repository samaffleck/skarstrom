# Vysyaraju et al. (2024) — He/N₂ PSA on binderless Zeolite 5A

## Citation

> Vysyaraju, R., Sawada, J.A., Estupiñan Perez, L., Rajendran, A. *Simulation and
> experimental demonstration of helium purification from He/N₂ mixtures by
> pressure swing adsorption.* Chemical Engineering Journal (2024).
> DOI: [10.1016/j.cej.2024.155929](https://doi.org/10.1016/j.cej.2024.155929)


## Model notes

- Helium is the **light, essentially non-adsorbing** product; N₂ is the adsorbed
  heavy impurity (opposite dynamic to the CO₂-capture papers). N₂ is fitted with a
  concentration-form **dual-site Langmuir** (Table 2); He has no isotherm/LDF.
- Mass transfer is macropore-diffusion controlled: constant LDF
  `k = 15·εp·Dm/(τ·rp²) ≈ 105 /s` (equilibrium-controlled separation).
- Energy balance uses the **ISOSTERIC** heat-of-adsorption convention — the N₂
  breakthrough excursion (~9 K) is only reproduced when the full isosteric heat
  `q_st ≈ 18.7 kJ/mol` is released; the internal-energy form (`dH − Cp,a·T`,
  used for the Ward/Pini CO₂ papers) releases only ~10 kJ/mol here and
  under-predicts the rise by ~4 K.

## Files

| File | Purpose |
|------|---------|
| `vysyarajua_2024.cpp` | GoogleTest cases |

## Reproduced cases

### Dynamic column breakthrough (Section 4.2) — validates the N₂ isotherm + energy balance

| Test (`Vysyaraju2024Test.*`) | Condition | N₂ loading (model / paper) | Peak T (model / paper) |
|------------------------------|-----------|----------------------------|------------------------|
| `DCB_Ambient_N2`      | 100 % N₂, 0.95 bar, 295 K (Fig. 4)  | 0.53 / 0.56 mol/kg | 304 / 304 K |
| `DCB_HighPressure_N2` | 100 % N₂, 11.2 bar, 298 K (Fig. 5b) | 2.12 / 2.16 mol/kg | — |

### 5-step PSA cycle (Section 3.2) — validates He purity

`PSA_Exp1` reproduces Table 3 **Exp #1** (10 % He, t_ADS 15 s, P_H 15.1 bar):
PRES → ADS (He light product, top) → **PED** → BLO → **PER**, with a top-top
pressure equalisation (PED donor → PER consumer) via `MakeMutualConnection`.

| Metric | Model | Paper |
|--------|-------|-------|
| He purity  | 99.4 % | 99.9 % |
| He recovery | 63 % | 31 % |

Reaching a clean cyclic steady state required (see the header comment in the
`.cpp`): **valve-limited** PRES to `P_H` (a fixed-velocity feed under/overshoots;
a pressure ramp floods the bed), **seeding a He bank** at ~`P_H` (the cycle is
bistable), and a **flow-driven PER** (velocity-0, driven by the connection) so
it receives what PED released rather than over-pulling. Blow-down from `P_INT`
(not `P_H`) still vents enough N₂ to keep the front low.

**Purity is 0.5 pp below the paper model** (1.4 pp above the experiment).
**Recovery reads about twice the paper's**: in this single-bed deck the bed's
own equalisation gives back to it every mole PED released, whereas the
experiment and the paper's model are a 3-bed cycle. The cause of the gap is
not yet attributed. Recovery is computed by conservation
(`1 − blow-down He / fed He`); it agrees with product He / fed He to ~0.1 pp.

Until November 2026 the connection replayed the magnitude of PED's flow, so
the backflow through the PED valve for most of the step (the flow reverses
about a second in) reached PER as more inflow: PER took in
29 % more gas than PED released, and the conservation recovery read 17 %. With
mole-conserving equalisation the numbers above replaced the earlier 99.9 %
purity and 17 % recovery, and the tolerances were re-measured.

## Validation

The tests compare the simulation with the paper using the helpers in
`../../fixtures/validation_metrics.h`; the published numbers they use are in
[`reference_values.csv`](reference_values.csv) (each with its table / figure /
page). Every tolerance is the measured error rounded up, with the measured
value in a comment next to it. Results of the last run:
[validation page](../../../../docs/validation.md#vysyarajua-2024)
(generated from `../validation_summary.json`).

- The N2 breakthrough loadings match the reported value (ambient) and the DSL value (11.2 bar) to 2-4 %, and the peak bed temperature to 0.3 K.
- He purity is 0.54 pp below the paper model (1.4 pp above the experiment), and He recovery is 63 % against 31 %: the deck is a single bed with self-equalisation, the experiment and the paper's model are a 3-bed cycle.
