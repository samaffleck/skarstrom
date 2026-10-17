# Perez et al. (2019) — CO₂/N₂ VSA on Zeolite 13X

Validation cases reproducing the vacuum swing adsorption (VSA) cycles for
CO₂/N₂ separation reported by Perez et al. (2019).

## Citation

> Estupiñan Perez, L., Sarkar, P., Rajendran, A. *Experimental validation of
> multi-objective optimization techniques for design of vacuum swing adsorption
> processes.* Separation and Purification Technology 224:553–563 (2019).
> DOI: [10.1016/j.seppur.2019.05.039](https://doi.org/10.1016/j.seppur.2019.05.039)


## Files

| File | Purpose |
|------|---------|
| `perez_2019.cpp` | GoogleTest cases — one `TEST_F` per reproduced configuration |
| `perez_2019.h`   | `Perez2019Test` fixture + shared bed/isotherm setup (`SetUpPerez2019`) |

## Reproduced configurations

| Test (`Perez2019Test.*`) | Paper reference | Notes |
|--------------------------|-----------------|-------|
| `ConfigA`, `ConfigB`, `ConfigC1` | Table 3, A / B / C1 | Basic four-step VSA cycle, CO₂/N₂ |
| `ConfigD1`, `ConfigE`, `ConfigF` | Table 3, D1 / E / F | Four-step cycle with light-product pressurisation |
| `ConfigD1_ConnectedBC`   | Configuration D1 | D1 run with connected boundary streams; asserts BC consistency between the connected-stream and standalone formulations |

## Key parameters (from `SetUpPerez2019` in `perez_2019.h`)

- Adsorbent: Zeolite 13X, competitive (extended) dual-site Langmuir
  (concentration form; non-competitive until the validation change)
- Column: L = 0.3 m, Dᵢ = 0.032 m, ε_b = 0.41, 30 cells
- Adsorbent mass ≈ 0.162 kg; ρ_solid derived from bulk density
- Heat of adsorption: CO₂ 30.85 kJ/mol, N₂ 15.93 kJ/mol
- Effective axial thermal conductivity K_z = 0.0903 W/(m·K) (Table 1; was left
  at 0 behind a TODO)
- `ConfigD1_ConnectedBC` is a 5-cycle boundary-condition consistency check, not
  a comparison with the paper

> The authoritative parameter values live in `perez_2019.h`. This table is a
> human-readable summary — if the two ever disagree, the header wins.

## Validation

The tests compare the simulation with the paper using the helpers in
`../../fixtures/validation_metrics.h`; the published numbers they use are in
[`reference_values.csv`](reference_values.csv) (each with its table / figure /
page). Every tolerance is the measured error rounded up, with the measured
value in a comment next to it. Results of the last run:
[validation page](../../../../docs/validation.md#perez-2019)
(generated from `../validation_summary.json`).

- Four-step cycles: CO2 purity within 0.1-1.2 pp of the experiments, but recovery 8-10 pp high. LPP cycles: purity 1.9-3.4 pp low and recovery 2.7-4.1 pp high, at the cyclic steady state of the per-component mass balance (E and F need about 450 and 1000 cycles to reach it).
- adsim's macropore LDF coefficient omits the c/q* capacity factor used in the paper (much faster kinetics), and its Carman-Kozeny constant is 180 against the paper's 150. The decks now use the competitive DSL (Sect. 4.1) and the paper's K_z = 0.0903 W/(m K).
