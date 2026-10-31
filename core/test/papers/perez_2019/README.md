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

## Provenance (WP10 audit, 2026-10-31, from the paper text and SI)

| Deck value | Value | Source | Note |
|---|---|---|---|
| Column length, inner radius, outer radius | 0.3 m, 0.016 m, 0.0191 m | Table 1 | wall thickness 3.1 mm (derived) |
| Bed voidage, particle voidage, tortuosity | 0.41, 0.35, 3 | Table 1 ("assumed") | |
| Particle radius | 7.5 × 10⁻⁴ m | Table 1 ("measured") | |
| Adsorbent mass | 0.162 kg per column | abstract, Sect. 2 | particle density derived: m / (V (1 − ε)) = 1138 kg/m³ |
| Heats of adsorption | 30.85 (CO₂), 15.93 (N₂) kJ/mol | Table 1 ("measured", Clausius–Clapeyron average) | |
| Adsorbent C_p, wall ρ / C_p / k | 900 J/(kg K); 7800 kg/m³, 502 J/(kg K), 16 W/(m K) | Table 1 | |
| Gas and adsorbed-phase C_p | 30.7 J/(mol K), constant | Table 1 | the engine takes the gas C_p from CoolProp; the adsorbed-phase value is set to 30.7 |
| Viscosity, D_m | 1.72 × 10⁻⁵ Pa s, 1.60 × 10⁻⁵ m²/s, constant | Table 1 | the engine computes both per cell |
| K_z | 0.0903 W/(m K) | Table 1 | |
| h_in, h_out | 8.6 (fitted to the DCB temperature curves, Sect. 4.2), 10 ("assumed") W/(m² K) | Table 1 | |
| Isotherm | competitive DSL, concentration form | Table 2, Sect. 4.1 and Fig. 2c | q_sb 2.768, q_sd 2.758 mol/kg; CO₂ b₀ 3.934 × 10⁻⁷, d₀ 1.256 × 10⁻⁷ m³/mol, −ΔU_b 40.98, −ΔU_d 34.19 kJ/mol; N₂ b₀ = d₀ = 3.134 × 10⁻⁶, 15.71 / 15.72 kJ/mol |
| LDF coefficient | **macropore, no capacity factor** | SI Table S1, Eq. S4: k_i = 15 ε_p D_p / r_p² · c_i / q_i\* | the capacity factor is in the paper's model (see "Kinetics check") |
| Pressure drop | **Kozeny–Carman (180)** | SI Eq. S5: Blake–Kozeny with 150; Sect. 3 assumption (Darcy adequate) | |
| Wall end condition | **zero flux** | SI Table S2: T_w = T_a at z = 0 and z = L in every step | the engine has a fixed-end option |
| Pressurisation profile | linear from P_evac to P_ads over the step | Sect. 4.3 and Table S2 | as the deck |
| Blowdown profile | **exponential reaching 99 % of the change at the step end** (λ = ln 100 / 19 s = 0.242 s⁻¹) | Sect. 4.3: P = P_int + (P_ads − P_int) e^(−λ_BD t), λ_BD = 0.1742 s⁻¹ (Table 1, fitted to one experiment) | the deck's decay is 1.4× faster |
| Evacuation profile | **single exponential, 99 % at the step end** (λ = ln 100 / 57 s = 0.081 s⁻¹) | Sect. 4.3: two-term, P = P_L + (P_H1 − P_L) e^(−λ₁ t) + (P_H2 − P_L) e^(−λ₂ t), λ₁ = 0.019, λ₂ = 0.101 s⁻¹ (Table 1); the split P_H1 / P_H2 is not printed | the slow term (53 s time constant) keeps the experiment's pressure above P_low for most of the step; the deck reaches it |
| Step switching | fixed step times (19 s blowdown, 57 s evacuation, 30 s pressurisation) | Sect. 5.2: blowdown, evacuation and pressurisation end on a **pressure trigger**; Sect. 6.1: evacuation lengths differed by about 10 s between cycles | the deck's step times are not printed in the paper; their source is unrecorded |
| Idle | none | Sect. 6: a one-second hold with all valves closed between steps | |
| Feed | 15 % CO₂ / 85 % N₂, 298 K, interstitial velocity of Table 3 × ε | Sect. 2, Table 3 | the deck converts v_feed (Table 3, interstitial) to superficial with ε = 0.41 |
| Dead volume | none | Sect. 2: about 5 % of the column volume (12 cm³ of 240 cm³ measured for the DCB rig) | |
| CSS | component balance to 0.5 % | Sect. 3: < 1 % for each species over 5 consecutive cycles, 80–150 cycles | |
| Cells | 30 | not printed | |

Bold entries differ from the paper's model. Optimisation set-up (Sect. 5.2, for
WP11): decision variables t_ads 20–200 s, P_int 3–50 kPa, P_low 3–10 kPa,
v_feed 0.05–1.5 m/s; NSGA-II, population 96, 50 generations, 4800 evaluations
per front, each to CSS; objectives purity and recovery; P_low converged near
its 3 kPa bound for every point. Rig: Pfeiffer MVP 040-2 diaphragm pumps,
thermocouples at 8, 15 and 27 cm from the inlet.

### Attribution runs (2026-10-31)

One model-form change at a time on the bundled decks (Python, same engine), CO₂
purity / recovery at CSS in %. Experiment and paper model from Table 3.

| Config A (exp 99.7 / 58.7; paper model 99.3 / 60.0) | purity | recovery |
|---|---|---|
| deck before the audit (single exponential reaching 99 % of the evacuation pressure at the step end, Kozeny–Carman 180, zero-flux wall ends, plain macropore LDF) | 98.94 | 68.23 |
| + wall ends at T_a | 98.94 | 68.23 |
| + Ergun (Blake–Kozeny 150) | 98.96 | 68.56 |
| + blowdown at the paper's λ_BD = 0.1742 s⁻¹ (does not reach P_int in the 19 s step) | 96.89 | 70.80 |
| + capacity-factor LDF | 98.95 | 67.76 |
| + two-term evacuation, weight 0.1 / 0.3 / 0.5 / 0.7 / 1.0 in the slow term | 99.02 / 99.24 / 99.45 / — / 99.89 | 67.00 / 63.20 / 59.12 / — / 47.92 |
| two-term evacuation at 0.5, wall ends, capacity LDF: **the decks now** (Ergun left out, see below) | **99.49** | **59.01** |

| Two-term evacuation, weight 0.3 / 0.5 / 0.7 | purity | recovery | experiment | paper model |
|---|---|---|---|---|
| B (deck before: 96.01 / 70.28) | 96.45 / 96.75 / 97.05 | 64.99 / 61.01 / 56.99 | 96.1 / 61.8 | 95.2 / 64.1 |
| C1 (before: 86.89 / 72.76) | 87.37 / 87.67 / 87.96 | 65.42 / 60.30 / 55.27 | 88.1 / 64.6 (C2: 86.9 / 71.9) | 83.4 / 68.2 |
| D1, LPP (before: 93.72 / 90.22) | 93.73 / 93.83 / 93.94 | 90.14 / 89.85 / 89.55 | 97.1 / 86.1 (D2: 94.5 / 83.4) | 96.1 / 86.8 |

The four-step recovery excess was the evacuation pressure history, not the
kinetics: the paper's two-term profile keeps the column above P_low for most of
the 57 s step where the deck's single exponential reached it, so less CO₂ is
drawn. The split between the two terms (P_H1, P_H2) is not printed; the decks
use an equal split (weight 0.5), which was **chosen on A** (it reproduces A's
measured purity and recovery to 0.2 and 0.7 pp) and then **predicts B** to 0.6
and 1.1 pp without further change. C1 lands 4.5 pp below experiment C1 and 12
pp below its repeat C2, which disagree with each other by 7 pp; the paper
attributes its recovery uncertainty to the pressure trigger ending the
evacuation about 10 s earlier or later between cycles. The LPP cycles, whose
167 s evacuation reaches P_low whatever the profile, move by less than 0.5 pp;
their 3 pp recovery excess is a separate item (the composition-coupled LPP
connection, as in the Nguyen 2023 decks). The capacity-factor LDF, the wall
ends and the pressure-drop constant are each worth under 0.5 pp. Ergun (the
Blake–Kozeny constant of SI Eq. S5 plus the inertial term) is left out of the
decks: it is worth 0.3 pp, the paper calls the pressure drop negligible, and it
fails the initial-condition calculation of the fast-feed connected
configuration of `ConfigD1_ConnectedBC` (reproduced in Python with every other
form switched off), which is an engine follow-up. The paper's
blowdown rate hurts purity because the deck's 19 s blowdown does not reach
P_int at that rate (the experiment ended the step on a pressure trigger, whose
duration is not printed), so the deck keeps its blowdown profile.

Open: digitising Fig. 6(c) (the measured pressure history of C1 at CSS) would
pin the split of the evacuation terms and remove the one chosen value.

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

### Kinetics check (Config A, 2026-10-31)

The deck's macropore LDF coefficient has no capacity factor (the engine's
`MassTransferCoefficientMacropore` computes c and ρ_b q\* per cell but
multiplies by a hard-coded 1.0), so its CO₂ coefficient is of order 50 s⁻¹
where the Haghpanah et al. (2013) form, k = 15 ε_p D_m c / (τ r_p² ρ_p q\*),
gives about 0.1 s⁻¹. Replacing the model by constant coefficients of that
magnitude (0.1631 s⁻¹ CO₂, 0.2044 s⁻¹ N₂, as reported for 13X by Haghpanah et
al. (2013) and to be confirmed against Perez 2019 Table 1 once the paper is in
`assets/`) on Config A gives, at CSS:

| LDF coefficient | CO₂ purity | CO₂ recovery |
|---|---|---|
| macropore, no capacity factor (deck) | 98.9 % | 68.2 % |
| constant 0.1631 / 0.2044 s⁻¹ | 94.0 % | 66.6 % |
| experiment (Table 3) | 99.7 % | 58.7 % |
| paper's model (Table 3) | 99.3 % | 60.0 % |

Slower kinetics lower the recovery only by 1.6 pp and cost 5 pp of purity, so
the 8–10 pp recovery excess is not a kinetics effect; the paper's own
capacity-factor form (now the `macropore_capacity` coefficient, which the decks
use) moves it by 0.5 pp. See "Attribution runs" below for the cause.
