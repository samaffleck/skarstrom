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
- Column: L = 0.3 m, Dᵢ = 0.032 m, ε_b = 0.41, 80 cells (WENO advection, the engine default; see "Grid and advection scheme")
- Initial bed: the feed gas (15 % CO₂ / 85 % N₂) at 1 bar and 298 K, loadings in equilibrium with it
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
| Gas holdup | bed voids only, ε_t = ε_b = 0.41 | SI Table S1, Eqs. S1, S2, S6: (1 − ε)/ε with ε = 0.41; ε_p appears only in the LDF coefficient (Eq. S4) | `adsorbent[13X].particle_void_fraction` = 0; `kinetics[C].macropore_void_fraction` = 0.35 and tortuosity 3 keep the rate (total voidage 0.6165 before 2026-11-07) |
| Particle radius | 7.5 × 10⁻⁴ m | Table 1 ("measured") | |
| Adsorbent mass | 0.162 kg per column | abstract, Sect. 2 | particle density derived: m / (V (1 − ε)) = 1138 kg/m³ |
| Heats of adsorption | 30.85 (CO₂), 15.93 (N₂) kJ/mol | Table 1 ("measured", Clausius–Clapeyron average) | |
| Adsorbent C_p, wall ρ / C_p / k | 900 J/(kg K); 7800 kg/m³, 502 J/(kg K), 16 W/(m K) | Table 1 | |
| Gas and adsorbed-phase C_p | 30.7 J/(mol K), constant | Table 1 | `fluid.gas_heat_capacity_mode` = constant, `fluid.gas_heat_capacity[*]` = 30.7 and `fluid.adsorbed_heat_capacity[*]` = 30.7 J/(mol K) |
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
| CSS | component balance to 0.5 % and end-of-cycle state change ≤ 2×10⁻⁵, on 3 consecutive cycles | Sect. 3: < 1 % for each species over 5 consecutive cycles, 80–150 cycles | |
| Cells | 80, WENO advection (30 upwind cells before 2026-11-04) | not printed | grid-converged: recovery within 0.1 pp of its grid limit |

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
decks: it is worth 0.3 pp and the paper calls the pressure drop negligible.
(It used to fail the initial-condition calculation of the fast-feed connected
configuration of `ConfigD1_ConnectedBC`: the feed switch opens the vent valve
on a bed pressurised 1325 Pa above it, and the Ergun face could not take the
surge. The flow boundary now solves its face explicitly, regression test
`ErgunValveSurgeTest`, and `ConfigD1_ConnectedBC` runs on Ergun.) The paper's
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

- Four-step cycles (A, B, C1): CO2 purity within 0.15-0.86 pp of the experiments and recovery within 0.24-2.2 pp. LPP cycles (D1, E, F): purity 0.9-3.0 pp low and recovery 1.6-3.4 pp high (gas held in the bed voids only since 2026-11-07, see below), at cyclic steady state (per-component balance plus the end-of-cycle state test; from the feed-loaded start A to F stop after 49, 57, 65, 94, 207 and 517 cycles).
- 2026-11-04, grid and advection scheme: see the subsection below; the decks now run WENO on 80 cells, and every recovery rose by 0.6-2.6 pp.
- 2026-11-02, CSS rule: the balance-only test used before stopped 26 to 166 cycles short of CSS (A to F after 44, 55, 71, 155, 352 and 674 cycles), while purity and recovery were still moving by up to 0.17 pp (F). Every number on the validation page moved by that much; the only bound crossed was A's recovery against the experiment (0.294 → 0.417 pp, bound 0.4 → 0.5). The bound follows the measured error at CSS; nothing in the model changed.
- adsim's macropore LDF coefficient omits the c/q* capacity factor used in the paper (much faster kinetics), and its Carman-Kozeny constant is 180 against the paper's 150. The decks now use the competitive DSL (Sect. 4.1) and the paper's K_z = 0.0903 W/(m K).

### Gas holdup in the bed voids only (2026-11-07)

The paper's balances hold the gas in the bed voids only: SI Table S1, Eqs.
S1, S2 and S6 use (1 − ε)/ε with ε = 0.41, and the particle voidage ε_p = 0.35
enters only the macropore LDF coefficient (Eq. S4). The decks held gas in the
total voidage ε_t = ε_b + (1 − ε_b) ε_p = 0.6165 (the engine's convention,
[model conventions](../../../../docs/model-conventions.md#total-voidage)). They
now set the adsorbent's `particle_void_fraction` to 0, so ε_t = ε_b = 0.41,
and keep `macropore_void_fraction` = 0.35 and tortuosity 3 in the LDF
coefficient, so the rate is unchanged (as the Haghpanah 2013 decks do).

Every result moves by at most 1.13 pp (F's purity). Purity rises in every
case; recovery rises in the four-step cycles (A 0.50, B 0.37, C1 0.16 pp) and
falls in the LPP cycles (0.35-0.48 pp), which brings D1, E and F closer to the
experiments and to the paper's model. Five bounds were crossed and are
loosened to the measured error rounded up (bold); no bound was tightened.
Before: total voidage; after: bed voids only; WENO on 80 cells, at CSS in both.

| Case | Quantity | Reference | Before | After | Error before → after | Bound before → after |
|---|---|---|---|---|---|---|
| ConfigA | CO2 purity vs experiment | 99.70 | 99.51 | 99.55 | 0.192 → 0.149 | 0.3 |
| ConfigA | CO2 recovery vs experiment | 58.70 | 59.64 | 60.14 | 0.944 → 1.443 | **1.2 → 1.5** |
| ConfigA | CO2 purity vs paper model | 99.30 | 99.51 | 99.55 | 0.208 → 0.251 | 0.3 |
| ConfigA | CO2 recovery vs paper model | 60.00 | 59.64 | 60.14 | 0.356 → 0.143 | 0.43 |
| ConfigB | CO2 purity vs experiment | 96.10 | 96.79 | 96.96 | 0.693 → 0.861 | **0.7 → 0.9** |
| ConfigB | CO2 recovery vs experiment | 61.80 | 61.66 | 62.03 | 0.144 → 0.235 | **0.2 → 0.3** |
| ConfigB | CO2 purity vs paper model | 95.20 | 96.79 | 96.96 | 1.593 → 1.761 | 1.8 |
| ConfigB | CO2 recovery vs paper model | 64.10 | 61.66 | 62.03 | 2.444 → 2.065 | 4.3 |
| ConfigC1 | CO2 purity vs experiment | 88.10 | 87.95 | 88.43 | 0.147 → 0.325 | **0.2 → 0.4** |
| ConfigC1 | CO2 recovery vs experiment | 64.60 | 62.23 | 62.39 | 2.370 → 2.209 | 2.9 |
| ConfigC1 | CO2 purity vs paper model | 83.40 | 87.95 | 88.43 | 4.553 → 5.025 | **4.8 → 5.1** |
| ConfigC1 | CO2 recovery vs paper model | 68.20 | 62.23 | 62.39 | 5.970 → 5.809 | 9.7 |
| ConfigC1 | CO2 purity vs repeat experiment C2 | 86.90 | 87.95 | 88.43 | 1.053 → 1.525 | – |
| ConfigC1 | CO2 recovery vs repeat experiment C2 | 71.90 | 62.23 | 62.39 | 9.670 → 9.509 | – |
| ConfigD1 | CO2 purity vs experiment | 97.10 | 93.79 | 94.08 | 3.306 → 3.016 | 3.9 |
| ConfigD1 | CO2 recovery vs experiment | 86.10 | 89.83 | 89.49 | 3.732 → 3.385 | 4.5 |
| ConfigD1 | CO2 purity vs paper model | 96.10 | 93.79 | 94.08 | 2.306 → 2.016 | 2.7 |
| ConfigD1 | CO2 recovery vs paper model | 86.80 | 89.83 | 89.49 | 3.032 → 2.685 | 3.7 |
| ConfigD1 | CO2 purity vs repeat experiment D2 | 94.50 | 93.79 | 94.08 | 0.706 → 0.416 | – |
| ConfigD1 | CO2 recovery vs repeat experiment D2 | 83.40 | 89.83 | 89.49 | 6.432 → 6.085 | – |
| ConfigE | CO2 purity vs experiment | 88.50 | 86.65 | 87.25 | 1.852 → 1.250 | 2.1 |
| ConfigE | CO2 recovery vs experiment | 91.70 | 94.36 | 93.92 | 2.663 → 2.216 | 3.2 |
| ConfigE | CO2 purity vs paper model | 90.40 | 86.65 | 87.25 | 3.752 → 3.150 | 4.3 |
| ConfigE | CO2 recovery vs paper model | 91.10 | 94.36 | 93.92 | 3.263 → 2.816 | 4 |
| ConfigF | CO2 purity vs experiment | 77.70 | 75.65 | 76.78 | 2.051 → 0.916 | 2.6 |
| ConfigF | CO2 recovery vs experiment | 95.60 | 97.69 | 97.21 | 2.087 → 1.606 | 2.6 |
| ConfigF | CO2 purity vs paper model | 79.90 | 75.65 | 76.78 | 4.251 → 3.116 | 4.4 |
| ConfigF | CO2 recovery vs paper model | 94.60 | 97.69 | 97.21 | 3.087 → 2.606 | 3.8 |

Note (not applied here): with the new holdup several bounds are well above the
measured error and could be tightened in a separate change: A purity vs
experiment (0.149 pp, bound 0.3) and recovery vs paper model (0.143, 0.43);
B recovery vs paper model (2.07, 4.3); C1 recovery vs experiment (2.21, 2.9)
and vs paper model (5.81, 9.7); D1 all four (3.02 / 3.39 / 2.02 / 2.69 against
3.9 / 4.5 / 2.7 / 3.7); E all four (1.25 / 2.22 / 3.15 / 2.82 against
2.1 / 3.2 / 4.3 / 4.0); F all four (0.92 / 1.61 / 3.12 / 2.61 against
2.6 / 2.6 / 4.4 / 3.8).

Open: the four-step decks' equal split (weight 0.5) of the two evacuation
terms was chosen on A under the old (total-voidage) holdup (see "Attribution
runs"); it has not been re-checked with the gas in the bed voids only. The
decks keep 0.5 until that is done.

### Gas heat capacity as in the paper (2026-11)

The paper uses a constant gas heat capacity of 30.7 J/(mol K). The deck now uses the constant gas heat capacity mode with 30.7 J/(mol K) for CO₂ and N₂ (before: the default ideal-gas C_p0 at the local temperature). Every error moves by at most 0.005 pp: Config A recovery vs paper model 0.351 → 0.356 pp (bound 0.43), Config B purity vs experiment 0.691 → 0.693 pp (bound 0.7). No bound changed. All papers' settings: [../README.md](../README.md#gas-heat-capacity-per-paper).

### Grid and advection scheme (2026-11-04)

The decks ran first-order upwind advection on 30 cells, and a grid study
(each deck on 10-160 cells, Richardson extrapolation from the three finest
grids; [Numerics](../../../../docs/numerics.md#advection-schemes)) showed that
their recoveries were not grid-converged: on 30 upwind cells the CO2 recovery
was 0.65 (A), 1.4 (B), 2.8 (C1), 1.3 (D1) and 1.4 pp (E) below its grid limit,
and 0.1 pp needs 200 to over 1000 upwind cells. With WENO, now the engine's
default, 23-65 cells reach 0.1 pp (A 49, B 51, C1 65, D1 26, E 23; F within
0.1 pp from 20 cells), so every deck now runs WENO on 80 cells. Purity was
already within about 0.1 pp at 30 upwind cells except for C1 (0.4 pp).

The decks also start from a feed-loaded bed (owner's decision, 2026-11-04):
the feed gas (15 % CO2) with the loadings in equilibrium with it, instead of a
clean N2 bed. The cyclic steady state is the same (F from a clean bed after
2276 cycles: 75.58 % / 97.69 %; from the feed-loaded bed after 517: 75.65 % /
97.69 %, the difference being the CSS tolerance), but it is reached in far
fewer cycles.

Every comparison moved (before: upwind, 30 cells, clean start; after: WENO,
80 cells, feed-loaded start; at CSS in both). The recoveries rose by 0.6-2.6 pp, which brings the four-step
cycles B and C1 closer to the experiments (B 1.5 → 0.14 pp, C1 5.0 → 2.4 pp)
and moves A and the LPP cycles D1, E and F further above them (A 0.34 → 0.95,
D1 2.5 → 3.7, E 1.3 → 2.7, F 0.54 → 2.1 pp). F's purity fell by 0.4 pp. The
earlier, closer agreement of A, D1, E and F was partly the upwind scheme's
numerical dispersion; the bounds follow the grid-converged errors (bold:
changed; loosened where the error grew past the old bound, tightened where it
fell well below it). Cycles to CSS from the feed-loaded start: A 49, B 57,
C1 65, D1 94, E 207, F 517 (from a clean bed with WENO: 60, 81, 99, 233, 690,
2276; upwind from a clean bed: 70, 81, 97, 190, 426, 840: the sharper fronts
slow the approach, the feed-loaded start more than makes up for it). The caps
are about 30 % above the feed-start stop for E (270) and F (680); A-D1 keep
250.

| Case | Quantity | Reference | Before | After | Error before → after | Bound before → after |
|---|---|---|---|---|---|---|
| ConfigA | CO2 purity vs experiment | 99.70 | 99.48 | 99.51 | 0.22 → 0.19 | 0.3 |
| ConfigA | CO2 recovery vs experiment | 58.70 | 59.04 | 59.65 | 0.34 → 0.95 | **0.5 → 1.2** |
| ConfigA | CO2 purity vs paper model | 99.30 | 99.48 | 99.51 | 0.18 → 0.21 | 0.3 |
| ConfigA | CO2 recovery vs paper model | 60.00 | 59.04 | 59.65 | 0.96 → 0.35 | **1.2 → 0.43** |
| ConfigB | CO2 purity vs experiment | 96.10 | 96.69 | 96.79 | 0.59 → 0.69 | 0.7 |
| ConfigB | CO2 recovery vs experiment | 61.80 | 60.28 | 61.66 | 1.52 → 0.14 | **1.7 → 0.2** |
| ConfigB | CO2 purity vs paper model | 95.20 | 96.69 | 96.79 | 1.49 → 1.59 | 1.8 |
| ConfigB | CO2 recovery vs paper model | 64.10 | 60.28 | 61.66 | 3.82 → 2.44 | 4.3 |
| ConfigC1 | CO2 purity vs experiment | 88.10 | 87.49 | 87.95 | 0.61 → 0.15 | **0.62 → 0.2** |
| ConfigC1 | CO2 recovery vs experiment | 64.60 | 59.60 | 62.23 | 5.00 → 2.37 | **5.5 → 2.9** |
| ConfigC1 | CO2 purity vs paper model | 83.40 | 87.49 | 87.95 | 4.09 → 4.55 | 4.8 |
| ConfigC1 | CO2 recovery vs paper model | 68.20 | 59.60 | 62.23 | 8.60 → 5.97 | 9.7 |
| ConfigC1 | CO2 purity vs repeat experiment C2 | 86.90 | 87.49 | 87.95 | 0.59 → 1.05 | – |
| ConfigC1 | CO2 recovery vs repeat experiment C2 | 71.90 | 59.60 | 62.23 | 12.30 → 9.67 | – |
| ConfigD1 | CO2 purity vs experiment | 97.10 | 93.72 | 93.79 | 3.38 → 3.31 | 3.9 |
| ConfigD1 | CO2 recovery vs experiment | 86.10 | 88.55 | 89.84 | 2.45 → 3.74 | **3.5 → 4.5** |
| ConfigD1 | CO2 purity vs paper model | 96.10 | 93.72 | 93.79 | 2.38 → 2.31 | 2.7 |
| ConfigD1 | CO2 recovery vs paper model | 86.80 | 88.55 | 89.84 | 1.75 → 3.04 | **2.7 → 3.7** |
| ConfigD1 | CO2 purity vs repeat experiment D2 | 94.50 | 93.72 | 93.79 | 0.78 → 0.71 | – |
| ConfigD1 | CO2 recovery vs repeat experiment D2 | 83.40 | 88.55 | 89.84 | 5.15 → 6.44 | – |
| ConfigE | CO2 purity vs experiment | 88.50 | 86.61 | 86.65 | 1.89 → 1.85 | 2.1 |
| ConfigE | CO2 recovery vs experiment | 91.70 | 92.97 | 94.37 | 1.27 → 2.67 | **2.5 → 3.2** |
| ConfigE | CO2 purity vs paper model | 90.40 | 86.61 | 86.65 | 3.79 → 3.75 | 4.3 |
| ConfigE | CO2 recovery vs paper model | 91.10 | 92.97 | 94.37 | 1.87 → 3.27 | **3.2 → 4** |
| ConfigF | CO2 purity vs experiment | 77.70 | 76.06 | 75.65 | 1.64 → 2.05 | **1.9 → 2.6** |
| ConfigF | CO2 recovery vs experiment | 95.60 | 96.14 | 97.69 | 0.54 → 2.09 | **1.9 → 2.6** |
| ConfigF | CO2 purity vs paper model | 79.90 | 76.06 | 75.65 | 3.84 → 4.25 | 4.4 |
| ConfigF | CO2 recovery vs paper model | 94.60 | 96.14 | 97.69 | 1.54 → 3.09 | **3 → 3.8** |

### Paper protocol (report only, 2026-11-05)

The paper stops its runs once each species' balance closes to 1 % over 5
consecutive cycles (Sect. 3), with no test of the bed state. The LPP tests
(D1, E, F) therefore re-run each deck under that rule after the CSS run (Mass
Balance (Per-Component), tolerance 0.01, 5 consecutive cycles, `state_tolerance`
0, same feed-loaded start) and print the stop cycle and KPIs next to the CSS
result; nothing is asserted and the validation page keeps the CSS values.

| Case | Paper rule: stop / CO2 purity / recovery | CSS: stop / purity / recovery | Paper model | Experiment |
|---|---|---|---|---|
| D1 | 44 / 93.83 % / 89.94 % | 94 / 93.79 % / 89.84 % | 96.1 % / 86.8 % | 97.1 % / 86.1 % |
| E | 84 / 86.79 % / 94.47 % | 207 / 86.65 % / 94.37 % | 90.4 % / 91.1 % | 88.5 % / 91.7 % |
| F | 152 / 76.19 % / 97.75 % | 517 / 75.65 % / 97.69 % | 79.9 % / 94.6 % | 77.7 % / 95.6 % |

The paper's rule stops two to three times earlier; F's purity is then 0.5 pp
higher (still drifting down) and every recovery about 0.1 pp higher. So the
protocol does not explain the gap to the paper's model (D1-F recoveries 2-3 pp
above it, purities 2-4 pp below), and the CSS values stand.

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
