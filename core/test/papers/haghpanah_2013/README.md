# Haghpanah et al. (2013) — four-step CO₂/N₂ PVSA on zeolite 13X

The four-step pressure–vacuum swing cycle for post-combustion CO₂ capture
(pressurisation with feed, high-pressure adsorption, forward blowdown, reverse
evacuation) of Haghpanah et al. (2013). It is checked against the paper's own
results and against the independent re-implementation of the same model by
Glover et al. (2026), which also gives the dominant multiplier of the cycle map
at a benchmark design.

## Citation

> R. Haghpanah, A. Majumder, R. Nilam, A. Rajendran, S. Farooq, I.A. Karimi,
> M. Amanullah. *Multiobjective Optimization of a Four-Step Adsorption Process
> for Postcombustion CO₂ Capture Via Finite Volume Simulation.* Ind. Eng. Chem.
> Res. 52 (2013) 4249–4265. DOI: [10.1021/ie302658y](https://doi.org/10.1021/ie302658y)

> T. Glover, T. Papathanasiou, R. Pini. *Accelerating Simulation and
> Optimisation of Cyclic Adsorption Processes with Differentiable Programming.*
> arXiv:2606.21523 (2026).

Isotherm: the fit of Krishnamurthy, Haghpanah, Rajendran and Farooq (APCChE
2012 conference), as cited by Haghpanah et al. (2013), ref. 38, and printed in
their Table 2.

## Files

| File | Purpose |
|------|---------|
| `haghpanah_2013.h` | `Haghpanah2013Test` fixture; `BuildHaghpanah` builds the whole deck for one design |
| `haghpanah_2013.cpp` | the cases below |
| `reference_values.csv` | the printed numbers compared with (both papers) |

## Reproduced cases

| Test (`Haghpanah2013Test.*`) | Design | Compared with | Deck |
|---|---|---|---|
| `Sec61` | Sect. 6.1 / Table 4: t = 15/15/30/40 s, P_H 1 bar, P_I 0.2 bar, P_L 0.1 bar, v_feed 1 m/s; 30 cells | Haghpanah's purity 88.74 % and recovery 35.83 % (p. 4257) | `Haghpanah2013_Sec61` |
| `PL002` … `PL010` | the minimum-energy designs of Table 6 (= Glover Table S2) at P_L = 0.02, 0.03, 0.05, 0.07, 0.1 bar; t_press 20 s; 30 cells | Haghpanah's energies (Table 6), Glover's purity, recovery, productivity and energy (Table 4) | `Haghpanah2013_PL002` … `_PL010` |
| `Fig5Benchmark` | Glover's CSS benchmark (Fig. 5 caption): P_H 8, P_I 1.5, P_L 0.05 bar, v_feed 0.8 m/s, t = 20/50/100/100 s; 10 cells (Glover's grid) | the dominant multiplier ρ = 0.9509 (Glover p. 27); cycles to CSS from a feed- and an N₂-saturated bed (82 and 115, pp. 23–24) | `Haghpanah2013_Fig5Benchmark` |

Every case starts from a bed full of N₂ at P_L and 298.15 K (Sect. 6.1) and runs
to CSS (Mass Balance (Per-Component), default state tolerance, at most 500
cycles as in Sect. 6.2). Runtime: Sec61 1 min (465 cycles), the Table 6 cases
30–40 s each, Fig5Benchmark 45 s (two runs of 200 cycles on 10 cells).

## Provenance

Sources: **H** = Haghpanah et al. (2013), **G** = Glover et al. (2026). Paths are
the registry paths of the shipped decks (`skarstrom params <deck> --all`).

| Value | Deck parameter (registry path) | Source | Note |
|---|---|---|---|
| 1 m | `layer[1].length` | H Table 3 | |
| inner radius 0.1445 m → diameter 0.289 m | `bed.diameter` | H Table 3 | G prints 0.145 m (0.7 % in area); per-volume and per-tonne KPIs do not depend on it at a fixed v_feed |
| wall thickness 0.0175 m | `wall.thickness` | H Table 3: r_out 0.162 m | derived: 0.162 − 0.1445 |
| ε = 0.37 | `layer[1].bed_void_fraction` | H Table 3 | |
| particle density 1130 kg/m³ | `adsorbent[1].solid_density` | H Table 3 (ρ_s, the pellet including its pores) | ρ_s (1 − ε) = 711.9 kg/m³ = G's bed density 712 |
| **particle voidage in the holdup: 0** | `adsorbent[1].particle_void_fraction` | H Eqs. A.2–A.3: gas is held in ε only | ε_t = ε_b (docs/model-conventions.md "Total voidage") |
| ε_p = 0.35, τ = 3 in the rate | `adsorbent[1].kinetics[C].macropore_void_fraction`, `.tortuosity` | H Table 3 | |
| r_p = 1 mm | `layer[1].particle_diameter` = 2 mm | H Table 3 | |
| C_p,s 1070 J/(kg K) | `adsorbent[1].heat_capacity` | H Table 3 | |
| C_p,g = C_p,a = 30.7 J/(mol K), constant | `fluid.gas_heat_capacity_mode` = constant, `fluid.gas_heat_capacity[*]`, `fluid.adsorbed_heat_capacity[*]` | H Table 3 | |
| μ = 1.72 × 10⁻⁵ Pa s, constant | `fluid.viscosity_mode` = constant, `fluid.viscosity` | H Table 3; App. A (μ independent of T) | |
| D_m = 1.60 × 10⁻⁵ m²/s, constant | `adsorbent[1].kinetics[C].molecular_diffusion_coefficient` | H Table 3 | G prints 1.5 × 10⁻⁷, which is inconsistent with its own k (its Eq. S36 gives 26.25 1/s only for 1.5 × 10⁻⁵): not used |
| D_L = 0.7 D_m + v_feed r_p | `layer[1].dispersion_coefficient` (constant model) | H Eq. 42 (G Eq. 14) | derived per design, e.g. 1.0112 × 10⁻³ m²/s at v_feed 1 m/s |
| K_z = 0.09 W/(m K) | `layer[1].heat_dispersion_coefficient` | H Table 3 | |
| h_in 8.6, h_out 2.5 W/(m² K) | `layer[1].internal_htc`, `layer[1].external_htc` | H Table 3 | |
| wall ρ 7800 kg/m³, C_p 502 J/(kg K), k 16 W/(m K) | `wall.density`, `wall.heat_capacity`, `wall.thermal_conductivity` | H Table 3 | |
| wall ends at T_a = 298.15 K | `wall.boundary` = fixed, `wall.fixed_temperature` | H Eq. 12, SI Table 2 | |
| T_feed = T_a = 298.15 K | stream and boundary temperatures, `bed.ambient_temperature` | H Table 3 | |
| feed 15 % CO₂ / 85 % N₂ | stream `y` | H Sect. 6.1 | |
| isotherm: CO₂ q_sb 3.09, q_sd 2.54 mol/kg, b₀ 8.65 × 10⁻⁷, d₀ 2.63 × 10⁻⁸ m³/mol, ΔU_b −36 641.21, ΔU_d −35 690.66 J/mol; N₂ q_sb 5.84 mol/kg, b₀ 2.50 × 10⁻⁶ m³/mol, ΔU_b −15 800 J/mol, site d empty | `adsorbent[1].isotherm[C].m1/b1/q1/m2/b2/q2` | H Table 2 (Krishnamurthy et al. 2012) | competitive dual-site Langmuir in concentration form (H Eqs. 37–39), the engine's `extended_dual_site_langmuir_c` |
| LDF k_i = (c_i/q_i\*) 15 ε_p D_m / (τ r_p²) | `adsorbent[1].kinetics[C]` = `macropore_capacity` | H Eqs. 40–41 | q\* per particle volume (ρ_s q\*), the dimensionally consistent reading; the paper does not state the basis (see "Attribution") |
| **heat of adsorption, constant: CO₂ 35 983.8, N₂ 15 800 J/mol** | `adsorbent[1].heat_of_adsorption[C]` | **assumption** (not printed; checked: Tables 2–4, Sect. 2, App. A; G p. S-13 says the paper does not specify it) | CO₂: dual-site isosteric weighting (|ΔU_b| ∂q_b/∂c + |ΔU_d| ∂q_d/∂c)/(∂q_b/∂c + ∂q_d/∂c) at the feed partial pressure 0.15 bar, 298.15 K (as the Ferré 2026 decks); N₂: its single site. The engine has no dynamic Clausius–Clapeyron heat (G Eqs. 15–16) |
| heat-of-adsorption convention: internal energy | `fluid.heat_of_adsorption_convention` | H Eq. A.8 carries −C_p,a T Σ ∂q/∂t | |
| pressure drop: Darcy with the Blake–Kozeny constant 150 | `layer[1].momentum_model` = blake_kozeny | H Eq. 6 | |
| exponential pressure profiles, λ = 0.5 1/s | step boundary pressure: `double_exponential`, both rates 0.5 1/s, weight 1 | H Eqs. 57–59, Table 4 | exactly P_end + (P_start − P_end) e^(−λt), from the bed pressure at the start of the step |
| interstitial feed velocity v_feed | stream velocity (superficial) = ε v_feed | H Table 3 ("interstitial feed velocity") | |
| t_press = 20 s for the Table 6 designs | `step[1].time` | H Sect. 7.1, G p. 21 | not listed in Table 6 / Table S2 |
| P_L = 0.05 bar, t_press = 20 s for Fig5Benchmark | | G p. 21 (fixed in G's optimisation) | not repeated in the Fig. 5 caption |
| compression and evacuation: isentropic, γ = 1.4, η = 0.72, from / to 1 bar | stream and boundary `pump` (compressor on the pressurisation and feed inlets, vacuum pump on the blowdown and evacuation outlets), energy model isothermal | H Eqs. 66–69, p. 4260; G p. S-13 (constant 0.72 for the verification) | the outlet-temperature model only sets the temperature of the gas entering the bed; the power uses γ either way, so the feed enters at 298.15 K as in the paper |
| initial bed: N₂ at P_L, 298.15 K | `initial_conditions.*` | H Sect. 6.1 | |

## Model-form audit

| Item | Paper | Deck |
|---|---|---|
| Isotherm and competition | competitive DSL, shared denominators per site (H Eqs. 37–39) | same |
| Rate model | LDF, k from Eq. 40 with c_i/q_i\* | `macropore_capacity`, which caps c/(ρ_s q\*) at 1 (N₂ only reaches it; the cap changes no KPI by more than 0.01 pp: a constant k_N₂ = 42 1/s gives the same results) |
| Gas holdup | ε only (Eqs. A.2–A.3) | ε_t = ε_b (`particle_void_fraction` 0) |
| Axial dispersion | constant D_L per design (Eq. 42) | constant |
| Axial conduction | constant K_z | constant |
| Wall | wall energy balance, ends at T_a | same |
| Pressure drop | Darcy, Blake–Kozeny 150 | same |
| Energy balance | H Eq. A.8: no gas-phase accumulation term, −(C_p,g/R) ∂P/∂t on the right | the engine's column energy balance (gas, solid and adsorbed phases in local equilibrium). An isothermal run moves the KPIs by ≤ 0.15 pp and the energy by 1 % (Sect. 6.1), so the form is not what limits the agreement |
| Heat of adsorption | not specified | constant (above) |
| Inlet composition / temperature | Danckwerts at the feed end (SI Table 2) | the engine's inlet boundary |
| KPIs | boundary-flux integrals (H Eqs. 60–64, 66–71) | the engine's boundary-flux KPIs: purity and recovery of the evacuation product, CO₂ fed in pressurisation and adsorption, energy = time-integrated pump power per kg of CO₂ product (J/kg ÷ 3600 = kWh/t), productivity = CO₂ product per column volume per cycle time |
| Scheme | WENO, 30 cells (optimisation), 1000 (Sect. 6.1 reference) | WENO, 30 cells (10 for Fig5Benchmark) |

### Grid

30 against 60 cells (the deck, once, 2026-11-07):

| Case | Metric | 30 cells | 60 cells |
|---|---|---|---|
| Sec61 | purity / recovery | 77.04 % / 37.98 % | 77.05 % / 38.01 % |
| Sec61 | energy | 141.28 kWh/t | 141.20 kWh/t |
| PL002 | purity / recovery | 90.16 % / 87.75 % | 90.13 % / 87.82 % |
| PL002 | energy, productivity | 169.12 kWh/t, 0.2915 mol/(m³ s) | 168.73 kWh/t, 0.2918 mol/(m³ s) |

(Measured before the macropore holdup was removed.) Doubling the grid moves
every KPI by ≤ 0.07 pp or 0.25 %: 30 cells, Haghpanah's optimisation grid, is
converged for these comparisons.

## Measured errors (2026-11-07)

Purity and recovery in percentage points (Glover prints fractions: 0.906 is
stored as 90.6 %), the others relative. Bounds are regression bounds set from the
measured error; a metric is asserted where it is within the 5 % target, and
recorded without a bound (–) otherwise.

| Case | Metric | Paper | Engine | Error | Bound |
|---|---|---|---|---|---|
| Sec61 | CO₂ purity vs H | 88.74 % | 79.68 % | 9.06 pp | – |
| Sec61 | CO₂ recovery vs H | 35.83 % | 37.42 % | 1.59 pp | 1.9 |
| PL002 | energy vs H Table 6 | 148.96 kWh/t | 167.46 | 12.4 % | – |
| PL002 | CO₂ purity vs G | 90.6 % | 90.86 % | 0.26 pp | 0.35 |
| PL002 | CO₂ recovery vs G | 89.4 % | 87.86 % | 1.54 pp | 1.8 |
| PL002 | productivity vs G | 0.297 mol/(m³ s) | 0.2888 | 2.8 % | 3.3 % |
| PL002 | productivity vs G Ref. | 0.312 | 0.2888 | 7.4 % | – |
| PL002 | energy vs G | 143.6 kWh/t | 167.46 | 16.6 % | – |
| PL002 | productivity vs H (per adsorbent volume) | 0.49 mol/(m³ ads s) | 0.4584 | 6.5 % | – |
| PL003 | energy vs H Table 6 | 201.50 | 225.31 | 11.8 % | – |
| PL003 | CO₂ purity vs G | 90.6 % | 90.94 % | 0.34 pp | 0.45 |
| PL003 | CO₂ recovery vs G | 90.7 % | 89.17 % | 1.53 pp | 1.8 |
| PL003 | productivity vs G | 0.508 | 0.4910 | 3.4 % | 4 % |
| PL003 | productivity vs G Ref. | 0.512 | 0.4910 | 4.1 % | 4.8 % |
| PL003 | energy vs G | 194.5 | 225.31 | 15.8 % | – |
| PL005 | energy vs H Table 6 | 298.92 | 317.85 | 6.3 % | – |
| PL005 | CO₂ purity vs G | 90.3 % | 89.69 % | 0.61 pp | 0.75 |
| PL005 | CO₂ recovery vs G | 90.9 % | 90.73 % | 0.17 pp | 0.25 |
| PL005 | productivity vs G | 0.333 | 0.3305 | 0.7 % | 1.2 % |
| PL005 | productivity vs G Ref. | 0.348 | 0.3305 | 5.0 % | – |
| PL005 | energy vs G | 289.3 | 317.85 | 9.9 % | – |
| PL007 | energy vs H Table 6 | 388.68 | 402.55 | 3.6 % | 4.2 % |
| PL007 | CO₂ purity vs G | 91.0 % | 90.57 % | 0.43 pp | 0.55 |
| PL007 | CO₂ recovery vs G | 91.4 % | 90.78 % | 0.61 pp | 0.75 |
| PL007 | productivity vs G | 0.704 | 0.6910 | 1.8 % | 2.3 % |
| PL007 | productivity vs G Ref. | 0.688 | 0.6910 | 0.4 % | 1 % |
| PL007 | energy vs G | 373.2 | 402.55 | 7.9 % | – |
| PL010 | energy vs H Table 6 | 499.65 | 503.23 | 0.7 % | 1.2 % |
| PL010 | CO₂ purity vs G | 91.8 % | 90.65 % | 1.15 pp | 1.35 |
| PL010 | CO₂ recovery vs G | 91.1 % | 91.57 % | 0.47 pp | 0.6 |
| PL010 | productivity vs G | 0.721 | 0.7213 | 0.04 % | 0.5 % |
| PL010 | productivity vs G Ref. | 0.694 | 0.7213 | 3.9 % | 4.6 % |
| PL010 | energy vs G | 486.8 | 503.23 | 3.4 % | 4 % |
| Fig5Benchmark | dominant multiplier (fitted decay) | 0.9509 | 0.95131 | 0.00041 | 0.0015 |
| Fig5Benchmark | engine CSS rate, median over the window | 0.9509 | 0.95134 | 0.00044 | – |
| Fig5Benchmark | cycles to a 10⁻⁵ state change, N₂ start | 115 | 99 | 14 % | order of magnitude |
| Fig5Benchmark | cycles to a 10⁻⁵ state change, feed start | 82 | 72 | 12 % | order of magnitude |

Glover's "Ref." productivities are not printed in Haghpanah et al. (2013)
(presumably read from its Fig. 12), so they are soft references. Glover's
productivity is per column volume (its Eq. 55) and Haghpanah's per adsorbent
volume: 0.312 / (1 − ε) = 0.495 reproduces Haghpanah's 0.49, so the engine's
column-volume value is divided by (1 − ε) for that comparison.

### CSS multiplier (Fig5Benchmark)

The deck runs 200 cycles from each start with the stopping rule held off, and
the end-of-cycle state change of every cycle is read from the log. From the N₂
start, the change decays geometrically from cycle 52 (below 10⁻⁴) to cycle 143
(above 10⁻⁶, a few times the noise floor of the state change at reltol 10⁻⁶).
The least-squares slope of ln(change) over that window gives ρ = 0.95131; the
CSS check's own rate estimate (geometric mean over the last five changes) has a
median of 0.95134 and spans 0.95043–0.95166 over the window, so the bound of
0.0015 is the diagnostic's own scatter. From the feed start the window is
cycles 31–115 and the fit 0.94824: that start is closer to the CSS, and the
dominant mode is not yet alone at the start of its window. Glover's 0.9509 is the
spectral radius of the sensitivity matrix, from its own model (constant LDF
coefficient, its energy balance); the engine's value is that of the
Haghpanah-form deck above.

Cycles to a 10⁻⁵ state change are 99 (N₂ start) and 72 (feed start) against
Glover's 115 and 82. The criteria differ (Glover: the largest relative change of
any state; the engine: its norm of the end-of-cycle state change), so only the
order of magnitude and the order of the two starts are asserted. The engine's
own stopping rule stops at cycles 86 and 60.

## Attribution

### Specific energy (PL002, PL003)

The energy is 0.7–6 % above Table 6 for P_L = 0.05–0.1 bar and 12 % above for
0.02–0.03 bar. The productivity and purity at PL002 are within 3 % and 0.3 pp
of Glover, so about the same CO₂ is produced; the excess is in how deep a vacuum
it leaves at. At PL002 the evacuation pump does 95 % of the work (200 kJ per
cycle for 8.9 mol, 22.4 kJ/mol; the isentropic work per mole at 0.02 bar is
24.8 kJ/mol and at 0.05 bar 16.3). One-factor runs on PL002 (energy in kWh/t):

| Change | Energy | Recovery | Note |
|---|---|---|---|
| deck | 167.5 | 87.9 % | |
| 60 cells | 168.7 | 87.8 % | (before the holdup change, deck 169.1) |
| isosteric heat-of-adsorption convention | 167.0 | 87.1 % | (before the holdup change) |
| macropore holdup on (ε_p = 0.35 in ε_t) | 169.1 | 87.8 % | |
| constant LDF 26.25 1/s (Glover) | 135.4 | 91.6 % | (with the holdup on) |
| D_m in the rate × 10 | 147.2 | 90.9 % | |
| D_m in the rate × 163.6 | 134.2 | 91.3 % | |

The energy depends on the CO₂ rate: with the coefficient of Eq. 40 read per
particle volume, k_CO₂ ≈ 0.04 1/s at the feed, CO₂ desorbs late and leaves at
the deepest vacuum. Eq. 40 does not state the basis of q\*: the factor 163.6 is
R T₀ ρ_s q_s0 / P₀, what the ratio becomes if it is formed from the paper's
dimensionless variables (c̄ = c R T₀/P₀, x = q/q_s0). A factor of about 10
matches Table 6 at PL002, but every faster rate raises the Sect. 6.1 recovery to
46 % against the paper's 35.83 % (below), which the dimensional reading
reproduces. No single reading of Eq. 40 fits both of the paper's own results, so
the deck keeps the dimensionally consistent one and the energy at P_L ≤ 0.05 bar
is recorded, not asserted.

### Sect. 6.1 purity

Recovery matches (37.4 % against 35.83 %), purity does not (79.7 % against
88.74 %). One-factor runs on Sec61:

| Change | Purity | Recovery | Energy |
|---|---|---|---|
| deck | 79.68 % | 37.42 % | 131.8 kWh/t |
| isothermal | 79.83 % | 38.95 % | 130.7 |
| CO₂ heat of adsorption 36 641 J/mol (site b) | 79.68 % | 37.39 % | 131.8 |
| isosteric convention | 79.43 % | 36.47 % | 132.8 |
| no wall heat transfer | 79.70 % | 37.40 % | 131.7 |
| Ruthven dispersion (local velocity) | 79.84 % | 37.82 % | 131.1 |
| constant k_N₂ = 42 1/s (no cap) | 79.67 % | 37.41 % | 131.6 |
| D_m in the rate × 10 | 81.85 % | 45.55 % | 114.2 |
| D_m in the rate × 163.6 | 81.87 % | 45.83 % | 111.6 |
| constant LDF 26.25 1/s (with the holdup on) | 80.01 % | 47.72 % | 117.4 |
| macropore holdup on | 77.04 % | 37.98 % | 141.3 |
| 60 cells (holdup on) | 77.05 % | 38.01 % | 141.2 |

Nothing thermal, the dispersion closure, the grid or the N₂ rate moves the
purity by more than 0.3 pp, and the CO₂ rate by about 2 pp at the cost of the
recovery. The design has a small evacuation swing (0.2 → 0.1 bar) after a short
adsorption step whose CO₂ front reaches the product end (p. 4257), so the purity
is set by the N₂ left in the bed after blowdown; the cause of the 9 pp is open.
The Table 6 purities, from the same deck, agree with Glover within 1.2 pp.

### What is not modelled

- Glover's dynamic heat of adsorption (G Eqs. 15–16): the engine uses a
  constant heat per component.
- Haghpanah's energy-balance form (Eq. A.8): see the model-form audit; the
  isothermal run bounds its effect.

## Validation protocol

1. **Source**: both papers and their SI were read in full and every input and
   reference value transcribed; the numbers here are the printed ones (no value
   read off a figure).
2. **Provenance**: above; one assumption (the heat of adsorption), one
   misprint not used (G's D_m).
3. **Model form**: above; gas holdup corrected to the paper's (ε only).
4. **Data**: no experimental curves (both references are model results).
5. **Metrics**: the papers' KPI definitions (H Eqs. 60–71, G Eq. 55); purity
   and recovery in pp, the rest relative.
6. **Attribution**: above; two open disagreements (energy at P_L ≤ 0.05 bar,
   Sect. 6.1 purity).
7. **Report**: bounds from the measured errors; the validation page regenerated.
