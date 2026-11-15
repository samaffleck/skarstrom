# Ward & Pini (2022) — four-step CO₂/N₂ PVSA on zeolite 13X and ZIF-36-FRL

The four-step pressure–vacuum swing cycle with feed pressurisation (adsorption,
forward blowdown, reverse evacuation, feed pressurisation) that Ward and Pini
optimised for post-combustion capture from a 1000 MW coal plant. The paper's
process model is the Haghpanah et al. (2013) model, re-implemented
([../haghpanah_2013/](../haghpanah_2013/)), with a pressure-dependent pump
efficiency and a column that scales with its length. The SI prints the ten
cost-optimal designs completely, with the model's purity, recovery,
productivity and energy. The deck is run at each design and compared with these
numbers (a code-to-code comparison; there is no experiment).

This is not the Ward & Pini (2022) breakthrough paper in *Adsorption*
([../ward_2022/](../ward_2022/)).

## Citation

> A. Ward, R. Pini. *Efficient Bayesian Optimization of Industrial-Scale
> Pressure-Vacuum Swing Adsorption Processes for CO₂ Capture.* Ind. Eng. Chem.
> Res. 61 (2022) 13650–13668. DOI: [10.1021/acs.iecr.2c02313](https://doi.org/10.1021/acs.iecr.2c02313)

Sources below: **P** = the paper, **S** = its Supporting Information.

## Files

| File | Purpose |
|------|---------|
| `ward_pini_2022.h` | `WardPini2022Test` fixture; `BuildWardPini` builds the whole deck for one design |
| `ward_pini_2022.cpp` | the ten cases below |
| `reference_values.csv` | the printed numbers compared with (S Tables S2, S3) |

## Reproduced cases

| Test (`WardPini2022Test.*`) | Design | Deck |
|---|---|---|
| `Zeolite13X_PL002` … `_PL010` | S Table S2: the minimum-capture-cost design of zeolite 13X at p_L = 0.02, 0.03, 0.05, 0.07, 0.1 bar | `PL003` (= P Table 7): `WardPini2022_13X_PL003` |
| `ZIF36FRL_PL002` … `_PL010` | S Table S3: the same for ZIF-36-FRL | `PL007` (= P Table 7): `WardPini2022_ZIF36FRL_PL007` |

Each design fixes p_H, p_I, p_L, t_ads, t_bd, t_evac, the interstitial feed
velocity v_F and the column length L; t_pres is 20 s for every design (P Sect.
3.3). Every case starts from a bed full of N₂ at p_H and 298.15 K (P p. 13653)
and runs to CSS (Mass Balance (Per-Component), default state tolerance, at most
500 cycles): 120–170 cycles, 10–14 s each. The two P Table 7 designs are
shipped as example decks; the other eight are run only by the tests.

## Provenance

Paths are the registry paths of the shipped decks (`skarstrom params <deck> --all`).

| Value | Deck parameter (registry path) | Source | Note |
|---|---|---|---|
| L per design (1.57–4.18 m) | `layer[1].length` | S Tables S2, S3 | |
| **r_in = L / 6.43** | `bed.diameter` = 2 L / 6.43 | P Sect. 6.1: "the aspect ratio … is fixed … such that L/r_in = 6.43" | P Table 3 gives r_in = 0.145 m at L = 1 m (L / r_in = 6.9); the stated ratio is used. Per-volume and per-tonne KPIs depend on it only through the wall heat transfer |
| **wall thickness (0.162/0.145 − 1) r_in** | `wall.thickness` | **assumption** (not stated for the scaled columns; checked: P Sects. 2.3, 6.1, Table 3, S Sects. S.1, S.6) | the r_out / r_in of the Table 3 column, kept as the column scales |
| ε = 0.37 | `layer[1].bed_void_fraction` | P Table 3 | |
| ρ_b = 712 kg/m³ | `adsorbent[1].solid_density` = 712 / (1 − 0.37) = 1130.2 kg/m³ | P Table 3 | derived; the same for both materials (P Sect. 2.3) |
| **particle voidage in the holdup: 0** | `adsorbent[1].particle_void_fraction` | S Eq. 1 (ψ = R T₀ ρ_b q_s0 / (p₀ ε)): gas is held in ε only | ε_t = ε_b (docs/model-conventions.md "Total voidage") |
| ε_p = 0.35, τ_p = 3 in the rate | `adsorbent[1].kinetics[C].macropore_void_fraction`, `.tortuosity` | P Table 3 | |
| r_p = 1 mm | `layer[1].particle_diameter` = 2 mm | P Table 3 | |
| D_m = 1.5 × 10⁻⁵ m²/s | `adsorbent[1].kinetics[C].molecular_diffusion_coefficient` | P Table 3 | Haghpanah prints 1.6 × 10⁻⁵ |
| D_L = 0.7 D_m + v_F r_p | `layer[1].dispersion_coefficient` (constant model) | P Eq. 5 | derived per design |
| K_z = 0.09 W/(m K) | `layer[1].heat_dispersion_coefficient` | P Table 3 | |
| C_p,g = C_p,a = 30.7 J/(mol K), constant | `fluid.gas_heat_capacity_mode` = constant, `fluid.gas_heat_capacity[*]`, `fluid.adsorbed_heat_capacity[*]` | P Table 3 | |
| C_p,s = 1070 J/(kg K) | `adsorbent[1].heat_capacity` | P Table 3 | |
| μ = 1.72 × 10⁻⁵ Pa s, constant | `fluid.viscosity_mode` = constant, `fluid.viscosity` | P Table 3 | |
| h_in 8.6, h_out 2.5 W/(m² K) | `layer[1].internal_htc`, `layer[1].external_htc` | P Table 3 | |
| wall ρ 7800 kg/m³, C_p 502 J/(kg K), k 16 W/(m K) | `wall.density`, `wall.heat_capacity`, `wall.thermal_conductivity` | P Table 3 | |
| wall ends at T_a | `wall.boundary` = fixed, `wall.fixed_temperature` = 298.15 K | S Eq. 18 | |
| feed 15 % CO₂ / 85 % N₂, 1 bar, 298.15 K | stream `y`, temperatures, `bed.ambient_temperature` | P Sect. 2.1, Table 3 | |
| 13X isotherm: CO₂ q_b 3.09, q_d 2.54 mol/kg, b₀ 8.65 × 10⁻⁷, d₀ 2.63 × 10⁻⁸ m³/mol, ΔU_b −36 600, ΔU_d −35 700 J/mol; N₂ q_b 5.84 mol/kg, b₀ 2.50 × 10⁻⁶ m³/mol, ΔU_b −15 800 J/mol | `adsorbent[1].isotherm[C].m1/b1/q1/m2/b2/q2` | P Table 1 | extended dual-site Langmuir in concentration form (P Eqs. 1–4), `extended_dual_site_langmuir_c`; the ΔU are the rounded values of Haghpanah's Table 2 |
| ZIF-36-FRL isotherm: CO₂ q_b 2.01, q_d 1.01 mol/kg, b₀ 5.00 × 10⁻⁸, d₀ 9.10 × 10⁻⁷ m³/mol, ΔU_b −38 100, ΔU_d −21 900 J/mol; N₂ q_b 3.02 mol/kg, b₀ 6.29 × 10⁻⁷ m³/mol, ΔU_b −16 230 J/mol | the same | P Table 1 | the deck's adsorbent is inline; ZIF-36-FRL is not in the adsorbent library |
| LDF k_i = (c_i / q_i\*) 15 ε_p D_m / (τ_p r_p²) | `adsorbent[1].kinetics[C]` = `macropore_capacity` | P Eq. 6 | q\* per particle volume (ρ_s q\*), the dimensionally consistent reading, as the Haghpanah decks |
| heat of adsorption −ΔU_b,i + R T_F: CO₂ 39 079 (13X), 40 579 (ZIF) J/mol; N₂ 18 279 (13X), 18 709 (ZIF) J/mol | `adsorbent[1].heat_of_adsorption[C]` | S Eq. 8 (σ_i carries −ΔU_b,i + R T₀ T̄) | site b only, as the paper; R T at T_F instead of the local T (≤ 0.1 kJ/mol over the run) |
| heat-of-adsorption convention: internal energy | `fluid.heat_of_adsorption_convention` | P Table 2: the source is Σ (σ_i − Ω₃ T̄) ∂x_i/∂τ, Ω₃ ∝ C_p,a | |
| pressure drop: Darcy, 150 (1 − ε)² / ε² | `layer[1].momentum_model` = blake_kozeny | P Table 2 | |
| exponential pressure profiles, λ = 0.5 1/s | step boundary pressure: `double_exponential`, both rates 0.5 1/s, weight 1 | S Eqs. 24, 30, 36; P Table 3 | from the bed pressure at the start of the step |
| **adsorption at p_H** | `step[2].top` pressure = p_H | P Sect. 2.3; S Eq. 13 | S Eq. 13 writes p̄ = 1 with p₀ = p_F (S Sect. S.1), but S Eqs. 24, 30, 36 normalise by p_H: the text (adsorption at p_H) is used |
| interstitial feed velocity v_F | stream velocity (superficial) = ε v_F | S Eq. 45 (ε π r_in² v₀) | |
| compression and evacuation: isentropic, γ = 1.4, from p_F / to p_atm = 1 bar | stream and boundary `pump`; energy model isothermal (the gas enters the bed at T_F) | S Eqs. 45–48 | no work when the pressure ratio is below 1 (S Eqs. 46, 48) |
| η(p) = 0.8 × 19.55 p / (1 + 19.55 p), p in bar, every pump | pump efficiency `vacuum_pump` (eta_max 0.8, alpha 19.55) | P Eq. 11, S Eq. 49; P Table 3 (η = η(p)); Ward & Pini 2024 nomenclature ("adiabatic vacuum/compression efficiency") | evaluated at the suction pressure: the bed for the vacuum pumps, p_F for the compressors (η = 0.761) |
| initial bed: N₂ at p_H, 298.15 K | `initial_conditions.*` | P p. 13653 | |
| grid: WENO, 10 cells | `layer[1].n_cells` | P Sect. 2.2 (N = 10) | see "Grid" |

## Model-form audit

| Item | Paper | Deck |
|---|---|---|
| Isotherm and competition | extended DSL, shared denominators per site (P Eq. 1) | same |
| Rate model | LDF, k from P Eq. 6 with c_i/q_i\* | `macropore_capacity` |
| Gas holdup | ε only (S Eq. 1) | ε_t = ε_b |
| Axial dispersion, conduction | constant D_L (P Eq. 5), constant K_z | same |
| Wall | wall energy balance, ends at T_a | same |
| Pressure drop | Darcy, 150 | same |
| Energy balance | P Table 2 (Haghpanah's form: no gas accumulation term, Ω₂ ∂p̄/∂τ) | the engine's column energy balance (see the Haghpanah README: an isothermal run moves the KPIs by ≤ 0.15 pp) |
| KPIs | P Eqs. 7–10 (S Eqs. 37–48, from column inventories and boundary integrals) | the engine's boundary-flux KPIs: purity and recovery of the evacuation product, CO₂ fed in pressurisation and adsorption, energy = time-integrated pump power per kg of CO₂ product (J/kg ÷ 3600 = kWh/t), productivity = CO₂ product per π r_in² L per cycle time |
| CSS | relative change of purity, recovery, productivity and energy < 0.5 % over 10 cycles | Mass Balance (Per-Component), tighter |

## Measured errors (2026-11-15)

Purity and recovery in percentage points, productivity and energy relative.
Bounds are regression bounds set from the measured error; a metric is asserted
where it is within the 5 % target (5 pp), and recorded without a bound (–)
otherwise.

| Case | Purity paper / engine | Recovery paper / engine | Productivity [mol/(m³ s)] paper / engine | Energy [kWh/t] paper / engine |
|---|---|---|---|---|
| 13X_PL002 | 96.4 / 96.77 (0.37 pp, bound 0.45) | 90.2 / 84.66 (5.54 pp, –) | 1.64 / 1.522 (7.2 %, –) | 647 / 642.5 (0.7 %, bound 1.2 %) |
| 13X_PL003 | 95.8 / 95.85 (0.05 pp, 0.1) | 89.3 / 87.04 (2.26 pp, 2.6) | 1.13 / 1.092 (3.4 %, 4 %) | 599 / 555.5 (7.3 %, –) |
| 13X_PL005 | 96.5 / 96.24 (0.26 pp, 0.35) | 90.4 / 89.63 (0.77 pp, 0.9) | 1.27 / 1.256 (1.1 %, 1.5 %) | 693 / 573.0 (17.3 %, –) |
| 13X_PL007 | 95.9 / 95.48 (0.42 pp, 0.5) | 92.0 / 90.22 (1.78 pp, 2.0) | 2.14 / 2.097 (2.0 %, 2.5 %) | 865 / 547.6 (36.7 %, –) |
| 13X_PL010 | 94.9 / 94.50 (0.40 pp, 0.5) | 89.9 / 89.52 (0.38 pp, 0.45) | 0.994 / 0.9975 (0.4 %, 0.8 %) | 666 / 500.0 (24.9 %, –) |
| ZIF_PL002 | 96.6 / 96.38 (0.22 pp, 0.3) | 89.5 / 89.16 (0.34 pp, 0.45) | 0.782 / 0.7852 (0.4 %, 0.8 %) | 641 / 596.5 (6.9 %, –) |
| ZIF_PL003 | 95.8 / 95.64 (0.16 pp, 0.25) | 89.9 / 89.14 (0.76 pp, 0.9) | 0.586 / 0.5830 (0.5 %, 0.9 %) | 523 / 455.4 (12.9 %, –) |
| ZIF_PL005 | 95.6 / 95.07 (0.53 pp, 0.65) | 90.5 / 90.39 (0.11 pp, 0.2) | 0.811 / 0.8090 (0.2 %, 0.6 %) | 562 / 422.2 (24.9 %, –) |
| ZIF_PL007 | 94.2 / 93.82 (0.38 pp, 0.45) | 90.3 / 90.47 (0.17 pp, 0.25) | 0.750 / 0.7604 (1.4 %, 1.8 %) | 531 / 417.1 (21.4 %, –) |
| ZIF_PL010 | 96.3 / 95.69 (0.61 pp, 0.75) | 89.3 / 90.54 (1.24 pp, 1.45) | 1.01 / 1.0298 (2.0 %, 2.4 %) | 692 / 514.5 (25.6 %, –) |

The paper's conclusion holds in the deck: at the two P Table 7 designs, the
ZIF-36-FRL process uses less energy than the 13X process (417 against 556
kWh/t; paper 531 against 599) and has the lower productivity (0.76 against 1.09
mol/(m³ s); paper 0.750 against 1.13).

### Energy by step (engine, kWh/t)

| Case | Pressurisation | Adsorption | Blowdown | Evacuation |
|---|---|---|---|---|
| 13X_PL002 | 7.1 | 252.2 | 0.0 | 383.2 |
| 13X_PL003 | 13.9 | 287.7 | 0.1 | 253.9 |
| 13X_PL005 | 25.8 | 405.0 | 0.0 | 142.2 |
| 13X_PL007 | 25.1 | 424.5 | 0.0 | 97.9 |
| 13X_PL010 | 28.8 | 398.5 | 0.0 | 72.7 |
| ZIF_PL002 | 7.0 | 228.6 | 0.1 | 360.9 |
| ZIF_PL003 | 7.9 | 196.7 | 0.7 | 250.1 |
| ZIF_PL005 | 11.6 | 272.7 | 0.0 | 137.9 |
| ZIF_PL007 | 18.5 | 300.3 | 0.0 | 98.3 |
| ZIF_PL010 | 29.2 | 419.8 | 0.0 | 65.4 |

The blowdown costs nearly nothing: it ends at p_I ≥ 0.48 bar, and no work is
done above 1 bar (S Eq. 46).

### Grid

10 against 30 cells (the deck, once, 2026-11-15):

| Case | Recovery 10 / 30 cells | Productivity 10 / 30 | Energy 10 / 30 |
|---|---|---|---|
| 13X_PL002 | 84.66 / 85.55 % | 1.522 / 1.540 | 642.5 / 637.8 |
| 13X_PL003 | 87.04 / 88.48 % | 1.092 / 1.115 | 555.5 / 546.9 |
| 13X_PL007 | 90.22 / 90.64 % | 2.097 / 2.108 | 547.6 / 544.5 |
| ZIF_PL002 | 89.16 / 90.69 % | 0.7852 / 0.8000 | 596.5 / 583.0 |
| ZIF_PL007 | 90.47 / 91.08 % | 0.7604 / 0.7663 | 417.1 / 413.3 |

Over all ten designs, purity moves by ≤ 0.3 pp, energy by ≤ 2.3 %. Recovery rises by
0.4–1.5 pp at 30 cells: 10 cells is not grid-converged, but it is the paper's
grid, and the comparison is code to code. The deck keeps 10 cells.

## Attribution

One-factor runs on all ten designs (the deck, 10 cells).

| Change | Purity | Recovery | Productivity | Energy |
|---|---|---|---|---|
| compressor efficiency 0.72 (Haghpanah) instead of η(1 bar) = 0.761 | 0 | 0 | 0 | +2.3 to +5.0 % |
| isosteric heat −ΔU + R T (`isosteric_rt`) instead of the internal-energy form | +0.4 to +0.7 pp | −1.0 to −2.2 pp | −1.4 to −2.7 % | −0.4 to +2.1 % |
| D_m in the rate × 10 | −1.0 to +0.04 pp | +0.4 to +6.2 pp | +0.4 to +8.6 % | −2.8 to −18.7 % |
| 30 cells | −0.3 to 0 pp | +0.4 to +1.5 pp | +0.5 to +2.1 % | −0.5 to −2.3 % |

### Recovery of 13X at p_L = 0.02 and 0.03 bar

The deck is 5.5 and 2.3 pp low (with productivity 7 % and 3 % low). The
ZIF-36-FRL designs at the same pressures agree within 0.8 pp. A ten times faster
CO₂ rate gives 90.8 % and 90.2 % (paper 90.2 % and 89.3 %): the shortfall is CO₂
that desorbs too slowly in the deep vacuum. The rate coefficient of P Eq. 6 is
proportional to c_i / q_i\*, which is smallest for the steep 13X isotherm at low
pressure. The paper does not state the basis of q\* in Eq. 6 (the same open
question as Haghpanah's Eq. 40, see the Haghpanah README), so the deck keeps
the dimensionally consistent reading.

### Specific energy

The energy agrees at p_L = 0.02 bar (0.7 % for 13X, 6.9 % for ZIF-36-FRL) and
for 13X at 0.03 bar (7.3 %). In the six designs where feed compression is most
of the work, the deck is 13–37 % low. Checks:

- The engine's feed-compression work matches the printed formula (S Eq. 45).
  By hand, for 13X_PL007: ε π r_in² v_F p_H / (R T) t_ads moles at
  γ/(γ − 1) R T [(9.01)^0.286 − 1] / 0.761 per mole give 420 kWh/t; the engine
  logs 424.5.
- The one-factor runs above move the energy by at most 5 %, except the faster
  rate, which lowers it.
- The paper's number cannot be matched by one factor on the compression work
  (it would be ×1.0 to ×1.7 between designs), nor on the evacuation work
  (×1.0 to ×4.2). A least-squares fit of both factors over the ten designs
  still leaves errors of up to 14 %.
- The paper does not give its energy by step, so the gap cannot be put on one
  step. P Table 7's electricity cost (599 kWh/t × 0.06 $/kWh = 35.94 $/t)
  confirms that the printed energies are the ones the paper used.

The cause is open. The energy is asserted only where it is within 5 % (13X_PL002).

### What is not modelled

- The paper's CSS rule (KPIs settled to 0.5 % over 10 cycles): the engine runs
  to its own, tighter CSS.
- P Table 2's energy-balance form: see the model-form audit.

## Validation protocol

1. **Source**: the paper and its SI were read in full and every input and
   reference value transcribed. The numbers here are the printed ones (no value
   read off a figure).
2. **Provenance**: above. Two assumptions: the wall of the scaled columns and
   the adsorption pressure (S Eq. 13 against the text).
3. **Model form**: above; gas holdup in the bed voids only, the paper's grid.
4. **Data**: no experimental curves (the references are model results).
5. **Metrics**: the paper's KPI definitions (P Eqs. 7–10); purity and recovery
   in pp, the rest relative.
6. **Attribution**: above. Two open disagreements: the 13X recovery at
   p_L ≤ 0.03 bar, and the energy where compression dominates.
7. **Report**: bounds from the measured errors; the validation page regenerated.
