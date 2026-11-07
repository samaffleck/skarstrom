# Streb et al. (2021) — CO₂/CH₄ separation breakthrough

## Citation

> Streb, A., Mazzotti, M. *Adsorption for efficient low carbon hydrogen
> production: part 1 — adsorption equilibrium and breakthrough studies for
> H₂/CO₂/CH₄ on zeolite 13X.* Adsorption 27:541–558 (2021).
> DOI: [10.1007/s10450-021-00306-y](https://doi.org/10.1007/s10450-021-00306-y)


## Files

| File | Purpose |
|------|---------|
| `streb_2021.cpp` | `Streb2021Test` fixture + reproduced breakthrough cases |

## Reproduced cases

| Test (`Streb2021Test.*`) | Paper reference | Notes |
|--------------------------|-----------------|-------|
| `ExperimentA1` | Table 2, A1 (25 °C, 1 bar, 20 cm³/s) | Binary CO₂/CH₄ 50:50, bundled example |
| `ExperimentA2` | A2 (25 °C, 10 bar, 10 cm³/s) | Binary |
| `ExperimentA3` | A3 (25 °C, 10 bar, 20 cm³/s) | Binary, reference experiment (Fig. 6), bundled example |
| `ExperimentA4` | A4 (25 °C, 10 bar, 40 cm³/s) | Binary |
| `ExperimentA5` | A5 (25 °C, 25 bar, 20 cm³/s) | Binary |
| `ExperimentA6` | A6 (45 °C, 10 bar, 20 cm³/s) | Binary |
| `ExperimentA7` | A7 (65 °C, 10 bar, 20 cm³/s) | Binary |
| `ExperimentB1` | B1 (25 °C, 1 bar, 20 cm³/s) | Ternary H₂:CO₂:CH₄ 75:20:5 |
| `ExperimentB2` | B2 (25 °C, 10 bar, 10 cm³/s) | Ternary |

## Provenance (WP10 audit, 2026-10-31, from the paper text)

| Deck value | Value | Source | Note |
|---|---|---|---|
| Column length, inner diameter | 1.2 m, 0.025 m | Table 1 | |
| Wall thickness | 2.5 mm | Table 1, d_o 0.03 m | |
| Wall heat capacity | 4 × 10⁶ J/(m³ K) as ρ_w C_p,w = 8000 × 500 | Table 1, C_w | the paper lists it per kg; the value is volumetric |
| Wall axial conduction | 0 | Sect. 3.2, assumption "axial conductivity of the column wall is neglected" | |
| Particle diameter | 1.8 mm | Table 1 (pellets 1.6–2 mm, Sect. 2.1) | |
| Particle density, skeletal density, bed density | 1085, 2359, 595 kg/m³ | Table 1 | ε_p = 1 − 1085/2359 = 0.540, ε_b = 1 − 595/1085 = 0.452 (derived) |
| Adsorbent C_p | 1100 J/(kg K) | Table 1 | |
| Adsorbent thermal conductivity | 0.15 W/(m K) | not in the paper | enters only the Dixon coefficient's static term |
| Isotherm | extended Sips (Eqs. 1–4, 6) | Table 3 (CO₂ from ref. [14]) | B_1 converted from 10⁻⁴ bar⁻¹ to Pa⁻¹, B_2 from kJ/mol to J/mol; T_ref 25 °C |
| Competition | **extended Sips** | Sect. 4.2.1: the k values in Table 4 were estimated with IAST; IAST and RAST "predict the time when CH₄ breakthrough occurs and the shape of the CO₂ breakthrough front better than extended isotherms" | the engine has no IAST for Sips |
| LDF coefficients | CO₂ 0.06, CH₄ 0.2, H₂ 1.0 s⁻¹ | Table 4 (a, b) | the ternary-fitted values (0.09, 0.4, 1.2; footnote c) were not used by the paper for prediction either |
| Heats of adsorption | 36.3, 18.8, 10.7 kJ/mol | Sect. 3.1.1 (temperature-averaged isosteric heats) | constant, as the paper assumes |
| Axial dispersion, axial conduction | 0, 0 | Sect. 3.2 assumptions ("diffusion is neglected"; Eq. 18 has no conduction term) | |
| Bed-to-wall HTC | Dixon (1996) | Eq. 21 and ref. [39] | h_L lumps λ_e and α_w |
| Wall-to-thermofluid HTC | 220 W/(m² K) | Sect. 3.2.1, h_w from ref. [23] | the column is double-jacketed and thermostated, so the ambient is the experiment temperature |
| Pressure drop | Ergun | Eq. 20 | |
| Feed | Table 2 flows "at P and T", 50:50 CO₂:CH₄ or 75:20:5 H₂:CO₂:CH₄ | Table 2, Sect. 2.1 | **the paper fits the feed velocity per experiment within the MFC's ±2 % of full scale (Sect. 3.2.1); the deck uses the set-point** |
| Initial gas | helium | Sect. 2.2 | |
| Thermocouples | z = 10, 35, 60, 85, 110 cm, centre line (r = 0) | Sect. 2.2 | the paper notes the centre line cools slower than its 1D model |
| Sampling line | **not modelled** | Table 1: 10.51 m × 6 mm downstream, modelled by the paper as plug flow at column pressure up to the BPR and at ambient after | about 0.3 L; a delay of order seconds to tens of seconds depending on pressure |

Bold entries differ from the paper's model.

## Validation

The tests compare the simulation with the paper using the helpers in
`../../fixtures/validation_metrics.h`; the published numbers they use are in
[`reference_values.csv`](reference_values.csv) (each with its table / figure /
page). Every tolerance is the measured error rounded up, with the measured
value in a comment next to it. Results of the last run:
[validation page](../../../../docs/validation.md#streb-2021)
(generated from `../validation_summary.json`).

- CO2 front positions agree with the digitised experiments to 0.5-9 %, CH4 fronts to 0.3-15 %. Curve NRMSE is large (13-83 % for CH4) because the fronts are nearly vertical: a 10-20 s shift on a 20 s front dominates the RMSE. Thermocouple RMSE is 11-33 K, driven by the timing and height of the 30-100 K temperature fronts at the centre line.
- adsim uses extended Sips competition and the MFC set-point velocity; the paper shows IAST/RAST predict the binary cases better than extended Sips and fits the feed velocity per experiment within the MFC uncertainty. The ~10 m of sampling line to the MS is not modelled.
- A1/A3 previously compared simulated CO2 against the digitised CH4 curve (alphabetical file order) with a bound of 0.8 in mole fraction; fixed.
- **Flow-inlet and inlet-energy fixes (2026-11).** On a Kozeny–Carman bed the flow boundary condition used to compute the boundary-face flow with the ghost cell's viscosity alone, while the face velocity that carries gas into the bed uses the mean of the ghost and first-cell values. The ghost holds the CO₂/CH₄ feed (μ ≈ 1.3 × 10⁻⁵ Pa s) and the first cell starts in helium (μ ≈ 2 × 10⁻⁵ Pa s) and is heated by the front, so the bed received a few per cent less gas than the set point the stream log records. The bed now takes in exactly the set flow, and the CH₄ fronts arrive about 3 % earlier: CH₄ front-time errors A2 7.6 → 10.5 %, A3 0.3 → 3.1 %, A5 12.1 → 14.2 %, A6 8.7 → 11.7 %, A7 7.3 → 10.2 %, B1 0.7 → 2.1 %; A3 y_CH₄ NRMSE 0.21 → 0.50, A6 y_CO₂ NRMSE 0.048 → 0.10, A2 thermocouple RMSE 14.1 → 18.6 K. The Danckwerts inlet change for the energy balance (no conduction across the inlet face, upstream value the feed temperature) moves these cases by less than 1 % once the flow is right, and the species change does not move them. The affected tolerances were re-measured (see `streb_2021.cpp`).
- **Constant gas heat capacity, as the paper assumes (2026-11).** The engine's default gas heat capacity used to be CoolProp's real-gas value of each pure component at the initial pressure and temperature, held for the run; it is now the ideal-gas C_p0(T) of the NASA TM-4513 polynomials at the local temperature (`docs/equation-of-state.md`). Streb et al. state in their model assumptions that constant, average values are used for the heat capacity of the gas phase, and print no value. So these decks use the constant gas heat capacity mode, with each component's ideal-gas C_p0 at the case's feed (= initial) temperature as the average: our choice, since the paper gives no number. This is the lowest value of the run (the bed heats by 15–40 K), and lower than the old held value, which at these pressures was the real-gas C_p of pure CO₂ at the total pressure (9 % above ideal at 10 bar, 29 % at 25 bar). Errors as old held value → default C_p0(T) → constant C_p0(T_feed), the last being what the tests now measure:
  - A5 CO₂ front 5.94 → 6.21 → 6.48 % (tolerance 0.061 → **0.066**), A5 y_CO₂ NRMSE 0.173 → 0.179 → 0.185 (0.18 → **0.19**);
  - A6 CO₂ front 2.73 → 2.68 → 2.92 % (0.028 → **0.030**), A6 y_CO₂ NRMSE 0.161 → 0.158 → 0.170 (0.17 → **0.175**);
  - the four tolerances were raised to the measured values rounded up, with the owner's approval; every other bound holds;
  - thermocouple RMSE improves in every case: A3 12.1 → 12.1 → 11.4 K, A5 16.2 → 14.9 → 14.3 K, A7 9.0 → 9.1 → 8.7 K; A3 and A7 CO₂ fronts improve (A7 0.49 → 0.57 → 0.31 %);
  - the 1 bar cases (A1, B1, B2) move by less than 0.1 %.
- **Molar-basis gas heat capacity (2026-11).** The gas-mixture heat capacity used to mix the per-kilogram pure-component values by mole fraction, which is far off when the components' molar masses differ: here the bed starts in helium (4 g/mol) and is fed CO₂/CH₄ (44/16 g/mol), so a He-rich cell had several times its true heat capacity per mole. With the mass-consistent mixture the temperature fronts fit the thermocouples better in every case (pooled RMSE A3 14.0 → 11.4 K, A5 18.3 → 16.1 K, A7 11.1 → 9.3 K, B2 32.9 → 27.0 K) and the composition fronts arrive 0.5–1.6 % earlier, which moves the front-time and curve errors that were close to their bounds: A2 y_CO₂ NRMSE 0.338 → 0.350; A3 CH₄ front 3.1 → 3.8 %; A4 CO₂ front 8.7 → 10.0 %, CH₄ front 17.8 → 18.5 %; A5 y_CO₂ NRMSE 0.155 → 0.168, CO₂ front 5.2 → 5.7 %; A6 y_CO₂ NRMSE 0.102 → 0.157, CO₂ front 1.5 → 2.6 %. Those eight tolerances were re-measured (see `streb_2021.cpp`); A1, A7, B1 and B2 stay within theirs.
- **Competitive-Sips trace fix (WP19-1, 2026-11).** `extended_sips_2` used to multiply each competitive term by y²/(y² + 10⁻⁶), which halved a component's loading at y = 10⁻³ and cut it by 1 % at 10⁻²; it now uses the published isotherm down to y = 10⁻⁵ (see `docs/isotherms.md`). Only the toes of the fronts pass through those mole fractions, so the metrics barely move: the largest changes are the ternary CH₄ fronts (B1 front-time error 2.13 → 2.38 %, B2 13.7 → 14.0 %, B2 y_CH₄ NRMSE 0.319 → 0.320); every other metric moves by less than 0.002 in NRMSE, 0.1 % in front time and 0.03 K in thermocouple RMSE. No tolerance changed.

### Grid and advection scheme (2026-11-04)

The decks ran first-order upwind advection on 60 cells. A grid study of A1
and A3 (10-160 cells, Richardson extrapolation;
[Numerics](../../../../docs/numerics.md#advection-schemes)) found the CO2
front converged on 20-40 cells, while the CH4 half-feed time needs 120-170
upwind or 110-150 WENO cells to be within 0.5 % of its grid limit (the
thermal front and the light component's dispersion keep first-order
behaviour, so the schemes differ little here). All cases now run WENO, the
engine's default, on 160 cells.

Before: upwind, 60 cells; after: WENO, 160 cells (bold: changed bound). The
CH4 fronts come 0.6-2.7 % later and the CO2 fronts move by less than 0.2 %; the
outlet curves are sharper, so the y_CH4 NRMSE of A3 (0.56 → 0.60) and the
y_CO2 NRMSE of A4 (0.287 → 0.291) passed their bounds (0.58 → 0.73, 0.29 →
0.35), and B1's CH4 front time came to within 0.2 % (bound 0.024 → 0.004).
Both runs include the competitive-Sips trace-loading fix.

| Case | Quantity | Reference | Before | After | Error before → after | Bound before → after |
|---|---|---|---|---|---|---|
| ExperimentA1 | outlet y_CO2 | 0.5 | 0.2191 | 0.2237 | 0.4382 → 0.4474 | 0.53 |
| ExperimentA1 | outlet y_CH4 | 0.5 | 0.2274 | 0.2621 | 0.4548 → 0.5243 | 0.83 |
| ExperimentA1 | CO2 front time (y_CO2 = y_feed/2) | 3662.94 | 3932.80 | 3929.81 | 0.07367 → 0.07286 | 0.098 |
| ExperimentA1 | CH4 front time (y_CH4 = y_feed/2) | 458.57 | 470.74 | 475.72 | 0.02653 → 0.03741 | 0.12 |
| ExperimentA1 | thermocouple T (5 positions, pooled) | 0 | 19.47 | 20.30 | 19.47 → 20.3 | 25 |
| ExperimentA2 | outlet y_CO2 | 0.5 | 0.1753 | 0.1784 | 0.3505 → 0.3569 | 0.37 |
| ExperimentA2 | outlet y_CH4 | 0.5 | 0.4195 | 0.4202 | 0.839 → 0.8404 | 0.96 |
| ExperimentA2 | CO2 front time (y_CO2 = y_feed/2) | 965.66 | 912.83 | 912.28 | 0.05471 → 0.05528 | 0.058 |
| ExperimentA2 | CH4 front time (y_CH4 = y_feed/2) | 364.98 | 324.67 | 326.97 | 0.1104 → 0.1041 | 0.13 |
| ExperimentA2 | thermocouple T (5 positions, pooled) | 0 | 18.59 | 19.88 | 18.59 → 19.88 | 22 |
| ExperimentA3 | outlet y_CO2 | 0.5 | 0.08861 | 0.09217 | 0.1772 → 0.1843 | 0.27 |
| ExperimentA3 | outlet y_CH4 | 0.5 | 0.2796 | 0.3016 | 0.5591 → 0.6031 | **0.58 → 0.73** |
| ExperimentA3 | CO2 front time (y_CO2 = y_feed/2) | 414.86 | 434.75 | 434.30 | 0.04795 → 0.04686 | 0.08 |
| ExperimentA3 | CH4 front time (y_CH4 = y_feed/2) | 156.53 | 150.53 | 151.67 | 0.03832 → 0.031 | 0.041 |
| ExperimentA3 | thermocouple T (5 positions, pooled) | 0 | 11.42 | 12.07 | 11.42 → 12.07 | 19 |
| ExperimentA3 | CH4 breakthrough time (y > 0.01) | 140.00 | 145.99 | 150.11 | 5.99 → 10.11 | – |
| ExperimentA3 | CO2 breakthrough time (y > 0.01) | 400.00 | 416.19 | 420.64 | 16.19 → 20.64 | – |
| ExperimentA4 | outlet y_CO2 | 0.5 | 0.1432 | 0.1453 | 0.2865 → 0.2906 | **0.29 → 0.35** |
| ExperimentA4 | outlet y_CH4 | 0.5 | 0.396 | 0.3964 | 0.7921 → 0.7928 | 0.91 |
| ExperimentA4 | CO2 front time (y_CO2 = y_feed/2) | 226.18 | 203.61 | 203.36 | 0.09977 → 0.1009 | 0.105 |
| ExperimentA4 | CH4 front time (y_CH4 = y_feed/2) | 85.66 | 69.83 | 70.42 | 0.1849 → 0.1779 | 0.195 |
| ExperimentA4 | thermocouple T (5 positions, pooled) | 0 | 16.39 | 16.68 | 16.39 → 16.68 | 20 |
| ExperimentA5 | outlet y_CO2 | 0.5 | 0.08388 | 0.08582 | 0.1678 → 0.1716 | 0.18 |
| ExperimentA5 | outlet y_CH4 | 0.5 | 0.3299 | 0.3302 | 0.6598 → 0.6604 | 0.76 |
| ExperimentA5 | CO2 front time (y_CO2 = y_feed/2) | 189.45 | 178.59 | 178.29 | 0.05737 → 0.05893 | 0.061 |
| ExperimentA5 | CH4 front time (y_CH4 = y_feed/2) | 94.98 | 80.90 | 81.39 | 0.1482 → 0.1431 | 0.17 |
| ExperimentA5 | thermocouple T (5 positions, pooled) | 0 | 16.08 | 16.08 | 16.08 → 16.08 | 21 |
| ExperimentA6 | outlet y_CO2 | 0.5 | 0.07849 | 0.07962 | 0.157 → 0.1592 | 0.17 |
| ExperimentA6 | outlet y_CH4 | 0.5 | 0.375 | 0.3751 | 0.75 → 0.7503 | 0.85 |
| ExperimentA6 | CO2 front time (y_CO2 = y_feed/2) | 447.88 | 436.28 | 435.78 | 0.02588 → 0.02701 | 0.028 |
| ExperimentA6 | CH4 front time (y_CH4 = y_feed/2) | 160.57 | 140.85 | 142.13 | 0.1228 → 0.1149 | 0.14 |
| ExperimentA6 | thermocouple T (5 positions, pooled) | 0 | 11.76 | 11.76 | 11.76 → 11.76 | 13 |
| ExperimentA7 | outlet y_CO2 | 0.5 | 0.01408 | 0.01714 | 0.02815 → 0.03428 | 0.14 |
| ExperimentA7 | outlet y_CH4 | 0.5 | 0.3451 | 0.3458 | 0.6902 → 0.6917 | 0.72 |
| ExperimentA7 | CO2 front time (y_CO2 = y_feed/2) | 426.45 | 429.25 | 428.67 | 0.006571 → 0.005207 | 0.04 |
| ExperimentA7 | CH4 front time (y_CH4 = y_feed/2) | 146.22 | 130.47 | 131.84 | 0.1077 → 0.09834 | 0.12 |
| ExperimentA7 | thermocouple T (5 positions, pooled) | 0 | 9.25 | 9.14 | 9.25 → 9.143 | 14 |
| ExperimentB1 | outlet y_CO2 | 0.2 | 0.07541 | 0.07535 | 0.377 → 0.3768 | 0.45 |
| ExperimentB1 | outlet y_CH4 | 0.05 | 0.006493 | 0.00632 | 0.1299 → 0.1264 | 0.15 |
| ExperimentB1 | CO2 front time (y_CO2 = y_feed/2) | 8474.18 | 9020.71 | 9018.64 | 0.06449 → 0.06425 | 0.081 |
| ExperimentB1 | CH4 front time (y_CH4 = y_feed/2) | 337.67 | 329.62 | 337.01 | 0.02382 → 0.001939 | **0.024 → 0.004** |
| ExperimentB1 | outlet y_H2 | 0.75 | 0.2651 | 0.2793 | 0.3535 → 0.3724 | 0.42 |
| ExperimentB1 | thermocouple T (5 positions, pooled) | 0 | 16.11 | 17.69 | 16.11 → 17.69 | 20 |
| ExperimentB2 | outlet y_CO2 | 0.2 | 0.08957 | 0.09007 | 0.4478 → 0.4503 | 0.57 |
| ExperimentB2 | outlet y_CH4 | 0.05 | 0.01602 | 0.01645 | 0.3204 → 0.329 | 0.35 |
| ExperimentB2 | CO2 front time (y_CO2 = y_feed/2) | 2694.28 | 2405.45 | 2405.20 | 0.1072 → 0.1073 | 0.11 |
| ExperimentB2 | CH4 front time (y_CH4 = y_feed/2) | 713.25 | 613.58 | 630.36 | 0.1397 → 0.1162 | 0.15 |
| ExperimentB2 | outlet y_H2 | 0.75 | 0.2983 | 0.2984 | 0.3978 → 0.3979 | 0.46 |
| ExperimentB2 | thermocouple T (5 positions, pooled) | 0 | 27.01 | 28.27 | 27.01 → 28.27 | 38 |

### Attribution runs (A3, 2026-10-31)

One change at a time on the bundled A3 deck, front times at half the feed
fraction (y = 0.25) against the digitised Fig. 6:

| Variant | CO₂ NRMSE | CO₂ front sim / exp [s] | CH₄ NRMSE | CH₄ front sim / exp [s] |
|---|---|---|---|---|
| deck (set-point velocity, no line) | 0.229 | 444 / 415 | 0.212 | 157 / 157 |
| feed velocity −2 % | 0.264 | 453 / 415 | 0.306 | 160 / 157 |
| feed velocity +2 % | 0.174 | 434 / 415 | 0.349 | 154 / 157 |
| feed velocity +5 % | 0.064 | 421 / 415 | 0.588 | 149 / 157 |
| 0.3 L top void (the sampling line at column pressure) | 0.272 | 462 / 415 | 0.549 | 168 / 157 |

A velocity shift moves the two fronts together, so it cannot close the CO₂
error (7 %) without opening a CH₄ error: with extended Sips the CO₂ front is
late *relative to* the CH₄ front, which is what the paper reports for extended
isotherms against IAST/RAST (Sect. 4.2.1). The sampling line as a void at
column pressure delays both fronts by 3–7 % and is not what the paper models
(plug flow at column pressure to the BPR, ambient after it). The remaining
disagreement is therefore attributed to the competition model, not to the
column model; an IAST option for the Sips isotherm is the open item (WP10-E3).
