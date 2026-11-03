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
- **Molar-basis gas heat capacity (2026-11).** The gas-mixture heat capacity used to mix the per-kilogram pure-component values by mole fraction, which is far off when the components' molar masses differ: here the bed starts in helium (4 g/mol) and is fed CO₂/CH₄ (44/16 g/mol), so a He-rich cell had several times its true heat capacity per mole. With the mass-consistent mixture the temperature fronts fit the thermocouples better in every case (pooled RMSE A3 14.0 → 11.4 K, A5 18.3 → 16.1 K, A7 11.1 → 9.3 K, B2 32.9 → 27.0 K) and the composition fronts arrive 0.5–1.6 % earlier, which moves the front-time and curve errors that were close to their bounds: A2 y_CO₂ NRMSE 0.338 → 0.350; A3 CH₄ front 3.1 → 3.8 %; A4 CO₂ front 8.7 → 10.0 %, CH₄ front 17.8 → 18.5 %; A5 y_CO₂ NRMSE 0.155 → 0.168, CO₂ front 5.2 → 5.7 %; A6 y_CO₂ NRMSE 0.102 → 0.157, CO₂ front 1.5 → 2.6 %. Those eight tolerances were re-measured (see `streb_2021.cpp`); A1, A7, B1 and B2 stay within theirs.

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
