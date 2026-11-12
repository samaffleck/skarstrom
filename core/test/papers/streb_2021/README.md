# Streb et al. (2021) — CO₂/CH₄ separation breakthrough

## Citation

> Streb, A., Mazzotti, M. *Adsorption for efficient low carbon hydrogen
> production: part 1 — adsorption equilibrium and breakthrough studies for
> H₂/CO₂/CH₄ on zeolite 13X.* Adsorption 27:541–558 (2021).
> DOI: [10.1007/s10450-021-00306-y](https://doi.org/10.1007/s10450-021-00306-y)


## Files

| File | Purpose |
|------|---------|
| `streb_2021.cpp` | `Streb2021Test` fixture + reproduced breakthrough cases, the IAST variant of each (recorded only) and the IAST pure-limit check |
| `streb_calibration.py` | The A3 calibration: flow factors from the thermocouples, the least-squares fit and the predictions ([Calibration on A3](#calibration-on-a3-2026-11)) |

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
| Competition | **extended Sips** | Sect. 4.2.1: the k values in Table 4 were estimated with IAST; IAST and RAST "predict the time when CH₄ breakthrough occurs and the shape of the CO₂ breakthrough front better than extended isotherms" | IAST over the same pure Sips terms is recorded for every case (`*_IAST_RecordOnly`); it fits worse overall, see "IAST versus extended Sips" |
| LDF coefficients | **CO₂ 0.0592, CH₄ 0.137 s⁻¹ (calibrated on A3)**; H₂ 1.0 s⁻¹ | Table 4 (a, b): CO₂ 0.06, CH₄ 0.2 s⁻¹, themselves fitted to A3 | see [Calibration on A3](#calibration-on-a3-2026-11); the ternary-fitted values (0.09, 0.4, 1.2; footnote c) were not used by the paper for prediction either |
| Heats of adsorption | 36.3, 18.8, 10.7 kJ/mol | Sect. 3.1.1 (temperature-averaged isosteric heats) | constant, as the paper assumes |
| Axial dispersion, axial conduction | 0, 0 | Sect. 3.2 assumptions ("diffusion is neglected"; Eq. 18 has no conduction term) | |
| Bed-to-wall HTC | Dixon (1996) | Eq. 21 and ref. [39] | h_L lumps λ_e and α_w |
| Wall-to-thermofluid HTC | **49.7 W/(m² K) (calibrated on A3)** | Sect. 3.2.1: 220 W/(m² K), h_w from ref. [23] | the column is double-jacketed and thermostated, so the ambient is the experiment temperature |
| Pressure drop | Ergun | Eq. 20 | |
| Feed | Table 2 flows "at P and T" × a flow factor per case, 50:50 CO₂:CH₄ or 75:20:5 H₂:CO₂:CH₄ | Table 2, Sect. 2.1 | the paper fits the feed velocity per experiment within the MFC's ±2 % of full scale (Sect. 3.2.1); **the flow factors (0.887–1.057) are fitted to each case's thermocouples** ([Calibration on A3](#calibration-on-a3-2026-11)) |
| Initial gas | helium | Sect. 2.2 | |
| Thermocouples | z = 10, 35, 60, 85, 110 cm, centre line (r = 0) | Sect. 2.2 | the paper notes the centre line cools slower than its 1D model |
| Sampling line | **0.142 L of gas at column pressure, plug flow (calibrated on A3)**: an inert layer above the bed, bed voidage 0.99, 40 cells, adiabatic, no dispersion | Table 1: 10.51 m × 6 mm downstream, modelled by the paper as plug flow at column pressure up to the BPR and at ambient after | 0.13 L for a 4 mm bore, 0.30 L for 6 mm. **The deck keeps the whole line at column pressure and the gas at the bed-outlet temperature.** The engine's void region is one well-mixed cell and smears the CH₄ front; the engine has no plug-flow void option yet |

Bold entries differ from the paper's model.

## Validation

The tests compare the simulation with the paper using the helpers in
`../../fixtures/validation_metrics.h`; the published numbers they use are in
[`reference_values.csv`](reference_values.csv) (each with its table / figure /
page). Every tolerance is the measured error rounded up, with the measured
value in a comment next to it. Results of the last run:
[validation page](../../../../docs/validation.md#streb-2021)
(generated from `../validation_summary.json`).

- With the decks calibrated on A3 ([Calibration on A3](#calibration-on-a3-2026-11)), the CO₂ front (half the feed fraction) is within 1.4 % of the measurement in every case and the CH₄ front within 2.4 %, except B1 (4.7 % early). The outlet NRMSE is 1.3–16 % for CO₂ and 4.7–31 % for CH₄, and the thermocouple RMSE 5.5–11.4 K (before: 0.3–11 %, 0.2–18 %, 2.4–45 %, 13–84 % and 8.7–28.5 K).
- Only A3 informs the calibrated values. The other eight cases are predictions, except their feed flows: each case's flow factor is fitted to the speed of its own CO₂ temperature front, which uses none of the outlet data. Six cases need more than 2 % of the set point (up to −11.3 % in B2). All but A4 are within the paper's ±2 % of the MFC full scale if the full scale is at least 31 NL/min; A4 needs 102 NL/min.
- The decks keep extended Sips competition; the paper shows IAST/RAST predict the binary cases better than extended Sips ([IAST versus extended Sips](#iast-versus-extended-sips-2026-11-10), measured before the calibration).
- The sections below record the history of the decks, oldest first. Their numbers are the results at the time.
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

### Calibration on A3 (2026-11)

Before this calibration the decks used the paper's values and the Table 2 set-point flows. They had no sampling line. The CO₂ fronts were 0.3–11 % off and the CH₄ fronts 0.2–18 %, and the thermocouple RMSE was 8.7–28.5 K.

**The feed flow differs from the set point from run to run.** The thermocouples see the CO₂ front inside the column, so the time it takes to pass from 35 to 110 cm measures the front speed. The sampling line does not affect it. A2, A3 and A4 differ only in their flow (10, 20 and 40 cm³/s at 10 bar and 25 °C), so they have the same equilibrium. But their simulated fronts were 4.6 % fast, 6.8 % slow and 8.2 % fast relative to the measured ones. No model parameter can move these three cases in opposite directions; a flow error can. The paper fits the feed velocity per experiment for the same reason (Sect. 3.2.1).

**The rest of the outlet error is the sampling line.** With the flow corrected, the outlet CO₂ fronts at 10–25 bar came 1–3 % early and the CH₄ fronts 6–16 % early. A plug-flow line delays the CH₄ front by more than the CO₂ front, because the outlet flow before CO₂ breakthrough is about half the feed. That is what the residuals need.

**Method** ([`streb_calibration.py`](streb_calibration.py) repeats each stage).

1. The A3 deck was calibrated by weighted least squares (scipy `least_squares`). The residuals were the digitised outlet y_CO₂ and y_CH₄ (weight 1/0.03) and the five thermocouples (weight 1/(3 K)). The fitted values were the CO₂ and CH₄ LDF coefficients, the wall-to-thermofluid coefficient h_w and the gas volume of the sampling line. A3's flow factor was held at 1.02, the +2 % bound: its thermocouples ask for more (1.033 with the calibrated values).
2. Each other case's flow factor was then set so that the simulated CO₂ temperature front crosses the thermocouples at 35, 60, 85 and 110 cm at the measured speed: the slope of front-arrival time against position, where the front arrives when the temperature reaches half its rise. The other values are A3's.
3. Every case was then run and scored with the test's metrics.

The sampling line is an inert layer above the bed (bed voidage 0.99, 40 cells, no dispersion, no heat exchange), so it carries the outlet gas as plug flow. The engine's void region is a single well-mixed cell. On A3, at an earlier step of the fit (k_CO₂ 0.084, k_CH₄ 0.39 s⁻¹, h_w 108 W/(m² K), line 0.18 L, flow factor 1.067), a well-mixed top void of 0.18 L delays the CO₂ front correctly (−0.7 %) but leaves the CH₄ front 5.2 % early with an NRMSE of 0.29. A plug-flow delay of the same volume gives −0.5 %, +0.9 % and 0.12. A well-mixed tank reaches half the step at 0.69 of its residence time and adds an exponential tail, so it cannot match the sharp CH₄ front.

| Value | Calibrated (± 1 σ of the fit) | Paper |
|---|---|---|
| LDF coefficient, CO₂ | 0.0592 s⁻¹ (± 4.7 %) | 0.06 s⁻¹ (Table 4, fitted to A3) |
| LDF coefficient, CH₄ | 0.137 s⁻¹ (± 13 %) | 0.2 s⁻¹ (Table 4, fitted to A3) |
| Wall-to-thermofluid h_w | 49.7 W/(m² K) (± 6.7 %) | 220 W/(m² K) (Sect. 3.2.1, from ref. [23]) |
| Sampling-line gas volume | 0.142 ± 0.003 L at column pressure | 10.51 m × 6 mm (Table 1): 0.13 L for a 4 mm bore, 0.30 L for 6 mm |
| A3 flow factor | 1.020 (on the bound) | the paper fits the velocity per experiment within ±2 % of the MFC full scale |

The low h_w is what brings the thermocouples closer. The 1D model cools the bed behind the front faster than the centre-line thermocouples do (the paper notes the same), and a smaller h_w slows the cooling. h_w is therefore an effective value for the centre line, not the jacket coefficient. The fitted LDF coefficients are close to the paper's, which were fitted to the same experiment.

**Flow factors.** Set point (normal litres per minute at 0 °C, 1.01325 bar), the factor from the thermocouples, and the full scale at which the deviation is 2 % of it:

| Case | Set point [NL/min] | Flow factor | Deviation [NL/min] | Full scale for 2 % [NL/min] |
|---|---|---|---|---|
| A1 | 1.09 | 1.051 (+5.1 %) | +0.06 | 2.8 |
| A2 | 5.43 | 0.939 (−6.1 %) | −0.33 | 16.6 |
| A3 | 10.86 | 1.020 (+2.0 %, calibration bound; the thermocouples give 1.033) | +0.22 | 10.9 |
| A4 | 21.71 | 0.906 (−9.4 %) | −2.04 | 102 |
| A5 | 27.14 | 0.978 (−2.2 %) | −0.60 | 29.9 |
| A6 | 10.17 | 0.965 (−3.5 %) | −0.36 | 17.8 |
| A7 | 9.57 | 0.995 (−0.5 %) | −0.05 | 2.4 |
| B1 | 1.09 | 1.057 (+5.7 %) | +0.06 | 3.1 |
| B2 | 5.43 | 0.887 (−11.3 %) | −0.61 | 30.7 |

The paper's bound is ±2 % of the MFC full scale, and the full scale is not printed. A5's set point means it is at least 27 NL/min if one controller delivers every case. With a full scale of 31 NL/min every case except A4 is within the bound. A4 needs 2.04 NL/min. We first held every case to ±2 % of the *set point*, the strictest reading. That fits A3, A5, A6 and A7 but not the other five (Python, the same engine and decks):

| Case | ±2 % of the set point: CO₂ / CH₄ front, T RMSE | Flow from the thermocouples (the decks) |
|---|---|---|
| A1 | +3.9 / +4.6 %, 15.9 K | +0.6 / +1.3 %, 6.9 K |
| A2 | −3.7 / −2.9 %, 19.1 K | +0.8 / +1.8 %, 8.5 K |
| A4 | −8.5 / −7.5 %, 14.3 K | −0.3 / +0.7 %, 7.1 K |
| B1 | +3.6 / −1.1 %, 14.2 K | −0.1 / −4.7 %, 7.4 K |
| B2 | −9.9 / −8.7 %, 25.1 K | −0.1 / +1.7 %, 9.1 K |

So the decks use the thermocouple flow factors.

Results of the paper test (the decks as shipped; front-time errors are signed, simulation minus experiment):

| Case | CO₂ front, before → after | CH₄ front, before → after | NRMSE y_CO₂ / y_CH₄, before → after | Thermocouple RMSE [K], before → after |
|---|---|---|---|---|
| A3 (calibration) | +4.4 → +1.4 % | −3.3 → +2.4 % | 0.18 / 0.61 → 0.083 / 0.18 | 11.4 → 5.5 |
| A1 | +7.3 → +0.6 % | +3.7 → +1.3 % | 0.45 / 0.52 → 0.16 / 0.20 | 20.3 → 6.9 |
| A2 | −5.6 → +0.8 % | −10.6 → +1.8 % | 0.36 / 0.84 → 0.10 / 0.26 | 19.9 → 8.5 |
| A4 | −10.4 → −0.3 % | −18.0 → +0.7 % | 0.29 / 0.79 → 0.023 / 0.29 | 16.4 → 7.1 |
| A5 | −6.5 → −1.4 % | −14.9 → −1.0 % | 0.18 / 0.66 → 0.050 / 0.31 | 14.3 → 11.4 |
| A6 | −2.9 → −0.2 % | −11.6 → +0.2 % | 0.17 / 0.75 → 0.013 / 0.12 | 11.6 → 7.5 |
| A7 | +0.3 → −0.4 % | −9.9 → −0.3 % | 0.024 / 0.69 → 0.014 / 0.26 | 8.7 → 7.2 |
| B1 | +6.4 → −0.1 % | −0.2 → −4.7 % | 0.38 / 0.13 → 0.045 / 0.065 (y_H₂ 0.37 → 0.25) | 17.8 → 7.4 |
| B2 | −10.7 → −0.1 % | −11.6 → +1.7 % | 0.45 / 0.33 → 0.050 / 0.047 (y_H₂ 0.40 → 0.22) | 28.5 → 9.1 |

All tolerances were re-measured: the measured error × 1.1, rounded up, with a floor of 0.2 % for the front times. B1's CH₄ front-time bound rises from 0.004 to 0.055. Before the calibration that front was within 0.2 %, but its neighbours were 6–11 % off. Every other bound falls.

Limits:

- The flow factors use each case's own thermocouples, so the thermocouple RMSE of A1, A2 and A4–B2 is not an independent check of the front speed. It still checks the height and shape of the temperature fronts. The outlet curves of those eight cases are predictions.
- The fit holds the whole sampling line at column pressure. The paper puts a back-pressure regulator in the line, with the gas at ambient pressure after it. The fitted volume is therefore the line's equivalent volume at column pressure. One volume fits 1, 10 and 25 bar, which suggests that most of the line is upstream of the regulator.
- The line gas leaves the bed at the bed-outlet temperature and stays there, because the layer is adiabatic. A real sampling line is at room temperature.
- The CH₄ outlet NRMSE stays at 12–31 % in the binary cases. Most of it comes from the CH₄ roll-up and its fall at CO₂ breakthrough, whose height depends on the competition model (extended Sips here).
- A5 (25 bar) keeps the largest thermocouple RMSE, 11.4 K.

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
disagreement was therefore attributed to the competition model, not to the
column model. The IAST comparison below tests that attribution: IAST moves the
CO₂ front later, not earlier, so the competition model does not explain it.

### IAST versus extended Sips (2026-11-10)

Measured on the decks before the calibration on A3. On the calibrated decks (the same flow factors, line and transport values) IAST puts the CO₂ fronts 0.3–11 % late (extended Sips: within 1.4 %), the CH₄ fronts within 0.2–4.7 % and the thermocouple RMSE at 6.6–22 K (extended Sips: 5.5–11.4 K), so extended Sips stays the decks' mixing rule.

Every case was run a second time with the adsorbent's mixing rule set to IAST
over the pure Sips 2 terms of the same Table 3 parameters
(`Experiment<case>_IAST_RecordOnly`, recorded on the validation page as
`Experiment<case>_IAST`, without assertion). The pure-component limits are
identical: `IastPureLimitsMatchExtendedSips` checks that IAST at y_i = 1
gives the extended-Sips loading to 1e-10 for every component at 25, 45 and
65 °C from 1 mbar to 25 bar. Only the competition changes. It changes a lot
for this pair: at 50:50 CO₂/CH₄, 10 bar, 25 °C, extended Sips gives 5.34 mol/kg
CO₂ and 0.58 mol/kg CH₄, IAST gives 6.13 and 0.005 mol/kg. CO₂'s small Sips
exponent (n = 0.425) gives it a spreading pressure that excludes CH₄ almost
completely.

Errors (fractions for the curves and front times, K for the thermocouples).
Before = extended Sips, the decks' model and its bound; after = IAST.
Bold: IAST better.

| Case | Quantity | Extended Sips | IAST | Bound (extended Sips) |
|---|---|---|---|---|
| A1 | outlet y_CO2 NRMSE | 0.447 | 0.507 | 0.53 |
| A1 | outlet y_CH4 NRMSE | 0.524 | **0.427** | 0.83 |
| A1 | CO2 front time | 0.0728 (3929 s) | 0.0986 (4024 s) | 0.098 |
| A1 | CH4 front time | 0.0374 (476 s) | **0.0203** (468 s) | 0.12 |
| A1 | thermocouple RMSE | 20.3 | 21.6 | 25 |
| A2 | outlet y_CO2 NRMSE | 0.360 | **0.358** | 0.37 |
| A2 | outlet y_CH4 NRMSE | 0.841 | **0.829** | 0.96 |
| A2 | CO2 front time | 0.0563 (911 s) | **0.0497** (1014 s) | 0.058 |
| A2 | CH4 front time | 0.106 (326 s) | 0.125 (319 s) | 0.13 |
| A2 | thermocouple RMSE | 19.9 | **19.5** | 22 |
| A3 | outlet y_CO2 NRMSE | 0.178 | 0.318 | 0.27 |
| A3 | outlet y_CH4 NRMSE | 0.609 | 0.657 | 0.73 |
| A3 | CO2 front time | 0.0444 (433 s) | 0.148 (476 s) | 0.08 |
| A3 | CH4 front time | 0.0333 (151 s) | 0.0551 (148 s) | 0.041 |
| A3 | thermocouple RMSE | 11.4 | 17.9 | 19 |
| A3 | CH4 breakthrough (y > 0.01), text: about 140 s | 150 s | **146 s** | – |
| A3 | CO2 breakthrough (y > 0.01), text: about 400 s | 420 s | 465 s | – |
| A4 | outlet y_CO2 NRMSE | 0.295 | **0.104** | 0.35 |
| A4 | outlet y_CH4 NRMSE | 0.794 | **0.770** | 0.91 |
| A4 | CO2 front time | 0.104 (203 s) | **0.0268** (220 s) | 0.105 |
| A4 | CH4 front time | 0.180 (70.3 s) | 0.202 (68.4 s) | 0.195 |
| A4 | thermocouple RMSE | 16.4 | **10.5** | 20 |
| A5 | outlet y_CO2 NRMSE | 0.185 | 0.211 | 0.19 |
| A5 | outlet y_CH4 NRMSE | 0.663 | 0.670 | 0.76 |
| A5 | CO2 front time | 0.0648 (177 s) | 0.0700 (203 s) | 0.066 |
| A5 | CH4 front time | 0.149 (80.8 s) | 0.167 (79.2 s) | 0.17 |
| A5 | thermocouple RMSE | 14.3 | 14.6 | 21 |
| A6 | outlet y_CO2 NRMSE | 0.171 | 0.248 | 0.175 |
| A6 | outlet y_CH4 NRMSE | 0.752 | 0.766 | 0.85 |
| A6 | CO2 front time | 0.0292 (435 s) | 0.0468 (469 s) | 0.030 |
| A6 | CH4 front time | 0.116 (142 s) | 0.137 (139 s) | 0.14 |
| A6 | thermocouple RMSE | 11.6 | 14.0 | 13 |
| A7 | outlet y_CO2 NRMSE | 0.0244 | 0.240 | 0.14 |
| A7 | outlet y_CH4 NRMSE | 0.691 | 0.726 | 0.72 |
| A7 | CO2 front time | 0.0031 (428 s) | 0.0670 (455 s) | 0.04 |
| A7 | CH4 front time | 0.0994 (132 s) | 0.120 (129 s) | 0.12 |
| A7 | thermocouple RMSE | 8.7 | 13.0 | 14 |
| B1 | outlet y_CO2 NRMSE | 0.377 | 0.386 | 0.45 |
| B1 | outlet y_CH4 NRMSE | 0.127 | 0.136 | 0.15 |
| B1 | outlet y_H2 NRMSE | 0.372 | 0.372 | 0.42 |
| B1 | CO2 front time | 0.0643 (9019 s) | 0.0685 (9055 s) | 0.081 |
| B1 | CH4 front time | 0.0018 (337 s) | 0.0021 (337 s) | 0.004 |
| B1 | thermocouple RMSE | 17.9 | 18.2 | 20 |
| B2 | outlet y_CO2 NRMSE | 0.450 | **0.419** | 0.57 |
| B2 | outlet y_CH4 NRMSE | 0.328 | 0.346 | 0.35 |
| B2 | outlet y_H2 NRMSE | 0.398 | **0.397** | 0.46 |
| B2 | CO2 front time | 0.108 (2405 s) | **0.0861** (2462 s) | 0.11 |
| B2 | CH4 front time | 0.116 (631 s) | 0.128 (622 s) | 0.15 |
| B2 | thermocouple RMSE | 28.5 | **27.5** | 38 |

IAST is better on 14 of the 47 asserted metrics and worse on 33. It would
exceed 14 of the current bounds (A1, A3, A4, A5, A6, A7; none in A2, B1 or B2).
So extended Sips stays the decks' mixing rule, no bound changed, and IAST is
recorded only.

**What the competition model explains.** IAST holds more CO₂ and almost no
CH₄, so every CO₂ front comes later (A1 +2.4 %, A2 +11 %, A3 +9.9 %, A4 +8.6 %,
A5 +14 %, A6 +7.8 %, A7 +6.4 %, B1 +0.4 %, B2 +2.4 %) and every CH₄ half-feed
front 1.4–2.7 % earlier. This helps where extended Sips was early (the CO₂
fronts of A2, A4 and B2) and hurts where it was on time or late (A1, A3, A5–A7).
The CH₄ breakthrough of the reference experiment A3 moves from 150 to 146 s,
nearer the "after 140 s" of the text. That agrees with the paper's statement
about the CH₄ breakthrough time. With the Table 3 pure isotherms and the
Table 4 LDF coefficients, we could not reproduce the paper's statement about
the CO₂ front.

**What remains.**

- *Feed velocity.* A2, A3 and A4 differ only in flow (10, 20 and 40 cm³/s at
  10 bar and 25 °C), so their equilibrium is the same. Their extended-Sips CO₂
  front errors still have opposite signs (−5.6 %, +4.4 %, −10.4 %). A
  competition model shifts all three the same way, so it cannot remove this
  pattern. A per-experiment flow error can, and that is what the paper's
  velocity fit (within ±2 % of the MFC full scale, Sect. 3.2.1) absorbs.
- *Sampling line.* The CH₄ half-feed fronts are 10–20 % early at 10–25 bar
  with either model, and IAST makes them earlier still. The unmodelled
  sampling line delays both fronts by the same seconds. That delay is a larger
  fraction of a 70–160 s CH₄ front than of a 200–4000 s CO₂ front.
- *Centre-line temperature.* The thermocouple RMSE follows the CO₂ front
  timing (IAST −6 K in A4, +6 K in A3). The remaining 9–28 K is the timing
  and height of the temperature fronts at the centre line, which the paper
  also notes its 1D model does not capture.

**Engine changes this needed (`core/src/equations/iast.cpp`).** Under IAST
the decks first failed within 10⁻¹⁰ s. There were two independent causes:

- The front toes reach subnormal mole fractions (y = 5 × 10⁻³²¹), where the
  spreading pressure underflows and the solve returned NaN. Mole fractions
  at or below 10⁻¹⁵⁰ now take IAST's first-order dilute limit. With equal
  capacities and unit exponents (IAST = extended Langmuir) the A3 run now
  reproduces the extended run (CO₂ front 542.57 s against 542.57 s); before
  the fix it failed although the residuals agreed to 10⁻²⁰.
- CO₂'s Sips term (n < 1) has an infinite Henry slope, so the dilute loading
  of an absent component diverged (−1.3 × 10⁷ mol/kg at y = −10⁻¹⁴). Inside
  the IAST solve the pure term is continued below p₀ = 10⁻⁴ P by a monotone C¹
  form whose spreading pressure at p₀ is exact. IAST is therefore unchanged
  for every component with y ≥ 10⁻⁴. At 10⁻⁵ the continuation is too steep
  (A3 covered 0.1 s of 4100 s in a minute). At 10⁻⁴ and at 10⁻³ every metric
  above agrees to 0.001.

**Run time.** Wall time per case, measured with two cases at a time on four
cores (extended Sips → IAST): A1 21 → 162 s, A2 10 → 46 s, A3 9 → 39 s,
A4 6 → 32 s, A5 6 → 28 s, A6 7 → 34 s, A7 6 → 31 s, B1 60 → 214 s,
B2 29 → 86 s. That is 3–8 times slower (IAST cost per residual, issue #162).
The nine recorded IAST cases add about 11 CPU-minutes to the paper tier,
whose per-test limit is 3600 s.
