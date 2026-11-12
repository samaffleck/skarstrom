# Ward et al. (2022) — CO₂ breakthrough on activated carbon (lab-scale rig)

## Citation

> Ward, A., Pini, R. *Integrated uncertainty quantification and sensitivity
> analysis of single-component dynamic column breakthrough experiments.*
> Adsorption 28:161–183 (2022).
> DOI: [10.1007/s10450-022-00361-z](https://doi.org/10.1007/s10450-022-00361-z)


## Files

| File | Purpose |
|------|---------|
| `ward_2022.cpp` | `Ward2022Test` fixture + reproduced breakthrough cases |

## Reproduced cases

| Test (`Ward2022Test.*`) | Notes |
|-------------------------|-------|
| `Q50ccm`  | CO₂ breakthrough at 50 cm³/min feed |
| `Q100ccm` | CO₂ breakthrough at 100 cm³/min feed |
| `Q200ccm` | CO₂ breakthrough at 200 cm³/min feed |
| `MTC_Invariance` | Sanity check: equilibrium loading must be invariant to the mass-transfer coefficient (breakthrough may shift, equilibrium must not) |

## Provenance (WP10 audit, 2026-10-31, from the paper text)

Every deck value, with where it comes from. "Paper model" is what Ward & Pini's
own 1D model used (Sect. 3, Appendix A); the experiment is that of their ref. [5]
(Pini, Joss, Hejazi 2021).

This table records the audit. Since the calibration (2026-11, see
[Calibration](#calibration-2026-11) below) the decks replace four of these values
with calibrated ones: the ambient and wall-end temperature, q_s,2, the factor on
Eq. 4 and the inside HTC. They also use the particle density 435 / 0.544
(799.6 kg/m³) in place of 799.

| Deck value | Value | Source | Note |
|---|---|---|---|
| Column length | 0.286 m | Table 1, L | |
| Column inner diameter | 0.0189 m | Table 1, r_in = 0.00945 m | |
| Wall thickness | 3.25 mm | Table 1, r_out − r_in = 0.0127 − 0.00945 | aluminium: ρ_w 2700, C_p,w 902, K_w 205 (Table 1) |
| Bed voidage | 0.456 | Table 1, ε | |
| Particle diameter | 3 mm | Table 1, r_p = 1.5 × 10⁻³ m | |
| Particle density | 799 kg/m³ | derived: ρ_b / (1 − ε) = 435 / 0.544 = 799.6 (Table 1) | |
| Particle voidage | 0 in the deck (0.59 in the LDF) | Table 1, ε_p = 0.59 | the paper's balances (A1, A2) hold up gas in the bed voidage only, so the deck gives the particle no gas holdup; 0.59 is used in the Haynes coefficient as the paper does (Eq. 6) |
| Adsorbent C_p, k_s | 1040 J/(kg K), 0.63 W/(m K) | Table 1 | |
| Adsorbed-phase C_p | 36.96 J/(mol K) | Table 1, C_p,a = 840 J/(kg K) × 0.044 kg/mol | |
| Ambient temperature | 291.4 K | Table 1, T_a | |
| Feed temperature, pressure | 293 K, 1 bar | Table 1, T_F, p_F | |
| Feed | pure CO₂ at 50 / 100 / 200 cm³/min | Sect. 2 | |
| Initial gas / carrier | **air in the deck; helium in the experiment** | Sect. 2: "the CO₂/He system" | the paper's model uses constant CO₂ gas properties (C_p,g 840 J/(kg K), μ 1.72 × 10⁻⁵ Pa s, k_g 0.0162 W/(m K), D_m 0.605 × 10⁻⁵ m²/s; Table 1); the engine takes the other mixture properties from CoolProp; the gas heat capacity is the paper's: constant 840 J/(kg K) per component, i.e. 840 J/(kg K) for every CO₂/He mixture (see below) |
| Isotherm | DSL, concentration form | Table 2, Eqs. 1–3 | q_s,1 0.418, q_s,2 6.04 mol/kg; b_0,1 2.05 × 10⁻⁶, b_0,2 3.32 × 10⁻⁶ m³/mol; ΔU_1 −29 980, ΔU_2 −20 367 J/mol |
| Heat of adsorption | 29 980 J/mol | Table 2 (−ΔU_1) | the paper's energy balance (A4) releases −ΔU_1 + RT ≈ 32.4 kJ/mol per mole adsorbed, for all loading |
| LDF coefficient | Haynes (film + macropore), τ_p 3, ε_p 0.59 | Eqs. 6–7, Table 1 | the engine evaluates D_m per cell from the mixture; the paper uses the constant 0.605 × 10⁻⁵ m²/s |
| Axial dispersion | **Ruthven, D_m = 1.6 × 10⁻⁵ m²/s, no factor** | Eq. 4: D_L = 0.7 D_m + v_F r_p, **× 4** (Sect. 3.2) | the paper multiplies Eq. 4 by 4 "for an optimal match to the experimental outlet composition profiles"; at 100 ccm the deck's D_L is about a third of the paper's |
| Axial thermal conductivity | Yagi–Kunii | Eqs. 8–10 (K_z0 = k_g (k_s/k_g)^n, n = 0.28 − 0.757 log ε − 0.057 log(k_s/k_g); K_z = K_z0 + 0.75 k_g Pr Re_p) | the engine's "Yagi and Kunii" model is this correlation |
| Inside HTC | 77.5 / 57.9 / 45.2 W/(m² K) | Table 6 (fitted per flow) | |
| Outside HTC | Holman, c 0.027, m 0.8, n 1/3, u 0.5 m/s | Eq. 5 and text | |
| Pressure drop | **Kozeny–Carman (180)** | Eq. A3: Blake–Kozeny with 150 | negligible at 1 bar and these flows |
| Wall end condition | **zero flux** | Table 4: T_w(0) = T_w(L) = T_a ("tight thermal contact") | the engine has a fixed-end option |
| Cells | 60 | Table 5: 30 | |
| Thermocouple | z = 0.0715 m | Sect. 3.2: z = L/4 | |
| Indicators | t_br at y_out = 0.5; ΔT_ads = max T_tc − T_F | Table 9 | the test used the initial 291.4 K instead of T_F until this audit |

Bold entries differ from the paper's model; the attribution runs below measure
what each one does.

### Attribution runs (2026-10-31)

One model-form change at a time on the bundled decks (Python, same engine),
against the digitised experiment. NRMSE of the outlet CO₂ curve; the front
time where y_out = 0.5 (Table 9's t_br); the thermocouple RMSE; ΔT_ads per
Table 9 (measured: 13.1 K at 100 ccm, 6.4 K at 50, 20.3 K at 200).

| Variant (100 ccm) | y NRMSE | t_br sim / exp [s] | T RMSE [K] | ΔT_ads sim [K] |
|---|---|---|---|---|
| deck before the audit (air, Ruthven D_m 1.6e-5, 29.98 kJ/mol, zero-flux wall ends, Kozeny–Carman) | 0.379 | 1182 / 1259 (−6.1 %) | 4.24 | 24.9 |
| + paper's dispersion (Eq. 4 × 4, D_m 0.605e-5) | 0.400 | 1150 / 1259 (−8.6 %) | 3.66 | 15.1 |
| + −ΔU₁ + RT | 0.383 | 1169 / 1259 (−7.1 %) | 4.98 | 26.3 |
| + helium carrier | 0.378 | 1183 / 1259 (−6.1 %) | 3.97 | 23.9 |
| + wall ends at T_a (Table 4) | 0.368 | 1324 / 1259 (+5.1 %) | 3.30 | 22.8 |
| + Ergun | 0.383 | 1169 / 1259 (−7.1 %) | 4.55 | 25.0 |
| all five (the paper's model) | **0.129** | **1280 / 1259 (+1.6 %)** | **2.65** | **13.2** |

| Variant | 50 ccm: y NRMSE, t_br error, T RMSE, ΔT_ads | 200 ccm: y NRMSE, t_br error, T RMSE, ΔT_ads |
|---|---|---|
| deck before the audit | 0.382, −4.3 %, 2.64 K, 14.4 K | 0.336, −4.6 %, 5.24 K, 31.6 K |
| paper's dispersion + wall ends at T_a | 0.185, +1.9 %, 1.42 K, 5.4 K | 0.132, +1.7 %, 3.87 K, 21.6 K |
| all five | 0.184, +1.9 %, 1.38 K, 5.6 K | 0.029, +0.1 %, 5.50 K, 22.9 K |

Two forms carry the result: the paper's (fitted) dispersion coefficient halves
the temperature rise and brings it inside the paper's 95 % interval, and the
wall ends held at ambient move the front from 4–6 % early to 1–2 % late. The
carrier gas, the RT term in the heat and the pressure-drop constant are
each worth less than 1 % in the front time. The heat-of-adsorption
explanation recorded before the audit (that a single site-1 value "heats the
bed more") was wrong: the paper's own model uses that value plus RT.

The decks now use all five forms, i.e. the paper's model; the remaining
error (the table's last rows) is the engine-vs-paper-model difference plus
the digitising error of the curves.

## Validation

The tests compare the simulation with the paper using the helpers in
`../../fixtures/validation_metrics.h`; the published numbers they use are in
[`reference_values.csv`](reference_values.csv) (each with its table / figure /
page). Every tolerance is the measured error rounded up, with the measured
value in a comment next to it. Results of the last run:
[validation page](../../../../docs/validation.md#ward-2022)
(generated from `../validation_summary.json`).

- With the calibrated decks ([Calibration](#calibration-2026-11)), the CO2 front (y = 0.5) is within 0.4 % of the measured time at all three flows, the outlet-curve NRMSE is 2.7-3.9 % and the thermocouple RMSE 0.18-1.40 K. The temperature rise per Table 9 is 6.7 / 13.1 / 20.9 K against the measured 6.4 / 13.1 / 20.3 K.
- The sections below record the history of the decks, oldest first. Their numbers are the results at the time.

### Danckwerts inlet fix (2026-11)

The inlet face of the bed used to let axial dispersion carry CO₂ into the bed
on top of the feed, eb·D·C·(1 − y₁)/Δz with the ghost cell at the feed's
y = 1. With the paper's (×4) dispersion coefficient that is about 2 % of the
fed CO₂ at 50 cm³/min (the per-step component balance against the stream log
was off by 2.1 %), which the stream log and the breakthrough indicators never
saw. The inlet is now Danckwerts (the face carries exactly the feed's
u·C·y_in), so the bed holds the same capacity but is fed less CO₂ and the
front arrives later:

| Case | front time sim / exp [s], before → after | y NRMSE, before → after | T RMSE [K], before → after |
|---|---|---|---|
| 50 ccm | 2623 → 2682 / 2574 (+1.9 % → +4.2 %) | 0.184 → 0.332 | 1.19 → 1.18 |
| 100 ccm | 1280 → 1301 / 1259 (+1.7 % → +3.3 %) | 0.129 → 0.244 | 2.06 → 1.18 |
| 200 ccm | 591.7 → 599.8 / 591.1 (+0.1 % → +1.5 %) | 0.029 → 0.120 | 4.01 → 2.76 |

The temperature fit improves (the inlet energy balance now carries exactly
the feed enthalpy, with no conduction across the inlet face), the front-time
and curve fits get worse: part of the previous agreement came from the
unaccounted CO₂. The front-time and NRMSE tolerances were re-measured
accordingly (see `ward_2022.cpp`). At 100 ccm ΔT_ads is now 12.89 K, just
below the paper model's 95 % interval [12.9, 17.2] K (informational metric).

### Gas heat capacity as in the paper (2026-11)

Table 1 gives one constant heat capacity, 840 J/(kg K), for the CO₂/He gas and the adsorbed phase. The deck now uses the constant gas heat capacity mode with 840 J/(kg K) for each component, converted with its own molar mass (CO₂ 36.97, He 3.362 J/(mol K)). The two differ by a factor of 11 per mole, but the engine forms the mixture's heat capacity per kilogram as Σ y_i C_p,i / Σ y_i M_i, so every CO₂/He mixture has exactly the paper's 840 J/(kg K). Before, the default ideal-gas C_p0(T) (He 20.8 J/(mol K), i.e. 5193 J/(kg K)). Changes: Q200 front-time error 1.724 → 1.735 % (bound 1.8 %), outlet NRMSE 0.1402 → 0.1410 (bound 0.17), thermocouple RMSE 2.58 → 2.56 K; Q50 and Q100 move by less than 0.002 K and 0.0001 in NRMSE. No bound changed. All papers' settings: [../README.md](../README.md#gas-heat-capacity-per-paper).

### Grid and advection scheme (2026-11-04)

The decks keep their 60 cells and now run WENO, the engine's default, instead
of first-order upwind. A grid study (10-160 cells, Richardson extrapolation;
[Numerics](../../../../docs/numerics.md#advection-schemes)) found the CO2
half-feed time within 0.5 % of its grid limit on about 23 WENO (40 upwind)
cells, so 60 are converged. The outlet curves are sharper than the digitised
measurements: every outlet NRMSE rose by about 0.02 and Q200's passed its
bound (0.14 → 0.17); the front times moved by 0.1-0.3 % (later), and the
temperature rises at the thermocouple moved 0.4-0.9 K towards the
experiment. Before: upwind; after: WENO; both on 60 cells (bold: changed
bound).

| Case | Quantity | Reference | Before | After | Error before → after | Bound before → after |
|---|---|---|---|---|---|---|
| Q100ccm | outlet y_CO2 | 1.00 | 0.2427 | 0.2645 | 0.2427 → 0.2645 | 0.29 |
| Q100ccm | CO2 front time (y_CO2 = 0.5) | 1258.95 | 1300.35 | 1303.77 | 0.03289 → 0.0356 | 0.039 |
| Q100ccm | thermocouple T at z = 0.0715 m | 0 | 1.18 | 1.15 | 1.178 → 1.147 | 2.4 |
| Q100ccm | breakthrough time (y_CO2,out = 0.5) | 1249.50 | 1300.35 | 1303.77 | 50.85 → 54.27 | – |
| Q100ccm | temperature rise at z = 0.0715 m | 15.05 | 12.90 | 13.65 | 2.15 → 1.40 | – |
| Q200ccm | outlet y_CO2 | 1.00 | 0.1194 | 0.1404 | 0.1194 → 0.1404 | **0.14 → 0.17** |
| Q200ccm | CO2 front time (y_CO2 = 0.5) | 591.08 | 599.74 | 601.28 | 0.01466 → 0.01726 | 0.018 |
| Q200ccm | thermocouple T at z = 0.0715 m | 0 | 2.76 | 2.57 | 2.759 → 2.574 | 4.6 |
| Q200ccm | breakthrough time (y_CO2,out = 0.5) | 575.00 | 599.74 | 601.28 | 24.74 → 26.28 | – |
| Q200ccm | temperature rise at z = 0.0715 m | 23.60 | 22.54 | 23.43 | 1.06 → 0.17 | – |
| Q50ccm | outlet y_CO2 | 1.00 | 0.332 | 0.3521 | 0.332 → 0.3521 | 0.39 |
| Q50ccm | CO2 front time (y_CO2 = 0.5) | 2573.62 | 2681.75 | 2689.05 | 0.04201 → 0.04485 | 0.049 |
| Q50ccm | thermocouple T at z = 0.0715 m | 0 | 1.18 | 1.24 | 1.177 → 1.243 | 1.4 |
| Q50ccm | breakthrough time (y_CO2,out = 0.5) | 2540.50 | 2681.75 | 2689.05 | 141.25 → 148.55 | – |
| Q50ccm | temperature rise at z = 0.0715 m | 8.45 | 5.42 | 5.85 | 3.03 → 2.60 | – |

### Calibration (2026-11)

After the steps above the decks used the paper's model with the paper's values,
and the CO2 front still arrived 1.7-4.5 % late (most at 50 ccm), with an
outlet NRMSE of 14-35 % and a thermocouple RMSE of 1.15-2.56 K. A recheck of every
deck value against Tables 1, 2, 4, 5 and 6 found no transcription error. The
equations of the paper (Eqs. 4-10, Appendix A) were not available in readable
form for this recheck. The only remaining difference was the particle density,
which the decks rounded to 799 kg/m³ (Table 1 gives 435 / 0.544 = 799.6).

The measured thermocouple gives the cause. It reads 291.5 K at the start (Table
1's ambient, 291.4 K), but after the front it settles at about 293 K, the feed
temperature, and stays there to the end of the run. With the wall ends held at
291.4 K (Table 4), the simulated bed cools to 291.4 K behind the front. A cooler
bed holds more CO2, so the front arrives late, and most at the lowest flow,
where the bed has the most time to cool. Insulated wall ends (zero flux) do not
fix this: the bed then keeps its heat and the front arrives 5-10 % early.

The decks are now calibrated. A weighted least-squares fit (scipy
`least_squares`; residual weights 0.03 in y and 0.5 K, the assumed digitising
error) to the outlet y_CO2 and thermocouple points of all three flows together
gave these values. The ambient temperature was set to 293 K from the measured
thermocouple, not fitted.

| Value | Calibrated (± 1 σ of the fit) | Paper |
|---|---|---|
| Ambient and wall-end temperature | 293 K (set) | 291.4 K (Table 1) |
| q_s,2 | 5.89 ± 0.01 mol/kg | 6.04 mol/kg (Table 2); 95 % interval 5.64-6.42 (Table 7) |
| Factor on Eq. 4 | 3.75 ± 0.12 | 4 (Sect. 3.2, itself fitted to the outlet curves) |
| h_in, 50 / 100 / 200 ccm | 97.3 ± 6.6 / 81.5 ± 2.4 / 69.0 ± 1.1 W/(m² K) | 77.5 / 57.9 / 45.2 W/(m² K) (Table 6, fitted) |

The rate model keeps the engine's Haynes coefficient (no capacity factor, the
mixture's molecular diffusivity, tortuosity 3). The coefficient is about 40 1/s,
so the bed is close to equilibrium. A fit with the capacity factor rho_p q*/c
and the constant D_m of Table 1 in the coefficient (about 0.15 1/s) gave an
almost equal fit (total cost 267 against 314). Eq. 6 of the paper could not be
checked for the capacity factor, so the decks do not use it. A constant gas
viscosity (Table 1: 1.72e-5 Pa s) had no measurable effect on the fit; the
decks keep the engine's viscosity model, which keeps the 320-cell outlet in the
regime that `FineGridBreakthrough` tests.

One change at a time on the decks before the calibration (Python, the same
runs as the paper test; front-time error, outlet NRMSE and ΔT_ads per Table 9;
measured ΔT_ads 6.4 / 13.1 / 20.3 K):

| Variant | 50 ccm | 100 ccm | 200 ccm |
|---|---|---|---|
| decks before the calibration | +4.5 %, 35.2 %, 5.9 K | +3.6 %, 26.5 %, 13.7 K | +1.7 %, 14.1 %, 23.4 K |
| particle density 799.6 only | +4.6 %, 35.5 %, 5.9 K | +3.6 %, 26.9 %, 13.7 K | +1.8 %, 14.7 %, 23.4 K |
| ambient and wall ends at 293 K only | +1.6 %, 16.1 %, 7.4 K | +0.7 %, 5.1 %, 15.1 K | −1.0 %, 9.0 %, 24.3 K |
| q_s,2 5.89 only | +2.3 %, 22.2 %, 5.9 K | +1.4 %, 11.1 %, 13.6 K | −0.3 %, 4.2 %, 23.3 K |
| factor 3.75 on Eq. 4 only | +4.6 %, 36.2 %, 6.1 K | +3.7 %, 27.5 %, 14.1 K | +1.8 %, 15.1 %, 23.8 K |
| h_in 97.3 / 81.5 / 69.0 only | +4.8 %, 36.5 %, 4.8 K | +4.8 %, 31.9 %, 11.3 K | +5.2 %, 29.3 %, 19.5 K |
| ambient + q_s,2 | −0.5 %, 5.4 %, 7.4 K | −1.4 %, 12.4 %, 15.1 K | −3.0 %, 21.3 %, 24.2 K |
| all (the calibrated decks) | +0.0 %, 2.8 %, 6.7 K | −0.0 %, 2.7 %, 13.1 K | +0.4 %, 3.9 %, 20.9 K |

The ambient temperature carries most of the front-time correction, and q_s,2
the rest. The larger h_in then takes the excess heat out of the bed at the
higher ambient, and the dispersion factor sets the front width.

The fit is reproducible with the calibration API: `examples/calibrate_ward_2022.py`
(`skarstrom.calibrate`, the same cases, datasets, weights and parameters, from the
paper's values) gives q_s,2 5.89 ± 0.01 mol/kg, the factor 3.74 ± 0.12 and h_in
96.7 ± 6.7 / 81.6 ± 2.4 / 69.0 ± 1.1 W/(m² K), within 0.1 σ of the table above and
with the same standard errors. Its estimate has a slightly lower objective than the
rounded values above (S = 639.9 against 642.9), mostly from q_s,2 = 5.894 before
rounding.

Results of the paper test (the decks as shipped):

| Case | Front-time error, before → after | Outlet NRMSE, before → after | Thermocouple RMSE [K], before → after |
|---|---|---|---|
| 50 ccm | +4.49 % → +0.021 % | 0.352 → 0.0276 | 1.243 → 0.181 |
| 100 ccm | +3.56 % → +0.019 % | 0.265 → 0.0272 | 1.147 → 0.675 |
| 200 ccm | +1.73 % → +0.40 % | 0.141 → 0.0389 | 2.555 → 1.403 |

The tolerances are re-measured (front-time bounds with a floor of 0.2 %).

Limits:

- The calibration is fitted to the same digitised curves that the test compares
  with, so this is a calibration, not an independent validation. The raw data of
  ref. [5] would allow an independent check.
- At 200 ccm the simulated bed cools slightly too fast just after the
  temperature peak. ΔT_ads at 50 ccm (6.7 K) is below the paper model's 95 %
  interval [7.36, 9.53] K, as is the measured rise (6.4 K).
- The adsorbent library entry `activated_carbon_ward2022` keeps the published
  Table 2 isotherm (q_s,2 6.04 mol/kg). It is no longer deck-checked against
  `Ward2022_Q50ccm`.
