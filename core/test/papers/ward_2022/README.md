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
| Initial gas / carrier | **air in the deck; helium in the experiment** | Sect. 2: "the CO₂/He system" | the paper's model uses constant CO₂ gas properties (C_p,g 840 J/(kg K), μ 1.72 × 10⁻⁵ Pa s, k_g 0.0162 W/(m K), D_m 0.605 × 10⁻⁵ m²/s; Table 1); the engine takes mixture properties from CoolProp |
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

- The CO2 front (y = 0.5) arrives 4-6 % earlier than measured. The simulated front is much sharper than the measured one, so the curve NRMSE is 34-38 %; the thermocouple RMSE is 2.6-4.2 K, and the simulated temperature rise (15-33 K) is above the paper model's 95 % interval (7-27 K).
- No transcription error was found (column length, radius, voidage, bed density, flows and isotherm match Tables 1-2). Likely causes: a single constant heat of adsorption of 29.98 kJ/mol (the site-1 internal energy, while most of the loading sits on site 2 at 20.4 kJ/mol) heats the bed more and lowers its capacity; the deck uses air as carrier where the paper uses He; the sampling line is not modelled. The ambient temperature was corrected to Table 1's 291.4 K.
