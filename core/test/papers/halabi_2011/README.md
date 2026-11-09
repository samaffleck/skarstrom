# Halabi (2011) — sorption-enhanced steam methane reforming (Rh catalyst + K-HTC)

## Citation

M. H. Halabi, *Sorption enhanced catalytic reforming of methane for pure
hydrogen production: experimental and modeling*, PhD thesis, Eindhoven
University of Technology, 2011.
[doi:10.6100/IR709035](https://doi.org/10.6100/IR709035) (open access:
<https://pure.tue.nl/ws/files/3339515/709035.pdf>).

Chapters used: 8 (the Rh/Ce0.6Zr0.4O2 rate model, Tables 1–4; the uptake model
Eq. 6; particle data, Table 12) and 9 (the fixed-bed experiments: Tables 1–4
and 6, Fig. 7, pp. 191–205). Page numbers are the thesis's printed page numbers.

## Files

- `halabi_2011.cpp`: fixture `Halabi2011Test`. Builds every deck from kind names
  (`deck_builder.h`) and the parameter registry, and holds the thesis's constants
  with their kPa → bar conversions.
- `reference_values.csv`: the numbers the thesis prints (Fig. 7 legend, Table 6,
  p. 203), each with its source.
- `assets/` (git-ignored): the thesis PDF and a local text extraction.

No bundled example deck: the Python package has no reaction builder yet, and the
decks need the reaction structure (three LHHW reactions of 14 adsorption terms each).

## Reproduced cases

| Test | What it does | Asserted |
|---|---|---|
| `GasPhaseEquilibrium` | Ch. 8 rate model on an **inert** sorbent bed (no CO₂ uptake), 5 s of feed. Outlet CH₄ conversion against the equilibria the thesis prints: Fig. 7 legend (3 vol% CH₄, S/C 4; 1.5 / 4.5 / 6 bar) and the "CH₄ equil. conversion" column of Table 6 (3.94 vol% CH₄, S/C 3–6; 1.5 / 4.5 / 6 bar). 15 runs, ~17 s. | yes |
| `SorptionEnhanced_RecordOnly` | All 12 Table 6 runs with the sorbent active, 900 s of feed. Records the conversion and dry Ar-free H₂ at 60 s, the maximum conversion after start-up, and the CO₂ breakthrough time. ~35 s. | no ([below](#sorption-enhanced-runs-recorded-not-asserted)) |

Fig. 7's caption and the text on p. 205 give its feed as 3 vol% CH₄, 0.1 vol% H₂
and S/C 4, at 75 nml/min. The thesis's own K_I and K_II give 57.15 / 40.14 / 36.34 %
for that feed, against the printed 57.1 / 40.1 / 36.8 %. S/C 6 at 3.94 % would
give 61.2 / 43.0 / 39.0 %, so the caption's feed is the one used.

## Provenance

Registry paths use the deck's names (`adsorbent[K-HTC]`, `adsorbent[Rh catalyst]`,
`layer[Bed]`, `step[Feed]`).
`rx = adsorbent[Rh catalyst].reaction[R1|R2|R3]`;
`kin = adsorbent[K-HTC].kinetics[CarbonDioxide]`.

| Value | Deck parameter (registry path) | Source (table / equation / page / sentence) | Note (derived arithmetic, assumption status) |
|---|---|---|---|
| 9.0 mm | `bed.diameter` | Ch. 9, p. 193: "quartz tube with an internal diameter of 9.0 mm" | |
| 4.2 cm | `layer[Bed].length` | Ch. 9, p. 195: "The final bed height inside the reactor was 4.2 cm" | |
| 0.3 mm | `layer[Bed].particle_diameter` | Ch. 9, p. 195: "average size of the catalyst and the sorbent pellets … 0.3 mm" | |
| 2.0 g K-HTC + 0.120 g catalyst, one mixed layer | `layer[Bed]` composition 0.9434 / 0.0566 by mass | Ch. 9, p. 195: "loaded with 120 mg of catalyst mixed with 2 g sorbent" | w_cat = 0.120 / 2.120 |
| 2000 kg/m³ | `adsorbent[K-HTC].solid_density`, `adsorbent[Rh catalyst].solid_density` | Ch. 8, Table 12, "Pellet density" | **assumption** for both solids. Ch. 9 gives no pellet density, only bed densities (Table 2: 748.6 kg/m³ for the sorbent bed; p. 201: 50.4 and 748.6 kg/m³ in the mixed bed) |
| 0.603 | `layer[Bed].bed_void_fraction` | derived from the masses, the bed volume and the pellet density | ε_b = 1 − 2.12 g / (π/4 · 0.009² · 0.042 m³) / 2000 = 0.603. Moves the results by several points ([attribution](#attribution-runs-final-deck)) |
| 0.074 | `adsorbent[K-HTC].particle_void_fraction` | Ch. 9, Table 2: pore volume 0.037 cm³/g | 0.037 cm³/g × 2000 kg/m³ |
| 0.40 | `adsorbent[Rh catalyst].particle_void_fraction` | Ch. 9, Table 1: pore volume 0.20 cm³/g | 0.20 cm³/g × 2000 kg/m³ |
| 450 °C, isothermal | `initial.temperature`, `step[Feed].bottom[Feed inlet].temperature`, isothermal energy balance | Ch. 9, Table 6 heading (450 °C); p. 197: "temperature gradient along the bed … was ± 3 °C" (the 4.0 cm breakthrough bed) | isothermal assumed from the ±3 °C |
| Feed, Table 6 | `step[Feed].bottom[Feed inlet].y[*]` | Ch. 9, p. 203: "3.94 vol% CH4, 0.1vol% H2, variable S/C ratio, and balance Ar … 75 nml/min" | y_H2O = S/C × 0.0394, Ar the balance |
| Feed, Fig. 7 | as above | Ch. 9, p. 205 and the Fig. 7 caption: 3 vol% CH₄, 0.1 vol% H₂, S/C 4, 75 nml/min | |
| 5.577×10⁻⁵ mol/s | `step[Feed].bottom[Feed inlet].molar_flow` | 75 nml/min (above) | "nml" taken as 0 °C, 1 atm: 101325 × 1.25×10⁻⁶ / (R × 273.15) |
| 1.5 / 4.5 / 6.0 bar | `step[Feed].top.pressure`, `initial.pressure` | Ch. 9, Table 6 | back-pressure at the outlet |
| R1, R2, R3 | `rx` (kind `lhhw`, partial-pressure basis, reversible) | Ch. 8, Table 1, Eqs. 1–3 | orders: R1 CH₄ 1, H₂O 1, H₂ −2.5; R2 CO 1, H₂O 1, H₂ −1; R3 CH₄ 1, H₂O 2, H₂ −3.5 |
| k₁, k₂, k₃ | `rx.rate.pre_exp_factor`, `rx.rate.activation_energy` | Ch. 8, Table 3: 1.62×10⁸ kPa^0.5, 2.34×10⁵ kPa⁻¹, 4.55×10⁸ kPa^0.5 mol/(kg_cat s); 83.8, 15.1, 89.2 kJ/mol | to bar: × 100^s, s the total order (−0.5, 1, −0.5): 1.62×10⁷, 2.34×10⁷, 4.55×10⁷ |
| K_I, K_II, K_III | `rx.rate.equilibrium_pre_exp_factor`, `rx.rate.equilibrium_energy` | Ch. 8, Table 2: 1.198×10¹⁷ exp(−26830/T) kPa², 1.767×10⁻² exp(4400/T), 2.117×10¹⁵ exp(−22430/T) kPa² | kPa² → bar²: × 10⁻⁴; E = θ R |
| Ω_l Ω_S | `rx.rate.adsorption[1..14].*`, `rx.rate.denominator_exponent` = 1 | Ch. 8, Table 1, Eq. 4; Table 4 (K_i,o, ΔH_i) | **form**: the engine's LHHW has one denominator. The product (1 + Σ_l)(1 + Σ_S) is expanded exactly into 1 + 4 + 2 + 8 terms; a product term has the product of the K_o, the sum of the ΔH and the sum of the exponents. Van 't Hoff K = K_o exp(−ΔH/RT). To bar: × 100^s with s the term's degree (CH₄ term p_CH₄/p_H₂^0.5: s = 0.5; CO, CO₂, H₂: s = 1; H₂O term p_H₂O/p_H₂: s = 0). The table labels every K_i,o "kPa⁻¹"; the conversion follows the degree of each term |
| Effectiveness 1 | `rx.effectiveness` | Ch. 9, Table 1 (300 µm pellets); ch. 8, p. 163 (intrinsic kinetics, no diffusion limitation) | **assumption**: no intraparticle limitation in 300 µm catalyst pellets |
| Heat of reaction 0 | `rx.heat_of_reaction` | | isothermal run: does not enter |
| Freundlich q = k p^(1/n), wet, k = 1.11 mol/(kg bar^(1/n)), n = 1.76 | `adsorbent[K-HTC].isotherm[CarbonDioxide].a` = k^n × 10⁻⁵ Pa⁻¹, `.v` = n | Ch. 9, Table 4 ("Wet conditions"), Eq. 2 | the engine's Freundlich is q = (a p)^(1/v) with p in Pa. **Measured at 400 °C, used at 450 °C as is** (the thesis gives no temperature dependence). Steam, H₂, CH₄, CO and Ar are inert on the sorbent (ch. 8, assumption 7, p. 166) |
| Uptake, Eq. 6 | `kin` series resistance: macropore term only, capacity factor `slope`, `kin.tortuosity` = 1, `kin.macropore_void_fraction` = 0.4, `kin.molecular_diffusion_coefficient` = 6.6×10⁻⁷ m²/s | Ch. 8, Table 5, Eq. 6 (D_p = 1.1×10⁻⁶ m²/s); Table 12 (tortuosity 1, pellet voidage 0.4) | see [Uptake](#uptake) |
| Ergun momentum, Ruthven dispersion, 40 cells | `layer[Bed]` | Ch. 8, Table 11 (Ergun) | Ruthven with the engine default D_m = 1.6×10⁻⁵ m²/s (**assumption**; no value in the thesis for this bed). 80 cells moves every recorded value by ≤ 0.05 point |
| Start-up | `initial.y[*]` | | see [Start-up](#start-up) |

### Uptake

Ch. 8, Eq. 6: k = 15 ε_p D_p / (r_p² [ε_p + (1 − ε_p) ρ_p R T ∂q*/∂p_CO₂]).
The engine's macropore term with the isotherm-slope capacity factor is
k = 15 ε_p D / (τ r_p² K), with K = ρ_p ∂q*/∂c = ρ_p R T ∂q*/∂p. With D = (1 − ε_p) D_p
the two agree up to a factor 1 + ε_p / ((1 − ε_p) K). For this isotherm K is above 70
at every CO₂ partial pressure below 1 bar (about 380 at 0.024 bar), so the difference
is below 1 %. The
deck uses this form. With D = D_p instead, the uptake is 1.67× faster and the
60 s conversions at 1.5 bar rise by 0.6–1.9 points
([attribution](#attribution-runs-final-deck)).

### Start-up

From the experiment's starting bed (H₂O / H₂ / Ar after the reduction step,
pp. 192, 195), the integrator fails at the first steps: "too much work". The time
scale of the reactions is microseconds, and the rates hold p_H₂^−2.5 to p_H₂^−3.5
at the H₂-poor front. This is skarstrom#231. The deck therefore starts the bed
in the gas-phase equilibrium of the feed at the run pressure (computed in the
test from the thesis's K_I and K_II), with its CO₂ replaced by Ar, so the
sorbent starts clean. A start without CH₄ as well moves the recorded 60 s
conversions by ≤ 0.17 point and the breakthrough times not at all. The outlet
carries the start-up gas for about one residence time (~0.7 s at 1.5 bar), and
the first second shows a CO₂ pulse of up to 0.16 % from the shift reaction. The
test therefore takes the maximum conversion over t ≥ 5 s, and the breakthrough
as the last upward crossing of the 5 % level.

## Model-form audit

- **Rate model.** Ch. 8 Tables 1–4 in the engine's `lhhw` kind, partial pressures in bar.
  The product of the two site balances is expanded into one denominator; this is
  exact. Reversible, with the thesis's K_I, K_II, K_III. Table 2's K_III is
  exactly K_I K_II, so the test's equilibrium solver uses R1 and R2.
- **Isotherm.** Freundlich, not the Langmuir of ch. 8 (Table 5: m = 0.65 mol/kg,
  b = 23.6 bar⁻¹). Ch. 9 fits Freundlich to its own sorbent, and the experiments
  are ch. 9's. No temperature dependence (fitted at 400 °C).
- **Uptake.** LDF with Eq. 6 (above). Only CO₂ adsorbs.
- **Momentum.** Ergun. **Dispersion.** Ruthven (default D_m).
- **Energy.** Isothermal, from the ±3 °C the thesis measured in the 4.0 cm
  breakthrough bed (p. 197).
- **Gas.** Ideal gas.
- **Not modelled.** The sampling line to the mass spectrometer (the thesis
  corrects breakthrough times for its dead time, p. 196), and the reactor's
  heating zones above the bed.

## Validation

### Gas-phase equilibrium (asserted)

Inert sorbent: the outlet of the 4.2 cm bed is at the equilibrium of the rate
model within a residence time. These are code-to-code comparisons. The thesis's
values are its own equilibrium calculations, printed to one decimal.

| Feed | P [bar] | Thesis [%] | Model [%] | Error [pp] |
|---|---|---|---|---|
| Fig. 7 (3 %, S/C 4) | 1.5 | 57.1 | 57.15 | 0.05 |
| | 4.5 | 40.1 | 40.14 | 0.04 |
| | 6.0 | 36.8 | 36.34 | 0.46 |
| Table 6, S/C 3 | 1.5 / 4.5 / 6.0 | 45.3 / 32.1 / 29.2 | 47.25 / 32.71 / 29.55 | 1.95 / 0.61 / 0.35 |
| Table 6, S/C 4 | 1.5 / 4.5 / 6.0 | 51.2 / 36.1 / 31.8 | 52.86 / 36.76 / 33.22 | 1.66 / 0.66 / 1.42 |
| Table 6, S/C 5 | 1.5 / 4.5 / 6.0 | 55.6 / 39.6 / 35.9 | 57.41 / 40.13 / 36.30 | 1.81 / 0.53 / 0.40 |
| Table 6, S/C 6 | 1.5 / 4.5 / 6.0 | 59.4 / 42.4 / 38.5 | 61.21 / 43.04 / 38.97 | 1.81 / 0.64 / 0.47 |

The model is identical to a direct solution of the two equilibria with the
thesis's K_I and K_II. So the remaining differences are in the thesis's
numbers, not in the transport. The Table 6 column lies 0.3–2.0 points below the
equilibrium of Table 6's stated feed at every row. Its basis (feed composition
or K values) is not stated.

Tolerances (regression bounds):
- Fig. 7: **±1 point**, measured 0.46 (6 bar). Covers the three Fig. 7 values.
- Table 6 equilibrium column: **±2.5 points**, measured 1.95 (S/C 3, 1.5 bar).
  Covers its twelve values.

### Sorption-enhanced runs (recorded, not asserted)

The Table 6 runs with the sorbent active. The thesis gives one conversion and
one H₂ concentration per run "during the sorption enhanced reaction" (p. 203),
without a read time or an averaging rule.

| S/C | P [bar] | Thesis X [%] | Model X, 60 s | Model X max (t ≥ 5 s) | Thesis H₂ [%] | Model H₂, 60 s | Model breakthrough [min] |
|---|---|---|---|---|---|---|---|
| 3 | 1.5 | 82.0 | 81.25 | 94.15 | 94.7 | 94.35 | 1.6 |
| 4 | 1.5 | 92.7 | 90.52 | 99.15 | 98.0 | 97.30 | 1.7 |
| 5 | 1.5 | 98.0 | 96.59 | 99.91 | 99.5 | 99.05 | 1.8 |
| 6 | 1.5 | 99.6 | 99.28 | 99.97 | 99.8 | 99.80 | 1.9 |
| 3 | 4.5 | 74.9 | 72.10 | 83.67 | 91.9 | 91.18 | 2.9 |
| 4 | 4.5 | 85.2 | 80.27 | 91.04 | 95.6 | 94.19 | 2.9 |
| 5 | 4.5 | 93.4 | 86.39 | 95.28 | 98.2 | 96.18 | 2.9 |
| 6 | 4.5 | 97.5 | 90.99 | 97.60 | 99.2 | 97.55 | 2.9 |
| 3 | 6.0 | 72.0 | 69.12 | 79.82 | 91.3 | 89.98 | 3.4 |
| 4 | 6.0 | 83.4 | 76.91 | 87.28 | 95.3 | 93.02 | 3.3 |
| 5 | 6.0 | 91.4 | 82.80 | 92.00 | 97.5 | 95.05 | 3.3 |
| 6 | 6.0 | 95.2 | 87.36 | 95.01 | 98.8 | 96.49 | 3.3 |

The thesis gives the breakthrough times of these runs as a range: "4–15 min"
(p. 203). H₂ is on the dry, Ar-free basis. The model's breakthrough is the last
time the outlet CO₂ reaches 5 % of its post-breakthrough value (the thesis's
definition for its breakthrough experiments, p. 194).

Why these rows are not asserted:

1. **No plateau.** In the model the conversion is highest when the start-up gas
   has left (5 s) and falls steadily until CO₂ breaks through. At S/C 3 and 1.5 bar
   it is 94.2 % at 5 s and 81.3 % at 60 s, and after breakthrough it is the
   gas-phase equilibrium (47.3 %): 13 points between 5 s and 60 s. So no single read time stands for the thesis's value. A reading at
   60 s appeared to agree within 1 point at 1.5 bar, but that agreement is an
   artefact of the read time, not a validation.
2. **Early breakthrough.** The model breaks through after 1.6–3.4 min; the
   thesis reports 4–15 min. A stoichiometric estimate agrees with the model.
   The wet Freundlich gives 0.13 mol/kg at the CO₂ partial pressure in front of
   the sorption front (~0.024 bar at 1.5 bar). With 2 g of sorbent and the
   feed's CO₂ production of ~2.2 µmol/s, that is about 2 min.
3. **Maximum conversion too high at low S/C.** At 1.5 bar the maximum is 12.2 and
   6.5 points above the thesis at S/C 3 and 4. At 4.5 and 6 bar the maximum is
   0.2 points below to 8.8 points above the thesis, and the 60 s value is 2.8–8.6
   points below it.

The input that disagrees is the **sorbent capacity at low CO₂ partial
pressure**. The thesis's own data disagree with each other there:

- **Table 3** (p. 197) gives 0.38 mol/kg at 0.05 bar (dry).
- **The Henry constant** the thesis fits over 0.025–0.15 bar (3.243 dry, 3.516 wet
  mol/(kg bar), p. 200) gives 0.16–0.18 mol/kg at 0.05 bar, a factor ~2 lower.
- **The Freundlich fits** lie 22–54 % (dry) and 33 % (wet) below the Table 3
  points at 0.05–0.09 bar. Above 0.14 bar they agree with them within about 15 %.

Capacity is the most sensitive input (below), but a uniform scaling does not
reconcile the table. q × 1.5 moves the 4.5 and 6 bar rows up by 3.0–4.1 points at
60 s and delays breakthrough to 4.3 min or later. But it also moves the 1.5 bar
rows up by 0.6–4.9 points, to 4.1 points above the thesis at S/C 3. The open question and
what would resolve it are in
[skarstrom#245](https://github.com/samaffleck/skarstrom/issues/245). Until then
the rows stay record-only. Breakthrough times and CO+CO₂ levels are not
compared row by row; the thesis prints neither per run in a form the model can
reproduce.

### Attribution runs (final deck)

One input changed at a time, all twelve runs, 2026-11-09:

| Change | X at 60 s, 1.5 bar (S/C 3 → 6) | X at 60 s, 4.5 / 6 bar | Breakthrough [min] |
|---|---|---|---|
| Final deck | 81.25 / 90.52 / 96.59 / 99.28 | 72.1–91.0 / 69.1–87.4 | 1.6–1.9 / 2.9 / 3.3–3.4 |
| 80 cells | +0.04 / +0.05 / +0.04 / +0.02 | ≤ +0.05 | unchanged |
| q × 1.5 (Freundlich k × 1.5) | +4.9 / +4.1 / +2.2 / +0.6 | +3.0–4.1 / +3.2–3.9 | 2.4–2.9 / 4.3–4.4 / 4.9 or later |
| Uptake 1.67× faster (D = D_p, Eq. 6 without its (1 − ε_p) term) | +1.5 / +1.9 / +1.6 / +0.6 | +0.9–1.0 / +0.8–0.9 | +0.1 or unchanged |
| ε_b 0.50 (more sorbent in the bed at the same pellet density) | +4.0 / +3.8 / +2.3 / +0.7 | +3.1–3.6 / +3.3–3.5 | 2.1–2.5 / 3.6–3.7 / 4.1–4.3 |
| ε_b 0.70 (less sorbent) | −5.7 / −6.1 / −5.3 / −3.1 | −4.3–4.7 / −4.1–4.5 | 1.2–1.4 / 2.1–2.2 / 2.5 |
| Start-up without CH₄ (S/C 3, 6 at 1.5 bar) | +0.17 / – / – / +0.10 | – | unchanged |

The bed voidage moves the results mainly through the sorbent mass in the bed:
at a fixed pellet density, ρ_b = (1 − ε_b) ρ_p. The thesis fixes the sorbent mass
(2 g), and ε_b = 0.603 reproduces it with the assumed 2000 kg/m³.

## Open items

- The sorbent's low-pressure capacity (skarstrom#245); then assert the Table 6
  rows.
- The experiment's own start-up (product-free bed) fails the integrator
  (skarstrom#231).
- Not done: the temperature runs (Figs. 5, 6), the cyclic runs (section 9.3.4)
  and the breakthrough experiments (Table 3) as their own case.
