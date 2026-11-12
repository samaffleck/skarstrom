# Nguyen et al. (2022) — CO₂ breakthrough on CALF-20

## Citation

> Nguyen, T.T.T., Lin, J.-B., Shimizu, G.K.H., Rajendran, A. *Separation of CO₂
> and N₂ on a hydrophobic metal organic framework CALF-20.* Chemical Engineering
> Journal 442:136263 (2022).
> DOI: [10.1016/j.cej.2022.136263](https://doi.org/10.1016/j.cej.2022.136263)

The experimental conditions of every breakthrough run are in Table A2.1 of the
first author's thesis: T.T.T. Nguyen, *Development of adsorption processes for
post-combustion CO₂ capture using CALF-20* (PhD thesis, University of Alberta).
Chapter 2 of the thesis is the paper. Both were read in full (2026-11-12); the
supporting information was not available, and the thesis's appendix holds its
tables (A2.1, A2.2).

This is a different publication from `../nguyen_2023` (the VSA-cycle paper,
DOI 10.1016/j.cej.2022.139550), which cites it as its ref. [10].

## Files

| File | Purpose |
|------|---------|
| `nguyen_2022.cpp` | `Nguyen2022Test`: the six breakthrough experiments, each saved as an example deck |
| `nguyen_2022.h`   | `Nguyen2022Fixture`: the column, CALF-20, the transport and the feed of one case |

The measured curves and the paper's model curves are in
[`core/data/nguyen_2022/`](../../../data/nguyen_2022/) (its README records how they were digitised).

## Reproduced cases

Fig. 4 is single-component CO₂ (balance helium) at 50 cm³/min; Fig. 5 is
CO₂/N₂ (no helium) at 100 cm³/min. The column starts full of helium (Sect. 3.2).
Conditions: thesis Table A2.1.

| Test (`Nguyen2022Test.*`) | Example deck | Feed | Flow [cm³/min] | T [°C] | P [bar] |
|---|---|---|---|---|---|
| `Fig4_100pct` | `Nguyen2022_Fig4_100pct` | 100 % CO₂ | 50 | 22.4 | 0.97 |
| `Fig4_70pct` | `Nguyen2022_Fig4` | 70 % CO₂, 30 % He | 50 | 22.3 | 0.97 |
| `Fig4_50pct` | `Nguyen2022_Fig4_50pct` | 50 % CO₂, 50 % He | 50 | 22.4 | 0.97 |
| `Fig5_30pct` | `Nguyen2022_Fig5_30pct` | 30 % CO₂, 70 % N₂ | 100 | 21.8 | 0.98 |
| `Fig5_15pct` | `Nguyen2022_Fig5` | 15 % CO₂, 85 % N₂ | 100 | 20.3 | 0.97 |
| `Fig5_5pct` | `Nguyen2022_Fig5_5pct` | 5 % CO₂, 95 % N₂ | 100 | 23.4 | 0.97 |

## Provenance (2026-11-12, from the paper and the thesis)

| Deck value | Value | Source | Note |
|---|---|---|---|
| Column inner / outer diameter | 0.0282 / 0.0318 m | Table 2 | wall 1.8 mm |
| Bed length | **0.0782 m** | derived: m / (ρ_p (1 − ε) A) with the 16.7 g of Table 2 | Table 2 gives 0.0786 m; the Table 2 values are 0.5 % inconsistent, and the deck keeps the mass |
| Bed voidage, particle density | 0.4, 570 kg/m³ | Table 2 | |
| Particle voidage, tortuosity | 0.35, 3 | Table 2 ("assumed") | used in the LDF coefficient only |
| Gas holdup | bed voidage only | Eqs. 1–3 | the adsorbent's particle voidage is 0 in the deck |
| Particle diameter | **1.5 mm (assumption)** | Sect. 2: pellets 1.0–3.0 mm | not in Table 2 or the thesis; enters the LDF coefficient, the dispersion and the pressure drop |
| Isotherm | extended dual-site Langmuir, concentration form, equal-energy-site N₂ | Table 3, Eq. 20 | |
| LDF coefficient | macropore with the capacity factor, D_m 1.496 × 10⁻⁵ m²/s | Eq. 5, Table 2 | engine preset `macropore_capacity` (k = 15 ε_p D_m c / (τ r_p² ρ_p q*)) |
| Heats of adsorption | 38, 15 kJ/mol | Fig. 2d (isosteric, about 37.5 and 15 kJ/mol) | read from the figure; Table 2 does not list them |
| Gas heat capacity | constant: CO₂ 841, N₂ 1040, He 5200 J/(kg K) per mole | Table 2 (NIST) | |
| Adsorbed-phase heat capacity | equal to the gas's | Table 2 ("assumed") | |
| Adsorbent heat capacity | 1371 J/(kg K) | Table 2 (fitted) | |
| Wall | 7800 kg/m³, 502 J/(kg K), 16 W/(m K) | Table 2 | |
| Wall-end condition | **zero flux (assumption)** | not stated | ends at ambient: see the attribution |
| h_in, h_out | 6.963, 14 W/(m² K) | Table 2 (fitted) | |
| Axial conduction K_z | 0.4 W/(m K) | Table 2 (fitted) | |
| Axial dispersion | D_L = 0.7 D_m + v r_p, D_m 1.496 × 10⁻⁵ | Eq. 6 | the engine's Ruthven form |
| Pressure drop | Blake–Kozeny (Darcy, 150) | Eq. 7 | |
| Cells | 30 | Sect. 4.1 | grid-converged (below) |
| Flow, T, P, feed gas | per case (table above) | thesis Table A2.1 | the flow is taken at the feed conditions |
| Thermocouple | 1.82 cm from the outlet (z = L − 0.0182 m) | Sect. 3.2 | the figure captions say z = 0.75 L |

Bold entries are assumptions or differ from the paper.

## Validation

### The data were the paper's model curves (2026-11-12)

The curves in `core/data/nguyen_2022/Fig{4,5}/` that the tests used until
2026-11-12 trace the paper's **simulation lines**, not its measured markers:
overlaid on the 300 dpi figures they follow the solid lines exactly, also where
the markers lie far from them (the 50 % CO₂ front of Fig. 4). The tests therefore
compared Skarstrom with Nguyen's model, under the label "experiment". Those files
are now in `paper_model/` and give the "vs paper model" (code-to-code) metrics.
The measured markers were re-digitised from the figures (`Y_CO2`, `T_CO2` and,
new, the outlet `Y_N2` of Fig. 5), with the time axis converted from the paper's
dimensionless time t̄ = t v / L with each run's flow; the data README records
the method and its limits.

### Results

Every case is a prediction with the paper's values; nothing is fitted. Signed
front-time errors (simulation minus measurement); thermocouple RMSE in K; the
N₂ curve as an RMSE in mole fraction.

| Case | CO₂ front vs measured | y_CO₂ NRMSE | Thermocouple RMSE [K] | y_N₂ RMSE | CO₂ front vs paper model |
|---|---|---|---|---|---|
| Fig4_100pct | +6.5 % | 0.29 | 8.0 | – | −1.8 % |
| Fig4_70pct | +3.2 % | 0.12 | 5.2 | – | −1.7 % |
| Fig4_50pct | +10.9 % | 0.19 | 6.3 | – | −3.9 % |
| Fig5_30pct | −7.2 % | 0.15 | 2.7 | 0.28 | −0.1 % |
| Fig5_15pct | +5.6 % | 0.13 | 1.0 | 0.14 | −4.0 % |
| Fig5_5pct | −9.7 % | 0.09 | 1.7 | 0.036 | −4.1 % |

The paper's own model against the same measured markers (its digitised lines,
the same metrics): CO₂ fronts +7.6, +5.0, +15.4, −9.6, +10.0, −5.8 %; y_CO₂
NRMSE 0.10–0.31; thermocouple RMSE 1.5–5.3 K, with the temperature peaks within
0.1–1.9 K. So Skarstrom places the fronts as well as the published model does
(closer in five of six cases), and both miss the measured fronts by up to 10–15 %
with mixed signs. Skarstrom's Fig. 4 temperature peaks are 6–8 K above the
measurement and 5–6 K above the paper's model; its Fig. 5 temperatures are as
close as the paper's (peaks within 1.5 K).

Tolerances: the measured error × 1.1, rounded up, with a floor of 0.2 % for
the front times (`nguyen_2022.cpp`).

### Attribution (Python, one change at a time, 2026-11-12)

From the deck shipped before this audit (CO₂ in helium for Fig. 5 too, 295 K for
every case, the engine's ideal-gas heat capacities, the macropore LDF without the
capacity factor and the mixture's D_m, Ruthven with D_m 1.6 × 10⁻⁵,
Kozeny–Carman, the particle voidage in the gas holdup), against the measured
markers (front error; thermocouple RMSE):

| Change | Fig. 4 (100 / 70 / 50 %) | Fig. 5 (30 / 15 / 5 %) |
|---|---|---|
| deck before the audit | +6.9 / +4.2 / +11.7 %; 7.3 / 5.3 / 6.1 K | −9.4 / +4.7 / −5.7 %; 2.7 / 1.8 / 0.8 K; y_N₂ NRMSE 0.72 / 0.85 / 0.95 |
| CO₂ in N₂ for Fig. 5 (Sect. 6, Table A2.1) | – | −9.4 / +4.2 / −7.1 %; y_N₂ 0.34 / 0.17 / 0.035 |
| per-case temperature (Table A2.1) | +6.7 / +4.0 / +11.4 % | −9.3 / +6.1 / −8.3 % |
| paper gas heat capacity | < 0.1 % | < 0.1 % |
| capacity-factor LDF, D_m 1.496 × 10⁻⁵ | < 0.1 %, curves within 0.01 NRMSE | < 0.4 % |
| no gas holdup in the pores (on the paper's other forms) | −1.1 / −0.8 / −0.6 % | < 0.3 % |
| bed length 0.0786 m (Table 2) | +0.7 / +0.6 / +0.7 % | +0.6 / +0.6 / +0.5 % |

Tested and not adopted (all on the paper's forms):

- **Wall ends held at ambient:** Fig. 4 peaks drop to +1.4–4.2 K above the
  measurement but the fronts move 4–5 % later, and Fig. 5 peaks fall 2–5 K below
  it.
- **CO₂ capacity × 0.97** (the breakthrough loadings of Table A2.1 are 2–5 %
  below the isotherm for Fig. 4 and within ±2.4 % for Fig. 5): Fig. 4 fronts
  improve by 2–3 %, the 30 % and 5 % Fig. 5 fronts get 2.5 % worse.
- **A calibration of C_p,s, h_in, h_out and K_z on the Fig. 4 temperatures** (the
  paper's own procedure, Sect. 5.4): Fig. 4 thermocouple RMSE 2.5–4.0 K, but the
  predicted 30 % and 15 % Fig. 5 peaks fall 5–7 K below the measurement. A joint fit of the
  capacity, h_in, h_out and the particle diameter on Fig. 4 leaves the Fig. 5
  fronts at −11 / +3 / −12 %. No calibration is shipped.
- **Unit slips in the paper's heat capacities** (the per-kg values of Table 2 used
  per mole in the adsorbed phase, or in the gas): neither reproduces the paper's
  model temperatures.
- **Heat of adsorption 35.06 kJ/mol** (−ΔU_b): the Fig. 4 peaks drop by about 2.4 K, the Fig. 5 peaks by 1–2 K (too low then), and the fronts move up to 2 % later.
- **Thermocouple at 0.75 L** (the figure captions): the peaks move by about 0.3 K.

### Open items

- The 6–8 K excess of the Fig. 4 temperature peaks (CO₂/He). With the same
  inputs the paper's model is 5–6 K cooler; the cause is not found.
- The run-to-run scatter of the measured fronts (±10 %, mixed signs in Fig. 5)
  that neither model captures. The curves are in dimensionless time, so a flow
  error mostly cancels in them.
- The particle diameter and the wall-end condition are not printed.

### Grid and advection scheme (2026-11-04)

History: this study compared the decks of the time with the curves then in `core/data/nguyen_2022/Fig{4,5}/`, which are the paper's model lines (now `paper_model/`), not measurements, and its case names are the old ones. The grid result holds: 30 cells are converged.

The decks keep their 30 cells and now run WENO, the engine's default, instead
of first-order upwind. A grid study (10-160 cells, Richardson extrapolation;
[Numerics](../../../../docs/numerics.md#advection-schemes)) found the CO2
half-feed time within 0.5 % of its grid limit on 10 (Fig. 4) and 21 (Fig. 5)
WENO cells (11 and 26 upwind), so 30 cells are converged. The outlet curves
are sharper: every NRMSE rose by 0.001-0.017, the front times moved by less
than 0.3 % (19 s of 7400 for the 5 % feed, whose error fell from 0.31 % to
0.06 %; that bound was tightened). Before: upwind; after: WENO; both on 30
cells (bold: changed bound).

| Case | Quantity | Reference | Before | After | Error before → after | Bound before → after |
|---|---|---|---|---|---|---|
| Fig4_BinaryBreakthrough | outlet y_CO2 | 0.7 | 0.05689 | 0.06231 | 0.08128 → 0.08901 | 0.1 |
| Fig4_BinaryBreakthrough | CO2 front time (y_CO2 = y_feed/2) | 1747.30 | 1730.66 | 1733.98 | 0.009525 → 0.007626 | 0.027 |
| Fig4_Binary_100pct | outlet y_CO2 | 1.00 | 0.1021 | 0.1053 | 0.1021 → 0.1053 | 0.12 |
| Fig4_Binary_100pct | CO2 front time (y_CO2 = y_feed/2) | 1258.45 | 1248.50 | 1250.91 | 0.007906 → 0.00599 | 0.017 |
| Fig4_Binary_50pct | outlet y_CO2 | 0.5 | 0.04642 | 0.0486 | 0.09283 → 0.0972 | 0.14 |
| Fig4_Binary_50pct | CO2 front time (y_CO2 = y_feed/2) | 2441.30 | 2360.61 | 2364.22 | 0.03305 → 0.03158 | 0.06 |
| Fig5_TernaryBreakthrough | outlet y_CO2 | 0.15 | 0.01132 | 0.01275 | 0.0755 → 0.08497 | 0.097 |
| Fig5_TernaryBreakthrough | CO2 front time (y_CO2 = y_feed/2) | 3364.00 | 3206.85 | 3202.04 | 0.04672 → 0.04814 | 0.059 |
| Fig5_Ternary_30pct | outlet y_CO2 | 0.3 | 0.02911 | 0.03033 | 0.09704 → 0.1011 | 0.13 |
| Fig5_Ternary_30pct | CO2 front time (y_CO2 = y_feed/2) | 1708.31 | 1727.23 | 1727.87 | 0.01108 → 0.01145 | 0.039 |
| Fig5_Ternary_5pct | outlet y_CO2 | 0.05 | 0.006039 | 0.006875 | 0.1208 → 0.1375 | 0.15 |
| Fig5_Ternary_5pct | CO2 front time (y_CO2 = y_feed/2) | 7397.26 | 7420.14 | 7401.34 | 0.003093 → 0.0005517 | **0.026 → 0.0026** |
