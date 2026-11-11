# Nguyen et al. (2025) — humid CO₂/N₂ VSA on CALF-20: the dry cases A0 and B0

## Citation

> T.T.T. Nguyen, G.K.H. Shimizu, A. Rajendran, *Experimental demonstration of
> humid post-combustion CO₂ capture by vacuum swing adsorption using CALF-20*,
> Adsorption 31 (2025) 45. DOI: [10.1007/s10450-025-00600-z](https://doi.org/10.1007/s10450-025-00600-z).
> Publisher correction: [10.1007/s10450-025-00628-1](https://doi.org/10.1007/s10450-025-00628-1)
> (it adds the missing panels b of Figs. 3–9 and a reference number; no table changes).

The folder is named `nguyen_2024` after the year of submission (the paper was
received in September 2024); the example decks are `Nguyen2024_A0` and
`Nguyen2024_B0`.

The paper is experimental: it runs a two-column, four-step VSA cycle with light
product pressurisation (LPP) on the rig, columns and structured CALF-20 of
Nguyen et al. (2023) ([`../nguyen_2023/`](../nguyen_2023/README.md)), with a dry
and a humidified 15 % CO₂ / N₂ feed, and reports the measured purity, recovery
and productivity (Table 1). It has no process model of its own.

## Files

| File | Purpose |
|------|---------|
| `nguyen_2024.cpp` | `Nguyen2024Test.CaseA0`, `CaseB0`; `DISABLED_Attribution` (the one-factor runs below) |
| `reference_values.csv` | every Table 1 and SI Table A.1 number used, with table and page |

The fixture derives from `Nguyen2023Test` (`../nguyen_2023/nguyen_2023.h`) and
builds the `Nguyen2023_Exp_VII` deck structure (`RunNguyen2023LPP` in
`../nguyen_2023/nguyen_2023.cpp`) at this paper's conditions.

## Reproduced cases

| Test (`Nguyen2024Test.*`) | Table 1 row | t_ads [s] | v_feed [m/s] | P_I [bar] | P_L [bar] | Example deck |
|---|---|---|---|---|---|---|
| `CaseA0` | Case 1, 0 % RH (A0) | 180 | 0.081 | 0.119 | 0.032 | `Nguyen2024_A0` |
| `CaseB0` | Case 2, 0 % RH (B0) | 157 | 0.024 | 0.046 | 0.031 | `Nguyen2024_B0` |

B0 is the operating point of Nguyen 2023 Exp. VII (that paper measured
94.3 / 88.4 % and, in a repeat, 93.8 / 87.5 %; this paper 92.4 / 83.3 %).

The humid cases A1–A4 and B1–B4 are not reproduced (open item 1).

## Validation protocol

### 1. Source

- **Partial.** The supporting information (SI, open access) and the
  publisher correction were read in full. The main text (paywalled) was not
  available here; the Table 1 values, the KPI definitions (p. 6, Eqs. 2–4), the
  cycle description (Fig. 2, pp. 4–5), P_H = 0.96 bar (p. 5) and the
  thermocouple positions (p. 3) were transcribed from it by the maintainers
  and are cited by table and page. Step 1 closes when the main text is read.
- SI Table A.1 (column and transport parameters) was checked against the SI
  directly.

### 2. Provenance audit

The deck is the `Nguyen2023_Exp_VII` deck; only the rows marked "this paper"
were changed. Deviations from this paper's SI are in bold.

| Value | Deck parameter (registry path) | Source (table / equation / page / sentence) | Note (derived arithmetic, assumption status) |
|---|---|---|---|
| P_H = 0.96 bar | raffinate vents of `Feed 1` and `Feed 2`, LPP end pressure, `initial_conditions.pressure` | this paper: Table 1 caption; p. 5 "atmospheric pressure P_H = 0.96 bar" | VII deck: 1.0 bar |
| T = 295.15 K | feed and LPP stream temperatures, `initial_conditions.temperature`, `initial_conditions.wall_temperature`, `bed.ambient_temperature` | this paper: Table 1 caption (22 °C) | VII deck: 293 K (an assumption there) |
| Feed 15 % CO₂ / 85 % N₂ | `endpoint[...].y[CO2]` of both feed inlets | this paper p. 2 (premixed cylinder) | dry cases only |
| t_ads, v_feed, P_I, P_L | `step[Feed 2].time` = t_ads − 50 s; feed velocity; blowdown and evacuation end pressures; LPP start pressure | this paper: Table 1 | v_feed read as interstitial: superficial = 0.4 v_feed, as in the Nguyen 2023 decks (see "Feed velocity" below) |
| Step times: Feed 1 50 s, blowdown 60 s, evacuation 40 s, pressurisation 50 s | `step[...].time` | VII deck | **the paper ends blowdown, evacuation and pressurisation on pressure triggers and prints no step times** (open item 2, #271) |
| Blowdown and evacuation profiles | exponential to P_I and P_L, F = 100 (99 % of the change at the step end: 0.077 and 0.115 s⁻¹) | VII deck | **SI Table A.1: α_BD = 0.153 s⁻¹, α_EVAC = 0.075 s⁻¹, fitted to the measured pressure histories**; attribution run (1) |
| Column inner / outer diameter 0.021 / 0.025 m | `bed.diameter`, `wall.thickness` 2 mm | SI Table A.1 (measured) | |
| Column length | `layer[1].length` = 0.3039 m | **SI Table A.1: 0.267 m (measured; the text p. 2 says 26.5 cm)** | derived in the VII deck from m_ads, ρ_p and ε = 0.4: L = 0.036 / (570 · 0.6 · π/4 · 0.021²) = 0.3039 m, so the 36.0 g inventory is kept; 36.0 g in 0.267 m needs ε = 1 − 0.036 / (570 · 9.248e-5 m³) = 0.317; attribution run (2) |
| Bed voidage 0.4 | `layer[1].bed_void_fraction` | VII deck (Nguyen 2023) | not in SI Table A.1 |
| Adsorbent mass 36.0 g, particle density 570 kg/m³ | `adsorbent[1].solid_density` = 570 | SI Table A.1 (measured) | adsorbent volume m/ρ_p = 63.2 cm³ (productivity basis) |
| Particle voidage 0.35, tortuosity 3, pellet diameter 1 mm | `adsorbent[1].particle_void_fraction`, `kinetics[C].*`, `layer[1].particle_diameter` | VII deck (Nguyen 2023) | not in SI Table A.1 |
| Isotherm: competitive DSL (CO₂ and N₂) | `adsorbent[1].isotherm[C].*` | VII deck (Nguyen 2023, SI Table S1) | this paper prints no isotherm for CO₂ or N₂ |
| Heats of adsorption 38.0 / 15.0 kJ/mol | `adsorbent[1].heat_of_adsorption[C]` | VII deck | |
| LDF: macropore, D_m from the state | `adsorbent[1].kinetics[C].*` | VII deck | |
| Axial dispersion: Ruthven, D_m = 1.6e-5 m²/s | `layer[1]` dispersion | VII deck | **SI Table A.1: D_m = 1.496e-5 m²/s (printed unit "m s⁻²")**; attribution run (3) |
| Viscosity | polynomial mixture value (CoolProp fit) | VII deck | **SI Table A.1: 1.812e-5 kg/(m s), constant**; attribution run (3) |
| Pressure drop | Kozeny–Carman (180) | VII deck | **the Rajendran group's model uses Blake–Kozeny (150) (Nguyen 2023 and Estupiñan Perez 2019); this paper prints no form**; attribution run (4) |
| Gas C_p N₂ 1040, CO₂ 841 J/(kg K), constant; adsorbed phase = gas | `fluid.gas_heat_capacity[C]` 29.13 / 37.01 J/(mol K), `fluid.adsorbed_heat_capacity_mode` = gas_phase | SI Table A.1 (NIST; adsorbed "assumed") | H₂O 4165 J/(kg K) is a liquid value; not used here |
| Adsorbent C_p 1350 J/(kg K) | `adsorbent[1].heat_capacity` | SI Table A.1 (fitted) | |
| Wall ρ 7800 kg/m³, C_p 502 J/(kg K), k 10 W/(m K) | `wall.*` | SI Table A.1 (standard) | |
| K_z = 0.4 W/(m K) | `layer[1]` heat dispersion (constant) | SI Table A.1 (fitted) | |
| h_in 8, h_out 10 W/(m² K) | `layer[1]` internal / external heat transfer | SI Table A.1 (fitted) | |
| LPP connection, N₂ vents | `Feed 1` light-product stream connected to the pressurisation inlet; vents back-feed N₂ | VII deck (conserving since #287) | unchanged |
| KPIs | feed = `Feed 1` + `Feed 2` bottom; product = evacuation bottom | this paper Eqs. 2–4 (p. 6) | purity = CO₂ / (CO₂ + N₂) of the evacuation product; recovery = CO₂ in it / CO₂ fed in the adsorption steps; productivity = CO₂ in it / (m_ads / ρ_p · t_cycle). The engine's own KPI records give the same purity and recovery |
| CSS | component mass balance, ≥ 5 cycles; caps 50 (A0) and 430 (B0) | VII deck | A0 stops after 39 cycles, B0 after 326 |
| Grid | WENO, 80 cells | VII deck (grid-converged, Nguyen 2023 README) | |

**Feed velocity.** Table 1 prints a "feed velocity" without saying whether it
is superficial or interstitial. The Nguyen 2023 decks read it as interstitial,
and the decks here do the same (superficial = ε v with ε = 0.4). With the
printed recovery and the 307 s cycle of B0, mass conservation gives a
productivity of 0.132 mol/(m³ s), 10 % below the printed 0.147; with the
superficial reading it would be 0.33, so the interstitial reading is the
consistent one. For A0 the same arithmetic gives 0.405 at the 330 s cycle
against the printed 0.577: that value needs a 232 s cycle or a feed velocity of
0.115 m/s. A0's printed productivity is inconsistent with its own cycle time
under either reading of v_feed, so it is recorded, not asserted. (The paper
ends the evacuation, blowdown and pressurisation steps on pressure triggers,
so its real cycle may be shorter than the deck's; that is open item 2.)

### 3. Model-form audit (against the VII deck)

The paper has no model, so the audit is of the deck against the experiment as
described:

- **Cycle.** Two columns in the experiment; one column in the deck, which
  emulates the partner column by connecting the Feed 1 light-product draw to
  its own pressurisation step (the Nguyen 2023 single-bed emulation). The
  connection conserves moles (#287).
- **Pressure histories.** The experiment switches steps on pressure triggers;
  the deck uses fixed step times with an exponential approach to P_I and P_L.
  The SI's fitted decay rates differ from the deck's (run (1)).
- **Column.** The deck keeps the 36.0 g inventory in a 0.304 m column (ε 0.4)
  rather than the measured 0.267 m (run (2)).
- **Transport.** D_m and μ come from the state rather than the SI's constants
  (run (3)); the pressure drop is Kozeny–Carman, not Blake–Kozeny (run (4)).
- **Not modelled.** Dead volumes and piping, the valve hold times, and the
  humid feed.

### 4. Data

No curves are compared. The paper's Fig. 3 (temperatures at 1, 7.8 and 19 cm),
Fig. 4 (outlet CO₂ mole fractions) and SI Figs. A.2–A.11 (pressure histories)
have not been digitised; no measured data are published ("available on
request").

### 5. Metrics

Purity and recovery as absolute errors in percentage points, productivity as
a relative error, against Table 1 (experiment).

### 6. Attribution

One change at a time from the deck, both cases
(`Nguyen2024Test.DISABLED_Attribution`; run with
`--gtest_also_run_disabled_tests --gtest_filter='Nguyen2024Test.DISABLED_*'`
from `core/`).

ATTRIBUTION_TABLE

### 7. Measured errors

| Case | Quantity | Paper (± printed) | Engine | Error | Bound |
|---|---|---|---|---|---|
| A0 | CO₂ purity | 98.3 ± 0.88 % | 94.89 % | −3.41 pp | info |
| A0 | CO₂ recovery | 71.2 ± 1.46 % | 81.98 % | +10.78 pp | info |
| A0 | CO₂ productivity | 0.577 ± 0.014 mol/(m³ s) | 0.467 | −19.1 % | info |
| B0 | CO₂ purity | 92.4 ± 1.31 % | 93.33 % | +0.93 pp | 1.0 pp |
| B0 | CO₂ recovery | 83.3 ± 1.97 % | 94.87 % | +11.57 pp | info |
| B0 | CO₂ productivity | 0.147 ± 0.005 mol/(m³ s) | 0.1500 | +2.0 % | 2.5 % |

Cycles to CSS and run time (one core, Linux release build): A0 39 cycles,
21 s; B0 326 cycles, 100 s.

**Reading the agreement.** B0's purity and productivity reproduce within the
printed uncertainty, and they are asserted. B0's productivity agreement is
partly a compensation: the engine recovers 11.6 points more CO₂ than measured
(14 % more product), which offsets the 10 % gap between the printed
productivity and the one the printed recovery and the deck's cycle time give
(see "Feed velocity"). A0's purity is 3.4 points low, outside the printed
uncertainty, and is recorded.

**Recovery.** Both cases over-recover by about 11 points. This is the known
LPP over-recovery of the Nguyen 2023 LPP decks (Exp. VII: 95.25 % against
88.4 %), and it is larger here because the measured recoveries of this paper
are lower than those of Nguyen 2023 at the same operating point (B0 83.3 %
against VII 88.4 %). None of the printed
parameters closes it (attribution); in the model the CO₂ front stays inside
the bed, so the raffinate and the blowdown gas carry less CO₂ than the
outlet mole fractions of the paper's Fig. 4 show.

## Open items

1. **Humid cases (A1–A4, B1–B4).** The paper prints no water isotherm for
   CALF-20; it is in its reference [18], Nguyen et al., Ind. Eng. Chem. Res. 63
   (2024) 3265–3281 (competitive CO₂/H₂O adsorption on CALF-20). The engine has
   no competitive CO₂/H₂O isotherm kind for a physisorbent yet (#272).
2. **Pressure-triggered steps.** The experiment ends blowdown, evacuation and
   pressurisation on pressure triggers and prints no step times; the deck
   keeps the VII deck's 60 / 40 / 50 s (#271). A measured cycle time would
   also settle A0's productivity.
3. **Recovery gap** of about 11 points in both cases (see above), shared with
   the Nguyen 2023 LPP cases.
4. **Main text** not read here (step 1).
5. Temperatures (Fig. 3) and outlet compositions (Fig. 4) are not compared (no
   digitised or published data).
