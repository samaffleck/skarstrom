# Ferré et al. (2026) — CO₂/CH₄ co-production PVSA (biogas upgrading)

Validation of the single-bed **4-step pressure–vacuum swing adsorption (PVSA)**
cycle for simultaneous high-purity CO₂ (heavy product) and CH₄ (light product)
recovery from biogas, reported by Ferré et al. (2026).

## Citation

> Ferré, A.; Wedler, C.; Ward, A.; Burger, J.; Pini, R.; Petit, C. (2026),
> *"Co-production of high-purity carbon dioxide and methane via a single-bed
> pressure-vacuum swing adsorption process"*, Separation and Purification
> Technology **406**, 138821. DOI: `10.1016/j.seppur.2026.138821`.


## Model notes

- Detailed dynamic column model is the **Ward & Pini 2022** framework (ref [29]).
- Isotherm: extended (competitive) dual-site Langmuir in **concentration form**
  (`ExtendedDualSiteLangmuirIsotherm_Concentration`; the deck used the
  non-competitive form until the validation change); energies are
  the **internal energy of adsorption ΔU** (not isosteric ΔH), so the fixture
  uses `HeatOfAdsorptionConvention::INTERNAL_ENERGY`.
- Cycle: (1) adsorption → CH₄ raffinate (light product, top); (2) forward
  blowdown (top → pᵢ, not collected); (3) reverse evacuation (bottom → p_L) →
  CO₂ (heavy product); (4) feed pressurisation with feed gas.
- Feed: 40 % CO₂ / 60 % CH₄, 298.15 K, 1 bar. Bed/model parameters: Table 2.
  Optimised operating points and target KPIs: Table 5.

## Files

| File | Purpose |
|------|---------|
| `ferre_2026.cpp` | GoogleTest cases — one `TEST_F` per reproduced operating point |

## Reproduced configurations

Each case reproduces the Table 5 **recovery–recovery, max Re_CO₂** operating point
for one of the paper's four top-performing adsorbents (a shared parameterised
`RunFerreCase` helper drives them). Isotherm parameters from SI Tables S1/S2,
densities from S3 (`ρ_b = ρ_ads·(1−ε)`), operating points from Table 5.

| Test (`Ferre2026Test.*`) | Adsorbent | Paper targets (Pu_CO₂ / Re_CO₂ / Pu_CH₄ / Re_CH₄, %) |
|--------------------------|-----------|------------------------------------------------------|
| `NaY_F`      | NaY (Feng)   | 95.1 / 97.8 / 98.8 / 66.0 |
| `NaY_W`      | NaY (Wedler) | 95.0 / 95.9 / 97.6 / 65.4 |
| `Zeolite13X` | 13X          | 95.0 / 95.9 / 98.0 / 48.6 |
| `MIL160Al`   | MIL-160(Al)  | 95.0 / 94.1 / 97.6 / 35.8 |

_Other paper materials (LiX, 5A, AC, ZIF-8, ZIF-71) were screened out and have no
Table 5 operating point to validate against; their isotherm fits are in the SI._

## Validation

The tests compare the simulation with the paper using the helpers in
`../../fixtures/validation_metrics.h`; the published numbers they use are in
[`reference_values.csv`](reference_values.csv) (each with its table / figure /
page). Every tolerance is the measured error rounded up, with the measured
value in a comment next to it. Results of the last run:
[validation page](../../../../docs/validation.md#ferre-2026)
(generated from `../validation_summary.json`).

- Product purities agree with Table 5 within 0.1-1.8 percentage points and CO2 recovery within 0.6-3 pp; CH4 recovery is 5-7 pp low for the zeolites and 11 pp high for MIL-160(Al), whose CH4 productivity is 59 % high.
- The decks now use the extended (competitive) DSL isotherm the paper states (Sect. 3.1); they previously used the non-competitive form. adsim uses a constant LDF coefficient and a single heat of adsorption per component instead of the paper's loading-dependent forms.
- Specific energy is not compared: the decks drive the pressure steps with pressure boundary conditions, not pump models.

### Gas heat capacity as in the paper (2026-11)

The paper uses a constant 30.70 J/(mol K) for the gas and the adsorbed phase. The deck now uses the constant gas heat capacity mode with 30.7 J/(mol K) for CO₂ and CH₄ (before: the default ideal-gas C_p0 at the local temperature; the adsorbed phase was already 30.7). Changes: MIL-160(Al) CO₂ purity error 2.152 → 2.186 pp (bound 2.2), CO₂ recovery 2.287 → 2.333 pp; NaY(F) CO₂ purity 0.392 → 0.335 pp, CH₄ recovery 3.62 → 3.53 pp; NaY(W) CO₂ purity 1.265 → 1.214 pp; 13X CH₄ purity 0.005 → 0.014 pp (bound 0.076), CO₂ recovery 0.244 → 0.261 pp. No bound changed. All papers' settings: [../README.md](../README.md#gas-heat-capacity-per-paper).

### Grid and advection scheme (2026-11-04)

The decks ran first-order upwind advection on 30 cells. A grid study of
`Ferre2026_NaYF` (10-160 cells, Richardson extrapolation;
[Numerics](../../../../docs/numerics.md#advection-schemes)) found that the CO2
recovery needs about 610 upwind cells, or 110 WENO cells, to be within 0.1 pp
of its grid limit. All four decks now run WENO, the engine's default, on 120
cells (CSS after 141-159 cycles).

Before: upwind, 30 cells; after: WENO, 120 cells, both at CSS (bold: changed
bound). CO2 recovery rose by 0.4-0.8 pp and CH4 purity by 0.2-0.4 pp; that
brings NaY(F), 13X and MIL-160(Al) closer to Table 5 (bounds tightened where
the error fell well below them) and moves NaY(W)'s CH4 purity 0.3 pp further
above the paper (bound 0.45 → 0.65 pp). The CH4 recovery, the largest
disagreement, hardly moves (≤ 0.26 pp).

| Case | Quantity | Reference | Before | After | Error before → after | Bound before → after |
|---|---|---|---|---|---|---|
| MIL160Al | CO2 purity (heavy product) | 95.00 | 97.02 | 97.12 | 2.02 → 2.12 | 2.2 |
| MIL160Al | CO2 recovery | 94.10 | 91.13 | 91.82 | 2.97 → 2.28 | 3.5 |
| MIL160Al | CH4 purity (light product) | 97.60 | 96.43 | 96.83 | 1.17 → 0.77 | 1.4 |
| MIL160Al | CH4 recovery | 35.80 | 47.14 | 47.23 | 11.34 → 11.43 | 13 |
| MIL160Al | CH4 productivity | 0.9 | 1.43 | 1.43 | 0.5938 → 0.5932 | 0.68 |
| NaY_F | CO2 purity (heavy product) | 95.10 | 94.64 | 94.65 | 0.46 → 0.45 | 1.1 |
| NaY_F | CO2 recovery | 97.80 | 97.18 | 97.60 | 0.62 → 0.20 | **1.3 → 0.25** |
| NaY_F | CH4 purity (light product) | 98.80 | 98.54 | 98.75 | 0.26 → 0.05 | **0.59 → 0.097** |
| NaY_F | CH4 recovery | 66.00 | 62.25 | 62.21 | 3.75 → 3.79 | 5.6 |
| NaY_F | CH4 productivity | 3.20 | 3.10 | 3.10 | 0.03135 → 0.03086 | 0.046 |
| NaY_W | CO2 purity (heavy product) | 95.00 | 93.73 | 93.70 | 1.27 → 1.30 | 1.8 |
| NaY_W | CO2 recovery | 95.90 | 95.83 | 96.46 | 0.07 → 0.56 | 0.89 |
| NaY_W | CH4 purity (light product) | 97.60 | 97.82 | 98.14 | 0.22 → 0.54 | **0.45 → 0.65** |
| NaY_W | CH4 recovery | 65.40 | 59.85 | 59.77 | 5.55 → 5.63 | 7.6 |
| NaY_W | CH4 productivity | 3.10 | 2.89 | 2.89 | 0.06667 → 0.06812 | 0.085 |
| Zeolite13X | CO2 purity (heavy product) | 95.00 | 93.80 | 93.67 | 1.20 → 1.33 | 1.5 |
| Zeolite13X | CO2 recovery | 95.90 | 94.85 | 95.62 | 1.05 → 0.28 | **1.6 → 0.35** |
| Zeolite13X | CH4 purity (light product) | 98.00 | 97.64 | 97.97 | 0.36 → 0.03 | **0.73 → 0.076** |
| Zeolite13X | CH4 recovery | 48.60 | 43.88 | 43.61 | 4.72 → 4.99 | 5.7 |
| Zeolite13X | CH4 productivity | 1.30 | 1.19 | 1.19 | 0.08467 → 0.08582 | 0.11 |
