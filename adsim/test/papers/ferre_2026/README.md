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
