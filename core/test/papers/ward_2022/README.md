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

## Key parameters

- Adsorbent: activated carbon, lab-scale test rig (Tables 1, 2, 6). Setup is
  inline per test.
- Ambient temperature 291.4 K (Table 1; the decks had 293 K, the feed
  temperature, until the validation change).
- The paper's experiments use He as the carrier; the decks use a
  non-adsorbing "Air" component.
- Digitised data in `core/data/ward_2022/Qin=<flow>/`.
- `MTC_Invariance` runs the Q=80 cm³/min setup at high vs. low MTC and asserts
  the accumulated CO₂ at equilibrium matches.

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
