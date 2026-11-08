### Added — humid DAC breakthrough example
- Example deck `Illustrative_Humid_DAC_Breakthrough`: a non-isothermal breakthrough of 400 ppm CO₂ in humid
  N₂ (y_H₂O = 0.0115, 49 % RH at 20 °C, 1 bar) on the APDES-NFC library adsorbent, with the GAB water
  isotherm, the humid Toth CO₂ isotherm that reads the adsorbed water loading, LDF kinetics and the water
  heat of sorption. An integration test checks that the run is healthy and the mole balances of all three
  components close, that water breaks through first and the outlet loadings reach the GAB and humid Toth
  values at the local state, that a dry feed reaches the lower dry CO₂ capacity, and that the temperature
  peaks before the CO₂ front.
- `docs/isotherms.md` evaluates the APDES-NFC dry and humid loadings at the feed state; the TSA tutorial
  (`docs/tsa-regeneration.md`) has a "Humid feeds" section.
