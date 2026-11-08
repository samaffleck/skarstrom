### Added — water and humid CO2 parameters for three more library entries
- Lewatit VP OC 1065 and Purolite A110 (Low et al. 2023) now include water: the Young-type GAB isotherm (`gab_young`, Table 9, Eqs. 17-22), the limiting water heat of adsorption (Table 10) and the measured water isotherm points (288.25, 298.15 and 308.15 K).
- Triamine-grafted alumina (Grossmann et al. 2025) now includes the GAB water isotherm (Table 5), the water heat of adsorption and the CO2 weighted-average dual-site Toth (`toth_wadst`, Tables 4 and 6, referenced to T0 = 308 K).

### Changed — Lewatit VP OC 1065 CO2 isotherm is now the humid WADST fit
- The CO2 isotherm of `lewatit_vp_oc_1065_low2023` is the `toth_wadst` fit of Young et al. (2021) instead of the dry td-Toth fit of Low et al. (2023). At zero water loading it is the dry site of that fit, which differs from the earlier dry isotherm. The dry td-Toth parameters stay in the entry notes.

### Fixed — removing a component left an empty isotherm object
- Removing a component from an adsorbent (including loading a library adsorbent into a deck that lacks one of its components) now also drops an isotherm object that no component uses any more. Before, a library adsorbent with a Water isotherm loaded into a deck without Water left an empty isotherm in the deck, which the deck-to-script round trip could not reproduce. Results do not change: the empty isotherm was never evaluated.
