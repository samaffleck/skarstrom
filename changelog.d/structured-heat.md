### Added — heat in monolith and laminate layers
- Heat-dispersion model `channel` for structured layers: the gas conducts over the open frontal area, λ_f = OFA·k_g, under LTNE; under LTE the single temperature sees OFA·k_g + (1 − OFA)·k_s.
- Axial conduction through the solid of a structured layer under LTNE, λ_s = (1 − OFA)·k_s, with the face conductance of the two half cells in series. It is zero on faces that touch a packed layer or a void and at the bed ends; packed beds are unchanged.
- Verification tests `StructuredHeat.*`: wall-to-solid and gas-to-solid relaxation against the analytic exponentials, steady Fourier conduction through the solid, energy closure, and packed beds free of solid conduction.

### Changed — structured layers take the column-wall flux through the solid
- Under LTNE, the whole bed-to-wall flux of a structured layer goes into the solid; the gas share is 0 instead of the open frontal area. Packed layers keep the ε_b split.
- The deck checks refuse Ruthven and Yagi–Kunii heat dispersion and Dixon and Leva wall coefficients on a structured layer, and the `channel` model on a packed one. The structured-layer warning about these closures is gone; the LTE warning stays.
- `build::SetLayerKind` to a structured kind installs `channel` heat dispersion and a constant h_i (replacing Dixon or Leva; a constant coefficient keeps its value). Switching back to packed restores constant heat dispersion.
