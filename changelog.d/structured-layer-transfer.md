### Added — structured layers: channel heat and mass transfer
- Reynolds number of a structured layer on the channel velocity and hydraulic
  diameter, `rho (|u| / OFA) D_h / mu`.
- Channel film coefficients: the fully developed laminar Nu = Sh of Shah and
  London (2.976 square, 3.657 circle, 7.541 plates), or Hawthorn's
  developing-flow form over the layer length with the new per-layer switch
  `layer[ID].entrance_effects` (`StructuredGeometry` v2).
- The macropore rate coefficients use Glueckauf's slab LDF `3 D_eff / l^2`
  (half the wall or sheet, or the washcoat) on a structured layer; Haynes uses
  the slab's film and pore resistances with the channel Sh. Packed beds keep
  `60 D_eff / d_p^2`, unchanged.
- New gas-solid model `film` (bed-wide): Wakao-Funazkri on packed layers,
  the channel Nu on structured layers. The constant model uses the layer's
  specific surface (`6 (1 - eb) / dp` for packed beds, unchanged).
- New dispersion model `taylor_aris`: `D_m + u_ch^2 D_h^2 / (kappa D_m)`,
  kappa 192 (circle; square as an approximation) and 840 (plates).
- `build::SetLayerKind` installs Taylor-Aris and `film` for a structured
  layer and restores the packed defaults when switching back.
- Deck checks reject Ruthven dispersion on a structured layer and Taylor-Aris
  on a packed one.
- Operating limits: a warning when a structured layer's channel Reynolds
  number passes 2000 (`settings.warn_channel_reynolds`, `WarningSettings` v1).
- Verification tests `StructuredTransfer.*` against each source's own result
  and a parallel-passage breakthrough against the LDF moment analysis.
