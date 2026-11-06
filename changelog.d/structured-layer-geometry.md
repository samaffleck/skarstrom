### Added — structured layers: monolith and laminate geometry
- A layer can now be a monolith (square or circular channels) or a laminate
  (sheets with a gap) as well as a packed bed: `MixedAdsorbentLayer::kind_`
  with `StructuredGeometry` (cell density, wall thickness, sheet thickness,
  gap, washcoat thickness). Builder: `build::SetLayerKind`,
  `build::SetChannelShape`, `build::LayerKinds`, `build::ChannelShapeKinds`.
- The flow in a structured layer is laminar Hagen-Poiseuille,
  `K = OFA D_h^2 / (2 fRe)` with fRe 14.23 (square), 16 (circle), 24 (plates),
  and no inertial term; the layer's momentum model is not used. The open
  frontal area takes the bed void fraction's place (gas holdup, bulk density,
  dispersion area, heat split).
- Derived geometry for front ends: `LayerGeometrySummary` (open fraction,
  hydraulic diameter, specific surface, diffusion length, channel count,
  washcoat fraction).
- Parameter registry: a structured layer lists its kind's fields
  (`channel_shape`, `cell_density`, `wall_thickness` for a monolith;
  `sheet_thickness`, `gap` for a laminate; `washcoat_thickness`) and keeps
  `particle_diameter` while the packed-bed transfer closures still read it;
  `bed_void_fraction` and `momentum_model` are not listed for it.
- Deck checks refuse a geometry without an open channel. Run-start warnings:
  LTE with a structured layer (LTNE recommended), and the packed-bed closures
  a structured layer still uses.
- `MixedAdsorbentLayer` is serialised at version 4 (`Kind`,
  `Structured_Geometry`), `StructuredGeometry` at version 1. Older decks load
  as packed beds with unchanged results.

### Removed
- The unused `StructuredLayer` class.
