### Added — inert support solids in the adsorbent library
- Five non-adsorbing support solids for structured contactors, each a
  library material of the new family `support`: cordierite (dense crystal),
  pure aluminium, stainless steel 304L and 316L, and alumina-silica ceramic
  fibre paper. Mix one by mass with a sorbent in a monolith or laminate
  layer's composition (`sim.add_library_adsorbent("cordierite")`). Every
  value is cited (USGS Bulletin 2131 and Open-File Report 88-690, Ho, Powell
  and Liley 1972, ANL-75-55, Headley et al. 2019); the ceramic paper's
  thermal conductivity has no citable value at adsorption temperatures and
  keeps the engine default, as its entry says.
- Library catalogue: a property may carry a `temperature_table` (values
  against temperature, transcribed from the source); the engine uses the
  constant value at 300 K. Support entries need density, heat capacity and
  thermal conductivity instead of a particle size, and their particle
  voidage is 0 by convention (provenance `convention`).
