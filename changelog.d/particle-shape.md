### Added — cylindrical particles (extrudates, fibres) in packed layers
- A packed layer has a particle shape: `sphere` (the default, and the shape of every existing deck) or
  `cylinder`, an extrudate or fibre of diameter `particle_diameter` = 2 r_f, taken as infinitely long.
  The linear-driving-force diffusion geometry follows the shape: an LDF constant of 8 on the radius
  (k = 8 D / r_f^2) and V/S = r_f / 2 for the film term, from Glueckauf's parabolic profile. The film,
  macropore, Knudsen and capacity-factor forms of the generic coefficient model use it; the micropore
  term stays that of a spherical crystal. Pressure drop, dispersion, the film correlations and the
  specific surface area keep d_p as an equivalent-sphere diameter.
- Registry leaf `layer[L].particle_shape`, deck edit `SetParticleShape` (Python `Layer.particle_shape`,
  `sim.edit` op `set_particle_shape`, `build_kinds()["particle_shape"]`), deck-view layer field
  `particle_shape`, and the `cylinder` LaTeX fragments of the resistance descriptors (the display
  equation takes the layer's shape). `MixedAdsorbentLayer` is at class version 7; older decks load
  as spheres and give bitwise the same results.
