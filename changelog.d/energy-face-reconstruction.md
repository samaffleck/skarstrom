### Changed — face temperature follows the advection scheme
- The energy balance's convective face temperature is now reconstructed with
  the selected advection scheme (WENO or van Leer), with the same stencils and
  the same upwind faces (bed ends, void cells) as the mole fractions. It used
  to be the upwind value whatever the scheme, which kept the temperature first
  order. In the manufactured-solution test the temperature now converges at
  order 2 with WENO and van Leer, and the other fields of a thermal run keep
  their second order on fine grids. Upwind runs are unchanged. Results of
  thermal runs with WENO (the default) or van Leer move; the re-saved example
  decks and the paper comparisons record by how much.
