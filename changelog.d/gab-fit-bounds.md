### Fixed — GAB fit set-up and its test
- The default fit set-up for the GAB kinds keeps K between 0.3 and 1.2 and
  the monolayer capacity below ten times the largest loading, the physical
  range for water, so the search stays out of the near-linear corner
  (K -> 0 with a very large q_m) where the parameters are not identifiable.
- The `gab_young` fit test used a non-physical target (K = 0.14 at 298 K),
  where the fit is ill-conditioned and ended in a different corner on macOS
  arm64; it now uses the consistent Young et al. (2021) Lewatit set, and a
  new test checks that the GAB fits recover the same parameters for every
  global-search budget and a last-digit perturbation of the data.
