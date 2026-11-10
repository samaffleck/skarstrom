### Added — validation: Halabi (2011), sorption-enhanced steam methane reforming
- New paper suite `Halabi2011Test` (`core/test/papers/halabi_2011/`, ctest label
  `paper`): Rh/Ce0.6Zr0.4O2 catalyst mixed with K2CO3-promoted hydrotalcite
  (M. H. Halabi, PhD thesis, TU Eindhoven 2011, doi:10.6100/IR709035). The
  decks are built in C++ from kind names. The thesis's rate model goes into
  three `lhhw` reactions; its two site balances are expanded exactly into one
  denominator of 14 terms.
- Asserted: with an inert sorbent, the outlet CH₄ conversion matches the
  equilibria the thesis prints (Fig. 7 within 0.46 points, bound 1; Table 6
  within 1.95 points, bound 2.5).
- Recorded only: the twelve sorption-enhanced runs of Table 6. CO₂ breaks
  through after 1.6–3.4 min in the model against 4–15 min in the thesis, and
  the conversion has no plateau before breakthrough. The open question is the
  sorbent capacity below 0.05 bar (#245). See the suite's README and the
  Validation page.
