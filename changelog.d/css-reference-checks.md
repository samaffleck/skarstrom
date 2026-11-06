### Added — CSS stops verified against a long reference
- A slow-tier test (`CSSReference`) checks plain, Aitken and Anderson stops
  and a perturbed warm start against a long reference run, on eight example
  decks. It compares the full end-of-cycle state (loadings, mole fractions,
  temperatures and pressure), the endpoint histories and the KPIs. Every
  distance must be within c times the stop's state tolerance, with one bound
  (c <= 100) for all decks and methods. The measured c is at most 37.
- `docs/cyclic-steady-state.md` tabulates c per deck and method, gives the
  reference drift, and lists the decks left out and why.
