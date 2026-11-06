### Added — CSS start independence verified
- A slow-tier test (`CSSStartIndependence`) runs eight example decks to CSS
  from the shipped start, a clean bed, a feed-saturated bed, a hot and a
  cold bed (+-40 K, through the warm-start API) and, where one exists, a
  neighbouring deck's CSS. Every start must reach the shipped start's CSS
  within c <= 100 times the state tolerance, in the full state, the
  histories and the KPIs. The measured c is at most 74 (Perez 2019 E), and
  no deck has a second steady state.
- `docs/cyclic-steady-state.md` tabulates the result, and explains why the
  hot and cold starts do not edit the initial temperature: with CoolProp
  gas heat capacities that also changes the model for the whole run.
