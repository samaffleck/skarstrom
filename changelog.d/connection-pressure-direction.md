### Added — connection pressure warning
- A new end-of-run warning, `connection_pressure` (`settings.warn_connection_pressure`, on by
  default), reports gas that enters a bed through a connection while that bed's end pressure is
  more than 1 kPa above the partner's logged pressure (more than 0.1 % of the cycle feed in one
  step, last cycle; pump streams are not checked). It is a diagnostic: a connected endpoint still
  replays its partner's previous run and does not enforce the order of the pressures. (A cap on
  connected pressure profiles and a pressure-difference scaling of replayed flows were tried and
  withdrawn before release: they made the light-product decks settle more slowly and noisily.)
