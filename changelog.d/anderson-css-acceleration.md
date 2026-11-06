### Added — Anderson acceleration of the approach to cyclic steady state
- New CSS acceleration method *Anderson* (`AccelerationMethod::ANDERSON`).
  It moves the slow fields (loadings, solid and wall temperatures, mole
  fractions) along a combination of up to four earlier cycle-to-cycle steps,
  so it can remove several slowly decaying modes in one jump, where Aitken
  removes one. It uses the same projection and the same re-seed check as
  Aitken. Measured CSS stops (plain / Aitken / Anderson): Perez 2019 E
  205 / 65 / 45, Nguyen 2023 X over 800 / 548 / 190, Perez 2019 A 48 / 32 / 35,
  at KPIs as close to a long reference run as the plain stop's. Decks that do
  not select it are unchanged.
