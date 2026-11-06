### Changed — TSA tutorial heated by a heating medium
- The illustrative TSA deck (`Illustrative_TSA_13X`) and its tutorial now
  heat with a 150 C medium (steam or hot oil) and cool with 25 C water.
  Both pass through a jacket on the wall and an internal heat exchanger in
  the packing, replacing the 30 W heater and the -12 W cooler. The deck's
  results change: CSS at cycle 27, CO2 purity 86.3 %, recovery 96.0 %,
  thermal energy 10.9 MJ/kg, all of it from the medium duty.
- The tutorial explains the choice, with a table showing that a 100 C
  jacket alone is limited by the bed-to-wall heat transfer on this column.

### Added — heating-medium verification
- Analytic heat-up tests against closed forms: the jacket heating the wall
  alone, the jacket heating wall and bed (exact two-node solution), and the
  internal exchanger heating the bed alone.
- `sk.to_script` writes a step's heating medium (its temperature, jacket
  coefficient and internal exchanger). A ramped medium temperature is left
  as a comment, since the builder cannot create the ramp yet.
