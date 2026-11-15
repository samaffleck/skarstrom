### Added — natural convection from a slender vertical column (Yang 1985)
- A new wall-to-ambient heat transfer kind, `yang_vertical_cylinder`: natural convection from a vertical cylinder of the bed height and the outer wall diameter, valid at any height-to-diameter ratio. It is Churchill and Chu's plate correlation with the conduction term replaced by 0.36 L/d, so a slender laboratory column loses more heat than the vertical-plate kind gives (about 44 % more for a 1 m column 30 mm across). It has no parameters.
- The slender-column warning of `churchill_chu_vertical` now names this kind.
