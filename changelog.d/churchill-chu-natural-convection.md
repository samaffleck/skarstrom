### Added — natural convection from the column wall (Churchill and Chu 1975)
- Two new wall-to-ambient heat transfer kinds: `churchill_chu_vertical` (the column as a vertical plate of the bed height) and `churchill_chu_horizontal` (a horizontal cylinder of the outer wall diameter). They compute h_o from Churchill and Chu's (1975) natural-convection correlations in each cell, from that cell's wall temperature and the ambient temperature, with the air properties at the film temperature. They have no parameters.
- A run-start warning when a column with the vertical kind is too slender to be a vertical plate (D/L < 35/Gr_L^(1/4) for a wall 10 K above ambient). The plate correlation then under-predicts the heat loss.

