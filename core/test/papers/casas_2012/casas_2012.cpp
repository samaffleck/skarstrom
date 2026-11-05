// STL includes
#define _USE_MATH_DEFINES
#include <cmath>
#include <filesystem>

// AdSim includes
#include "../../../src/pch.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/solver/kpis.h"
#include "adsim/math/constants.h"
#include "adsim/equations/isotherm.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"

// GTest includes
#include <gtest/gtest.h>

// Paper fixture (Casas et al. 2012) — see README.md in this folder.
#include "casas_2012.h"

namespace adsim {


  // Reference experiment (Casas 2012, Fig. 4): 50/50 CO2/H2, 15 bar, 25 °C.
  // This is the single experiment the paper's four dynamic parameters (h_L, h_W,
  // k_CO2, k_H2) were fitted to; every other case is a prediction with the same
  // parameter set.
  TEST_F(Casas2012Test, Fig4_Reference_50_50_15bar_25C) {
    data.info.title = "Casas 2012 - Fig. 4 (reference, 50/50, 15 bar, 25 C)";
    data.info.description =
        "Binary CO2/H2 breakthrough on AP3-60 activated carbon, reproducing the "
        "reference experiment of Casas et al. (2012), Fig. 4: an equimolar CO2/H2 "
        "feed at 15 bar and 25 C. The bed is initially pressurised with inert "
        "helium and then fed the CO2/H2 mixture at a constant superficial flow of "
        "10 cm3/s. A competitive Sips isotherm (parameters from the companion "
        "equilibrium study, Schell et al. 2012) describes the equilibrium; mass "
        "transfer uses lumped LDF coefficients fitted in the paper.";
    data.info.author = "Sam Affleck";

    RunBreakthroughCase(15e5, 298.15, 0.50, 1800.0, "Casas2012_Fig4_50_50_15bar_25C",
                        "Fig4_Reference_50_50_15bar_25C",
                        0.18, 0.024);  // measured rel. error: H2 0.16443, CO2 0.007326 (before the WP19-1 competitive-Sips trace fix: H2 0.167584, CO2 0.0092)
  }

  // Prediction (Casas 2012, Fig. 6 family): CO2-lean feed, 25/75 CO2/H2, 15 bar, 25 C.
  TEST_F(Casas2012Test, Fig6_25_75_15bar_25C) {
    data.info.title = "Casas 2012 - Fig. 6 (25/75, 15 bar, 25 C)";
    data.info.description =
        "Binary CO2/H2 breakthrough on AP3-60 activated carbon (Casas et al. 2012). "
        "CO2-lean 25/75 CO2/H2 feed at 15 bar and 25 C, predicted with the same "
        "isotherm and transport parameters fitted to the 50/50 reference case.";
    data.info.author = "Sam Affleck";

    RunBreakthroughCase(15e5, 298.15, 0.25, 1800.0, "Casas2012_Fig6_25_75_15bar_25C",
                        "Fig6_25_75_15bar_25C",
                        0.22, 0.016);  // measured rel. error: H2 0.194276, CO2 0.004447 (before WP19-1: H2 0.197555, CO2 0.006444)
  }

  // Prediction (Casas 2012, Fig. 6 family): CO2-rich feed, 75/25 CO2/H2, 15 bar, 25 C.
  TEST_F(Casas2012Test, Fig6_75_25_15bar_25C) {
    data.info.title = "Casas 2012 - Fig. 6 (75/25, 15 bar, 25 C)";
    data.info.description =
        "Binary CO2/H2 breakthrough on AP3-60 activated carbon (Casas et al. 2012). "
        "CO2-rich 75/25 CO2/H2 feed at 15 bar and 25 C, predicted with the same "
        "isotherm and transport parameters fitted to the 50/50 reference case.";
    data.info.author = "Sam Affleck";

    RunBreakthroughCase(15e5, 298.15, 0.75, 1800.0, "Casas2012_Fig6_75_25_15bar_25C",
                        "Fig6_75_25_15bar_25C",
                        0.047, 0.15);  // measured rel. error: H2 0.015323, CO2 0.099855 (before WP19-1: H2 0.011687, CO2 0.098061)
  }

} // namespace adsim
