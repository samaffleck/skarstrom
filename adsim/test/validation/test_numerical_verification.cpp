/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Verification, as distinct from validation.
//
// The other files in this directory ask "does the solver solve the right
// equations?". These ask "does it solve the equations right?" -- that the
// discretisation is consistent and converges at its design order, that the
// answer does not depend on knobs it should not depend on (integrator
// tolerance, choice of units of scale), and that symmetries of the continuous
// problem survive discretisation.
//
// These catch a class of defect the physics tests cannot: an error that is
// present on every grid, and therefore invisible to any single fixed-grid
// comparison, but which shows up as the wrong order of accuracy or as a broken
// symmetry.

#include "breakthrough_case.h"

namespace adsim {

  using namespace adsim::validation;


  // ===========================================================================
  // 1. Order of accuracy
  // ===========================================================================
  //
  // With physical dispersion switched off and no adsorption, the exact solution
  // is a travelling step: the breakthrough curve is a shock, with zero variance.
  // What the solver produces instead is a front smeared by the upwind scheme's
  // numerical diffusion, D_num = u*dx/2, giving a response variance proportional
  // to dx.
  //
  // So the variance of the outlet curve is a direct, dimensional measure of the
  // discretisation error with a known exact value of zero, and halving dx must
  // halve it. Measuring that ratio is a measurement of the *observed order of
  // accuracy* -- the standard way to demonstrate a discretisation is consistent,
  // and something no fixed-grid comparison can show.
  //
  // Both moments come from integrals of (1 - F), which avoids differentiating
  // the breakthrough curve numerically:
  //     m1 = int (1 - F) dt,   m2 = 2 int t (1 - F) dt,   var = m2 - m1^2

  TEST(NumericalVerification, UpwindSchemeConvergesAtFirstOrder) {
    std::vector<double> variance;
    std::vector<double> first_moment;

    for (int n_cells : {40, 80, 160, 320}) {
      BreakthroughCase c;
      c.n_cells = n_cells;
      // Inert tracer: with no isotherm and no dispersion the exact response is a
      // pure step, so every bit of the measured variance is numerical diffusion.
      // (A finite LDF coefficient would add a grid-independent kinetic variance
      // that the refinement sequence would converge to instead of to zero.)
      c.adsorbing = false;
      c.run_time = 250.0;
      c.max_step = 0.1;

      UserData data;
      Build(data, c);

      CoupledSolver solver{data};
      std::string error;
      ASSERT_TRUE(data.PerformChecks(error)) << error;
      ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0)
        << "Solver did not complete at n_cells = " << n_cells;

      auto [times, y_out] = OutletSeries(solver.GetLog(), c.feed_from_top);
      ASSERT_FALSE(times.empty());
      ASSERT_NEAR(y_out.back(), c.y_feed, 1e-3 * c.y_feed)
        << "Breakthrough incomplete at n_cells = " << n_cells;

      std::vector<double> unbroken(times.size()), t_unbroken(times.size());
      for (size_t n = 0; n < times.size(); ++n) {
        unbroken[n] = 1.0 - y_out[n] / c.y_feed;
        t_unbroken[n] = times[n] * unbroken[n];
      }
      const double m1 = TimeIntegral(times, unbroken);
      const double m2 = 2.0 * TimeIntegral(times, t_unbroken);
      first_moment.push_back(m1);
      variance.push_back(m2 - m1 * m1);
    }

    ASSERT_EQ(variance.size(), 4u);

    // The first moment is grid-invariant (it is fixed by the column capacity),
    // so it must not drift as the grid changes.
    for (size_t i = 1; i < first_moment.size(); ++i) {
      EXPECT_NEAR(first_moment[i], first_moment[0], 0.01 * first_moment[0])
        << "Stoichiometric time drifted with grid refinement";
    }

    // ... while the variance, which is pure numerical diffusion, must halve.
    for (size_t i = 1; i < variance.size(); ++i) {
      ASSERT_GT(variance[i], 0.0);
      const double ratio = variance[i - 1] / variance[i];
      const double observed_order = std::log2(ratio);
      EXPECT_NEAR(observed_order, 1.0, 0.2)
        << "Observed spatial order of accuracy is " << observed_order
        << ", expected 1 for upwind (variance " << variance[i - 1] << " -> "
        << variance[i] << ")";
    }
  }

  // ===========================================================================
  // 2. Integrator tolerance independence
  // ===========================================================================
  //
  // The spatial discretisation fixes the answer; the time integrator should only
  // control how accurately that answer is reached. Loosening reltol by five
  // orders of magnitude must therefore change the result by far less than the
  // discretisation error. If it does not, the defaults are under-resolved and
  // every result in the suite is quietly tolerance-dependent.

  TEST(NumericalVerification, SolutionIsIndependentOfIntegratorTolerance) {
    std::vector<std::vector<double>> outlet;
    std::vector<std::vector<double>> times;

    for (double reltol : {1e-4, 1e-9}) {
      BreakthroughCase c;
      c.adsorbing = true;
      c.reltol = reltol;
      c.run_time = 600.0;

      UserData data;
      Build(data, c);

      CoupledSolver solver{data};
      std::string error;
      ASSERT_TRUE(data.PerformChecks(error)) << error;
      ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

      auto [t, y] = OutletSeries(solver.GetLog(), c.feed_from_top);
      times.push_back(t);
      outlet.push_back(y);
    }

    // Compare on a common time grid, since the two runs take different steps.
    double max_difference = 0.0;
    for (double t = 10.0; t <= 590.0; t += 10.0) {
      const double loose = interp_linear(times[0], outlet[0], t);
      const double tight = interp_linear(times[1], outlet[1], t);
      max_difference = std::max(max_difference, std::abs(loose - tight));
    }

    EXPECT_LT(max_difference, 1e-4 * 0.05)
      << "Outlet composition moved by " << max_difference
      << " when reltol was tightened from 1e-4 to 1e-9; the default integrator "
         "tolerance is not resolving the spatial solution";
  }

  // ===========================================================================
  // 3. Scale invariance
  // ===========================================================================
  //
  // The governing equations depend on the column cross-section only through the
  // superficial velocity, u = n_dot R T / (A P). Doubling the diameter while
  // quadrupling the molar feed rate leaves u unchanged, so every intensive field
  // -- pressure, temperature, composition, loading -- must be bit-for-bit the
  // same function of axial position and time. Any place where an area, a volume
  // or a diameter has leaked into the physics with the wrong power shows up
  // immediately.

  TEST(NumericalVerification, ResultsAreInvariantUnderColumnScaling) {
    std::vector<std::vector<double>> final_composition;
    std::vector<std::vector<double>> final_pressure;

    for (int scale : {1, 2}) {
      BreakthroughCase c;
      c.adsorbing = true;
      c.diameter = 0.1 * scale;
      c.molar_flow = 0.005 * scale * scale;  // same superficial velocity
      c.run_time = 75.0;  // front roughly mid-column

      UserData data;
      Build(data, c);

      CoupledSolver solver{data};
      std::string error;
      ASSERT_TRUE(data.PerformChecks(error)) << error;
      ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

      const auto& log = solver.GetLog();
      final_composition.push_back(
        log.GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION, "Oxygen")
           .values.back());
      final_pressure.push_back(log.GetVariable(0, 0, VariableType::PRESSURE).values.back());
    }

    ASSERT_EQ(final_composition[0].size(), final_composition[1].size());

    // The front must be somewhere in the middle of the bed, or the comparison is
    // between two uniform fields and proves nothing.
    const auto& y = final_composition[0];
    ASSERT_GT(*std::max_element(y.begin(), y.end()), 0.9 * 0.05);
    ASSERT_LT(*std::min_element(y.begin(), y.end()), 0.1 * 0.05);

    for (size_t i = 0; i < final_composition[0].size(); ++i) {
      EXPECT_NEAR(final_composition[1][i], final_composition[0][i], 1e-6)
        << "Composition profile changed under column scaling at cell " << i;
      EXPECT_NEAR(final_pressure[1][i], final_pressure[0][i], 1e-3)
        << "Pressure profile changed under column scaling at cell " << i;
    }
  }

  // ===========================================================================
  // 4. Flow-direction symmetry
  // ===========================================================================
  //
  // A uniform bed has no preferred direction. Feeding the same stream in at the
  // top instead of the bottom must produce the mirror image of the same profile.
  // This is the sharpest available test of the upwinding: the scheme switches
  // branch when the velocity changes sign, and an asymmetry there (a wrong
  // index, a face/cell mix-up, a boundary condition that only handles one
  // orientation) is invisible in every co-current test.

  TEST(NumericalVerification, TopFedBedMirrorsBottomFedBed) {
    std::vector<std::vector<double>> profile;

    for (bool from_top : {false, true}) {
      BreakthroughCase c;
      c.adsorbing = true;
      c.feed_from_top = from_top;
      c.run_time = 75.0;  // front roughly mid-column

      UserData data;
      Build(data, c);

      CoupledSolver solver{data};
      std::string error;
      ASSERT_TRUE(data.PerformChecks(error)) << error;
      ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0)
        << "Solver did not complete (feed from top = " << from_top << ")";

      profile.push_back(
        solver.GetLog()
          .GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION, "Oxygen")
          .values.back());
    }

    ASSERT_EQ(profile[0].size(), profile[1].size());
    const size_t n = profile[0].size();

    // Front must be partway through, otherwise both profiles are flat.
    ASSERT_GT(*std::max_element(profile[0].begin(), profile[0].end()), 0.9 * 0.05);
    ASSERT_LT(*std::min_element(profile[0].begin(), profile[0].end()), 0.1 * 0.05);

    // Compare interior cells against their mirror image.
    for (size_t i = 1; i + 1 < n; ++i) {
      const size_t mirrored = n - 1 - i;
      EXPECT_NEAR(profile[1][mirrored], profile[0][i], 1e-6)
        << "Top-fed profile is not the mirror of the bottom-fed profile at cell " << i
        << " (bottom-fed " << profile[0][i] << " vs top-fed mirrored "
        << profile[1][mirrored] << ")";
    }
  }

} // namespace adsim
