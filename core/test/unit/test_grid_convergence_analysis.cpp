/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// The grid convergence analysis (adsim/solver/grid_convergence.h): observed
// order, Richardson extrapolation, GCI and cells for a target error on exact
// power-law sequences (any constant ratio), the status of sequences that are
// not in the asymptotic range, and the two published grid studies: the table
// in core/test/validation/test_grid_convergence.cpp and the example-deck table
// in docs/numerics.md ("Convergence of the results").

#include "adsim/solver/grid_convergence.h"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <vector>

namespace adsim {

  namespace {

    // f(h) = f* + a h^p on grids h_k = r^-k, coarse to fine.
    std::vector<double> PowerLaw(double f_star, double a, double p, double r, int levels) {
      std::vector<double> f;
      for (int k = 0; k < levels; ++k) f.push_back(f_star + a * std::pow(r, -p * k));
      return f;
    }

  } // namespace

  // ===========================================================================
  // Exact power laws: every estimate is recovered to rounding
  // ===========================================================================

  TEST(GridConvergenceAnalysis, RecoversOrderAndLimitOfAPowerLawForAnyRatio) {
    for (double r : {2.0, 1.5, 3.0, std::sqrt(2.0)}) {
      for (double p : {1.0, 1.5, 2.4}) {
        for (int levels : {3, 5}) {
          const double f_star = 0.7, a = -0.05;
          const auto f = PowerLaw(f_star, a, p, r, levels);
          const auto res = AnalyseGridConvergence(f, r);
          SCOPED_TRACE(testing::Message() << "r = " << r << ", p = " << p << ", " << levels << " levels");
          ASSERT_EQ(res.status, GridConvergenceStatus::Monotone);
          EXPECT_TRUE(res.InAsymptoticRange());
          EXPECT_TRUE(res.monotone);
          EXPECT_TRUE(res.shrinking);
          EXPECT_FALSE(res.order_assumed);
          EXPECT_NEAR(res.order, p, 1e-9);
          EXPECT_NEAR(res.extrapolated, f_star, 1e-12);
          ASSERT_EQ(res.differences.size(), size_t(levels - 1));
          ASSERT_EQ(res.orders.size(), size_t(levels - 2));
          for (double q : res.orders) EXPECT_NEAR(q, p, 1e-9);
          // Exact Richardson: |d_fine| / (r^p - 1) is the finest grid's error.
          const double error = std::abs(f.back() - f_star);
          EXPECT_NEAR(res.gci, 1.25 * error / std::abs(f.back()), 1e-12);
          EXPECT_NEAR(res.error_band, 1.25 * error, 1e-12);
          EXPECT_TRUE(std::isnan(res.cells_for_target));  // no cell count given
        }
      }
    }
  }

  TEST(GridConvergenceAnalysis, CellsForTargetErrorFollowsThePowerLaw) {
    // Cells N = N0 r^k; f(N) = f* + a (N0 / N)^p. The cell count at which
    // |f(N) - f*| = target is N0 (|a| / target)^(1/p).
    const double f_star = 0.9, a = 0.02, p = 1.2, n0 = 10.0;
    for (double r : {2.0, 1.5}) {
      const auto f = PowerLaw(f_star, a, p, r, 4);
      const double n_fine = n0 * std::pow(r, 3);
      GridConvergenceOptions opt;
      opt.target_error = 1e-3;
      const auto res = AnalyseGridConvergence(f, r, opt, n_fine);
      ASSERT_EQ(res.status, GridConvergenceStatus::Monotone);
      const double expected = n0 * std::pow(std::abs(a) / opt.target_error, 1.0 / p);
      EXPECT_NEAR(res.cells_for_target, expected, 1e-9 * expected) << "r = " << r;
      // A loose target is met below the finest grid; the estimate is not clamped.
      opt.target_error = 1e-2;
      EXPECT_NEAR(AnalyseGridConvergence(f, r, opt, n_fine).cells_for_target,
                  n0 * std::pow(std::abs(a) / 1e-2, 1.0 / p), 1e-9);
    }
  }

  TEST(GridConvergenceAnalysis, CellCountsGiveTheRatio) {
    const auto f2 = PowerLaw(0.5, 0.1, 1.0, 2.0, 5);
    const auto res = AnalyseGridConvergence(std::vector<int>{10, 20, 40, 80, 160}, f2);
    EXPECT_EQ(res.refinement_ratio, 2.0);  // exact when every pair doubles
    ASSERT_EQ(res.status, GridConvergenceStatus::Monotone);
    EXPECT_NEAR(res.order, 1.0, 1e-9);
    // 1e-3 error at N = 10 * (0.1 / 1e-3) = 1000 cells.
    EXPECT_NEAR(res.cells_for_target, 1000.0, 1e-6);

    const auto f15 = PowerLaw(0.5, 0.1, 2.0, 1.5, 3);
    const auto r15 = AnalyseGridConvergence(std::vector<int>{12, 18, 27}, f15);
    EXPECT_EQ(r15.refinement_ratio, 1.5);
    EXPECT_NEAR(r15.order, 2.0, 1e-9);

    // Cell counts rounded per layer: the ratio is the geometric mean.
    const auto rounded = AnalyseGridConvergence(std::vector<int>{10, 15, 22, 33}, PowerLaw(0.5, 0.1, 1.0, 1.5, 4));
    EXPECT_EQ(rounded.status, GridConvergenceStatus::Monotone);
    EXPECT_NEAR(rounded.refinement_ratio, std::cbrt(3.3), 1e-12);

    // Not a constant ratio, or mismatched inputs.
    EXPECT_EQ(AnalyseGridConvergence(std::vector<int>{10, 20, 30}, f15).status, GridConvergenceStatus::Invalid);
    EXPECT_EQ(AnalyseGridConvergence(std::vector<int>{10, 20}, f15).status, GridConvergenceStatus::Invalid);
    EXPECT_EQ(AnalyseGridConvergence(std::vector<int>{0, 20, 40}, f15).status, GridConvergenceStatus::Invalid);
  }

  // ===========================================================================
  // Not in the asymptotic range
  // ===========================================================================

  TEST(GridConvergenceAnalysis, OscillatoryConvergenceHasNoOrder) {
    const auto res = AnalyseGridConvergence({0.60, 0.62, 0.615}, 2.0);
    EXPECT_EQ(res.status, GridConvergenceStatus::Oscillatory);
    EXPECT_FALSE(res.InAsymptoticRange());
    EXPECT_FALSE(res.monotone);
    EXPECT_TRUE(std::isnan(res.order));
    EXPECT_TRUE(std::isnan(res.orders[0]));
    EXPECT_TRUE(std::isnan(res.extrapolated));
    EXPECT_TRUE(std::isnan(res.cells_for_target));
    // GCI with the assumed order (1): 1.25 |d| / |f| / (2 - 1).
    EXPECT_TRUE(res.order_assumed);
    EXPECT_NEAR(res.gci, 1.25 * 0.005 / 0.615, 1e-12);

    // A sign change on the coarse grids makes the whole study oscillatory, even
    // when the three finest grids alone would give an order.
    const auto coarse = AnalyseGridConvergence({0.50, 0.70, 0.60, 0.64, 0.66}, 2.0, {}, 160);
    EXPECT_EQ(coarse.status, GridConvergenceStatus::Oscillatory);
    ASSERT_EQ(coarse.orders.size(), 3u);
    EXPECT_TRUE(std::isnan(coarse.orders[0]));
    EXPECT_NEAR(coarse.orders[2], 1.0, 1e-12);
    EXPECT_TRUE(std::isnan(coarse.order));
  }

  TEST(GridConvergenceAnalysis, GrowingDifferencesAreDivergent) {
    // Growing, equal (exact in binary), and equal on the finest pair only.
    for (const auto& f : {std::vector<double>{1.0, 1.1, 1.3}, std::vector<double>{1.0, 1.25, 1.5},
                          std::vector<double>{1.0, 1.5, 1.75, 2.0}}) {
      const auto res = AnalyseGridConvergence(f, 2.0, {}, 80);
      EXPECT_EQ(res.status, GridConvergenceStatus::Divergent);
      EXPECT_TRUE(res.monotone);
      EXPECT_FALSE(res.shrinking);
      EXPECT_TRUE(std::isnan(res.order));
      EXPECT_TRUE(std::isnan(res.gci));
      EXPECT_TRUE(std::isnan(res.extrapolated));
      EXPECT_TRUE(std::isnan(res.cells_for_target));
    }
  }

  TEST(GridConvergenceAnalysis, EdgeCases) {
    EXPECT_EQ(AnalyseGridConvergence({}, 2.0).status, GridConvergenceStatus::TooFewLevels);
    const auto two = AnalyseGridConvergence({0.5, 0.6}, 2.0);
    EXPECT_EQ(two.status, GridConvergenceStatus::TooFewLevels);
    EXPECT_EQ(two.differences.size(), 1u);
    EXPECT_TRUE(std::isnan(two.order));

    const double nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(AnalyseGridConvergence({0.5, nan, 0.6}, 2.0).status, GridConvergenceStatus::Invalid);
    EXPECT_EQ(AnalyseGridConvergence({0.5, 0.55, 0.57}, 1.0).status, GridConvergenceStatus::Invalid);
    EXPECT_EQ(AnalyseGridConvergence({0.5, 0.55, 0.57}, nan).status, GridConvergenceStatus::Invalid);

    // The two finest values equal: grid independent at the finest grid.
    for (const auto& f : {std::vector<double>{0.5, 0.6, 0.6}, std::vector<double>{0.6, 0.6, 0.6}}) {
      const auto res = AnalyseGridConvergence(f, 2.0);
      EXPECT_EQ(res.status, GridConvergenceStatus::GridIndependent);
      EXPECT_EQ(res.extrapolated, 0.6);
      EXPECT_EQ(res.gci, 0.0);
      EXPECT_TRUE(std::isnan(res.order));
    }

    EXPECT_EQ(ToString(GridConvergenceStatus::Monotone), "monotone");
    EXPECT_EQ(ToString(GridConvergenceStatus::Oscillatory), "oscillatory");
    EXPECT_EQ(ToString(GridConvergenceStatus::TooFewLevels), "too few levels");

    // Building blocks: NaN where undefined.
    EXPECT_TRUE(std::isnan(ObservedOrder(1.0, 1.1, 1.05, 2.0)));
    EXPECT_TRUE(std::isnan(ObservedOrder(1.0, 1.1, 1.1, 2.0)));
    EXPECT_TRUE(std::isnan(GridConvergenceIndex(0.1, 0.0, 2.0, 1.0)));
    EXPECT_TRUE(std::isnan(RichardsonExtrapolation(0.1, 0.2, 2.0, 0.0)));
    EXPECT_TRUE(std::isnan(CellsForTargetError(0.5, 0.51, 80, 1.0, 0.0)));
  }

  // ===========================================================================
  // The published grid studies
  // ===========================================================================

  // The measured table of core/test/validation/test_grid_convergence.cpp
  // (O2 purity, recovery and specific energy on 10, 20, 40, 80 cells), printed
  // to five significant figures, with p, GCI and Richardson value computed
  // from the unrounded KPIs. Recomputed from the rounded values the order
  // moves by up to 0.12 (WENO purity, whose finest difference 9e-5 has only
  // two printed digits) and the GCI by up to 18 %; the extrapolated value
  // stays within 1.6e-5 relative of the printed one (itself rounded to five
  // figures, i.e. to 1.3e-5 relative).
  struct Row {
    const char* name;
    std::vector<double> f;
    double order, gci, extrapolated;
  };

  TEST(GridConvergenceAnalysis, ReproducesTheCycleKpiGridStudy) {
    const std::vector<Row> rows = {
      {"upwind purity", {0.38801, 0.39304, 0.39503, 0.39572}, 1.54, 1.2e-3, 0.39608},
      {"WENO purity", {0.39339, 0.39555, 0.39602, 0.39611}, 2.50, 5.7e-5, 0.39613},
      {"van Leer purity", {0.39420, 0.39577, 0.39606, 0.39611}, 2.47, 3.6e-5, 0.39612},
      {"upwind recovery", {0.67868, 0.68348, 0.68533, 0.68597}, 1.55, 6.0e-4, 0.68630},
      {"WENO recovery", {0.68377, 0.68581, 0.68625, 0.68633}, 2.48, 3.1e-5, 0.68634},
      {"van Leer recovery", {0.68453, 0.68601, 0.68628, 0.68633}, 2.40, 2.2e-5, 0.68634},
      {"upwind energy", {854701, 848717, 846428, 845650}, 1.56, 5.9e-4, 845250},
      {"WENO energy", {848354, 845841, 845308, 845213}, 2.48, 3.1e-5, 845192},
      {"van Leer energy", {847419, 845599, 845271, 845209}, 2.41, 2.1e-5, 845194},
    };
    for (const auto& row : rows) {
      SCOPED_TRACE(row.name);
      const auto res = AnalyseGridConvergence(std::vector<int>{10, 20, 40, 80}, row.f);
      ASSERT_EQ(res.status, GridConvergenceStatus::Monotone);
      EXPECT_NEAR(res.order, row.order, 0.15);
      EXPECT_NEAR(res.gci, row.gci, 0.25 * row.gci);
      EXPECT_NEAR(res.extrapolated, row.extrapolated, 2.5e-5 * row.extrapolated);
      EXPECT_EQ(res.finest, row.f.back());
    }
  }

  // docs/numerics.md, "Convergence of the results": the example decks on 10 to
  // 160 cells, printed to four digits. The rounding is too coarse to recompute
  // the orders (a finest difference of 0.0002 has one digit), so the order,
  // GCI and cells columns are checked through the formulas with the printed
  // p, f80, f160 and f*: GCI and cells within 15 % (the printed values carry
  // the unrounded KPIs). From the five values themselves, each recovery
  // sequence is monotone and its Richardson value within 3e-4 of the printed
  // f*. (The rows whose two finest printed values are equal are left out.)
  struct DeckRow {
    const char* name;
    std::vector<double> f;
    double order, gci, extrapolated, cells;
  };

  TEST(GridConvergenceAnalysis, ReproducesTheExampleDeckGridStudy) {
    const std::vector<DeckRow> rows = {
      {"Perez2019_A upwind recovery", {0.5772, 0.5874, 0.5921, 0.5945, 0.5957}, 0.96, 0.26e-2, 0.5969, 200},
      {"Perez2019_A WENO recovery", {0.5872, 0.5934, 0.5956, 0.5965, 0.5968}, 1.49, 0.036e-2, 0.5970, 50},
      {"Perez2019_A van Leer recovery", {0.5897, 0.5951, 0.5963, 0.5967, 0.5969}, 1.40, 0.023e-2, 0.5970, 35},
      {"Ward2024 upwind recovery", {0.8936, 0.9054, 0.9120, 0.9149, 0.9165}, 0.86, 0.27e-2, 0.9185, 360},
      {"Ward2024 WENO recovery", {0.9040, 0.9129, 0.9166, 0.9180, 0.9187}, 1.14, 0.074e-2, 0.9192, 95},
      {"Ward2024 van Leer recovery", {0.9097, 0.9157, 0.9176, 0.9185, 0.9189}, 0.96, 0.062e-2, 0.9193, 70},
      {"Illustrative_TSA_13X upwind recovery", {0.6602, 0.6634, 0.6649, 0.6655, 0.6658}, 1.12, 0.05e-2, 0.6661, 50},
    };
    const std::vector<int> cells = {10, 20, 40, 80, 160};
    for (const auto& row : rows) {
      SCOPED_TRACE(row.name);
      const double f80 = row.f[3], f160 = row.f[4];
      EXPECT_NEAR(GridConvergenceIndex(f80, f160, 2.0, row.order), row.gci, 0.15 * row.gci);
      EXPECT_NEAR(CellsForTargetError(f160, row.extrapolated, 160, row.order, 1e-3), row.cells, 0.15 * row.cells);
      EXPECT_NEAR(RichardsonExtrapolation(f80, f160, 2.0, row.order), row.extrapolated, 1e-4);

      const auto res = AnalyseGridConvergence(cells, row.f);
      EXPECT_EQ(res.status, GridConvergenceStatus::Monotone);
      EXPECT_NEAR(res.extrapolated, row.extrapolated, 3e-4);
    }
    // The upwind purities of Ward2024 (falling) and the TSA deck are monotone too.
    EXPECT_EQ(AnalyseGridConvergence(cells, {0.9643, 0.9628, 0.9618, 0.9613, 0.9610}).status,
              GridConvergenceStatus::Monotone);
    EXPECT_EQ(AnalyseGridConvergence(cells, {0.8333, 0.8343, 0.8348, 0.8350, 0.8351}).status,
              GridConvergenceStatus::Monotone);
  }

} // namespace adsim
