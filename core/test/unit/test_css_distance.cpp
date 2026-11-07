/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// The remaining-distance diagnostic of the CSS check
// (CyclicSteadyState::EstimateDistance) on synthetic histories of the
// end-of-cycle state change: a geometric decay with a known rate must give that
// rate and the geometric remaining distance, and a history that is too short,
// not decaying, or broken by an acceleration jump must give "not estimated".

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/solver/cyclic_steady_state.h"

namespace adsim {

  namespace {
    // Delta_c = a * rho^c for c = 0..n-1 (oldest first).
    std::vector<double> Geometric(double a, double rho, int n) {
      std::vector<double> d;
      for (int c = 0; c < n; ++c) d.push_back(a * std::pow(rho, c));
      return d;
    }
  }  // namespace

  TEST(CSSDistance, GeometricDecayRecoversRateAndDistance) {
    for (double rho : {0.5, 0.9, 0.97, 0.99}) {
      const auto changes = Geometric(1e-3, rho, 12);
      const double tol = 2e-5;
      const auto d = CyclicSteadyState::EstimateDistance(changes, tol);
      ASSERT_TRUE(d.estimated) << rho << ": " << d.note;
      EXPECT_NEAR(d.rate, rho, 1e-12) << rho;
      // The true remaining distance of a geometric sequence: the sum of all
      // the steps still to come, last * rho / (1 - rho).
      const double expected = changes.back() * rho / (1.0 - rho);
      EXPECT_NEAR(d.distance, expected, 1e-12 * expected) << rho;
      EXPECT_NEAR(d.ratio, expected / tol, 1e-9 * expected / tol) << rho;
      EXPECT_EQ(d.changes_used, CyclicSteadyState::kDistanceMaxChanges);
      EXPECT_TRUE(d.note.empty());
    }
  }

  TEST(CSSDistance, ShortestWindowIsEnough) {
    const auto changes = Geometric(1e-4, 0.8, CyclicSteadyState::kDistanceMinChanges);
    const auto d = CyclicSteadyState::EstimateDistance(changes, 1e-5);
    ASSERT_TRUE(d.estimated) << d.note;
    EXPECT_NEAR(d.rate, 0.8, 1e-12);
    EXPECT_EQ(d.changes_used, CyclicSteadyState::kDistanceMinChanges);
  }

  TEST(CSSDistance, TooShortHistoryIsNotEstimated) {
    for (int n = 0; n < CyclicSteadyState::kDistanceMinChanges; ++n) {
      const auto d = CyclicSteadyState::EstimateDistance(Geometric(1e-4, 0.9, n), 2e-5);
      EXPECT_FALSE(d.estimated) << n;
      EXPECT_EQ(d.note, "history too short") << n;
      EXPECT_TRUE(std::isnan(d.rate) && std::isnan(d.distance) && std::isnan(d.ratio)) << n;
    }
  }

  TEST(CSSDistance, NonMonotoneWindowIsNotEstimated) {
    // An acceleration jump: the change collapses, then the next natural cycle
    // shows the remaining slow mode again. Neither the trailing rise nor a
    // two-change run is a decay to fit.
    const std::vector<double> after_jump{4e-4, 3.6e-4, 3.2e-4, 2e-5, 6e-5};
    auto d = CyclicSteadyState::EstimateDistance(after_jump, 2e-5);
    EXPECT_FALSE(d.estimated);
    EXPECT_EQ(d.note, "not a monotone decay");

    const std::vector<double> rise_then_fall{1e-5, 3e-5, 2e-5};
    d = CyclicSteadyState::EstimateDistance(rise_then_fall, 2e-5);
    EXPECT_FALSE(d.estimated);

    // Noise floor: the changes wander instead of decaying.
    const std::vector<double> noise{1.1e-6, 0.9e-6, 1.2e-6, 1.0e-6, 1.05e-6, 0.95e-6};
    d = CyclicSteadyState::EstimateDistance(noise, 2e-5);
    EXPECT_FALSE(d.estimated);

    // A constant change has no decay: rate 1 would be an infinite distance.
    d = CyclicSteadyState::EstimateDistance({1e-5, 1e-5, 1e-5, 1e-5}, 2e-5);
    EXPECT_FALSE(d.estimated);
  }

  TEST(CSSDistance, UsesOnlyTheTrailingDecay) {
    // A bump earlier in the history does not stop the estimate when the last
    // kDistanceMinChanges or more changes decay; only they are fitted.
    std::vector<double> changes{1e-4, 5e-6, 3e-5};
    for (double v : Geometric(8e-5, 0.9, 4)) changes.push_back(v);
    const auto d = CyclicSteadyState::EstimateDistance(changes, 2e-5);
    ASSERT_TRUE(d.estimated) << d.note;
    EXPECT_EQ(d.changes_used, 4);
    EXPECT_NEAR(d.rate, 0.9, 1e-12);
  }

  TEST(CSSDistance, UnavailableChangesBreakTheWindow) {
    // NaN (no previous cycle) and non-positive values are not decay steps.
    const double nan = std::nan("");
    auto d = CyclicSteadyState::EstimateDistance({1e-4, 8e-5, nan, 5e-5, 4e-5}, 2e-5);
    EXPECT_FALSE(d.estimated);
    EXPECT_EQ(d.changes_used, 2);
    d = CyclicSteadyState::EstimateDistance({1e-4, 8e-5, 6e-5, 0.0}, 2e-5);
    EXPECT_FALSE(d.estimated);
    EXPECT_EQ(d.changes_used, 0);
  }

  TEST(CSSDistance, RateIsClipped) {
    // A decay slower than kDistanceMaxRate is reported at that rate, so the
    // distance stays finite (and large).
    const auto changes = Geometric(1e-5, 0.99995, 5);
    const auto d = CyclicSteadyState::EstimateDistance(changes, 2e-5);
    ASSERT_TRUE(d.estimated);
    EXPECT_EQ(d.rate, CyclicSteadyState::kDistanceMaxRate);
    EXPECT_NEAR(d.distance, changes.back() * 999.0, 1e-9);
    EXPECT_GT(d.ratio, CyclicSteadyState::kDistanceWarningRatio);
  }

  TEST(CSSDistance, NoToleranceGivesNoRatio) {
    // The Molar Profile method has no state tolerance: the distance is still
    // reported in the state-change norm, the ratio is not.
    const auto d = CyclicSteadyState::EstimateDistance(Geometric(1e-4, 0.9, 6), 0.0);
    ASSERT_TRUE(d.estimated);
    EXPECT_TRUE(std::isfinite(d.distance));
    EXPECT_TRUE(std::isnan(d.ratio));
  }

}  // namespace adsim
