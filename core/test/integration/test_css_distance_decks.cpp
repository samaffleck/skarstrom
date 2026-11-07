/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// The remaining-distance diagnostic of the CSS check on bundled decks
// (CyclicSteadyState::EstimateDistance; synthetic cases are in
// unit/test_css_distance.cpp).
//
//  * Illustrative_TSA_13X with Aitken acceleration (fast tier): the cycles
//    that follow a committed jump report "not estimated", the solver's stop
//    line carries the estimate, and css_history() agrees with the solver.
//  * Perez2019_E as shipped (slow tier, ~90 s): the slowest deck of the
//    long-reference verification (test_css_reference.cpp), whose plain stop at
//    cycle 205 was measured 37 tolerances from the reference CSS in the state
//    norm. The estimate at the stop must be of that size.

#include <cmath>
#include <cstdio>
#include <set>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/simulation.h"

namespace adsim {

  namespace {
    // Runs the deck at MEDIUM log level, recording the 0-based cycles that
    // started from a committed acceleration jump (from the order of the
    // solver's messages: one CSS check line per cycle, then the acceleration
    // outcome) and the "CSS reached" line.
    struct Observed {
      RunResult run;
      std::set<int> jump_cycles;
      std::string stop_line;
      std::vector<std::string> warnings;
    };

    Observed RunObserved(Simulation& sim) {
      Observed o;
      int checked = 0;  // CSS check lines seen = index of the next cycle
      RunOptions options;
      options.log_level = LogLevel::MEDIUM;
      options.on_message = [&](LogLevel, const std::string& m) {
        if (m.find("CSS check:") != std::string::npos) ++checked;
        if (m.rfind("CSS reached at cycle", 0) == 0) { ++checked; o.stop_line = m; }
        if (m.find("CSS acceleration applied") != std::string::npos) o.jump_cycles.insert(checked);
        if (m.rfind("Warning: the run stopped at CSS", 0) == 0) o.warnings.push_back(m);
      };
      o.run = sim.Run(options);
      return o;
    }
  }  // namespace

  TEST(CSSDistanceDecks, NotEstimatedRightAfterAJump) {
    Simulation sim = Simulation::FromFile("../examples/decks/Illustrative_TSA_13X.json");
    ASSERT_EQ(sim.Data().settings.css.method_, CSSMethod::MASS_BALANCE_COMPONENT);  // GUARD
    sim.Data().settings.acc_css.method_ = AccelerationMethod::AITKEN;
    sim.Data().settings.max_cycles_to_keep_ = 40;  // keep every cycle for the history
    sim.Data().cycle.number_of_cycles_ = 40;
    const Observed o = RunObserved(sim);
    ASSERT_TRUE(o.run.Ok()) << o.run.message;
    ASSERT_TRUE(o.run.css_reached);
    ASSERT_FALSE(o.jump_cycles.empty()) << "GUARD: Aitken made no jump";
    EXPECT_NE(o.stop_line.find("est. distance to CSS"), std::string::npos) << o.stop_line;

    const auto history = sim.CSSHistory();
    int checked_after_jump = 0;
    for (const CSSRecord& r : history) {
      // A cycle within kDistanceMinChanges - 1 cycles of the cycle that
      // started from a jump has too few natural changes since it.
      for (int j : o.jump_cycles) {
        if (r.cycle >= j && r.cycle < j + CyclicSteadyState::kDistanceMinChanges) {
          EXPECT_FALSE(r.distance_estimated) << "cycle " << r.cycle << " after the jump into " << j;
          EXPECT_EQ(r.distance_note, "too few cycles since the last acceleration jump") << r.cycle;
          EXPECT_TRUE(std::isnan(r.distance) && std::isnan(r.distance_ratio)) << r.cycle;
          ++checked_after_jump;
        }
      }
      if (r.distance_estimated) {
        EXPECT_TRUE(std::isfinite(r.distance) && r.distance > 0.0) << r.cycle;
        EXPECT_GE(r.convergence_rate, 0.0);
        EXPECT_LE(r.convergence_rate, CyclicSteadyState::kDistanceMaxRate);
        EXPECT_NEAR(r.distance_ratio, r.distance / r.state_tolerance, 1e-12 * r.distance_ratio);
      } else {
        EXPECT_FALSE(r.distance_note.empty()) << r.cycle;
      }
    }
    EXPECT_GT(checked_after_jump, 0);

    const CSSRecord& stop = history.back();
    ASSERT_EQ(stop.cycle, o.run.css_cycle);
    std::printf("  Illustrative_TSA_13X + Aitken: stop at cycle %d, %zu jumps, %s\n", stop.cycle,
                o.jump_cycles.size(),
                stop.distance_estimated ? ("d/eps " + std::to_string(stop.distance_ratio)).c_str()
                                        : ("not estimated: " + stop.distance_note).c_str());
    // A fast TSA deck ends well within the warning threshold.
    EXPECT_TRUE(o.warnings.empty());
  }

  TEST(CSSDistanceDecks, Perez2019EReportsTheRemainingDistance) {
    Simulation sim = Simulation::FromFile("../examples/decks/Perez2019_E.json");
    ASSERT_EQ(sim.Data().settings.acc_css.method_, AccelerationMethod::NONE);  // GUARD: plain run
    const Observed o = RunObserved(sim);
    ASSERT_TRUE(o.run.Ok()) << o.run.message;
    ASSERT_TRUE(o.run.css_reached);
    // The diagnostic does not move the stop (the plain stop of the reference
    // verification).
    EXPECT_EQ(o.run.css_cycle, 205);

    const CSSRecord stop = sim.CSSHistory().back();
    ASSERT_EQ(stop.cycle, o.run.css_cycle);
    ASSERT_TRUE(stop.distance_estimated) << stop.distance_note;
    std::printf("  Perez2019_E: stop at cycle %d, rate %.4f, distance %.3e, d/eps %.1f "
                "(measured against the long reference: 37)\n",
                stop.cycle, stop.convergence_rate, stop.distance, stop.distance_ratio);
    EXPECT_NE(o.stop_line.find("x tolerance"), std::string::npos) << o.stop_line;
    // Same size as the measured 37: within a factor of about 3 either way.
    EXPECT_GT(stop.distance_ratio, 12.0);
    EXPECT_LT(stop.distance_ratio, 100.0);
    EXPECT_GT(stop.convergence_rate, 0.8);
    EXPECT_TRUE(o.warnings.empty());
  }

}  // namespace adsim
