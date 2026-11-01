// test_css_status.cpp
//
// Regression tests for L-13 (and a null-boundary guard found alongside it):
//
//  * CoupledSolver::IntegrateCycle stops at the cycle where CSS is reached,
//    which can be before number_of_cycles_ - 1. The desktop GUI used to decide
//    "At CSS" vs "Completed" by re-checking CSS at the LAST configured cycle,
//    which has no log data when the run stopped early, so a run that reached
//    CSS was reported as "Completed". The solver now exposes the CSS outcome
//    of its last run (GetCSSCycleIndex / ReachedCSS).
//  * Cycle::PerformChecks must reject a step with a null top/bottom boundary
//    condition with an error instead of dereferencing it.

#include <string>

#include "../src/pch.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"

namespace adsim {

  // The shared feed/hold deck (UserDataTest::AddFeedHoldDeck).
  class CSSStatusTest : public UserDataTest {
  protected:
    static constexpr int kNumberOfCycles = 5;

    void SetUp() override {
      UserDataTest::SetUp();

      AddFeedHoldDeck(kNumberOfCycles);
    }

    // Highest cycle index present in the solver's equation log.
    static int LastLoggedCycle(const CoupledSolver& solver) {
      const auto& cycles = solver.GetLog().sys_data.cycle_index_set_;
      return cycles.empty() ? -1 : *cycles.rbegin();
    }
  };

  // A run that reaches CSS before the last configured cycle must report the
  // cycle it stopped at, even though the final configured cycle never ran.
  TEST_F(CSSStatusTest, ReportsCSSReachedBeforeLastCycle) {
    auto& css = data.settings.css;
    css.method_ = CSSMethod::MOLAR_PROFILE;
    css.min_cycles_ = 2;
    css.molar_absolute_tolerance = 1e12; // any profile change counts as converged
    css.molar_relative_tolerance = 0.0;

    CoupledSolver solver{data};
    int flag = solver.IntegrateCycle(data.cycle, false);
    ASSERT_GE(flag, 0) << "Solver did not complete";

    EXPECT_TRUE(solver.ReachedCSS());
    EXPECT_EQ(solver.GetCSSCycleIndex(), 1);    // min_cycles_ = 2 -> 0-based cycle 1
    EXPECT_EQ(LastLoggedCycle(solver), 1);      // run stopped there

    // The check the GUI used to make (CSS at the last configured cycle) says
    // "not at CSS" for this run, because that cycle never ran.
    EXPECT_FALSE(css.IsAtCSS(solver.GetLog(), kNumberOfCycles - 1));
  }

  // Without a CSS method the run goes to the last cycle and reports no CSS.
  TEST_F(CSSStatusTest, ReportsNoCSSWhenCheckDisabled) {
    data.settings.css.method_ = CSSMethod::NONE;

    CoupledSolver solver{data};
    int flag = solver.IntegrateCycle(data.cycle, false);
    ASSERT_GE(flag, 0) << "Solver did not complete";

    EXPECT_FALSE(solver.ReachedCSS());
    EXPECT_EQ(solver.GetCSSCycleIndex(), -1);
    EXPECT_EQ(LastLoggedCycle(solver), kNumberOfCycles - 1);
  }

  TEST_F(CSSStatusTest, NullTopBoundaryConditionRejectedByChecks) {
    data.cycle.GetStep("Hold").top_boundary_condition_.reset();
    std::string error;
    EXPECT_FALSE(data.cycle.PerformChecks(data.fluid, error));
    EXPECT_NE(error.find("Top boundary condition of step 'Hold' is null"), std::string::npos) << error;
  }

  TEST_F(CSSStatusTest, NullBottomBoundaryConditionRejectedByChecks) {
    data.cycle.GetStep("Feed").bottom_boundary_condition_.reset();
    std::string error;
    EXPECT_FALSE(data.PerformChecks(error));
    EXPECT_NE(error.find("Bottom boundary condition of step 'Feed' is null"), std::string::npos) << error;
  }

} // namespace adsim
