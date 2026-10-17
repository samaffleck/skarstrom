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

  // Small two-step deck: a short feed step followed by a closed hold step.
  class CSSStatusTest : public UserDataTest {
  protected:
    static constexpr int kNumberOfCycles = 5;

    void SetUp() override {
      UserDataTest::SetUp();

      AddComponents({ "Oxygen", "Nitrogen" });

      auto [layer_id, ads_id] = AddLayer("13X");
      GetLayer(layer_id).domain_.n_cells = 10;
      AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0, layer_id);

      AddInitialConditions(101325, 298, 298, { {"Oxygen", 0.0}, {"Nitrogen", 1.0} });

      AddStep("Feed", 5.0);
      Step& feed = data.cycle.GetStep("Feed");
      auto feed_flow_bc = std::make_unique<FlowBoundaryCondition>();
      feed_flow_bc->direction_ = BoundaryDirection::INLET;
      feed_flow_bc->location_ = BoundaryLocation::BOTTOM;
      feed_flow_bc->streams_.push_back(CreateNoValveInletStream(
        CreateMolarFlowBC(0.1), { {"Oxygen", 0.21}, {"Nitrogen", 0.79} }, 298));
      feed.bottom_boundary_condition_ = std::move(feed_flow_bc);

      auto feed_pressure_out = std::make_unique<PressureBoundaryCondition>();
      feed_pressure_out->direction_ = BoundaryDirection::OUTLET;
      feed_pressure_out->location_ = BoundaryLocation::TOP;
      feed_pressure_out->pressure_ = CreateConstantBoundaryValue(101325);
      feed_pressure_out->mole_fraction_["Oxygen"] = 0.0;
      feed_pressure_out->mole_fraction_["Nitrogen"] = 1.0;
      feed.top_boundary_condition_ = std::move(feed_pressure_out);

      AddStep("Hold", 5.0);
      AddClosedBC(BoundaryLocation::BOTTOM, "Hold");
      AddClosedBC(BoundaryLocation::TOP, "Hold");

      ConfigureCycle(kNumberOfCycles);
      data.settings.max_step = 1.0;
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
