// test_profile_runtime_values.cpp
//
// Regression tests for L-05: pressure-driven profile parents (PressureBoundaryCondition,
// ValveBoundaryStream) call InitialiseStartValue / InitialiseStepTime at every step
// start (start = current bed-end pressure, ramp = whole step). Those used to assign
// through ExpressionParam::operator=, which rewrote the user's expression strings
// (start_value_ / ramp_time_ / step_time_) and persisted the runtime values into the
// saved deck. The runtime values now live in transient members, so:
//
//  * after a run the user's expressions are untouched (formulas survive);
//  * the profile values the solver sees are unchanged — the user's start/ramp are
//    still ignored for pressure/valve parents, exactly as before.

#include <string>
#include <vector>

#include "../src/pch.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/boundary_value.h"
#include "adsim/boundary_conditions/step.h"

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"

namespace adsim {

  // Three-step deck: pressurise through a linear pressure BC, blow down through a
  // valve with an exponential downstream-pressure profile, then hold.
  class ProfileRuntimeValuesTest : public UserDataTest {
  protected:
    LinearBoundaryValue* press_profile = nullptr;       // pressure-BC parent
    ExponentialBoundaryValue* valve_profile = nullptr;  // valve-stream parent

    void SetUp() override {
      UserDataTest::SetUp();

      AddComponents({ "Oxygen", "Nitrogen" });

      auto [layer_id, ads_id] = AddLayer("13X");
      GetLayer(layer_id).domain_.n_cells = 10;
      AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0, layer_id);

      AddInitialConditions(101325, 298, 298, { {"Oxygen", 0.0}, {"Nitrogen", 1.0} });

      AddStep("Press", 5.0);
      {
        Step& s = data.cycle.GetStep("Press");
        auto profile = CreateLinearBoundaryValue(2e5, 1.0);
        press_profile = profile.get();
        auto bc = CreatePressureInletBoundary(std::move(profile),
          { {"Oxygen", 0.21}, {"Nitrogen", 0.79} }, 298);
        bc->direction_ = BoundaryDirection::INLET;
        bc->location_ = BoundaryLocation::BOTTOM;
        s.bottom_boundary_condition_ = std::move(bc);
      }
      AddClosedBC(BoundaryLocation::TOP, "Press");

      AddStep("Blowdown", 5.0);
      AddClosedBC(BoundaryLocation::BOTTOM, "Blowdown");
      {
        Step& s = data.cycle.GetStep("Blowdown");
        auto profile = CreateExponentialBoundaryValue(101325, 1.0);
        valve_profile = profile.get();
        auto bc = std::make_unique<FlowBoundaryCondition>();
        bc->direction_ = BoundaryDirection::OUTLET;
        bc->location_ = BoundaryLocation::TOP;
        bc->streams_.push_back(CreateValveOutletStream(std::move(profile), 0.03));
        s.top_boundary_condition_ = std::move(bc);
      }

      AddStep("Hold", 5.0);
      AddClosedBC(BoundaryLocation::BOTTOM, "Hold");
      AddClosedBC(BoundaryLocation::TOP, "Hold");

      ConfigureCycle(2);
      data.settings.max_step = 1.0;
    }

    // Give the user-facing start/ramp fields formulas that would visibly change the
    // profile if they were (wrongly) consumed for a pressure/valve parent.
    void SetUserFormulas(const std::string& start, const std::string& ramp) {
      EvalContext ctx;
      press_profile->start_value_.SetExpr(start);
      press_profile->ramp_time_.SetExpr(ramp);
      valve_profile->start_value_.SetExpr(start);
      valve_profile->step_time_.SetExpr(ramp);
      ASSERT_TRUE(press_profile->start_value_.Evaluate(ctx));
      ASSERT_TRUE(press_profile->ramp_time_.Evaluate(ctx));
      ASSERT_TRUE(valve_profile->start_value_.Evaluate(ctx));
      ASSERT_TRUE(valve_profile->step_time_.Evaluate(ctx));
    }

    std::vector<double> Run(CoupledSolver& solver) {
      const int retval = solver.IntegrateCycle(data.cycle, false);
      EXPECT_GE(retval, 0) << "Solver did not complete";
      const int n = data.state.GetSystemSize();
      std::vector<double> y(n, 0.0), yp(n, 0.0);
      data.state.PackBanded(y.data(), yp.data());
      return y;
    }
  };

  TEST_F(ProfileRuntimeValuesTest, RunLeavesUserExpressionsUntouched) {
    SetUserFormulas("3e5 / 2", "2 * 1.5");
    std::string error;
    ASSERT_TRUE(data.cycle.PerformChecks(data.fluid, error)) << error;

    CoupledSolver solver{data};
    Run(solver);

    EXPECT_EQ(press_profile->start_value_.GetExpr(), "3e5 / 2");
    EXPECT_EQ(press_profile->ramp_time_.GetExpr(), "2 * 1.5");
    EXPECT_EQ(valve_profile->start_value_.GetExpr(), "3e5 / 2");
    EXPECT_EQ(valve_profile->step_time_.GetExpr(), "2 * 1.5");
    EXPECT_DOUBLE_EQ(press_profile->start_value_.Get(), 1.5e5);
    EXPECT_DOUBLE_EQ(press_profile->ramp_time_.Get(), 3.0);
    EXPECT_DOUBLE_EQ(valve_profile->start_value_.Get(), 1.5e5);
    EXPECT_DOUBLE_EQ(valve_profile->step_time_.Get(), 3.0);
  }

  // The user's start/ramp have never driven a pressure/valve profile (they were
  // overwritten at step start); they still must not. Two runs that differ only in
  // those fields produce bit-identical states.
  TEST_F(ProfileRuntimeValuesTest, UserStartAndRampDoNotAffectPressureProfiles) {
    std::string error;
    ASSERT_TRUE(data.cycle.PerformChecks(data.fluid, error)) << error;
    CoupledSolver solver{data};

    SetUserFormulas("1", "1");
    const std::vector<double> first = Run(solver);

    SetUserFormulas("9e5", "0.25");
    const std::vector<double> second = Run(solver);

    ASSERT_EQ(first.size(), second.size());
    ASSERT_FALSE(first.empty());
    for (size_t i = 0; i < first.size(); ++i)
      ASSERT_EQ(first[i], second[i]) << "y[" << i << "]";
  }

} // namespace adsim
