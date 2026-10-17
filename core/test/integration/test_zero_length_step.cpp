// test_zero_length_step.cpp
//
// Regression tests for E-01: a step with zero duration used to pass
// UserData::PerformChecks, after which the chunked (web) solver spun forever —
// its integration loop `while (t_curr < t_final)` never ran, so the step never
// completed and every Integrate() chunk returned IN_PROGRESS.
//
//  * Cycle::PerformChecks must reject a step time <= 0 (ABSOLUTE and FRACTION
//    mode).
//  * As a belt-and-braces guard, the chunked solver must complete a zero-length
//    step immediately rather than spin, if one ever reaches it.

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
  class ZeroLengthStepTest : public UserDataTest {
  protected:
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

      ConfigureCycle(2);
      data.settings.max_step = 1.0;
    }
  };

  TEST_F(ZeroLengthStepTest, PositiveStepTimesPassChecks) {
    std::string error;
    EXPECT_TRUE(data.PerformChecks(error)) << error;
  }

  TEST_F(ZeroLengthStepTest, ZeroStepTimeRejectedByChecks) {
    data.cycle.GetStep("Hold").step_time_ = 0.0;
    std::string error;
    EXPECT_FALSE(data.PerformChecks(error));
    EXPECT_NE(error.find("Hold"), std::string::npos) << error;
    EXPECT_NE(error.find("greater than zero"), std::string::npos) << error;
  }

  TEST_F(ZeroLengthStepTest, NegativeStepTimeRejectedByChecks) {
    data.cycle.GetStep("Feed").step_time_ = -1.0;
    std::string error;
    EXPECT_FALSE(data.PerformChecks(error));
    EXPECT_NE(error.find("Feed"), std::string::npos) << error;
  }

  TEST_F(ZeroLengthStepTest, ZeroFractionRejectedInFractionMode) {
    // Fractions sum to 1, so only the zero-duration check can catch this.
    data.cycle.step_time_mode_ = StepTimeMode::FRACTION;
    data.cycle.total_cycle_time_ = 10.0;
    data.cycle.GetStep("Feed").step_fraction_ = 1.0;
    data.cycle.GetStep("Hold").step_fraction_ = 0.0;
    std::string error;
    EXPECT_FALSE(data.PerformChecks(error));
    EXPECT_NE(error.find("Hold"), std::string::npos) << error;
    EXPECT_NE(error.find("greater than zero"), std::string::npos) << error;
  }

  // Defensive solver guard: bypass the checks by zeroing the step after Start()
  // (step times are read when each step initialises). Before the guard, every
  // chunk returned IN_PROGRESS forever once the zero-length step was reached; the
  // chunk cap turns that hang into a test failure.
  TEST_F(ZeroLengthStepTest, ChunkedSolverCompletesZeroLengthStep) {
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;

    CoupledSolverChuncked solver{data};
    solver.Start();
    ASSERT_TRUE(solver.IsActive()) << solver.GetStatus().GetMessage();

    data.cycle.GetStep("Hold").step_time_ = 0.0;

    constexpr int kMaxChunks = 20000;
    int chunks = 0;
    while (solver.IsActive() && chunks < kMaxChunks) {
      solver.Integrate(50);
      ++chunks;
    }

    const ChunkStatus status = solver.GetStatus();
    ASSERT_FALSE(solver.IsActive())
      << "Chunked solver still active after " << chunks << " chunks (zero-length step never completed)";
    EXPECT_NE(status.GetStatus(), SolverStatus::FAILED) << status.GetMessage();
    EXPECT_TRUE(status.GetStatus() == SolverStatus::COMPLETED
                || status.GetStatus() == SolverStatus::AT_CSS) << status.GetMessage();
  }

} // namespace adsim
