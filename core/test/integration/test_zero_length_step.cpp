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

  // The shared feed/hold deck (UserDataTest::AddFeedHoldDeck).
  class ZeroLengthStepTest : public UserDataTest {
  protected:
    void SetUp() override {
      UserDataTest::SetUp();

      AddFeedHoldDeck(2);
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
