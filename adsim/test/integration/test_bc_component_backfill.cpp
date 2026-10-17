// test_bc_component_backfill.cpp
//
// Regression tests for E-15: inlet mole-fraction maps that do not cover every
// fluid component.
//
//  * PressureBoundaryCondition::PerformChecks only summed the entries present,
//    so {N2: 1.0} with an N2+O2 fluid passed and Apply() then threw
//    std::out_of_range from mole_fraction_.at("Oxygen") inside the IDA residual.
//  * IBoundaryStream::PerformChecks only back-filled when the component COUNT
//    differed, so a map with one component swapped for another ({N2, Argon}
//    against {O2, N2}) skipped the back-fill and threw the same way.
//  * An exception escaping the residual unwound through SUNDIALS C frames.
//    UpdateResiduals now catches it, returns an IDA error code and the solver
//    reports the message; the run fails cleanly instead of terminating.

#include <stdexcept>
#include <string>

#include "../src/pch.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"

namespace adsim {

  namespace {
    // Constant pressure that throws once the step-relative time reaches t_throw.
    class ThrowingBoundaryValue : public ConstantBoundaryValue {
    public:
      ThrowingBoundaryValue(double value, double t_throw)
        : ConstantBoundaryValue(value), t_throw_(t_throw) {}
      double GetValue(double time) override {
        if (time >= t_throw_) throw std::out_of_range("injected residual failure");
        return ConstantBoundaryValue::GetValue(time);
      }
    private:
      double t_throw_;
    };
  }

  // Two-step deck on an O2 + N2 fluid. "Feed" has a flow inlet at the bottom and
  // a pressure outlet at the top; "Press" has a pressure inlet at the bottom and
  // a closed top. Tests overwrite the inlet compositions.
  class BCComponentBackfillTest : public UserDataTest {
  protected:
    void SetUp() override {
      UserDataTest::SetUp();

      AddComponents({ "Oxygen", "Nitrogen" });

      auto [layer_id, ads_id] = AddLayer("13X");
      GetLayer(layer_id).domain_.n_cells = 10;
      AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0, layer_id);

      AddInitialConditions(101325, 298, 298, { {"Oxygen", 0.21}, {"Nitrogen", 0.79} });

      AddStep("Feed", 5.0);
      Step& feed = data.cycle.GetStep("Feed");
      auto feed_flow_bc = std::make_unique<FlowBoundaryCondition>();
      feed_flow_bc->direction_ = BoundaryDirection::INLET;
      feed_flow_bc->location_ = BoundaryLocation::BOTTOM;
      feed_flow_bc->streams_.push_back(CreateNoValveInletStream(
        CreateMolarFlowBC(0.1), { {"Oxygen", 0.21}, {"Nitrogen", 0.79} }, 298));
      feed.bottom_boundary_condition_ = std::move(feed_flow_bc);

      auto feed_out = std::make_unique<PressureBoundaryCondition>();
      feed_out->direction_ = BoundaryDirection::OUTLET;
      feed_out->location_ = BoundaryLocation::TOP;
      feed_out->pressure_ = CreateConstantBoundaryValue(101325);
      feed.top_boundary_condition_ = std::move(feed_out);

      AddStep("Press", 5.0);
      Step& press = data.cycle.GetStep("Press");
      auto press_in = std::make_unique<PressureBoundaryCondition>();
      press_in->direction_ = BoundaryDirection::INLET;
      press_in->location_ = BoundaryLocation::BOTTOM;
      press_in->pressure_ = CreateLinearBoundaryValue(101325, 150000, 5.0);
      press_in->temperature_ = 298;
      press_in->mole_fraction_ = { {"Oxygen", 0.21}, {"Nitrogen", 0.79} };
      press.bottom_boundary_condition_ = std::move(press_in);
      AddClosedBC(BoundaryLocation::TOP, "Press");

      ConfigureCycle(1);
    }

    PressureBoundaryCondition& PressInlet() {
      return *dynamic_cast<PressureBoundaryCondition*>(
        data.cycle.GetStep("Press").bottom_boundary_condition_.get());
    }

    PressureBoundaryCondition& FeedOutlet() {
      return *dynamic_cast<PressureBoundaryCondition*>(
        data.cycle.GetStep("Feed").top_boundary_condition_.get());
    }

    IBoundaryStream& FeedStream() {
      auto* fbc = dynamic_cast<FlowBoundaryCondition*>(
        data.cycle.GetStep("Feed").bottom_boundary_condition_.get());
      return *fbc->streams_.front();
    }

    // Runs the non-chunked solver, collecting its log. Returns the IDA flag.
    int Run(std::string& log) {
      CoupledSolver solver{data};
      solver.SetMessageCallback([&log](LogLevel, const std::string& msg) { log += msg + "\n"; });
      return solver.IntegrateCycle(data.cycle);
    }
  };

  // Desktop scenario from the bug report: the step's map was filled before O2
  // was added to the fluid. The check must back-fill O2 = 0 and the run work.
  TEST_F(BCComponentBackfillTest, PressureInletMissingComponentIsBackFilled) {
    PressInlet().mole_fraction_ = { {"Nitrogen", 1.0} };

    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;

    const auto& mf = PressInlet().mole_fraction_;
    ASSERT_EQ(mf.count("Oxygen"), 1u);
    EXPECT_DOUBLE_EQ(mf.at("Oxygen"), 0.0);
    EXPECT_DOUBLE_EQ(mf.at("Nitrogen"), 1.0);

    std::string log;
    int flag = 0;
    EXPECT_NO_THROW(flag = Run(log));
    EXPECT_GE(flag, 0) << log;
  }

  // Entries for components no longer in the fluid are dropped, so the map is
  // keyed exactly by fluid.components_ and the sum check sees the real total.
  TEST_F(BCComponentBackfillTest, PressureInletStaleComponentIsErased) {
    PressInlet().mole_fraction_ = { {"Nitrogen", 1.0}, {"Argon", 0.0} };

    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;

    const auto& mf = PressInlet().mole_fraction_;
    EXPECT_EQ(mf.size(), 2u);
    EXPECT_EQ(mf.count("Argon"), 0u);
    EXPECT_EQ(mf.count("Oxygen"), 1u);
  }

  TEST_F(BCComponentBackfillTest, PressureInletSumStillChecked) {
    PressInlet().mole_fraction_ = { {"Nitrogen", 0.7}, {"Argon", 0.3} };

    std::string error;
    EXPECT_FALSE(data.PerformChecks(error));
    EXPECT_NE(error.find("Sum of mole fractions"), std::string::npos) << error;
  }

  // A stream map with one component swapped for another has the right COUNT, so
  // the old back-fill (only run on a count mismatch) skipped it and the residual
  // threw on the missing O2.
  TEST_F(BCComponentBackfillTest, FlowInletSwappedComponentIsBackFilled) {
    FeedStream().mole_fraction = { {"Nitrogen", 1.0}, {"Argon", 0.0} };

    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;

    const auto& mf = FeedStream().mole_fraction;
    EXPECT_EQ(mf.size(), 2u);
    EXPECT_EQ(mf.count("Argon"), 0u);
    ASSERT_EQ(mf.count("Oxygen"), 1u);
    EXPECT_DOUBLE_EQ(mf.at("Oxygen"), 0.0);
    EXPECT_DOUBLE_EQ(mf.at("Nitrogen"), 1.0);

    std::string log;
    int flag = 0;
    EXPECT_NO_THROW(flag = Run(log));
    EXPECT_GE(flag, 0) << log;
  }

  // Residual-time lookups no longer throw even if a map is edited after the
  // checks ran (e.g. between chunks): a missing component reads as 0. The edit
  // is on the second step so the composition change lands on a step boundary.
  TEST_F(BCComponentBackfillTest, MissingKeyAfterChecksDoesNotThrow) {
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;

    CoupledSolverChuncked solver{data};
    solver.Start();
    ASSERT_TRUE(solver.IsActive()) << solver.GetStatus().GetMessage();

    PressInlet().mole_fraction_ = { {"Nitrogen", 1.0} };

    constexpr int kMaxChunks = 20000;
    int chunks = 0;
    EXPECT_NO_THROW({
      while (solver.IsActive() && chunks < kMaxChunks) {
        solver.Integrate(50);
        ++chunks;
      }
    });
    EXPECT_NE(solver.GetStatus().GetStatus(), SolverStatus::FAILED) << solver.GetStatus().GetMessage();
  }

  // An exception thrown mid-step inside the residual must not escape through
  // IDA: the run fails with IDA_RES_FAIL and the message reaches the log.
  TEST_F(BCComponentBackfillTest, ResidualExceptionFailsRunCleanly) {
    FeedOutlet().pressure_ = std::make_unique<ThrowingBoundaryValue>(101325, 1.0);

    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;

    std::string log;
    int flag = 0;
    EXPECT_NO_THROW(flag = Run(log));
    EXPECT_LT(flag, 0) << log;
    EXPECT_NE(log.find("Residual evaluation threw: injected residual failure"), std::string::npos) << log;
    EXPECT_TRUE(data.residual_error_.empty()) << "message should be consumed by the failure log";
  }

  // Same, but the residual throws during the consistent-IC solve.
  TEST_F(BCComponentBackfillTest, ResidualExceptionDuringInitialisationFailsCleanly) {
    FeedOutlet().pressure_ = std::make_unique<ThrowingBoundaryValue>(101325, 0.0);

    std::string log;
    int flag = 0;
    EXPECT_NO_THROW(flag = Run(log));
    EXPECT_LT(flag, 0) << log;
    EXPECT_NE(log.find("Initialisation failed"), std::string::npos) << log;
    EXPECT_NE(log.find("injected residual failure"), std::string::npos) << log;
  }

  // Chunked (web/GUI) solver: the failure status carries the message.
  TEST_F(BCComponentBackfillTest, ResidualExceptionFailsChunkedRunCleanly) {
    FeedOutlet().pressure_ = std::make_unique<ThrowingBoundaryValue>(101325, 1.0);

    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;

    CoupledSolverChuncked solver{data};
    solver.Start();
    ASSERT_TRUE(solver.IsActive()) << solver.GetStatus().GetMessage();

    constexpr int kMaxChunks = 20000;
    int chunks = 0;
    EXPECT_NO_THROW({
      while (solver.IsActive() && chunks < kMaxChunks) {
        solver.Integrate(50);
        ++chunks;
      }
    });
    const ChunkStatus status = solver.GetStatus();
    EXPECT_EQ(status.GetStatus(), SolverStatus::FAILED);
    EXPECT_NE(status.GetMessage().find("injected residual failure"), std::string::npos)
      << status.GetMessage();
  }

} // namespace adsim
