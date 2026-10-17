// test_step_retry.cpp
//
// Regression tests for L-02: the "reduced-step retry" after a recoverable IDA
// failure called IDASetInitStep(h_retry). IDA only reads hin before its first
// accepted step (nst == 0), so after that the retry did not reduce anything
// and its log line was misleading; worse, the tiny hin (clamped to min_step,
// 1e-15 s by default) stayed in IDA's memory and was applied by the next
// IDAReInit, i.e. the first step after a CSS-acceleration re-seed.
//
//  * A re-seed (ValidateAndCommitReseed / Reinitialise) after a retry starts
//    with the configured init_step, not the retry's reduced step.
//  * Before the first accepted step a retry restarts at half the failed step.
//  * After an accepted step a retry leaves hin alone and continues from IDA's
//    own step; once that step sits at the min_step floor it restarts from a
//    consistent initial condition instead (F-01).

#include <cmath>
#include <string>

#include "../src/pch.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"

namespace adsim {

  // Single feed step on a small 13X bed, driven directly through SolverData.
  class StepRetryTest : public UserDataTest {
  protected:
    static constexpr double kStepTime = 5.0;

    void SetUp() override {
      UserDataTest::SetUp();

      AddComponents({ "Oxygen", "Nitrogen" });

      auto [layer_id, ads_id] = AddLayer("13X");
      GetLayer(layer_id).domain_.n_cells = 10;
      AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0, layer_id);

      AddInitialConditions(101325, 298, 298, { {"Oxygen", 0.0}, {"Nitrogen", 1.0} });

      AddStep("Feed", kStepTime);
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

      ConfigureCycle(1);
      data.settings.max_step = 1.0;
    }

    // Mirrors CoupledSolver::IntegrateCycle / Initialise / IntegrateStep up to
    // the first IDASolve call.
    void Start(SolverData& sd) {
      std::string err;
      ASSERT_TRUE(data.PerformChecks(err)) << err;
      Step& first = data.cycle.step_list_[0];
      first.ClearCycleSummary();
      data.Initialise();
      data.state.SetCurrentStep(&first);
      data.step_start_time_ = 0.0;
      first.Initialise(data.state);
      ASSERT_EQ(sd.Initialise(data, err), IDA_SUCCESS) << err;
      first.ClearOutletConditions();
      data.cycle.ResolveBufferTankSources();
      sd.SetStopTime(kStepTime);
      sd.OnStepStart(0);
    }

    // Takes n accepted internal steps; fails the test on an IDA error.
    void Advance(SolverData& sd, int n) {
      double t = sd.GetTime();
      for (int i = 0; i < n; ++i) {
        const int flag = sd.OneStep(kStepTime, t);
        ASSERT_GE(flag, 0) << "IDA failed at internal step " << i << " (t=" << t << ")";
        if (flag == IDA_TSTOP_RETURN) break;
      }
    }
  };

  // The defect itself: a hin left behind by a retry is consumed by the next
  // IDAReInit. Simulate the old retry (SetInitStep(min_step) after accepted
  // steps), then re-seed: the first step after the re-seed must use init_step.
  TEST_F(StepRetryTest, ReseedAfterRetryStartsAtConfiguredInitStep) {
    const double init_step = data.settings.init_step.Get();
    const double min_step = data.settings.min_step.Get();
    ASSERT_GT(init_step, min_step);

    SolverData sd;
    Start(sd);
    Advance(sd, 1);
    EXPECT_DOUBLE_EQ(sd.GetActualInitStep(), init_step);
    Advance(sd, 20);
    ASSERT_GT(sd.GetNumSteps(), 0);

    // What the old reduced-step retry left in IDA's memory.
    ASSERT_EQ(sd.SetInitStep(min_step), IDA_SUCCESS);

    // CSS-acceleration commit path: re-seed from the current state.
    sd.Unpack(data);
    ASSERT_TRUE(sd.ValidateAndCommitReseed(data, 0.01));
    ASSERT_EQ(sd.GetNumSteps(), 0);
    Advance(sd, 1);
    EXPECT_DOUBLE_EQ(sd.GetActualInitStep(), init_step)
        << "re-seed started from a stale retry step";

    // Reinitialise path.
    Advance(sd, 20);
    ASSERT_EQ(sd.SetInitStep(min_step), IDA_SUCCESS);
    sd.Unpack(data);
    ASSERT_EQ(sd.Reinitialise(data), IDA_SUCCESS);
    ASSERT_EQ(sd.GetNumSteps(), 0);
    Advance(sd, 1);
    EXPECT_DOUBLE_EQ(sd.GetActualInitStep(), init_step)
        << "re-initialisation started from a stale retry step";
  }

  // A genuine first-step failure: an oversized first step with a tiny
  // error-test budget fails before any step is accepted. The retry restarts
  // at half the failed step, and a later re-seed is back at init_step.
  TEST_F(StepRetryTest, FirstStepFailureRestartsAtHalfStepAndDoesNotLeak) {
    data.settings.init_step = 1.0;
    data.settings.reltol = 1e-10;
    data.settings.max_err_test_fails = 3;
    const double min_step = data.settings.min_step.Get();

    SolverData sd;
    Start(sd);

    double t = sd.GetTime();
    const int flag = sd.OneStep(kStepTime, t);
    ASSERT_EQ(flag, IDA_ERR_FAIL) << "deck no longer fails its first step";
    ASSERT_EQ(sd.GetNumSteps(), 0);
    EXPECT_DOUBLE_EQ(sd.GetActualInitStep(), 1.0);

    const SolverData::StepRetry retry = sd.PrepareStepRetry(flag, min_step);
    ASSERT_TRUE(retry.retry);
    EXPECT_TRUE(retry.restarted);
    EXPECT_GT(retry.h, min_step);
    EXPECT_LT(retry.h, 0.5);

    // The retried first step really starts from the reduced step.
    sd.OneStep(kStepTime, t);
    EXPECT_DOUBLE_EQ(sd.GetActualInitStep(), retry.h);

    // Configure a sane deck for the rest of the run, then re-seed: the retry's
    // hin must not be reused.
    data.settings.init_step = 1e-10;
    data.settings.reltol = 1e-6;
    sd.Unpack(data);
    ASSERT_TRUE(sd.ValidateAndCommitReseed(data, 0.01));
    double t_after = sd.GetTime();
    sd.OneStep(kStepTime, t_after);
    EXPECT_DOUBLE_EQ(sd.GetActualInitStep(), 1e-10);
  }

  // After an accepted step IDA ignores hin: a retry continues from IDA's own
  // step without re-initialising, and is refused at the min_step floor.
  TEST_F(StepRetryTest, RetryAfterAcceptedStepContinuesFromIdaStep) {
    SolverData sd;
    Start(sd);
    Advance(sd, 20);
    const long int nst = sd.GetNumSteps();
    ASSERT_GT(nst, 0);
    const double h_ida = std::abs(sd.GetLastStep());
    ASSERT_GT(h_ida, 0.0);

    const SolverData::StepRetry retry = sd.PrepareStepRetry(IDA_ERR_FAIL, data.settings.min_step.Get());
    EXPECT_TRUE(retry.retry);
    EXPECT_FALSE(retry.restarted);
    EXPECT_GT(retry.h, 0.0);
    EXPECT_EQ(sd.GetNumSteps(), nst) << "a continuation retry must not re-initialise";

    // A linear-solver failure is not a step-size failure; still retried by
    // continuing.
    const SolverData::StepRetry lsetup = sd.PrepareStepRetry(IDA_LSETUP_FAIL, retry.h);
    EXPECT_TRUE(lsetup.retry);
    EXPECT_FALSE(lsetup.restarted);
    EXPECT_EQ(sd.GetNumSteps(), nst);

    // A step-size failure with the step already at the floor: continuing would
    // repeat the failed attempt, so the integrator restarts from a consistent
    // IC at the last accepted point instead (F-01), and integrates on.
    const double t_fail = sd.GetTime();
    const SolverData::StepRetry restart = sd.PrepareStepRetry(IDA_ERR_FAIL, retry.h);
    EXPECT_TRUE(restart.retry);
    EXPECT_TRUE(restart.restarted);
    EXPECT_TRUE(restart.consistent_ic);
    EXPECT_EQ(sd.GetNumSteps(), 0) << "the restart must re-initialise";
    EXPECT_DOUBLE_EQ(sd.GetTime(), t_fail);
    Advance(sd, 5);
    EXPECT_GT(sd.GetNumSteps(), 0);
    EXPECT_DOUBLE_EQ(sd.GetActualInitStep(), data.settings.init_step.Get());
  }

  // F-01: cycle steps run back to back without re-initialising IDA, so a step
  // that switches a boundary condition makes the algebraic ghost-cell values
  // jump. IDA's error test only passes a step across the jump that is ~1e-8 of
  // the last step before the boundary, and the next one has to be smaller
  // still (the history reads the jump as a huge algebraic derivative). When
  // the last step before the boundary is a short remainder clipped by tstop,
  // that falls below min_step and the run failed with "error test failed
  // repeatedly or with |h| = hmin" (CSSAcceleration.ReducesCyclesToCSS at
  // t = 5010 s). The retry must restart from a consistent IC and carry on.
  TEST_F(StepRetryTest, BoundarySwitchAfterShortStepIsCrossed) {
    // Second step: bottom closed, top switched to a 3 bar pressure inlet.
    AddStep("Pressurise", kStepTime);
    Step& press = data.cycle.GetStep("Pressurise");
    auto top = std::make_unique<PressureBoundaryCondition>();
    top->direction_ = BoundaryDirection::INLET;
    top->location_ = BoundaryLocation::TOP;
    top->pressure_ = CreateConstantBoundaryValue(3.0e5);
    top->temperature_ = 298.0;
    top->mole_fraction_["Oxygen"] = 0.21;
    top->mole_fraction_["Nitrogen"] = 0.79;
    press.top_boundary_condition_ = std::move(top);
    AddClosedBC(BoundaryLocation::BOTTOM, "Pressurise");

    const double min_step = data.settings.min_step.Get();
    SolverData sd;
    Start(sd);
    Advance(sd, 30);

    // End the feed step with a 1e-6 s remainder clipped by tstop. With it the
    // first step across the switch (~1e-8 of it) is still accepted, but leaves a
    // ghost-pressure derivative of ~jump / h in the history, and the next step
    // collapses to min_step. A restart that kept that derivative would
    // collapse again, so this also pins the zeroed algebraic derivatives.
    double t = sd.GetTime();
    ASSERT_LT(t, kStepTime - 0.5);
    const double t_boundary = t + 1e-6;
    ASSERT_EQ(sd.SetStopTime(t_boundary), IDA_SUCCESS);
    ASSERT_EQ(sd.OneStep(kStepTime, t), IDA_TSTOP_RETURN);
    ASSERT_LT(std::abs(sd.GetLastStep()), 2e-6) << "the last step is not a short remainder";

    // Switch to the second step, as CoupledSolver::IntegrateStep does.
    data.step_start_time_ = t;
    press.Initialise(data.state);
    const double t_final = t + kStepTime;
    ASSERT_EQ(sd.SetStopTime(t_final), IDA_SUCCESS);
    sd.OnStepStart(0);

    const double t_target = t + 0.5;
    int consistent_restarts = 0;
    int retries = 0;
    while (t < t_target) {
      const int flag = sd.OneStep(t_final, t);
      if (flag >= 0) {
        retries = 0;
        continue;
      }
      ASSERT_LT(++retries, 5) << "retry budget exhausted at t = " << t;
      const SolverData::StepRetry retry = sd.PrepareStepRetry(flag, min_step);
      ASSERT_TRUE(retry.retry) << "IDA could not cross the boundary switch (flag " << flag
                               << ", t = " << t << ")";
      if (retry.consistent_ic) ++consistent_restarts;
    }
    EXPECT_GE(t, t_target);
    // GUARD: the switch must actually have collapsed the step to the floor, or
    // this test does not exercise the restart.
    EXPECT_GE(consistent_restarts, 1);
  }

} // namespace adsim
