// test_recovered_steps.cpp
//
// Recovered integrator failures are quiet and counted: a run whose steps IDA
// fails (an initial step far too large, with a small error-test budget) but
// the engine's step retry recovers writes nothing to stderr, logs each
// recovery at MEDIUM with IDA's own message (nothing at LOW), counts it in
// IntegrationStats::recovered_steps, and the STEP_RETRIES measure reads the
// same count. A run whose retry budget runs out logs the failure at LOW with
// IDA's last message, still without writing to stderr.

#include <string>
#include <vector>

#include "../src/pch.h"
#include "adsim/solver/run_health.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/simulation.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"

namespace adsim {

  class RecoveredStepsTest : public UserDataTest {
  protected:
    static constexpr double kAtm = 101325.0;

    // The run-health fixture's two-step cycle (test_run_health.cpp): an
    // upward feed, then a downward purge at a slightly lower pressure.
    void SetUp() override {
      UserDataTest::SetUp();
      const std::unordered_map<std::string, double> feed{ {"CarbonDioxide", 0.15}, {"Nitrogen", 0.85} };
      AddComponents({ "CarbonDioxide", "Nitrogen" });
      AddInitialConditions(kAtm, 293, 293, feed);
      auto [layer_id, ads_id] = AddLayer("Adsorbent");
      GetLayer(layer_id).domain_.n_cells = 10;
      AddIsotherm({ "Nitrogen" }, ads_id, std::make_unique<HenryIsotherm>("Nitrogen"));
      AddIsotherm({ "CarbonDioxide" }, ads_id, std::make_unique<LangmuirIsotherm>("CarbonDioxide"));
      AddLinearDrivingForceModel(ads_id, "CarbonDioxide", 0.5);
      AddLinearDrivingForceModel(ads_id, "Nitrogen", 0.5);
      AddConstantInternalHeatTransferCoefficientModel(0.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);
      AddStep("Feed", 20.0);
      {
        std::vector<std::unique_ptr<IBoundaryStream>> streams;
        streams.push_back(CreateNoValveInletStream(CreateVelocityBC(0.05), feed, 293));
        AddFlowBC(std::move(streams), BoundaryLocation::BOTTOM, BoundaryDirection::INLET, "Feed");
      }
      AddPressureOutletBC(CreateConstantBoundaryValue(kAtm), BoundaryLocation::TOP, "Feed");
      AddStep("Purge", 10.0);
      {
        std::vector<std::unique_ptr<IBoundaryStream>> streams;
        streams.push_back(CreateNoValveInletStream(CreateVelocityBC(0.02), feed, 293));
        AddFlowBC(std::move(streams), BoundaryLocation::TOP, BoundaryDirection::INLET, "Purge");
      }
      AddPressureOutletBC(CreateConstantBoundaryValue(0.95 * kAtm), BoundaryLocation::BOTTOM, "Purge");
      ConfigureCycle(1);
      data.settings.css.method_ = CSSMethod::NONE;
    }

    struct Message { LogLevel level; std::string text; };

    static bool Contains(const std::string& s, const std::string& part) { return s.find(part) != std::string::npos; }
  };

  TEST_F(RecoveredStepsTest, RecoveredFailuresAreCountedLoggedAtMediumAndQuiet) {
    // Every step starts with a 1 s initial step (far too large for the switch
    // of boundary conditions), and IDA may fail the error test only three
    // times per step: IDA gives up on the first step of each cycle step, and
    // the engine's retry (a restart at half the failed step) gets past it.
    data.settings.init_step = ExpressionParam{1.0};
    data.settings.max_err_test_fails = 3;
    Simulation sim;
    sim.SetData(std::move(data));
    std::vector<Message> messages;
    RunOptions options;
    options.log_level = LogLevel::HIGH;
    options.on_message = [&messages](LogLevel l, const std::string& m) { messages.push_back({ l, m }); };

    testing::internal::CaptureStderr();
    const RunResult r = sim.Run(options);
    const std::string err = testing::internal::GetCapturedStderr();

    EXPECT_TRUE(r.Ok()) << r.message;
    EXPECT_EQ(err, "");  // nothing from SUNDIALS reaches stderr
    EXPECT_GT(r.stats.recovered_steps, 0);

    long medium = 0;
    for (const auto& m : messages) {
      if (!Contains(m.text, "recovered step")) continue;
      EXPECT_EQ(m.level, LogLevel::MEDIUM) << m.text;
      EXPECT_TRUE(Contains(m.text, "IDA: ")) << m.text;  // IDA's own message is quoted
      ++medium;
    }
    EXPECT_EQ(medium, r.stats.recovered_steps);
    long routed = 0;  // IDA's messages themselves, at HIGH
    for (const auto& m : messages) {
      if (m.level == LogLevel::LOW) EXPECT_FALSE(Contains(m.text, "IDA")) << m.text;
      if (m.level == LogLevel::HIGH && Contains(m.text, "  IDA: ")) ++routed;
    }
    EXPECT_GE(routed, r.stats.recovered_steps);

    // STEP_RETRIES reads the same count: its measure is the IDA failures plus
    // the recovered steps (one cycle, so the cycle's counters are the run's).
    const auto& report = sim.RunHealth();
    ASSERT_EQ(report.cycles.size(), 1u);
    EXPECT_EQ(report.cycles[0].rejected,
              r.stats.num_err_test_fails + r.stats.num_conv_fails + r.stats.recovered_steps);
    EXPECT_EQ(report.cycles[0].steps, r.stats.num_steps);
  }

  TEST_F(RecoveredStepsTest, ChunkedRunLogsTheSameRecoveries) {
    data.settings.init_step = ExpressionParam{1.0};
    data.settings.max_err_test_fails = 3;
    Simulation sim;
    sim.SetData(std::move(data));
    const RunResult r = sim.Run();  // blocking, for the count
    ASSERT_TRUE(r.Ok()) << r.message;
    ASSERT_GT(r.stats.recovered_steps, 0);

    std::vector<Message> messages;
    testing::internal::CaptureStderr();
    ASSERT_TRUE(sim.Start(LogLevel::MEDIUM, [&messages](LogLevel l, const std::string& m) { messages.push_back({ l, m }); }));
    RunStatus st = RunStatus::Running;
    while (st == RunStatus::Running) st = sim.Advance(1000);
    const std::string err = testing::internal::GetCapturedStderr();
    EXPECT_EQ(st, RunStatus::Completed);
    EXPECT_EQ(err, "");
    long medium = 0;
    for (const auto& m : messages)
      if (Contains(m.text, "recovered step")) {
        EXPECT_EQ(m.level, LogLevel::MEDIUM) << m.text;
        ++medium;
      }
    EXPECT_EQ(medium, r.stats.recovered_steps);
  }

  TEST_F(RecoveredStepsTest, HealthyRunRecoversNothing) {
    Simulation sim;
    sim.SetData(std::move(data));
    testing::internal::CaptureStderr();
    const RunResult r = sim.Run();
    const std::string err = testing::internal::GetCapturedStderr();
    EXPECT_TRUE(r.Ok()) << r.message;
    EXPECT_EQ(err, "");
    EXPECT_EQ(r.stats.recovered_steps, 0);
  }

  TEST_F(RecoveredStepsTest, ExhaustedRetryBudgetLogsAtLowWithIdasMessage) {
    // One Newton iteration at an absurd tolerance and one convergence failure
    // allowed per step: no step converges, so the retries run out and the run
    // fails.
    data.settings.nonlin_conv_coef = ExpressionParam{1e-12};
    data.settings.max_nonlin_iters = 1;
    data.settings.max_conv_fails = 1;
    Simulation sim;
    sim.SetData(std::move(data));
    std::vector<Message> messages;
    RunOptions options;
    options.log_level = LogLevel::MEDIUM;
    options.on_message = [&messages](LogLevel l, const std::string& m) { messages.push_back({ l, m }); };

    testing::internal::CaptureStderr();
    const RunResult r = sim.Run(options);
    const std::string err = testing::internal::GetCapturedStderr();

    EXPECT_FALSE(r.Ok());
    EXPECT_EQ(err, "");
    bool failure_logged = false;
    for (const auto& m : messages)
      if (m.level == LogLevel::LOW && Contains(m.text, "failed to converge after")) {
        failure_logged = true;
        EXPECT_TRUE(Contains(m.text, "after 4 retries")) << m.text;
        EXPECT_TRUE(Contains(m.text, "IDA: ")) << m.text;  // IDA's last message is quoted
      }
    EXPECT_TRUE(failure_logged);
  }

}  // namespace adsim
