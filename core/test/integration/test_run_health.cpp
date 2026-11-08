// test_run_health.cpp
//
// Run-health checks (run_health.h): a healthy two-step cycle is silent, and
// deliberately broken runs raise each warning: a component source with no
// matching overall source (MOLE_FRACTION_SUM), a leak of one component
// (MASS_BALANCE_STUCK), an absurd nonlinear tolerance (STEP_RETRIES) and a log shorter
// than the diagnostics read (HISTORY_TRUNCATED). The broken runs use the
// manufactured-source verification hook (manufactured_sources.h), which no
// deck can set.

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "../src/pch.h"
#include "adsim/solver/run_health.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/objects/manifests.h"
#include "adsim/simulation.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"

namespace adsim {

  class RunHealthTest : public UserDataTest {
  protected:
    static constexpr double kAtm = 101325.0;
    // Manufactured CO2 sink with no overall sink [mol/(m3 s)]: strong
    // enough to break the sum, weak enough for the run to finish (a source of the
    // same size fails the run).
    static constexpr double kSumSource = -0.02;
    const std::unordered_map<std::string, double> feed{ {"CarbonDioxide", 0.15}, {"Nitrogen", 0.85} };

    void SetUp() override {
      UserDataTest::SetUp();
      AddComponents({ "CarbonDioxide", "Nitrogen" });
      AddInitialConditions(kAtm, 293, 293, { {"CarbonDioxide", 0.15}, {"Nitrogen", 0.85} });
      auto [layer_id, ads_id] = AddLayer("Adsorbent");
      GetLayer(layer_id).domain_.n_cells = 10;
      AddIsotherm({ "Nitrogen" }, ads_id, std::make_unique<HenryIsotherm>("Nitrogen"));
      AddIsotherm({ "CarbonDioxide" }, ads_id, std::make_unique<LangmuirIsotherm>("CarbonDioxide"));
      AddLinearDrivingForceModel(ads_id, "CarbonDioxide", 0.5);
      AddLinearDrivingForceModel(ads_id, "Nitrogen", 0.5);
      AddConstantInternalHeatTransferCoefficientModel(0.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);

      // Upward feed, then a downward purge at a slightly lower pressure.
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
      ConfigureCycle(8);
      data.settings.css.method_ = CSSMethod::NONE;
    }

    // Runs the cycle; returns the solver's LOW-level messages.
    std::vector<std::string> Run(CoupledSolver& solver) {
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      std::vector<std::string> messages;
      solver.SetLogLevel(LogLevel::LOW);
      solver.SetMessageCallback([&messages](LogLevel, const std::string& msg) { messages.push_back(msg); });
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0);
      return messages;
    }

    static const RunHealthWarning* Find(const RunHealthReport& r, RunHealthCheck check) {
      for (const auto& w : r.warnings)
        if (w.check == check) return &w;
      return nullptr;
    }

    static std::string Codes(const RunHealthReport& r) {
      std::string out;
      for (const auto& w : r.warnings) out += w.code + " (" + std::to_string(w.value) + ") ";
      return out;
    }
  };

  TEST_F(RunHealthTest, HealthyCycleIsSilent) {
    CoupledSolver solver{ data };
    const auto messages = Run(solver);
    const auto& report = solver.GetRunHealth();
    EXPECT_TRUE(report.Healthy()) << Codes(report);
    ASSERT_EQ(report.cycles.size(), 8u);
    ASSERT_EQ(report.components.size(), 2u);
    for (const auto& c : report.cycles) {
      EXPECT_LT(c.mole_fraction_sum_error, 1e-6) << "cycle " << c.cycle;
      EXPECT_GT(c.steps, 0) << "cycle " << c.cycle;
      ASSERT_EQ(c.mass_balance_error.size(), 2u);
      for (double e : c.mass_balance_error) EXPECT_LT(std::abs(e), 1e-3) << "cycle " << c.cycle;
    }
    EXPECT_TRUE(std::isnan(report.cycles.front().state_change));  // no previous cycle
    EXPECT_GE(report.cycles.back().state_change, 0.0);
    for (const auto& m : messages) EXPECT_EQ(m.find("mole fractions sum"), std::string::npos) << m;
  }

  TEST_F(RunHealthTest, ComponentSourceBreaksTheMoleFractionSum) {
    // CO2 is removed from its own balance but not from the overall one: the mole
    // fractions no longer sum to 1.
    data.manufactured_sources_.species["CarbonDioxide"] = [](double, double) { return kSumSource; };
    ConfigureCycle(2);
    CoupledSolver solver{ data };
    const auto messages = Run(solver);
    const auto& report = solver.GetRunHealth();
    const RunHealthWarning* w = Find(report, RunHealthCheck::MoleFractionSum);
    ASSERT_NE(w, nullptr) << Codes(report);
    EXPECT_EQ(w->code, "MOLE_FRACTION_SUM");
    EXPECT_GT(w->value, RunHealthThresholds::kMoleFractionSum);
    EXPECT_EQ(w->threshold, RunHealthThresholds::kMoleFractionSum);
    EXPECT_GE(w->cycle, 0);
    EXPECT_GE(w->step, 0);
    EXPECT_FALSE(w->step_name.empty());
    ASSERT_FALSE(w->suspects.empty());
    EXPECT_NE(w->suspects.front().find("advection scheme"), std::string::npos);
    EXPECT_NE(w->hint.find("upwind"), std::string::npos);
    // Logged once at the end of the run, at LOW.
    EXPECT_EQ(std::count(messages.begin(), messages.end(), w->message), 1);
  }

  TEST_F(RunHealthTest, SwitchedOffCheckIsNotReported) {
    data.manufactured_sources_.species["CarbonDioxide"] = [](double, double) { return kSumSource; };
    data.settings.warnings_.mole_fraction_sum = false;
    ConfigureCycle(2);
    CoupledSolver solver{ data };
    Run(solver);
    EXPECT_EQ(Find(solver.GetRunHealth(), RunHealthCheck::MoleFractionSum), nullptr);
    // The measure is still recorded.
    ASSERT_FALSE(solver.GetRunHealth().cycles.empty());
    EXPECT_GT(solver.GetRunHealth().cycles.back().mole_fraction_sum_error, RunHealthThresholds::kMoleFractionSum);
  }

  TEST_F(RunHealthTest, LeakKeepsTheMassBalanceOpen) {
    // CO2 disappears at a constant rate from its own and the overall balance:
    // the composition stays consistent, but the cycle balance never closes.
    data.manufactured_sources_.species["CarbonDioxide"] = [](double, double) { return -0.02; };
    data.manufactured_sources_.total = [](double, double) { return -0.02; };
    ConfigureCycle(12);
    CoupledSolver solver{ data };
    Run(solver);
    const auto& report = solver.GetRunHealth();
    EXPECT_EQ(Find(report, RunHealthCheck::MoleFractionSum), nullptr) << Codes(report);
    const RunHealthWarning* w = Find(report, RunHealthCheck::MassBalanceStuck);
    ASSERT_NE(w, nullptr) << Codes(report);
    EXPECT_EQ(w->code, "MASS_BALANCE_STUCK");
    EXPECT_EQ(w->component, "CarbonDioxide");
    EXPECT_GT(w->value, RunHealthThresholds::kStuckBalance);
    EXPECT_GE(w->cycle, RunHealthThresholds::kStuckWindow - 1);
    EXPECT_LE(report.cycles[w->cycle].state_change, RunHealthThresholds::kStuckStateChange);
  }

  TEST_F(RunHealthTest, AbsurdNonlinearToleranceRaisesStepRetries) {
    // An absurd nonlinear tolerance (1e-6 of the error-test tolerance) with
    // at most two Newton iterations per step: most attempts fail to
    // converge, yet the run finishes. (A tight or loose relative tolerance
    // alone does not do it: IDA adapts its step and stays at or below 0.15.)
    data.settings.nonlin_conv_coef = ExpressionParam{1e-6};
    data.settings.max_nonlin_iters = 2;
    ConfigureCycle(1);
    CoupledSolver solver{ data };
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    solver.IntegrateCycle(data.cycle, false);
    const auto& report = solver.GetRunHealth();
    const RunHealthWarning* w = Find(report, RunHealthCheck::StepRetries);
    ASSERT_NE(w, nullptr) << Codes(report);
    EXPECT_EQ(w->code, "STEP_RETRIES");
    EXPECT_GT(w->value, RunHealthThresholds::kStepRetryRatio);
  }

  TEST_F(RunHealthTest, ShortLogIsReportedAndResultIsUnhealthy) {
    data.settings.max_cycles_to_keep_ = 2;
    Simulation sim;
    sim.SetData(std::move(data));
    const RunResult r = sim.Run();
    EXPECT_TRUE(r.Ok()) << r.message;
    EXPECT_EQ(r.status, RunStatus::Completed);  // the status is unchanged
    EXPECT_FALSE(r.healthy);
    ASSERT_EQ(r.health_warnings.size(), 1u);
    const auto& w = r.health_warnings.front();
    EXPECT_EQ(w.code, "HISTORY_TRUNCATED");
    EXPECT_LT(w.value, w.threshold);
    EXPECT_EQ(w.cycle, -1);
    EXPECT_EQ(sim.RunHealth().cycles.size(), 8u);
    // KPIs, profiles and the other checks are unaffected; switching the check
    // off makes the same run healthy.
    sim.Data().settings.warnings_.history_truncated = false;
    const RunResult again = sim.Run();
    EXPECT_TRUE(again.healthy);
    EXPECT_TRUE(again.health_warnings.empty());
  }

  TEST_F(RunHealthTest, ChunkedRunReportsTheSameWarnings) {
    data.settings.max_cycles_to_keep_ = 2;
    Simulation sim;
    sim.SetData(std::move(data));
    ASSERT_TRUE(sim.Start(LogLevel::LOW));
    RunStatus st = RunStatus::Running;
    while (st == RunStatus::Running) st = sim.Advance(1000);
    EXPECT_EQ(st, RunStatus::Completed);
    const RunResult& r = sim.LastResult();
    EXPECT_FALSE(r.healthy);
    ASSERT_EQ(r.health_warnings.size(), 1u);
    EXPECT_EQ(r.health_warnings.front().code, "HISTORY_TRUNCATED");
  }

  TEST(RunHealth, EveryCodeHasASwitchedManifestEntry) {
    for (int i = 0; i < 4; ++i) {
      const auto check = static_cast<RunHealthCheck>(i);
      EXPECT_STREQ(RunHealthCode(check), kRunHealthCodeList[i]);
      const auto& defs = WarningDefinitions();
      const auto it = std::find_if(defs.begin(), defs.end(),
                                   [&](const WarningDefinition& w) { return w.code == RunHealthCode(check); });
      ASSERT_NE(it, defs.end()) << RunHealthCode(check);
      EXPECT_TRUE(it->default_on);
    }
  }

} // namespace adsim
