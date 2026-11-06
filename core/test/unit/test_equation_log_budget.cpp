/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// The equation log's memory budget per (cycle, step): LogState records a full
// snapshot at every accepted solver step, so a step that takes very many
// solver steps (a stalled run) used to grow without bound (skarstrom#74). With
// a budget, the step's rows stay bounded by the budget, not by the number of
// solver steps; the first and the newest row are kept exactly; the spacing
// stays roughly even; a step under the budget is logged exactly as before; and
// the solver is told once per step.

#include <algorithm>
#include <string>
#include <vector>

#include "../fixtures/user_data_test.h"
#include "adsim/solver/log.h"

namespace adsim {

  class EquationLogBudget : public UserDataTest {
  protected:
    void SetUp() override {
      AddComponents({"CarbonDioxide", "Nitrogen"});
      auto [layer_id, ads_id] = AddLayer("Layer1");
      GetLayer(layer_id).domain_.n_cells = 10;
      auto henry = std::make_unique<HenryIsotherm>();
      for (const auto& c : {"CarbonDioxide", "Nitrogen"}) {
        henry->params_[c].k = 0.0;
        henry->params_[c].Q = 0.0;
      }
      AddIsotherm({"CarbonDioxide", "Nitrogen"}, ads_id, std::move(henry));
      for (const auto& c : {"CarbonDioxide", "Nitrogen"}) {
        AddLinearDrivingForceModel(ads_id, c, 1.0);
        SetHOA(ads_id, c, 0.0);
      }
      AddConstantInternalHeatTransferCoefficientModel(10.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(10.0, layer_id);
      std::get<CylinderBed>(data.bed.shape_).diameter = 0.05;
      AddInitialConditions(1.0e5, 298.15, 298.15, {{"CarbonDioxide", 0.0}, {"Nitrogen", 1.0}});
      AddStep("Feed", 100.0);
      data.Initialise();
      log_.Initialise(data);
      log_.message_sink_ = [this](const std::string& m) { messages_.push_back(m); };
    }

    // Log `n` rows of the same state at times 0, dt, 2 dt, ... in (cycle, step).
    void LogRows(int n, int cycle = 0, int step = 0, double dt = 1e-3) {
      for (int i = 0; i < n; ++i) log_.LogState(data.state, i * dt, step, cycle);
    }

    const Log& StepLog(int cycle = 0, int step = 0) const { return log_.GetLog(cycle, step); }

    EquationLog log_;
    std::vector<std::string> messages_;
  };

  // Under the budget (the default 2 GiB here) every row is kept, as before.
  TEST_F(EquationLogBudget, UnderTheBudgetEveryRowIsKept) {
    LogRows(1500);
    const Log& l = StepLog();
    EXPECT_EQ(l.step_data.times.size(), 1500u);
    EXPECT_EQ(l.step_data.thinnings, 0);
    EXPECT_GT(l.step_data.row_bytes, 0u);
    for (const auto& [type, data] : l.variables) EXPECT_EQ(data.values.size(), 1500u);
    EXPECT_TRUE(messages_.empty());
  }

  // Over the budget the rows stay between half the budget and the budget,
  // however many solver steps are logged; the first and the newest row are
  // exact, the times stay increasing and roughly even, and every variable
  // keeps one value per time.
  TEST_F(EquationLogBudget, OverTheBudgetRowsStayBounded) {
    LogRows(1);
    const std::size_t row_bytes = StepLog().step_data.row_bytes;
    ASSERT_GT(row_bytes, 0u);
    const StepHeatDuty duty = StepLog().heat_duty;  // per-step value: thinning must not touch it
    const std::size_t max_rows = 400;
    log_.max_bytes_per_step_ = max_rows * row_bytes;

    const int n = 20000;  // 50x the budget
    LogRows(n);           // a new (cycle, step) entry was not opened: rows continue from the first one
    const Log& l = StepLog();
    const auto& t = l.step_data.times;
    EXPECT_LE(t.size(), max_rows);
    EXPECT_GE(t.size(), max_rows / 2);
    EXPECT_GT(l.step_data.thinnings, 0);
    EXPECT_DOUBLE_EQ(t.front(), 0.0);
    EXPECT_DOUBLE_EQ(t.back(), (n - 1) * 1e-3);  // the newest row: the end-of-step state
    EXPECT_TRUE(std::is_sorted(t.begin(), t.end()));
    double gap_min = 1e300, gap_max = 0.0;
    for (std::size_t i = 1; i + 1 < t.size(); ++i) {  // interior gaps (the last one can be short)
      gap_min = std::min(gap_min, t[i] - t[i - 1]);
      gap_max = std::max(gap_max, t[i] - t[i - 1]);
    }
    EXPECT_LE(gap_max, 2.0 * gap_min + 1e-12) << "spacing should stay roughly even";
    for (const auto& [type, data] : l.variables) EXPECT_EQ(data.values.size(), t.size());
    for (const auto& [component, vars] : l.component_variables)
      for (const auto& [type, data] : vars) EXPECT_EQ(data.values.size(), t.size()) << component;

    EXPECT_EQ(l.heat_duty.heating_power, duty.heating_power);
    EXPECT_EQ(l.heat_duty.cooling_power, duty.cooling_power);

    // One message for the step, however many thinnings it took.
    ASSERT_EQ(messages_.size(), 1u);
    EXPECT_NE(messages_[0].find("Feed"), std::string::npos) << messages_[0];
    EXPECT_NE(messages_[0].find("memory budget"), std::string::npos) << messages_[0];
  }

  // The budget is per (cycle, step): a new step starts empty, and its own first
  // overrun reports again.
  TEST_F(EquationLogBudget, TheBudgetIsPerStep) {
    LogRows(1);
    log_.max_bytes_per_step_ = 100 * StepLog().step_data.row_bytes;
    LogRows(1000, 0);
    LogRows(50, 1);  // next cycle, under the budget
    EXPECT_EQ(StepLog(1).step_data.times.size(), 50u);
    EXPECT_EQ(StepLog(1).step_data.thinnings, 0);
    LogRows(1000, 1);
    EXPECT_LE(StepLog(1).step_data.times.size(), 100u);
    EXPECT_EQ(messages_.size(), 2u);
  }

  // 0 disables the budget.
  TEST_F(EquationLogBudget, ZeroDisablesTheBudget) {
    log_.max_bytes_per_step_ = 0;
    LogRows(3000);
    EXPECT_EQ(StepLog().step_data.times.size(), 3000u);
    EXPECT_TRUE(messages_.empty());
  }

}  // namespace adsim
