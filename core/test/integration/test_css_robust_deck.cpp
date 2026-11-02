// End-to-end check of the cyclic-steady-state stopping rule on a bundled deck.
//
// Perez2019_A (examples/decks) uses Mass Balance (Per-Component) with the
// default state tolerance and consecutive cycles. Under the old balance-only,
// single-cycle rule it stopped at cycle 43 (0-based) with the CO2 recovery
// 1.3e-3 short of its converged value; with the robust rule it must stop close
// to the 1000-cycle values below (a run with the CSS check off, measured for
// this test). Loading pre-rule decks is covered in test_css_stopping_rule.cpp.

#include <cmath>
#include <string>

#include <gtest/gtest.h>

#include "adsim/simulation.h"

namespace adsim {

  namespace {
    constexpr double kConvergedPurity = 0.9948831;    // CO2, cycle 1000
    constexpr double kConvergedRecovery = 0.5911967;  // CO2, cycle 1000

    const KPIRecord* CarbonDioxide(const std::vector<KPIRecord>& kpis) {
      for (const auto& k : kpis)
        if (k.component == "CarbonDioxide") return &k;
      return nullptr;
    }
  }  // namespace

  TEST(CSSRobustDeck, Perez2019AStopsNearTheConvergedKPIs) {
    Simulation sim = Simulation::FromFile("../examples/decks/Perez2019_A.json");
    auto& css = sim.Data().settings.css;
    ASSERT_EQ(css.method_, CSSMethod::MASS_BALANCE_COMPONENT);
    EXPECT_DOUBLE_EQ(css.state_tolerance_.Get(), CyclicSteadyState::kDefaultStateTolerance);
    EXPECT_EQ(css.consecutive_cycles_, CyclicSteadyState::kDefaultConsecutiveCycles);

    sim.Data().cycle.number_of_cycles_ = 300;
    RunResult run = sim.Run();
    ASSERT_TRUE(run.Ok()) << run.message;
    ASSERT_TRUE(run.css_reached);
    EXPECT_GT(run.css_cycle, 50);

    // On failure, say which cycles are logged and what each one's KPIs are.
    std::string diag = "css_cycle " + std::to_string(run.css_cycle) + ", last logged " +
                       std::to_string(sim.LastCycle()) + "; KPI config:";
    for (const auto& cfg : sim.Data().kpi_config.per_component)
      diag += " [" + cfg.product_component + " tracked " + std::to_string(cfg.tracked) + " feed " +
              std::to_string(cfg.feed_endpoints.size()) + " product " +
              std::to_string(cfg.product_endpoints.size()) + "]";
    diag += "; log molar masses:";
    for (const auto& [name, m] : sim.Log().sys_data.molar_mass) diag += " " + name + "=" + std::to_string(m);
    diag += "; logged cycles:";
    for (int c : sim.Cycles()) {
      const KPIRecord* k = CarbonDioxide(sim.KPIs(c));
      diag += " " + std::to_string(c);
      if (k) diag += "(valid " + std::to_string(k->valid) + ", purity " + std::to_string(k->purity) +
                     ", recovery " + std::to_string(k->recovery) + ")";
    }

    const KPIRecord* co2 = CarbonDioxide(sim.KPIs(-1));
    ASSERT_NE(co2, nullptr) << diag;
    EXPECT_NEAR(co2->purity, kConvergedPurity, 1e-4) << diag;
    EXPECT_NEAR(co2->recovery, kConvergedRecovery, 1e-4) << diag;

    // The last record is the stopping cycle: the balance and the state change
    // within their tolerances on every cycle of the window, which the log holds.
    const auto history = sim.CSSHistory();
    ASSERT_GE(static_cast<int>(history.size()), css.CyclesToLog());
    const CSSRecord& last = history.back();
    EXPECT_EQ(last.cycle, run.css_cycle);
    EXPECT_TRUE(last.converged);
    EXPECT_EQ(last.passing_cycles, CyclicSteadyState::kDefaultConsecutiveCycles);
    EXPECT_LE(last.error, last.tolerance);
    EXPECT_LE(last.state_change, last.state_tolerance);
    for (size_t i = 0; i + 1 < history.size(); ++i) EXPECT_FALSE(history[i].converged) << i;
  }

  TEST(CSSRobustDeck, Perez2019AWithTheOldRuleStopsEarly) {
    // state_tolerance 0 and one cycle restore the balance-only test, which
    // passes long before the KPIs settle (kept to document why it changed).
    Simulation sim = Simulation::FromFile("../examples/decks/Perez2019_A.json");
    auto& css = sim.Data().settings.css;
    css.state_tolerance_ = 0.0;
    css.consecutive_cycles_ = 1;
    sim.Data().settings.max_cycles_to_keep_ = 1;  // no floor for the balance-only test
    sim.Data().cycle.number_of_cycles_ = 300;
    RunResult run = sim.Run();
    ASSERT_TRUE(run.Ok()) << run.message;
    ASSERT_TRUE(run.css_reached);
    EXPECT_LT(run.css_cycle, 50);
    const KPIRecord* co2 = CarbonDioxide(sim.KPIs(-1));
    ASSERT_NE(co2, nullptr);
    EXPECT_GT(std::abs(co2->recovery - kConvergedRecovery), 5e-4);
    EXPECT_EQ(sim.Cycles().size(), 2u);  // the last cycle and one earlier
  }

}  // namespace adsim
