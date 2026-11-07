/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Anderson acceleration of the approach to cyclic steady state on two bundled
// decks, as shipped (their own grid, start and CSS settings).
//
// Measured on these decks (cycle of the CSS stop, 0-based):
//
//   deck          plain  Aitken  Anderson
//   Perez2019_A      47      32        35
//   Perez2019_E     206      65        43
//
// The reference KPIs are those of a long plain run with the CSS check off
// (cycle 499 for A, 999 for E). Over its last few hundred cycles the
// reference itself wanders by about 1e-7 (A) and 3e-6 (E) in purity and 2e-6
// (A) and 2e-5 (E) in recovery. The stops, references and errors were
// measured on the decks with their constant gas heat capacity (30.7 J/(mol K));
// re-measure them whenever the decks' inputs change. Each test runs Anderson
// only (a plain run of E takes 70 s) and checks that it stops well before the plain run, at KPIs no
// further from the reference than the plain stop's or within 5e-5 of it (half
// the 1e-4 to which CSSRobustDeck holds a plain run).
// Slow tier: about 30 s (A) and 20 s (E).

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include <gtest/gtest.h>

#include "adsim/simulation.h"

namespace adsim {

  namespace {
    struct DeckCase {
      const char* deck;
      int plain_stop;           // measured, plain run
      int max_anderson_stop;    // bound for the Anderson stop
      double ref_purity;        // CO2, long plain run
      double ref_recovery;
      double plain_purity_error;    // |plain stop - reference|, measured
      double plain_recovery_error;
    };

    constexpr double kKpiAllowance = 5e-5;

    void CheckAnderson(const DeckCase& c) {
      Simulation sim = Simulation::FromFile(std::string("../examples/decks/") + c.deck + ".json");
      ASSERT_NE(sim.Data().settings.css.method_, CSSMethod::NONE);  // GUARD: the deck stops at CSS
      sim.Data().settings.acc_css.method_ = AccelerationMethod::ANDERSON;
      sim.Data().cycle.number_of_cycles_ = 300;
      int applied = 0;
      RunOptions options;
      options.log_level = LogLevel::MEDIUM;
      options.on_message = [&applied](LogLevel, const std::string& m) {
        if (m.find("CSS acceleration applied") != std::string::npos) ++applied;
      };
      const RunResult run = sim.Run(options);
      ASSERT_TRUE(run.Ok()) << run.message;
      ASSERT_TRUE(run.css_reached) << c.deck << ": no CSS in 300 cycles with Anderson";

      const KPIRecord* co2 = nullptr;
      const auto kpis = sim.KPIs(-1);  // keep alive: co2 points into it
      for (const auto& k : kpis)
        if (k.component == "CarbonDioxide") co2 = &k;
      ASSERT_NE(co2, nullptr);
      const double purity_error = std::abs(co2->purity - c.ref_purity);
      const double recovery_error = std::abs(co2->recovery - c.ref_recovery);
      std::printf("  %s: Anderson stop at cycle %d (plain %d), %d jumps; purity error %.2e"
                  " (plain %.2e), recovery error %.2e (plain %.2e)\n",
                  c.deck, run.css_cycle, c.plain_stop, applied, purity_error, c.plain_purity_error,
                  recovery_error, c.plain_recovery_error);

      EXPECT_GT(applied, 0);
      EXPECT_LE(run.css_cycle, c.max_anderson_stop) << c.deck << ": Anderson no longer shortens the run";
      // At the same CSS.
      EXPECT_LE(purity_error, std::max(c.plain_purity_error, kKpiAllowance)) << c.deck;
      EXPECT_LE(recovery_error, std::max(c.plain_recovery_error, kKpiAllowance)) << c.deck;
    }
  }  // namespace

  TEST(CSSAndersonDecks, Perez2019AStopsSoonerAtTheSameKpis) {
    // Measured: stop at 35, purity error 6.6e-7, recovery error 2.5e-5.
    CheckAnderson({"Perez2019_A", 47, 40, 0.99508079, 0.59641507, 3.4e-7, 2.5e-5});
  }

  TEST(CSSAndersonDecks, Perez2019EStopsSoonerAtTheSameKpis) {
    // Measured: stop at 43, purity error 7.3e-6, recovery error 6.3e-6.
    CheckAnderson({"Perez2019_E", 206, 70, 0.86641524, 0.94362659, 6.9e-5, 1.9e-6});
  }

}  // namespace adsim
