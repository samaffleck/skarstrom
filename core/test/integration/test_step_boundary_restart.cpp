// test_step_boundary_restart.cpp
//
// Every cycle step starts IDA afresh (SolverData::RestartForNewStep): fresh
// history and hin = init_step at the new boundary conditions, instead of the
// previous step's step size.
//
// Before, steps ran back to back and the first step of a cycle step inherited
// the size of the previous step's last, tstop-clipped, step. With the
// boundary condition's algebraic jump out of the error test (suppress_alg_),
// a step whose forcing is a short pulse at its start could then be stepped
// over in one step. On Vysyaraju2024_PSA_exp1 the PER step replays the PED
// blowdown (2.8 mol/s decaying in ~0.1 s): when the BLO step happened to end
// on a long remainder (0.2 s instead of 0.015 s), PER's first step jumped the
// whole pulse, the bed missed its equalisation gas, and the cycle restarted
// its approach to CSS (PRES uptake 0.1497 -> 0.284 -> -0.002 -> 0.071 mol...).
// Which runs tripped hung on round-off: a reltol 1e-4 away from the deck's did
// natively, and the WebAssembly build did at the deck's own settings, so the
// web solver disagreed with the native one.

#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "adsim/simulation.h"
#include "adsim/boundary_conditions/step.h"

#include <gtest/gtest.h>

namespace adsim {

  namespace {
    const char* kDeck = "../examples/decks/Vysyaraju2024_PSA_exp1.json";
  }

  // The first internal step of every cycle step is at most init_step.
  TEST(StepBoundaryRestart, EveryStepStartsAtTheInitialStep) {
    Simulation sim = Simulation::FromFile(kDeck);
    const double h0 = sim.Data().settings.init_step.Get();
    ASSERT_GT(h0, 0.0);

    std::map<std::pair<int, int>, double> first_step;  // (cycle, step) -> first dt
    RunOptions run;
    run.cycles = 2;
    run.on_step = [&](double, double dt, double, int cycle, int step) {
      first_step.emplace(std::make_pair(cycle, step), std::abs(dt));
    };
    const RunResult result = sim.Run(run);
    ASSERT_TRUE(result.Ok()) << result.message;

    const int steps = static_cast<int>(sim.Data().cycle.step_list_.size());
    ASSERT_EQ(first_step.size(), static_cast<size_t>(2 * steps));
    for (const auto& [key, dt] : first_step) {
      EXPECT_LE(dt, h0 * (1.0 + 1e-12))
          << "cycle " << key.first << " step " << key.second << " started with a step of " << dt << " s";
    }
  }

  // With a reltol 1e-4 away from the deck's (which, before the restart, made
  // PER step over its pulse at cycle 47), the pressurisation uptake still
  // approaches CSS smoothly over the last cycles.
  TEST(StepBoundaryRestart, VysyarajuPressurisationApproachesCSSSmoothly) {
    Simulation sim = Simulation::FromFile(kDeck);
    UserData& data = sim.Data();
    data.settings.reltol = 1.0001e-6;

    std::vector<double> pres_moles;
    RunOptions run;
    run.progress = [&](const RunProgress& p) {
      if (p.step == 0) {
        const CycleEndpointSummary* s = data.cycle.step_list_[0].FindSummaryForEndpoint("");
        pres_moles.push_back(s ? s->total_moles : std::nan(""));
      }
      return true;
    };
    const RunResult result = sim.Run(run);
    ASSERT_TRUE(result.Ok()) << result.message;
    ASSERT_GE(pres_moles.size(), 40u);

    // Cycle to cycle the uptake creeps up by < 0.3 % from cycle 20 on (the
    // step over the pulse moved it by 90 % and more).
    for (size_t c = 20; c < pres_moles.size(); ++c) {
      EXPECT_NEAR(pres_moles[c], pres_moles[c - 1], 3e-3 * std::abs(pres_moles[c - 1]))
          << "PRES uptake jumped at cycle " << c;
    }
    // 0.1681 since connected streams keep the sign of the donor's flow (with
    // the magnitude replayed, PER took in gas PED never released: 0.1498).
    EXPECT_NEAR(pres_moles.back(), 0.1681, 0.0005);
  }

} // namespace adsim
