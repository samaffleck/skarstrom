// test_chunked_warm_start.cpp
//
// The chunked run (Simulation::Start / Advance) and the cycle-state API, as
// for blocking runs:
//
//  * Without a start state the chunked run is unchanged: it ends in the same
//    state as a blocking cold run, bit for bit, and its final state is kept.
//  * Start from a state at cyclic steady state (the argument, or
//    SetInitialState) reaches the same CSS as the cold chunked run in fewer
//    cycles; the state's cycles count toward the CSS minimum number of
//    cycles, and the deck's minimum is unchanged between chunks.
//  * A state that does not fit the deck makes Start return false with the
//    reason; a stopped chunked run has no final state.
//  * A state on another grid is adapted first, and the run says so.
//
// Tolerance: a warm and a cold run stop at the same CSS within the bound the
// CSS stopping rule leaves (Per-Component balance 0.5 %): the end states
// differ by kCssStateTol relative to each field's RMS (measured on Linux/GCC:
// 2.4e-4 on the buffer-tank deck).

#include <algorithm>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/simulation.h"
#include "adsim/solver/cycle_state.h"
#include "cycle_state_decks.h"

namespace adsim {

namespace {

using namespace test_decks;

constexpr double kCssStateTol = 2e-3;

// Largest relative difference between two states over the differential
// entries of every field: max |a - b| / max(RMS(b), floor) per field.
double StateDifference(const CycleState& a, const CycleState& b) {
  EXPECT_EQ(a.y.size(), b.y.size());
  if (a.y.size() != b.y.size()) return 1e300;
  double worst = 0.0;
  for (const auto& f : b.fields) {
    if (!f.differential) continue;
    double ss = 0.0;
    for (size_t g : f.index) ss += b.y[g] * b.y[g];
    const double scale = std::max(std::sqrt(ss / std::max<size_t>(f.index.size(), 1)), 1e-8);
    for (size_t g : f.index) worst = std::max(worst, std::abs(a.y[g] - b.y[g]) / scale);
  }
  return worst;
}

Simulation CssDeck() {
  Simulation sim = BufferTankSimulation();
  auto& css = sim.Data().settings.css;
  css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
  css.min_cycles_ = 5;
  css.consecutive_cycles_ = 3;
  sim.Data().cycle.number_of_cycles_ = 200;
  return sim;
}

// Drives a started chunked run to its end in short chunks. `between` runs
// after every chunk while the run is still going.
RunStatus Drive(Simulation& sim, const std::function<void()>& between = {}) {
  RunStatus st = RunStatus::Running;
  while (st == RunStatus::Running) {
    st = sim.Advance(50);
    if (st == RunStatus::Running && between) between();
  }
  return st;
}

}  // namespace

TEST(ChunkedWarmStart, ColdChunkedRunIsUnchangedAndKeepsItsFinalState) {
  Simulation blocking = BufferTankSimulation();
  ASSERT_TRUE(blocking.Run({.cycles = 3}).Ok());
  const CycleState reference = blocking.FinalState();

  Simulation chunked = BufferTankSimulation();
  chunked.Data().cycle.number_of_cycles_ = 3;
  ASSERT_TRUE(chunked.Start(LogLevel::LOW));
  EXPECT_FALSE(chunked.HasFinalState());  // not before the run ends
  ASSERT_EQ(Drive(chunked), RunStatus::Completed) << chunked.LastResult().message;
  ASSERT_TRUE(chunked.HasFinalState());
  const CycleState& x = chunked.FinalState();
  EXPECT_EQ(x.cycles, 3);
  EXPECT_EQ(x.y, reference.y);    // same integrator path: bit for bit
  EXPECT_EQ(x.yp, reference.yp);
  EXPECT_EQ(chunked.Messages().end(),
            std::find_if(chunked.Messages().begin(), chunked.Messages().end(),
                         [](const SolverMessage& m) { return m.text.rfind("Initial state:", 0) == 0; }));
}

TEST(ChunkedWarmStart, WarmChunkedRunReachesTheSameCssInFewerCycles) {
  Simulation cold = CssDeck();
  ASSERT_TRUE(cold.Start(LogLevel::LOW));
  ASSERT_EQ(Drive(cold), RunStatus::ReachedCSS) << cold.LastResult().message;
  const int cold_cycles = cold.LastResult().cycles_run;
  ASSERT_TRUE(cold.HasFinalState());
  const CycleState x = cold.FinalState();
  EXPECT_EQ(x.cycles, cold_cycles);

  // The state as the argument of Start.
  Simulation warm = CssDeck();
  ASSERT_TRUE(warm.Start(LogLevel::LOW, {}, &x));
  // The deck's CSS minimum is not changed between chunks (only for the check).
  ASSERT_EQ(Drive(warm, [&] { EXPECT_EQ(warm.Data().settings.css.min_cycles_, 5); }), RunStatus::ReachedCSS)
      << warm.LastResult().message;
  const int warm_cycles = warm.LastResult().cycles_run;
  std::cout << "[ chunked warm start ] buffer tank: cold " << cold_cycles << " cycles to CSS, warm "
            << warm_cycles << ", state difference " << StateDifference(warm.FinalState(), x) << "\n";
  EXPECT_LE(warm_cycles, 4);
  EXPECT_LT(warm_cycles, cold_cycles);
  EXPECT_LT(StateDifference(warm.FinalState(), x), kCssStateTol);
  EXPECT_EQ(warm.FinalState().cycles, cold_cycles + warm_cycles);
  EXPECT_EQ(warm.Data().settings.css.min_cycles_, 5);
  EXPECT_FALSE(warm.HasInitialState());  // the argument applies to one run

  // The state set with SetInitialState: the same run.
  Simulation set = CssDeck();
  set.SetInitialState(x);
  ASSERT_TRUE(set.Start(LogLevel::LOW));
  ASSERT_EQ(Drive(set), RunStatus::ReachedCSS) << set.LastResult().message;
  EXPECT_EQ(set.LastResult().cycles_run, warm_cycles);
  EXPECT_EQ(set.FinalState().y, warm.FinalState().y);

  // The final state itself can be the start (Start copies it first).
  ASSERT_TRUE(warm.Start(LogLevel::LOW, {}, &warm.FinalState()));
  ASSERT_EQ(Drive(warm), RunStatus::ReachedCSS) << warm.LastResult().message;
  EXPECT_EQ(warm.FinalState().cycles, cold_cycles + warm_cycles + warm.LastResult().cycles_run);
}

TEST(ChunkedWarmStart, StateThatDoesNotFitIsRefused) {
  Simulation sim = BufferTankSimulation();
  ASSERT_TRUE(sim.Run({.cycles = 1}).Ok());
  const CycleState x = sim.FinalState();

  Simulation wrong = BufferTankSimulation();
  wrong.Data().fluid.AddComponent("Argon");
  std::vector<std::string> said;
  EXPECT_FALSE(wrong.Start(LogLevel::LOW, [&](LogLevel, const std::string& m) { said.push_back(m); }, &x));
  EXPECT_EQ(wrong.LastResult().status, RunStatus::Failed);
  EXPECT_NE(wrong.LastResult().message.find("components"), std::string::npos) << wrong.LastResult().message;
  EXPECT_NE(wrong.LastResult().message.find(wrong.CheckState(x)), std::string::npos);
  EXPECT_FALSE(said.empty());
  EXPECT_FALSE(wrong.IsRunning());
  EXPECT_FALSE(wrong.HasFinalState());
  EXPECT_EQ(wrong.Advance(10), RunStatus::Failed);
}

TEST(ChunkedWarmStart, StoppedRunHasNoFinalState) {
  Simulation sim = BufferTankSimulation();
  sim.Data().cycle.number_of_cycles_ = 2;
  ASSERT_TRUE(sim.Start(LogLevel::LOW));
  ASSERT_EQ(Drive(sim), RunStatus::Completed);
  ASSERT_TRUE(sim.HasFinalState());

  // The next run clears it; stopping that run leaves none.
  ASSERT_TRUE(sim.Start(LogLevel::LOW, {}, &sim.FinalState()));
  EXPECT_FALSE(sim.HasFinalState());
  sim.Advance(1);
  sim.Stop();
  EXPECT_EQ(sim.LastResult().status, RunStatus::Cancelled);
  EXPECT_FALSE(sim.HasFinalState());
  EXPECT_EQ(sim.Advance(10), RunStatus::Cancelled);
  EXPECT_FALSE(sim.HasFinalState());
}

TEST(ChunkedWarmStart, StateOnAnotherGridIsAdaptedAndTheRunSaysSo) {
  Simulation coarse = TwoLayerSimulation(6, 4);
  ASSERT_TRUE(coarse.Run({.cycles = 1}).Ok());
  const CycleState x = coarse.FinalState();

  Simulation fine = TwoLayerSimulation(12, 8);
  fine.Data().cycle.number_of_cycles_ = 1;
  std::vector<std::pair<LogLevel, std::string>> said;
  ASSERT_TRUE(fine.Start(LogLevel::LOW, [&](LogLevel l, const std::string& m) { said.push_back({l, m}); }, &x))
      << fine.LastResult().message;
  bool noted = false;
  for (const auto& [level, text] : said)
    noted = noted || (level == LogLevel::LOW && text.rfind("Initial state:", 0) == 0);
  EXPECT_TRUE(noted);
  ASSERT_TRUE(RunStatus::Completed == Drive(fine)) << fine.LastResult().message;
  EXPECT_EQ(fine.FinalState().cycles, x.cycles + 1);
}

}  // namespace adsim
