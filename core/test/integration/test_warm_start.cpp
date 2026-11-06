// test_warm_start.cpp
//
// Warm starts from a previous run (RunOptions::initial_state,
// Simulation::SetInitialState) after the edits a study makes between runs:
//
//  * the endpoint histories are stretched onto changed step times
//    (ResampleCycleStateHistories), and a zero-length step clears them;
//  * a warm start after a step-time, feed-composition or cell-count edit
//    reaches the same cyclic steady state as a cold start of the edited deck
//    (purity and recovery within kKpiTol) in fewer cycles;
//  * RunOptions::initial_state applies to one run only, and the run says how
//    it adapted the state.

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/simulation.h"
#include "adsim/solver/cycle_state.h"
#include "cycle_state_decks.h"

namespace adsim {

namespace {

using namespace test_decks;

// Purity and recovery of a warm and a cold start agree to the error the CSS
// stopping rule leaves (Per-Component balance 0.5 %, state change 2e-5).
constexpr double kKpiTol = 2e-3;

const CycleStateEndpoint& EndpointOfStep(const CycleState& s, size_t step_index) {
  for (const auto& e : s.endpoints)
    if (e.step_id == s.steps[step_index].id && !e.history.time_series_.empty()) return e;
  throw std::runtime_error("no endpoint with history on that step");
}

}  // namespace

TEST(WarmStart, HistoriesAreStretchedOntoChangedStepTimes) {
  Simulation sim = BufferTankSimulation();
  ASSERT_TRUE(sim.Run({.cycles = 2}).Ok());
  const CycleState x = sim.FinalState();

  // Feed (step 2) from 20 s to 25 s.
  Simulation edited = BufferTankSimulation();
  edited.Data().cycle.step_list_[1].step_time_ = 25.0;
  EXPECT_EQ(edited.CheckState(x), "");
  std::string note;
  const CycleState y = edited.PrepareState(x, &note);
  EXPECT_NE(note.find("\"Feed\" 20 s -> 25 s"), std::string::npos) << note;
  EXPECT_DOUBLE_EQ(y.steps[1].time, 25.0);
  EXPECT_EQ(y.y, x.y);  // same grid: the DAE vector is untouched

  const auto& before = EndpointOfStep(x, 1);
  const auto& after = EndpointOfStep(y, 1);
  ASSERT_EQ(after.history.time_series_.size(), before.history.time_series_.size());
  EXPECT_NEAR(after.history.time_series_.back(), 1.25 * before.history.time_series_.back(), 1e-12);
  EXPECT_EQ(after.history.molar_flow_, before.history.molar_flow_);  // same shape
  ASSERT_TRUE(before.summary.is_valid);
  EXPECT_NEAR(after.summary.total_moles, 1.25 * before.summary.total_moles, 1e-12 * std::abs(before.summary.total_moles));
  // Purge (step 1) kept its time and its history.
  EXPECT_EQ(EndpointOfStep(y, 0).history.time_series_, EndpointOfStep(x, 0).history.time_series_);

  // A zero-length step has nothing to stretch: its histories are cleared.
  CycleState zero = x;
  zero.steps[1].time = 0.0;
  const CycleState z = ResampleCycleStateHistories(zero, edited.Data());
  for (const auto& e : z.endpoints)
    if (e.step_id == z.steps[1].id) {
      EXPECT_TRUE(e.history.time_series_.empty());
      EXPECT_FALSE(e.summary.is_valid);
    }

  // ... and the run from the stretched state starts and says so.
  const RunResult r = edited.Run({.cycles = 1, .initial_state = &x});
  ASSERT_TRUE(r.Ok()) << r.message;
  bool said = false;
  for (const auto& m : edited.Messages()) said = said || m.text.find("Initial state: histories stretched") == 0;
  EXPECT_TRUE(said);
}

TEST(WarmStart, InitialStateOptionAppliesToOneRun) {
  Simulation sim = BufferTankSimulation();
  ASSERT_TRUE(sim.Run({.cycles = 3}).Ok());
  const CycleState x = sim.FinalState();

  Simulation other = BufferTankSimulation();
  ASSERT_TRUE(other.Run({.cycles = 1, .initial_state = &x}).Ok());
  EXPECT_EQ(other.FinalState().cycles, 4);
  EXPECT_FALSE(other.HasInitialState());
  ASSERT_TRUE(other.Run({.cycles = 1}).Ok());
  EXPECT_EQ(other.FinalState().cycles, 1);  // the next run starts cold again

  // The final state itself can be passed back in (it is replaced by the run).
  ASSERT_TRUE(other.Run({.cycles = 1, .initial_state = &other.FinalState()}).Ok());
  EXPECT_EQ(other.FinalState().cycles, 2);

  // A state that does not fit fails the run with the reason.
  Simulation wrong = BufferTankSimulation();
  wrong.Data().fluid.AddComponent("Argon");
  const RunResult r = wrong.Run({.cycles = 1, .initial_state = &x});
  EXPECT_EQ(r.status, RunStatus::Failed);
  EXPECT_NE(r.message.find("components"), std::string::npos) << r.message;
}

namespace {

struct Edit {
  const char* label;
  void (*apply)(Simulation&);
};

// Cold and warm (from the base deck's CSS) runs of the edited deck: same
// purity and recovery, fewer cycles warm.
template <class Make>
void ExpectWarmMatchesCold(Make make, const std::vector<Edit>& edits, const char* deck) {
  Simulation base = make();
  ASSERT_EQ(base.Run().status, RunStatus::ReachedCSS) << deck << ": " << base.LastResult().message;
  const CycleState x = base.FinalState();
  for (const auto& edit : edits) {
    Simulation cold = make();
    edit.apply(cold);
    const RunResult rc = cold.Run();
    ASSERT_EQ(rc.status, RunStatus::ReachedCSS) << deck << " / " << edit.label << ": " << rc.message;
    Simulation warm = make();
    edit.apply(warm);
    const RunResult rw = warm.Run({.initial_state = &x});
    ASSERT_EQ(rw.status, RunStatus::ReachedCSS) << deck << " / " << edit.label << ": " << rw.message;
    const auto kc = cold.KPIs(-1), kw = warm.KPIs(-1);
    ASSERT_EQ(kc.size(), kw.size());
    ASSERT_FALSE(kc.empty());
    for (size_t i = 0; i < kc.size(); ++i) {
      EXPECT_NEAR(kw[i].purity, kc[i].purity, kKpiTol) << deck << " / " << edit.label << " " << kc[i].component;
      EXPECT_NEAR(kw[i].recovery, kc[i].recovery, kKpiTol) << deck << " / " << edit.label << " " << kc[i].component;
    }
    std::cout << "[ warm start ] " << deck << " / " << edit.label << ": cold " << rc.cycles_run << " cycles, warm "
              << rw.cycles_run << "; purity " << kc[0].purity << " / " << kw[0].purity << ", recovery "
              << kc[0].recovery << " / " << kw[0].recovery << "\n";
    EXPECT_LT(rw.cycles_run, rc.cycles_run) << deck << " / " << edit.label;
  }
}

}  // namespace

TEST(WarmStart, ReachesTheColdCssAfterTypicalEdits) {
  // The TSA example on a coarse grid (CSS in about 15 cycles).
  auto make = [] {
    Simulation sim = Simulation::FromFile("../examples/decks/Illustrative_TSA_13X.json");
    for (auto& layer : sim.Data().bed.layers_.layers_) layer.domain_.n_cells = 6;
    sim.Data().cycle.number_of_cycles_ = 200;
    return sim;
  };
  ExpectWarmMatchesCold(make,
                        {{"adsorption step +10 %",
                          [](Simulation& s) { s.SetParameter("step[1].time", 1.1 * s.GetParameter("step[1].time")); }},
                         {"feed CO2 0.15 -> 0.20",
                          [](Simulation& s) {
                            s.SetParameter("endpoint[1].y[CarbonDioxide]", 0.20);
                            s.SetParameter("endpoint[1].y[Nitrogen]", 0.80);
                          }},
                         {"cells 6 -> 9",
                          [](Simulation& s) {
                            for (auto& layer : s.Data().bed.layers_.layers_) layer.domain_.n_cells = 9;
                          }}},
                        "Illustrative_TSA_13X");
}

TEST(WarmStart, ReachesTheColdCssAfterAStepTimeEditOnABufferTankDeck) {
  // The buffer tank reads the previous cycle's summary of a step whose time
  // changes, so the warm start goes through the stretched histories.
  auto make = [] {
    Simulation sim = BufferTankSimulation();
    auto& css = sim.Data().settings.css;
    css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
    css.min_cycles_ = 5;
    css.consecutive_cycles_ = 3;
    sim.Data().cycle.number_of_cycles_ = 200;
    return sim;
  };
  // This deck tracks no KPIs; compare the end states instead.
  Simulation base = make();
  ASSERT_EQ(base.Run().status, RunStatus::ReachedCSS);
  const CycleState x = base.FinalState();
  auto edit = [](Simulation& s) { s.Data().cycle.step_list_[1].step_time_ = 24.0; };
  Simulation cold = make();
  edit(cold);
  const RunResult rc = cold.Run();
  ASSERT_EQ(rc.status, RunStatus::ReachedCSS) << rc.message;
  Simulation warm = make();
  edit(warm);
  const RunResult rw = warm.Run({.initial_state = &x});
  ASSERT_EQ(rw.status, RunStatus::ReachedCSS) << rw.message;
  const CycleState a = warm.FinalState(), b = cold.FinalState();
  ASSERT_EQ(a.y.size(), b.y.size());
  double worst = 0.0;
  for (const auto& f : b.fields) {
    if (!f.differential) continue;
    double ss = 0.0;
    for (size_t g : f.index) ss += b.y[g] * b.y[g];
    const double scale = std::max(std::sqrt(ss / f.index.size()), 1e-8);
    for (size_t g : f.index) worst = std::max(worst, std::abs(a.y[g] - b.y[g]) / scale);
  }
  std::cout << "[ warm start ] buffer tank / Feed 20 s -> 24 s: cold " << rc.cycles_run << " cycles, warm "
            << rw.cycles_run << "; end-state difference " << worst << "\n";
  EXPECT_LT(worst, 1e-2);
  EXPECT_LT(rw.cycles_run, rc.cycles_run);
}

}  // namespace adsim
