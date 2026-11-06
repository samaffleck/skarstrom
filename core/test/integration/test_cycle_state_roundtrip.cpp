// test_cycle_state_roundtrip.cpp
//
// The cycle-state API (adsim/solver/cycle_state.h, Simulation::FinalState /
// SetInitialState / RunCycle):
//
//  * Round trip: N cycles, take the state, run one cycle from it on a fresh
//    Simulation of the same deck, and compare with N + 1 cycles run straight.
//    Checked on a deck with inter-bed connections (FiveStepProductPurge,
//    coarsened) and on a deck with a buffer tank whose consumer runs before its
//    source (so it reads the PREVIOUS cycle's summary). Dropping the histories
//    from the state must break the match on both, which shows they are needed.
//  * Warm start from a cyclic-steady-state state reaches CSS in a few cycles.
//  * A state that does not fit the deck is refused with a reason.
//
// Tolerance: the restarted cycle differs from the straight one only by the
// integrator's path (a fresh IDA start instead of the restart at the cycle
// boundary, a band probe at another state), so the end states agree to the
// integrator tolerance: max |dx| / max(RMS(field), abstol) over each field of
// the differential state is below kRoundTripTol (measured on Linux/GCC:
// 3e-12 buffer tank, 1.3e-5 Vysyaraju2024, 4e-7 FiveStepProductPurge; the
// bound leaves room for platform differences). Without the histories the
// first two miss by 0.23 and 0.25.

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/boundary_value.h"
#include "adsim/boundary_conditions/buffer_tank_source.h"
#include "adsim/boundary_conditions/flow_boundary_value.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/equations/isotherm.h"
#include "adsim/objects/user_data.h"
#include "adsim/simulation.h"
#include "adsim/solver/cycle_state.h"
#include "cycle_state_decks.h"

namespace adsim {

namespace {

using namespace test_decks;

constexpr double kRoundTripTol = 1e-4;
// A state without its histories must miss by far more than the tolerance.
constexpr double kMissingHistoryMin = 1e-2;

// Largest relative difference between two states over the differential
// entries of every field: max |a - b| / max(RMS(b), floor) per field.
double StateDifference(const CycleState& a, const CycleState& b) {
  EXPECT_EQ(a.y.size(), b.y.size());
  if (a.y.size() != b.y.size()) return 1e300;
  double worst = 0.0;
  for (const auto& f : b.fields) {
    if (!f.differential) continue;
    double ss = 0.0;
    size_t n = 0;
    for (size_t g : f.index) { ss += b.y[g] * b.y[g]; ++n; }
    const double scale = std::max(std::sqrt(ss / std::max<size_t>(n, 1)), 1e-8);
    for (size_t g : f.index) worst = std::max(worst, std::abs(a.y[g] - b.y[g]) / scale);
  }
  return worst;
}

// An example deck with inter-bed connections, on a coarse grid with the CSS
// check off so every run does its cycles.
Simulation ExampleSimulation(const std::string& deck) {
  Simulation sim = Simulation::FromFile("../examples/decks/" + deck + ".json");
  for (auto& layer : sim.Data().bed.layers_.layers_) layer.domain_.n_cells = 20;
  sim.Data().settings.css.method_ = CSSMethod::NONE;
  return sim;
}

// check_history: also check that the state without its histories misses (on
// decks where a connection or tank reads a step that runs later in the cycle,
// so the previous cycle's history matters).
template <class Make>
void ExpectRoundTrip(Make make, int n, const char* label, bool check_history = true) {
  Simulation straight = make();
  ASSERT_TRUE(straight.Run({.cycles = n + 1}).Ok()) << label << ": " << straight.LastResult().message;
  const CycleState reference = straight.FinalState();
  EXPECT_EQ(reference.cycles, n + 1);

  Simulation first = make();
  ASSERT_TRUE(first.Run({.cycles = n}).Ok()) << label << ": " << first.LastResult().message;
  const CycleState x = first.FinalState();
  EXPECT_EQ(x.cycles, n);
  EXPECT_TRUE(x.HasHistory()) << label;

  Simulation restarted = make();
  const RunResult r = restarted.RunCycle(x);
  ASSERT_TRUE(r.Ok()) << label << ": " << r.message;
  EXPECT_EQ(r.cycles_run, 1);
  EXPECT_EQ(restarted.Data().cycle.number_of_cycles_, make().Data().cycle.number_of_cycles_)
      << "RunCycle must not change the deck's number of cycles";
  const CycleState gx = restarted.FinalState();
  EXPECT_EQ(gx.cycles, n + 1);
  const double diff = StateDifference(gx, reference);
  std::cout << "[ cycle state ] " << label << ": round trip difference " << diff << "\n";
  EXPECT_LT(diff, kRoundTripTol) << label;

  if (!check_history) return;
  // Without the histories the restart replays nothing (connections) or the
  // first-cycle assumption (buffer tank): the cycle differs clearly.
  CycleState bare = x;
  for (auto& e : bare.endpoints) { e.history.Clear(); e.summary.Clear(); }
  Simulation no_history = make();
  ASSERT_TRUE(no_history.RunCycle(bare).Ok()) << label;
  const double miss = StateDifference(no_history.FinalState(), reference);
  std::cout << "[ cycle state ] " << label << ": without histories " << miss << "\n";
  EXPECT_GT(miss, kMissingHistoryMin) << label;
}

}  // namespace

TEST(CycleStateRoundTrip, BufferTankDeck) {
  ExpectRoundTrip(BufferTankSimulation, 2, "buffer tank");
}

TEST(CycleStateRoundTrip, ConnectedBedsDeck) {
  // Vysyaraju2024: a pressure-equalisation step replays a partner step that
  // runs later in the cycle, i.e. the previous cycle's logged outlet.
  ExpectRoundTrip([] { return ExampleSimulation("Vysyaraju2024_PSA_exp1"); }, 2, "Vysyaraju2024_PSA_exp1");
  // Four connections; each reads a step that ran earlier in the same cycle.
  ExpectRoundTrip([] { return ExampleSimulation("FiveStepProductPurge"); }, 2, "FiveStepProductPurge",
                  /*check_history=*/false);
}

TEST(CycleStateRoundTrip, RepeatedRunCycleStepsTheCycleMap) {
  // RunCycle() without a state continues from the last final state, so three
  // calls are three cycles of one run.
  Simulation straight = BufferTankSimulation();
  ASSERT_TRUE(straight.Run({.cycles = 3}).Ok());
  Simulation stepped = BufferTankSimulation();
  for (int k = 0; k < 3; ++k) ASSERT_TRUE(stepped.RunCycle().Ok()) << "cycle " << k;
  EXPECT_EQ(stepped.FinalState().cycles, 3);
  EXPECT_LT(StateDifference(stepped.FinalState(), straight.FinalState()), kRoundTripTol);
}

TEST(CycleStateRoundTrip, StateLayoutAndDifferentialVector) {
  Simulation sim = BufferTankSimulation();
  ASSERT_TRUE(sim.RunCycle().Ok());
  CycleState x = sim.FinalState();
  EXPECT_EQ(x.components, (std::vector<std::string>{"Oxygen", "Nitrogen"}));
  ASSERT_EQ(x.layers.size(), 1u);
  EXPECT_EQ(x.layers[0].cells, 10);
  ASSERT_EQ(x.steps.size(), 2u);
  EXPECT_EQ(x.steps[0].name, "Purge");
  EXPECT_EQ(x.y.size(), x.yp.size());
  size_t entries = 0;
  for (const auto& f : x.fields) entries += f.index.size();
  EXPECT_EQ(entries, x.y.size());

  const auto* u = x.FindField("Velocity");
  ASSERT_NE(u, nullptr);
  EXPECT_FALSE(u->differential);
  EXPECT_NE(x.FindField("Pressure"), nullptr);
  const auto y_o2 = x.FieldValues("Mole Fraction", "Oxygen");
  EXPECT_FALSE(y_o2.empty());
  EXPECT_NE(x.FindField("Adsorbed Concentration", "Nitrogen"), nullptr);
  EXPECT_THROW(x.FieldValues("No Such Field"), std::out_of_range);

  // Differential round trip, and a size check.
  std::vector<double> d = x.Differential();
  ASSERT_FALSE(d.empty());
  EXPECT_LT(d.size(), x.y.size());  // velocity and the ghost cells are algebraic
  CycleState copy = x;
  for (double& v : d) v *= 1.0;
  copy.SetDifferential(d);
  EXPECT_EQ(copy.y, x.y);
  EXPECT_THROW(copy.SetDifferential(std::vector<double>(d.size() + 1, 0.0)), std::invalid_argument);
}

TEST(CycleStateRoundTrip, IncompatibleStatesAreRefusedWithAReason) {
  Simulation sim = BufferTankSimulation();
  ASSERT_TRUE(sim.Run({.cycles = 2}).Ok());
  const CycleState x = sim.FinalState();
  EXPECT_EQ(sim.CheckState(x), "");

  {  // a different grid is not a refusal: the run remaps the state
    Simulation other = BufferTankSimulation();
    other.Data().bed.layers_.layers_[0].domain_.n_cells = 12;
    EXPECT_EQ(other.CheckState(x), "");
    EXPECT_NO_THROW(other.SetInitialState(x));
  }
  {  // a different component list
    Simulation other = BufferTankSimulation();
    other.Data().fluid.AddComponent("Argon");
    EXPECT_NE(other.CheckState(x).find("components"), std::string::npos);
  }
  {  // a changed step time is not a refusal: the histories are stretched
    Simulation other = BufferTankSimulation();
    other.Data().cycle.step_list_[1].step_time_ = 25.0;
    EXPECT_EQ(other.CheckState(x), "");
  }
  {  // a removed step
    Simulation other = BufferTankSimulation();
    other.Data().cycle.step_list_.pop_back();
    EXPECT_NE(other.CheckState(x).find("step structure"), std::string::npos);
  }
  {  // a different DAE layout (isothermal): found when the run starts
    Simulation other = BufferTankSimulation();
    other.Data().bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
    ASSERT_EQ(other.CheckState(x), "");
    const RunResult r = other.RunCycle(x);
    EXPECT_EQ(r.status, RunStatus::Failed);
    EXPECT_NE(r.message.find("state layout"), std::string::npos) << r.message;
    EXPECT_FALSE(other.HasFinalState());
  }
}

TEST(CycleStateRoundTrip, WarmStartFromCssReachesCssInAFewCycles) {
  // Buffer-tank deck with the Per-Component mass-balance CSS check (min 5
  // cycles, 3 consecutive passes): a cold run needs many cycles; a warm start
  // from its CSS state needs only the consecutive-cycles rule's minimum (the
  // first warm cycle has no previous cycle in the log to compare with).
  auto make = [] {
    Simulation sim = BufferTankSimulation();
    auto& css = sim.Data().settings.css;
    css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
    css.min_cycles_ = 5;
    css.consecutive_cycles_ = 3;
    sim.Data().cycle.number_of_cycles_ = 200;
    return sim;
  };
  Simulation cold = make();
  const RunResult rc = cold.Run();
  ASSERT_EQ(rc.status, RunStatus::ReachedCSS) << rc.message;

  Simulation warm = make();
  warm.SetInitialState(cold.FinalState());
  const RunResult rw = warm.Run();
  std::cout << "[ cycle state ] buffer tank: cold " << rc.cycles_run << " cycles to CSS, warm "
            << rw.cycles_run << "\n";
  ASSERT_EQ(rw.status, RunStatus::ReachedCSS) << rw.message;
  EXPECT_LE(rw.cycles_run, 4);
  EXPECT_LT(rw.cycles_run, rc.cycles_run);
  EXPECT_EQ(warm.FinalState().cycles, rc.cycles_run + rw.cycles_run);
  EXPECT_EQ(warm.Data().settings.css.min_cycles_, 5) << "the CSS floor must be restored after the run";
  EXPECT_TRUE(warm.HasInitialState());
  warm.ClearInitialState();
  EXPECT_FALSE(warm.HasInitialState());
}

}  // namespace adsim
