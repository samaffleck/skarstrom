// test_warm_start_cache.cpp
//
// Warm-started optimisation and studies (adsim/solver/state_cache.h):
//
//  * RunWithStateCache starts from the nearest cached state, caches the final
//    state of a finished run, and skips (with a note) a state that does not
//    fit the deck;
//  * OptimizeSimulation warm-starts every candidate after the first, re-runs
//    the optimum cold and reports the cold value; with warm_start off no run
//    is warm;
//  * a warm-started sensitivity sample reports the same KPIs as a cold one
//    within the CSS stopping rule's error.
//
// Coarse grids (fast tier). The measured savings on full decks
// are in docs/cyclic-steady-state.md (Warm-started optimisation and studies).

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/simulation.h"
#include "adsim/solver/optimisation.h"
#include "adsim/solver/sensitivity.h"
#include "adsim/solver/state_cache.h"
#include "cycle_state_decks.h"

namespace adsim {

namespace {

using namespace test_decks;

std::string ReadDeck(const std::string& name) {
  std::ifstream f("../examples/decks/" + name);
  std::stringstream s;
  s << f.rdbuf();
  return s.str();
}

}  // namespace

TEST(WarmStartCache, RunsStartFromTheNearestCachedState) {
  // The TSA example on 6 cells, run to CSS: the Adsorption step time (step 1)
  // is the decision variable.
  const std::string deck = ReadDeck("Illustrative_TSA_13X.json");
  ASSERT_FALSE(deck.empty());
  auto Tsa = [&](double t_ads) {
    Simulation sim = Simulation::FromJson(deck);
    sim.SetParameter("layer[1].cells", 6);
    sim.SetParameter("step[1].time", t_ads);
    return sim;
  };
  CycleStateCache cache({500.0}, {700.0});

  Simulation a = Tsa(600.0);
  CachedRunInfo info;
  RunResult cold = RunWithStateCache(a, &cache, {600.0}, {}, &info);
  ASSERT_TRUE(cold.css_reached) << cold.message;
  EXPECT_FALSE(info.warm);  // empty cache: cold
  EXPECT_TRUE(info.cached);
  EXPECT_EQ(cache.Size(), 1u);

  Simulation b = Tsa(640.0);
  RunResult r = RunWithStateCache(b, &cache, {640.0}, {}, &info);
  ASSERT_TRUE(r.css_reached) << r.message;
  EXPECT_TRUE(info.warm);
  EXPECT_NEAR(info.distance, 0.2, 1e-12);
  EXPECT_EQ(info.start_cycles, cold.cycles_run);
  EXPECT_TRUE(info.note.empty()) << info.note;
  EXPECT_LT(r.cycles_run, cold.cycles_run);
  ASSERT_TRUE(b.HasFinalState());
  EXPECT_EQ(b.FinalState().cycles, cold.cycles_run + r.cycles_run);  // the cycles behind the start count
  EXPECT_EQ(cache.Size(), 2u);

  // A run that stops at the cycle limit is not cached: its end state depends
  // on its start.
  Simulation c = Tsa(560.0);
  r = RunWithStateCache(c, &cache, {560.0}, {.cycles = 1}, &info);
  ASSERT_TRUE(r.Ok()) << r.message;
  EXPECT_FALSE(r.css_reached);
  EXPECT_TRUE(info.warm);
  EXPECT_FALSE(info.cached);
  EXPECT_EQ(cache.Size(), 2u);

  // A state of another deck structure is refused: cold run, with a note.
  Simulation two = TwoLayerSimulation(4, 4);
  ASSERT_TRUE(two.Run({.cycles = 1}).Ok());
  CycleStateCache other({500.0}, {700.0});
  other.Add({600.0}, two.FinalState());
  Simulation d = Tsa(600.0);
  std::vector<std::string> messages;
  RunOptions options;
  options.cycles = 1;
  options.on_message = [&](LogLevel, const std::string& m) { messages.push_back(m); };
  r = RunWithStateCache(d, &other, {600.0}, options, &info);
  ASSERT_TRUE(r.Ok()) << r.message;
  EXPECT_FALSE(info.warm);
  EXPECT_TRUE(info.skipped);
  EXPECT_NE(info.note.find("warm start skipped"), std::string::npos) << info.note;
  ASSERT_FALSE(messages.empty());
  EXPECT_EQ(messages.front(), info.note);

  // The cache chooses the start.
  CycleState x = a.FinalState();
  options.initial_state = &x;
  EXPECT_THROW(RunWithStateCache(d, &cache, {600.0}, options), std::invalid_argument);
  // No cache: a plain cold run.
  r = RunWithStateCache(d, nullptr, {600.0}, {.cycles = 1}, &info);
  EXPECT_TRUE(r.Ok());
  EXPECT_FALSE(info.warm);
  EXPECT_FALSE(info.cached);
}

namespace {

// The TSA example on a coarse grid, run to CSS: the test is about the
// optimisation loop, not the optimum.
SimulationOptimizerSettings SmallOptimisation(bool warm) {
  SimulationOptimizerSettings s;
  s.n_init_samples = 2;
  s.n_iterations = 1;
  s.warm_start = warm;
  s.settings = {{"layer[1].cells", "6"}};
  return s;
}

}  // namespace

TEST(WarmStartCache, OptimisationWarmStartsCandidatesAndReportsTheColdOptimum) {
  const std::string deck = ReadDeck("Illustrative_TSA_13X.json");
  ASSERT_FALSE(deck.empty());
  Simulation sim = Simulation::FromJson(deck);
  ASSERT_EQ(sim.Data().kpi_config.per_component.empty(), false);
  std::string component;
  for (const auto& c : sim.Data().kpi_config.per_component)
    if (c.tracked) component = c.product_component;
  ASSERT_FALSE(component.empty());

  const std::string path = "step[1].time";
  const double t = sim.GetParameter(path);
  const std::vector<OptimizationVariable> vars = {{path, 0.8 * t, 1.2 * t}};
  OptimizationCostFunction cost;
  cost.AddObjective(component, KPIType::RECOVERY, 1.0, ObjectiveDirection::MAXIMIZE);

  int calls = 0;
  SimulationOptimizerSettings warm = SmallOptimisation(true);
  warm.on_evaluation = [&](const OptimizationEvaluation&) { ++calls; };
  const SimulationOptimizationResult r = OptimizeSimulation(deck, vars, cost, warm);
  ASSERT_EQ(r.evaluations.size(), 3u);
  EXPECT_EQ(calls, 4);  // three evaluations and the cold re-run
  EXPECT_FALSE(r.evaluations[0].warm);
  for (size_t i = 1; i < r.evaluations.size(); ++i) {
    EXPECT_TRUE(r.evaluations[i].warm) << i << ": " << r.evaluations[i].note;
    EXPECT_TRUE(r.evaluations[i].accepted) << r.evaluations[i].message;
    EXPECT_TRUE(r.evaluations[i].css_reached);
    EXPECT_LT(r.evaluations[i].cycles_run, r.evaluations[0].cycles_run);
    EXPECT_TRUE(std::isfinite(r.evaluations[i].distance));
  }
  ASSERT_TRUE(r.rerun_cold);
  EXPECT_FALSE(r.cold.warm);
  EXPECT_EQ(r.cold.index, -1);
  EXPECT_EQ(r.cold.x, r.x);
  EXPECT_EQ(r.cost, r.cold.cost);  // the cold value is the reported one
  ASSERT_GE(r.best_index, 0);
  EXPECT_EQ(r.search_cost, r.evaluations[r.best_index].cost);
  EXPECT_GE(r.x[0], vars[0].lower);
  EXPECT_LE(r.x[0], vars[0].upper);
  EXPECT_LT(r.cost, 0.0);  // -recovery
  EXPECT_FALSE(r.cancelled);
  // The cold re-run agrees with the search's warm value within the CSS bound
  // (c <= 100 state tolerances, docs/cyclic-steady-state.md) -- unless the
  // optimum is the first candidate, which ran cold too.
  EXPECT_NEAR(r.cost, r.search_cost, 100.0 * sim.GetParameter("css.state_tolerance"));

  SimulationOptimizerSettings no_warm = SmallOptimisation(false);
  no_warm.n_iterations = 0;
  no_warm.rerun_cold = false;
  const SimulationOptimizationResult cold = OptimizeSimulation(deck, vars, cost, no_warm);
  EXPECT_FALSE(cold.rerun_cold);
  EXPECT_EQ(cold.cost, cold.search_cost);
  for (const auto& e : cold.evaluations) EXPECT_FALSE(e.warm);
  // The first candidate is the same and runs cold either way.
  EXPECT_EQ(cold.evaluations[0].x, r.evaluations[0].x);
  EXPECT_EQ(cold.evaluations[0].cost, r.evaluations[0].cost);

  // The set-up is checked before any run.
  EXPECT_THROW(CheckOptimization(deck, {{"step[999].time", 1, 2}}, cost), std::invalid_argument);
  EXPECT_THROW(CheckOptimization(deck, {{path, 2, 1}}, cost), std::invalid_argument);
  EXPECT_THROW(CheckOptimization(deck, {vars[0], vars[0]}, cost), std::invalid_argument);
  EXPECT_THROW(CheckOptimization(deck, vars, OptimizationCostFunction{}), std::invalid_argument);
  OptimizationCostFunction untracked;
  untracked.AddObjective("Helium", KPIType::PURITY, 1.0, ObjectiveDirection::MAXIMIZE);
  EXPECT_THROW(CheckOptimization(deck, vars, untracked), std::invalid_argument);
  OptimizationCostFunction lcoc;
  lcoc.AddObjective(component, KPIType::LCOC, 1.0, ObjectiveDirection::MINIMIZE);
  EXPECT_THROW(CheckOptimization(deck, vars, lcoc), std::invalid_argument);
}

TEST(WarmStartCache, AWarmSensitivitySampleMatchesTheColdOne) {
  const std::string deck = ReadDeck("Illustrative_TSA_13X.json");
  Simulation sim = Simulation::FromJson(deck);
  const std::string path = "step[1].time";
  const std::vector<sensitivity::Input> inputs = {sensitivity::Relative(path, sim.GetParameter(path), 0.1)};
  const std::vector<sensitivity::Output> outputs = sensitivity::DefaultOutputs(sim.Data());
  ASSERT_FALSE(outputs.empty());
  sensitivity::RunPolicy policy;
  policy.settings = {{"layer[1].cells", "6"}};

  CycleStateCache cache;
  const auto nominal = sensitivity::EvaluateSample(deck, inputs, {}, outputs, policy, {}, &cache, {0.5});
  ASSERT_TRUE(nominal.accepted) << nominal.message;
  EXPECT_FALSE(nominal.warm);
  ASSERT_EQ(cache.Size(), 1u);

  const std::vector<double> value = {inputs[0].Quantile(0.75)};
  const auto warm = sensitivity::EvaluateSample(deck, inputs, value, outputs, policy, {}, &cache, {0.75});
  const auto cold = sensitivity::EvaluateSample(deck, inputs, value, outputs, policy);
  ASSERT_TRUE(warm.accepted) << warm.message;
  ASSERT_TRUE(cold.accepted) << cold.message;
  EXPECT_TRUE(warm.warm);
  EXPECT_FALSE(cold.warm);
  EXPECT_LT(warm.cycles_run, cold.cycles_run);
  // Purity and recovery agree within the CSS bound (c <= 100 state
  // tolerances, docs/cyclic-steady-state.md).
  const double tol = 100.0 * sim.GetParameter("css.state_tolerance");
  for (size_t o = 0; o < outputs.size(); ++o) {
    if (outputs[o].kpi != "purity" && outputs[o].kpi != "recovery") continue;
    EXPECT_NEAR(warm.outputs[o], cold.outputs[o], tol) << outputs[o].Name();
  }
  EXPECT_THROW(sensitivity::EvaluateSample(deck, inputs, value, outputs, policy, {}, &cache, {}),
               std::invalid_argument);
}

}  // namespace adsim
