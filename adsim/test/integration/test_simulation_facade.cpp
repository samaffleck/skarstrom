// test_simulation_facade.cpp
//
// The Simulation facade (adsim/simulation.h): load/save, parameters through
// the registry, blocking and chunked runs, progress and cancel, and the typed
// results (KPIs, CSS history, series, profiles, stream logs). The KPI and
// series checks compare against the same quantities computed straight from a
// CoupledSolver run, so the facade is pinned to the engine's own numbers.

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/simulation.h"
#include "adsim/solver/kpis.h"
#include "adsim/solver/parameterisation.h"

namespace adsim {

namespace {

const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";

std::string Deck(const std::string& name) { return (kExamples / name).string(); }

std::string ReadFile(const std::filesystem::path& p) {
  std::ifstream in(p);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

// Perez 2019 A (four-step CO2 VSA, CO2 tracked) cut to three cycles: < 1 s.
Simulation ShortPerez() {
  Simulation sim = Simulation::FromFile(Deck("Perez2019_A.json"));
  sim.SetParameter("cycle.number_of_cycles", 3);
  return sim;
}

}  // namespace

TEST(SimulationFacade, LoadsBothJsonRootsAndRoundTrips) {
  Simulation a = Simulation::FromFile(Deck("Ward_2022_Q=100[ccm].json"));  // "value0" root
  const std::string json = a.ToJson();                                      // "simulation" root
  ASSERT_NE(json.find("\"simulation\""), std::string::npos);

  Simulation b = Simulation::FromJson(json);
  Simulation c = Simulation::FromJson(ReadFile(Deck("Ward_2022_Q=100[ccm].json")));
  for (Simulation* s : {&b, &c}) {
    EXPECT_EQ(s->Components(), a.Components());
    EXPECT_EQ(s->StepNames(), a.StepNames());
    EXPECT_EQ(s->Parameters().Size(), a.Parameters().Size());
    for (const auto& p : a.Parameters().All()) {
      const Parameter* q = s->Parameters().Find(p.Path());
      ASSERT_NE(q, nullptr) << p.Path();
      EXPECT_EQ(q->GetText(), p.GetText()) << p.Path();
    }
  }

  EXPECT_THROW(Simulation::FromJson(""), std::invalid_argument);
  EXPECT_ANY_THROW(Simulation::FromJson("{\"nope\": 1}"));
  EXPECT_ANY_THROW(Simulation::FromFile(Deck("does_not_exist.json")));
}

TEST(SimulationFacade, SaveWritesALoadableDeck) {
  Simulation a = Simulation::FromFile(Deck("Ward_2022_Q=100[ccm].json"));
  a.SetParameter("settings.reltol", 2e-6);
  const auto path = std::filesystem::temp_directory_path() / "adsim_facade_save.json";
  a.Save(path.string());
  Simulation b = Simulation::FromFile(path.string());
  EXPECT_DOUBLE_EQ(b.GetParameter("settings.reltol"), 2e-6);
  std::filesystem::remove(path);
}

TEST(SimulationFacade, ParametersGetSetAndValidate) {
  Simulation sim = Simulation::FromFile(Deck("Perez2019_A.json"));
  // Static, integer and string fields.
  sim.SetParameter("cycle.number_of_cycles", 4);
  EXPECT_EQ(sim.Data().cycle.number_of_cycles_, 4);
  sim.SetParameterText("info.title", "renamed");
  EXPECT_EQ(sim.Data().info.title, "renamed");
  sim.SetParameter("css.min_cycles", 2);
  EXPECT_EQ(sim.Data().settings.css.min_cycles_, 2);

  // A model parameter by ID path and by alias.
  const Step& first = sim.Data().cycle.step_list_.front();
  const std::string time_path = "step[1].time";
  ASSERT_NE(sim.Parameters().Find(time_path), nullptr);
  sim.SetParameter(time_path, 33.0);
  EXPECT_DOUBLE_EQ(first.step_time_.Get(), 33.0);
  sim.SetParameter("step[" + first.step_name_ + "].time", 34.0);
  EXPECT_DOUBLE_EQ(first.step_time_.Get(), 34.0);

  // Expressions with inline units evaluate; a bad expression is rejected and
  // the old value restored.
  sim.SetParameterText(time_path, "1 [min]");
  EXPECT_DOUBLE_EQ(sim.GetParameter(time_path), 60.0);
  EXPECT_THROW(sim.SetParameterText(time_path, "1 +* 2"), std::invalid_argument);
  EXPECT_DOUBLE_EQ(sim.GetParameter(time_path), 60.0);
  EXPECT_EQ(sim.GetParameterText(time_path), "1 [min]");

  EXPECT_THROW(sim.GetParameter("no.such.parameter"), std::out_of_range);
  EXPECT_THROW(sim.SetParameterText("cycle.number_of_cycles", "three"), std::invalid_argument);
}

TEST(SimulationFacade, BlockingRunMatchesTheEngine) {
  Simulation sim = ShortPerez();
  std::vector<RunProgress> progress;
  RunOptions opt;
  opt.progress = [&](const RunProgress& p) { progress.push_back(p); return true; };
  const RunResult r = sim.Run(opt);
  ASSERT_TRUE(r.Ok()) << r.message;
  EXPECT_EQ(r.cycles_run, 3);
  EXPECT_GT(r.stats.num_steps, 0);
  EXPECT_GT(r.wall_time_s, 0.0);

  const int steps = static_cast<int>(sim.Data().cycle.step_list_.size());
  ASSERT_EQ(static_cast<int>(progress.size()), 3 * steps);
  EXPECT_DOUBLE_EQ(progress.back().fraction, 1.0);
  EXPECT_EQ(sim.Cycles(), (std::vector<int>{0, 1, 2}));
  EXPECT_EQ(sim.LastCycle(), 2);

  // Same deck straight through CoupledSolver.
  UserData ref;
  LoadFromJson(ref, Deck("Perez2019_A.json"));
  ref.cycle.number_of_cycles_ = 3;
  std::string eval_err;
  ASSERT_TRUE(ref.EvaluateAllExpressions(eval_err)) << eval_err;
  CoupledSolver solver(ref);
  ASSERT_GE(solver.IntegrateCycle(ref.cycle), 0);

  const auto kpis = sim.KPIs();
  ASSERT_EQ(kpis.size(), 3u);  // CO2 tracked, three cycles
  for (const auto& k : kpis) {
    EXPECT_EQ(k.component, "CarbonDioxide");
    ASSERT_TRUE(k.valid);
    const auto* cfg = ref.kpi_config.FindByComponent(k.component);
    ASSERT_NE(cfg, nullptr);
    KPIParams p(solver.GetLog(), solver.GetBoundaryLog());
    p.cycle_index_ = k.cycle;
    p.product_comp_ = k.component;
    p.feed_endpoints = cfg->feed_endpoints;
    p.product_endpoints = cfg->product_endpoints;
    EXPECT_DOUBLE_EQ(k.purity, CalculatePurity(p));
    EXPECT_DOUBLE_EQ(k.recovery, CalculateRecovery(p));
    EXPECT_DOUBLE_EQ(k.productivity_mass, CalculateProductivityMass(p));
    EXPECT_GT(k.purity, 0.0);
    EXPECT_LE(k.purity, 1.0);
  }
  EXPECT_EQ(sim.KPIs(-1).front().cycle, 2);

  std::vector<double> t, v;
  ASSERT_TRUE(ExtractComponentTimeSeries(solver.GetLog(), ComponentVariableType::MOLE_FRACTION, "CarbonDioxide",
                                         ParamSpatialLoc::Top, 0.0, t, v));
  const TimeSeries s = sim.Series("mole_fraction", "CarbonDioxide", SeriesLocation::Top);
  EXPECT_EQ(s.time, t);
  EXPECT_EQ(s.values, v);
  EXPECT_EQ(s.unit, "mol/mol");

  // Earlier cycle, interior position, bulk variable.
  const TimeSeries p0 = sim.Series("pressure", "", SeriesLocation::Position, 0.1, 0);
  ASSERT_FALSE(p0.time.empty());
  EXPECT_DOUBLE_EQ(p0.time.front(), 0.0);
  EXPECT_EQ(p0.unit, "Pa");
  EXPECT_THROW(sim.Series("pressure", "CarbonDioxide"), std::invalid_argument);
  EXPECT_THROW(sim.Series("mole_fraction"), std::invalid_argument);
  EXPECT_THROW(sim.Series("bogus"), std::invalid_argument);
  EXPECT_THROW(sim.Series("pressure", "", SeriesLocation::Top, 0.0, 7), std::out_of_range);
}

TEST(SimulationFacade, ProfilesStreamsAndCSSHistory) {
  Simulation sim = ShortPerez();
  ASSERT_TRUE(sim.Run().Ok());

  const AxialProfile prof = sim.Profile("temperature", "", -1, 1);
  EXPECT_EQ(prof.cycle, 2);
  EXPECT_EQ(prof.step, 1);
  EXPECT_EQ(prof.cols, sim.Positions().size());
  EXPECT_EQ(prof.rows, prof.time.size());
  EXPECT_EQ(prof.values.size(), prof.rows * prof.cols);
  for (double x : prof.values) EXPECT_GT(x, 200.0);
  const AxialProfile q = sim.Profile("adsorbed_concentration", "CarbonDioxide", 0, 0);
  EXPECT_EQ(q.unit, "mol/kg");

  const auto streams = sim.Streams();
  ASSERT_FALSE(streams.empty());
  bool any = false;
  for (const auto& st : streams) {
    const TimeSeries f = sim.StreamSeries(st.step, st.name, "molar_flow_rate");
    EXPECT_EQ(f.time.size(), f.values.size());
    const TimeSeries y = sim.StreamSeries(st.step, st.name, "mole_fraction", "CarbonDioxide", 1);
    EXPECT_EQ(y.unit, "mol/mol");
    any = any || !f.time.empty();
  }
  EXPECT_TRUE(any);
  EXPECT_THROW(sim.StreamSeries(0, "no such stream", "molar_flow_rate"), std::out_of_range);

  const auto css = sim.CSSHistory();
  ASSERT_EQ(css.size(), 3u);
  EXPECT_EQ(css[2].cycle, 2);
}

TEST(SimulationFacade, ProgressCallbackAndCancelPollCancel) {
  {
    Simulation sim = ShortPerez();
    RunOptions opt;
    int calls = 0;
    opt.progress = [&](const RunProgress&) { return ++calls < 2; };
    const RunResult r = sim.Run(opt);
    EXPECT_EQ(r.status, RunStatus::Cancelled);
    EXPECT_EQ(calls, 2);
  }
  {
    Simulation sim = ShortPerez();
    RunOptions opt;
    int polls = 0;
    opt.cancel = [&] { return ++polls >= 50; };
    const RunResult r = sim.Run(opt);
    EXPECT_EQ(r.status, RunStatus::Cancelled);
    EXPECT_EQ(polls, 50);
  }
}

TEST(SimulationFacade, CyclesOptionAndFailedChecks) {
  Simulation sim = Simulation::FromFile(Deck("Perez2019_A.json"));
  RunOptions opt;
  opt.cycles = 1;
  const RunResult r = sim.Run(opt);
  ASSERT_TRUE(r.Ok()) << r.message;
  EXPECT_EQ(r.cycles_run, 1);

  // A deck that fails its checks reports Failed with the reason.
  Simulation bad = Simulation::FromFile(Deck("Perez2019_A.json"));
  bad.Data().cycle.step_list_.clear();
  const RunResult rb = bad.Run();
  EXPECT_EQ(rb.status, RunStatus::Failed);
  EXPECT_FALSE(rb.message.empty());
  std::string why;
  EXPECT_FALSE(bad.Check(why));
  EXPECT_FALSE(why.empty());
}

TEST(SimulationFacade, ChunkedRunCompletes) {
  Simulation sim = Simulation::FromFile(Deck("Ward_2022_Q=100[ccm].json"));
  ASSERT_TRUE(sim.Start());
  EXPECT_TRUE(sim.IsRunning());
  EXPECT_THROW(sim.SetParameter("settings.reltol", 1e-5), std::logic_error);
  RunStatus st = RunStatus::Running;
  for (int i = 0; i < 100000 && st == RunStatus::Running; ++i) st = sim.Advance(50);
  EXPECT_TRUE(st == RunStatus::Completed || st == RunStatus::ReachedCSS) << sim.LastResult().message;
  EXPECT_FALSE(sim.IsRunning());
  EXPECT_DOUBLE_EQ(sim.ProgressPercent(), 100.0);
  EXPECT_TRUE(sim.HasResults());
  EXPECT_FALSE(sim.Series("mole_fraction", "CarbonDioxide").values.empty());
  EXPECT_FALSE(sim.Messages().empty());
}

}  // namespace adsim
