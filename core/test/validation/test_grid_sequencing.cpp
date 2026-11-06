// test_grid_sequencing.cpp
//
// Grid sequencing: does running a coarse grid to cyclic steady state and
// warm-starting a finer one from it reach the fine-grid CSS faster than the
// fine grid from cold? Three paths per deck, each ending in a run to CSS on
// the deck's own (fine) grid of N cells per layer:
//   (a) cold:  N cells from the deck's initial conditions;
//   (b) 2-level: N/2 cells to CSS, then N cells warm-started from that state;
//   (c) 3-level: N/4 -> N/2 -> N cells, each level warm-started from the last.
// A warm start onto another grid goes through the warm-start API, which
// remaps the state per layer (PrepareCycleState / RemapCycleState); the
// state's cycle count carries over, as it does for any warm start.
//
// Asserted (no wall-clock assertions, which would be flaky):
//   * every path reaches CSS, and its final fine-grid state, histories and
//     KPIs agree with path (a) within c <= 100 times the state tolerance,
//     the bound of CSSReference and CSSStartIndependence;
//   * a warm-started fine level needs fewer fine-grid cycles than path (a).
// The cycles per level and the wall time of each path are printed; the table
// in docs/cyclic-steady-state.md comes from this output.
//
// Decks: six of the CSSReference decks (TSA, VPSA, PSA and two slower
// cycles). Slow tier: about 2-8 minutes per deck.

#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/simulation.h"
#include "css_state_compare.h"

namespace adsim {

  using namespace adsim::validation;

  namespace {

    constexpr double kBound = 100.0;  // as CSSReference
    constexpr int kMaxCycles = 1000;

    struct Level {
      int cells = 0;       // cells of the first layer (the decks here have one)
      int css_cycle = -1;  // CSS stop, counted from the start of this level's run
      int cycles = 0;      // cycles run at this level
      double seconds = 0.0;
    };

    struct Path {
      std::vector<Level> levels;
      CycleState state;
      std::vector<KPIRecord> kpis;
      double eps = 0.0;
      double Seconds() const {
        double s = 0.0;
        for (const auto& l : levels) s += l.seconds;
        return s;
      }
      int Cycles() const {
        int c = 0;
        for (const auto& l : levels) c += l.cycles;
        return c;
      }
    };

    std::string DeckPath(const std::string& deck) { return "../examples/decks/" + deck + ".json"; }

    // One level: `deck` with every layer's cell count scaled by `fraction`
    // (rounded), run to CSS from `from` when given.
    Level RunLevel(const std::string& deck, double fraction, const CycleState* from, Path& path) {
      Simulation sim = Simulation::FromFile(DeckPath(deck));
      auto& data = sim.Data();
      if (data.settings.css.method_ == CSSMethod::NONE) data.settings.css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
      data.settings.acc_css.method_ = AccelerationMethod::NONE;
      Level level;
      for (auto& layer : data.bed.layers_.layers_) {
        layer.domain_.n_cells = std::max(4, static_cast<int>(std::lround(layer.domain_.n_cells * fraction)));
        if (level.cells == 0) level.cells = layer.domain_.n_cells;
      }
      RunOptions options;
      options.cycles = kMaxCycles;
      options.initial_state = from;
      options.on_message = [](LogLevel, const std::string&) {};
      const auto start = std::chrono::steady_clock::now();
      const RunResult run = sim.Run(options);
      level.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
      EXPECT_TRUE(run.Ok()) << deck << ": " << run.message;
      path.eps = data.settings.css.state_tolerance_.Get();
      if (!run.Ok()) return level;
      level.css_cycle = run.css_reached ? run.css_cycle : -1;
      level.cycles = run.cycles_run;
      path.state = sim.FinalState();
      path.kpis = sim.KPIs(-1);
      return level;
    }

    Path RunPath(const std::string& deck, const std::vector<double>& fractions) {
      Path path;
      for (double f : fractions) {
        const CycleState from = path.state;  // copy: RunLevel replaces path.state
        path.levels.push_back(RunLevel(deck, f, from.Empty() ? nullptr : &from, path));
      }
      return path;
    }

    void Print(const std::string& deck, const char* name, const Path& p) {
      std::string levels;
      for (const auto& l : p.levels)
        levels += std::to_string(l.cells) + ":" + std::to_string(l.cycles) + "cyc/" +
                  std::to_string(static_cast<int>(std::lround(l.seconds))) + "s ";
      std::printf("[GRIDSEQ] %-24s %-8s total %4d cycles %7.1f s  | %s\n", deck.c_str(), name, p.Cycles(),
                  p.Seconds(), levels.c_str());
    }

  }  // namespace

  class GridSequencing : public ::testing::TestWithParam<const char*> {};

  TEST_P(GridSequencing, CoarseToFineReachesTheSameCss) {
    const std::string deck = GetParam();
    const Path cold = RunPath(deck, {1.0});
    const Path two = RunPath(deck, {0.5, 1.0});
    const Path three = RunPath(deck, {0.25, 0.5, 1.0});
    Print(deck, "cold", cold);
    Print(deck, "2-level", two);
    Print(deck, "3-level", three);
    ASSERT_GE(cold.levels.back().css_cycle, 0) << deck << ": the cold fine run did not reach CSS";
    const double eps = cold.eps;

    for (const auto* p : {&two, &three}) {
      const char* name = p == &two ? "2-level" : "3-level";
      for (const auto& l : p->levels)
        EXPECT_GE(l.css_cycle, 0) << deck << " " << name << ": no CSS at " << l.cells << " cells";
      const CssDistance d = Distance(p->state, p->kpis, cold.state, cold.kpis);
      std::printf("[GRIDSEQ] %-24s %-8s vs cold: c_state %.3f (%s) c_history %.3f c_kpi %.3f\n", deck.c_str(), name,
                  d.state / eps, d.state_field.c_str(), d.history / eps, d.kpi / eps);
      EXPECT_LE(d.state, kBound * eps) << deck << " " << name << ": bed state (" << d.state_field << ")";
      EXPECT_LE(d.history, kBound * eps) << deck << " " << name << ": endpoint histories";
      EXPECT_LE(d.kpi, kBound * eps) << deck << " " << name << ": KPIs";
      // The point of sequencing: the warm-started fine level runs fewer
      // fine-grid cycles than the cold fine run.
      EXPECT_LT(p->levels.back().cycles, cold.levels.back().cycles) << deck << " " << name;
    }
  }

  INSTANTIATE_TEST_SUITE_P(ExampleDecks, GridSequencing,
                           ::testing::Values("Illustrative_TSA_13X", "Perez2019_A", "Perez2019_E", "Nguyen2023_Exp_I",
                                             "Vysyaraju2024_PSA_exp1", "Ward2024"),
                           [](const ::testing::TestParamInfo<const char*>& info) { return std::string(info.param); });

}  // namespace adsim
