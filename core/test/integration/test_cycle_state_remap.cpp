// test_cycle_state_remap.cpp
//
// Transfer of a CycleState onto other cell counts (RemapCycleState,
// Simulation::RemapState; a run started from a state on another grid remaps it
// itself):
//
//  * every component's gas and adsorbed inventory in every layer is kept to
//    rounding, fine -> coarse and coarse -> fine, ratios 2 and 1.5, on a deck
//    with two layers of different cell counts, a second adsorbent present in
//    one layer only and a void region;
//  * refining by an integer ratio and coarsening back gives the state back;
//  * the remapped state starts a run cleanly (two-layer deck and example decks);
//  * a warm start on a refined grid from the coarse grid's CSS reaches CSS in
//    fewer cycles than a cold start;
//  * histories and summaries are carried over unchanged; a state that does not
//    fit (components, layers, void regions, steps) is still refused.

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

constexpr double kGasConstant = 8.31446261815324;

// The state after two cycles on (cells1, cells2), computed once per grid.
const CycleState& TwoLayerState(int cells1, int cells2) {
  static std::map<std::pair<int, int>, CycleState> cache;
  auto it = cache.find({cells1, cells2});
  if (it != cache.end()) return it->second;
  Simulation sim = TwoLayerSimulation(cells1, cells2);
  const RunResult r = sim.Run();
  EXPECT_TRUE(r.Ok()) << r.message;
  return cache[{cells1, cells2}] = sim.FinalState();
}

// Inventory of each layer (per unit area, up to each layer's constant
// porosity and solid density): per component, the gas sum of c_i dz and, per
// adsorbent, the adsorbed sum of q dz. dz = layer length / cells.
struct LayerInventory {
  std::map<std::string, double> gas;                          // component -> sum c_i dz
  std::map<std::pair<size_t, std::string>, double> adsorbed;  // (adsorbent, component) -> sum q dz
};

std::vector<LayerInventory> Inventory(const CycleState& s, const std::vector<double>& lengths) {
  std::vector<LayerInventory> out(s.layers.size());
  // First cell of each layer.
  std::vector<size_t> first(s.layers.size());
  size_t c = 1 + (s.bottom_void ? 1 : 0);
  for (size_t l = 0; l < s.layers.size(); ++l) { first[l] = c; c += static_cast<size_t>(s.layers[l].cells); }
  auto layer_of = [&](size_t cell) -> int {
    for (size_t l = 0; l < s.layers.size(); ++l)
      if (cell >= first[l] && cell < first[l] + static_cast<size_t>(s.layers[l].cells)) return static_cast<int>(l);
    return -1;
  };
  auto cell_values = [&](const CycleStateField* f) {
    std::map<size_t, double> m;
    if (f) for (size_t k = 0; k < f->index.size(); ++k) m[f->cell[k]] = s.y[f->index[k]];
    return m;
  };
  const auto P = cell_values(s.FindField("Pressure"));
  const auto T = cell_values(s.FindField("Fluid Temperature"));
  const auto C = cell_values(s.FindField("Concentration"));
  for (const auto& comp : s.components) {
    for (const auto& [cell, y] : cell_values(s.FindField("Mole Fraction", comp))) {
      const int l = layer_of(cell);
      if (l < 0) continue;
      const double total = C.empty() ? P.at(cell) / (kGasConstant * T.at(cell)) : C.at(cell);
      out[l].gas[comp] += total * y * lengths[l] / s.layers[l].cells;
    }
  }
  for (const auto& f : s.fields) {
    if (f.adsorbent == 0) continue;
    for (size_t k = 0; k < f.index.size(); ++k) {
      const int l = layer_of(f.cell[k]);
      if (l < 0) continue;
      out[l].adsorbed[{f.adsorbent, f.component}] += s.y[f.index[k]] * lengths[l] / s.layers[l].cells;
    }
  }
  return out;
}

void ExpectSameInventory(const CycleState& a, const CycleState& b, const char* label) {
  const std::vector<double> lengths = {0.6, 0.4};
  const auto ia = Inventory(a, lengths), ib = Inventory(b, lengths);
  ASSERT_EQ(ia.size(), ib.size());
  for (size_t l = 0; l < ia.size(); ++l) {
    ASSERT_EQ(ia[l].gas.size(), 2u) << label;
    for (const auto& [comp, v] : ia[l].gas)
      EXPECT_NEAR(ib[l].gas.at(comp), v, 1e-12 * std::abs(v)) << label << ": gas " << comp << " in layer " << l;
    ASSERT_EQ(ia[l].adsorbed.size(), l == 0 ? 2u : 4u) << label;  // B only in the top layer
    for (const auto& [key, v] : ia[l].adsorbed)
      EXPECT_NEAR(ib[l].adsorbed.at(key), v, 1e-12 * std::abs(v))
          << label << ": adsorbed " << key.second << " on adsorbent " << key.first << " in layer " << l;
  }
}

}  // namespace

TEST(CycleStateRemap, KeepsEveryLayersInventory) {
  struct Case { int from1, from2, to1, to2; const char* label; };
  const Case cases[] = {
      {12, 8, 6, 4, "fine -> coarse, ratio 2"},
      {6, 4, 12, 8, "coarse -> fine, ratio 2"},
      {12, 8, 8, 12, "uneven: 1.5 down in layer 1, 1.5 up in layer 2"},
      {8, 12, 12, 8, "uneven: 1.5 up in layer 1, 1.5 down in layer 2"},
  };
  for (const auto& k : cases) {
    const CycleState& x = TwoLayerState(k.from1, k.from2);
    Simulation target = TwoLayerSimulation(k.to1, k.to2);
    EXPECT_EQ(target.CheckState(x), "") << k.label;
    const CycleState r = target.RemapState(x);
    EXPECT_EQ(r.layers[0].cells, k.to1);
    EXPECT_EQ(r.layers[1].cells, k.to2);
    EXPECT_EQ(r.cycles, x.cycles);
    EXPECT_EQ(r.y.size(), r.yp.size());
    EXPECT_TRUE(std::all_of(r.y.begin(), r.y.end(), [](double v) { return std::isfinite(v); })) << k.label;
    ExpectSameInventory(x, r, k.label);

    // Histories and summaries are per step: carried over unchanged.
    ASSERT_EQ(r.endpoints.size(), x.endpoints.size());
    for (size_t e = 0; e < r.endpoints.size(); ++e) {
      EXPECT_EQ(r.endpoints[e].endpoint_id, x.endpoints[e].endpoint_id);
      EXPECT_EQ(r.endpoints[e].history.time_series_, x.endpoints[e].history.time_series_);
      EXPECT_EQ(r.endpoints[e].summary.total_moles, x.endpoints[e].summary.total_moles);
    }
    // Mole fractions stay in [0, 1] and sum to one in every bed cell (the remap is monotone).
    const auto yo = r.FieldValues("Mole Fraction", "Oxygen");
    const auto yn = r.FieldValues("Mole Fraction", "Nitrogen");
    for (size_t i = 0; i < yo.size(); ++i) {
      EXPECT_GE(yo[i], -1e-12);
      EXPECT_LE(yo[i], 1.0 + 1e-12);
      EXPECT_NEAR(yo[i] + yn[i], 1.0, 1e-9) << k.label << " cell " << i;
    }
  }
}

TEST(CycleStateRemap, PressureIsTheMonotoneCellAverage) {
  // Illustrative_TSA_13X at the end of a cycle has a sharp temperature front
  // near the feed end. Remapping C (= P / R T) and T separately overshot the
  // pressure there (by about 20 Pa above both neighbours on 30 -> 45 cells),
  // and IDA's consistent-IC solve then failed; with the concentration-weighted
  // temperature P is the cell average of P.
  auto make = [](int cells) {
    Simulation sim = Simulation::FromFile("../examples/decks/Illustrative_TSA_13X.json");
    for (auto& layer : sim.Data().bed.layers_.layers_) layer.domain_.n_cells = cells;
    sim.Data().settings.css.method_ = CSSMethod::NONE;
    return sim;
  };
  Simulation coarse = make(30);
  ASSERT_TRUE(coarse.Run({.cycles = 1}).Ok());
  const CycleState x = coarse.FinalState();
  for (int cells : {45, 20}) {
    Simulation fine = make(cells);
    const CycleState r = fine.RemapState(x);
    // Each target cell's pressure lies between the pressures of the source
    // cells it overlaps (bed cells 1..n; ghosts excluded).
    const auto p_src = x.FieldValues("Pressure"), p_dst = r.FieldValues("Pressure");
    const int n = 30, m = cells;
    for (int j = 0; j < m; ++j) {
      const int i0 = j * n / m, i1 = std::min(n - 1, ((j + 1) * n - 1) / m);
      double lo = 1e300, hi = -1e300;
      for (int i = i0; i <= i1; ++i) { lo = std::min(lo, p_src[1 + i]); hi = std::max(hi, p_src[1 + i]); }
      EXPECT_GE(p_dst[1 + j], lo - 1e-9 * hi) << cells << " cells, cell " << j;
      EXPECT_LE(p_dst[1 + j], hi + 1e-9 * hi) << cells << " cells, cell " << j;
    }
    const RunResult run = fine.RunCycle(x);
    EXPECT_TRUE(run.Ok()) << "30 -> " << cells << " cells: " << run.message;
  }
}

TEST(CycleStateRemap, RefineThenCoarsenIsTheIdentity) {
  for (int ratio : {2, 3}) {
    const CycleState& x = TwoLayerState(6, 4);
    Simulation fine = TwoLayerSimulation(6 * ratio, 4 * ratio);
    const CycleState up = fine.RemapState(x);
    Simulation coarse = TwoLayerSimulation(6, 4);
    const CycleState back = coarse.RemapState(up);
    ASSERT_EQ(back.y.size(), x.y.size());
    double worst = 0.0;
    for (size_t g = 0; g < x.y.size(); ++g)
      worst = std::max(worst, std::abs(back.y[g] - x.y[g]) / std::max(1e-8, std::abs(x.y[g])));
    std::cout << "[ remap ] ratio " << ratio << ": refine + coarsen, largest relative change " << worst << "\n";
    EXPECT_LT(worst, 1e-12) << "ratio " << ratio;
  }
}

TEST(CycleStateRemap, SameGridIsACopy) {
  const CycleState& x = TwoLayerState(6, 4);
  Simulation same = TwoLayerSimulation(6, 4);
  const CycleState r = same.RemapState(x);
  EXPECT_EQ(r.y, x.y);
}

TEST(CycleStateRemap, RemappedStateStartsARun) {
  const CycleState& x = TwoLayerState(12, 8);
  for (auto [c1, c2] : {std::pair{6, 4}, std::pair{18, 12}, std::pair{8, 12}}) {
    Simulation sim = TwoLayerSimulation(c1, c2);
    const RunResult r = sim.RunCycle(x);  // the run remaps the state itself
    ASSERT_TRUE(r.Ok()) << c1 << "/" << c2 << ": " << r.message;
    EXPECT_EQ(sim.FinalState().layers[0].cells, c1);
    EXPECT_EQ(sim.FinalState().cycles, x.cycles + 1);
  }
}

TEST(CycleStateRemap, RemappedStateStartsARunOnExampleDecks) {
  for (const char* deck : {"Vysyaraju2024_PSA_exp1", "Illustrative_TSA_13X", "Perez2019_A"}) {
    auto make = [&](int cells) {
      Simulation sim = Simulation::FromFile(std::string("../examples/decks/") + deck + ".json");
      for (auto& layer : sim.Data().bed.layers_.layers_) layer.domain_.n_cells = cells;
      sim.Data().settings.css.method_ = CSSMethod::NONE;
      return sim;
    };
    Simulation coarse = make(10);
    ASSERT_TRUE(coarse.Run({.cycles = 1}).Ok()) << deck << ": " << coarse.LastResult().message;
    const CycleState x = coarse.FinalState();
    for (int cells : {15, 20}) {
      Simulation fine = make(cells);
      const RunResult r = fine.RunCycle(x);
      EXPECT_TRUE(r.Ok()) << deck << " 10 -> " << cells << " cells: " << r.message;
    }
    // ... and back down.
    Simulation fine = make(20);
    ASSERT_TRUE(fine.RunCycle(x).Ok());
    Simulation back = make(10);
    EXPECT_TRUE(back.RunCycle(fine.FinalState()).Ok()) << deck << " 20 -> 10 cells";
  }
}

TEST(CycleStateRemap, WarmStartOnARefinedGridFromCoarseCss) {
  // Buffer-tank deck with the Per-Component mass-balance CSS check: CSS on
  // 6/4 cells, then 12/8 cells cold and warm-started from the coarse CSS.
  auto make = [](int c1, int c2) {
    Simulation sim = TwoLayerSimulation(c1, c2, /*top_void=*/false);
    auto& css = sim.Data().settings.css;
    css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
    css.min_cycles_ = 5;
    css.consecutive_cycles_ = 3;
    sim.Data().cycle.number_of_cycles_ = 200;
    return sim;
  };
  Simulation coarse = make(6, 4);
  ASSERT_EQ(coarse.Run().status, RunStatus::ReachedCSS) << coarse.LastResult().message;
  Simulation cold = make(12, 8);
  const RunResult rc = cold.Run();
  ASSERT_EQ(rc.status, RunStatus::ReachedCSS) << rc.message;
  Simulation warm = make(12, 8);
  warm.SetInitialState(coarse.FinalState());
  const RunResult rw = warm.Run();
  ASSERT_EQ(rw.status, RunStatus::ReachedCSS) << rw.message;
  std::cout << "[ remap ] 6/4 -> 12/8 cells: coarse CSS after " << coarse.LastResult().cycles_run
            << " cycles; fine cold " << rc.cycles_run << ", warm from coarse " << rw.cycles_run << "\n";
  EXPECT_LT(rw.cycles_run, rc.cycles_run);
}

TEST(CycleStateRemap, IncompatibleStatesAreStillRefused) {
  const CycleState& x = TwoLayerState(6, 4);
  {
    Simulation other = TwoLayerSimulation(12, 8, /*top_void=*/false);
    const std::string why = other.CheckState(x);
    EXPECT_NE(why.find("void"), std::string::npos) << why;
    EXPECT_THROW(other.RemapState(x), std::invalid_argument);
  }
  {
    Simulation other = TwoLayerSimulation(12, 8);
    other.Data().bed.layers_.layers_.pop_back();
    EXPECT_NE(other.CheckState(x).find("different layers"), std::string::npos);
  }
  {
    Simulation other = TwoLayerSimulation(12, 8);
    other.Data().fluid.AddComponent("Argon");
    EXPECT_NE(other.CheckState(x).find("components"), std::string::npos);
  }
  {
    Simulation other = TwoLayerSimulation(12, 8);
    other.Data().cycle.step_list_.pop_back();
    EXPECT_NE(other.CheckState(x).find("step structure"), std::string::npos);
  }
}

}  // namespace adsim
