/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// A vacuum pump stream's ultimate (blank-off) pressure, and the step-end
// pressure output (step_end_pressure.h) that a study or the optimiser reads.

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/simulation.h"
#include "adsim/solver/optimisation.h"
#include "adsim/solver/sensitivity.h"
#include "adsim/solver/step_end_pressure.h"

namespace adsim {
namespace {

const std::filesystem::path kDecks = std::filesystem::current_path() / ".." / "examples" / "decks";

// The Perez 2019 VSA column (feed, blowdown, evacuation, pressurisation),
// coarse (20 cells), one cycle, with the evacuation drawn through the bottom
// by a constant-power vacuum pump instead of the deck's pressure profile. It
// opens after the blowdown, below the 1 bar discharge pressure: a constant-
// power pump opened at a pressure ratio of about 1 draws a near-unbounded
// flow (issue #250).
struct Column {
  Simulation sim;
  size_t evacuation = 0;  // step ID
  int evac_index = 0;     // its index in the cycle
  size_t pump = 0;        // endpoint ID
};

void Build(Column& c) {
  c.sim.Load((kDecks / "Perez2019_A.json").string());
  UserData& d = c.sim.Data();
  const auto& steps = d.cycle.step_list_;
  for (size_t s = 0; s < steps.size(); ++s)
    if (steps[s].step_name_ == "Evacuation") {
      c.evacuation = steps[s].id_;
      c.evac_index = static_cast<int>(s);
    }
  build::SetFlowBoundary(d, c.evacuation, "bottom", "outlet");
  c.pump = build::AddPumpStream(d, c.evacuation, "bottom", "Vacuum pump", "constant", "isothermal", "vacuum_pump");
  c.sim.RefreshParameters();
  const std::vector<std::pair<std::string, std::string>> values = {
      {"layer[1].cells", "20"},
      {"css.method", "0"},
      {"step[Evacuation].time", "300"},
      {"step[Evacuation].bottom[Vacuum pump].power", "100"},
      {"step[Evacuation].bottom[Vacuum pump].pump.eta_max", "0.7"},
      {"step[Evacuation].bottom[Vacuum pump].pump.discharge_pressure", "1e5"},
  };
  for (const auto& [path, text] : values) c.sim.SetParameterText(path, text);
}

// Runs one cycle; returns the evacuation's end pressure at the pumped end.
double RunColumn(Column& c) {
  std::string error;
  EXPECT_TRUE(c.sim.Check(error)) << error;
  RunOptions options;
  options.cycles = 1;
  const RunResult r = c.sim.Run(options);
  EXPECT_TRUE(r.Ok()) << r.message;
  const std::vector<StepEndPressure> ends = c.sim.StepEndPressures(-1);
  EXPECT_EQ(ends.size(), 4u);
  return static_cast<size_t>(c.evac_index) < ends.size() ? ends[static_cast<size_t>(c.evac_index)].bottom
                                                          : std::nan("");
}

}  // namespace

// The factor is 1 without an ultimate pressure, 0 at and below it, the speed
// law 1 - p_ult/P well above it, and C1 where the two meet.
TEST(PumpUltimatePressure, FactorIsSmoothAndVanishesAtTheUltimatePressure) {
  const double p_ult = 1e4;
  EXPECT_EQ(UltimatePressureFactor(5e3, 0.0), 1.0);
  EXPECT_EQ(UltimatePressureFactor(5e3, p_ult), 0.0);
  EXPECT_EQ(UltimatePressureFactor(p_ult, p_ult), 0.0);
  EXPECT_DOUBLE_EQ(UltimatePressureFactor(5e4, p_ult), 1.0 - p_ult / 5e4);

  // Value and slope continuous across the blend's ends (u = 0 and u = d).
  const double d = kUltimatePressureBlend;
  auto f_of_u = [&](double u) { return UltimatePressureFactor(p_ult / (1.0 - u), p_ult); };
  const double h = 1e-7;
  for (double u : {0.0, d}) {
    const double below = f_of_u(u - h), at = f_of_u(u), above = f_of_u(u + h);
    EXPECT_NEAR(below, at, 1e-6);
    EXPECT_NEAR(above, at, 1e-6);
    const double slope_below = (at - below) / h, slope_above = (above - at) / h;
    EXPECT_NEAR(slope_below, slope_above, 1e-4) << "u = " << u;
  }
  // Non-decreasing in the suction pressure.
  double last = 0.0;
  for (double P = 0.9 * p_ult; P < 3 * p_ult; P += 10.0) {
    const double f = UltimatePressureFactor(P, p_ult);
    EXPECT_GE(f, last - 1e-15);
    last = f;
  }
}

// The leaf is listed on a vacuum pump stream, defaults to 0 (off), and a
// negative value fails the checks.
TEST(PumpUltimatePressure, IsARegistryLeafThatDefaultsToOff) {
  Column c;
  Build(c);
  const std::string path = "endpoint[" + std::to_string(c.pump) + "].pump.ultimate_pressure";
  ASSERT_NE(c.sim.Parameters().Find(path), nullptr);
  EXPECT_EQ(c.sim.GetParameter(path), 0.0);
  EXPECT_NE(c.sim.Parameters().Find("step[Evacuation].bottom[Vacuum pump].pump.ultimate_pressure"), nullptr);
  c.sim.SetParameterText(path, "-1");
  std::string error;
  EXPECT_FALSE(c.sim.Check(error));
  EXPECT_NE(error.find("ultimate pressure"), std::string::npos) << error;
}

// The pump draws the column down to its ultimate pressure and holds it there,
// where without one it goes far below.
TEST(PumpUltimatePressure, IsReachedAndHeld) {
  Column free;
  Build(free);
  const double p_free = RunColumn(free);

  const double p_ult = 4e3;  // below the blowdown's end (~7.7 kPa)
  Column floor;
  Build(floor);
  floor.sim.SetParameterText("step[Evacuation].bottom[Vacuum pump].pump.ultimate_pressure", std::to_string(p_ult));
  const double p_end = RunColumn(floor);

  EXPECT_LT(p_free, 0.5 * p_ult);
  EXPECT_GE(p_end, p_ult);
  EXPECT_LT(p_end, 1.02 * p_ult);

  // Never below the ultimate pressure at the pumped end during the
  // evacuation, and settled there: its last fifth moves the pressure by less
  // than 1 % of p_ult.
  const std::vector<StepEndPressure> ends = floor.sim.StepEndPressures(-1);
  const size_t e = static_cast<size_t>(floor.evac_index);
  ASSERT_GT(e, 0u);
  const double t0 = ends[e - 1].time, t1 = ends[e].time;
  const TimeSeries s = floor.sim.Series("pressure", {}, SeriesLocation::Bottom, 0.0, -1);
  double p_at_80 = 0.0;
  for (size_t i = 0; i < s.time.size(); ++i) {
    if (s.time[i] <= t0 || s.time[i] > t1) continue;
    EXPECT_GE(s.values[i], p_ult * (1.0 - 1e-9)) << "t = " << s.time[i];
    if (s.time[i] <= t0 + 0.8 * (t1 - t0)) p_at_80 = s.values[i];
  }
  EXPECT_LT(p_at_80 - p_end, 0.01 * p_ult);
}

// The step-end pressures equal the last logged values of the bed-end series
// and of the axial profile.
TEST(StepEndPressure, EqualsTheLastValueOfTheSeriesAndTheProfile) {
  Column c;
  Build(c);
  RunColumn(c);
  const std::vector<StepEndPressure> ends = c.sim.StepEndPressures(-1);
  ASSERT_EQ(ends.size(), 4u);
  const auto& steps = c.sim.Data().cycle.step_list_;
  for (size_t s = 0; s < ends.size(); ++s) {
    EXPECT_EQ(ends[s].step, static_cast<int>(s));
    EXPECT_EQ(ends[s].step_name, steps[s].step_name_);
    EXPECT_EQ(ends[s].step_id, steps[s].id_);
  }

  const TimeSeries bottom = c.sim.Series("pressure", {}, SeriesLocation::Bottom, 0.0, -1);
  const TimeSeries top = c.sim.Series("pressure", {}, SeriesLocation::Top, 0.0, -1);
  ASSERT_FALSE(bottom.time.empty());
  // The last step's end is the last row of the cycle.
  EXPECT_EQ(ends.back().bottom, bottom.values.back());
  EXPECT_EQ(ends.back().top, top.values.back());
  EXPECT_EQ(ends.back().time, bottom.time.back());
  // Every other step's end is its last row: the last row at its end time
  // (the next step's first row repeats that time).
  for (size_t s = 0; s + 1 < ends.size(); ++s) {
    size_t last = 0;
    for (size_t i = 0; i < bottom.time.size(); ++i)
      if (bottom.time[i] < ends[s].time || (bottom.time[i] == ends[s].time && last < i && bottom.time[last] < ends[s].time)) last = i;
    EXPECT_EQ(ends[s].bottom, bottom.values[last]) << ends[s].step_name;
    EXPECT_EQ(ends[s].top, top.values[last]) << ends[s].step_name;
  }

  // Bed statistics from each step's profile, last row.
  for (int step = 0; step < 4; ++step) {
    const AxialProfile p = c.sim.Profile("pressure", {}, -1, step);
    ASSERT_GT(p.rows, 0u);
    const std::vector<double> row(p.values.end() - static_cast<std::ptrdiff_t>(p.cols), p.values.end());
    const auto [lo, hi] = std::minmax_element(row.begin(), row.end());
    EXPECT_DOUBLE_EQ(ends[step].bed_min, *lo);
    EXPECT_DOUBLE_EQ(ends[step].bed_max, *hi);
    EXPECT_GE(ends[step].bed_mean, *lo);
    EXPECT_LE(ends[step].bed_mean, *hi);
  }
  // The pumped end is the lowest point at the end of the evacuation.
  EXPECT_LE(ends[static_cast<size_t>(c.evac_index)].bottom, ends[static_cast<size_t>(c.evac_index)].bed_min);
}

// Step references resolve by deck ID or name; the outputs and the optimiser
// terms read the record.
TEST(StepEndPressure, StepReferencesAndOutputs) {
  Column c;
  Build(c);
  const auto& steps = c.sim.Data().cycle.step_list_;
  EXPECT_TRUE(IsStepReference("step[Evacuation]"));
  EXPECT_FALSE(IsStepReference("CarbonDioxide"));
  EXPECT_FALSE(IsStepReference("step[]"));
  EXPECT_EQ(ResolveStepReference(steps, "step[Evacuation]"), c.evac_index);
  EXPECT_EQ(ResolveStepReference(steps, "step[" + std::to_string(c.evacuation) + "]"), c.evac_index);
  EXPECT_THROW(ResolveStepReference(steps, "step[Hold]"), std::invalid_argument);

  const sensitivity::Output o = sensitivity::ParseOutput("step[Evacuation].end_pressure_bottom");
  EXPECT_EQ(o.component, "step[Evacuation]");
  EXPECT_EQ(o.kpi, "end_pressure_bottom");
  EXPECT_THROW(sensitivity::ParseOutput("step[Evacuation].purity"), std::invalid_argument);

  // A sample reads the step output from the run.
  const std::string deck = c.sim.ToJson();
  sensitivity::RunPolicy policy;
  policy.require_css = false;
  policy.cycles = 1;
  const std::vector<sensitivity::Output> outputs = {o, sensitivity::ParseOutput("step[Feed].end_pressure_mean")};
  const std::vector<sensitivity::Input> inputs = {
      sensitivity::Uniform("step[Evacuation].bottom[Vacuum pump].power", 80.0, 120.0)};
  sensitivity::CheckStudy(deck, inputs, outputs, policy);
  EXPECT_THROW(sensitivity::CheckStudy(deck, inputs, {sensitivity::ParseOutput("step[Hold].end_pressure_top")}, policy),
               std::invalid_argument);
  const sensitivity::SampleOutcome s = sensitivity::EvaluateSample(deck, inputs, {}, outputs, policy);
  ASSERT_TRUE(s.accepted) << s.message;
  const double p_end = RunColumn(c);
  EXPECT_NEAR(s.outputs[0], p_end, 1e-6 * p_end);
  EXPECT_TRUE(std::isfinite(s.outputs[1]));

  // An optimiser constraint on a step-end pressure: the cost adds its
  // violation, and a term mixing a step with a component KPI is rejected.
  OptimizationCostFunction cost;
  cost.AddObjective("step[Evacuation]", KPIType::END_PRESSURE_BOTTOM, 1.0, ObjectiveDirection::MINIMIZE);
  cost.AddConstraint("step[Evacuation]", KPIType::END_PRESSURE_BOTTOM, 1e4, 2.0, ConstraintSense::GREATER_OR_EQUAL);
  PerComponentKPIValues values;
  values["step[Evacuation]"][KPIType::END_PRESSURE_BOTTOM] = 8e3;
  EXPECT_DOUBLE_EQ(cost.ComputeCost(values), 8e3 + 2.0 * 2e3);
  const std::vector<OptimizationVariable> variables = {{"step[Evacuation].bottom[Vacuum pump].power", 80.0, 120.0}};
  SimulationOptimizerSettings settings;
  settings.cycles = 1;
  CheckOptimization(deck, variables, cost, settings);
  OptimizationCostFunction mixed;
  mixed.AddObjective("step[Evacuation]", KPIType::PURITY, 1.0, ObjectiveDirection::MAXIMIZE);
  EXPECT_THROW(CheckOptimization(deck, variables, mixed, settings), std::invalid_argument);
  OptimizationCostFunction wrong;
  wrong.AddObjective("CarbonDioxide", KPIType::END_PRESSURE_TOP, 1.0, ObjectiveDirection::MINIMIZE);
  EXPECT_THROW(CheckOptimization(deck, variables, wrong, settings), std::invalid_argument);
}

}  // namespace adsim
