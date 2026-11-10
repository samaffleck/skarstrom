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

// A 0.3 m 13X column filled with N2 at 1 bar, closed at the top, evacuated
// from the bottom by a constant-power vacuum pump (20 cells, short steps).
struct Column {
  Simulation sim;
  size_t evacuation = 0, hold = 0, pump = 0;
};

void Build(Column& c, double evacuation_time = 60.0) {
  UserData& d = c.sim.Data();
  build::AddComponent(d, "CarbonDioxide");
  build::AddComponent(d, "Nitrogen");
  const AdsorbentID zeolite = build::AddAdsorbent(d, "13X");
  build::SetIsotherm(d, zeolite, "extended_dual_site_langmuir_c", {"CarbonDioxide", "Nitrogen"});
  build::SetRateModel(d, zeolite, "CarbonDioxide", "ldf", "constant");
  build::SetRateModel(d, zeolite, "Nitrogen", "ldf", "constant");
  build::SetBedShape(d, "cylinder");
  build::SetWallMaterial(d, "ss316l");
  build::SetWallBoundary(d, "zero_flux");
  build::SetThermalEquilibrium(d, "lte");
  const size_t layer = build::AddLayer(d, "13X");
  build::SetLayerComposition(d, layer, {{zeolite, 1.0}});
  build::SetMomentum(d, layer, "kozeny_carman");
  build::SetDispersion(d, layer, "ruthven");
  build::SetHeatDispersion(d, layer, "constant");
  build::SetInternalHeatTransfer(d, layer, "constant");
  build::SetExternalHeatTransfer(d, layer, "constant");

  c.evacuation = build::AddStep(d, "Evacuation");
  build::SetFlowBoundary(d, c.evacuation, "bottom", "outlet");
  c.pump = build::AddPumpStream(d, c.evacuation, "bottom", "Vacuum pump", "constant", "isothermal", "vacuum_pump");
  // A closed hold step after it: its end pressure is the evacuation's.
  c.hold = build::AddStep(d, "Hold");
  c.sim.RefreshParameters();

  const std::vector<std::pair<std::string, std::string>> values = {
      {"bed.diameter", "0.032"},
      {"wall.thickness", "0.0031"},
      {"adsorbent[13X].solid_density", "1138"},
      {"adsorbent[13X].heat_capacity", "900"},
      {"adsorbent[13X].thermal_conductivity", "0.63"},
      {"adsorbent[13X].particle_void_fraction", "0.35"},
      {"adsorbent[13X].isotherm[CarbonDioxide].m1", "2.768"},
      {"adsorbent[13X].isotherm[CarbonDioxide].b1", "3.934e-7"},
      {"adsorbent[13X].isotherm[CarbonDioxide].q1", "40980"},
      {"adsorbent[13X].isotherm[CarbonDioxide].m2", "2.758"},
      {"adsorbent[13X].isotherm[CarbonDioxide].b2", "1.256e-7"},
      {"adsorbent[13X].isotherm[CarbonDioxide].q2", "34190"},
      {"adsorbent[13X].isotherm[Nitrogen].m1", "2.768"},
      {"adsorbent[13X].isotherm[Nitrogen].b1", "3.134e-6"},
      {"adsorbent[13X].isotherm[Nitrogen].q1", "15710"},
      {"adsorbent[13X].isotherm[Nitrogen].m2", "2.758"},
      {"adsorbent[13X].isotherm[Nitrogen].b2", "3.134e-6"},
      {"adsorbent[13X].isotherm[Nitrogen].q2", "15720"},
      {"adsorbent[13X].kinetics[CarbonDioxide].k", "0.5"},
      {"adsorbent[13X].kinetics[Nitrogen].k", "0.5"},
      {"adsorbent[13X].heat_of_adsorption[CarbonDioxide]", "30850"},
      {"adsorbent[13X].heat_of_adsorption[Nitrogen]", "15930"},
      {"layer[13X].length", "0.3"},
      {"layer[13X].cells", "20"},
      {"layer[13X].bed_void_fraction", "0.41"},
      {"layer[13X].particle_diameter", "0.0015"},
      {"layer[13X].molecular_diffusion_coefficient", "1.6e-5"},
      {"initial.pressure", "1e5"},
      {"initial.temperature", "298"},
      {"initial.y[Nitrogen]", "1"},
      {"initial.y[CarbonDioxide]", "0"},
      {"step[Evacuation].time", std::to_string(evacuation_time)},
      {"step[Evacuation].bottom[Vacuum pump].power", "5"},
      {"step[Evacuation].bottom[Vacuum pump].pump.eta_max", "0.7"},
      {"step[Evacuation].bottom[Vacuum pump].pump.discharge_pressure", "1e5"},
      {"step[Hold].time", "5"},
  };
  for (const auto& [path, text] : values) c.sim.SetParameterText(path, text);
}

double Run(Column& c) {
  std::string error;
  EXPECT_TRUE(c.sim.Check(error)) << error;
  RunOptions options;
  options.cycles = 1;
  const RunResult r = c.sim.Run(options);
  EXPECT_TRUE(r.Ok()) << r.message;
  const std::vector<StepEndPressure> ends = c.sim.StepEndPressures(-1);
  EXPECT_EQ(ends.size(), 2u);
  return ends.empty() ? std::nan("") : ends.front().bottom;
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
  const double p_free = Run(free);

  const double p_ult = 3e4;
  Column floor;
  Build(floor);
  floor.sim.SetParameterText("step[Evacuation].bottom[Vacuum pump].pump.ultimate_pressure", std::to_string(p_ult));
  const double p_end = Run(floor);

  EXPECT_LT(p_free, 0.5 * p_ult);
  EXPECT_GE(p_end, p_ult);
  EXPECT_LT(p_end, 1.02 * p_ult);

  // Never below the ultimate pressure at the pumped end, and settled there:
  // the last fifth of the step moves the pressure by less than 1 %.
  const TimeSeries s = floor.sim.Series("pressure", {}, SeriesLocation::Bottom, 0.0, -1);
  const double t_end = floor.sim.StepEndPressures(-1).front().time;
  double p_at_80 = 0.0;
  for (size_t i = 0; i < s.time.size() && s.time[i] <= t_end; ++i) {
    EXPECT_GE(s.values[i], p_ult * (1.0 - 1e-9)) << "t = " << s.time[i];
    if (s.time[i] <= 0.8 * t_end) p_at_80 = s.values[i];
  }
  EXPECT_LT(p_at_80 - p_end, 0.01 * p_ult);
}

// The step-end pressures equal the last logged values of the bed-end series
// and of the axial profile.
TEST(StepEndPressure, EqualsTheLastValueOfTheSeriesAndTheProfile) {
  Column c;
  Build(c);
  Run(c);
  const std::vector<StepEndPressure> ends = c.sim.StepEndPressures(-1);
  ASSERT_EQ(ends.size(), 2u);
  EXPECT_EQ(ends[0].step_name, "Evacuation");
  EXPECT_EQ(ends[1].step_name, "Hold");

  const TimeSeries bottom = c.sim.Series("pressure", {}, SeriesLocation::Bottom, 0.0, -1);
  const TimeSeries top = c.sim.Series("pressure", {}, SeriesLocation::Top, 0.0, -1);
  ASSERT_FALSE(bottom.time.empty());
  // The last step's end is the last row of the cycle.
  EXPECT_EQ(ends[1].bottom, bottom.values.back());
  EXPECT_EQ(ends[1].top, top.values.back());
  EXPECT_EQ(ends[1].time, bottom.time.back());
  // The evacuation's end is the last row at its end time.
  size_t last_evac = 0;
  for (size_t i = 0; i < bottom.time.size(); ++i)
    if (bottom.time[i] <= ends[0].time) last_evac = i;
  EXPECT_EQ(ends[0].bottom, bottom.values[last_evac]);
  EXPECT_NEAR(ends[0].time, 60.0, 1e-9);

  // Bed statistics from the profile's last row.
  for (int step = 0; step < 2; ++step) {
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
  EXPECT_LE(ends[0].bottom, ends[0].bed_min);
}

// Step references resolve by deck ID or name; the outputs and the optimiser
// terms read the record.
TEST(StepEndPressure, StepReferencesAndOutputs) {
  Column c;
  Build(c);
  const auto& steps = c.sim.Data().cycle.step_list_;
  EXPECT_TRUE(IsStepReference("step[Evacuation]"));
  EXPECT_FALSE(IsStepReference("CarbonDioxide"));
  EXPECT_EQ(ResolveStepReference(steps, "step[Hold]"), 1);
  EXPECT_EQ(ResolveStepReference(steps, "step[" + std::to_string(steps[1].id_) + "]"), 1);
  EXPECT_THROW(ResolveStepReference(steps, "step[Feed]"), std::invalid_argument);

  const sensitivity::Output o = sensitivity::ParseOutput("step[Evacuation].end_pressure_bottom");
  EXPECT_EQ(o.component, "step[Evacuation]");
  EXPECT_EQ(o.kpi, "end_pressure_bottom");
  EXPECT_THROW(sensitivity::ParseOutput("step[Evacuation].purity"), std::invalid_argument);

  // A sample reads the step output from the run.
  const std::string deck = c.sim.ToJson();
  sensitivity::RunPolicy policy;
  policy.require_css = false;
  policy.cycles = 1;
  const std::vector<sensitivity::Output> outputs = {o, sensitivity::ParseOutput("step[Hold].end_pressure_mean")};
  sensitivity::CheckStudy(deck, {}, outputs, policy);
  EXPECT_THROW(sensitivity::CheckStudy(deck, {}, {sensitivity::ParseOutput("step[Feed].end_pressure_top")}, policy),
               std::invalid_argument);
  const sensitivity::SampleOutcome s = sensitivity::EvaluateSample(deck, {}, {}, outputs, policy);
  ASSERT_TRUE(s.accepted) << s.message;
  const double p_end = Run(c);
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
  const std::vector<OptimizationVariable> variables = {
      {"step[Evacuation].bottom[Vacuum pump].power", 2.0, 8.0}};
  SimulationOptimizerSettings settings;
  settings.cycles = 1;
  CheckOptimization(deck, variables, cost, settings);
  OptimizationCostFunction mixed;
  mixed.AddObjective("step[Evacuation]", KPIType::PURITY, 1.0, ObjectiveDirection::MAXIMIZE);
  EXPECT_THROW(CheckOptimization(deck, variables, mixed, settings), std::invalid_argument);
  OptimizationCostFunction wrong;
  wrong.AddObjective("Nitrogen", KPIType::END_PRESSURE_TOP, 1.0, ObjectiveDirection::MINIMIZE);
  EXPECT_THROW(CheckOptimization(deck, variables, wrong, settings), std::invalid_argument);
}

}  // namespace adsim
