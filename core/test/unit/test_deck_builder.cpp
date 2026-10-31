/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// Tests for building a deck from scratch (adsim/objects/deck_builder.h): the
// builder creates structure by kind name, the parameter registry sets values.

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/objects/deck_builder.h"
#include "adsim/simulation.h"

namespace adsim {
namespace {

bool Contains(const std::vector<std::string>& v, const std::string& s) {
  return std::find(v.begin(), v.end(), s) != v.end();
}

std::string ThrowMessage(const std::function<void()>& f) {
  try {
    f();
  } catch (const std::invalid_argument& e) {
    return e.what();
  }
  return "";
}

// A one-bed CO2/N2 cycle: feed (stream in, valve out), blowdown, evacuation
// with a vacuum pump, pressurisation; a heater on the feed step; KPIs on CO2.
struct Built {
  Simulation sim;
  AdsorbentID zeolite = 0;
  size_t layer = 0, feed = 0, blowdown = 0, evacuation = 0, pressurisation = 0;
  size_t feed_in = 0, product = 0, purge_out = 0, purge_in = 0, evac_end = 0, press_end = 0;
};

void Build(Built& b) {
  UserData& d = b.sim.Data();
  build::AddComponent(d, "CarbonDioxide");
  build::AddComponent(d, "Nitrogen");

  b.zeolite = build::AddAdsorbent(d, "13X");
  build::SetIsotherm(d, b.zeolite, "extended_dual_site_langmuir_c", {"CarbonDioxide", "Nitrogen"});
  build::SetRateModel(d, b.zeolite, "CarbonDioxide", "ldf", "macropore");
  build::SetRateModel(d, b.zeolite, "Nitrogen", "ldf", "constant");

  build::SetBedShape(d, "cylinder");
  build::SetWallMaterial(d, "ss316l");
  build::SetWallBoundary(d, "zero_flux");
  build::SetThermalEquilibrium(d, "lte");

  b.layer = build::AddLayer(d, "13X");
  build::SetLayerComposition(d, b.layer, {{b.zeolite, 1.0}});
  build::SetMomentum(d, b.layer, "kozeny_carman");
  build::SetDispersion(d, b.layer, "ruthven");
  build::SetHeatDispersion(d, b.layer, "constant");
  build::SetInternalHeatTransfer(d, b.layer, "constant");
  build::SetExternalHeatTransfer(d, b.layer, "constant");

  b.feed = build::AddStep(d, "Feed");
  build::SetFlowBoundary(d, b.feed, "bottom", "inlet");
  b.feed_in = build::AddStream(d, b.feed, "bottom", "Feed inlet", "velocity", "constant");
  build::SetTemperatureProfile(d, b.feed_in, "linear");
  build::SetFlowBoundary(d, b.feed, "top", "outlet");
  b.product = build::AddValveStream(d, b.feed, "top", "Product", "constant");
  b.purge_out = build::AddStream(d, b.feed, "top", "Purge out", "molar_flow", "constant");
  build::SetEnergySource(d, b.feed, b.layer, true);

  b.blowdown = build::AddStep(d, "Blowdown");
  build::SetPressureBoundary(d, b.blowdown, "top", "outlet", "exponential");

  b.evacuation = build::AddStep(d, "Evacuation");
  b.evac_end = build::SetPressureBoundary(d, b.evacuation, "bottom", "outlet", "exponential");
  build::SetPump(d, b.evac_end, "vacuum_pump", "adiabatic", "vacuum_pump");
  build::SetFlowBoundary(d, b.evacuation, "top", "inlet");
  b.purge_in = build::AddStream(d, b.evacuation, "top", "Purge in", "molar_flow", "constant");

  b.pressurisation = build::AddStep(d, "Pressurisation");
  b.press_end = build::SetPressureBoundary(d, b.pressurisation, "bottom", "inlet", "linear");

  build::Connect(d, b.purge_out, b.purge_in);
  build::SetKPI(d, "CarbonDioxide", {b.feed_in, b.press_end}, {b.evac_end});
  b.sim.RefreshParameters();

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
      {"adsorbent[13X].kinetics[CarbonDioxide].tortuosity", "3"},
      {"adsorbent[13X].kinetics[CarbonDioxide].macropore_void_fraction", "0.35"},
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
      {"step[Feed].time", "60"},
      {"step[Feed].bottom[Feed inlet].velocity", "0.08"},
      {"step[Feed].bottom[Feed inlet].temperature_linear_start", "298"},
      {"step[Feed].bottom[Feed inlet].temperature_linear_end", "305"},
      {"step[Feed].bottom[Feed inlet].temperature_linear_time", "10"},
      {"step[Feed].bottom[Feed inlet].y[CarbonDioxide]", "0.15"},
      {"step[Feed].bottom[Feed inlet].y[Nitrogen]", "0.85"},
      {"step[Feed].top[Product].pressure", "1e5"},
      {"step[Feed].top[Product].valve.cv", "10"},
      {"step[Feed].top[Purge out].molar_flow", "1e-4"},
      {"step[Feed].energy[13X].power", "-1"},
      {"step[Blowdown].time", "20"},
      {"step[Blowdown].top.pressure_exponential_end", "2e4"},
      {"step[Evacuation].time", "60"},
      {"step[Evacuation].bottom.pressure_exponential_end", "5e3"},
      {"step[Evacuation].bottom.pump.eta_max", "0.7"},
      {"step[Evacuation].bottom.pump.discharge_pressure", "1e5"},
      {"step[Evacuation].top[Purge in].molar_flow", "1e-4"},
      {"step[Evacuation].top[Purge in].y[Nitrogen]", "1"},
      {"step[Pressurisation].time", "20"},
      {"step[Pressurisation].bottom.pressure_linear_end", "1e5"},
      {"step[Pressurisation].bottom.y[CarbonDioxide]", "0.15"},
      {"step[Pressurisation].bottom.y[Nitrogen]", "0.85"},
  };
  for (const auto& [path, text] : values) b.sim.SetParameterText(path, text);
}

}  // namespace

TEST(DeckBuilder, KindListsMatchTheEngineEnums) {
  EXPECT_TRUE(Contains(build::IsothermKinds(), "inert"));
  EXPECT_TRUE(Contains(build::IsothermKinds(), "extended_dual_site_langmuir_c"));
  EXPECT_EQ(build::FlowKinds(), (std::vector<std::string>{"velocity", "volumetric_flow", "molar_flow", "mass_flow"}));
  EXPECT_EQ(build::ProfileKinds(), (std::vector<std::string>{"constant", "linear", "exponential", "log_linear", "double_exponential"}));
  EXPECT_EQ(build::BedShapeKinds(), (std::vector<std::string>{"cylinder", "square", "rectangle"}));
  EXPECT_EQ(build::MomentumKinds(), (std::vector<std::string>{"kozeny_carman", "ergun"}));
  EXPECT_EQ(build::WallMaterialKinds().front(), "custom");
}

TEST(DeckBuilder, UnknownKindsAndObjectsThrowWithTheValidChoices) {
  Simulation sim;
  UserData& d = sim.Data();
  build::AddComponent(d, "Nitrogen");
  build::AddComponent(d, "Nitrogen");  // adding a component twice is a no-op
  EXPECT_EQ(d.fluid.components_.size(), 1u);
  EXPECT_THROW(build::AddComponent(d, "Unobtainium"), std::invalid_argument);
  const AdsorbentID a = build::AddAdsorbent(d, "5A");
  EXPECT_NE(ThrowMessage([&] { build::SetIsotherm(d, a, "nope", {"Nitrogen"}); }).find("henry"), std::string::npos);
  EXPECT_THROW(build::SetIsotherm(d, a, "henry", {"Argon"}), std::invalid_argument);
  EXPECT_THROW(build::SetIsotherm(d, a + 7, "henry", {"Nitrogen"}), std::invalid_argument);
  EXPECT_NE(ThrowMessage([&] { build::SetRateModel(d, a, "Nitrogen", "fick", "constant"); }).find("ldf2"),
            std::string::npos);
  EXPECT_NE(ThrowMessage([&] { build::SetBedShape(d, "sphere"); }).find("cylinder"), std::string::npos);
  const size_t layer = build::AddLayer(d, "bed");
  EXPECT_THROW(build::SetMomentum(d, layer, "darcy"), std::invalid_argument);
  EXPECT_THROW(build::SetLayerComposition(d, layer, {{a + 7, 1.0}}), std::invalid_argument);
  const size_t step = build::AddStep(d, "Feed");
  EXPECT_THROW(build::AddStep(d, "Feed"), std::invalid_argument);
  EXPECT_THROW(build::SetPressureBoundary(d, step, "middle", "inlet", "constant"), std::invalid_argument);
  EXPECT_THROW(build::SetPressureBoundary(d, step, "top", "sideways", "constant"), std::invalid_argument);
  // Streams need a flow boundary first.
  EXPECT_THROW(build::AddStream(d, step, "bottom", "feed", "velocity", "constant"), std::invalid_argument);
  build::SetFlowBoundary(d, step, "bottom", "inlet");
  EXPECT_THROW(build::AddStream(d, step, "bottom", "feed", "litres", "constant"), std::invalid_argument);
  const size_t valve = build::AddValveStream(d, step, "bottom", "valve", "constant");
  EXPECT_THROW(build::SetPump(d, valve, "compressor", "isothermal", "constant"), std::invalid_argument);
  EXPECT_THROW(build::SetKPI(d, "Argon", {valve}, {valve}), std::invalid_argument);
  EXPECT_THROW(build::Connect(d, valve, 999), std::invalid_argument);
  EXPECT_THROW(build::AddVariable(d, "T", "temperature", "", ""), std::invalid_argument);    // reserved
  EXPECT_THROW(build::AddVariable(d, "2x", "time", "", ""), std::invalid_argument);          // not a name
  EXPECT_THROW(build::AddVariable(d, "t_feed", "minutes", "", ""), std::invalid_argument);   // category
}

TEST(DeckBuilder, UserVariablesDriveOtherValues) {
  Built b;
  Build(b);
  UserData& d = b.sim.Data();
  build::AddVariable(d, "t_feed", "time", "", "feed time");
  build::AddVariable(d, "t_half", "time", "t_feed / 2", "");
  b.sim.RefreshParameters();
  b.sim.SetParameter("variable[t_feed]", 90.0);
  b.sim.SetParameterText("step[Blowdown].time", "t_half");
  EXPECT_EQ(b.sim.GetParameter("step[Blowdown].time"), 45.0);
  EXPECT_EQ(build::VariableCategoryKind(d.variable_registry.Find("t_feed")->unit_category), "time");
  EXPECT_TRUE(Contains(build::VariableCategoryKinds(), "custom"));
}

TEST(DeckBuilder, BuildsADeckThatPassesTheChecksAndRuns) {
  Built b;
  Build(b);
  std::string error;
  ASSERT_TRUE(b.sim.Check(error)) << error;

  const UserData& d = b.sim.Data();
  EXPECT_EQ(d.fluid.components_, (std::vector<std::string>{"CarbonDioxide", "Nitrogen"}));
  ASSERT_EQ(d.cycle.step_list_.size(), 4u);
  EXPECT_EQ(b.sim.GetParameter("adsorbent[13X].isotherm[CarbonDioxide].q1"), 40980.0);
  EXPECT_EQ(b.sim.GetParameter("layer[13X].cells"), 20.0);
  EXPECT_EQ(b.sim.GetParameter("step[Feed].energy[13X].power"), -1.0);
  EXPECT_EQ(b.sim.GetParameter("step[Evacuation].bottom.pump.eta_max"), 0.7);
  EXPECT_EQ(b.sim.GetParameter("cycle.number_of_cycles"), 1.0);

  RunOptions options;
  options.cycles = 1;
  const RunResult r = b.sim.Run(options);
  ASSERT_TRUE(r.Ok()) << r.message;
  const auto kpis = b.sim.KPIs(-1);
  ASSERT_EQ(kpis.size(), 1u);
  EXPECT_EQ(kpis[0].component, "CarbonDioxide");
  EXPECT_TRUE(kpis[0].valid);
}

TEST(DeckBuilder, ReplacingAnEndDropsItsStreamsAndConnections) {
  Built b;
  Build(b);
  UserData& d = b.sim.Data();
  build::SetClosed(d, b.feed, "top");
  b.sim.RefreshParameters();
  EXPECT_EQ(b.sim.Parameters().Find("step[Feed].top[Purge out].molar_flow"), nullptr);
  std::string error;
  EXPECT_TRUE(b.sim.Check(error)) << error;  // the purge inlet is no longer connected to anything
  // The KPI product was the evacuation boundary; replacing that end drops the ref.
  ASSERT_EQ(d.kpi_config.FindByComponent("CarbonDioxide")->product_endpoints.size(), 1u);
  build::SetFlowBoundary(d, b.evacuation, "bottom", "outlet");
  EXPECT_TRUE(d.kpi_config.FindByComponent("CarbonDioxide")->product_endpoints.empty());
  EXPECT_EQ(d.kpi_config.FindByComponent("CarbonDioxide")->feed_endpoints.size(), 2u);
  // An isotherm reassignment moves components between isotherm objects.
  build::SetIsotherm(d, b.zeolite, "henry", {"Nitrogen"});
  const auto& map = d.adsorbent_directory.Get(b.zeolite).adsorption_source_.IsothermMap();
  EXPECT_NE(map.at("Nitrogen"), map.at("CarbonDioxide"));
  b.sim.RefreshParameters();
  EXPECT_NE(b.sim.Parameters().Find("adsorbent[13X].isotherm[Nitrogen].k"), nullptr);
  EXPECT_EQ(b.sim.Parameters().Find("adsorbent[13X].isotherm[Nitrogen].m1"), nullptr);
}

// A pump-driven stream takes its flow from a shaft power: a vacuum pump on an
// outlet, a compressor on an inlet; its power model and pump are in the registry.
TEST(DeckBuilder, PumpStreamsAreDrivenByAPower) {
  Built b;
  Build(b);
  UserData& d = b.sim.Data();
  build::SetFlowBoundary(d, b.evacuation, "bottom", "outlet");
  const size_t vac = build::AddPumpStream(d, b.evacuation, "bottom", "Vacuum pump", "constant", "isothermal",
                                          "vacuum_pump");
  build::SetKPI(d, "CarbonDioxide", {b.feed_in, b.press_end}, {vac});
  b.sim.RefreshParameters();
  b.sim.SetParameterText("step[Evacuation].bottom[Vacuum pump].power", "40");
  b.sim.SetParameterText("step[Evacuation].bottom[Vacuum pump].pump.eta_max", "0.7");
  b.sim.SetParameterText("step[Evacuation].bottom[Vacuum pump].pump.discharge_pressure", "1e5");
  std::string error;
  ASSERT_TRUE(b.sim.Check(error)) << error;
  RunOptions options;
  options.cycles = 1;
  const RunResult r = b.sim.Run(options);
  ASSERT_TRUE(r.Ok()) << r.message;
  EXPECT_GT(b.sim.KPIs(-1).at(0).energy_consumption, 0.0);

  // A power polynomial of the boundary pressure, and the pump cannot be removed.
  build::SetFlowBoundary(d, b.pressurisation, "bottom", "inlet");
  const size_t comp = build::AddPumpStream(d, b.pressurisation, "bottom", "Compressor", "polynomial", "adiabatic",
                                           "constant");
  b.sim.RefreshParameters();
  EXPECT_NE(b.sim.Parameters().Find("endpoint[" + std::to_string(comp) + "].power_poly_f"), nullptr);
  EXPECT_NE(b.sim.Parameters().Find("endpoint[" + std::to_string(comp) + "].pump.suction_pressure"), nullptr);
  EXPECT_NE(ThrowMessage([&] { build::SetPump(d, comp, "none", "isothermal", "constant"); }).find("always has a pump"),
            std::string::npos);
}

// A buffer tank feeds an inlet with a share of another endpoint's product.
TEST(DeckBuilder, BufferTanksDrawFromAnotherEndpoint) {
  Built b;
  Build(b);
  UserData& d = b.sim.Data();
  build::SetBufferTank(d, b.press_end, b.product);
  b.sim.RefreshParameters();
  b.sim.SetParameterText("step[Pressurisation].bottom.buffer_tank.percentage", "40 [%]");
  EXPECT_DOUBLE_EQ(b.sim.GetParameter("step[Pressurisation].bottom.buffer_tank.percentage"), 40.0);
  std::string error;
  ASSERT_TRUE(b.sim.Check(error)) << error;
  RunOptions options;
  options.cycles = 2;
  ASSERT_TRUE(b.sim.Run(options).Ok());

  EXPECT_NE(ThrowMessage([&] { build::SetBufferTank(d, b.evac_end, b.product); }).find("only an inlet"),
            std::string::npos);
  EXPECT_NE(ThrowMessage([&] { build::SetBufferTank(d, b.purge_in, b.product); }).find("connected"),
            std::string::npos);
  build::SetBufferTank(d, b.press_end, 0);
  b.sim.RefreshParameters();
  EXPECT_EQ(b.sim.Parameters().Find("step[Pressurisation].bottom.buffer_tank.percentage"), nullptr);
}

}  // namespace adsim
