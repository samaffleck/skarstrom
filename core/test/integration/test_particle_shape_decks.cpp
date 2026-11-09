/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// The particle shape (MixedAdsorbentLayer v7) moves no result of a deck that
// does not select a cylinder. For every example deck: every layer loads as a
// sphere, and a run is bitwise identical to the same deck written as a v6
// archive (no Particle_Shape field, the layer class at version 6). Short runs
// (2 cycles, 10 cells per layer, CSS off): the assertion is identity, which a
// short run shows as well as a long one.

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <cereal/external/rapidjson/document.h>
#include <cereal/external/rapidjson/stringbuffer.h>
#include <cereal/external/rapidjson/writer.h>
#include <gtest/gtest.h>

#include "adsim/objects/porous_media.h"
#include "adsim/simulation.h"

namespace adsim {
namespace {

const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";

bool SameBits(const std::vector<double>& a, const std::vector<double>& b) {
  return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
}
bool SameBits(const std::vector<std::vector<double>>& a, const std::vector<std::vector<double>>& b) {
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); ++i)
    if (!SameBits(a[i], b[i])) return false;
  return true;
}
bool SameBits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

void ExpectSameLogs(const EquationLog& a, const EquationLog& b) {
  const auto& la = a.GetSimulationLog();
  const auto& lb = b.GetSimulationLog();
  ASSERT_EQ(la.size(), lb.size());
  for (const auto& [cycle, steps] : la) {
    ASSERT_TRUE(lb.count(cycle));
    for (const auto& [step, log] : steps) {
      SCOPED_TRACE("cycle " + std::to_string(cycle) + " step " + std::to_string(step));
      const Log& other = lb.at(cycle).at(step);
      ASSERT_TRUE(SameBits(log.step_data.times, other.step_data.times));
      ASSERT_EQ(log.variables.size(), other.variables.size());
      for (const auto& [type, v] : log.variables)
        EXPECT_TRUE(SameBits(v.values, other.variables.at(type).values)) << v.variable_name;
      ASSERT_EQ(log.component_variables.size(), other.component_variables.size());
      for (const auto& [comp, vars] : log.component_variables) {
        ASSERT_EQ(vars.size(), other.component_variables.at(comp).size()) << comp;
        for (const auto& [type, v] : vars)
          EXPECT_TRUE(SameBits(v.values, other.component_variables.at(comp).at(type).values)) << comp << " " << v.variable_name;
      }
    }
  }
}

// Runs both and compares the integrator path, the KPIs, the logs and the
// final state.
void ExpectSameRun(Simulation& generic, Simulation& legacy) {
  RunOptions options;
  options.on_message = [](LogLevel, const std::string&) {};
  const RunResult a = generic.Run(options);
  const RunResult b = legacy.Run(options);
  ASSERT_TRUE(a.Ok()) << a.message;
  ASSERT_TRUE(b.Ok()) << b.message;
  EXPECT_EQ(a.cycles_run, b.cycles_run);
  EXPECT_EQ(a.stats.num_steps, b.stats.num_steps);
  EXPECT_EQ(a.stats.num_res_evals, b.stats.num_res_evals);
  EXPECT_EQ(a.stats.num_nonlin_iters, b.stats.num_nonlin_iters);
  EXPECT_EQ(a.stats.num_err_test_fails, b.stats.num_err_test_fails);
  const auto ka = generic.KPIs();
  const auto kb = legacy.KPIs();
  ASSERT_EQ(ka.size(), kb.size());
  for (size_t i = 0; i < ka.size(); ++i) {
    EXPECT_TRUE(SameBits(ka[i].purity, kb[i].purity));
    EXPECT_TRUE(SameBits(ka[i].recovery, kb[i].recovery));
    EXPECT_TRUE(SameBits(ka[i].productivity_mass, kb[i].productivity_mass));
    EXPECT_TRUE(SameBits(ka[i].energy_consumption, kb[i].energy_consumption));
    EXPECT_TRUE(SameBits(ka[i].mass_balance_error, kb[i].mass_balance_error));
  }
  ExpectSameLogs(generic.Log(), legacy.Log());
  EXPECT_TRUE(SameBits(generic.FinalState().y, legacy.FinalState().y));
  EXPECT_TRUE(SameBits(generic.FinalState().yp, legacy.FinalState().yp));
}

// The deck as a MixedAdsorbentLayer v6 archive: the layer class version (on
// the first layer, where cereal writes it) set to 6 and every Particle_Shape
// field removed. Returns the number of fields removed.
int ToV6(const std::string& v7, std::string& v6) {
  CEREAL_RAPIDJSON_NAMESPACE::Document doc;
  doc.Parse(v7.c_str());
  EXPECT_FALSE(doc.HasParseError());
  auto& layers = doc.MemberBegin()->value["Bed"]["Layers"]["Layers"];
  int removed = 0;
  for (auto& layer : layers.GetArray()) {
    if (layer.HasMember("cereal_class_version")) {
      EXPECT_EQ(layer["cereal_class_version"].GetInt(), 7);
      layer["cereal_class_version"].SetInt(6);
    }
    if (layer.RemoveMember("Particle_Shape")) ++removed;
  }
  CEREAL_RAPIDJSON_NAMESPACE::StringBuffer buffer;
  CEREAL_RAPIDJSON_NAMESPACE::Writer<CEREAL_RAPIDJSON_NAMESPACE::StringBuffer> writer(buffer);
  doc.Accept(writer);
  v6 = buffer.GetString();
  return removed;
}

void Shorten(Simulation& sim) {
  sim.SetParameter("css.method", 0);
  sim.SetParameter("cycle.number_of_cycles", std::min(2.0, sim.GetParameter("cycle.number_of_cycles")));
  for (auto& layer : sim.Data().bed.layers_.layers_) layer.domain_.n_cells = 10;
}

}  // namespace

TEST(ParticleShapeDecks, ExampleDecksAreSpheresAndRunBitwiseLikeV6Archives) {
  std::vector<std::string> names;
  for (const auto& e : std::filesystem::directory_iterator(kExamples))
    if (e.path().extension() == ".json") names.push_back(e.path().stem().string());
  std::sort(names.begin(), names.end());
  ASSERT_GE(names.size(), 43u);
  for (const auto& name : names) {
    SCOPED_TRACE(name);
    Simulation current = Simulation::FromFile((kExamples / (name + ".json")).string());
    for (const auto& layer : current.Data().bed.layers_.layers_)
      EXPECT_EQ(layer.particle_shape_, ParticleShape::Sphere) << layer.name_;
    std::string v6;
    const int removed = ToV6(current.ToJson(), v6);
    EXPECT_EQ(removed, static_cast<int>(current.Data().bed.layers_.layers_.size()));
    Simulation old = Simulation::FromJson(v6);
    Shorten(current);
    Shorten(old);
    ExpectSameRun(current, old);
  }
}

}  // namespace adsim
