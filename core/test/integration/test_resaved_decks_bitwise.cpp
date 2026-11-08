/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// The example decks are stored on the generic mass-transfer coefficient model
// (series_resistance with explicit switches), not on the old coefficient
// classes (mass_transfer_legacy.h). For every example deck: the file names no
// old class, every coefficient loads as the generic model, and a run is
// bitwise identical to the same deck with the equivalent old class put back in
// every slot (so the stored form changes no result). Short runs (2 cycles, 10
// cells per layer, CSS off): the assertion is identity, which a short run shows
// as well as a long one.

#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/equations/mass_transfer_model.h"
#include "adsim/simulation.h"

namespace adsim {
namespace {

using Generic = MassTransferCoefficientSeriesResistance;
using T = MassTransferCoefficientType;

const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";

// The old class equivalent to a generic model's preset (the reverse of
// ToResistances). A model that matches no old class fails the test.
std::unique_ptr<IMassTransferCoefficientModel> ToLegacy(IMassTransferCoefficientModel* m) {
  auto* g = dynamic_cast<Generic*>(m);
  EXPECT_NE(g, nullptr);
  if (!g) return nullptr;
  switch (g->PresetType()) {
    case T::CONSTANT: {
      auto c = std::make_unique<MassTransferCoefficientConstant>();
      c->constant_mtc_ = g->k_;
      return c;
    }
    case T::ARRHENIUS: {
      auto c = std::make_unique<MassTransferCoefficientArrhenius>();
      c->pre_exp_factor_ = g->k_;
      c->activation_energy_ = g->activation_energy_;
      return c;
    }
    case T::MACROPORE: {
      auto c = std::make_unique<MassTransferCoefficientMacropore>();
      c->tortuosity_ = g->tortuosity_;
      c->macropore_void_fraction_ = g->macropore_void_fraction_;
      c->molecular_diffusion_coefficient_ = g->molecular_diffusion_coefficient_;
      return c;
    }
    case T::MACROPORE_KNUDSEN: {
      auto c = std::make_unique<MassTransferCoefficientMacroporeKnudsen>();
      c->pore_radius_ = g->pore_radius_;
      c->tortuosity_ = g->tortuosity_;
      c->macropore_void_fraction_ = g->macropore_void_fraction_;
      c->molecular_diffusion_coefficient_ = g->molecular_diffusion_coefficient_;
      return c;
    }
    case T::HAYNES_1973: {
      auto c = std::make_unique<MassTransferCoefficient_Haynes>();
      c->tortuosity_ = g->tortuosity_;
      c->macropore_void_fraction_ = g->macropore_void_fraction_;
      return c;
    }
    case T::MACROPORE_CAPACITY: {
      auto c = std::make_unique<MassTransferCoefficientMacroporeCapacity>();
      c->tortuosity_ = g->tortuosity_;
      c->macropore_void_fraction_ = g->macropore_void_fraction_;
      c->concentration_floor_ = g->concentration_floor_;
      c->molecular_diffusion_coefficient_ = g->molecular_diffusion_coefficient_;
      return c;
    }
    case T::SERIES_RESISTANCE: {
      // Only a switch set the version-1 zeros can express.
      const auto r = g->Resistances();
      if (r.knudsen || r.lumped) return nullptr;  // no old class expresses it: keep the generic model
      EXPECT_EQ(g->CapacityMode(),
                g->capacity_factor_.Get() > 0.0 ? CapacityFactorMode::FIXED : CapacityFactorMode::SLOPE);
      auto c = std::make_unique<LegacySeriesResistanceV1>();
      c->film_resistance_factor_ = r.film ? g->film_resistance_factor_ : ExpressionParam(0.0);
      c->tortuosity_ = r.macropore ? g->tortuosity_ : ExpressionParam(0.0);
      c->macropore_void_fraction_ = g->macropore_void_fraction_;
      c->molecular_diffusion_coefficient_ = g->molecular_diffusion_coefficient_;
      c->crystal_radius_ = r.micropore ? g->crystal_radius_ : ExpressionParam(0.0);
      c->micropore_diffusivity_ = g->micropore_diffusivity_;
      c->micropore_activation_energy_ = g->micropore_activation_energy_;
      c->capacity_factor_ = g->capacity_factor_;
      c->concentration_floor_ = g->concentration_floor_;
      return c;
    }
  }
  return nullptr;
}

// Puts the old classes back in every coefficient slot; returns how many.
int ToLegacy(UserData& d) {
  int n = 0;
  for (auto& [id, ads] : d.adsorbent_directory.GetAll()) {
    for (auto& [c, model] : ads.adsorption_source_.mass_transfer_models_) {
      if (!model) continue;
      if (auto legacy = ToLegacy(model->mtc_model_.get())) {
        model->mtc_model_ = std::move(legacy);
        ++n;
      }
      if (auto* dual = dynamic_cast<DualRateLinearDrivingForce*>(model.get()); dual && dual->slow_mtc_model_) {
        if (auto legacy = ToLegacy(dual->slow_mtc_model_.get())) {
          dual->slow_mtc_model_ = std::move(legacy);
          ++n;
        }
      }
    }
  }
  return n;
}

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

Simulation ShortDeck(const std::string& name) {
  Simulation sim = Simulation::FromFile((kExamples / (name + ".json")).string());
  sim.SetParameter("css.method", 0);
  sim.SetParameter("cycle.number_of_cycles", std::min(2.0, sim.GetParameter("cycle.number_of_cycles")));
  for (auto& layer : sim.Data().bed.layers_.layers_) layer.domain_.n_cells = 10;
  return sim;
}

}  // namespace

TEST(ResavedDecks, StoredOnGenericModelAndRunsBitwiseLikeLegacy) {
  std::vector<std::string> names;
  for (const auto& e : std::filesystem::directory_iterator(kExamples))
    if (e.path().extension() == ".json") names.push_back(e.path().stem().string());
  std::sort(names.begin(), names.end());
  ASSERT_GE(names.size(), 43u);
  for (const auto& name : names) {
    SCOPED_TRACE(name);
    std::ifstream in(kExamples / (name + ".json"));
    std::stringstream text;
    text << in.rdbuf();
    for (const char* old : {"MassTransferCoefficientConstant\"", "MassTransferCoefficientArrhenius\"",
                            "MassTransferCoefficientMacropore\"", "MassTransferCoefficientMacroporeKnudsen\"",
                            "MassTransferCoefficient_Haynes\"", "MassTransferCoefficientMacroporeCapacity\""})
      EXPECT_EQ(text.str().find(old), std::string::npos) << "stored on old class " << old;
    Simulation generic = ShortDeck(name);
    Simulation legacy = ShortDeck(name);
    for (auto& [id, ads] : generic.Data().adsorbent_directory.GetAll())
      for (auto& [c, m] : ads.adsorption_source_.mass_transfer_models_) {
        if (!m) continue;
        EXPECT_NE(dynamic_cast<Generic*>(m->mtc_model_.get()), nullptr) << c;
      }
    ToLegacy(legacy.Data());
    ExpectSameRun(generic, legacy);
  }
}

}  // namespace adsim
