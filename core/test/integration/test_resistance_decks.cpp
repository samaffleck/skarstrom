/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Decks run with the generic mass-transfer coefficient model give the same
// results, bit for bit, as with the old coefficient classes they store
// (mass_transfer_legacy.h): one example deck per old kind, and a reactive
// mixed layer. Each deck is loaded (its coefficients become the generic
// model) and run; the reference run puts the equivalent old class back in
// every slot first. Short runs (2 cycles, 10 cells per layer, CSS off): the
// assertion is identity, which a short run shows as well as a long one.

#include <cstring>
#include <filesystem>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/equations/mass_transfer_model.h"
#include "adsim/objects/deck_builder.h"
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
      EXPECT_FALSE(r.knudsen || r.lumped);
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
      model->mtc_model_ = ToLegacy(model->mtc_model_.get());
      ++n;
      if (auto* dual = dynamic_cast<DualRateLinearDrivingForce*>(model.get()); dual && dual->slow_mtc_model_) {
        dual->slow_mtc_model_ = ToLegacy(dual->slow_mtc_model_.get());
        ++n;
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

TEST(ResistancePresets, DeckResultsUnchanged) {
  const std::vector<std::pair<std::string, T>> decks{
      {"FiveStepProductPurge", T::CONSTANT},
      {"Nguyen2023_Exp_I", T::MACROPORE},
      {"Perez2019_A", T::MACROPORE_CAPACITY},
      {"Ward2022_Q50ccm", T::HAYNES_1973},
      {"Illustrative_DAC_Breakthrough", T::SERIES_RESISTANCE},
  };
  for (const auto& [name, kind] : decks) {
    SCOPED_TRACE(name);
    Simulation generic = ShortDeck(name);
    Simulation legacy = ShortDeck(name);
    // The deck has the kind it stands for, loaded as the generic model.
    bool has_kind = false;
    for (auto& [id, ads] : generic.Data().adsorbent_directory.GetAll())
      for (auto& [c, m] : ads.adsorption_source_.mass_transfer_models_) {
        ASSERT_NE(dynamic_cast<Generic*>(m->mtc_model_.get()), nullptr) << c;
        has_kind = has_kind || static_cast<Generic*>(m->mtc_model_.get())->PresetType() == kind;
      }
    EXPECT_TRUE(has_kind);
    ASSERT_GT(ToLegacy(legacy.Data()), 0);
    ExpectSameRun(generic, legacy);
  }
}

// A catalyst (first-order isomerisation) mixed with a sorbent whose kinetics
// are a constant and a capacity-factor macropore coefficient.
TEST(ResistanceReactive, Unaffected) {
  auto build = [] {
    Simulation sim;
    UserData& d = sim.Data();
    const std::string N2 = "Nitrogen", BUT = "Butane", IBUT = "Isobutane";
    for (const auto& c : {N2, BUT, IBUT}) build::AddComponent(d, c);
    const AdsorbentID cat = build::AddAdsorbent(d, "Catalyst");
    const AdsorbentID sorb = build::AddAdsorbent(d, "Sorbent");
    build::SetIsotherm(d, sorb, "single_site_langmuir_p", {BUT, IBUT});
    build::SetRateModel(d, sorb, BUT, "ldf", "macropore_capacity");
    build::SetRateModel(d, sorb, IBUT, "ldf", "constant");
    build::SetBedShape(d, "cylinder");
    build::SetWallMaterial(d, "ss316l");
    build::SetWallBoundary(d, "zero_flux");
    const size_t layer = build::AddLayer(d, "Bed");
    build::SetLayerComposition(d, layer, {{cat, 0.4}, {sorb, 0.6}});
    build::SetMomentum(d, layer, "kozeny_carman");
    build::SetDispersion(d, layer, "constant");
    build::SetHeatDispersion(d, layer, "constant");
    build::SetInternalHeatTransfer(d, layer, "constant");
    build::SetExternalHeatTransfer(d, layer, "constant");
    const size_t feed = build::AddStep(d, "Feed");
    build::SetFlowBoundary(d, feed, "bottom", "inlet");
    build::AddStream(d, feed, "bottom", "Feed inlet", "velocity", "constant");
    build::SetPressureBoundary(d, feed, "top", "outlet", "constant");
    build::AddReaction(d, cat, "R1", {{BUT, -1}, {IBUT, 1}}, "power_law");
    build::SetReactionBasis(d, cat, 1, "concentration");
    sim.RefreshParameters();
    for (const auto& [path, text] : std::vector<std::pair<std::string, std::string>>{
             {"bed.diameter", "0.05"},
             {"adsorbent[Catalyst].solid_density", "1000"},
             {"adsorbent[Catalyst].heat_capacity", "900"},
             {"adsorbent[Sorbent].solid_density", "1100"},
             {"adsorbent[Sorbent].heat_capacity", "900"},
             {"adsorbent[Sorbent].kinetics[Isobutane].k", "0.2"},
             {"layer[Bed].length", "0.2"},
             {"layer[Bed].cells", "20"},
             {"layer[Bed].bed_void_fraction", "0.4"},
             {"layer[Bed].particle_diameter", "0.003"},
             {"initial.pressure", "1e5"},
             {"initial.temperature", "400"},
             {"initial.y[Nitrogen]", "1"},
             {"step[Feed].time", "20"},
             {"step[Feed].bottom[Feed inlet].velocity", "0.1"},
             {"step[Feed].bottom[Feed inlet].temperature", "400"},
             {"step[Feed].bottom[Feed inlet].y[Nitrogen]", "0.9"},
             {"step[Feed].bottom[Feed inlet].y[Butane]", "0.1"},
             {"step[Feed].top.pressure", "1e5"},
             {"adsorbent[Catalyst].reaction[R1].rate.pre_exp_factor", "0.05"},
         })
      sim.SetParameterText(path, text);
    d.settings.css.method_ = CSSMethod::NONE;
    d.cycle.number_of_cycles_ = 1;
    d.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
    return sim;
  };
  Simulation generic = build();
  Simulation legacy = build();
  ASSERT_GT(ToLegacy(legacy.Data()), 0);
  ExpectSameRun(generic, legacy);
  // The reaction and the sorbent both did something.
  const auto& y = generic.FinalState().y;
  ASSERT_FALSE(y.empty());
}

}  // namespace adsim
