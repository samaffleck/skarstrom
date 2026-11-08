/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Dual-rate linear driving force (MassTransferModelType::DUAL_RATE_LDF): two
// parallel LDF sites sharing the component's equilibrium loading q*,
//
//   dq1/dt = k1 (phi q* - q1),   dq2/dt = k2 ((1 - phi) q* - q2),   q = q1 + q2,
//
// with the total q and the slow-site q2 in the state. Checked here:
//   * the descriptor, the registry paths (split, slow_<name>) and the kind
//     switch that keeps the two sites' coefficient kinds equal;
//   * split = 1 is the first-order LDF bit for bit (no extra unknown);
//   * the two-rate relaxation q / q* = phi (1 - e^(-k1 t)) + (1 - phi)(1 - e^(-k2 t))
//     for a step in q*;
//   * the cycle mole balance closes with the slow-site state;
//   * the deck round trip, and the warm-start remap with the extra state;
//   * a 400 ppm breakthrough shows the slow tail.

#include "validation_helpers.h"

#include <filesystem>
#include <set>

#include "adsim/objects/deck_edit.h"
#include "adsim/objects/model_kinds.h"
#include "adsim/simulation.h"
#include "adsim/solver/cycle_state.h"

namespace adsim {

using namespace adsim::validation;

namespace {

const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";

// Perez 2019 A (four-step CO2 VSA, CO2 tracked, macropore-capacity LDF for
// both components) on 20 cells for a few cycles: under a second per run.
Simulation ShortPerez(int cycles = 2, int cells = 20) {
  Simulation sim = Simulation::FromFile((kExamples / "Perez2019_A.json").string());
  sim.SetParameter("cycle.number_of_cycles", cycles);
  sim.SetParameter("css.method", 0);
  for (auto& layer : sim.Data().bed.layers_.layers_) layer.domain_.n_cells = cells;
  return sim;
}

// Every rate model of the deck switched to dual_rate_ldf, keeping its
// coefficient kind (and so its values: same names, same 1/s unit).
void ConvertToDualRate(Simulation& sim) {
  UserData& d = sim.Data();
  for (const auto& [id, ads] : d.adsorbent_directory.GetAll()) {
    for (const auto& c : d.fluid.components_) {
      auto& models = d.adsorbent_directory.Get(id).adsorption_source_.mass_transfer_models_;
      auto it = models.find(c);
      ASSERT_TRUE(it != models.end() && it->second && it->second->mtc_model_);
      const auto coeff = it->second->mtc_model_->GetModelType();
      std::string err;
      ASSERT_TRUE(SetMassTransferKind(d, id, c, MassTransferModelType::DUAL_RATE_LDF, coeff, err)) << err;
    }
  }
  sim.RefreshParameters();
}

DualRateLinearDrivingForce* Dual(Simulation& sim, AdsorbentID id, const std::string& c) {
  auto& m = sim.Data().adsorbent_directory.Get(id).adsorption_source_.mass_transfer_models_.at(c);
  return dynamic_cast<DualRateLinearDrivingForce*>(m.get());
}

AdsorbentID FirstAdsorbent(Simulation& sim) {
  return sim.Data().bed.layers_.layers_.front().composition_.begin()->first;
}

// Every logged series of two runs, compared exactly.
void ExpectSameLogs(const EquationLog& a, const EquationLog& b) {
  const auto& la = a.GetSimulationLog();
  const auto& lb = b.GetSimulationLog();
  ASSERT_EQ(la.size(), lb.size());
  for (const auto& [cycle, steps] : la) {
    ASSERT_TRUE(lb.count(cycle));
    for (const auto& [step, log] : steps) {
      SCOPED_TRACE("cycle " + std::to_string(cycle) + " step " + std::to_string(step));
      const Log& other = lb.at(cycle).at(step);
      ASSERT_EQ(log.step_data.times, other.step_data.times);
      ASSERT_EQ(log.variables.size(), other.variables.size());
      for (const auto& [type, v] : log.variables) EXPECT_EQ(v.values, other.variables.at(type).values) << v.variable_name;
      ASSERT_EQ(log.component_variables.size(), other.component_variables.size());
      for (const auto& [comp, vars] : log.component_variables) {
        ASSERT_EQ(vars.size(), other.component_variables.at(comp).size()) << comp;
        for (const auto& [type, v] : vars)
          EXPECT_EQ(v.values, other.component_variables.at(comp).at(type).values) << v.variable_name;
      }
    }
  }
}

}  // namespace

// ---------------------------------------------------------------------------
// The kind, its descriptor and the registry
// ---------------------------------------------------------------------------

TEST(DualRateLdf, DescriptorAppendsTheKind) {
  // Appended to the enum: ldf and ldf2 keep their stored values.
  EXPECT_EQ(static_cast<int>(MassTransferModelType::LDF), 0);
  EXPECT_EQ(static_cast<int>(MassTransferModelType::LDF2), 1);
  const KindDescriptor* k = FindKind(ModelFamily::RateModel, "dual_rate_ldf");
  ASSERT_NE(k, nullptr);
  EXPECT_EQ(k, KindAt(ModelFamily::RateModel, static_cast<size_t>(MassTransferModelType::DUAL_RATE_LDF)));
  EXPECT_EQ(k->label, RateModelTypeList[static_cast<size_t>(MassTransferModelType::DUAL_RATE_LDF)]);
  const KindParameter* split = k->Find("split");
  ASSERT_NE(split, nullptr);
  EXPECT_EQ(split->unit, "-");
  EXPECT_EQ(split->default_value, 1.0);
  EXPECT_GT(split->lower, 0.0);
  EXPECT_EQ(split->upper, 1.0);
  EXPECT_TRUE(split->tunable);

  // The model: split 1 (the LDF), a constant coefficient for both sites.
  DualRateLinearDrivingForce m;
  EXPECT_EQ(m.split_.Get(), 1.0);
  EXPECT_FALSE(m.HasSlowSite());
  ASSERT_NE(m.SlowCoefficientModel(), nullptr);
  EXPECT_DOUBLE_EQ(m.Rate(1.0, 0.25, 2.0), LinearDrivingForce().Rate(1.0, 0.25, 2.0));
  m.split_ = 0.6;
  EXPECT_TRUE(m.HasSlowSite());
  // Rates of the two sites at q = 0.5, q2 = 0.1 (q1 = 0.4), q* = 1.
  EXPECT_DOUBLE_EQ(DualRateLinearDrivingForce::FastRate(0.6, 1.0, 0.5, 0.1, 2.0), 2.0 * (0.6 - 0.4));
  EXPECT_DOUBLE_EQ(DualRateLinearDrivingForce::SlowRate(0.6, 1.0, 0.1, 0.5), 0.5 * (0.4 - 0.1));
}

TEST(DualRateLdf, RegistryHasSplitAndTheSlowSiteCoefficient) {
  Simulation sim = ShortPerez();
  const AdsorbentID id = FirstAdsorbent(sim);
  const std::string kp = "adsorbent[" + std::to_string(id) + "].kinetics[CarbonDioxide].";
  EXPECT_THROW(sim.GetParameter(kp + "split"), std::exception);  // LDF: no split

  ConvertToDualRate(sim);
  DualRateLinearDrivingForce* m = Dual(sim, id, "CarbonDioxide");
  ASSERT_NE(m, nullptr);
  // The slow site has the fast site's kind and, after the switch, its values.
  ASSERT_EQ(m->slow_mtc_model_->GetModelType(), m->mtc_model_->GetModelType());
  const auto fast = ModelFields(m->mtc_model_.get());
  ASSERT_FALSE(fast.empty());
  for (const auto& [name, p] : fast) {
    SCOPED_TRACE(name);
    EXPECT_EQ(sim.GetParameter(kp + "slow_" + name), p->Get());
    const ParameterInfo* info = nullptr;
    for (const auto& q : sim.Parameters().All())
      if (q.Info().path == kp + "slow_" + name) info = &q.Info();
    ASSERT_NE(info, nullptr);
    EXPECT_NE(info->name.find("slow-site"), std::string::npos);
  }
  EXPECT_EQ(sim.GetParameter(kp + "split"), 1.0);

  // Setting the slow site leaves the fast site alone.
  const std::string first = fast.front().first;
  const double before = fast.front().second->Get();
  sim.SetParameter(kp + "slow_" + first, 2.0 * before + 1.0);
  EXPECT_EQ(fast.front().second->Get(), before);
  EXPECT_EQ(ModelFields(m->slow_mtc_model_.get()).front().second->Get(), 2.0 * before + 1.0);
  sim.SetParameter(kp + "split", 0.7);
  EXPECT_EQ(m->split_.Get(), 0.7);

  // Switching the coefficient kind moves both sites; the slow site keeps
  // what the two kinds share (tortuosity, void fraction, D_m for the
  // macropore kinds).
  std::string err;
  ASSERT_TRUE(SetMassTransferKind(sim.Data(), id, "CarbonDioxide", MassTransferModelType::DUAL_RATE_LDF,
                                  MassTransferCoefficientType::MACROPORE, err)) << err;
  m = Dual(sim, id, "CarbonDioxide");
  ASSERT_NE(m, nullptr);
  EXPECT_EQ(m->mtc_model_->GetModelType(), MassTransferCoefficientType::MACROPORE);
  EXPECT_EQ(m->slow_mtc_model_->GetModelType(), MassTransferCoefficientType::MACROPORE);
  EXPECT_EQ(m->split_.Get(), 0.7);
  ASSERT_TRUE(SetMassTransferKind(sim.Data(), id, "CarbonDioxide", MassTransferModelType::DUAL_RATE_LDF,
                                  MassTransferCoefficientType::CONSTANT, err)) << err;
  m = Dual(sim, id, "CarbonDioxide");
  EXPECT_EQ(m->slow_mtc_model_->GetModelType(), MassTransferCoefficientType::CONSTANT);
}

// ---------------------------------------------------------------------------
// split = 1 is the LDF, bit for bit
// ---------------------------------------------------------------------------

TEST(DualRateLdf, SplitOneReproducesTheLdfBitForBit) {
  Simulation ldf = ShortPerez();
  Simulation dual = ShortPerez();
  ConvertToDualRate(dual);
  const RunResult a = ldf.Run();
  const RunResult b = dual.Run();
  ASSERT_TRUE(a.Ok()) << a.message;
  ASSERT_TRUE(b.Ok()) << b.message;
  EXPECT_EQ(a.stats.num_steps, b.stats.num_steps);
  EXPECT_EQ(a.stats.num_res_evals, b.stats.num_res_evals);
  EXPECT_EQ(a.stats.num_nonlin_iters, b.stats.num_nonlin_iters);
  EXPECT_FALSE(dual.Data().state.HasSlowSites());

  const auto ka = ldf.KPIs();
  const auto kb = dual.KPIs();
  ASSERT_EQ(ka.size(), kb.size());
  ASSERT_FALSE(ka.empty());
  for (size_t i = 0; i < ka.size(); ++i) {
    EXPECT_EQ(ka[i].purity, kb[i].purity);
    EXPECT_EQ(ka[i].recovery, kb[i].recovery);
    EXPECT_EQ(ka[i].productivity_mass, kb[i].productivity_mass);
    EXPECT_EQ(ka[i].energy_consumption, kb[i].energy_consumption);
    EXPECT_EQ(ka[i].mass_balance_error, kb[i].mass_balance_error);
  }
  ExpectSameLogs(ldf.Log(), dual.Log());
  EXPECT_EQ(ldf.FinalState().y, dual.FinalState().y);
}

// ---------------------------------------------------------------------------
// A slow site: an extra unknown per cell, the balance closes
// ---------------------------------------------------------------------------

namespace {

// Perez 2019 A with dual-rate CO2: 70 % of q* on the fast site (the deck's
// coefficient), 30 % on a slow site ten times slower (macropore capacity
// kind: the tortuosity scales k as 1 / tau).
Simulation SlowSitePerez(int cycles = 2, int cells = 20, double split = 0.7) {
  Simulation sim = ShortPerez(cycles, cells);
  ConvertToDualRate(sim);
  const AdsorbentID id = FirstAdsorbent(sim);
  const std::string kp = "adsorbent[" + std::to_string(id) + "].kinetics[CarbonDioxide].";
  sim.SetParameter(kp + "split", split);
  sim.SetParameter(kp + "slow_tortuosity", 10.0 * sim.GetParameter(kp + "tortuosity"));
  return sim;
}

}  // namespace

TEST(DualRateLdf, SlowSiteIsAStateAndTheMoleBalanceCloses) {
  Simulation ldf = ShortPerez();
  Simulation sim = SlowSitePerez();
  const RunResult a = ldf.Run();
  const RunResult r = sim.Run();
  ASSERT_TRUE(a.Ok()) << a.message;
  ASSERT_TRUE(r.Ok()) << r.message;

  // One more unknown in every bed cell (and the two ghost cells), for CO2 only.
  const Equations& s = sim.Data().state;
  ASSERT_TRUE(s.HasSlowSites());
  ASSERT_EQ(s.ads_concentration_slow.size(), 1u);
  EXPECT_EQ(s.ads_concentration_slow.begin()->second.count("CarbonDioxide"), 1u);
  EXPECT_EQ(s.ads_concentration_slow.begin()->second.count("Nitrogen"), 0u);
  const size_t extra = sim.FinalState().y.size() - ldf.FinalState().y.size();
  EXPECT_EQ(extra, 20u + 2u);
  const CycleStateField* f = sim.FinalState().FindField("CarbonDioxide Slow Adsorbed Concentration");
  ASSERT_NE(f, nullptr);
  EXPECT_TRUE(f->differential);

  // The slow-site loading is logged ("q slow") and lies between 0 and the total.
  const TimeSeries q2 = sim.Series("slow_adsorbed_concentration", "CarbonDioxide", SeriesLocation::Bottom);
  const TimeSeries q = sim.Series("adsorbed_concentration", "CarbonDioxide", SeriesLocation::Bottom);
  ASSERT_FALSE(q2.values.empty());
  ASSERT_EQ(q2.values.size(), q.values.size());
  for (size_t i = 0; i < q.values.size(); ++i) {
    EXPECT_GE(q2.values[i], 0.0);
    EXPECT_LT(q2.values[i], q.values[i]);
  }
  // Nitrogen has no slow site: 0.
  const TimeSeries n2 = sim.Series("slow_adsorbed_concentration", "Nitrogen", SeriesLocation::Bottom);
  for (double v : n2.values) EXPECT_EQ(v, 0.0);

  // The cycle mole balance closes (the inventory holds the total loading).
  const auto kpis = sim.KPIs();
  ASSERT_FALSE(kpis.empty());
  for (const auto& k : kpis) {
    SCOPED_TRACE("cycle " + std::to_string(k.cycle));
    ASSERT_TRUE(std::isfinite(k.quadrature_mass_balance_error));
    EXPECT_LT(std::abs(k.quadrature_mass_balance_error), 1e-4);
    EXPECT_LT(std::abs(k.mass_balance_error), 1e-3);
  }
  const BoundaryLog& b_log = sim.Boundaries();
  for (int st = 0; st < static_cast<int>(sim.Data().cycle.step_list_.size()); ++st) {
    const StepBalance* b = b_log.GetStepBalance(sim.LastCycle(), st);
    ASSERT_NE(b, nullptr);
    ASSERT_TRUE(b->valid);
    EXPECT_LT(std::abs(b->RelativeEndpointDiscrepancy()), 1e-5) << "step " << st;
  }
  // And the slower uptake shows in the KPIs.
  EXPECT_NE(kpis.back().recovery, ldf.KPIs().back().recovery);
}

TEST(DualRateLdf, SplitOutsideTheUnitIntervalIsAnError) {
  Simulation sim = ShortPerez(1);
  ConvertToDualRate(sim);
  const std::string kp = "adsorbent[" + std::to_string(FirstAdsorbent(sim)) + "].kinetics[CarbonDioxide].";
  for (double bad : {0.0, 1.5}) {
    sim.SetParameter(kp + "split", bad);
    std::string err;
    EXPECT_FALSE(sim.Check(err)) << bad;
    EXPECT_NE(err.find("split"), std::string::npos) << err;
  }
}

// ---------------------------------------------------------------------------
// Deck round trip
// ---------------------------------------------------------------------------

TEST(DualRateLdf, DeckRoundTrip) {
  Simulation sim = SlowSitePerez();
  const AdsorbentID id = FirstAdsorbent(sim);
  const std::string kp = "adsorbent[" + std::to_string(id) + "].kinetics[CarbonDioxide].";
  sim.SetParameterText(kp + "slow_concentration_floor", "2e-3");
  const std::string json = sim.ToJson();
  Simulation back = Simulation::FromJson(json);
  DualRateLinearDrivingForce* m = Dual(back, id, "CarbonDioxide");
  ASSERT_NE(m, nullptr);
  EXPECT_EQ(m->split_.Get(), 0.7);
  ASSERT_NE(m->slow_mtc_model_, nullptr);
  EXPECT_EQ(m->slow_mtc_model_->GetModelType(), Dual(sim, id, "CarbonDioxide")->slow_mtc_model_->GetModelType());
  for (const auto& p : sim.Parameters().All()) {
    if (p.Info().path.rfind(kp, 0) != 0) continue;
    SCOPED_TRACE(p.Info().path);
    EXPECT_EQ(back.GetParameter(p.Info().path), sim.GetParameter(p.Info().path));
    EXPECT_EQ(back.GetParameterText(p.Info().path), sim.GetParameterText(p.Info().path));
  }
  // The slow site survives a second round trip.
  Simulation again = Simulation::FromJson(back.ToJson());
  ASSERT_NE(Dual(again, id, "CarbonDioxide"), nullptr);
  EXPECT_EQ(again.GetParameterText(kp + "slow_concentration_floor"), "2e-3");
}

// ---------------------------------------------------------------------------
// Warm starts with the extra state
// ---------------------------------------------------------------------------

TEST(DualRateLdf, WarmStartRemapsTheSlowSite) {
  Simulation coarse = SlowSitePerez(1, 20);
  ASSERT_TRUE(coarse.Run().Ok());
  const CycleState x = coarse.FinalState();
  const AdsorbentID id = FirstAdsorbent(coarse);

  // Onto a finer grid: the slow-site field is remapped like the loading and
  // keeps its layer inventory (sum of q2 x cell width, uniform cells).
  Simulation fine = SlowSitePerez(1, 40);
  const CycleState y = fine.PrepareState(x);
  EXPECT_TRUE(fine.CheckState(y).empty());
  auto inventory = [](const CycleState& s, size_t cells, const std::string& name) {
    const CycleStateField* f = s.FindField(name);
    EXPECT_NE(f, nullptr);
    if (!f) return 0.0;
    double sum = 0.0;
    for (size_t k = 0; k < f->index.size(); ++k)
      if (f->cell[k] >= 1 && f->cell[k] <= cells) sum += s.y[f->index[k]];
    return sum / static_cast<double>(cells);
  };
  const std::string slow = "CarbonDioxide Slow Adsorbed Concentration";
  EXPECT_NEAR(inventory(y, 40, slow), inventory(x, 20, slow), 1e-12 * std::abs(inventory(x, 20, slow)));
  ASSERT_TRUE(fine.Run().Ok());

  // From an LDF state (no slow site): the slow site takes (1 - phi) of the
  // total loading, so the run starts with the same inventory.
  Simulation ldf = ShortPerez(1);
  ASSERT_TRUE(ldf.Run().Ok());
  const CycleState z = ldf.FinalState();
  Simulation dual = SlowSitePerez(1);
  std::string note;
  const CycleState w = dual.PrepareState(z, &note);
  EXPECT_NE(note.find("slow-site"), std::string::npos) << note;
  EXPECT_TRUE(dual.CheckState(w).empty());
  const auto total = w.FieldValues("Adsorbed Concentration", "CarbonDioxide", id);
  const auto q2 = w.FieldValues(slow);
  ASSERT_EQ(total.size(), q2.size());
  EXPECT_EQ(total, z.FieldValues("Adsorbed Concentration", "CarbonDioxide", id));
  for (size_t i = 0; i < q2.size(); ++i) EXPECT_DOUBLE_EQ(q2[i], (1.0 - 0.7) * total[i]);
  dual.SetInitialState(w);
  const RunResult warm = dual.Run();
  ASSERT_TRUE(warm.Ok()) << warm.message;

  // Back onto the LDF deck (or split 1): the slow site is dropped.
  Simulation again = ShortPerez(1);
  const CycleState v = again.PrepareState(x);
  EXPECT_TRUE(again.CheckState(v).empty());
  EXPECT_EQ(v.FindField(slow), nullptr);
  again.SetInitialState(v);
  ASSERT_TRUE(again.Run().Ok());

  // The Run path adapts the state itself (an optimiser's warm start).
  Simulation direct = SlowSitePerez(1);
  direct.SetInitialState(z);
  const RunResult d = direct.Run();
  EXPECT_TRUE(d.Ok()) << d.message;
}

// ---------------------------------------------------------------------------
// Two-rate relaxation after a step in q*
// ---------------------------------------------------------------------------
//
// The bed of test_ldf_kinetics.cpp: a clean bed, isothermal, a Henry isotherm
// with a tiny capacity and a fast, strongly dispersed feed, so every cell sees
// a constant q* within a few seconds while the loading relaxes over minutes.
// Each site then relaxes on its own exponential from q = 0:
//
//   q(t) / q* = phi (1 - e^(-k1 t)) + (1 - phi) (1 - e^(-k2 t)),
//   q2(t) / q* = (1 - phi) (1 - e^(-k2 t)).

class DualRateRelaxationTest : public UserDataTest {};

TEST_F(DualRateRelaxationTest, LoadingFollowsTheTwoExponentials) {
  AddComponents({"Nitrogen", "Oxygen"});
  auto [layer_id, ads_id] = AddLayer("Packing");
  auto& layer = GetLayer(layer_id);
  layer.domain_.length = 1.0;
  layer.domain_.n_cells = 20;
  layer.eb_ = 0.4;
  layer.particle_diameter_ = 0.002;
  std::get<CylinderBed>(data.bed.shape_).diameter = 0.1;
  auto& adsorbent = GetAdsorbent(ads_id);
  adsorbent.ep_ = 0.0;
  adsorbent.solid_density_ = 1000.0;
  adsorbent.heat_capacity_ = 900.0;
  adsorbent.thermal_conductivity_ = 0.0;
  AddConstantInternalHeatTransferCoefficientModel(0.0, layer_id);
  AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);

  const double henry_k = 1.0e-11, Qiso = 15000.0, T = 298.0, P = 1.0e5, y_ads = 0.5;
  const double k1 = 0.02, k2 = 0.004, phi = 0.6;
  auto isotherm = std::make_unique<HenryIsotherm>();
  isotherm->params_["Oxygen"].k = henry_k;
  isotherm->params_["Oxygen"].Q = Qiso;
  AddIsotherm({"Oxygen"}, ads_id, std::move(isotherm));
  {
    auto m = std::make_unique<DualRateLinearDrivingForce>();
    static_cast<MassTransferCoefficientConstant&>(*m->mtc_model_).constant_mtc_ = k1;
    static_cast<MassTransferCoefficientConstant&>(*m->slow_mtc_model_).constant_mtc_ = k2;
    m->split_ = phi;
    GetAdsorbent(ads_id).adsorption_source_.mass_transfer_models_["Oxygen"] = std::move(m);
  }
  SetHOA(ads_id, "Oxygen", 0.0);
  auto disp = std::make_unique<DispersionModelConstant>();
  disp->dispersion_coefficient_ = 0.05;
  GetLayer(layer_id).equations_.dispersion_model_ = std::move(disp);

  AddInitialConditions(P, T, T, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}}, 0.0);
  const double run_time = 600.0;
  AddStep("Fill", run_time);
  Step& step = data.cycle.GetStep("Fill");
  auto inlet = CreateNoValveInletStream(CreateMolarFlowBC(0.2), {{"Nitrogen", 1.0 - y_ads}, {"Oxygen", y_ads}}, T);
  auto inlet_bc = std::make_unique<FlowBoundaryCondition>();
  inlet_bc->direction_ = BoundaryDirection::INLET;
  inlet_bc->location_ = BoundaryLocation::BOTTOM;
  inlet_bc->streams_.push_back(std::move(inlet));
  step.bottom_boundary_condition_ = std::move(inlet_bc);
  auto outlet_bc = std::make_unique<PressureBoundaryCondition>();
  outlet_bc->direction_ = BoundaryDirection::OUTLET;
  outlet_bc->location_ = BoundaryLocation::TOP;
  outlet_bc->pressure_ = CreateConstantBoundaryValue(P);
  outlet_bc->mole_fraction_["Nitrogen"] = 1.0 - y_ads;
  outlet_bc->mole_fraction_["Oxygen"] = y_ads;
  step.top_boundary_condition_ = std::move(outlet_bc);
  ConfigureCycle(1);
  data.settings.reltol = 1e-9;
  data.settings.max_step = 0.5;

  CoupledSolver solver{data};
  std::string error;
  ASSERT_TRUE(data.PerformChecks(error)) << error;
  ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
  const auto& log = solver.GetLog();
  const auto& dom = log.GetDomain();

  const double q_star = henry_k * std::exp(Qiso / (kR * T)) * P * y_ads;
  const auto q = log.GetComponentVariable(0, 0, ComponentVariableType::ADSORBED_CONCENTRATION, "Oxygen");
  const auto q2 = log.GetComponentVariable(0, 0, ComponentVariableType::SLOW_ADSORBED_CONCENTRATION, "Oxygen");
  ASSERT_FALSE(q.values.empty());
  ASSERT_EQ(q2.values.size(), q.values.size());
  std::vector<double> q_mean(q.values.size()), q2_mean(q.values.size());
  for (size_t n = 0; n < q.values.size(); ++n) {
    q_mean[n] = InteriorMean(q.values[n], dom.volume) / q_star;
    q2_mean[n] = InteriorMean(q2.values[n], dom.volume) / q_star;
  }
  ASSERT_LT(q_mean.front(), 0.01);  // starts clean

  // The tolerance is the gas-fill transient (a few seconds against 1/k1 = 50
  // s), as for the single-rate test; the two rates are a factor 5 apart, so a
  // single exponential misses these points by 0.05 to 0.2.
  for (double t : {25.0, 50.0, 100.0, 200.0, 400.0, 600.0}) {
    SCOPED_TRACE("t = " + std::to_string(t));
    const double fast = phi * (1.0 - std::exp(-k1 * t));
    const double slow = (1.0 - phi) * (1.0 - std::exp(-k2 * t));
    EXPECT_NEAR(interp_linear(q.times, q_mean, t), fast + slow, 0.02);
    EXPECT_NEAR(interp_linear(q2.times, q2_mean, t), slow, 0.01);
  }
}

// ---------------------------------------------------------------------------
// A 400 ppm breakthrough: the slow tail
// ---------------------------------------------------------------------------
//
// CO2 at 400 ppm in nitrogen fed to a clean isothermal bed (Henry isotherm,
// q* = 0.05 mol/kg at the feed). With a single LDF the outlet rises to the
// feed in one front; with the same fast coefficient on 60 % of q* and a
// twenty times slower site holding the rest, the front comes earlier and the
// outlet then creeps towards the feed while the slow site keeps loading: the
// slow tail that direct-air-capture breakthroughs show (Stampi-Bombelli et
// al. 2024, doi:10.1021/acs.iecr.4c01392). Qualitative: logged, and checked
// only for its shape.

// One breakthrough run on its own deck (a fixture object per run).
struct BreakthroughDeck : UserDataTest {
  void TestBody() override {}
  // Outlet CO2 mole fraction over the feed step, sampled at `times`.
  std::vector<double> Run(bool dual, const std::vector<double>& times) {
    SetUp();
    AddComponents({"Nitrogen", "CarbonDioxide"});
    auto [layer_id, ads_id] = AddLayer("Sorbent");
    auto& layer = GetLayer(layer_id);
    layer.domain_.length = 0.1;
    layer.domain_.n_cells = 20;
    layer.eb_ = 0.4;
    std::get<CylinderBed>(data.bed.shape_).diameter = 0.02;
    auto& adsorbent = GetAdsorbent(ads_id);
    adsorbent.solid_density_ = 1000.0;
    AddConstantInternalHeatTransferCoefficientModel(0.0, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);
    auto isotherm = std::make_unique<HenryIsotherm>();
    isotherm->params_["CarbonDioxide"].k = 0.05 / 40.0;  // [mol/(kg Pa)]: 0.05 mol/kg at 40 Pa
    isotherm->params_["CarbonDioxide"].Q = 0.0;
    AddIsotherm({"CarbonDioxide"}, ads_id, std::move(isotherm));
    SetHOA(ads_id, "CarbonDioxide", 0.0);
    const double k1 = 0.01;
    if (dual) {
      auto m = std::make_unique<DualRateLinearDrivingForce>();
      static_cast<MassTransferCoefficientConstant&>(*m->mtc_model_).constant_mtc_ = k1;
      static_cast<MassTransferCoefficientConstant&>(*m->slow_mtc_model_).constant_mtc_ = k1 / 20.0;
      m->split_ = 0.6;
      GetAdsorbent(ads_id).adsorption_source_.mass_transfer_models_["CarbonDioxide"] = std::move(m);
    } else {
      AddLinearDrivingForceModel(ads_id, "CarbonDioxide", k1);
    }
    AddInitialConditions(1e5, 298.0, 298.0, {{"Nitrogen", 1.0}, {"CarbonDioxide", 0.0}}, 0.0);
    AddStep("Feed", times.back());
    std::vector<std::unique_ptr<IBoundaryStream>> in;
    in.push_back(CreateNoValveInletStream(CreateVelocityBC(0.1), {{"Nitrogen", 1.0 - 4e-4}, {"CarbonDioxide", 4e-4}},
                                          298.0));
    AddFlowBC(std::move(in), BoundaryLocation::BOTTOM, BoundaryDirection::INLET, "Feed");
    AddPressureOutletBC(CreateConstantBoundaryValue(1e5), BoundaryLocation::TOP, "Feed");
    ConfigureCycle(1);
    CoupledSolver solver{data};
    std::string error;
    EXPECT_TRUE(data.PerformChecks(error)) << error;
    EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0);
    const auto y = solver.GetLog().GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION, "CarbonDioxide");
    std::vector<double> top(y.values.size());
    for (size_t n = 0; n < y.values.size(); ++n) top[n] = y.values[n][y.values[n].size() - 2] / 4e-4;
    std::vector<double> out;
    for (double t : times) out.push_back(interp_linear(y.times, top, t));
    return out;
  }
};

TEST(DualRateBreakthrough, SlowSiteGivesTheTail) {
  const std::vector<double> times{500, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000};
  const auto single = BreakthroughDeck().Run(false, times);
  const auto dual = BreakthroughDeck().Run(true, times);
  for (size_t i = 0; i < times.size(); ++i)
    std::printf("  t = %5.0f s   c/c0 single LDF %.4f   dual-rate %.4f\n", times[i], single[i], dual[i]);
  // Earlier front: the fast site holds only 60 % of the capacity.
  EXPECT_GT(dual[1], single[1] + 0.05);
  // The tail: the single LDF is through, the dual-rate outlet still several
  // per cent below the feed and rising.
  EXPECT_GT(single.back(), 0.99);
  EXPECT_LT(dual.back(), single.back() - 0.02);
  EXPECT_GT(dual.back(), dual[times.size() - 3] + 0.01);
  // Monotone approach to the feed.
  for (size_t i = 1; i < times.size(); ++i) EXPECT_GE(dual[i], dual[i - 1] - 1e-6);
}

}  // namespace adsim
