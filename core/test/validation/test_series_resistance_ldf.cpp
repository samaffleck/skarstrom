/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Series-resistance LDF coefficient (MassTransferCoefficientType::SERIES_RESISTANCE):
//
//   1/k = f K r_p / (3 k_f) + K tau r_p^2 / (15 eps_p D_m) + r_c^2 / (15 D_c),
//   D_c = D_c0 exp(-E_a / (R T)),   K = rho_p dq*/dc (floored at 1).
//
// Checked here:
//   * the descriptor, the registry leaves and the deck round trip;
//   * each term and their sum against the closed forms at a few states, with
//     K from the isotherm slope (Henry and Langmuir), fixed, and the floor;
//   * the Arrhenius temperature dependence of the micropore term;
//   * with K = 1 and no micropore term the coefficient is haynes bit for bit,
//     per cell and over a whole run of a packed-bed example deck;
//   * a short breakthrough with all three terms (qualitative, logged).

#include "validation_helpers.h"

#include <cstdio>
#include <filesystem>

#include "adsim/math/constants.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/deck_edit.h"
#include "adsim/objects/model_kinds.h"
#include "adsim/simulation.h"

namespace adsim {

using namespace adsim::validation;

namespace {

const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";

constexpr int kCells = 5;
constexpr AdsorbentID kAds = 1;
const std::string kCO2 = "CarbonDioxide";
constexpr double kDp = 2.0e-3;   // particle diameter [m]
constexpr double kRp = kDp / 2;  // particle radius [m]
constexpr double kEp = 0.35;     // macropore void fraction [-]
constexpr double kTau = 3.0;     // tortuosity [-]
constexpr double kDm = 1.6e-5;   // state molecular diffusivity [m2/s]
constexpr double kSh = 4.0;      // state Sherwood number [-]
constexpr double kRho = 1100.0;  // particle density [kg/m3]
const double R = constants::kIdealGasConstant;

// A uniform state at (T, P, y) with the transport properties the coefficient reads.
struct State {
  Equations eq;
  MixedAdsorbentLayer layer{"bed"};
  Fluid fluid;

  State(double T, double P, double y) {
    eq.fluid_temperature.var.assign(kCells, T);
    eq.solid_temperature.var.assign(kCells, T);
    eq.pressure.var.assign(kCells, P);
    eq.sherwood_number.var.assign(kCells, kSh);
    Equation yi{"Mole Fraction", CellType::CellCentre, EquationType::Algebraic, nullptr};
    yi.var.assign(kCells, y);
    eq.mole_fraction.emplace(kCO2, std::move(yi));
    IndependentEquation dm{"Molecular Diffusion", CellType::CellCentre, nullptr};
    dm.var.assign(kCells, kDm);
    eq.molecular_diffusion.emplace(kCO2, std::move(dm));
    IndependentEquation k{"MTC", CellType::CellCentre, nullptr};
    k.var.assign(kCells, -1.0);
    eq.mass_transfer_coefficient[kAds].emplace(kCO2, std::move(k));
    IndependentEquation q{"q*", CellType::CellCentre, nullptr};
    q.var.assign(kCells, -7.0);  // sentinel: the coefficient must leave it as it was
    eq.eq_ads_concentration[kAds].emplace(kCO2, std::move(q));
    layer.domain_.start_index = 1;
    layer.domain_.end_index = kCells - 2;
    layer.particle_diameter_ = kDp;
    fluid.components_ = {kCO2};
    fluid.molar_mass_[kCO2] = 0.04401;
  }

  double Update(IMassTransferCoefficientModel& m) {
    const auto y0 = eq.mole_fraction.at(kCO2).var;
    m.UpdateMassTransferCoefficient(eq, kAds, fluid, kCO2, layer, nullptr);
    EXPECT_EQ(eq.mole_fraction.at(kCO2).var, y0);  // restored
    for (double q : eq.eq_ads_concentration.at(kAds).at(kCO2).var) EXPECT_EQ(q, -7.0);  // restored
    const auto& k = eq.mass_transfer_coefficient.at(kAds).at(kCO2).var;
    for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n)
      EXPECT_EQ(k[n], k[layer.domain_.start_index]);  // uniform state -> uniform k
    return k[layer.domain_.start_index];
  }
};

// The three resistances [s], written out independently of the model.
double Film(double K, double f = 1.0) { return f * K * kRp / (3.0 * (kSh * kDm / kDp)); }
double Macro(double K, double tau = kTau) { return K * tau * kRp * kRp / (15.0 * kEp * kDm); }
double Micro(double rc, double Dc0, double Ea, double T) { return rc * rc / (15.0 * Dc0 * std::exp(-Ea / (R * T))); }

MassTransferCoefficientSeriesResistance Series() {
  MassTransferCoefficientSeriesResistance m;
  m.tortuosity_ = kTau;
  m.macropore_void_fraction_ = kEp;
  return m;
}

}  // namespace

// ---------------------------------------------------------------------------
// Descriptor, registry and round trip
// ---------------------------------------------------------------------------

TEST(SeriesResistanceLdf, DescriptorAppendsTheKind) {
  // Appended: the stored values of the older kinds do not move.
  EXPECT_EQ(static_cast<int>(MassTransferCoefficientType::HAYNES_1973), 3);
  EXPECT_EQ(static_cast<int>(MassTransferCoefficientType::MACROPORE_CAPACITY), 5);
  EXPECT_EQ(static_cast<int>(MassTransferCoefficientType::SERIES_RESISTANCE), 6);
  const KindDescriptor* k = FindKind(ModelFamily::RateCoefficient, "series_resistance");
  ASSERT_NE(k, nullptr);
  EXPECT_EQ(k, KindAt(ModelFamily::RateCoefficient, static_cast<size_t>(MassTransferCoefficientType::SERIES_RESISTANCE)));
  EXPECT_EQ(k->label, MassTransferCoefficientTypeList[static_cast<size_t>(MassTransferCoefficientType::SERIES_RESISTANCE)]);
  for (const char* p : {"film_resistance_factor", "tortuosity", "macropore_void_fraction",
                        "molecular_diffusion_coefficient", "crystal_radius", "micropore_diffusivity",
                        "micropore_activation_energy", "capacity_factor", "concentration_floor"})
    EXPECT_NE(k->Find(p), nullptr) << p;
  EXPECT_EQ(k->Find("crystal_radius")->unit, "m");
  EXPECT_EQ(k->Find("micropore_diffusivity")->unit, "m2/s");
  EXPECT_EQ(k->Find("micropore_activation_energy")->unit, "J/mol");
  // Micropore term off and K from the isotherm by default.
  EXPECT_EQ(k->Find("crystal_radius")->default_value, 0.0);
  EXPECT_EQ(k->Find("micropore_diffusivity")->default_value, 0.0);
  EXPECT_EQ(k->Find("capacity_factor")->default_value, 0.0);
  EXPECT_EQ(k->Find("film_resistance_factor")->default_value, 1.0);
}

TEST(SeriesResistanceLdf, RegistryLeavesAndRoundTrip) {
  Simulation sim;
  UserData& d = sim.Data();
  build::AddComponent(d, kCO2);
  const AdsorbentID a = build::AddAdsorbent(d, "13X");
  build::SetIsotherm(d, a, "henry", {kCO2});
  build::SetRateModel(d, a, kCO2, "ldf", "series_resistance");
  sim.RefreshParameters();
  const std::string kp = "adsorbent[13X].kinetics[" + kCO2 + "].";
  sim.SetParameterText(kp + "crystal_radius", "1.5 [um]");
  sim.SetParameterText(kp + "micropore_diffusivity", "2e-12 [m2/s]");
  sim.SetParameterText(kp + "micropore_activation_energy", "20 [kJ/mol]");
  sim.SetParameter(kp + "capacity_factor", 250.0);
  sim.SetParameter(kp + "film_resistance_factor", 0.5);
  EXPECT_DOUBLE_EQ(sim.GetParameter(kp + "crystal_radius"), 1.5e-6);
  EXPECT_DOUBLE_EQ(sim.GetParameter(kp + "micropore_activation_energy"), 2e4);
  EXPECT_THROW(sim.SetParameterText(kp + "crystal_radius", "2 [bar]"), std::invalid_argument);

  const std::string json = sim.ToJson();
  EXPECT_NE(json.find("\"Crystal_Radius\""), std::string::npos);
  Simulation back;
  back.LoadJson(json);
  auto& m = back.Data().adsorbent_directory.Get(a).adsorption_source_.mass_transfer_models_.at(kCO2)->mtc_model_;
  ASSERT_EQ(m->GetModelType(), MassTransferCoefficientType::SERIES_RESISTANCE);
  for (const char* p : {"crystal_radius", "micropore_diffusivity", "micropore_activation_energy", "capacity_factor",
                        "film_resistance_factor", "tortuosity", "macropore_void_fraction", "concentration_floor"})
    EXPECT_EQ(back.GetParameter(kp + p), sim.GetParameter(kp + p)) << p;
}

TEST(SeriesResistanceLdf, KindSwitchFromHaynesKeepsTheSharedValues) {
  Simulation sim;
  UserData& d = sim.Data();
  build::AddComponent(d, kCO2);
  const AdsorbentID a = build::AddAdsorbent(d, "13X");
  build::SetIsotherm(d, a, "henry", {kCO2});
  build::SetRateModel(d, a, kCO2, "ldf", "haynes");
  sim.RefreshParameters();
  const std::string kp = "adsorbent[13X].kinetics[" + kCO2 + "].";
  sim.SetParameter(kp + "tortuosity", 4.5);
  sim.SetParameter(kp + "macropore_void_fraction", 0.27);
  std::string err;
  ASSERT_TRUE(SetMassTransferKind(d, a, kCO2, MassTransferModelType::LDF, MassTransferCoefficientType::SERIES_RESISTANCE, err))
      << err;
  sim.RefreshParameters();
  EXPECT_EQ(sim.GetParameter(kp + "tortuosity"), 4.5);
  EXPECT_EQ(sim.GetParameter(kp + "macropore_void_fraction"), 0.27);
}

// ---------------------------------------------------------------------------
// The series sum against the closed forms
// ---------------------------------------------------------------------------

TEST(SeriesResistanceLdf, TermsAndSumMatchTheHandCalculation) {
  // Henry isotherm q* = k e^(Q/RT) y P: K = rho_p dq*/dc = rho_p k e^(Q/RT) R T,
  // independent of P and y.
  HenryIsotherm iso(kCO2);
  iso.params_[kCO2].k = 2.0e-9;
  iso.params_[kCO2].Q = 2.5e4;
  const double rc = 1.0e-6, Dc0 = 2.0e-12, Ea = 2.0e4;
  for (double T : {298.15, 333.15}) {
    for (double P : {1.0e5, 2.0e4}) {
      for (double y : {4e-4, 0.15}) {
        SCOPED_TRACE("T " + std::to_string(T) + " P " + std::to_string(P) + " y " + std::to_string(y));
        const double K = kRho * 2.0e-9 * std::exp(2.5e4 / (R * T)) * R * T;
        ASSERT_GT(K, 1.0);
        State s(T, P, y);
        auto m = Series();
        m.SetEquilibriumContext(&iso, kRho, nullptr);
        const double tol = 1e-5;  // the forward-difference slope

        // Film only.
        m.tortuosity_ = 0.0;
        EXPECT_NEAR(s.Update(m) * Film(K), 1.0, tol);
        // Film at half its correlation value.
        m.film_resistance_factor_ = 0.5;
        EXPECT_NEAR(s.Update(m) * Film(K, 0.5), 1.0, tol);
        // Macropore only.
        m.film_resistance_factor_ = 0.0;
        m.tortuosity_ = kTau;
        EXPECT_NEAR(s.Update(m) * Macro(K), 1.0, tol);
        // Micropore only: exact (no K).
        m.tortuosity_ = 0.0;
        m.crystal_radius_ = rc;
        m.micropore_diffusivity_ = Dc0;
        m.micropore_activation_energy_ = Ea;
        EXPECT_NEAR(s.Update(m) * Micro(rc, Dc0, Ea, T), 1.0, 1e-14);
        // All three in series.
        m.film_resistance_factor_ = 1.0;
        m.tortuosity_ = kTau;
        const double sum = Film(K) + Macro(K) + Micro(rc, Dc0, Ea, T);
        EXPECT_NEAR(s.Update(m) * sum, 1.0, tol);
        // A fixed capacity factor replaces the slope: exact.
        m.capacity_factor_ = 50.0;
        EXPECT_NEAR(s.Update(m) * (Film(50.0) + Macro(50.0) + Micro(rc, Dc0, Ea, T)), 1.0, 1e-14);
      }
    }
  }
}

TEST(SeriesResistanceLdf, CapacityFactorIsTheLangmuirSlope) {
  // q* = m b' p / (1 + b' p), b' = b e^(Q/RT): dq*/dc = m b' R T / (1 + b' p)^2,
  // which is below the secant m b' R T / (1 + b' p).
  LangmuirIsotherm iso(kCO2);
  iso.params_[kCO2].m1 = 3.0;
  iso.params_[kCO2].b1 = 1e-10;
  iso.params_[kCO2].q1 = 3.0e4;
  const double T = 300.0;
  const double bp = 1e-10 * std::exp(3.0e4 / (R * T));
  for (double p : {10.0, 1e3, 1e4, 5e4}) {
    State s(T, 1e5, p / 1e5);
    auto m = Series();
    m.film_resistance_factor_ = 0.0;  // macropore only: k = 1 / Macro(K)
    m.SetEquilibriumContext(&iso, kRho, nullptr);
    const double slope = kRho * 3.0 * bp * R * T / ((1 + bp * p) * (1 + bp * p));
    ASSERT_GT(slope, 1.0) << p;
    EXPECT_NEAR(s.Update(m) * Macro(slope), 1.0, 1e-5) << p;
  }
}

TEST(SeriesResistanceLdf, CapacityFactorFloorAndConcentrationFloor) {
  const double T = 298.15;
  // A non-adsorbing component (zero slope and zero secant): K floors at 1.
  HenryIsotherm none(kCO2);
  none.params_[kCO2].k = 0.0;
  {
    State s(T, 1e5, 0.1);
    auto m = Series();
    m.SetEquilibriumContext(&none, kRho, nullptr);
    EXPECT_NEAR(s.Update(m) * (Film(1.0) + Macro(1.0)), 1.0, 1e-14);
  }
  // No isotherm context at all: K = 1 as well.
  {
    State s(T, 1e5, 0.1);
    auto m = Series();
    EXPECT_NEAR(s.Update(m) * (Film(1.0) + Macro(1.0)), 1.0, 1e-14);
  }
  // At y = 0 the slope is taken at the concentration floor: finite, and for
  // a Henry isotherm the same as anywhere else.
  HenryIsotherm iso(kCO2);
  iso.params_[kCO2].k = 2.0e-9;
  iso.params_[kCO2].Q = 2.5e4;
  const double K = kRho * 2.0e-9 * std::exp(2.5e4 / (R * T)) * R * T;
  State s(T, 1e5, 0.0);
  auto m = Series();
  m.SetEquilibriumContext(&iso, kRho, nullptr);
  EXPECT_NEAR(s.Update(m) * (Film(K) + Macro(K)), 1.0, 1e-5);
}

TEST(SeriesResistanceLdf, MicroporeTermFollowsArrhenius) {
  // Micropore only: k = 15 D_c0 e^(-E_a / R T) / r_c^2, so
  // ln(k(T2) / k(T1)) = -(E_a / R)(1/T2 - 1/T1).
  const double rc = 2.0e-6, Dc0 = 5.0e-11, Ea = 3.0e4;
  auto m = Series();
  m.film_resistance_factor_ = 0.0;
  m.tortuosity_ = 0.0;
  m.crystal_radius_ = rc;
  m.micropore_diffusivity_ = Dc0;
  m.micropore_activation_energy_ = Ea;
  std::vector<double> T{280.0, 300.0, 330.0, 373.15};
  std::vector<double> k;
  for (double t : T) {
    State s(t, 1e5, 0.1);
    k.push_back(s.Update(m));
    EXPECT_NEAR(k.back(), 15.0 * Dc0 * std::exp(-Ea / (R * t)) / (rc * rc), 1e-13 * k.back()) << t;
  }
  for (size_t i = 1; i < T.size(); ++i) {
    EXPECT_GT(k[i], k[i - 1]);
    EXPECT_NEAR(std::log(k[i] / k[0]), -(Ea / R) * (1.0 / T[i] - 1.0 / T[0]), 1e-12);
  }
  // E_a = 0: temperature independent.
  m.micropore_activation_energy_ = 0.0;
  State a(280.0, 1e5, 0.1), b(373.15, 1e5, 0.1);
  EXPECT_EQ(a.Update(m), b.Update(m));
  // r_c = 0 or D_c0 = 0 drops the term: back to film + macropore.
  m.film_resistance_factor_ = 1.0;
  m.tortuosity_ = kTau;
  m.capacity_factor_ = 1.0;
  m.crystal_radius_ = 0.0;
  EXPECT_NEAR(a.Update(m) * (Film(1.0) + Macro(1.0)), 1.0, 1e-14);
  m.crystal_radius_ = rc;
  m.micropore_diffusivity_ = 0.0;
  EXPECT_NEAR(a.Update(m) * (Film(1.0) + Macro(1.0)), 1.0, 1e-14);
}

// ---------------------------------------------------------------------------
// Reduction to haynes, bit for bit
// ---------------------------------------------------------------------------

TEST(SeriesResistanceLdf, UnitCapacityAndNoMicroporeIsHaynesPerCell) {
  HenryIsotherm iso(kCO2);
  iso.params_[kCO2].k = 2.0e-9;
  for (bool structured : {false, true}) {
    for (double T : {280.0, 298.15, 350.0}) {
      State s(T, 1e5, 0.15);
      if (structured) {  // a laminate: slab of half-thickness l, channel k_f
        s.layer.kind_ = LayerKind::Laminate;
        s.layer.structured_.sheet_thickness_ = 4e-4;
        s.layer.structured_.gap_ = 6e-4;
      }
      MassTransferCoefficient_Haynes h;
      h.tortuosity_ = kTau;
      h.macropore_void_fraction_ = 0.2;
      auto m = Series();
      m.macropore_void_fraction_ = 0.2;
      m.capacity_factor_ = 1.0;
      m.SetEquilibriumContext(&iso, kRho, nullptr);
      EXPECT_EQ(s.Update(m), s.Update(h)) << T << " structured " << structured;
    }
  }
}

namespace {

// Ward 2022 (CO2/He breakthrough on activated carbon, haynes LDF for every
// component) on 15 cells over a shortened feed step.
Simulation ShortWard() {
  Simulation sim = Simulation::FromFile((kExamples / "Ward2022_Q100ccm.json").string());
  sim.SetParameter("cycle.number_of_cycles", 1);
  sim.SetParameter("css.method", 0);
  sim.SetParameter("step[1].time", std::min(sim.GetParameter("step[1].time"), 600.0));
  for (auto& layer : sim.Data().bed.layers_.layers_) layer.domain_.n_cells = 15;
  return sim;
}

// Every haynes coefficient switched to series_resistance (tortuosity and
// macropore void fraction carry over), with the given capacity factor.
int ConvertHaynes(Simulation& sim, double capacity_factor) {
  UserData& d = sim.Data();
  int converted = 0;
  for (const auto& [id, ads] : d.adsorbent_directory.GetAll()) {
    auto& models = d.adsorbent_directory.Get(id).adsorption_source_.mass_transfer_models_;
    for (auto& [c, model] : models) {
      if (!model || !model->mtc_model_ || model->mtc_model_->GetModelType() != MassTransferCoefficientType::HAYNES_1973)
        continue;
      std::string err;
      EXPECT_TRUE(SetMassTransferKind(d, id, c, model->GetModelType(), MassTransferCoefficientType::SERIES_RESISTANCE, err))
          << err;
      auto* s = dynamic_cast<MassTransferCoefficientSeriesResistance*>(models.at(c)->mtc_model_.get());
      EXPECT_NE(s, nullptr);
      if (s) s->capacity_factor_ = capacity_factor;
      ++converted;
    }
  }
  sim.RefreshParameters();
  return converted;
}

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

TEST(SeriesResistanceLdf, UnitCapacityAndNoMicroporeIsHaynesOverARun) {
  Simulation haynes = ShortWard();
  Simulation series = ShortWard();
  ASSERT_GT(ConvertHaynes(series, 1.0), 0);
  const RunResult a = haynes.Run();
  const RunResult b = series.Run();
  ASSERT_TRUE(a.Ok()) << a.message;
  ASSERT_TRUE(b.Ok()) << b.message;
  EXPECT_EQ(a.stats.num_steps, b.stats.num_steps);
  EXPECT_EQ(a.stats.num_res_evals, b.stats.num_res_evals);
  EXPECT_EQ(a.stats.num_nonlin_iters, b.stats.num_nonlin_iters);
  const auto ka = haynes.KPIs();
  const auto kb = series.KPIs();
  ASSERT_EQ(ka.size(), kb.size());
  for (size_t i = 0; i < ka.size(); ++i) {
    EXPECT_EQ(ka[i].purity, kb[i].purity);
    EXPECT_EQ(ka[i].recovery, kb[i].recovery);
    EXPECT_EQ(ka[i].productivity_mass, kb[i].productivity_mass);
    EXPECT_EQ(ka[i].energy_consumption, kb[i].energy_consumption);
    EXPECT_EQ(ka[i].mass_balance_error, kb[i].mass_balance_error);
  }
  ExpectSameLogs(haynes.Log(), series.Log());
  EXPECT_EQ(haynes.FinalState().y, series.FinalState().y);
}

// ---------------------------------------------------------------------------
// A breakthrough with all three terms (qualitative)
// ---------------------------------------------------------------------------
//
// CO2 at 400 ppm in nitrogen fed to a clean isothermal bed (Henry isotherm,
// q* = 0.05 mol/kg at the feed, K ~ 3e3). Film and macropore resistances
// (scaled by K) give k ~ 7e-3 1/s; adding a micropore resistance of the same
// order (r_c = 1 um, D_c ~ 7e-16 m2/s at 298 K) roughly halves k, which
// shows as an earlier, more spread front. Logged, and checked only for its
// shape: monotone, the slower kinetics breaking through first, both reaching
// the feed.

namespace {

struct SeriesBreakthroughDeck : UserDataTest {
  void TestBody() override {}
  std::vector<double> Run(bool micropore, const std::vector<double>& times) {
    SetUp();
    AddComponents({"Nitrogen", "CarbonDioxide"});
    auto [layer_id, ads_id] = AddLayer("Sorbent");
    auto& layer = GetLayer(layer_id);
    layer.domain_.length = 0.1;
    layer.domain_.n_cells = 20;
    layer.eb_ = 0.4;
    layer.particle_diameter_ = 2e-3;
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
    auto ldf = std::make_unique<LinearDrivingForce>();
    auto mtc = std::make_unique<MassTransferCoefficientSeriesResistance>();
    mtc->tortuosity_ = 3.0;
    mtc->macropore_void_fraction_ = 0.35;
    mtc->molecular_diffusion_coefficient_ = 1.6e-5;
    if (micropore) {
      mtc->crystal_radius_ = 1e-6;
      mtc->micropore_diffusivity_ = 2.1e-12;
      mtc->micropore_activation_energy_ = 2.0e4;
    }
    ldf->mtc_model_ = std::move(mtc);
    GetAdsorbent(ads_id).adsorption_source_.mass_transfer_models_["CarbonDioxide"] = std::move(ldf);
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
    const auto& log = solver.GetLog();
    const auto y = log.GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION, "CarbonDioxide");
    std::vector<double> top(y.values.size());
    for (size_t n = 0; n < y.values.size(); ++n) top[n] = y.values[n][y.values[n].size() - 2] / 4e-4;
    std::vector<double> out;
    for (double t : times) out.push_back(interp_linear(y.times, top, t));
    return out;
  }
};

}  // namespace

TEST(SeriesResistanceBreakthrough, AllThreeTermsRunAndSpreadTheFront) {
  const std::vector<double> times{500, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000};
  const auto two = SeriesBreakthroughDeck().Run(false, times);
  const auto three = SeriesBreakthroughDeck().Run(true, times);
  for (size_t i = 0; i < times.size(); ++i)
    std::printf("  t = %5.0f s   c/c0 film+macropore %.4f   film+macropore+micropore %.4f\n", times[i], two[i],
                three[i]);
  for (size_t i = 1; i < times.size(); ++i) {
    EXPECT_GE(two[i], two[i - 1] - 1e-6);
    EXPECT_GE(three[i], three[i - 1] - 1e-6);
  }
  // Slower kinetics: an earlier breakthrough ...
  EXPECT_GT(three[1], two[1] + 0.02);
  // ... and both reach the feed.
  EXPECT_GT(two.back(), 0.98);
  EXPECT_GT(three.back(), 0.95);
}

}  // namespace adsim
