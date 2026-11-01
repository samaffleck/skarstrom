/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// The macropore LDF coefficients (macropore, macropore_capacity,
// macropore_knudsen) with the optional constant molecular diffusivity:
// > 0 replaces the state's (T, P)-dependent D_m, 0 keeps it. The expected
// values are the closed forms, written out independently of the models.

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>

#include "adsim/equations/equation.h"
#include "adsim/equations/isotherm.h"
#include "adsim/equations/mass_transfer_model.h"
#include "adsim/math/constants.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/fluid.h"
#include "adsim/objects/porous_media.h"
#include "adsim/simulation.h"

namespace adsim {
namespace {

constexpr int kCells = 5;
constexpr AdsorbentID kAds = 1;
const std::string kCO2 = "CarbonDioxide";
constexpr double kT = 298.15;   // [K]
constexpr double kDp = 2.0e-3;  // particle diameter [m]
constexpr double kEp = 0.35;    // macropore void fraction [-]
constexpr double kTau = 3.0;    // tortuosity [-]
constexpr double kDm = 1.6e-5;  // constant D_m of Estupinan Perez et al. (2019) [m2/s]

// A uniform state at pressure P whose binary-diffusion D_m follows 1/P like the
// Chapman-Enskog value (D_m = 1.6e-5 m2/s at 1 bar).
struct State {
  Equations eq;
  MixedAdsorbentLayer layer{"bed"};
  Fluid fluid;

  explicit State(double P, double y = 0.15) {
    eq.fluid_temperature.var.assign(kCells, kT);
    eq.solid_temperature.var.assign(kCells, kT);
    eq.pressure.var.assign(kCells, P);
    Equation yi{"Mole Fraction", CellType::CellCentre, EquationType::Algebraic, nullptr};
    yi.var.assign(kCells, y);
    eq.mole_fraction.emplace(kCO2, std::move(yi));
    IndependentEquation dm{"Molecular Diffusion", CellType::CellCentre, nullptr};
    dm.var.assign(kCells, StateDm(P));
    eq.molecular_diffusion.emplace(kCO2, std::move(dm));
    IndependentEquation k{"MTC", CellType::CellCentre, nullptr};
    k.var.assign(kCells, -1.0);
    eq.mass_transfer_coefficient[kAds].emplace(kCO2, std::move(k));
    IndependentEquation q{"q*", CellType::CellCentre, nullptr};
    q.var.assign(kCells, 0.0);
    eq.eq_ads_concentration[kAds].emplace(kCO2, std::move(q));
    layer.domain_.start_index = 1;
    layer.domain_.end_index = kCells - 2;
    layer.particle_diameter_ = kDp;
    fluid.molar_mass_[kCO2] = 0.04401;
  }

  static double StateDm(double P) { return 1.6e-5 * 1.0e5 / P; }

  double Update(IMassTransferCoefficientModel& m) {
    m.UpdateMassTransferCoefficient(eq, kAds, fluid, kCO2, layer, nullptr);
    const auto& k = eq.mass_transfer_coefficient.at(kAds).at(kCO2).var;
    for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n)
      EXPECT_EQ(k[n], k[layer.domain_.start_index]);  // uniform state -> uniform k
    return k[layer.domain_.start_index];
  }
};

// 15 eps_p D / (tau r_p^2) = 60 eps_p D / (tau d_p^2).
double Macropore(double D) { return 60.0 * D * kEp / (kTau * kDp * kDp); }

TEST(MacroporeCoefficient, ConstantDiffusivityReplacesTheStateValue) {
  MassTransferCoefficientMacropore m;
  m.tortuosity_ = kTau;
  m.macropore_void_fraction_ = kEp;
  EXPECT_EQ(m.molecular_diffusion_coefficient_.Get(), 0.0);  // default: the state's D_m

  for (double P : {1.0e5, 1.0e4, 3.0e3}) {
    State s(P);
    m.molecular_diffusion_coefficient_ = 0.0;
    EXPECT_DOUBLE_EQ(s.Update(m), Macropore(State::StateDm(P))) << P;  // follows 1/P
    m.molecular_diffusion_coefficient_ = kDm;
    EXPECT_DOUBLE_EQ(s.Update(m), Macropore(kDm)) << P;  // independent of P
  }
}

TEST(MacroporeCoefficient, CapacityFactorWithConstantDiffusivity) {
  // Henry isotherm q* = k e^(Q/RT) P y: the capacity ratio c / (rho_p q*) =
  // 1 / (R T rho_p k e^(Q/RT)) is independent of P and y, so with a constant
  // D_m the whole coefficient is.
  HenryIsotherm iso(kCO2);
  iso.params_[kCO2].k = 2.0e-9;
  iso.params_[kCO2].Q = 2.5e4;
  const double rho_p = 1100.0;
  const double R = constants::kIdealGasConstant;
  const double ratio = 1.0 / (R * kT * rho_p * 2.0e-9 * std::exp(2.5e4 / (R * kT)));
  ASSERT_LT(ratio, 1.0);  // not capped

  MassTransferCoefficientMacroporeCapacity m;
  m.tortuosity_ = kTau;
  m.macropore_void_fraction_ = kEp;
  m.SetEquilibriumContext(&iso, rho_p);
  EXPECT_EQ(m.molecular_diffusion_coefficient_.Get(), 0.0);

  for (double P : {1.0e5, 1.0e4, 3.0e3}) {
    State s(P);
    m.molecular_diffusion_coefficient_ = 0.0;
    EXPECT_NEAR(s.Update(m), Macropore(State::StateDm(P)) * ratio, 1e-12 * Macropore(State::StateDm(P)) * ratio) << P;
    m.molecular_diffusion_coefficient_ = kDm;
    EXPECT_NEAR(s.Update(m), Macropore(kDm) * ratio, 1e-12 * Macropore(kDm) * ratio) << P;
  }
}

TEST(MacroporeCoefficient, KnudsenWithConstantDiffusivity) {
  MassTransferCoefficientMacroporeKnudsen m;
  m.tortuosity_ = kTau;
  m.macropore_void_fraction_ = kEp;
  m.pore_radius_ = 5e-8;
  const double R = constants::kIdealGasConstant;
  const double Dk = (2.0 / 3.0) * 5e-8 * std::sqrt(8.0 * R * kT / (constants::kPI * 0.04401));
  const auto expected = [&](double D) {
    const double Deff = 1.0 / (kTau / (Dk * kEp) + kTau / (D * kEp));  // Bosanquet, in series
    return 60.0 * Deff / (kDp * kDp);
  };
  for (double P : {1.0e5, 3.0e3}) {
    State s(P);
    m.molecular_diffusion_coefficient_ = 0.0;
    EXPECT_NEAR(s.Update(m), expected(State::StateDm(P)), 1e-12 * expected(State::StateDm(P))) << P;
    m.molecular_diffusion_coefficient_ = kDm;
    EXPECT_NEAR(s.Update(m), expected(kDm), 1e-12 * expected(kDm)) << P;
  }
}

// Registry leaf, builder kind and a JSON round trip.
TEST(MacroporeCoefficient, RegistryLeafAndRoundTrip) {
  for (const char* kind : {"macropore", "macropore_capacity", "macropore_knudsen"}) {
    Simulation sim;
    UserData& d = sim.Data();
    build::AddComponent(d, kCO2);
    const AdsorbentID a = build::AddAdsorbent(d, "13X");
    build::SetIsotherm(d, a, "henry", {kCO2});
    build::SetRateModel(d, a, kCO2, "ldf", kind);
    sim.RefreshParameters();
    const std::string path = "adsorbent[13X].kinetics[" + kCO2 + "].molecular_diffusion_coefficient";
    EXPECT_EQ(sim.GetParameter(path), 0.0) << kind;
    sim.SetParameterText(path, "1.6e-5 [m2/s]");
    EXPECT_DOUBLE_EQ(sim.GetParameter(path), 1.6e-5) << kind;
    EXPECT_THROW(sim.SetParameterText(path, "2 [bar]"), std::invalid_argument) << kind;

    const std::string json = sim.ToJson();
    EXPECT_NE(json.find("\"Molecular_Diffusion_Coefficient\""), std::string::npos) << kind;
    Simulation back;
    back.LoadJson(json);
    EXPECT_DOUBLE_EQ(back.GetParameter(path), 1.6e-5) << kind;
  }
}

// Decks written before the field existed (macropore v2, macropore_capacity v1)
// load with 0, i.e. the state's D_m, exactly as before.
TEST(MacroporeCoefficient, OlderDecksLoadWithTheStateDiffusivity) {
  for (const char* deck : {"Perez2019_A", "Nguyen2022_Fig4"}) {
    const std::string file = std::string("../examples/decks/") + deck + ".json";
    std::ifstream in(file);
    ASSERT_TRUE(in) << file;
    std::stringstream text;
    text << in.rdbuf();
    ASSERT_EQ(text.str().find("Molecular_Diffusion_Coefficient"), std::string::npos) << deck;
    Simulation sim;
    sim.Load(file);
    int found = 0;
    const std::string leaf = ".molecular_diffusion_coefficient";
    for (const auto& p : sim.Parameters().All()) {
      const std::string& path = p.Info().path;
      if (path.find(".kinetics[") != std::string::npos && path.size() > leaf.size() &&
          path.compare(path.size() - leaf.size(), leaf.size(), leaf) == 0) {
        EXPECT_EQ(sim.GetParameter(path), 0.0) << deck << " " << path;
        ++found;
      }
    }
    EXPECT_GT(found, 0) << deck;
  }
}

}  // namespace
}  // namespace adsim
