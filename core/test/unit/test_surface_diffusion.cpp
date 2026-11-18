/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Surface diffusion in the macropore term of the generic coefficient model
// (MassTransferCoefficientSeriesResistance, switch "surface"). With pore and
// surface fluxes in parallel at local equilibrium,
//   K R_macro = K tau l^2 / (c (eps_p D + tau K D_s)),
//   i.e. k_macro = (c / l^2) (eps_p D / (tau K) + D_s),
// with D_s = D_s0 exp(-E_s / (R T)).
//
// Checked here:
//   * D_s = 0 (switch on) gives today's coefficient bit for bit (memcmp) in
//     every evaluation form, on packed spheres, cylinders and slabs; and a
//     deck run with surface on and D_s = 0 gives the same final state;
//   * k_macro against the closed form at several states (each form);
//   * a large D_s dominates; the Arrhenius temperature dependence;
//   * the switch rules, the cereal versions, the deck check and warning, the
//     display equation and the registry leaves.

#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <limits>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <cereal/archives/json.hpp>
#include <cereal/types/memory.hpp>

#include "adsim/equations/equation.h"
#include "adsim/equations/isotherm.h"
#include "adsim/equations/mass_transfer_model.h"
#include "adsim/math/constants.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/deck_edit.h"
#include "adsim/objects/fluid.h"
#include "adsim/objects/model_kinds.h"
#include "adsim/objects/porous_media.h"
#include "adsim/simulation.h"
#include "adsim/solver/deck_warnings.h"

namespace adsim {
namespace {

using Generic = MassTransferCoefficientSeriesResistance;

const std::string kCO2 = "CarbonDioxide";
const std::string kN2 = "Nitrogen";
constexpr AdsorbentID kAds = 1;
constexpr int kCells = 9;  // ghost cells 0 and 8, the layer 1..7
const double kR = constants::kIdealGasConstant;
const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";

enum class Geometry { Sphere, Cylinder, Monolith, Laminate };

const char* Name(Geometry g) {
  switch (g) {
    case Geometry::Sphere: return "sphere";
    case Geometry::Cylinder: return "cylinder";
    case Geometry::Monolith: return "monolith";
    case Geometry::Laminate: return "laminate";
  }
  return "";
}
const std::vector<Geometry> kGeometries{Geometry::Sphere, Geometry::Cylinder, Geometry::Monolith, Geometry::Laminate};

// A non-uniform two-component state; with `odd` also a cell far below the
// concentration floor, c = 0, T < 0 and T = 0.
struct CellState {
  Equations eq;
  MixedAdsorbentLayer layer{"bed"};
  Fluid fluid;

  CellState(Geometry g, bool odd) {
    const std::vector<double> T = odd ? std::vector<double>{300, 300, 298, 320, -5, 0, 350, 280, 280}
                                      : std::vector<double>{300, 300, 298, 320, 305, 310, 350, 280, 280};
    const std::vector<double> P{1e5, 1e5, 2e5, 1e5, 1e5, 1e5, 5e4, 1.2e5, 1.2e5};
    const std::vector<double> y = odd ? std::vector<double>{0.15, 0.15, 1e-9, 0.0, 0.1, 0.1, 0.5, 0.05, 0.05}
                                      : std::vector<double>{0.15, 0.15, 0.02, 0.3, 0.1, 0.1, 0.5, 0.05, 0.05};
    eq.fluid_temperature.var = T;
    eq.solid_temperature.var = T;
    eq.pressure.var = P;
    eq.sherwood_number.var.resize(kCells);
    for (int n = 0; n < kCells; ++n) eq.sherwood_number.var[n] = 3.0 + 0.1 * n;
    for (const auto& c : {kCO2, kN2}) {
      Equation yi{"Mole Fraction", CellType::CellCentre, EquationType::Algebraic, nullptr};
      yi.var.resize(kCells);
      for (int n = 0; n < kCells; ++n) yi.var[n] = c == kCO2 ? y[n] : 1.0 - y[n];
      eq.mole_fraction.emplace(c, std::move(yi));
      IndependentEquation dm{"Molecular Diffusion", CellType::CellCentre, nullptr};
      dm.var.resize(kCells);
      for (int n = 0; n < kCells; ++n) dm.var[n] = (c == kCO2 ? 1.6e-5 : 2.0e-5) * (1.0 + 0.01 * n);
      eq.molecular_diffusion.emplace(c, std::move(dm));
      IndependentEquation k{"MTC", CellType::CellCentre, nullptr};
      k.var.assign(kCells, -1.0);
      eq.mass_transfer_coefficient[kAds].emplace(c, std::move(k));
      IndependentEquation q{"q*", CellType::CellCentre, nullptr};
      q.var.assign(kCells, -7.0);
      eq.eq_ads_concentration[kAds].emplace(c, std::move(q));
    }
    layer.domain_.start_index = 1;
    layer.domain_.end_index = kCells - 2;
    layer.particle_diameter_ = 2e-3;
    switch (g) {
      case Geometry::Sphere: break;
      case Geometry::Cylinder: layer.particle_shape_ = ParticleShape::Cylinder; break;
      case Geometry::Monolith: layer.kind_ = LayerKind::Monolith; break;
      case Geometry::Laminate:
        layer.kind_ = LayerKind::Laminate;
        layer.structured_.sheet_thickness_ = 4e-4;
        layer.structured_.gap_ = 6e-4;
        break;
    }
    fluid.components_ = {kCO2, kN2};
    fluid.molar_mass_[kCO2] = 0.04401;
    fluid.molar_mass_[kN2] = 0.028013;
  }

  const std::vector<double>& K() { return eq.mass_transfer_coefficient.at(kAds).at(kCO2).var; }
  void Update(IMassTransferCoefficientModel& m) { m.UpdateMassTransferCoefficient(eq, kAds, fluid, kCO2, layer, nullptr); }
};

bool SameBits(const std::vector<double>& a, const std::vector<double>& b) {
  return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

std::string Dump(const std::vector<double>& v) {
  std::ostringstream os;
  os.precision(17);
  for (double x : v) os << x << " ";
  return os.str();
}

// A two-component Langmuir adsorbent.
struct Sorbent {
  UserData d;
  AdsorbentID id = 0;
  Sorbent() {
    build::AddComponent(d, kCO2);
    build::AddComponent(d, kN2);
    id = build::AddAdsorbent(d, "Sorbent");
    build::SetIsotherm(d, id, "single_site_langmuir_p", {kCO2, kN2});
  }
  AdsorptionSource& Source() { return d.adsorbent_directory.Get(id).adsorption_source_; }
  IIsotherm* Isotherm() { return Source().GetIsotherm(kCO2).get(); }
};

// A model with exactly `tags` on (and the macropore term's pore diffusion as
// given by them) and capacity mode `mode`.
Generic Make(const std::vector<std::string>& tags, CapacityFactorMode mode) {
  Generic m;
  for (const auto& t : tags) m.Switch(t) = true;
  m.resistances_.lumped = std::find(tags.begin(), tags.end(), "lumped") != tags.end();
  m.SetCapacityMode(mode);
  if (mode == CapacityFactorMode::FIXED) m.capacity_factor_ = 40.0;
  m.tortuosity_ = 2.7;
  m.macropore_void_fraction_ = 0.31;
  m.pore_radius_ = 5e-9;
  m.crystal_radius_ = 1e-6;
  m.micropore_diffusivity_ = 2e-12;
  m.micropore_activation_energy_ = 1.5e4;
  m.k_ = 0.4;
  m.film_resistance_factor_ = 0.9;
  return m;
}

struct Case {
  const char* name;
  std::vector<std::string> tags;
  CapacityFactorMode mode;
};
// Every evaluation form: macropore alone with mode none (with and without
// Knudsen), with the secant mode, and the series sum.
const std::vector<Case> kCases{
    {"macropore, none", {"macropore", "molecular"}, CapacityFactorMode::NONE},
    {"macropore + knudsen, none", {"macropore", "molecular", "knudsen"}, CapacityFactorMode::NONE},
    {"knudsen only, none", {"macropore", "knudsen"}, CapacityFactorMode::NONE},
    {"macropore, secant", {"macropore", "molecular"}, CapacityFactorMode::SECANT},
    {"macropore, slope", {"macropore", "molecular"}, CapacityFactorMode::SLOPE},
    {"macropore, fixed", {"macropore", "molecular"}, CapacityFactorMode::FIXED},
    {"film + macropore, none", {"film", "macropore", "molecular"}, CapacityFactorMode::NONE},
    {"film + macropore + micropore, slope", {"film", "macropore", "molecular", "micropore"}, CapacityFactorMode::SLOPE},
    {"macropore + knudsen + lumped, secant", {"macropore", "molecular", "knudsen", "lumped"}, CapacityFactorMode::SECANT},
};

std::string Save(const std::unique_ptr<IMassTransferCoefficientModel>& m) {
  std::ostringstream os;
  {
    cereal::JSONOutputArchive ar(os);
    ar(cereal::make_nvp("model", m));
  }
  return os.str();
}

std::unique_ptr<IMassTransferCoefficientModel> Load(const std::string& text) {
  std::istringstream is(text);
  cereal::JSONInputArchive ar(is);
  std::unique_ptr<IMassTransferCoefficientModel> m;
  ar(cereal::make_nvp("model", m));
  return ToResistances(std::move(m));
}

}  // namespace

// ---------------------------------------------------------------------------
// D_s = 0 is today's coefficient, bit for bit
// ---------------------------------------------------------------------------

TEST(SurfaceDiffusion, ZeroDiffusivityIsTodaysCoefficientBitForBit) {
  Sorbent sorbent;
  for (const Case& c : kCases)
    for (Geometry g : kGeometries)
      for (bool with_isotherm : {false, true})
        for (double Es : {0.0, 2.0e4}) {
          // E_s != 0 skips cells with T <= 0 (the exponent divides by T): the
          // odd cells only with E_s = 0.
          const bool odd = Es == 0.0;
          SCOPED_TRACE(std::string(c.name) + ", " + Name(g) + (with_isotherm ? ", isotherm" : "") +
                       ", E_s = " + std::to_string(Es));
          Generic off = Make(c.tags, c.mode);
          Generic on = Make(c.tags, c.mode);
          std::string err;
          ASSERT_TRUE(on.SetResistance("surface", true, err)) << err;
          ASSERT_TRUE(on.SurfaceDiffusion());
          on.surface_diffusivity_ = 0.0;
          on.surface_activation_energy_ = Es;
          if (with_isotherm)
            for (Generic* m : {&off, &on}) m->SetEquilibriumContext(sorbent.Isotherm(), 1100.0, &sorbent.Source());
          CellState a(g, odd), b(g, odd);
          a.Update(off);
          b.Update(on);
          EXPECT_TRUE(SameBits(a.K(), b.K())) << "off " << Dump(a.K()) << "\non  " << Dump(b.K());
          for (const auto& comp : {kCO2, kN2}) {
            EXPECT_TRUE(SameBits(a.eq.eq_ads_concentration.at(kAds).at(comp).var,
                                 b.eq.eq_ads_concentration.at(kAds).at(comp).var));
            EXPECT_TRUE(SameBits(a.eq.mole_fraction.at(comp).var, b.eq.mole_fraction.at(comp).var));
          }
        }
}

// A deck with film + macropore (the series form) runs to the same final
// state with surface diffusion on at D_s = 0.
TEST(SurfaceDiffusion, ZeroDiffusivityLeavesADeckRunUnchanged) {
  auto run = [](bool surface) {
    Simulation sim = Simulation::FromFile((kExamples / "Illustrative_Diluted_Bed.json").string());
    for (auto& layer : sim.Data().bed.layers_.layers_) layer.domain_.n_cells = 12;
    int switched = 0;
    if (surface) {
      for (auto& [id, ads] : sim.Data().adsorbent_directory.GetAll())
        for (auto& [c, m] : ads.adsorption_source_.mass_transfer_models_) {
          auto* g = dynamic_cast<Generic*>(m->mtc_model_.get());
          if (!g || !g->Enabled("macropore")) continue;
          std::string err;
          EXPECT_TRUE(SetMassTransferResistance(sim.Data(), id, c, MassTransferSite::Fast, "surface", true, err)) << err;
          g->surface_diffusivity_ = 0.0;
          ++switched;
        }
      EXPECT_GT(switched, 0) << "the deck has a macropore term";
    }
    const RunResult r = sim.Run();
    EXPECT_TRUE(r.Ok()) << r.message;
    return sim.FinalState().y;
  };
  const std::vector<double> before = run(false);
  const std::vector<double> after = run(true);
  ASSERT_FALSE(before.empty());
  EXPECT_TRUE(SameBits(before, after));
}

// ---------------------------------------------------------------------------
// k_macro against the closed form
// ---------------------------------------------------------------------------

TEST(SurfaceDiffusion, CoefficientMatchesTheParallelFluxFormula) {
  const double Ds0 = 3e-9, Es = 1.2e4;
  for (Geometry g : kGeometries) {
    SCOPED_TRACE(Name(g));
    CellState s(g, false);
    const auto geo = s.layer.GetDiffusionGeometry();
    const double c_over_l2 = geo.c / (geo.length * geo.length);
    for (CapacityFactorMode mode : {CapacityFactorMode::NONE, CapacityFactorMode::FIXED}) {
      SCOPED_TRACE(static_cast<int>(mode));
      Generic m = Make({"macropore", "molecular"}, mode);
      std::string err;
      ASSERT_TRUE(m.SetResistance("surface", true, err)) << err;
      m.surface_diffusivity_ = Ds0;
      m.surface_activation_energy_ = Es;
      s.Update(m);
      const double K = mode == CapacityFactorMode::FIXED ? 40.0 : 1.0;
      for (int n = 1; n <= kCells - 2; ++n) {
        const double T = s.eq.fluid_temperature.var[n];
        const double D = s.eq.molecular_diffusion.at(kCO2).var[n];
        const double Ds = Ds0 * std::exp(-Es / (kR * T));
        // k = (c / l^2) (eps_p D / (tau K) + D_s), and 1/k = K tau l^2 / (c (eps_p D + tau K D_s)).
        const double expected = c_over_l2 * (0.31 * D / (2.7 * K) + Ds);
        EXPECT_NEAR(s.K()[n], expected, 1e-13 * expected) << n;
        EXPECT_NEAR(1.0 / s.K()[n], K * 2.7 / (c_over_l2 * (0.31 * D + 2.7 * K * Ds)), 1e-13 / expected) << n;
      }
    }
  }
}

// The secant form: K = 1 / min(c / (rho_p q*), 1) from a Henry isotherm.
TEST(SurfaceDiffusion, SecantFormMatchesTheFormula) {
  HenryIsotherm iso(kCO2);
  iso.params_[kCO2].k = 2.0e-9;
  iso.params_[kCO2].Q = 2.5e4;
  const double rho_p = 1100.0, Ds0 = 5e-10;
  Generic m = Make({"macropore", "molecular"}, CapacityFactorMode::SECANT);
  std::string err;
  ASSERT_TRUE(m.SetResistance("surface", true, err)) << err;
  m.surface_diffusivity_ = Ds0;
  m.SetEquilibriumContext(&iso, rho_p, nullptr);
  CellState s(Geometry::Sphere, false);
  s.Update(m);
  const double c_over_l2 = 60.0 / (2e-3 * 2e-3);
  for (int n = 1; n <= kCells - 2; ++n) {
    const double T = s.eq.fluid_temperature.var[n];
    const double K = kR * T * rho_p * 2.0e-9 * std::exp(2.5e4 / (kR * T));  // rho_p q* / c
    ASSERT_GT(K, 1.0);
    const double D = s.eq.molecular_diffusion.at(kCO2).var[n];
    const double expected = c_over_l2 * (0.31 * D / (2.7 * K) + Ds0);
    EXPECT_NEAR(s.K()[n], expected, 1e-12 * expected) << n;
  }
}

TEST(SurfaceDiffusion, LargeSurfaceDiffusivityDominates) {
  for (Geometry g : kGeometries) {
    SCOPED_TRACE(Name(g));
    CellState s(g, false);
    const auto geo = s.layer.GetDiffusionGeometry();
    Generic m = Make({"macropore", "molecular"}, CapacityFactorMode::FIXED);
    std::string err;
    ASSERT_TRUE(m.SetResistance("surface", true, err)) << err;
    double previous = std::numeric_limits<double>::infinity();
    for (double Ds : {1e-8, 1e-6, 1e-4, 1e-2}) {
      m.surface_diffusivity_ = Ds;
      s.Update(m);
      const double surface_only = geo.c * Ds / (geo.length * geo.length);
      const double ratio = s.K()[1] / surface_only;
      EXPECT_GT(ratio, 1.0);  // the pore term adds
      EXPECT_LT(ratio, previous);
      previous = ratio;
    }
    EXPECT_NEAR(previous, 1.0, 1e-5);
    // Without surface diffusion the coefficient is far lower.
    CellState t(g, false);
    ASSERT_TRUE(m.SetResistance("surface", false, err)) << err;
    t.Update(m);
    EXPECT_LT(t.K()[1], 1e-3 * s.K()[1]);
  }
}

TEST(SurfaceDiffusion, ArrheniusTemperatureDependence) {
  Generic m;
  m.surface_diffusivity_ = 2e-8;
  m.surface_activation_energy_ = 0.0;
  EXPECT_EQ(m.SurfaceDiffusivity(300.0), 2e-8);
  EXPECT_EQ(m.SurfaceDiffusivity(0.0), 2e-8) << "E_s = 0: constant, also at T = 0";
  m.surface_activation_energy_ = 1.8e4;
  for (double T : {250.0, 300.0, 400.0}) EXPECT_DOUBLE_EQ(m.SurfaceDiffusivity(T), 2e-8 * std::exp(-1.8e4 / (kR * T)));
  // ln(D_s(T1) / D_s(T2)) = (E_s / R) (1/T2 - 1/T1).
  EXPECT_NEAR(std::log(m.SurfaceDiffusivity(400.0) / m.SurfaceDiffusivity(300.0)), 1.8e4 / kR * (1.0 / 300.0 - 1.0 / 400.0),
              1e-12);
  // A surface-dominated coefficient follows it.
  Generic k = Make({"macropore", "molecular"}, CapacityFactorMode::FIXED);
  std::string err;
  ASSERT_TRUE(k.SetResistance("surface", true, err)) << err;
  k.surface_diffusivity_ = 10.0;  // the pore term is then < 1e-5 of the surface term
  k.surface_activation_energy_ = 1.8e4;
  CellState s(Geometry::Sphere, false);
  s.Update(k);
  const double T3 = s.eq.fluid_temperature.var[3], T6 = s.eq.fluid_temperature.var[6];  // 320 K, 350 K
  const double expected = std::exp(1.8e4 / kR * (1.0 / T3 - 1.0 / T6));
  EXPECT_NEAR(s.K()[6] / s.K()[3], expected, 1e-4 * expected);
}

// ---------------------------------------------------------------------------
// Switch rules, archives, checks, equation, registry
// ---------------------------------------------------------------------------

TEST(SurfaceDiffusion, SwitchRules) {
  UserData d;
  build::AddComponent(d, kCO2);
  const AdsorbentID a = build::AddAdsorbent(d, "AC");
  build::SetRateModel(d, a, kCO2, "ldf", "constant");
  auto model = [&]() -> Generic& {
    return dynamic_cast<Generic&>(*d.adsorbent_directory.Get(a).adsorption_source_.mass_transfer_models_.at(kCO2)->mtc_model_);
  };
  auto set = [&](const char* tag, bool on) {
    std::string err;
    const bool ok = SetMassTransferResistance(d, a, kCO2, MassTransferSite::Fast, tag, on, err);
    EXPECT_EQ(ok, err.empty()) << err;
    return ok;
  };
  EXPECT_FALSE(set("surface", true)) << "surface needs macropore";
  EXPECT_TRUE(set("surface", false)) << "already off with macropore";
  ASSERT_TRUE(set("macropore", true));
  ASSERT_TRUE(set("lumped", false));
  EXPECT_EQ(model().PresetTag(), "macropore");
  ASSERT_TRUE(set("surface", true));
  EXPECT_TRUE(model().SurfaceDiffusion());
  EXPECT_TRUE(model().MolecularDiffusion()) << "surface does not replace the pore diffusion";
  EXPECT_EQ(model().PresetTag(), "custom");
  EXPECT_EQ(model().PresetType(), MassTransferCoefficientType::SERIES_RESISTANCE);
  EXPECT_FALSE(set("molecular", false)) << "the pore diffusion stays on with surface";
  model().surface_diffusivity_ = 4e-9;
  ASSERT_TRUE(set("surface", false));
  EXPECT_EQ(model().PresetTag(), "macropore");
  ASSERT_TRUE(set("surface", true));
  EXPECT_EQ(model().surface_diffusivity_.Get(), 4e-9) << "values are kept";
  // Macropore off clears surface; on again does not restore it.
  ASSERT_TRUE(set("lumped", true));
  ASSERT_TRUE(set("macropore", false));
  EXPECT_FALSE(model().resistances_.surface);
  ASSERT_TRUE(set("macropore", true));
  EXPECT_FALSE(model().SurfaceDiffusion());
  // A stored surface switch without macropore (written field by field) is off.
  Generic g;
  g.resistances_.surface = true;
  EXPECT_FALSE(g.SurfaceDiffusion());
  EXPECT_FALSE(g.Enabled("surface"));
  // The parameters of the term.
  EXPECT_EQ(Generic::ResistanceParameters("surface"),
            (std::vector<std::string>{"surface_diffusivity", "surface_activation_energy"}));
  // A preset switches it off.
  ASSERT_TRUE(set("surface", true));
  std::string err;
  ASSERT_TRUE(SetMassTransferKind(d, a, kCO2, MassTransferModelType::LDF, MassTransferCoefficientType::MACROPORE, err))
      << err;
  EXPECT_FALSE(model().SurfaceDiffusion());
}

TEST(SurfaceDiffusion, OlderArchivesLoadWithSurfaceOff) {
  Generic m = Make({"macropore", "molecular"}, CapacityFactorMode::SLOPE);
  std::string err;
  ASSERT_TRUE(m.SetResistance("surface", true, err)) << err;
  m.surface_diffusivity_ = 7e-9;
  m.surface_activation_energy_ = 1.1e4;
  std::unique_ptr<IMassTransferCoefficientModel> base = std::make_unique<Generic>(m);
  const std::string v3 = Save(base);
  // Round trip.
  auto back = Load(v3);
  auto* g = dynamic_cast<Generic*>(back.get());
  ASSERT_NE(g, nullptr);
  EXPECT_TRUE(g->SurfaceDiffusion());
  EXPECT_EQ(g->surface_diffusivity_.Get(), 7e-9);
  EXPECT_EQ(g->surface_activation_energy_.Get(), 1.1e4);
  EXPECT_EQ(Save(back), v3);

  // The archive an engine before the surface term wrote: class versions 2,
  // no "Surface" switch and no surface parameters.
  std::string v2 = v3;
  v2 = std::regex_replace(v2, std::regex(",\\s*\"Surface\": (true|false)"), "");
  for (const char* field : {"Surface_Diffusivity", "Surface_Activation_Energy"})
    v2 = std::regex_replace(v2, std::regex(std::string(",\\s*\"") + field + "\": \\{[^}]*\\}"), "");
  const std::string before = v2;
  v2 = std::regex_replace(v2, std::regex("(\"data\": \\{\\s*\"cereal_class_version\": )3"), "$012");
  ASSERT_NE(v2, before) << "model class version not found";
  const std::string before_switches = v2;
  v2 = std::regex_replace(v2, std::regex("(\"Resistances\": \\{\\s*\"cereal_class_version\": )3"), "$012");
  ASSERT_NE(v2, before_switches) << "switch class version not found";
  ASSERT_EQ(v2.find("Surface"), std::string::npos);
  auto old = Load(v2);
  g = dynamic_cast<Generic*>(old.get());
  ASSERT_NE(g, nullptr);
  EXPECT_FALSE(g->resistances_.surface);
  EXPECT_FALSE(g->SurfaceDiffusion());
  EXPECT_EQ(g->surface_diffusivity_.Get(), 1e-10);  // the default
  EXPECT_EQ(g->surface_activation_energy_.Get(), 0.0);
  EXPECT_EQ(g->PresetTag(), "custom");  // macropore, slope
  EXPECT_TRUE(g->MolecularDiffusion());
}

TEST(SurfaceDiffusion, DeckCheckAndWarning) {
  auto check = [](CapacityFactorMode mode, double Ds) {
    auto model = CreateMassTransferModel(MassTransferModelType::LDF);
    auto g = std::make_unique<Generic>(Make({"macropore", "molecular"}, mode));
    std::string err;
    EXPECT_TRUE(g->SetResistance("surface", true, err)) << err;
    g->surface_diffusivity_ = Ds;
    model->mtc_model_ = std::move(g);
    std::vector<std::string> errors, warnings;
    CheckMassTransferModel(*model, "the CO2 kinetics of AC", errors, warnings);
    std::string e, w;
    for (const auto& s : errors) e += s + "\n";
    for (const auto& s : warnings) w += s + "\n";
    return std::make_pair(e, w);
  };
  const auto none = check(CapacityFactorMode::NONE, 1e-9);
  EXPECT_EQ(none.first, "");
  EXPECT_NE(none.second.find("surface diffusion is on with the capacity-factor mode none"), std::string::npos)
      << none.second;
  for (CapacityFactorMode mode : {CapacityFactorMode::SLOPE, CapacityFactorMode::SECANT, CapacityFactorMode::FIXED}) {
    const auto r = check(mode, 1e-9);
    EXPECT_EQ(r.first, "");
    EXPECT_EQ(r.second, "") << static_cast<int>(mode);
  }
  EXPECT_EQ(check(CapacityFactorMode::SLOPE, 0.0).first, "") << "D_s = 0 is allowed";
  EXPECT_NE(check(CapacityFactorMode::SLOPE, -1e-9).first.find("surface diffusivity"), std::string::npos);

  // In a deck: a run-start deck warning.
  Simulation sim = Simulation::FromFile((kExamples / "Illustrative_Diluted_Bed.json").string());
  bool switched = false;
  for (auto& [id, ads] : sim.Data().adsorbent_directory.GetAll())
    for (auto& [c, m] : ads.adsorption_source_.mass_transfer_models_) {
      auto* g = dynamic_cast<Generic*>(m->mtc_model_.get());
      if (!g || !g->Enabled("macropore") || g->CapacityMode() != CapacityFactorMode::NONE) continue;
      std::string err;
      ASSERT_TRUE(SetMassTransferResistance(sim.Data(), id, c, MassTransferSite::Fast, "surface", true, err)) << err;
      switched = true;
    }
  ASSERT_TRUE(switched);
  std::string error;
  EXPECT_TRUE(sim.Data().PerformChecks(error)) << error;
  const auto warnings = MassTransferWarnings(sim.Data().adsorbent_directory, sim.Data().bed);
  bool found = false;
  for (const auto& w : warnings) found = found || w.find("surface diffusion") != std::string::npos;
  EXPECT_TRUE(found);
}

TEST(SurfaceDiffusion, DisplayEquation) {
  Generic m = Make({"macropore", "molecular"}, CapacityFactorMode::SLOPE);
  const std::string plain = MassTransferEquationLatex(m, DiffusionShape::Sphere);
  EXPECT_EQ(plain.find("D_{s,i}"), std::string::npos);
  std::string err;
  ASSERT_TRUE(m.SetResistance("surface", true, err)) << err;
  const std::string sphere = MassTransferEquationLatex(m, DiffusionShape::Sphere);
  EXPECT_NE(sphere.find(R"(\frac{K\,\tau\,r_p^{2}}{15\left(\varepsilon_p\,D_{m,i} + \tau\,K\,D_{s,i}\right)})"),
            std::string::npos)
      << sphere;
  EXPECT_NE(sphere.find(R"(D_{s,i} &= D_{s0})"), std::string::npos) << sphere;
  EXPECT_EQ(sphere.find("E_s"), std::string::npos) << "E_s = 0: the constant form";
  m.surface_activation_energy_ = 1e4;
  const std::string arrhenius = MassTransferEquationLatex(m, DiffusionShape::Sphere);
  EXPECT_NE(arrhenius.find(R"(D_{s,i} &= D_{s0}\,e^{-E_s/RT})"), std::string::npos) << arrhenius;
  EXPECT_NE(MassTransferEquationLatex(m, DiffusionShape::Slab).find(R"({3\left(\varepsilon_p)"), std::string::npos);
  EXPECT_NE(MassTransferEquationLatex(m, DiffusionShape::Cylinder).find(R"({8\left(\varepsilon_p)"), std::string::npos);
  // Mode none: no K in the term.
  m.SetCapacityMode(CapacityFactorMode::NONE);
  EXPECT_NE(MassTransferEquationLatex(m, DiffusionShape::Sphere).find(R"(\tau\,D_{s,i})"), std::string::npos);
  // The descriptor.
  const KindDescriptor& kind = ResistanceKind();
  const ResistanceDescriptor* surface = nullptr;
  for (const auto& r : kind.resistances)
    if (r.tag == "surface") surface = &r;
  ASSERT_NE(surface, nullptr);
  EXPECT_EQ(surface->parent, "macropore");
  EXPECT_EQ(surface->parameters, Generic::ResistanceParameters("surface"));
  for (const auto& name : surface->parameters) EXPECT_NE(kind.Find(name), nullptr) << name;
  for (DiffusionShape shape : {DiffusionShape::Sphere, DiffusionShape::Slab, DiffusionShape::Cylinder})
    EXPECT_FALSE(surface->Fragments(shape).parent_term.empty());
}

TEST(SurfaceDiffusion, RegistryLeaves) {
  Simulation sim;
  UserData& d = sim.Data();
  build::AddComponent(d, kCO2);
  const AdsorbentID a = build::AddAdsorbent(d, "AC");
  build::SetIsotherm(d, a, "henry", {kCO2});
  build::SetRateModel(d, a, kCO2, "ldf", "macropore");
  sim.RefreshParameters();
  const std::string head = "adsorbent[AC].kinetics[" + kCO2 + "].";
  ASSERT_NE(sim.Parameters().Find(head + "resistances.surface"), nullptr);
  EXPECT_EQ(sim.GetParameter(head + "resistances.surface"), 0.0);
  EXPECT_EQ(sim.Parameters().Find(head + "surface_diffusivity"), nullptr) << "off: not listed";
  sim.SetParameterText(head + "resistances.surface", "true");
  sim.RefreshParameters();
  EXPECT_EQ(sim.GetParameter(head + "resistances.surface"), 1.0);
  sim.SetParameterText(head + "surface_diffusivity", "2e-9 [m2/s]");
  EXPECT_DOUBLE_EQ(sim.GetParameter(head + "surface_diffusivity"), 2e-9);
  sim.SetParameterText(head + "surface_activation_energy", "15 [kJ/mol]");
  EXPECT_DOUBLE_EQ(sim.GetParameter(head + "surface_activation_energy"), 1.5e4);
  EXPECT_THROW(sim.SetParameterText(head + "surface_diffusivity", "2 [bar]"), std::invalid_argument);
  // Macropore off: the switch and the parameters go.
  sim.SetParameterText(head + "resistances.lumped", "true");
  sim.SetParameterText(head + "resistances.macropore", "false");
  sim.RefreshParameters();
  EXPECT_EQ(sim.Parameters().Find(head + "resistances.surface"), nullptr);
  sim.SetParameterText(head + "resistances.macropore", "true");
  sim.RefreshParameters();
  sim.SetParameterText(head + "resistances.surface", "true");
  sim.RefreshParameters();

  const std::string json = sim.ToJson();
  Simulation back;
  back.LoadJson(json);
  EXPECT_EQ(back.GetParameter(head + "resistances.surface"), 1.0);
  EXPECT_DOUBLE_EQ(back.GetParameter(head + "surface_diffusivity"), 2e-9);
}

}  // namespace adsim
