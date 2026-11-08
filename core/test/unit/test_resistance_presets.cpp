/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// The generic mass-transfer coefficient model (MassTransferCoefficientSeriesResistance,
// mass_transfer_model.h): every older coefficient kind is a preset of it.
//
// Checked here:
//   * each old class and its conversion give the same coefficient bit for bit
//     (memcmp, not a tolerance) on packed and structured layers, at cells
//     below the concentration floor, at c = 0, T < 0 and T = 0, with a
//     single-component and an IAST two-component isotherm;
//   * version-1 series_resistance archives (each zero-sentinel pattern) load
//     with the switches they meant and the same coefficient;
//   * saves are stable after the first load (v1 -> v2 -> v2), and every old
//     tag is still accepted as a preset;
//   * the switch cascades of deck_edit.h and the derived preset tag;
//   * the slab (structured) series value against a hand calculation;
//   * the deck checks and warnings.

#include <cmath>
#include <cstring>
#include <filesystem>
#include <functional>
#include <map>
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
using T = MassTransferCoefficientType;

const std::string kCO2 = "CarbonDioxide";
const std::string kN2 = "Nitrogen";
constexpr AdsorbentID kAds = 1;
constexpr int kCells = 9;  // ghost cells 0 and 8, the layer 1..7
const double kR = constants::kIdealGasConstant;
const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";

enum class Geometry { Packed, Monolith, MonolithWashcoat, Laminate };

const char* Name(Geometry g) {
  switch (g) {
    case Geometry::Packed: return "packed";
    case Geometry::Monolith: return "monolith";
    case Geometry::MonolithWashcoat: return "monolith with washcoat";
    case Geometry::Laminate: return "laminate";
  }
  return "";
}

// A non-uniform two-component state: ordinary cells, a cell far below the
// concentration floor, c = 0, T < 0 and (optionally) T = 0.
struct CellState {
  Equations eq;
  MixedAdsorbentLayer layer{"bed"};
  Fluid fluid;

  CellState(Geometry g, bool t_zero) {
    const std::vector<double> T{300.0, 300.0, 298.0, 320.0, -5.0, t_zero ? 0.0 : 310.0, 350.0, 280.0, 280.0};
    const std::vector<double> P{1e5, 1e5, 2e5, 1e5, 1e5, 1e5, 5e4, 1.2e5, 1.2e5};
    const std::vector<double> y{0.15, 0.15, 1e-9, 0.0, 0.1, 0.1, 0.5, 0.05, 0.05};
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
      k.var.assign(kCells, -1.0);  // a skipped cell keeps this
      eq.mass_transfer_coefficient[kAds].emplace(c, std::move(k));
      IndependentEquation q{"q*", CellType::CellCentre, nullptr};
      q.var.assign(kCells, -7.0);  // the coefficient must restore it
      eq.eq_ads_concentration[kAds].emplace(c, std::move(q));
    }
    layer.domain_.start_index = 1;
    layer.domain_.end_index = kCells - 2;
    layer.particle_diameter_ = 2e-3;
    switch (g) {
      case Geometry::Packed: break;
      case Geometry::Monolith: layer.kind_ = LayerKind::Monolith; break;
      case Geometry::MonolithWashcoat:
        layer.kind_ = LayerKind::Monolith;
        layer.structured_.washcoat_thickness_ = 3e-5;
        break;
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
  void Update(IMassTransferCoefficientModel& m) {
    m.UpdateMassTransferCoefficient(eq, kAds, fluid, kCO2, layer, nullptr);
  }
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

// A two-component Langmuir adsorbent: none (each isotherm on its own) or IAST.
struct Sorbent {
  UserData d;
  AdsorbentID id = 0;
  explicit Sorbent(bool iast) {
    build::AddComponent(d, kCO2);
    build::AddComponent(d, kN2);
    id = build::AddAdsorbent(d, "Sorbent");
    build::SetIsotherm(d, id, "single_site_langmuir_p", {kCO2, kN2});
    if (iast) build::SetMixingRule(d, id, "iast");
  }
  AdsorptionSource& Source() { return d.adsorbent_directory.Get(id).adsorption_source_; }
  IIsotherm* Isotherm() { return Source().GetIsotherm(kCO2).get(); }
};

using Factory = std::function<std::unique_ptr<IMassTransferCoefficientModel>()>;

// The old class and its conversion on the same state: equal bit for bit, and
// both leave the state's mole fractions and q* as they were.
void ExpectSameCoefficient(const Factory& make, Geometry g, bool t_zero, Sorbent* sorbent, bool compile) {
  SCOPED_TRACE(std::string(Name(g)) + (t_zero ? ", T = 0 cell" : "") + (sorbent ? ", isotherm" : "") +
               (compile ? ", compiled" : ""));
  auto legacy = make();
  auto generic = ToResistances(make());
  ASSERT_NE(dynamic_cast<Generic*>(generic.get()), nullptr);
  // The preset matches the old kind (an arrhenius coefficient with E_a = 0 is the constant preset).
  auto* arrhenius = dynamic_cast<MassTransferCoefficientArrhenius*>(legacy.get());
  EXPECT_EQ(generic->GetModelType(), MassTransferCoefficientType::SERIES_RESISTANCE);
  EXPECT_EQ(CoefficientPreset(generic.get()), arrhenius && arrhenius->activation_energy_.Get() == 0.0
                                         ? MassTransferCoefficientType::CONSTANT
                                         : legacy->GetModelType());
  CellState a(g, t_zero), b(g, t_zero);
  for (auto* m : {legacy.get(), generic.get()}) {
    if (sorbent) m->SetEquilibriumContext(sorbent->Isotherm(), 1100.0, &sorbent->Source());
    if (compile) m->CompileForRun(nullptr);
  }
  a.Update(*legacy);
  b.Update(*generic);
  EXPECT_TRUE(SameBits(a.K(), b.K())) << "old " << Dump(a.K()) << "\nnew " << Dump(b.K());
  for (const auto& c : {kCO2, kN2}) {
    EXPECT_TRUE(SameBits(a.eq.eq_ads_concentration.at(kAds).at(c).var, b.eq.eq_ads_concentration.at(kAds).at(c).var));
    EXPECT_TRUE(SameBits(a.eq.mole_fraction.at(c).var, b.eq.mole_fraction.at(c).var));
  }
}

const std::vector<Geometry> kGeometries{Geometry::Packed, Geometry::Monolith, Geometry::MonolithWashcoat,
                                        Geometry::Laminate};

// ---- archives ---------------------------------------------------------------------------

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
  return ToResistances(std::move(m));  // as IMassTransferModel's loader does
}

// The version-1 archive of a series model with these values: the v2 save with
// the fields version 2 added removed and the class version set back to 1.
std::string VersionOne(const LegacySeriesResistanceV1& v1) {
  auto g = std::make_unique<Generic>();
  g->film_resistance_factor_ = v1.film_resistance_factor_;
  g->tortuosity_ = v1.tortuosity_;
  g->macropore_void_fraction_ = v1.macropore_void_fraction_;
  g->molecular_diffusion_coefficient_ = v1.molecular_diffusion_coefficient_;
  g->crystal_radius_ = v1.crystal_radius_;
  g->micropore_diffusivity_ = v1.micropore_diffusivity_;
  g->micropore_activation_energy_ = v1.micropore_activation_energy_;
  g->capacity_factor_ = v1.capacity_factor_;
  g->concentration_floor_ = v1.concentration_floor_;
  std::unique_ptr<IMassTransferCoefficientModel> base = std::move(g);
  std::string text = Save(base);
  for (const char* field : {"Resistances", "Mass_Transfer_Coefficient", "Activation_Energy", "Pore_Radius"})
    text = std::regex_replace(text, std::regex(std::string(",\\s*\"") + field + "\": \\{[^}]*\\}"), "");
  text = std::regex_replace(text, std::regex(",\\s*\"Capacity_Factor_Mode\": -?[0-9]+"), "");
  const std::string before = text;
  text = std::regex_replace(text, std::regex("(\"data\": \\{\\s*\"cereal_class_version\": )2"), "$011");  // group 01, then "1"
  EXPECT_NE(text, before) << "class version not found";
  EXPECT_EQ(text.find("Resistances"), std::string::npos);
  return text;
}

}  // namespace

// ---------------------------------------------------------------------------
// T1: each preset equals its old kind, bit for bit
// ---------------------------------------------------------------------------

TEST(ResistancePresets, EachPresetEqualsItsOldKind) {
  Sorbent plain(false), iast(true);

  // constant: a number (the run-constant fast path when compiled) and an expression of T.
  for (const char* expr : {"0.37", "0.1 + 1e-3 * T"}) {
    Factory make = [expr] {
      auto m = std::make_unique<MassTransferCoefficientConstant>();
      m->constant_mtc_.SetExpr(expr);
      m->constant_mtc_.Evaluate(EvalContext{});
      return m;
    };
    for (Geometry g : kGeometries)
      for (bool compile : {false, true}) ExpectSameCoefficient(make, g, true, nullptr, compile);
  }
  // arrhenius: E_a != 0 every cell; E_a = 0 without the T = 0 cell (the old
  // class gives exp(-0 / 0) = NaN there, the preset k).
  for (double Ea : {2.5e4, 0.0}) {
    Factory make = [Ea] {
      auto m = std::make_unique<MassTransferCoefficientArrhenius>();
      m->pre_exp_factor_ = 75.0;
      m->activation_energy_ = Ea;
      return m;
    };
    for (Geometry g : kGeometries)
      for (bool compile : {false, true}) ExpectSameCoefficient(make, g, Ea != 0.0, nullptr, compile);
  }
  // macropore and macropore_knudsen, with the state's D_m and a constant one.
  for (double Dm : {0.0, 1.6e-5}) {
    Factory macro = [Dm] {
      auto m = std::make_unique<MassTransferCoefficientMacropore>();
      m->tortuosity_ = 2.7;
      m->macropore_void_fraction_ = 0.31;
      m->molecular_diffusion_coefficient_ = Dm;
      return m;
    };
    Factory knudsen = [Dm] {
      auto m = std::make_unique<MassTransferCoefficientMacroporeKnudsen>();
      m->pore_radius_ = 5e-9;
      m->tortuosity_ = 2.7;
      m->macropore_void_fraction_ = 0.31;
      m->molecular_diffusion_coefficient_ = Dm;
      return m;
    };
    for (Geometry g : kGeometries) {
      ExpectSameCoefficient(macro, g, true, nullptr, false);
      ExpectSameCoefficient(knudsen, g, true, nullptr, false);
    }
  }
  // haynes (film + macropore).
  Factory haynes = [] {
    auto m = std::make_unique<MassTransferCoefficient_Haynes>();
    m->tortuosity_ = 2.5;
    m->macropore_void_fraction_ = 0.3;
    return m;
  };
  for (Geometry g : kGeometries) ExpectSameCoefficient(haynes, g, true, nullptr, false);
  // macropore_capacity: no isotherm, one isotherm, the IAST mixture.
  for (double Dm : {0.0, 1.6e-5}) {
    Factory capacity = [Dm] {
      auto m = std::make_unique<MassTransferCoefficientMacroporeCapacity>();
      m->tortuosity_ = 3.2;
      m->macropore_void_fraction_ = 0.36;
      m->concentration_floor_ = 2e-3;
      m->molecular_diffusion_coefficient_ = Dm;
      return m;
    };
    for (Geometry g : kGeometries) {
      ExpectSameCoefficient(capacity, g, true, nullptr, false);
      ExpectSameCoefficient(capacity, g, true, &plain, false);
      ExpectSameCoefficient(capacity, g, true, &iast, false);
    }
  }
  // series_resistance version 1 (every term, isotherm slope; and fixed K).
  for (double fixed : {0.0, 40.0}) {
    Factory series = [fixed] {
      auto m = std::make_unique<LegacySeriesResistanceV1>();
      m->film_resistance_factor_ = 0.8;
      m->macropore_void_fraction_ = 0.33;
      m->crystal_radius_ = 1e-6;
      m->micropore_diffusivity_ = 2e-12;
      m->micropore_activation_energy_ = 2e4;
      m->capacity_factor_ = fixed;
      return m;
    };
    for (Geometry g : kGeometries)
      for (Sorbent* s : {&plain, &iast}) ExpectSameCoefficient(series, g, true, s, false);
  }
}

// ---------------------------------------------------------------------------
// T2: version-1 series_resistance archives
// ---------------------------------------------------------------------------

TEST(ResistancePresets, SeriesV1LoadsUnchanged) {
  struct Pattern {
    const char* name;
    double f, tau, rc, Dc0, capacity;
    bool film, macro, micro;
    CapacityFactorMode mode;
  };
  const std::vector<Pattern> patterns{
      {"all on, slope", 1.0, 3.0, 1e-6, 2e-12, 0.0, true, true, true, CapacityFactorMode::SLOPE},
      {"film off", 0.0, 3.0, 1e-6, 2e-12, 0.0, false, true, true, CapacityFactorMode::SLOPE},
      {"macropore off", 1.0, 0.0, 1e-6, 2e-12, 0.0, true, false, true, CapacityFactorMode::SLOPE},
      {"micropore off (r_c = 0)", 1.0, 3.0, 0.0, 2e-12, 0.0, true, true, false, CapacityFactorMode::SLOPE},
      {"micropore off (D_c0 = 0)", 1.0, 3.0, 1e-6, 0.0, 0.0, true, true, false, CapacityFactorMode::SLOPE},
      {"default (micropore off)", 1.0, 3.0, 0.0, 0.0, 0.0, true, true, false, CapacityFactorMode::SLOPE},
      {"fixed K", 1.0, 3.0, 1e-6, 2e-12, 25.0, true, true, true, CapacityFactorMode::FIXED},
      {"micropore only, fixed K", 0.0, 0.0, 1e-6, 2e-12, 25.0, false, false, true, CapacityFactorMode::FIXED},
  };
  Sorbent sorbent(false);
  for (const auto& p : patterns) {
    SCOPED_TRACE(p.name);
    LegacySeriesResistanceV1 v1;
    v1.film_resistance_factor_ = p.f;
    v1.tortuosity_ = p.tau;
    v1.macropore_void_fraction_ = 0.3;
    v1.crystal_radius_ = p.rc;
    v1.micropore_diffusivity_ = p.Dc0;
    v1.micropore_activation_energy_ = 1.5e4;
    v1.capacity_factor_ = p.capacity;
    auto loaded = Load(VersionOne(v1));
    auto* g = dynamic_cast<Generic*>(loaded.get());
    ASSERT_NE(g, nullptr);
    EXPECT_EQ(g->resistances_.film, p.film);
    EXPECT_EQ(g->resistances_.macropore, p.macro);
    EXPECT_EQ(g->resistances_.micropore, p.micro);
    EXPECT_FALSE(g->resistances_.knudsen);
    EXPECT_FALSE(g->resistances_.lumped);
    EXPECT_EQ(g->CapacityMode(), p.mode);
    // A switched-off term's zero sentinels become its defaults.
    EXPECT_GT(g->film_resistance_factor_.Get(), 0.0);
    EXPECT_GT(g->tortuosity_.Get(), 0.0);
    EXPECT_GT(g->crystal_radius_.Get(), 0.0);
    EXPECT_GT(g->micropore_diffusivity_.Get(), 0.0);
    if (!p.film) EXPECT_EQ(g->film_resistance_factor_.Get(), 1.0);
    if (!p.micro && p.rc == 0.0) EXPECT_EQ(g->crystal_radius_.Get(), 1e-6);
    // The coefficient of the version-1 class, bit for bit.
    for (Geometry geo : kGeometries) {
      SCOPED_TRACE(Name(geo));
      CellState a(geo, true), b(geo, true);
      v1.SetEquilibriumContext(sorbent.Isotherm(), 1100.0, &sorbent.Source());
      g->SetEquilibriumContext(sorbent.Isotherm(), 1100.0, &sorbent.Source());
      a.Update(v1);
      b.Update(*g);
      EXPECT_TRUE(SameBits(a.K(), b.K())) << "v1 " << Dump(a.K()) << "\nv2 " << Dump(b.K());
    }
  }
}

// ---------------------------------------------------------------------------
// T4: round trip and the old tags
// ---------------------------------------------------------------------------

TEST(ResistancePresets, RoundTrip) {
  // A model alone: v1 -> load -> save (v2) -> load -> save.
  LegacySeriesResistanceV1 v1;
  v1.film_resistance_factor_ = 0.0;
  v1.crystal_radius_ = 2e-6;
  v1.micropore_diffusivity_ = 1e-12;
  const std::string second = Save(Load(VersionOne(v1)));
  const std::string third = Save(Load(second));
  EXPECT_EQ(second, third);
  EXPECT_NE(second.find("\"Resistances\""), std::string::npos);

  // Every old class: its archive loads as the generic model and saves stably.
  for (std::size_t i = 0; i < KindTags(ModelFamily::RateCoefficient).size(); ++i) {
    const auto type = static_cast<T>(i);
    if (type == T::SERIES_RESISTANCE) continue;  // not archived any more (covered above)
    SCOPED_TRACE(KindTags(ModelFamily::RateCoefficient)[i]);
    const std::string old = Save(CreateLegacyMassTransferCoefficientModel(type));
    EXPECT_EQ(old.find("SeriesResistance"), std::string::npos);
    auto loaded = Load(old);
    ASSERT_NE(dynamic_cast<Generic*>(loaded.get()), nullptr);
    EXPECT_EQ(CoefficientPreset(loaded.get()), type);
    const std::string a = Save(loaded);
    EXPECT_EQ(a, Save(Load(a)));
  }

  // A whole deck that stores old classes (haynes, constant) and the
  // version-1 series class with a dual-rate slow site.
  for (const char* deck : {"Ward2022_Q50ccm.json", "Illustrative_DAC_Breakthrough.json", "Perez2019_A.json"}) {
    SCOPED_TRACE(deck);
    Simulation sim = Simulation::FromFile((kExamples / deck).string());
    const std::string a = sim.ToJson();
    Simulation back;
    back.LoadJson(a);
    // The coefficient models (fast and slow sites) save the same after a reload.
    auto models = [](Simulation& s) {
      std::map<std::string, std::string> out;  // sorted: the deck's maps are unordered
      for (auto& [id, ads] : s.Data().adsorbent_directory.GetAll())
        for (auto& [c, m] : ads.adsorption_source_.mass_transfer_models_) {
          std::string& text = out[std::to_string(id) + " " + c];
          text = Save(m->mtc_model_);
          if (auto* dual = dynamic_cast<DualRateLinearDrivingForce*>(m.get())) text += Save(dual->slow_mtc_model_);
        }
      return out;
    };
    EXPECT_EQ(models(sim), models(back));
    EXPECT_NE(a.find("\"Resistances\""), std::string::npos);
    EXPECT_EQ(a.find("MassTransferCoefficientConstant"), std::string::npos);
    EXPECT_EQ(a.find("MassTransferCoefficient_Haynes"), std::string::npos);
    EXPECT_EQ(a.find("MassTransferCoefficientMacroporeCapacity"), std::string::npos);
  }

  // Every old tag is accepted as a preset, by the builder and by the kind switch.
  Simulation sim;
  UserData& d = sim.Data();
  build::AddComponent(d, kCO2);
  const AdsorbentID a = build::AddAdsorbent(d, "13X");
  for (const auto& tag : KindTags(ModelFamily::RateCoefficient)) {
    SCOPED_TRACE(tag);
    const auto type = static_cast<T>(FindKind(ModelFamily::RateCoefficient, tag) -
                                     KindDescriptors(ModelFamily::RateCoefficient).data());
    for (const char* rate : {"ldf", "ldf2", "dual_rate_ldf"}) {
      build::SetRateModel(d, a, kCO2, rate, tag);
      auto& m = d.adsorbent_directory.Get(a).adsorption_source_.mass_transfer_models_.at(kCO2);
      ASSERT_NE(dynamic_cast<Generic*>(m->mtc_model_.get()), nullptr);
      EXPECT_EQ(m->mtc_model_->GetModelType(), T::SERIES_RESISTANCE);
      // The tag of the generic model gives its defaults: the lumped term only.
      EXPECT_EQ(CoefficientPreset(m->mtc_model_.get()), type == T::SERIES_RESISTANCE ? T::CONSTANT : type);
    }
    std::string err;
    ASSERT_TRUE(SetMassTransferKind(d, a, kCO2, MassTransferModelType::LDF, type, err)) << err;
    EXPECT_EQ(CoefficientPreset(d.adsorbent_directory.Get(a).adsorption_source_.mass_transfer_models_.at(kCO2)->mtc_model_.get()),
              type == T::SERIES_RESISTANCE ? T::CONSTANT : type);
  }
}

// ---------------------------------------------------------------------------
// T5: the switch cascades
// ---------------------------------------------------------------------------

TEST(DeckEdit, ResistanceSwitchCascade) {
  UserData d;
  build::AddComponent(d, kCO2);
  const AdsorbentID a = build::AddAdsorbent(d, "13X");
  build::SetRateModel(d, a, kCO2, "ldf", "constant");
  auto model = [&]() -> Generic& {
    return dynamic_cast<Generic&>(*d.adsorbent_directory.Get(a).adsorption_source_.mass_transfer_models_.at(kCO2)->mtc_model_);
  };
  auto set = [&](const char* tag, bool on, MassTransferSite site = MassTransferSite::Fast) {
    std::string err;
    const bool ok = SetMassTransferResistance(d, a, kCO2, site, tag, on, err);
    EXPECT_EQ(ok, err.empty()) << err;
    return ok;
  };
  auto mode = [&](CapacityFactorMode m) {
    std::string err;
    EXPECT_TRUE(SetCapacityFactorMode(d, a, kCO2, MassTransferSite::Fast, m, err)) << err;
  };
  EXPECT_EQ(model().PresetTag(), "constant");

  // The last term cannot be switched off; knudsen needs macropore.
  EXPECT_FALSE(set("lumped", false));
  EXPECT_FALSE(set("knudsen", true));
  EXPECT_FALSE(set("barrier", true));
  EXPECT_TRUE(model().Enabled("lumped"));

  // The preset tag follows the switches.
  model().tortuosity_ = 4.2;
  ASSERT_TRUE(set("macropore", true));
  EXPECT_EQ(model().PresetTag(), "custom");
  ASSERT_TRUE(set("lumped", false));
  EXPECT_EQ(model().PresetTag(), "macropore");
  EXPECT_EQ(model().PresetType(), T::MACROPORE);
  EXPECT_EQ(model().GetModelType(), T::SERIES_RESISTANCE);
  ASSERT_TRUE(set("knudsen", true));
  EXPECT_EQ(model().PresetTag(), "macropore_knudsen");
  ASSERT_TRUE(set("knudsen", false));
  mode(CapacityFactorMode::SECANT);
  EXPECT_EQ(model().PresetTag(), "macropore_capacity");
  mode(CapacityFactorMode::NONE);
  ASSERT_TRUE(set("film", true));
  EXPECT_EQ(model().PresetTag(), "haynes");
  model().film_resistance_factor_ = 2.0;
  EXPECT_EQ(model().PresetTag(), "custom");
  model().film_resistance_factor_ = 1.0;
  EXPECT_EQ(model().PresetTag(), "haynes");
  mode(CapacityFactorMode::SLOPE);
  EXPECT_EQ(model().PresetTag(), "custom");
  mode(CapacityFactorMode::NONE);

  // Parameters survive a switch off and on.
  ASSERT_TRUE(set("macropore", false));
  EXPECT_EQ(model().tortuosity_.Get(), 4.2);
  ASSERT_TRUE(set("macropore", true));
  EXPECT_EQ(model().tortuosity_.Get(), 4.2);
  // Switching macropore off clears knudsen.
  ASSERT_TRUE(set("knudsen", true));
  ASSERT_TRUE(set("macropore", false));
  EXPECT_FALSE(model().Enabled("knudsen"));
  EXPECT_FALSE(model().resistances_.knudsen);
  // film is now the only term: refused off.
  EXPECT_FALSE(set("film", false));
  // The fixed mode with K <= 0 sets K to 1.
  model().capacity_factor_ = 0.0;
  mode(CapacityFactorMode::FIXED);
  EXPECT_EQ(model().capacity_factor_.Get(), 1.0);
  EXPECT_EQ(model().CapacityParameters(), std::vector<std::string>{"capacity_factor"});

  // A kind switch to a preset keeps the parameters it does not show.
  std::string err;
  model().crystal_radius_ = 3e-6;
  ASSERT_TRUE(SetMassTransferKind(d, a, kCO2, MassTransferModelType::LDF, T::CONSTANT, err)) << err;
  EXPECT_EQ(model().PresetTag(), "constant");
  EXPECT_EQ(model().crystal_radius_.Get(), 3e-6);
  EXPECT_EQ(model().tortuosity_.Get(), 4.2);

  // series_resistance (custom) keeps the switches, mode and values.
  ASSERT_TRUE(SetMassTransferKind(d, a, kCO2, MassTransferModelType::LDF, T::SERIES_RESISTANCE, err)) << err;
  EXPECT_EQ(model().PresetTag(), "constant");
  EXPECT_EQ(model().crystal_radius_.Get(), 3e-6);
  // A first enable replaces an out-of-range value by the term's default; the
  // values of a term switched off are kept.
  model().tortuosity_ = 0.0;
  EXPECT_FALSE(model().Enabled("macropore"));
  ASSERT_TRUE(set("macropore", true));
  EXPECT_EQ(model().tortuosity_.Get(), 3.0) << "a first enable replaces the zero by the default";
  ASSERT_TRUE(set("macropore", false));
  EXPECT_EQ(model().tortuosity_.Get(), 3.0);

  // Dual rate: the slow site starts as a copy, then its switches are its own.
  ASSERT_TRUE(SetMassTransferKind(d, a, kCO2, MassTransferModelType::DUAL_RATE_LDF, T::MACROPORE, err)) << err;
  auto& dual = dynamic_cast<DualRateLinearDrivingForce&>(*d.adsorbent_directory.Get(a).adsorption_source_.mass_transfer_models_.at(kCO2));
  auto& slow = dynamic_cast<Generic&>(*dual.slow_mtc_model_);
  EXPECT_EQ(slow.PresetTag(), "macropore");
  ASSERT_TRUE(set("micropore", true, MassTransferSite::Slow));
  EXPECT_TRUE(slow.Enabled("micropore"));
  EXPECT_FALSE(model().Enabled("micropore"));
  EXPECT_EQ(model().PresetTag(), "macropore");
  ASSERT_TRUE(set("film", true));
  EXPECT_FALSE(slow.Enabled("film"));
  // The deck sync leaves an existing slow site alone.
  EXPECT_FALSE(dual.SyncSlowSite());
  EXPECT_TRUE(slow.Enabled("micropore"));
  // A single-site model has no slow site.
  ASSERT_TRUE(SetMassTransferKind(d, a, kCO2, MassTransferModelType::LDF, T::MACROPORE, err)) << err;
  EXPECT_FALSE(SetMassTransferResistance(d, a, kCO2, MassTransferSite::Slow, "film", true, err));
  EXPECT_NE(err.find("slow site"), std::string::npos);
}

// ---------------------------------------------------------------------------
// T8: the slab series value
// ---------------------------------------------------------------------------

TEST(ResistanceStructured, SlabSanity) {
  for (Geometry g : {Geometry::Monolith, Geometry::MonolithWashcoat, Geometry::Laminate}) {
    SCOPED_TRACE(Name(g));
    CellState s(g, false);
    const auto geo = s.layer.GetDiffusionGeometry();
    const double l = geo.length;
    const double dh = s.layer.GetHydraulicDiameter();
    if (g == Geometry::Laminate) EXPECT_DOUBLE_EQ(l, 2e-4);          // half the sheet
    if (g == Geometry::MonolithWashcoat) EXPECT_DOUBLE_EQ(l, 3e-5);  // the washcoat
    EXPECT_DOUBLE_EQ(geo.volume_to_surface, l);

    Generic m;
    std::string err;
    ASSERT_TRUE(m.SetResistance("film", true, err) && m.SetResistance("macropore", true, err) &&
                m.SetResistance("micropore", true, err) && m.SetResistance("lumped", false, err))
        << err;
    m.tortuosity_ = 2.0;
    m.macropore_void_fraction_ = 0.4;
    m.crystal_radius_ = 1e-6;
    m.micropore_diffusivity_ = 3e-12;
    m.micropore_activation_energy_ = 0.0;
    s.Update(m);
    const int n = 1;
    const double Sh = s.eq.sherwood_number.var[n];
    const double Dm = s.eq.molecular_diffusion.at(kCO2).var[n];
    const double film = l / (Sh * Dm / dh);
    const double macro = 2.0 * l * l / (3.0 * 0.4 * Dm);
    const double micro = 1e-6 * 1e-6 / (15.0 * 3e-12);
    EXPECT_NEAR(s.K()[n] * (film + macro + micro), 1.0, 1e-14);

    // Film only: k proportional to Sh.
    ASSERT_TRUE(m.SetResistance("macropore", false, err) && m.SetResistance("micropore", false, err)) << err;
    s.Update(m);
    const double k1 = s.K()[n];
    EXPECT_NEAR(k1 * film, 1.0, 1e-14);
    for (double& sh : s.eq.sherwood_number.var) sh *= 2.0;
    s.Update(m);
    EXPECT_NEAR(s.K()[n] / k1, 2.0, 1e-14);
  }
}

// ---------------------------------------------------------------------------
// T10: deck checks and warnings
// ---------------------------------------------------------------------------

TEST(ResistanceChecks, EachErrorAndWarning) {
  auto check = [](const std::function<void(Generic&)>& edit, MassTransferModelType rate = MassTransferModelType::LDF) {
    auto model = CreateMassTransferModel(rate);
    model->mtc_model_ = std::make_unique<Generic>();
    edit(dynamic_cast<Generic&>(*model->mtc_model_));
    std::vector<std::string> errors, warnings;
    CheckMassTransferModel(*model, "the CO2 kinetics of 13X", errors, warnings);
    std::string e, w;
    for (const auto& s : errors) e += s + "\n";
    for (const auto& s : warnings) w += s + "\n";
    return std::make_pair(e, w);
  };
  auto on = [](Generic& m, std::initializer_list<const char*> tags) {
    for (const char* t : tags) m.resistances_ = [&] {
      MassTransferResistances r = m.resistances_;
      if (std::string(t) == "film") r.film = true;
      if (std::string(t) == "macropore") r.macropore = true;
      if (std::string(t) == "knudsen") r.knudsen = true;
      if (std::string(t) == "micropore") r.micropore = true;
      return r;
    }();
  };

  EXPECT_EQ(check([](Generic&) {}).first, "");
  EXPECT_EQ(check([](Generic& m) { m.k_ = 0.0; }).first, "") << "k = 0 is no uptake, allowed";
  EXPECT_NE(check([](Generic& m) { m.resistances_.lumped = false; }).first.find("no mass-transfer resistance"),
            std::string::npos);
  EXPECT_NE(check([](Generic& m) { m.k_ = -1.0; }).first.find("lumped mass transfer coefficient"), std::string::npos);
  EXPECT_NE(check([&](Generic& m) { on(m, {"macropore"}); m.tortuosity_ = 0.0; }).first.find("tortuosity"),
            std::string::npos);
  EXPECT_NE(check([&](Generic& m) { on(m, {"macropore"}); m.macropore_void_fraction_ = 1.2; }).first.find("void fraction"),
            std::string::npos);
  EXPECT_NE(check([&](Generic& m) { on(m, {"macropore"}); m.macropore_void_fraction_ = 0.0; }).first.find("void fraction"),
            std::string::npos);
  EXPECT_NE(check([&](Generic& m) { on(m, {"macropore", "knudsen"}); m.pore_radius_ = 0.0; }).first.find("pore radius"),
            std::string::npos);
  EXPECT_NE(check([&](Generic& m) { on(m, {"micropore"}); m.crystal_radius_ = 0.0; }).first.find("crystal radius"),
            std::string::npos);
  EXPECT_NE(check([&](Generic& m) { on(m, {"micropore"}); m.micropore_diffusivity_ = -1.0; }).first.find("micropore diffusivity"),
            std::string::npos);
  EXPECT_NE(check([&](Generic& m) { on(m, {"film"}); m.film_resistance_factor_ = -0.5; }).first.find("film resistance factor"),
            std::string::npos);
  EXPECT_NE(check([&](Generic& m) {
              on(m, {"film"});
              m.capacity_mode_ = static_cast<int>(CapacityFactorMode::FIXED);
              m.capacity_factor_ = 0.0;
            }).first.find("fixed capacity factor"),
            std::string::npos);
  // An off term's parameters are not checked.
  EXPECT_EQ(check([](Generic& m) { m.tortuosity_ = 0.0; m.crystal_radius_ = 0.0; }).first, "");

  // Warnings.
  const auto film0 = check([&](Generic& m) { on(m, {"film"}); m.film_resistance_factor_ = 0.0; });
  EXPECT_EQ(film0.first, "");
  EXPECT_NE(film0.second.find("film resistance factor of 0"), std::string::npos);
  const auto ldf2 = check([&](Generic& m) { on(m, {"macropore"}); }, MassTransferModelType::LDF2);
  EXPECT_NE(ldf2.second.find("second-order LDF"), std::string::npos);
  EXPECT_EQ(check([](Generic&) {}, MassTransferModelType::LDF2).second, "") << "lumped only: no warning";

  // In a deck: an error stops the checks; a warning is a run-start deck warning.
  Simulation sim = Simulation::FromFile((kExamples / "Ward2022_Q50ccm.json").string());
  std::string error;
  ASSERT_TRUE(sim.Data().PerformChecks(error)) << error;
  Generic* haynes = nullptr;
  for (auto& [id, ads] : sim.Data().adsorbent_directory.GetAll())
    for (auto& [c, m] : ads.adsorption_source_.mass_transfer_models_)
      if (auto* g = dynamic_cast<Generic*>(m->mtc_model_.get()); g && g->PresetType() == T::HAYNES_1973) haynes = g;
  ASSERT_NE(haynes, nullptr);
  haynes->film_resistance_factor_ = 0.0;
  EXPECT_TRUE(sim.Data().PerformChecks(error)) << error;
  const auto warnings = MassTransferWarnings(sim.Data().adsorbent_directory, sim.Data().bed);
  ASSERT_FALSE(warnings.empty());
  EXPECT_EQ(warnings.front().rfind("Warning: ", 0), 0u);
  haynes->tortuosity_ = -1.0;
  error.clear();
  EXPECT_FALSE(sim.Data().PerformChecks(error));
  EXPECT_NE(error.find("tortuosity"), std::string::npos) << error;
}

}  // namespace adsim
