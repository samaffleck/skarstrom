/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// One definition per model kind (adsim/objects/model_kinds.h): the
// descriptors agree with the model classes (parameter names, order and
// defaults), the builder's kind lists are their tags, the registry lists each
// kind's parameters with the descriptor's unit, label and bounds, and the
// isotherm fitter's kinds and default parameters follow them.

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/boundary_value.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/equations/isotherm_kinds.h"
#include "adsim/equations/mass_transfer_model.h"
#include "adsim/objects/compressor.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/deck_edit.h"
#include "adsim/objects/deck_io.h"
#include "adsim/objects/model_kinds.h"
#include "adsim/parameters/parameter_registry.h"
#include "adsim/solver/isotherm_fit.h"

namespace adsim {
namespace {

std::vector<std::string> Names(const ModelFieldList& fields) {
  std::vector<std::string> out;
  for (const auto& f : fields) out.push_back(f.first);
  return out;
}

std::vector<std::string> Names(const KindDescriptor& k) {
  std::vector<std::string> out;
  for (const auto& p : k.parameters) out.push_back(p.name);
  return out;
}

// The model's fields are the descriptor's, in order, at its defaults.
void ExpectMatches(const KindDescriptor& k, const ModelFieldList& fields, const std::string& what) {
  EXPECT_EQ(Names(fields), Names(k)) << what;
  for (const auto& [name, p] : fields) {
    const KindParameter* d = k.Find(name);
    if (!d) continue;
    EXPECT_EQ(p->Get(), d->default_value) << what << "." << name << " default";
    EXPECT_TRUE(d->lower <= d->upper) << what << "." << name << " bounds";
  }
}

UserData Perez() {
  std::ifstream in(std::filesystem::current_path() / ".." / "examples" / "decks" / "Perez2019_A.json",
                   std::ios::binary);
  std::ostringstream ss;
  ss << in.rdbuf();
  UserData d;
  LoadDeckFlexible(d, ss.str());
  return d;
}

// The registry entry at `path` has the descriptor parameter's label, unit and
// bounds (widened to the current value, as the registry does).
void ExpectEntry(ParameterRegistry& reg, const std::string& path, const KindParameter& d, const std::string& label,
                 const std::string& unit, double lo, double hi) {
  const Parameter* p = reg.Find(path);
  ASSERT_NE(p, nullptr) << path;
  const auto& i = p->Info();
  EXPECT_EQ(i.name, label) << path;
  EXPECT_EQ(i.unit, unit) << path;
  const double v = p->GetValue();
  EXPECT_DOUBLE_EQ(i.lower, std::min(lo, v)) << path;
  EXPECT_DOUBLE_EQ(i.upper, std::max(hi, v)) << path;
  if (!d.tunable) EXPECT_FALSE(i.tunable) << path;
}

void ExpectEntry(ParameterRegistry& reg, const std::string& path, const KindParameter& d) {
  ExpectEntry(reg, path, d, d.Label(), d.unit, d.lower, d.upper);
}

}  // namespace

TEST(ModelKinds, BuilderListsAreTheDescriptorTags) {
  using F = ModelFamily;
  const std::vector<std::pair<F, const std::vector<std::string>& (*)()>> lists = {
      {F::RateModel, &build::RateModelKinds},
      {F::RateCoefficient, &build::RateCoefficientKinds},
      {F::Momentum, &build::MomentumKinds},
      {F::Dispersion, &build::DispersionKinds},
      {F::HeatDispersion, &build::HeatDispersionKinds},
      {F::InternalHeatTransfer, &build::InternalHeatTransferKinds},
      {F::ExternalHeatTransfer, &build::ExternalHeatTransferKinds},
      {F::ParticleHeatTransfer, &build::ParticleHeatTransferKinds},
      {F::Profile, &build::ProfileKinds},
      {F::Flow, &build::FlowKinds},
      {F::Pump, &build::PumpKinds},
      {F::PumpEnergy, &build::PumpEnergyKinds},
      {F::PumpEfficiency, &build::PumpEfficiencyKinds},
      {F::PumpPower, &build::PumpPowerKinds},
  };
  for (const auto& [family, list] : lists) {
    // Legacy kinds (the former rate coefficients) are not offered.
    std::vector<std::string> offered;
    for (const auto& k : KindDescriptors(family))
      if (!k.legacy) offered.push_back(k.tag);
    EXPECT_EQ(list(), offered) << ModelFamilyName(family);
  }
  EXPECT_EQ(build::RateCoefficientKinds(), std::vector<std::string>{"series_resistance"});
  EXPECT_EQ(build::RateCoefficientPresets(),
            (std::vector<std::string>{"constant", "arrhenius", "macropore", "macropore_knudsen", "haynes",
                                      "macropore_capacity"}));
  // The engine's spellings (front ends must follow them).
  EXPECT_EQ(KindTags(F::Momentum).front(), "kozeny_carman");
  EXPECT_EQ(KindTags(F::Flow), (std::vector<std::string>{"velocity", "volumetric_flow", "molar_flow", "mass_flow"}));

  std::vector<std::string> iso;
  for (IsothermType t : AllIsothermTypes()) iso.push_back(IsothermKind(t).tag);
  EXPECT_EQ(build::IsothermKinds(), iso);
}

TEST(ModelKinds, DescriptorsMatchTheModelClasses) {
  using F = ModelFamily;
  auto each = [](F family, auto make) {
    const auto& kinds = KindDescriptors(family);
    for (size_t i = 0; i < kinds.size(); ++i) {
      auto model = make(i);
      ASSERT_TRUE(model) << ModelFamilyName(family) << " " << kinds[i].tag;
      ExpectMatches(kinds[i], ModelFields(model.get()), std::string(ModelFamilyName(family)) + "." + kinds[i].tag);
    }
  };
  // Rate coefficients: every kind is the generic model; a legacy kind (a
  // preset) binds its parameters under its own names (PresetFields), the
  // generic kind every parameter at the class defaults.
  for (size_t i = 0; i < KindDescriptors(F::RateCoefficient).size(); ++i) {
    const auto type = static_cast<MassTransferCoefficientType>(i);
    auto model = CreateMassTransferCoefficientModel(type);
    ASSERT_TRUE(model);
    const KindDescriptor& k = KindDescriptors(F::RateCoefficient)[i];
    EXPECT_EQ(k.legacy, type != MassTransferCoefficientType::SERIES_RESISTANCE) << k.tag;
    ExpectMatches(k, k.legacy ? PresetFields(model.get(), type) : ModelFields(model.get()), "rate_coefficient." + k.tag);
  }
  each(F::Dispersion, [](size_t i) { return CreateAxialDispersionModel(static_cast<AxialDispersionType>(i)); });
  each(F::HeatDispersion,
       [](size_t i) { return CreateAxialHeatDispersionModel(static_cast<AxialHeatDispersionType>(i)); });
  each(F::InternalHeatTransfer,
       [](size_t i) { return CreateInternalHeatTransferModel(static_cast<InternalHeatTransferType>(i)); });
  each(F::ExternalHeatTransfer,
       [](size_t i) { return CreateExternalHeatTransferModel(static_cast<ExternalHeatTransferType>(i)); });
  each(F::ParticleHeatTransfer,
       [](size_t i) { return CreateParticleHeatTransferModel(static_cast<ParticleHeatTransferType>(i)); });
  each(F::Profile, [](size_t i) { return CreateBoundaryValue(static_cast<BoundaryValueType>(i)); });
  each(F::PumpEfficiency, [](size_t i) { return CreatePumpEfficiencyModel(static_cast<PumpEfficiencyType>(i)); });
  each(F::PumpPower, [](size_t i) { return CreatePumpPowerModel(static_cast<PumpPowerType>(i)); });
  for (size_t k = 1; k < KindDescriptors(F::Pump).size(); ++k) {
    auto pump = CreatePump(static_cast<PumpType>(k - 1));
    ExpectMatches(*KindAt(F::Pump, k), ModelFields(pump.get()), "pump." + KindAt(F::Pump, k)->tag);
  }
  ValveBoundaryStream valve;
  ExpectMatches(*FindKind(F::Stream, "valve"), {{"cv", &valve.valve_.cv_}, {"xt", &valve.valve_.xt_}}, "valve");

  // Every family has one descriptor per enum value.
  EXPECT_EQ(KindDescriptors(F::RateModel).size(), static_cast<size_t>(MassTransferModelType::DUAL_RATE_LDF) + 1);
  // The rate models' own parameters bind to the model classes (dual_rate_ldf: split).
  for (size_t k = 0; k < KindDescriptors(F::RateModel).size(); ++k) {
    auto rate = CreateMassTransferModel(static_cast<MassTransferModelType>(k));
    ASSERT_NE(rate, nullptr);
    ExpectMatches(*KindAt(F::RateModel, k), ModelFields(rate.get()), "rate." + KindAt(F::RateModel, k)->tag);
  }
  EXPECT_EQ(KindDescriptors(F::RateCoefficient).size(), static_cast<size_t>(MassTransferCoefficientType::SERIES_RESISTANCE) + 1);
  EXPECT_EQ(KindDescriptors(F::Boundary).size(), static_cast<size_t>(BoundaryType::CLOSED) + 1);
  EXPECT_EQ(KindDescriptors(F::Stream).size(), static_cast<size_t>(BoundaryStreamType::PUMP) + 1);
  EXPECT_EQ(KindDescriptors(F::Profile).size(), static_cast<size_t>(BoundaryValueType::DOUBLE_EXPONENTIAL) + 1);
  EXPECT_EQ(KindDescriptors(F::Flow).size(), static_cast<size_t>(FlowType::MASS) + 1);
  EXPECT_EQ(KindDescriptors(F::Pump).size(), static_cast<size_t>(PumpType::VACUUM_PUMP) + 2);
  EXPECT_EQ(KindDescriptors(F::PumpEnergy).size(), static_cast<size_t>(CompressorEnergyType::ADIABATIC) + 1);
}

TEST(ModelKinds, IsothermDescriptorsMatchTheIsothermClasses) {
  const auto& all = IsothermKindDescriptors();
  ASSERT_EQ(all.size(), AllIsothermTypes().size());
  for (size_t i = 0; i < all.size(); ++i) {
    const auto& d = all[i];
    const auto t = static_cast<IsothermType>(i);
    EXPECT_EQ(d.type, t);
    EXPECT_EQ(d.label, isothermModelNames[i]) << d.tag;
    IsothermType back{};
    ASSERT_TRUE(KindFromString(d.tag, back)) << d.tag;
    EXPECT_EQ(back, t);
    auto iso = CreateIsotherm(t, "c");
    ASSERT_TRUE(iso) << d.tag;
    EnsureIsothermComponent(iso.get(), "c");
    ExpectMatches(d, IsothermFields(iso.get(), "c"), "isotherm." + d.tag);
    // A pair shares its layout (names and units) both ways.
    if (d.has_pair) {
      const auto& p = IsothermKind(d.pair);
      EXPECT_TRUE(p.has_pair && p.pair == t) << d.tag;
      ASSERT_EQ(Names(p), Names(d)) << d.tag;
      for (size_t j = 0; j < d.parameters.size(); ++j) EXPECT_EQ(p.parameters[j].unit, d.parameters[j].unit);
    }
    // IAST: every kind but the water-coupled ones has a pure form.
    EXPECT_EQ(d.iast == IastFamily::None, d.water_coupled) << d.tag;
  }
  // The fitter's list keeps its order.
  using T = IsothermType;
  const std::vector<T> fit = {T::HENRY, T::SINGLE_SITE_LANGMUIR_P, T::DUAL_SITE_LANGMUIR_P, T::SINGLE_SITE_LANGMUIR_C,
                              T::DUAL_SITE_LANGMUIR_C, T::EXTENDED_SINGLE_SITE_LANGMUIR_P,
                              T::EXTENDED_DUAL_SITE_LANGMUIR_P, T::FREUNDLICH, T::SIPS, T::DUAL_SITE_SIPS, T::TOTH,
                              T::SIPS_2, T::EXTENDED_SIPS_2, T::TOTH_2, T::GAB, T::GAB_YOUNG, T::TOTH_WADST,
                              T::TOTH_HUMID};
  EXPECT_EQ(isofit::FittableKinds(), fit);
}

TEST(ModelKinds, FitterDefaultSpecsNameTheDescriptorParameters) {
  // Langmuir-like points at three temperatures, dry and humid.
  std::vector<isofit::Point> pts;
  for (double T : {288.0, 298.0, 308.0})
    for (double P : {1e3, 5e3, 2e4, 5e4, 1e5})
      for (double qw : {0.0, 1.0}) {
        isofit::Point p;
        p.P = P;
        p.T = T;
        p.y = 1.0;
        p.q_w = qw;
        p.q = 3.0 * 1e-4 * P / (1.0 + 1e-4 * P) * (1.0 + 0.1 * qw) * (298.0 / T);
        pts.push_back(p);
      }
  for (IsothermType t : isofit::FittableKinds()) {
    std::vector<isofit::ParamSpec> specs;
    double t_ref = 0.0;
    std::string err;
    ASSERT_TRUE(isofit::DefaultSpecs(t, pts, specs, t_ref, err)) << IsothermKind(t).tag << ": " << err;
    std::set<std::string> spec_names, desc_names;
    for (const auto& s : specs) spec_names.insert(s.name);
    for (const auto& p : IsothermKind(t).parameters) desc_names.insert(p.name);
    EXPECT_EQ(spec_names, desc_names) << IsothermKind(t).tag;
  }
}

// The registry's entries for a kind are the descriptor's parameters.
TEST(ModelKinds, RegistryEntriesFollowTheDescriptors) {
  UserData d = Perez();
  std::string err;
  const AdsorbentID ads = d.bed.layers_.layers_.front().composition_.begin()->first;
  const std::string c = d.fluid.components_.front();
  const LayerID layer = d.bed.layers_.layers_.front().id_;
  const std::string ap = "adsorbent[" + std::to_string(ads) + "].";
  const std::string lp = "layer[" + std::to_string(layer) + "].";

  for (const auto& k : IsothermKindDescriptors()) {
    ASSERT_TRUE(SetIsothermKind(d, ads, c, k.type, err)) << err;
    ParameterRegistry reg(d);
    size_t n = 0;
    for (const auto& p : reg.All()) n += p.Path().rfind(ap + "isotherm[" + c + "].", 0) == 0;
    EXPECT_EQ(n, k.parameters.size()) << k.tag;
    for (const auto& p : k.parameters) ExpectEntry(reg, ap + "isotherm[" + c + "]." + p.name, p);
  }
  for (auto rate : {MassTransferModelType::LDF, MassTransferModelType::LDF2, MassTransferModelType::DUAL_RATE_LDF}) {
    for (size_t i = 0; i < KindDescriptors(ModelFamily::RateCoefficient).size(); ++i) {
      ASSERT_TRUE(SetMassTransferKind(d, ads, c, rate, static_cast<MassTransferCoefficientType>(i), err)) << err;
      ParameterRegistry reg(d);
      // The generic model's leaves, from its descriptor and current switches:
      // a switch per resistance (knudsen with macropore on), the mode while a
      // gas-side term is on, the parameters of the enabled terms and the mode.
      auto& mtm = d.adsorbent_directory.Get(ads).adsorption_source_.mass_transfer_models_.at(c);
      const std::string kp = ap + "kinetics[" + c + "].";
      auto check_site = [&](IMassTransferCoefficientModel* m, bool slow) {
        auto* g = dynamic_cast<MassTransferCoefficientSeriesResistance*>(m);
        ASSERT_NE(g, nullptr);
        const std::string pre = slow ? kSlowSitePrefix : "";
        auto entry = [&](const std::string& name) {
          const KindParameter& p0 = *ResistanceKind().Find(name);
          KindParameter p = rate == MassTransferModelType::LDF2 ? Ldf2RateConstant(p0) : p0;
          if (slow) p = SlowSiteParameter(p);
          ExpectEntry(reg, kp + p.name, p);
        };
        for (const auto& r : ResistanceKind().resistances) {
          const bool listed = r.parent.empty() || g->Enabled(r.parent);
          EXPECT_EQ(reg.Find(kp + pre + "resistances." + r.tag) != nullptr, listed) << r.tag;
          for (const auto& name : r.parameters) {
            if (g->Enabled(r.tag)) entry(name);
            else EXPECT_EQ(reg.Find(kp + pre + name), nullptr) << pre + name;
          }
        }
        EXPECT_EQ(reg.Find(kp + pre + "capacity_factor_mode") != nullptr, g->UsesCapacityFactor());
        for (const auto& name : g->CapacityParameters()) entry(name);
        // The former arrhenius leaf name resolves to the lumped k.
        if (g->Enabled("lumped")) EXPECT_EQ(reg.Find(kp + pre + "pre_exp_factor"), reg.Find(kp + pre + "k"));
      };
      check_site(mtm->mtc_model_.get(), false);
      if (rate == MassTransferModelType::DUAL_RATE_LDF) check_site(mtm->SlowCoefficientModel(), true);
      // The rate model's own parameters (dual_rate_ldf: split).
      for (const auto& p : KindAt(ModelFamily::RateModel, static_cast<size_t>(rate))->parameters)
        ExpectEntry(reg, ap + "kinetics[" + c + "]." + p.name, p);
    }
  }
  auto layer_family = [&](ModelFamily family, auto set) {
    for (size_t i = 0; i < KindDescriptors(family).size(); ++i) {
      ASSERT_TRUE(set(i)) << err;
      ParameterRegistry reg(d);
      for (const auto& p : KindAt(family, i)->parameters) ExpectEntry(reg, lp + p.name, p);
    }
  };
  layer_family(ModelFamily::Dispersion,
               [&](size_t i) { return SetDispersionKind(d, layer, static_cast<AxialDispersionType>(i), err); });
  layer_family(ModelFamily::HeatDispersion,
               [&](size_t i) { return SetHeatDispersionKind(d, layer, static_cast<AxialHeatDispersionType>(i), err); });
  layer_family(ModelFamily::InternalHeatTransfer, [&](size_t i) {
    return SetInternalHeatTransferKind(d, layer, static_cast<InternalHeatTransferType>(i), err);
  });
  layer_family(ModelFamily::ExternalHeatTransfer, [&](size_t i) {
    return SetExternalHeatTransferKind(d, layer, static_cast<ExternalHeatTransferType>(i), err);
  });
  for (size_t i = 0; i < KindDescriptors(ModelFamily::ParticleHeatTransfer).size(); ++i) {
    ASSERT_TRUE(SetParticleHeatTransferKind(d, static_cast<ParticleHeatTransferType>(i), err)) << err;
    ParameterRegistry reg(d);
    for (const auto& p : KindAt(ModelFamily::ParticleHeatTransfer, i)->parameters) ExpectEntry(reg, "bed." + p.name, p);
  }

  // Endpoints: pumps, pump efficiency and power, valves, flows and profiles.
  const size_t step = AddStep(d, "Kinds");
  const size_t top = dynamic_cast<PressureBoundaryCondition*>(d.cycle.FindStepById(step)->top_boundary_condition_.get())
                         ->endpoint_id_;
  const std::string tp = "endpoint[" + std::to_string(top) + "].";
  for (size_t k = 1; k < KindDescriptors(ModelFamily::Pump).size(); ++k) {
    ASSERT_TRUE(SetPumpKind(d, top, k, err)) << err;
    for (size_t e = 0; e < KindDescriptors(ModelFamily::PumpEfficiency).size(); ++e) {
      ASSERT_TRUE(SetPumpEfficiencyKind(d, top, static_cast<PumpEfficiencyType>(e), err)) << err;
      ParameterRegistry reg(d);
      for (const auto& p : KindAt(ModelFamily::Pump, k)->parameters) ExpectEntry(reg, tp + "pump." + p.name, p);
      for (const auto& p : KindAt(ModelFamily::PumpEfficiency, e)->parameters) ExpectEntry(reg, tp + "pump." + p.name, p);
    }
  }
  size_t pump_stream = 0, valve = 0;
  ASSERT_TRUE(AddStream(d, step, BoundaryLocation::BOTTOM, BoundaryStreamType::PUMP, pump_stream, err)) << err;
  ASSERT_TRUE(AddStream(d, step, BoundaryLocation::BOTTOM, BoundaryStreamType::VALVE, valve, err)) << err;
  for (size_t k = 0; k < KindDescriptors(ModelFamily::PumpPower).size(); ++k) {
    ASSERT_TRUE(SetPumpPowerKind(d, pump_stream, static_cast<PumpPowerType>(k), err)) << err;
    ParameterRegistry reg(d);
    for (const auto& p : KindAt(ModelFamily::PumpPower, k)->parameters)
      ExpectEntry(reg, "endpoint[" + std::to_string(pump_stream) + "]." + p.name, p);
  }
  {
    ParameterRegistry reg(d);
    for (const auto& p : FindKind(ModelFamily::Stream, "valve")->parameters)
      ExpectEntry(reg, "endpoint[" + std::to_string(valve) + "].valve." + p.name, p);
  }
  const size_t flow = dynamic_cast<FlowBoundaryCondition*>(d.cycle.FindStepById(step)->bottom_boundary_condition_.get())
                          ->streams_[0]->endpoint_id_;
  const std::string fp = "endpoint[" + std::to_string(flow) + "].";
  for (size_t f = 0; f < KindDescriptors(ModelFamily::Flow).size(); ++f) {
    ASSERT_TRUE(SetFlowKind(d, flow, static_cast<FlowType>(f), err)) << err;
    const KindParameter& q = KindAt(ModelFamily::Flow, f)->parameters.front();
    for (size_t k = 0; k < KindDescriptors(ModelFamily::Profile).size(); ++k) {
      ASSERT_TRUE(SetProfileKind(d, flow, ProfileTarget::Flow, static_cast<BoundaryValueType>(k), err)) << err;
      ParameterRegistry reg(d);
      const KindDescriptor& pk = *KindAt(ModelFamily::Profile, k);
      for (const auto& p : pk.parameters) {
        const bool quantity = p.unit == kProfileQuantityUnit;
        const std::string leaf = k == 0 ? q.name : q.name + "_" + pk.tag + "_" + p.name;
        ExpectEntry(reg, fp + leaf, p, k == 0 ? q.Label() : q.Label() + " (" + p.Label() + ")",
                    quantity ? q.unit : p.unit, quantity ? q.lower : p.lower, quantity ? q.upper : p.upper);
      }
    }
  }
}

}  // namespace adsim
