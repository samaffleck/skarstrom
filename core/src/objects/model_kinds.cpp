/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "adsim/objects/model_kinds.h"

#include <stdexcept>

#include "adsim/boundary_conditions/boundary_value.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/equations/mass_transfer_model.h"
#include "adsim/objects/compressor.h"

namespace adsim {

const KindParameter* KindDescriptor::Find(const std::string& name) const {
  for (const auto& p : parameters)
    if (p.name == name) return &p;
  return nullptr;
}

namespace {

using P = KindParameter;
using K = KindDescriptor;

// A rate-coefficient parameter's tortuosity / macropore void fraction /
// molecular diffusivity, shared by the macropore kinds.
P Tortuosity() { return {"tortuosity", "-", 3.0, 1.0, 10.0}; }
P MacroporeVoid(double d) { return {"macropore_void_fraction", "-", d, 0.01, 0.8, true, "macropore void fraction"}; }
P MolecularDiffusivity() {
  return {"molecular_diffusion_coefficient", "m2/s", 0.0, 0.0, 1e-3, true, "constant molecular diffusivity"};
}

const std::vector<K>& RateModels() {
  static const std::vector<K> k{{"ldf", "Linear Driving Force (1st Order)", {}},
                                {"ldf2", "Linear Driving Force (2nd Order)", {}}};
  return k;
}

const std::vector<K>& RateCoefficients() {
  static const std::vector<K> k{
      {"constant", "Constant", {{"k", "1/s", 0.5, 1e-4, 100.0, true, "mass transfer coefficient"}}},
      {"macropore", "Macropore Diffusion", {Tortuosity(), MacroporeVoid(0.35), MolecularDiffusivity()}},
      {"macropore_knudsen",
       "Macropore + Knudsen Diffusion",
       {{"pore_radius", "m", 1e-9, 1e-10, 1e-6, true, "pore radius"}, Tortuosity(), MacroporeVoid(0.2),
        MolecularDiffusivity()}},
      {"haynes", "Haynes (1973)", {Tortuosity(), MacroporeVoid(0.2)}},
      {"arrhenius",
       "Arrhenius",
       {{"pre_exp_factor", "1/s", 120.0, 0.01, 1e4, true, "MTC pre-exponential factor"},
        {"activation_energy", "J/mol", 1e4, 0.0, 1e5, true, "MTC activation energy"}}},
      {"macropore_capacity",
       "Macropore Diffusion (capacity factor)",
       {Tortuosity(), MacroporeVoid(0.35),
        {"concentration_floor", "mol/m3", 1e-3, 0.0, 10.0, true, "capacity-factor concentration floor"},
        MolecularDiffusivity()}},
  };
  return k;
}

const std::vector<K>& Momentum() {
  static const std::vector<K> k{{"kozeny_carman", "Darcy (Kozeny-Carman)", {}},
                                {"ergun", "Ergun", {}},
                                {"blake_kozeny", "Darcy (Blake-Kozeny)", {}}};
  return k;
}

const std::vector<K>& Dispersion() {
  static const std::vector<K> k{
      {"constant", "Constant", {{"dispersion_coefficient", "m2/s", 0.0, 1e-8, 1e-2, true, "dispersion coefficient"}}},
      {"ruthven",
       "Ruthven (1984)",
       {{"molecular_diffusion_coefficient", "m2/s", 1.6e-5, 1e-8, 1e-4, true, "molecular diffusion coefficient"}}},
      {"taylor_aris",
       "Taylor-Aris (channels)",
       {{"molecular_diffusion_coefficient", "m2/s", 1.6e-5, 1e-8, 1e-4, true, "molecular diffusion coefficient"}}},
  };
  return k;
}

const std::vector<K>& HeatDispersion() {
  // The constant model's coefficient is the effective axial thermal
  // conductivity of the energy balance, so W/(m K), not a diffusivity.
  static const std::vector<K> k{
      {"constant",
       "Constant",
       {{"heat_dispersion_coefficient", "W/(m-K)", 0.0, 1e-3, 100.0, true, "axial thermal conductivity"}}},
      {"ruthven", "Ruthven (1984)", {}},
      {"yagi_kunii", "Yagi and Kunii (1957)", {}},
      {"channel", "Channel conduction", {}},
  };
  return k;
}

const std::vector<K>& InternalHeatTransfer() {
  // Leva: Nu = eta1 Re^eta2 exp(-exp_coeff dp / D_t). The bounds span Leva's
  // canonical (0.813, 0.90) and Casas 2012's fitted (41.13, 0.32) values.
  static const std::vector<K> k{
      {"constant",
       "Constant",
       {{"internal_htc", "W/m2-K", 0.0, 1.0, 1e4, true, "internal heat transfer coefficient"}}},
      {"dixon", "Dixon (1996)", {}},
      {"leva",
       "Leva (1947)",
       {{"leva_eta1", "-", 0.813, 0.01, 100.0, true, "Leva Re coefficient (eta1)"},
        {"leva_eta2", "-", 0.90, 0.1, 1.5, true, "Leva Re exponent (eta2)"},
        {"leva_exp_coeff", "-", 6.0, 0.0, 20.0, true, "Leva wall-decay coefficient (exp_coeff)"}}},
  };
  return k;
}

const std::vector<K>& ExternalHeatTransfer() {
  static const std::vector<K> k{
      {"constant",
       "Constant",
       {{"external_htc", "W/m2-K", 0.0, 1.0, 1e4, true, "external heat transfer coefficient"}}},
      {"holman",
       "Holman (2008)",
       {{"holman_c", "-", 0.027, 0.001, 1.0, true, "Holman coefficient (c)"},
        {"holman_m", "-", 0.8, 0.1, 1.0, true, "Holman Re exponent (m)"},
        {"holman_n", "-", 0.3333333333, 0.1, 1.0, true, "Holman Pr exponent (n)"},
        {"external_velocity", "m/s", 0.5, 0.0, 10.0, true, "external velocity"}}},
  };
  return k;
}

const std::vector<K>& ParticleHeatTransfer() {
  static const std::vector<K> k{
      {"constant",
       "Constant",
       {{"particle_htc", "W/m2-K", 100.0, 1.0, 1e5, true, "gas-to-particle heat transfer coefficient"}}},
      {"film", "Film", {}},
  };
  return k;
}

const std::vector<K>& Boundaries() {
  static const std::vector<K> k{{"pressure", "Pressure", {}}, {"flow", "Flow", {}}, {"closed", "Closed", {}}};
  return k;
}

const std::vector<K>& Streams() {
  static const std::vector<K> k{
      {"flow", "Set flow", {}},
      {"valve",
       "Valve",
       {{"cv", "-", 100.0, 1e-3, 1e3, true, "valve flow coefficient (Cv, US gpm/psi^0.5)"},
        {"xt", "-", 0.53, 0.1, 1.0, true, "valve pressure-drop ratio factor (xT)"}}},
      {"pump", "Pump", {}},
  };
  return k;
}

const std::vector<K>& Profiles() {
  const char* q = kProfileQuantityUnit;
  static const std::vector<K> k{
      {"constant", "Constant", {{"value", q, 1e5, 0.0, 0.0}}},
      {"linear", "Linear", {{"start", q, 1.0, 0.0, 0.0}, {"end", q, 1e5, 0.0, 0.0}, {"time", "s", 1.0, 0.0, 1e4}}},
      {"exponential",
       "Exponential",
       {{"start", q, 1.0, 0.0, 0.0},
        {"end", q, 1e5, 0.0, 0.0},
        {"time", "s", 1.0, 0.0, 1e4},
        {"factor", "-", 100.0, 1.01, 1e6}}},
      {"log_linear",
       "Log-linear",
       {{"start", q, 1.0, 0.0, 0.0}, {"end", q, 1e5, 0.0, 0.0}, {"time", "s", 1.0, 0.0, 1e4}}},
      {"double_exponential",
       "Double exponential",
       {{"start", q, 1.0, 0.0, 0.0},
        {"end", q, 1e5, 0.0, 0.0},
        {"rate_1", "1/s", 0.02, 0.0, 10.0},
        {"rate_2", "1/s", 0.1, 0.0, 10.0},
        {"weight", "-", 0.5, 0.0, 1.0}}},
  };
  return k;
}

const std::vector<K>& Flows() {
  // One parameter each: the flow quantity, with the unit and bounds its
  // profile takes.
  static const std::vector<K> k{
      {"velocity", "Velocity", {{"velocity", "m/s", 0.0, 0.0, 10.0}}},
      {"volumetric_flow", "Volumetric flow", {{"volumetric_flow", "m3/s", 0.0, 0.0, 1.0}}},
      {"molar_flow", "Molar flow", {{"molar_flow", "mol/s", 0.0, 0.0, 100.0}}},
      {"mass_flow", "Mass flow", {{"mass_flow", "kg/s", 0.0, 0.0, 10.0}}},
  };
  return k;
}

const std::vector<K>& Pumps() {
  static const std::vector<K> k{
      {"none", "None", {}},
      {"compressor",
       "Compressor",
       {{"suction_pressure", "Pa", 1e5, 1e3, 1e7, true, "compressor suction pressure"}}},
      {"vacuum_pump",
       "Vacuum Pump",
       {{"discharge_pressure", "Pa", 1e5, 1e3, 1e7, true, "vacuum pump discharge pressure"}}},
  };
  return k;
}

const std::vector<K>& PumpEnergies() {
  static const std::vector<K> k{{"isothermal", "Isothermal", {}}, {"adiabatic", "Adiabatic", {}}};
  return k;
}

const std::vector<K>& PumpEfficiencies() {
  static const std::vector<K> k{
      {"constant", "Constant", {{"efficiency", "-", 0.7, 0.01, 1.0, true, "pump efficiency"}}},
      {"vacuum_pump",
       "Vacuum Pump",
       {{"eta_max", "-", 0.8, 0.01, 1.0, true, "vacuum pump peak efficiency"},
        {"alpha", "-", 19.55, 0.0, 100.0, true, "vacuum pump efficiency decay"}}},
      {"polynomial",
       "Polynomial",
       {{"poly_a", "-", 0.0, -1e3, 1e3, true, "pump efficiency polynomial a"},
        {"poly_b", "1/Pa", 0.0, -1e3, 1e3, true, "pump efficiency polynomial b (x P)"},
        {"poly_c", "1/Pa^2", 0.0, -1e3, 1e3, true, "pump efficiency polynomial c (x P^2)"},
        {"poly_d", "1/Pa^3", 0.0, -1e3, 1e3, true, "pump efficiency polynomial d (x P^3)"},
        {"poly_e", "1/Pa^4", 0.0, -1e3, 1e3, true, "pump efficiency polynomial e (x P^4)"},
        {"poly_f", "1/Pa^5", 0.0, -1e3, 1e3, true, "pump efficiency polynomial f (x P^5)"}}},
  };
  return k;
}

const std::vector<K>& PumpPowers() {
  static const std::vector<K> k{
      {"constant", "Constant", {{"power", "W", 0.0, 0.0, 1e6, true, "pump shaft power"}}},
      {"polynomial",
       "Polynomial",
       {{"power_poly_a", "W", 0.0, -1e6, 1e6, true, "pump power polynomial a"},
        {"power_poly_b", "W/Pa", 0.0, -1e3, 1e3, true, "pump power polynomial b (x P)"},
        {"power_poly_c", "W/Pa^2", 0.0, -1.0, 1.0, true, "pump power polynomial c (x P^2)"},
        {"power_poly_d", "W/Pa^3", 0.0, -1e-3, 1e-3, true, "pump power polynomial d (x P^3)"},
        {"power_poly_e", "W/Pa^4", 0.0, -1e-6, 1e-6, true, "pump power polynomial e (x P^4)"},
        {"power_poly_f", "W/Pa^5", 0.0, -1e-9, 1e-9, true, "pump power polynomial f (x P^5)"}}},
  };
  return k;
}

}  // namespace

const std::vector<ModelFamily>& AllModelFamilies() {
  using F = ModelFamily;
  static const std::vector<F> f{F::RateModel,           F::RateCoefficient,      F::Momentum,
                                F::Dispersion,          F::HeatDispersion,       F::InternalHeatTransfer,
                                F::ExternalHeatTransfer, F::ParticleHeatTransfer, F::Boundary,
                                F::Stream,              F::Profile,              F::Flow,
                                F::Pump,                F::PumpEnergy,           F::PumpEfficiency,
                                F::PumpPower};
  return f;
}

const char* ModelFamilyName(ModelFamily family) {
  switch (family) {
    case ModelFamily::RateModel: return "rate_model";
    case ModelFamily::RateCoefficient: return "rate_coefficient";
    case ModelFamily::Momentum: return "momentum";
    case ModelFamily::Dispersion: return "dispersion";
    case ModelFamily::HeatDispersion: return "heat_dispersion";
    case ModelFamily::InternalHeatTransfer: return "internal_heat_transfer";
    case ModelFamily::ExternalHeatTransfer: return "external_heat_transfer";
    case ModelFamily::ParticleHeatTransfer: return "particle_heat_transfer";
    case ModelFamily::Boundary: return "boundary";
    case ModelFamily::Stream: return "stream";
    case ModelFamily::Profile: return "profile";
    case ModelFamily::Flow: return "flow";
    case ModelFamily::Pump: return "pump";
    case ModelFamily::PumpEnergy: return "pump_energy";
    case ModelFamily::PumpEfficiency: return "pump_efficiency";
    case ModelFamily::PumpPower: return "pump_power";
  }
  return "unknown";
}

const std::vector<KindDescriptor>& KindDescriptors(ModelFamily family) {
  switch (family) {
    case ModelFamily::RateModel: return RateModels();
    case ModelFamily::RateCoefficient: return RateCoefficients();
    case ModelFamily::Momentum: return Momentum();
    case ModelFamily::Dispersion: return Dispersion();
    case ModelFamily::HeatDispersion: return HeatDispersion();
    case ModelFamily::InternalHeatTransfer: return InternalHeatTransfer();
    case ModelFamily::ExternalHeatTransfer: return ExternalHeatTransfer();
    case ModelFamily::ParticleHeatTransfer: return ParticleHeatTransfer();
    case ModelFamily::Boundary: return Boundaries();
    case ModelFamily::Stream: return Streams();
    case ModelFamily::Profile: return Profiles();
    case ModelFamily::Flow: return Flows();
    case ModelFamily::Pump: return Pumps();
    case ModelFamily::PumpEnergy: return PumpEnergies();
    case ModelFamily::PumpEfficiency: return PumpEfficiencies();
    case ModelFamily::PumpPower: return PumpPowers();
  }
  throw std::invalid_argument("unknown model family");
}

const std::vector<std::string>& KindTags(ModelFamily family) {
  // One cached list per family (built on first use; the tables are static).
  static const std::vector<std::vector<std::string>> tags = [] {
    std::vector<std::vector<std::string>> out;
    for (ModelFamily f : AllModelFamilies()) {
      std::vector<std::string> t;
      for (const auto& d : KindDescriptors(f)) t.push_back(d.tag);
      out.push_back(std::move(t));
    }
    return out;
  }();
  return tags.at(static_cast<std::size_t>(family));
}

const KindDescriptor* FindKind(ModelFamily family, const std::string& tag) {
  for (const auto& d : KindDescriptors(family))
    if (d.tag == tag) return &d;
  return nullptr;
}

const KindDescriptor* KindAt(ModelFamily family, std::size_t index) {
  const auto& all = KindDescriptors(family);
  return index < all.size() ? &all[index] : nullptr;
}

KindParameter Ldf2RateConstant(const KindParameter& p) {
  KindParameter out = p;
  if (p.unit != "1/s") return out;
  out.unit = kLdf2RateConstantUnit;
  out.lower = p.lower / 10.0;
  out.upper = p.upper * 10.0;
  return out;
}

// ---- isotherms -----------------------------------------------------------------------

namespace {

// Builder of one isotherm descriptor.
IsothermKindDescriptor Iso(IsothermType type, const char* tag, const char* label, IastFamily iast,
                           std::vector<KindParameter> params) {
  IsothermKindDescriptor d;
  d.type = type;
  d.tag = tag;
  d.label = label;
  d.iast = iast;
  d.parameters = std::move(params);
  return d;
}

std::vector<KindParameter> LangmuirP(double b1) {
  return {{"m1", "mol/kg", 3.0, 0.01, 50.0}, {"b1", "1/Pa", b1, 1e-15, 1e-3}, {"q1", "J/mol", 4e4, 0.0, 1e5}};
}
std::vector<KindParameter> LangmuirC(double b1) {
  return {{"m1", "mol/kg", 3.0, 0.01, 50.0}, {"b1", "m3/mol", b1, 1e-10, 1e-3}, {"q1", "J/mol", 4e4, 0.0, 1e5}};
}
std::vector<KindParameter> DualP() {
  return {{"m1", "mol/kg", 3.0, 0.01, 50.0}, {"b1", "1/Pa", 1e-12, 1e-15, 1e-3}, {"q1", "J/mol", 4.7e4, 0.0, 1e5},
          {"m2", "mol/kg", 3.0, 0.01, 50.0}, {"b2", "1/Pa", 1e-13, 1e-15, 1e-3}, {"q2", "J/mol", 4.2e4, 0.0, 1e5}};
}
std::vector<KindParameter> DualC() {
  return {{"m1", "mol/kg", 3.0, 0.01, 50.0}, {"b1", "m3/mol", 1e-8, 1e-10, 1e-3}, {"q1", "J/mol", 4e4, 0.0, 1e5},
          {"m2", "mol/kg", 3.0, 0.01, 50.0}, {"b2", "m3/mol", 1e-7, 1e-10, 1e-3}, {"q2", "J/mol", 4e4, 0.0, 1e5}};
}
std::vector<KindParameter> Sips2() {
  // b_i = B1 exp(B2/RT) multiplies the partial pressure in Pa.
  return {{"A1", "mol/kg", 4.0, 0.01, 50.0}, {"A2", "-", -0.6, -5.0, 5.0},  {"B1", "1/Pa", 3.6e-4, 1e-13, 1e-5},
          {"B2", "J/mol", 15.0e3, 0.0, 1e5}, {"C1", "-", 1.0, 0.1, 5.0},    {"C2", "-", 0.0, -5.0, 5.0},
          {"T_ref", "K", 298.0, 273.0, 473.0}};
}

std::vector<IsothermKindDescriptor> BuildIsothermKinds() {
  using T = IsothermType;
  using F = IastFamily;
  std::vector<IsothermKindDescriptor> k(static_cast<std::size_t>(T::TOTH_HUMID) + 1);
  auto put = [&](IsothermKindDescriptor d) { k[static_cast<std::size_t>(d.type)] = std::move(d); };

  put(Iso(T::SINGLE_SITE_LANGMUIR_P, "single_site_langmuir_p", "Single Site Langmuir", F::Langmuir, LangmuirP(1e-11)));
  put(Iso(T::DUAL_SITE_LANGMUIR_P, "dual_site_langmuir_p", "Dual Site Langmuir", F::DualSiteLangmuir, DualP()));
  put(Iso(T::INERT, "inert", "Inert", F::Inert, {}));
  put(Iso(T::HENRY, "henry", "Henry", F::Henry, {{"k", "-", 1.0e-12, 1e-15, 1e-3}, {"Q", "J/mol", 4e4, 0.0, 1e5}}));
  put(Iso(T::SINGLE_SITE_LANGMUIR_C, "single_site_langmuir_c", "Single Site Langmuir (Concentration)", F::Langmuir,
          LangmuirC(1e-7)));
  put(Iso(T::DUAL_SITE_LANGMUIR_C, "dual_site_langmuir_c", "Dual Site Langmuir (Concentration)",
          F::DualSiteLangmuir, DualC()));
  put(Iso(T::EXTENDED_SINGLE_SITE_LANGMUIR_P, "extended_single_site_langmuir_p", "Extended Single Site Langmuir",
          F::Langmuir, LangmuirP(1e-11)));
  put(Iso(T::EXTENDED_DUAL_SITE_LANGMUIR_P, "extended_dual_site_langmuir_p", "Extended Dual Site Langmuir",
          F::DualSiteLangmuir, DualP()));
  put(Iso(T::EXTENDED_SINGLE_SITE_LANGMUIR_C, "extended_single_site_langmuir_c",
          "Extended Single Site Langmuir (Concentration)", F::Langmuir, LangmuirC(1e-7)));
  put(Iso(T::EXTENDED_DUAL_SITE_LANGMUIR_C, "extended_dual_site_langmuir_c",
          "Extended Dual Site Langmuir (Concentration)", F::DualSiteLangmuir, DualC()));
  put(Iso(T::ANTI_LANGMUIR, "anti_langmuir", "Anti-Langmuir", F::AntiLangmuir,
          {{"a", "-", 1.0e-6, 1e-10, 1.0}, {"b", "-", 9.0e-6, 1e-10, 1.0}}));
  put(Iso(T::FREUNDLICH, "freundlich", "Freundlich", F::Freundlich,
          {{"a", "-", 1.0, 0.01, 100.0}, {"v", "-", 2.0, 0.1, 10.0}}));
  put(Iso(T::SIPS, "sips", "Sips", F::Sips,
          {{"q_sat", "mol/kg", 2.0, 0.01, 50.0}, {"b", "1/Pa", 1.0e-5, 1e-10, 1.0}, {"v", "-", 1.0, 0.1, 5.0}}));
  put(Iso(T::SIPS_2, "sips_2", "Sips 2", F::Sips, Sips2()));
  put(Iso(T::EXTENDED_SIPS_2, "extended_sips_2", "Extended Sips 2", F::Sips, Sips2()));
  put(Iso(T::DUAL_SITE_SIPS, "dual_site_sips", "Dual Site Sips", F::DualSiteSips,
          {{"q_sat1", "mol/kg", 2.0, 0.01, 50.0},
           {"b1", "1/Pa", 1.0e-5, 1e-10, 1.0},
           {"v1", "-", 1.0, 0.1, 5.0},
           {"q_sat2", "mol/kg", 2.0, 0.01, 50.0},
           {"b2", "1/Pa", 1.0e-5, 1e-10, 1.0},
           {"v2", "-", 1.0, 0.1, 5.0}}));
  put(Iso(T::REDLICH_PETERSON, "redlich_peterson", "Redlich-Peterson", F::RedlichPeterson,
          {{"a", "-", 2.0, 0.01, 100.0}, {"b", "1/Pa", 1.0e-5, 1e-10, 1.0}, {"v", "-", 1.0, 0.1, 5.0}}));
  put(Iso(T::TOTH, "toth", "Toth", F::Toth,
          {{"q_sat", "mol/kg", 2.0, 0.01, 50.0}, {"b", "1/Pa", 1.0e-5, 1e-10, 1.0}, {"v", "-", 1.0, 0.1, 5.0}}));
  put(Iso(T::TOTH_2, "toth_2", "Toth 2", F::Toth,
          {{"q_sat", "mol/kg", 3.715, 0.01, 50.0},
           {"x", "-", 0.0, -1.0, 1.0},
           {"T_0", "K", 324.0, 273.0, 473.0},
           {"b0", "1/Pa", 0.027, 1e-6, 1.0},
           {"dH0", "J/mol", 100000.0, 0.0, 1e6},
           {"w0", "-", 0.28, 0.01, 5.0},
           {"alpha", "-", 0.446, 0.0, 2.0}}));
  put(Iso(T::ARRHENIUS_SIPS, "arrhenius_sips", "Arrhenius Sips (Schell 2012)", F::Sips,
          {{"n_inf_a", "mol/kg", 1.0, 0.01, 50.0},
           {"n_inf_b", "J/mol", 0.0, 0.0, 1e5},
           {"k_A", "1/Pa", 1.0e-8, 1e-15, 1e-3},
           {"k_B", "J/mol", 0.0, 0.0, 1e5},
           {"c_ref", "-", 1.0, 0.1, 5.0},
           {"c_amp", "-", 0.0, -1.0, 1.0},
           {"c_beta", "1/K", 0.0, 0.0, 1.0},
           {"c_Tref", "K", 298.15, 273.0, 473.0}}));
  put(Iso(T::EXTENDED_LANGMUIR_FREUNDLICH, "extended_langmuir_freundlich", "Extended Langmuir-Freundlich", F::Sips,
          {{"k1", "mol/kg", 4.0, 0.01, 50.0},
           {"k2", "mol/(kg-K)", 0.0, -0.1, 0.1},
           {"k3", "-", 1.0, 1e-8, 1e3},
           {"k4", "K", 0.0, 0.0, 1e4},
           {"k5", "-", 1.0, 0.05, 5.0},
           {"k6", "K", 0.0, -1e3, 1e3},
           {"P_ref", "Pa", 1.0e5, 1.0, 1e7}}));
  put(Iso(T::GAB, "gab", "GAB", F::Gab,
          {{"q_m0", "mol/kg", 3.0, 0.01, 50.0},
           {"dH_m", "J/mol", 0.0, -1e5, 1e5},
           {"C0", "-", 1.0, 1e-6, 1e6},
           {"dH_C", "J/mol", 5000.0, -1e5, 1e5},
           {"K0", "-", 0.7, 0.01, 2.0},
           {"dH_K", "J/mol", 0.0, -1e5, 1e5}}));
  put(Iso(T::GAB_YOUNG, "gab_young", "GAB (Young 2021 energies)", F::Gab,
          {{"q_m", "mol/kg", 3.63, 0.01, 50.0},
           {"C_E", "J/mol", 47110.0, 0.0, 2e5},
           {"D", "1/K", 0.023744, 0.0, 0.1},
           {"F", "J/mol", 57706.0, 0.0, 2e5},
           {"G", "J/(mol-K)", -47.814, -500.0, 500.0},
           {"E10_A", "J/mol", 57220.0, 0.0, 2e5},
           {"E10_B", "J/(mol-K)", -44.38, -500.0, 500.0}}));
  put(Iso(T::TOTH_WADST, "toth_wadst", "Toth WADST (CO2 with water, Young 2021)", F::None,
          {{"T_0", "K", 298.15, 273.0, 473.0},
           {"q_inf_dry", "mol/kg", 4.86, 0.01, 50.0},
           {"chi_dry", "-", 0.0, -5.0, 5.0},
           {"b0_dry", "1/Pa", 2.85e-21, 1e-30, 1.0},
           {"dH_dry", "J/mol", 117798.0, 0.0, 3e5},
           {"tau0_dry", "-", 0.209, 0.01, 5.0},
           {"alpha_dry", "-", 0.523, -5.0, 5.0},
           {"q_inf_wet", "mol/kg", 9.035, 0.01, 50.0},
           {"chi_wet", "-", 0.0, -5.0, 5.0},
           {"b0_wet", "1/Pa", 1.230e-18, 1e-30, 1.0},
           {"dH_wet", "J/mol", 203687.0, 0.0, 3e5},
           {"tau0_wet", "-", 0.053, 0.01, 5.0},
           {"alpha_wet", "-", 0.053, -5.0, 5.0},
           {"A", "mol/kg", 1.532, 0.0, 50.0}}));
  put(Iso(T::TOTH_HUMID, "toth_humid", "Toth humid (CO2 with water, Stampi-Bombelli 2020)", F::None,
          {{"q_sat", "mol/kg", 2.38, 0.01, 50.0},
           {"x", "-", 0.0, -5.0, 5.0},
           {"T_0", "K", 296.0, 273.0, 473.0},
           {"b0", "1/Pa", 0.0707, 1e-12, 10.0},
           {"dH0", "J/mol", 57047.0, 0.0, 3e5},
           {"w0", "-", 0.4148, 0.01, 5.0},
           {"alpha", "-", -1.606, -5.0, 5.0},
           {"gamma", "kg/mol", 0.0016, -5.0, 5.0},
           {"beta", "kg/mol", 59.1, -10.0, 200.0}}));

  for (T t : {T::EXTENDED_SINGLE_SITE_LANGMUIR_P, T::EXTENDED_SINGLE_SITE_LANGMUIR_C, T::EXTENDED_DUAL_SITE_LANGMUIR_P,
              T::EXTENDED_DUAL_SITE_LANGMUIR_C, T::EXTENDED_SIPS_2, T::ARRHENIUS_SIPS, T::EXTENDED_LANGMUIR_FREUNDLICH})
    k[static_cast<std::size_t>(t)].competitive = true;
  for (T t : {T::TOTH_WADST, T::TOTH_HUMID}) k[static_cast<std::size_t>(t)].water_coupled = true;
  // The isotherm fitter's kinds, and those whose temperature dependence it fits.
  for (T t : {T::HENRY, T::SINGLE_SITE_LANGMUIR_P, T::DUAL_SITE_LANGMUIR_P, T::SINGLE_SITE_LANGMUIR_C,
              T::DUAL_SITE_LANGMUIR_C, T::EXTENDED_SINGLE_SITE_LANGMUIR_P, T::EXTENDED_DUAL_SITE_LANGMUIR_P,
              T::FREUNDLICH, T::SIPS, T::DUAL_SITE_SIPS, T::TOTH, T::SIPS_2, T::EXTENDED_SIPS_2, T::TOTH_2, T::GAB,
              T::GAB_YOUNG, T::TOTH_WADST, T::TOTH_HUMID})
    k[static_cast<std::size_t>(t)].fittable = true;
  for (T t : {T::HENRY, T::SINGLE_SITE_LANGMUIR_P, T::DUAL_SITE_LANGMUIR_P, T::SINGLE_SITE_LANGMUIR_C,
              T::DUAL_SITE_LANGMUIR_C, T::EXTENDED_SINGLE_SITE_LANGMUIR_P, T::EXTENDED_DUAL_SITE_LANGMUIR_P,
              T::SIPS_2, T::EXTENDED_SIPS_2, T::TOTH_2, T::GAB, T::GAB_YOUNG, T::TOTH_WADST, T::TOTH_HUMID})
    k[static_cast<std::size_t>(t)].temperature_dependent = true;
  // A pure kind and its competitive form share their parameter layout.
  const std::pair<T, T> pairs[] = {
      {T::SINGLE_SITE_LANGMUIR_P, T::EXTENDED_SINGLE_SITE_LANGMUIR_P},
      {T::SINGLE_SITE_LANGMUIR_C, T::EXTENDED_SINGLE_SITE_LANGMUIR_C},
      {T::DUAL_SITE_LANGMUIR_P, T::EXTENDED_DUAL_SITE_LANGMUIR_P},
      {T::DUAL_SITE_LANGMUIR_C, T::EXTENDED_DUAL_SITE_LANGMUIR_C},
      {T::SIPS_2, T::EXTENDED_SIPS_2},
  };
  for (const auto& [a, b] : pairs) {
    auto& da = k[static_cast<std::size_t>(a)];
    auto& db = k[static_cast<std::size_t>(b)];
    da.has_pair = db.has_pair = true;
    da.pair = b;
    db.pair = a;
  }
  return k;
}

}  // namespace

const std::vector<IsothermKindDescriptor>& IsothermKindDescriptors() {
  static const std::vector<IsothermKindDescriptor> k = BuildIsothermKinds();
  return k;
}

const IsothermKindDescriptor& IsothermKind(IsothermType type) {
  const auto& all = IsothermKindDescriptors();
  const auto i = static_cast<std::size_t>(type);
  if (i >= all.size()) throw std::invalid_argument("unknown isotherm kind " + std::to_string(i));
  return all[i];
}

// ---- field bindings ---------------------------------------------------------------------

ModelFieldList ModelFields(IMassTransferCoefficientModel* m) {
  ModelFieldList out;
  if (auto* c = dynamic_cast<MassTransferCoefficientConstant*>(m)) {
    out = {{"k", &c->constant_mtc_}};
  } else if (auto* c = dynamic_cast<MassTransferCoefficientMacropore*>(m)) {
    out = {{"tortuosity", &c->tortuosity_},
           {"macropore_void_fraction", &c->macropore_void_fraction_},
           {"molecular_diffusion_coefficient", &c->molecular_diffusion_coefficient_}};
  } else if (auto* c = dynamic_cast<MassTransferCoefficientMacroporeKnudsen*>(m)) {
    out = {{"pore_radius", &c->pore_radius_},
           {"tortuosity", &c->tortuosity_},
           {"macropore_void_fraction", &c->macropore_void_fraction_},
           {"molecular_diffusion_coefficient", &c->molecular_diffusion_coefficient_}};
  } else if (auto* c = dynamic_cast<MassTransferCoefficient_Haynes*>(m)) {
    out = {{"tortuosity", &c->tortuosity_}, {"macropore_void_fraction", &c->macropore_void_fraction_}};
  } else if (auto* c = dynamic_cast<MassTransferCoefficientArrhenius*>(m)) {
    out = {{"pre_exp_factor", &c->pre_exp_factor_}, {"activation_energy", &c->activation_energy_}};
  } else if (auto* c = dynamic_cast<MassTransferCoefficientMacroporeCapacity*>(m)) {
    out = {{"tortuosity", &c->tortuosity_},
           {"macropore_void_fraction", &c->macropore_void_fraction_},
           {"concentration_floor", &c->concentration_floor_},
           {"molecular_diffusion_coefficient", &c->molecular_diffusion_coefficient_}};
  }
  return out;
}

ModelFieldList ModelFields(IDispersionModel* m) {
  if (auto* d = dynamic_cast<DispersionModelConstant*>(m)) return {{"dispersion_coefficient", &d->dispersion_coefficient_}};
  if (auto* d = dynamic_cast<DispersionModelRuthven*>(m)) return {{"molecular_diffusion_coefficient", &d->Dm}};
  if (auto* d = dynamic_cast<DispersionModelTaylorAris*>(m)) return {{"molecular_diffusion_coefficient", &d->Dm}};
  return {};
}

ModelFieldList ModelFields(IHeatDispersionModel* m) {
  if (auto* d = dynamic_cast<HeatDispersionModelConstant*>(m))
    return {{"heat_dispersion_coefficient", &d->dispersion_coefficient_}};
  return {};
}

ModelFieldList ModelFields(IInternalHeatTransferModel* m) {
  if (auto* h = dynamic_cast<ConstantInternalHeatTransferModel*>(m))
    return {{"internal_htc", &h->heat_transfer_coefficient_}};
  if (auto* h = dynamic_cast<InternalHeatTransferModel_Leva*>(m))
    return {{"leva_eta1", &h->eta1}, {"leva_eta2", &h->eta2}, {"leva_exp_coeff", &h->exp_coeff}};
  return {};
}

ModelFieldList ModelFields(IExternalHeatTransferModel* m) {
  if (auto* h = dynamic_cast<ConstantExternalHeatTransferModel*>(m))
    return {{"external_htc", &h->heat_transfer_coefficient_}};
  if (auto* h = dynamic_cast<ExternalHeatTransferModel_Holman*>(m))
    return {{"holman_c", &h->c}, {"holman_m", &h->m}, {"holman_n", &h->n}, {"external_velocity", &h->external_velocity}};
  return {};
}

ModelFieldList ModelFields(IParticleHeatTransferModel* m) {
  if (auto* h = dynamic_cast<ConstantParticleHeatTransferModel*>(m))
    return {{"particle_htc", &h->heat_transfer_coefficient_}};
  return {};
}

ModelFieldList ModelFields(IBoundaryValue* m) {
  if (auto* v = dynamic_cast<ConstantBoundaryValue*>(m)) return {{"value", &v->constant_value_}};
  if (auto* v = dynamic_cast<LinearBoundaryValue*>(m))
    return {{"start", &v->start_value_}, {"end", &v->end_value_}, {"time", &v->ramp_time_}};
  if (auto* v = dynamic_cast<ExponentialBoundaryValue*>(m))
    return {{"start", &v->start_value_}, {"end", &v->end_value_}, {"time", &v->step_time_}, {"factor", &v->factor_}};
  if (auto* v = dynamic_cast<LogLinearBoundaryValue*>(m))
    return {{"start", &v->start_value_}, {"end", &v->end_value_}, {"time", &v->ramp_time_}};
  if (auto* v = dynamic_cast<DoubleExponentialBoundaryValue*>(m))
    return {{"start", &v->start_value_},
            {"end", &v->end_value_},
            {"rate_1", &v->rate_1_},
            {"rate_2", &v->rate_2_},
            {"weight", &v->weight_}};
  return {};
}

ModelFieldList ModelFields(IPumpEfficiencyModel* m) {
  if (auto* e = dynamic_cast<ConstantPumpEfficiencyModel*>(m)) return {{"efficiency", &e->efficiency_}};
  if (auto* e = dynamic_cast<VacuumPumpEfficiencyModel*>(m)) return {{"eta_max", &e->eta_max_}, {"alpha", &e->alpha_}};
  if (auto* e = dynamic_cast<PolyPumpEfficiencyModel*>(m))
    return {{"poly_a", &e->a_}, {"poly_b", &e->b_}, {"poly_c", &e->c_},
            {"poly_d", &e->d_}, {"poly_e", &e->e_}, {"poly_f", &e->f_}};
  return {};
}

ModelFieldList ModelFields(IPumpPowerModel* m) {
  if (auto* p = dynamic_cast<ConstantPumpPowerModel*>(m)) return {{"power", &p->power_}};
  if (auto* p = dynamic_cast<PolyPumpPowerModel*>(m))
    return {{"power_poly_a", &p->poly_.a}, {"power_poly_b", &p->poly_.b}, {"power_poly_c", &p->poly_.c},
            {"power_poly_d", &p->poly_.d}, {"power_poly_e", &p->poly_.e}, {"power_poly_f", &p->poly_.f}};
  return {};
}

ModelFieldList ModelFields(IPump* pump) {
  if (auto* c = dynamic_cast<Compressor*>(pump)) return {{"suction_pressure", &c->suction_pressure_}};
  if (auto* v = dynamic_cast<VacuumPump*>(pump)) return {{"discharge_pressure", &v->discharge_pressure_}};
  return {};
}

std::vector<std::string> CarryOverParameters(const KindDescriptor& from_kind, const ModelFieldList& from,
                                             const KindDescriptor& to_kind, const ModelFieldList& to) {
  std::vector<std::string> copied;
  for (const auto& [name, dst] : to) {
    const KindParameter* tp = to_kind.Find(name);
    const KindParameter* fp = from_kind.Find(name);
    if (!tp || !fp || tp->unit != fp->unit || !dst) continue;
    for (const auto& [from_name, src] : from) {
      if (from_name != name || !src) continue;
      *dst = *src;
      copied.push_back(name);
      break;
    }
  }
  return copied;
}

}  // namespace adsim
