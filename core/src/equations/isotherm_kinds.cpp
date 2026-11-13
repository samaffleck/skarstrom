/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "adsim/equations/isotherm_kinds.h"
#include "adsim/math/constants.h"
#include "adsim/objects/model_kinds.h"

#include <cmath>
#include <limits>

namespace adsim {

// Tags and flags come from the kind descriptors (objects/model_kinds.h).
const char* IsothermKindString(IsothermType t) {
  const auto& all = IsothermKindDescriptors();
  const auto i = static_cast<std::size_t>(t);
  return i < all.size() ? all[i].tag.c_str() : "unknown";
}

std::vector<std::pair<std::string, ExpressionParam*>>
IsothermFields(IIsotherm* iso, const std::string& component) {
  using T = IsothermType;
  std::vector<std::pair<std::string, ExpressionParam*>> out;
  if (!iso) return out;
  auto add = [&](const char* n, ExpressionParam& p) { out.emplace_back(n, &p); };
  // The dual-site Langmuir kinds' optional water factors (DslWaterFactors).
  auto add_water_factors = [&](auto& p) {
    add("lambda1", p.lambda1); add("lambda2", p.lambda2); add("lambda3", p.lambda3); add("gamma2", p.gamma2);
  };

  switch (iso->GetIsothermType()) {
    case T::INERT:
      break;  // supported, no params
    case T::HENRY: {
      auto* c = dynamic_cast<HenryIsotherm*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("k", it->second.k); add("Q", it->second.Q); } }
      break;
    }
    case T::SINGLE_SITE_LANGMUIR_P: {
      auto* c = dynamic_cast<LangmuirIsotherm*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("m1", it->second.m1); add("b1", it->second.b1); add("q1", it->second.q1); } }
      break;
    }
    case T::SINGLE_SITE_LANGMUIR_C: {
      auto* c = dynamic_cast<LangmuirIsotherm_Concentration*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("m1", it->second.m1); add("b1", it->second.b1); add("q1", it->second.q1); } }
      break;
    }
    case T::DUAL_SITE_LANGMUIR_P: {
      auto* c = dynamic_cast<DualSiteLangmuirIsotherm*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("m1", it->second.m1); add("b1", it->second.b1); add("q1", it->second.q1);
        add("m2", it->second.m2); add("b2", it->second.b2); add("q2", it->second.q2);
        add_water_factors(it->second); } }
      break;
    }
    case T::DUAL_SITE_LANGMUIR_C: {
      auto* c = dynamic_cast<DualSiteLangmuirIsotherm_Concentration*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("m1", it->second.m1); add("b1", it->second.b1); add("q1", it->second.q1);
        add("m2", it->second.m2); add("b2", it->second.b2); add("q2", it->second.q2);
        add_water_factors(it->second); } }
      break;
    }
    case T::FREUNDLICH: {
      auto* c = dynamic_cast<FreundlichIsotherm*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("a", it->second.a); add("v", it->second.v); } }
      break;
    }
    case T::SIPS: {
      auto* c = dynamic_cast<SipsIsotherm*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("q_sat", it->second.q_sat); add("b", it->second.b); add("v", it->second.v); } }
      break;
    }
    case T::TOTH: {
      auto* c = dynamic_cast<TothIsotherm*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("q_sat", it->second.q_sat); add("b", it->second.b); add("v", it->second.v); } }
      break;
    }
    case T::ANTI_LANGMUIR: {
      auto* c = dynamic_cast<AntiLangmuirIsotherm*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("a", it->second.a); add("b", it->second.b); } }
      break;
    }
    case T::REDLICH_PETERSON: {
      auto* c = dynamic_cast<RedlichPetersonIsotherm*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("a", it->second.a); add("b", it->second.b); add("v", it->second.v); } }
      break;
    }
    case T::EXTENDED_SINGLE_SITE_LANGMUIR_P: {
      auto* c = dynamic_cast<ExtendedLangmuirIsotherm*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("m1", it->second.m1); add("b1", it->second.b1); add("q1", it->second.q1); } }
      break;
    }
    case T::EXTENDED_SINGLE_SITE_LANGMUIR_C: {
      auto* c = dynamic_cast<ExtendedLangmuirIsotherm_Concentration*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("m1", it->second.m1); add("b1", it->second.b1); add("q1", it->second.q1); } }
      break;
    }
    case T::EXTENDED_DUAL_SITE_LANGMUIR_P: {
      auto* c = dynamic_cast<ExtendedDualSiteLangmuirIsotherm*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("m1", it->second.m1); add("b1", it->second.b1); add("q1", it->second.q1);
        add("m2", it->second.m2); add("b2", it->second.b2); add("q2", it->second.q2);
        add_water_factors(it->second); } }
      break;
    }
    case T::EXTENDED_DUAL_SITE_LANGMUIR_C: {
      auto* c = dynamic_cast<ExtendedDualSiteLangmuirIsotherm_Concentration*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("m1", it->second.m1); add("b1", it->second.b1); add("q1", it->second.q1);
        add("m2", it->second.m2); add("b2", it->second.b2); add("q2", it->second.q2);
        add_water_factors(it->second); } }
      break;
    }
    case T::SIPS_2: {
      auto* c = dynamic_cast<Sips2Isotherm*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("A1", it->second.A1); add("A2", it->second.A2); add("B1", it->second.B1);
        add("B2", it->second.B2); add("C1", it->second.C1); add("C2", it->second.C2);
        add("T_ref", it->second.T_ref); } }
      break;
    }
    case T::EXTENDED_SIPS_2: {
      auto* c = dynamic_cast<ExtendedSips2Isotherm*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("A1", it->second.A1); add("A2", it->second.A2); add("B1", it->second.B1);
        add("B2", it->second.B2); add("C1", it->second.C1); add("C2", it->second.C2);
        add("T_ref", it->second.T_ref); } }
      break;
    }
    case T::DUAL_SITE_SIPS: {
      auto* c = dynamic_cast<DualSiteSipsIsotherm*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("q_sat1", it->second.q_sat1); add("b1", it->second.b1); add("v1", it->second.v1);
        add("q_sat2", it->second.q_sat2); add("b2", it->second.b2); add("v2", it->second.v2); } }
      break;
    }
    case T::TOTH_2: {
      auto* c = dynamic_cast<Toth2Isotherm*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("q_sat", it->second.q_sat); add("x", it->second.x); add("T_0", it->second.T_0);
        add("b0", it->second.b0); add("dH0", it->second.dH0); add("w0", it->second.w0);
        add("alpha", it->second.alpha); } }
      break;
    }
    case T::ARRHENIUS_SIPS: {
      auto* c = dynamic_cast<ArrheniusSipsIsotherm*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("n_inf_a", it->second.n_inf_a); add("n_inf_b", it->second.n_inf_b);
        add("k_A", it->second.k_A); add("k_B", it->second.k_B);
        add("c_ref", it->second.c_ref); add("c_amp", it->second.c_amp);
        add("c_beta", it->second.c_beta); add("c_Tref", it->second.c_Tref); } }
      break;
    }
    case T::EXTENDED_LANGMUIR_FREUNDLICH: {
      auto* c = dynamic_cast<ExtendedLangmuirFreundlichIsotherm*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("k1", it->second.k1); add("k2", it->second.k2); add("k3", it->second.k3);
        add("k4", it->second.k4); add("k5", it->second.k5); add("k6", it->second.k6);
        add("P_ref", it->second.P_ref); } }
      break;
    }
    case T::GAB: {
      auto* c = dynamic_cast<GabIsotherm*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("q_m0", it->second.q_m0); add("dH_m", it->second.dH_m);
        add("C0", it->second.C0); add("dH_C", it->second.dH_C);
        add("K0", it->second.K0); add("dH_K", it->second.dH_K); } }
      break;
    }
    case T::GAB_YOUNG: {
      auto* c = dynamic_cast<GabYoungIsotherm*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("q_m", it->second.q_m); add("C_E", it->second.C_E); add("D", it->second.D);
        add("F", it->second.F); add("G", it->second.G);
        add("E10_A", it->second.E10_A); add("E10_B", it->second.E10_B); } }
      break;
    }
    case T::TOTH_WADST: {
      auto* c = dynamic_cast<TothWadstIsotherm*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        auto& p = it->second;
        add("T_0", p.T_0);
        add("q_inf_dry", p.q_inf_dry); add("chi_dry", p.chi_dry); add("b0_dry", p.b0_dry);
        add("dH_dry", p.dH_dry); add("tau0_dry", p.tau0_dry); add("alpha_dry", p.alpha_dry);
        add("q_inf_wet", p.q_inf_wet); add("chi_wet", p.chi_wet); add("b0_wet", p.b0_wet);
        add("dH_wet", p.dH_wet); add("tau0_wet", p.tau0_wet); add("alpha_wet", p.alpha_wet);
        add("A", p.A); } }
      break;
    }
    case T::TOTH_HUMID: {
      auto* c = dynamic_cast<TothHumidIsotherm*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        auto& p = it->second;
        add("q_sat", p.q_sat); add("x", p.x); add("T_0", p.T_0);
        add("b0", p.b0); add("dH0", p.dH0); add("w0", p.w0);
        add("alpha", p.alpha); add("gamma", p.gamma); add("beta", p.beta); } }
      break;
    }
    case T::CUBIC_LANGMUIR: {
      auto* c = dynamic_cast<CubicLangmuirIsotherm*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        auto& p = it->second;
        add("q_sb", p.q_sb); add("b0_1", p.b0_1); add("b0_2", p.b0_2); add("b0_3", p.b0_3);
        add("dU1", p.dU1); add("q_sd", p.q_sd); add("d0", p.d0); add("dU2", p.dU2); } }
      break;
    }
    // No `default`: -Wswitch flags a new IsothermType missing its param layout.
  }
  return out;
}

bool IsSupportedIsothermKind(IsothermType t) {
  using T = IsothermType;
  switch (t) {
    case T::INERT:
    case T::HENRY:
    case T::SINGLE_SITE_LANGMUIR_P:
    case T::SINGLE_SITE_LANGMUIR_C:
    case T::DUAL_SITE_LANGMUIR_P:
    case T::DUAL_SITE_LANGMUIR_C:
    case T::FREUNDLICH:
    case T::SIPS:
    case T::TOTH:
    case T::ANTI_LANGMUIR:
    case T::REDLICH_PETERSON:
    case T::EXTENDED_SINGLE_SITE_LANGMUIR_P:
    case T::EXTENDED_SINGLE_SITE_LANGMUIR_C:
    case T::EXTENDED_DUAL_SITE_LANGMUIR_P:
    case T::EXTENDED_DUAL_SITE_LANGMUIR_C:
    case T::SIPS_2:
    case T::EXTENDED_SIPS_2:
    case T::DUAL_SITE_SIPS:
    case T::TOTH_2:
    case T::ARRHENIUS_SIPS:
    case T::EXTENDED_LANGMUIR_FREUNDLICH:
    case T::GAB:
    case T::GAB_YOUNG:
    case T::TOTH_WADST:
    case T::TOTH_HUMID:
    case T::CUBIC_LANGMUIR:
      return true;
  }
  return false;  // out-of-range value; every enumerator is listed above (-Wswitch)
}

bool IsCompetitiveIsothermKind(IsothermType t) {
  const auto i = static_cast<std::size_t>(t);
  return i < IsothermKindDescriptors().size() && IsothermKindDescriptors()[i].competitive;
}

bool IsWaterCoupledIsothermKind(IsothermType t) {
  const auto i = static_cast<std::size_t>(t);
  return i < IsothermKindDescriptors().size() && IsothermKindDescriptors()[i].water_coupled;
}

bool IsothermUsesWaterLoading(IIsotherm* iso, const std::string& component) {
  using T = IsothermType;
  if (!iso) return false;
  const T t = iso->GetIsothermType();
  if (IsWaterCoupledIsothermKind(t)) return true;
  auto active = [&](auto* x) {
    if (!x) return false;
    auto it = x->params_.find(component);
    return it != x->params_.end() && DslWaterFactors::Active(it->second);
  };
  switch (t) {
    case T::DUAL_SITE_LANGMUIR_P: return active(dynamic_cast<DualSiteLangmuirIsotherm*>(iso));
    case T::EXTENDED_DUAL_SITE_LANGMUIR_P: return active(dynamic_cast<ExtendedDualSiteLangmuirIsotherm*>(iso));
    case T::DUAL_SITE_LANGMUIR_C: return active(dynamic_cast<DualSiteLangmuirIsotherm_Concentration*>(iso));
    case T::EXTENDED_DUAL_SITE_LANGMUIR_C:
      return active(dynamic_cast<ExtendedDualSiteLangmuirIsotherm_Concentration*>(iso));
    default: return false;
  }
}

bool PairedIsothermKind(IsothermType t, IsothermType& out) {
  const auto i = static_cast<std::size_t>(t);
  if (i >= IsothermKindDescriptors().size() || !IsothermKindDescriptors()[i].has_pair) return false;
  out = IsothermKindDescriptors()[i].pair;
  return true;
}

bool KindFromString(const std::string& kind, IsothermType& out) {
  for (auto t : AllIsothermTypes()) {
    if (kind == IsothermKindString(t)) {
      out = t;
      return true;
    }
  }
  return false;
}

void EnsureIsothermComponent(IIsotherm* iso, const std::string& c) {
  using T = IsothermType;
  if (!iso) return;
  switch (iso->GetIsothermType()) {
    case T::INERT: if (auto* p = dynamic_cast<InertIsotherm*>(iso)) (void)p->params_[c]; break;
    case T::HENRY: if (auto* p = dynamic_cast<HenryIsotherm*>(iso)) (void)p->params_[c]; break;
    case T::SINGLE_SITE_LANGMUIR_P: if (auto* p = dynamic_cast<LangmuirIsotherm*>(iso)) (void)p->params_[c]; break;
    case T::SINGLE_SITE_LANGMUIR_C: if (auto* p = dynamic_cast<LangmuirIsotherm_Concentration*>(iso)) (void)p->params_[c]; break;
    case T::DUAL_SITE_LANGMUIR_P: if (auto* p = dynamic_cast<DualSiteLangmuirIsotherm*>(iso)) (void)p->params_[c]; break;
    case T::DUAL_SITE_LANGMUIR_C: if (auto* p = dynamic_cast<DualSiteLangmuirIsotherm_Concentration*>(iso)) (void)p->params_[c]; break;
    case T::FREUNDLICH: if (auto* p = dynamic_cast<FreundlichIsotherm*>(iso)) (void)p->params_[c]; break;
    case T::SIPS: if (auto* p = dynamic_cast<SipsIsotherm*>(iso)) (void)p->params_[c]; break;
    case T::TOTH: if (auto* p = dynamic_cast<TothIsotherm*>(iso)) (void)p->params_[c]; break;
    case T::ANTI_LANGMUIR: if (auto* p = dynamic_cast<AntiLangmuirIsotherm*>(iso)) (void)p->params_[c]; break;
    case T::REDLICH_PETERSON: if (auto* p = dynamic_cast<RedlichPetersonIsotherm*>(iso)) (void)p->params_[c]; break;
    case T::EXTENDED_SINGLE_SITE_LANGMUIR_P: if (auto* p = dynamic_cast<ExtendedLangmuirIsotherm*>(iso)) (void)p->params_[c]; break;
    case T::EXTENDED_SINGLE_SITE_LANGMUIR_C: if (auto* p = dynamic_cast<ExtendedLangmuirIsotherm_Concentration*>(iso)) (void)p->params_[c]; break;
    case T::EXTENDED_DUAL_SITE_LANGMUIR_P: if (auto* p = dynamic_cast<ExtendedDualSiteLangmuirIsotherm*>(iso)) (void)p->params_[c]; break;
    case T::EXTENDED_DUAL_SITE_LANGMUIR_C: if (auto* p = dynamic_cast<ExtendedDualSiteLangmuirIsotherm_Concentration*>(iso)) (void)p->params_[c]; break;
    case T::SIPS_2: if (auto* p = dynamic_cast<Sips2Isotherm*>(iso)) (void)p->params_[c]; break;
    case T::EXTENDED_SIPS_2: if (auto* p = dynamic_cast<ExtendedSips2Isotherm*>(iso)) (void)p->params_[c]; break;
    case T::DUAL_SITE_SIPS: if (auto* p = dynamic_cast<DualSiteSipsIsotherm*>(iso)) (void)p->params_[c]; break;
    case T::TOTH_2: if (auto* p = dynamic_cast<Toth2Isotherm*>(iso)) (void)p->params_[c]; break;
    case T::ARRHENIUS_SIPS: if (auto* p = dynamic_cast<ArrheniusSipsIsotherm*>(iso)) (void)p->params_[c]; break;
    case T::EXTENDED_LANGMUIR_FREUNDLICH: if (auto* p = dynamic_cast<ExtendedLangmuirFreundlichIsotherm*>(iso)) (void)p->params_[c]; break;
    case T::GAB: if (auto* p = dynamic_cast<GabIsotherm*>(iso)) (void)p->params_[c]; break;
    case T::GAB_YOUNG: if (auto* p = dynamic_cast<GabYoungIsotherm*>(iso)) (void)p->params_[c]; break;
    case T::TOTH_WADST: if (auto* p = dynamic_cast<TothWadstIsotherm*>(iso)) (void)p->params_[c]; break;
    case T::TOTH_HUMID: if (auto* p = dynamic_cast<TothHumidIsotherm*>(iso)) (void)p->params_[c]; break;
    case T::CUBIC_LANGMUIR: if (auto* p = dynamic_cast<CubicLangmuirIsotherm*>(iso)) (void)p->params_[c]; break;
  }
}

double EvalIsotherm(IIsotherm* iso, const std::string& c, double T, double P, double yi, double q_water) {
  using K = IsothermType;
  double q = 0.0;
  if (!iso) return 0.0;
  namespace IF = IsothermFunction;
  switch (iso->GetIsothermType()) {
    case K::INERT: break;
    case K::HENRY: if (auto* x = dynamic_cast<HenryIsotherm*>(iso)) IF::Henry(q, T, P, yi, x->params_[c]); break;
    case K::SINGLE_SITE_LANGMUIR_P: if (auto* x = dynamic_cast<LangmuirIsotherm*>(iso)) IF::Langmuir(q, T, P, yi, x->params_[c]); break;
    case K::EXTENDED_SINGLE_SITE_LANGMUIR_P: if (auto* x = dynamic_cast<ExtendedLangmuirIsotherm*>(iso)) IF::Langmuir(q, T, P, yi, x->params_[c]); break;
    case K::SINGLE_SITE_LANGMUIR_C: if (auto* x = dynamic_cast<LangmuirIsotherm_Concentration*>(iso)) IF::Langmuir_Concentration(q, T, P, yi, x->params_[c]); break;
    case K::EXTENDED_SINGLE_SITE_LANGMUIR_C: if (auto* x = dynamic_cast<ExtendedLangmuirIsotherm_Concentration*>(iso)) IF::Langmuir_Concentration(q, T, P, yi, x->params_[c]); break;
    case K::DUAL_SITE_LANGMUIR_P: if (auto* x = dynamic_cast<DualSiteLangmuirIsotherm*>(iso)) IF::DualSiteLangmuir(q, T, P, yi, q_water, x->params_[c]); break;
    case K::EXTENDED_DUAL_SITE_LANGMUIR_P: if (auto* x = dynamic_cast<ExtendedDualSiteLangmuirIsotherm*>(iso)) IF::DualSiteLangmuir(q, T, P, yi, q_water, x->params_[c]); break;
    case K::DUAL_SITE_LANGMUIR_C: if (auto* x = dynamic_cast<DualSiteLangmuirIsotherm_Concentration*>(iso)) IF::DualSiteLangmuir_Concentration(q, T, P, yi, q_water, x->params_[c]); break;
    case K::EXTENDED_DUAL_SITE_LANGMUIR_C: if (auto* x = dynamic_cast<ExtendedDualSiteLangmuirIsotherm_Concentration*>(iso)) IF::DualSiteLangmuir_Concentration(q, T, P, yi, q_water, x->params_[c]); break;
    case K::ANTI_LANGMUIR: if (auto* x = dynamic_cast<AntiLangmuirIsotherm*>(iso)) IF::AntiLangmuir(q, T, P, yi, x->params_[c]); break;
    case K::FREUNDLICH: if (auto* x = dynamic_cast<FreundlichIsotherm*>(iso)) IF::Freundlich(q, T, P, yi, x->params_[c]); break;
    case K::SIPS: if (auto* x = dynamic_cast<SipsIsotherm*>(iso)) IF::Sips(q, T, P, yi, x->params_[c]); break;
    case K::SIPS_2: if (auto* x = dynamic_cast<Sips2Isotherm*>(iso)) IF::Sips2(q, T, P, yi, x->params_[c]); break;
    case K::EXTENDED_SIPS_2: if (auto* x = dynamic_cast<ExtendedSips2Isotherm*>(iso)) IF::ExtendedSips2(q, T, P, yi, x->params_[c]); break;
    case K::DUAL_SITE_SIPS: if (auto* x = dynamic_cast<DualSiteSipsIsotherm*>(iso)) IF::DualSiteSips(q, T, P, yi, x->params_[c]); break;
    case K::REDLICH_PETERSON: if (auto* x = dynamic_cast<RedlichPetersonIsotherm*>(iso)) IF::RedlichPeterson(q, T, P, yi, x->params_[c]); break;
    case K::TOTH: if (auto* x = dynamic_cast<TothIsotherm*>(iso)) IF::Toth(q, T, P, yi, x->params_[c]); break;
    case K::TOTH_2: if (auto* x = dynamic_cast<Toth2Isotherm*>(iso)) IF::Toth2(q, T, P, yi, x->params_[c]); break;
    case K::ARRHENIUS_SIPS: if (auto* x = dynamic_cast<ArrheniusSipsIsotherm*>(iso)) IF::ArrheniusSips(q, T, P, yi, x->params_[c]); break;
    case K::EXTENDED_LANGMUIR_FREUNDLICH: if (auto* x = dynamic_cast<ExtendedLangmuirFreundlichIsotherm*>(iso)) IF::ExtendedLangmuirFreundlich(q, T, P, yi, x->params_[c]); break;
    case K::GAB: if (auto* x = dynamic_cast<GabIsotherm*>(iso)) IF::Gab(q, T, P, yi, SaturationPressure(c, T), x->params_[c]); break;
    case K::GAB_YOUNG: if (auto* x = dynamic_cast<GabYoungIsotherm*>(iso)) IF::GabYoung(q, T, P, yi, SaturationPressure(c, T), x->params_[c]); break;
    case K::TOTH_WADST: if (auto* x = dynamic_cast<TothWadstIsotherm*>(iso)) IF::TothWadst(q, T, P, yi, q_water, x->params_[c]); break;
    case K::TOTH_HUMID: if (auto* x = dynamic_cast<TothHumidIsotherm*>(iso)) IF::TothHumid(q, T, P, yi, q_water, x->params_[c]); break;
    case K::CUBIC_LANGMUIR: if (auto* x = dynamic_cast<CubicLangmuirIsotherm*>(iso)) IF::CubicLangmuir(q, T, P, yi, x->params_[c]); break;
  }
  return q;
}

double IsostericHeat(IIsotherm* iso, const std::string& c, double T, double P, double yi, double q_water) {
  // q_st = R T^2 (d ln p / dT)_q = -R T^2 (dq/dT)_p / (dq/d ln p)_T, both
  // partials by central differences of EvalIsotherm (relative steps 1e-4).
  const double p = P * yi;
  if (!iso || !(p > 0.0) || !(T > 0.0)) return std::numeric_limits<double>::quiet_NaN();
  const double hT = 1e-4 * T, hl = 1e-4;
  const double qw = q_water;
  const double dq_dT = (EvalIsotherm(iso, c, T + hT, p, 1.0, qw) - EvalIsotherm(iso, c, T - hT, p, 1.0, qw)) / (2.0 * hT);
  const double dq_dlnp = (EvalIsotherm(iso, c, T, p * std::exp(hl), 1.0, qw) -
                          EvalIsotherm(iso, c, T, p * std::exp(-hl), 1.0, qw)) / (2.0 * hl);
  return -constants::kIdealGasConstant * T * T * dq_dT / dq_dlnp;
}

}  // namespace adsim
