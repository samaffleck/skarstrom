/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "adsim/equations/isotherm_kinds.h"

namespace adsim {

const char* IsothermKindString(IsothermType t) {
  using T = IsothermType;
  switch (t) {
    case T::INERT: return "inert";
    case T::HENRY: return "henry";
    case T::SINGLE_SITE_LANGMUIR_P: return "single_site_langmuir_p";
    case T::DUAL_SITE_LANGMUIR_P: return "dual_site_langmuir_p";
    case T::SINGLE_SITE_LANGMUIR_C: return "single_site_langmuir_c";
    case T::DUAL_SITE_LANGMUIR_C: return "dual_site_langmuir_c";
    case T::EXTENDED_SINGLE_SITE_LANGMUIR_P: return "extended_single_site_langmuir_p";
    case T::EXTENDED_DUAL_SITE_LANGMUIR_P: return "extended_dual_site_langmuir_p";
    case T::EXTENDED_SINGLE_SITE_LANGMUIR_C: return "extended_single_site_langmuir_c";
    case T::EXTENDED_DUAL_SITE_LANGMUIR_C: return "extended_dual_site_langmuir_c";
    case T::ANTI_LANGMUIR: return "anti_langmuir";
    case T::FREUNDLICH: return "freundlich";
    case T::SIPS: return "sips";
    case T::SIPS_2: return "sips_2";
    case T::EXTENDED_SIPS_2: return "extended_sips_2";
    case T::DUAL_SITE_SIPS: return "dual_site_sips";
    case T::REDLICH_PETERSON: return "redlich_peterson";
    case T::TOTH: return "toth";
    case T::TOTH_2: return "toth_2";
    case T::ARRHENIUS_SIPS: return "arrhenius_sips";
  }
  return "unknown";
}

std::vector<std::pair<std::string, ExpressionParam*>>
IsothermFields(IIsotherm* iso, const std::string& component) {
  using T = IsothermType;
  std::vector<std::pair<std::string, ExpressionParam*>> out;
  if (!iso) return out;
  auto add = [&](const char* n, ExpressionParam& p) { out.emplace_back(n, &p); };

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
        add("m2", it->second.m2); add("b2", it->second.b2); add("q2", it->second.q2); } }
      break;
    }
    case T::DUAL_SITE_LANGMUIR_C: {
      auto* c = dynamic_cast<DualSiteLangmuirIsotherm_Concentration*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("m1", it->second.m1); add("b1", it->second.b1); add("q1", it->second.q1);
        add("m2", it->second.m2); add("b2", it->second.b2); add("q2", it->second.q2); } }
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
        add("m2", it->second.m2); add("b2", it->second.b2); add("q2", it->second.q2); } }
      break;
    }
    case T::EXTENDED_DUAL_SITE_LANGMUIR_C: {
      auto* c = dynamic_cast<ExtendedDualSiteLangmuirIsotherm_Concentration*>(iso);
      if (c) { auto it = c->params_.find(component); if (it != c->params_.end()) {
        add("m1", it->second.m1); add("b1", it->second.b1); add("q1", it->second.q1);
        add("m2", it->second.m2); add("b2", it->second.b2); add("q2", it->second.q2); } }
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
      return true;
  }
  return false;  // out-of-range value; every enumerator is listed above (-Wswitch)
}

bool IsCompetitiveIsothermKind(IsothermType t) {
  using T = IsothermType;
  switch (t) {
    case T::EXTENDED_SINGLE_SITE_LANGMUIR_P:
    case T::EXTENDED_SINGLE_SITE_LANGMUIR_C:
    case T::EXTENDED_DUAL_SITE_LANGMUIR_P:
    case T::EXTENDED_DUAL_SITE_LANGMUIR_C:
    case T::EXTENDED_SIPS_2:
    case T::ARRHENIUS_SIPS:
      return true;
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
    case T::SIPS_2:
    case T::DUAL_SITE_SIPS:
    case T::TOTH_2:
      return false;
  }
  return false;  // out-of-range value; every enumerator is listed above (-Wswitch)
}

bool PairedIsothermKind(IsothermType t, IsothermType& out) {
  using T = IsothermType;
  static constexpr std::pair<T, T> kPairs[] = {
      {T::SINGLE_SITE_LANGMUIR_P, T::EXTENDED_SINGLE_SITE_LANGMUIR_P},
      {T::SINGLE_SITE_LANGMUIR_C, T::EXTENDED_SINGLE_SITE_LANGMUIR_C},
      {T::DUAL_SITE_LANGMUIR_P, T::EXTENDED_DUAL_SITE_LANGMUIR_P},
      {T::DUAL_SITE_LANGMUIR_C, T::EXTENDED_DUAL_SITE_LANGMUIR_C},
      {T::SIPS_2, T::EXTENDED_SIPS_2},
  };
  for (const auto& [pure, competitive] : kPairs) {
    if (t == pure) { out = competitive; return true; }
    if (t == competitive) { out = pure; return true; }
  }
  return false;
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
  }
}

double EvalIsotherm(IIsotherm* iso, const std::string& c, double T, double P, double yi) {
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
    case K::DUAL_SITE_LANGMUIR_P: if (auto* x = dynamic_cast<DualSiteLangmuirIsotherm*>(iso)) IF::DualSiteLangmuir(q, T, P, yi, x->params_[c]); break;
    case K::EXTENDED_DUAL_SITE_LANGMUIR_P: if (auto* x = dynamic_cast<ExtendedDualSiteLangmuirIsotherm*>(iso)) IF::DualSiteLangmuir(q, T, P, yi, x->params_[c]); break;
    case K::DUAL_SITE_LANGMUIR_C: if (auto* x = dynamic_cast<DualSiteLangmuirIsotherm_Concentration*>(iso)) IF::DualSiteLangmuir_Concentration(q, T, P, yi, x->params_[c]); break;
    case K::EXTENDED_DUAL_SITE_LANGMUIR_C: if (auto* x = dynamic_cast<ExtendedDualSiteLangmuirIsotherm_Concentration*>(iso)) IF::DualSiteLangmuir_Concentration(q, T, P, yi, x->params_[c]); break;
    case K::ANTI_LANGMUIR: if (auto* x = dynamic_cast<AntiLangmuirIsotherm*>(iso)) IF::AntiLangmuir(q, T, P, yi, x->params_[c]); break;
    case K::FREUNDLICH: if (auto* x = dynamic_cast<FreundlichIsotherm*>(iso)) IF::Freundlich(q, T, P, yi, x->params_[c]); break;
    case K::SIPS: if (auto* x = dynamic_cast<SipsIsotherm*>(iso)) IF::Sips(q, T, P, yi, x->params_[c]); break;
    case K::SIPS_2: if (auto* x = dynamic_cast<Sips2Isotherm*>(iso)) IF::Sips2(q, T, P, yi, x->params_[c]); break;
    case K::EXTENDED_SIPS_2: if (auto* x = dynamic_cast<ExtendedSips2Isotherm*>(iso)) IF::Sips2(q, T, P, yi, x->params_[c]); break;
    case K::DUAL_SITE_SIPS: if (auto* x = dynamic_cast<DualSiteSipsIsotherm*>(iso)) IF::DualSiteSips(q, T, P, yi, x->params_[c]); break;
    case K::REDLICH_PETERSON: if (auto* x = dynamic_cast<RedlichPetersonIsotherm*>(iso)) IF::RedlichPeterson(q, T, P, yi, x->params_[c]); break;
    case K::TOTH: if (auto* x = dynamic_cast<TothIsotherm*>(iso)) IF::Toth(q, T, P, yi, x->params_[c]); break;
    case K::TOTH_2: if (auto* x = dynamic_cast<Toth2Isotherm*>(iso)) IF::Toth2(q, T, P, yi, x->params_[c]); break;
    case K::ARRHENIUS_SIPS: if (auto* x = dynamic_cast<ArrheniusSipsIsotherm*>(iso)) IF::ArrheniusSips(q, T, P, yi, x->params_[c]); break;
  }
  return q;
}

}  // namespace adsim
