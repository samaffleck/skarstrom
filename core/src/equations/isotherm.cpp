/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/equations/isotherm.h"
#include "adsim/equations/equation.h"
#include "adsim/objects/porous_media.h"
#include "adsim/math/constants.h"

#include <cmath>


namespace adsim {

  inline double softplus(double x, double tau) {
    // numerically stable softplus: smooth approximation to max(0, x)
    if (x > 40.0*tau) return x;
    if (x < -40.0*tau) return 0.0;
    return tau * std::log1p(std::exp(x / tau));
  }

  inline double safe_pow(double base, double exponent) {
    return std::exp(exponent * std::log(std::max(base, 1e-300)));
  }

  constexpr double kSoftplusTau = 1e-6;

  namespace {
    inline double ExpOverRT(double T, double q) {
      return std::exp(q / (constants::kIdealGasConstant * T));
    }

    // The loading expressions with their temperature-only factors passed in,
    // shared by the per-cell functions below and the Update loops (which
    // compute each factor once per cell and component), so both paths do the
    // same arithmetic in the same order.
    inline double HenryQ(double e, double P, double yi, const HenryIsothermParameters& ips) {
      return ips.k * P * yi * e;
    }
    inline double LangmuirQ(double e1, double P, double yi, const LangmuirIsothermParameters& ips) {
      return ips.m1 * ips.b1 * e1 * P * yi / (1 + ips.b1 * e1 * P * yi);
    }
    inline double DualSiteLangmuirQ(double e1, double e2, double P, double yi, const DualSiteLangmuirIsothermParameters& ips) {
      return (ips.m1 * ips.b1 * e1 * P * yi / (1 + ips.b1 * e1 * P * yi)) +
        (ips.m2 * ips.b2 * e2 * P * yi / (1 + ips.b2 * e2 * P * yi));
    }
    inline double LangmuirConcentrationQ(double e1, double ci, const LangmuirIsothermParametersConcentration& ips) {
      return ips.m1 * ips.b1 * e1 * ci / (1 + ips.b1 * e1 * ci);
    }
    inline double DualSiteLangmuirConcentrationQ(double e1, double e2, double ci, const DualSiteLangmuirIsothermParametersConcentration& ips) {
      return (ips.m1 * ips.b1 * e1 * ci / (1 + ips.b1 * e1 * ci)) +
        (ips.m2 * ips.b2 * e2 * ci / (1 + ips.b2 * e2 * ci));
    }
    // eB = exp(B2 / (R T)), eA = exp(A2 ((T / T_ref) - 1)).
    inline double Sips2Q(double eB, double eA, double T, double P, double yi, const Sips2IsothermParameters& ips) {
      double yi_pos = softplus(yi, kSoftplusTau);
      double b_i = ips.B1 * eB;
      double s_i = ips.C1 + ips.C2 * ((T / ips.T_ref) - 1);
      double z = safe_pow(b_i * P * yi_pos, s_i);
      double q_max = ips.A1 * eA;
      return q_max * z / (1 + z);
    }

    // Competitive Sips term z = (bP y)^s of extended_sips_2 and arrhenius_sips
    // near zero mole fraction (docs/isotherms.md, "Competitive Sips kinds near
    // zero mole fraction"). For s < 1 the power has an infinite slope at y = 0,
    // which the Newton iteration cannot linearise. Below kSipsTraceY0 it is
    // replaced by the quadratic in y that is zero at y = 0 and matches the
    // power's value and slope at y0 (C1):
    //     z = (bP y0)^s t ((2 - s) - (1 - s) t),   t = y / y0.
    // Exact for y >= y0; concave and increasing on [0, y0] for s < 1 (no
    // inflection, so the isotherm stays favourable), with the finite Henry slope
    // (2 - s)(bP y0)^s / y0 at y = 0. Below 0 (integrator excursions ahead of a
    // front) the slope is continued by (2 - s)(bP y0)^s tanh(t), which stays
    // bounded. y0 is 1000x the mole-fraction absolute tolerance (1e-8), so a
    // Newton correction near y = 0 sees an almost linear term. For s >= 2 the
    // quadratic would not be monotone, so the power is used, clipped at 0
    // (C1 there since its slope at 0 is 0).
    //
    // This replaces the former smooth step y^2 / (y^2 + 1e-6), which halved the
    // loading at y = 1e-3 and cut it by 1 % at y = 1e-2 (WP19-1).
    constexpr double kSipsTraceY0 = 1e-5;
    inline double SipsTraceTerm(double bP, double y, double s) {
      if (y >= kSipsTraceY0) return safe_pow(bP * y, s);
      if (s >= 2.0) return y > 0.0 ? safe_pow(bP * y, s) : 0.0;
      const double z0 = safe_pow(bP * kSipsTraceY0, s);
      const double t = y / kSipsTraceY0;
      if (y >= 0.0) return z0 * t * ((2.0 - s) - (1.0 - s) * t);
      return z0 * (2.0 - s) * std::tanh(t);
    }
    // eX = exp(x (1 - T / T_0)), eH = exp(dH0 / (R T_0) (T_0 / T - 1)).
    inline double Toth2Q(double eX, double eH, double T, double P, double yi, const Toth2IsothermParameters& ips) {
      double yi_pos = softplus(yi, kSoftplusTau);
      double qs = ips.q_sat * eX;
      double b = ips.b0 * eH;
      double w = ips.w0 + ips.alpha * (1 - ips.T_0 / T);
      double bp = b * P * yi_pos;
      double bpv = safe_pow(bp, w);
      double denom = safe_pow(1 + bpv, (1 / w));
      return qs * bp / denom;
    }
  } // namespace

	void IsothermFunction::Henry(double& qi, double T, double P, double yi, const HenryIsothermParameters& ips) {
		qi = HenryQ(std::exp(ips.Q / (constants::kIdealGasConstant * T)), P, yi, ips);
	}

	void IsothermFunction::Langmuir(double& qi, double T, double P, double yi, const LangmuirIsothermParameters& ips) {
		const double e1 = std::exp(ips.q1 / (constants::kIdealGasConstant * T));
		qi = LangmuirQ(e1, P, yi, ips);
	}
	
	void IsothermFunction::AntiLangmuir(double& qi, double T, double P, double yi, const AntiLangmuirIsothermParameters& ips) {
		qi = ips.a * P * yi / (1 - ips.b * P * yi);
	}
	
	void IsothermFunction::Freundlich(double& qi, double T, double P, double yi, const FreundlichIsothermParameters& ips) {
		double yi_pos = softplus(yi, kSoftplusTau);
		qi = safe_pow(ips.a * P * yi_pos, (1 / ips.v));
	}

	void IsothermFunction::Sips(double& qi, double T, double P, double yi, const SipsIsothermParameters& ips) {
		double yi_pos = softplus(yi, kSoftplusTau);
		double z = safe_pow(ips.b * P * yi_pos, (1 / ips.v));
		qi = ips.q_sat * z / (1 + z);
	}

	void IsothermFunction::Sips2(double& qi, double T, double P, double yi, const Sips2IsothermParameters& ips) {
		const double& R = constants::kIdealGasConstant;
    const double eB = std::exp(ips.B2 / (R * T));
    const double eA = std::exp(ips.A2 * ((T / ips.T_ref) - 1));
    qi = Sips2Q(eB, eA, T, P, yi, ips);
	}

	void IsothermFunction::ExtendedSips2(double& qi, double T, double P, double yi, const Sips2IsothermParameters& ips) {
		// Same arithmetic as ExtendedSips2Isotherm::UpdateAll with one component.
		const double& R = constants::kIdealGasConstant;
		const double b_i = ips.B1 * std::exp(ips.B2 / (R * T));
		const double s_i = ips.C1 + ips.C2 * ((T / ips.T_ref) - 1);
		const double z = SipsTraceTerm(b_i * P, yi, s_i);
		const double q_max = ips.A1 * std::exp(ips.A2 * ((T / ips.T_ref) - 1));
		qi = q_max * z / (1 + z);
	}

	void IsothermFunction::ArrheniusSips(double& qi, double T, double P, double yi, const ArrheniusSipsIsothermParameters& ips) {
		const double R = constants::kIdealGasConstant;
		// Extract to plain doubles first so all arithmetic is numeric (not built
		// up as symbolic ExpressionParam operations).
		const double a = ips.n_inf_a, b = ips.n_inf_b;
		const double A = ips.k_A, B = ips.k_B;
		const double c_ref = ips.c_ref, c_amp = ips.c_amp, c_beta = ips.c_beta, c_Tref = ips.c_Tref;
		const double n_inf = a * std::exp(b / (R * T));
		const double k_i = A * std::exp(B / (R * T));
		const double c_i = c_ref + c_amp * (2.0 / M_PI) * std::atan(c_beta * (T - c_Tref));
		// The competitive term as ArrheniusSipsIsotherm::UpdateAll evaluates it.
		const double z = SipsTraceTerm(k_i * P, yi, c_i);
		qi = n_inf * z / (1 + z);
	}

	void IsothermFunction::DualSiteSips(double& qi, double T, double P, double yi, const DualSiteSipsIsothermParameters& ips) {
		double yi_pos = softplus(yi, kSoftplusTau);
		double z1 = safe_pow(ips.b1 * P * yi_pos, (1 / ips.v1));
		double z2 = safe_pow(ips.b2 * P * yi_pos, (1 / ips.v2));
		qi = ips.q_sat1 * z1 / (1 + z1) + ips.q_sat2 * z2 / (1 + z2);
	}
	
	void IsothermFunction::RedlichPeterson(double& qi, double T, double P, double yi, const RedlichPetersonIsothermParameters& ips) {
		double yi_pos = softplus(yi, kSoftplusTau);
		double z = safe_pow(ips.b * P * yi_pos, ips.v);
		qi = ips.a * P * yi_pos / (1 + z);
	}
	
	void IsothermFunction::Toth(double& qi, double T, double P, double yi, const TothIsothermParameters& ips) {
    double yi_pos = softplus(yi, kSoftplusTau);
    double bp = ips.b * P * yi_pos;
		double bpv = safe_pow(bp, ips.v);
		double denom = safe_pow(1 + bpv, (1 / ips.v));
		qi = ips.q_sat * bp / denom;
	}

	void IsothermFunction::Toth2(double& qi, double T, double P, double yi, const Toth2IsothermParameters& ips) {
		const double& R = constants::kIdealGasConstant;
    const double eX = std::exp(ips.x * (1 - T / ips.T_0));
    const double eH = std::exp(ips.dH0 / (R * ips.T_0) * (ips.T_0 / T - 1));
		qi = Toth2Q(eX, eH, T, P, yi, ips);
	}

	void IsothermFunction::DualSiteLangmuir(double& qi, double T, double P, double yi, const DualSiteLangmuirIsothermParameters& ips) {
		const double e1 = std::exp(ips.q1 / (constants::kIdealGasConstant * T));
		const double e2 = std::exp(ips.q2 / (constants::kIdealGasConstant * T));
		qi = DualSiteLangmuirQ(e1, e2, P, yi, ips);
	}

	void IsothermFunction::Langmuir_Concentration(double& qi, double T, double P, double yi, const LangmuirIsothermParametersConcentration& ips) {
		double ci = P * yi / (constants::kIdealGasConstant * T);
		const double e1 = std::exp(ips.q1 / (constants::kIdealGasConstant * T));
		qi = LangmuirConcentrationQ(e1, ci, ips);
	}

	void IsothermFunction::DualSiteLangmuir_Concentration(double& qi, double T, double P, double yi, const DualSiteLangmuirIsothermParametersConcentration& ips) {
		double ci = P * yi / (constants::kIdealGasConstant * T);
		const double e1 = std::exp(ips.q1 / (constants::kIdealGasConstant * T));
		const double e2 = std::exp(ips.q2 / (constants::kIdealGasConstant * T));
		qi = DualSiteLangmuirConcentrationQ(e1, e2, ci, ips);
	}

  void InertIsotherm::Update(IsothermUpdateParams& state) {
    auto& component = state.component;
    auto& qi_sat = state.qi_sat_map.at(component).var;
    int start_index = state.start_index;
    int end_index = state.end_index; 
    for (int n = start_index; n <= end_index; ++n) {
      qi_sat[n] = 0.0;
    }
	}

	void HenryIsotherm::Update(IsothermUpdateParams& state) {
    auto& component = state.component;
    auto& qi_sat = state.qi_sat_map.at(component).var;
    auto& yi = state.yi_map.at(component).var;
    auto& T = state.T;
    auto& P = state.P;
    auto& params = params_[component];
    int start_index = state.start_index;
    int end_index = state.end_index;

    const double Q = params.Q;
    for (int n = start_index; n <= end_index; ++n) {
      qi_sat[n] = HenryQ(ExpOverRT(T[n], Q), P[n], yi[n], params);
    }
	}

	void LangmuirIsotherm::Update(IsothermUpdateParams& state) {
    auto& component = state.component;
    auto& qi_sat = state.qi_sat_map.at(component).var;
    auto& yi = state.yi_map.at(component).var;
    auto& T = state.T;
    auto& P = state.P;
    auto& params = params_[component];
    int start_index = state.start_index;
    int end_index = state.end_index;

    const double q1 = params.q1;
    for (int n = start_index; n <= end_index; ++n) {
      qi_sat[n] = LangmuirQ(ExpOverRT(T[n], q1), P[n], yi[n], params);
    }
	}

  void AntiLangmuirIsotherm::Update(IsothermUpdateParams& state) {
    auto& component = state.component;
    auto& qi_sat = state.qi_sat_map.at(component).var;
    auto& yi = state.yi_map.at(component).var;
    auto& T = state.T;
    auto& P = state.P;
    auto& params = params_[component];
    int start_index = state.start_index;
    int end_index = state.end_index;

    for (int n = start_index; n <= end_index; ++n) {
      IsothermFunction::AntiLangmuir(qi_sat[n], T[n], P[n], yi[n], params);
    }
  }

  void FreundlichIsotherm::Update(IsothermUpdateParams& state) {
    auto& component = state.component;
    auto& qi_sat = state.qi_sat_map.at(component).var;
    auto& yi = state.yi_map.at(component).var;
    auto& T = state.T;
    auto& P = state.P;
    auto& params = params_[component];
    int start_index = state.start_index;
    int end_index = state.end_index;

    for (int n = start_index; n <= end_index; ++n) {
      IsothermFunction::Freundlich(qi_sat[n], T[n], P[n], yi[n], params);
    }
	}

  void SipsIsotherm::Update(IsothermUpdateParams& state) {
    auto& component = state.component;
    auto& qi_sat = state.qi_sat_map.at(component).var;
    auto& yi = state.yi_map.at(component).var;
    auto& T = state.T;
    auto& P = state.P;
    auto& params = params_[component];
    int start_index = state.start_index;
    int end_index = state.end_index;

    for (int n = start_index; n <= end_index; ++n) {
      IsothermFunction::Sips(qi_sat[n], T[n], P[n], yi[n], params);
    }
	}
	
  void Sips2Isotherm::Update(IsothermUpdateParams& state) {
    auto& component = state.component;
    auto& qi_sat = state.qi_sat_map.at(component).var;
    auto& yi = state.yi_map.at(component).var;
    auto& T = state.T;
    auto& P = state.P;
    auto& params = params_[component];
    int start_index = state.start_index;
    int end_index = state.end_index;

    const double R = constants::kIdealGasConstant;
    const double B2 = params.B2, A2 = params.A2, T_ref = params.T_ref;
    for (int n = start_index; n <= end_index; ++n) {
      const double eB = std::exp(B2 / (R * T[n]));
      const double eA = std::exp(A2 * ((T[n] / T_ref) - 1));
      qi_sat[n] = Sips2Q(eB, eA, T[n], P[n], yi[n], params);
    }
	}

  void DualSiteSipsIsotherm::Update(IsothermUpdateParams& state) {
    auto& component = state.component;
    auto& qi_sat = state.qi_sat_map.at(component).var;
    auto& yi = state.yi_map.at(component).var;
    auto& T = state.T;
    auto& P = state.P;
    auto& params = params_[component];
    int start_index = state.start_index;
    int end_index = state.end_index;

    for (int n = start_index; n <= end_index; ++n) {
      IsothermFunction::DualSiteSips(qi_sat[n], T[n], P[n], yi[n], params);
    }
  }

  void RedlichPetersonIsotherm::Update(IsothermUpdateParams& state) {
    auto& component = state.component;
    auto& qi_sat = state.qi_sat_map.at(component).var;
    auto& yi = state.yi_map.at(component).var;
    auto& T = state.T;
    auto& P = state.P;
    auto& params = params_[component];
    int start_index = state.start_index;
    int end_index = state.end_index;

    for (int n = start_index; n <= end_index; ++n) {
      IsothermFunction::RedlichPeterson(qi_sat[n], T[n], P[n], yi[n], params);
    }
	}

  void TothIsotherm::Update(IsothermUpdateParams& state) {
    auto& component = state.component;
    auto& qi_sat = state.qi_sat_map.at(component).var;
    auto& yi = state.yi_map.at(component).var;
    auto& T = state.T;
    auto& P = state.P;
    auto& params = params_[component];
    int start_index = state.start_index;
    int end_index = state.end_index;

    for (int n = start_index; n <= end_index; ++n) {
      IsothermFunction::Toth(qi_sat[n], T[n], P[n], yi[n], params);
    }
	}

  void Toth2Isotherm::Update(IsothermUpdateParams& state) {
    auto& component = state.component;
    auto& qi_sat = state.qi_sat_map.at(component).var;
    auto& yi = state.yi_map.at(component).var;
    auto& T = state.T;
    auto& P = state.P;
    auto& params = params_[component];
    int start_index = state.start_index;
    int end_index = state.end_index;

    const double R = constants::kIdealGasConstant;
    const double x = params.x, T_0 = params.T_0, dH0 = params.dH0;
    for (int n = start_index; n <= end_index; ++n) {
      const double eX = std::exp(x * (1 - T[n] / T_0));
      const double eH = std::exp(dH0 / (R * T_0) * (T_0 / T[n] - 1));
      qi_sat[n] = Toth2Q(eX, eH, T[n], P[n], yi[n], params);
    }
	}

  void DualSiteLangmuirIsotherm::Update(IsothermUpdateParams& state) {
    auto& component = state.component;
    auto& qi_sat = state.qi_sat_map.at(component).var;
    auto& yi = state.yi_map.at(component).var;
    auto& T = state.T;
    auto& P = state.P;
    auto& params = params_[component];
    int start_index = state.start_index;
    int end_index = state.end_index;

    const double q1 = params.q1, q2 = params.q2;
    for (int n = start_index; n <= end_index; ++n) {
      const double e1 = ExpOverRT(T[n], q1);
      const double e2 = ExpOverRT(T[n], q2);
      qi_sat[n] = DualSiteLangmuirQ(e1, e2, P[n], yi[n], params);
    }
	}

  void LangmuirIsotherm_Concentration::Update(IsothermUpdateParams& state) {
    auto& component = state.component;
    auto& qi_sat = state.qi_sat_map.at(component).var;
    auto& yi = state.yi_map.at(component).var;
    auto& T = state.T;
    auto& P = state.P;
    auto& params = params_[component];
    int start_index = state.start_index;
    int end_index = state.end_index;

    const double q1 = params.q1;
    for (int n = start_index; n <= end_index; ++n) {
      double ci = P[n] * yi[n] / (constants::kIdealGasConstant * T[n]);
      qi_sat[n] = LangmuirConcentrationQ(ExpOverRT(T[n], q1), ci, params);
    }
	}

  void DualSiteLangmuirIsotherm_Concentration::Update(IsothermUpdateParams& state) {
    auto& component = state.component;
    auto& qi_sat = state.qi_sat_map.at(component).var;
    auto& yi = state.yi_map.at(component).var;
    auto& T = state.T;
    auto& P = state.P;
    auto& params = params_[component];
    int start_index = state.start_index;
    int end_index = state.end_index;

    const double q1 = params.q1, q2 = params.q2;
    for (int n = start_index; n <= end_index; ++n) {
      double ci = P[n] * yi[n] / (constants::kIdealGasConstant * T[n]);
      const double e1 = ExpOverRT(T[n], q1);
      const double e2 = ExpOverRT(T[n], q2);
      qi_sat[n] = DualSiteLangmuirConcentrationQ(e1, e2, ci, params);
    }
	}

  // --- Competitive isotherms --------------------------------------------------
  // Each class's UpdateAll evaluates every requested component of the adsorbent
  // in one pass so the per-cell denominators (sums over all components, in
  // params_ map order) are shared; the per-component Update is the same pass
  // for a single component, so both give bit-identical loadings.

  void ExtendedLangmuirIsotherm::Update(IsothermUpdateParams& state) {
    UpdateAll(state, std::vector<std::string>{ state.component });
  }

  void ExtendedLangmuirIsotherm::UpdateAll(IsothermUpdateParams& state, const std::vector<std::string>& components) {
    auto& T = state.T;
    auto& P = state.P;
    int start_index = state.start_index;
    int end_index = state.end_index;

    // Materialise params for every requested component first (matches the
    // operator[] insertion behaviour of the former per-component Update).
    for (const auto& comp : components) params_[comp];

    // Hoist map lookups out of the cell loop. The denominator terms iterate
    // params_ in map order: summation order is observable in the result.
    struct Term { const LangmuirIsothermParameters* p; const std::vector<double>* yi; double q1; double e = 0.0; };
    std::vector<Term> terms;
    terms.reserve(params_.size());
    for (const auto& [comp, comp_params] : params_) {
      terms.push_back({ &comp_params, &state.yi_map.at(comp).var, comp_params.q1 });
    }

    struct Output { const LangmuirIsothermParameters* p; const std::vector<double>* yi; std::vector<double>* qi_sat; size_t term_index; };
    std::vector<Output> outputs;
    outputs.reserve(components.size());
    for (const auto& comp : components) {
      const auto* p = &params_.at(comp);
      size_t idx = 0;
      while (idx < terms.size() && terms[idx].p != p) ++idx;
      outputs.push_back({ p, &state.yi_map.at(comp).var, &state.qi_sat_map.at(comp).var, idx });
    }

    for (int n = start_index; n <= end_index; ++n) {
      // Denominator over all components, keeping each exponential so every
      // output component reuses it with the original operation order.
      double denominator = 0.0;
      for (auto& t : terms) {
        t.e = ExpOverRT(T[n], t.q1);
        denominator += t.p->b1 * t.e * P[n] * (*t.yi)[n];
      }

      for (const auto& out : outputs) {
        const double e = terms[out.term_index].e;
        (*out.qi_sat)[n] = out.p->m1 * out.p->b1 * e * P[n] * (*out.yi)[n] / (1 + denominator);
      }
    }
  }

	void ExtendedDualSiteLangmuirIsotherm::Update(IsothermUpdateParams& state) {
    UpdateAll(state, std::vector<std::string>{ state.component });
	}

  void ExtendedDualSiteLangmuirIsotherm::UpdateAll(IsothermUpdateParams& state, const std::vector<std::string>& components) {
    auto& T = state.T;
    auto& P = state.P;
    int start_index = state.start_index;
    int end_index = state.end_index;

    for (const auto& comp : components) params_[comp];

    struct Term { const DualSiteLangmuirIsothermParameters* p; const std::vector<double>* yi; double q1; double q2; double e1 = 0.0; double e2 = 0.0; };
    std::vector<Term> terms;
    terms.reserve(params_.size());
    for (const auto& [comp, comp_params] : params_) {
      terms.push_back({ &comp_params, &state.yi_map.at(comp).var, comp_params.q1, comp_params.q2 });
    }

    struct Output { const DualSiteLangmuirIsothermParameters* p; const std::vector<double>* yi; std::vector<double>* qi_sat; size_t term_index; };
    std::vector<Output> outputs;
    outputs.reserve(components.size());
    for (const auto& comp : components) {
      const auto* p = &params_.at(comp);
      size_t idx = 0;
      while (idx < terms.size() && terms[idx].p != p) ++idx;
      outputs.push_back({ p, &state.yi_map.at(comp).var, &state.qi_sat_map.at(comp).var, idx });
    }

    for (int n = start_index; n <= end_index; ++n) {
      double denominator1 = 0.0;
      double denominator2 = 0.0;
      for (auto& t : terms) {
        t.e1 = ExpOverRT(T[n], t.q1);
        t.e2 = ExpOverRT(T[n], t.q2);
        denominator1 += t.p->b1 * t.e1 * P[n] * (*t.yi)[n];
        denominator2 += t.p->b2 * t.e2 * P[n] * (*t.yi)[n];
      }

      for (const auto& out : outputs) {
        const auto& t = terms[out.term_index];
        (*out.qi_sat)[n] =
          out.p->m1 * out.p->b1 * t.e1 * P[n] * (*out.yi)[n] / (1 + denominator1) +
          out.p->m2 * out.p->b2 * t.e2 * P[n] * (*out.yi)[n] / (1 + denominator2);
      }
    }
  }

	void ExtendedLangmuirIsotherm_Concentration::Update(IsothermUpdateParams& state) {
    UpdateAll(state, std::vector<std::string>{ state.component });
	}

  void ExtendedLangmuirIsotherm_Concentration::UpdateAll(IsothermUpdateParams& state, const std::vector<std::string>& components) {
    auto& T = state.T;
    auto& P = state.P;
    int start_index = state.start_index;
    int end_index = state.end_index;
    const double R = constants::kIdealGasConstant;

    for (const auto& comp : components) params_[comp];

    struct Term { const LangmuirIsothermParametersConcentration* p; const std::vector<double>* yi; double q1; double e = 0.0; };
    std::vector<Term> terms;
    terms.reserve(params_.size());
    for (const auto& [comp, comp_params] : params_) {
      terms.push_back({ &comp_params, &state.yi_map.at(comp).var, comp_params.q1 });
    }

    struct Output { const LangmuirIsothermParametersConcentration* p; const std::vector<double>* yi; std::vector<double>* qi_sat; size_t term_index; };
    std::vector<Output> outputs;
    outputs.reserve(components.size());
    for (const auto& comp : components) {
      const auto* p = &params_.at(comp);
      size_t idx = 0;
      while (idx < terms.size() && terms[idx].p != p) ++idx;
      outputs.push_back({ p, &state.yi_map.at(comp).var, &state.qi_sat_map.at(comp).var, idx });
    }

    for (int n = start_index; n <= end_index; ++n) {
      // Denominator over all components (using concentration)
      double denominator = 0.0;
      for (auto& t : terms) {
        double ci = P[n] * (*t.yi)[n] / (R * T[n]);
        t.e = ExpOverRT(T[n], t.q1);
        denominator += t.p->b1 * t.e * ci;
      }

      for (const auto& out : outputs) {
        const double e = terms[out.term_index].e;
        double ci = P[n] * (*out.yi)[n] / (R * T[n]);
        (*out.qi_sat)[n] = out.p->m1 * out.p->b1 * e * ci / (1 + denominator);
      }
    }
  }

  void ExtendedDualSiteLangmuirIsotherm_Concentration::Update(IsothermUpdateParams& state) {
    UpdateAll(state, std::vector<std::string>{ state.component });
  }

  void ExtendedDualSiteLangmuirIsotherm_Concentration::UpdateAll(IsothermUpdateParams& state, const std::vector<std::string>& components) {
    auto& T = state.T;
    auto& P = state.P;
    int start_index = state.start_index;
    int end_index = state.end_index;
    const double R = constants::kIdealGasConstant;

    for (const auto& comp : components) params_[comp];

    struct Term { const DualSiteLangmuirIsothermParametersConcentration* p; const std::vector<double>* yi; double q1; double q2; double e1 = 0.0; double e2 = 0.0; };
    std::vector<Term> terms;
    terms.reserve(params_.size());
    for (const auto& [comp, comp_params] : params_) {
      terms.push_back({ &comp_params, &state.yi_map.at(comp).var, comp_params.q1, comp_params.q2 });
    }

    struct Output { const DualSiteLangmuirIsothermParametersConcentration* p; const std::vector<double>* yi; std::vector<double>* qi_sat; size_t term_index; };
    std::vector<Output> outputs;
    outputs.reserve(components.size());
    for (const auto& comp : components) {
      const auto* p = &params_.at(comp);
      size_t idx = 0;
      while (idx < terms.size() && terms[idx].p != p) ++idx;
      outputs.push_back({ p, &state.yi_map.at(comp).var, &state.qi_sat_map.at(comp).var, idx });
    }

    for (int n = start_index; n <= end_index; ++n) {
      // Denominators over all components (using concentration)
      double denominator1 = 0.0;
      double denominator2 = 0.0;
      for (auto& t : terms) {
        double ci = P[n] * (*t.yi)[n] / (R * T[n]);
        t.e1 = ExpOverRT(T[n], t.q1);
        t.e2 = ExpOverRT(T[n], t.q2);
        denominator1 += t.p->b1 * t.e1 * ci;
        denominator2 += t.p->b2 * t.e2 * ci;
      }

      for (const auto& out : outputs) {
        const auto& t = terms[out.term_index];
        double ci = P[n] * (*out.yi)[n] / (R * T[n]);
        (*out.qi_sat)[n] =
          out.p->m1 * out.p->b1 * t.e1 * ci / (1 + denominator1) +
          out.p->m2 * out.p->b2 * t.e2 * ci / (1 + denominator2);
      }
    }
  }

	void ExtendedSips2Isotherm::Update(IsothermUpdateParams& state) {
    UpdateAll(state, std::vector<std::string>{ state.component });
	}

  void ExtendedSips2Isotherm::UpdateAll(IsothermUpdateParams& state, const std::vector<std::string>& components) {
    const double& R = constants::kIdealGasConstant;
    auto& T = state.T;
    auto& P = state.P;
    int start_index = state.start_index;
    int end_index = state.end_index;

    for (const auto& comp : components) params_[comp];

    // Competitive Sips term z = (b_i P y_i)^{s_i}, continued below y = 1e-5
    // (SipsTraceTerm). eB = exp(B2 / (R T)) for this cell.
    auto get_z = [&](double temp, double press, double mf, double eB, const Sips2IsothermParameters& params) -> double {
      const double b_i = params.B1 * eB;
      const double s_i = params.C1 + params.C2 * ((temp / params.T_ref) - 1);
      return SipsTraceTerm(b_i * press, mf, s_i);
    };

    // Hoist map lookups out of the cell loop. The denominator terms iterate
    // params_ in map order, so the shared sum is bit-identical whichever
    // components are requested; each output component reuses its own z.
    struct Term { const Sips2IsothermParameters* p; const std::vector<double>* yi; double B2; double z = 0.0; };
    std::vector<Term> terms;
    terms.reserve(params_.size());
    for (const auto& [comp, comp_params] : params_) {
      terms.push_back({ &comp_params, &state.yi_map.at(comp).var, comp_params.B2 });
    }

    struct Output { const Sips2IsothermParameters* p; std::vector<double>* qi_sat; double A2; double T_ref; size_t term_index; };
    std::vector<Output> outputs;
    outputs.reserve(components.size());
    for (const auto& comp : components) {
      const auto* p = &params_.at(comp);
      size_t idx = 0;
      while (idx < terms.size() && terms[idx].p != p) ++idx;
      outputs.push_back({ p, &state.qi_sat_map.at(comp).var, p->A2, p->T_ref, idx });
    }

    for (int n = start_index; n <= end_index; ++n) {
      // Calculate denominator over all components
      double denominator = 0.0;
      for (auto& t : terms) {
        const double eB = std::exp(t.B2 / (R * T[n]));
        t.z = get_z(T[n], P[n], (*t.yi)[n], eB, *t.p);
        denominator += t.z;
      }

      for (const auto& out : outputs) {
        const double eA = std::exp(out.A2 * ((T[n] / out.T_ref) - 1));
        double q_max = out.p->A1 * eA;
        double z = terms[out.term_index].z;
        (*out.qi_sat)[n] = q_max * z / (1 + denominator);
      }
    }
  }

  // --- Arrhenius Sips (Schell/Casas 2012) -----------------------------------
  // Shared helpers: evaluate the temperature-dependent Sips parameters and the
  // competitive term z_i = (k_i·P·y_i)^{c_i}, continued below y = 1e-5 exactly
  // as in ExtendedSips2Isotherm (SipsTraceTerm), so the Jacobian stays finite as
  // a mole fraction passes through 0.
  namespace {

    inline double ArrSips_ninf(double T, const ArrheniusSipsIsothermParameters& p) {
      const double R = constants::kIdealGasConstant;
      const double a = p.n_inf_a, b = p.n_inf_b;
      return a * std::exp(b / (R * T));
    }
    inline double ArrSips_k(double T, const ArrheniusSipsIsothermParameters& p) {
      const double R = constants::kIdealGasConstant;
      const double A = p.k_A, B = p.k_B;
      return A * std::exp(B / (R * T));
    }
    inline double ArrSips_c(double T, const ArrheniusSipsIsothermParameters& p) {
      // c(T) = c_ref + c_amp·(2/π)·arctan(c_beta·(T − c_Tref)); c_amp = 0 → constant.
      const double c_ref = p.c_ref, c_amp = p.c_amp, c_beta = p.c_beta, c_Tref = p.c_Tref;
      return c_ref + c_amp * (2.0 / M_PI) * std::atan(c_beta * (T - c_Tref));
    }
    inline double ArrSips_z(double T, double P, double mf,
                            const ArrheniusSipsIsothermParameters& p) {
      const double k_i = ArrSips_k(T, p);
      const double c_i = ArrSips_c(T, p);
      return SipsTraceTerm(k_i * P, mf, c_i);
    }
  } // namespace

  // --- Extended Langmuir-Freundlich (Yang 1987; toPSAil) ----------------------
  // z_i = b_i (P y_i / P_ref)^{n_i}. For n < 1 the power has an infinite slope at
  // y = 0; it is continued below y = 1e-5 exactly as the competitive Sips kinds
  // are (SipsTraceTerm with bP = P / P_ref, s = n), so the isotherm is exact for
  // y >= 1e-5 and has a finite Henry slope below. SipsTraceTerm also covers
  // n >= 2, so the temperature-dependent n needs no special case.
  namespace {
    inline double Elf_qs(double T, const ExtendedLangmuirFreundlichIsothermParameters& p) {
      const double k1 = p.k1, k2 = p.k2;
      return k1 + k2 * T;
    }
    inline double Elf_z(double T, double P, double mf, const ExtendedLangmuirFreundlichIsothermParameters& p) {
      const double k3 = p.k3, k4 = p.k4, k5 = p.k5, k6 = p.k6, P_ref = p.P_ref;
      const double b = k3 * std::exp(k4 / T);
      const double n = 1.0 / (k5 + k6 / T);
      return b * SipsTraceTerm(P / P_ref, mf, n);
    }
  } // namespace

  void IsothermFunction::ExtendedLangmuirFreundlich(double& qi, double T, double P, double yi,
                                                    const ExtendedLangmuirFreundlichIsothermParameters& ips) {
    const double z = Elf_z(T, P, yi, ips);
    qi = Elf_qs(T, ips) * z / (1.0 + z);
  }

  void ExtendedLangmuirFreundlichIsotherm::Update(IsothermUpdateParams& state) {
    UpdateAll(state, std::vector<std::string>{ state.component });
  }

  void ExtendedLangmuirFreundlichIsotherm::UpdateAll(IsothermUpdateParams& state,
                                                     const std::vector<std::string>& components) {
    auto& T = state.T;
    auto& P = state.P;
    const int start_index = state.start_index;
    const int end_index = state.end_index;

    for (const auto& comp : components) params_[comp];

    // The denominator iterates params_ in map order (bit-identical whichever
    // components are requested), as the other competitive kinds do.
    struct Term { const ExtendedLangmuirFreundlichIsothermParameters* p; const std::vector<double>* yi; double z = 0.0; };
    std::vector<Term> terms;
    terms.reserve(params_.size());
    for (const auto& [comp, comp_params] : params_)
      terms.push_back({ &comp_params, &state.yi_map.at(comp).var });

    struct Output { const ExtendedLangmuirFreundlichIsothermParameters* p; std::vector<double>* qi_sat; size_t term_index; };
    std::vector<Output> outputs;
    outputs.reserve(components.size());
    for (const auto& comp : components) {
      const auto* p = &params_.at(comp);
      size_t idx = 0;
      while (idx < terms.size() && terms[idx].p != p) ++idx;
      outputs.push_back({ p, &state.qi_sat_map.at(comp).var, idx });
    }

    for (int n = start_index; n <= end_index; ++n) {
      double denominator = 0.0;
      for (auto& t : terms) {
        t.z = Elf_z(T[n], P[n], (*t.yi)[n], *t.p);
        denominator += t.z;
      }
      for (const auto& out : outputs)
        (*out.qi_sat)[n] = Elf_qs(T[n], *out.p) * terms[out.term_index].z / (1.0 + denominator);
    }
  }

  void ArrheniusSipsIsotherm::Update(IsothermUpdateParams& state) {
    UpdateAll(state, std::vector<std::string>{ state.component });
  }

  void ArrheniusSipsIsotherm::UpdateAll(IsothermUpdateParams& state,
                                        const std::vector<std::string>& components) {
    auto& T = state.T;
    auto& P = state.P;
    int start_index = state.start_index;
    int end_index = state.end_index;

    for (const auto& comp : components) params_[comp];

    // Hoist map lookups; iterate params_ in map order for a shared denominator
    // that is bit-identical whichever components are requested.
    struct Term { const ArrheniusSipsIsothermParameters* p; const std::vector<double>* yi; double z = 0.0; };
    std::vector<Term> terms;
    terms.reserve(params_.size());
    for (const auto& [comp, comp_params] : params_) {
      terms.push_back({ &comp_params, &state.yi_map.at(comp).var });
    }

    struct Output { const ArrheniusSipsIsothermParameters* p; std::vector<double>* qi_sat; size_t term_index; };
    std::vector<Output> outputs;
    outputs.reserve(components.size());
    for (const auto& comp : components) {
      const auto* p = &params_.at(comp);
      size_t idx = 0;
      while (idx < terms.size() && terms[idx].p != p) ++idx;
      outputs.push_back({ p, &state.qi_sat_map.at(comp).var, idx });
    }

    for (int n = start_index; n <= end_index; ++n) {
      double denominator = 0.0;
      for (auto& t : terms) {
        t.z = ArrSips_z(T[n], P[n], (*t.yi)[n], *t.p);
        denominator += t.z;
      }

      for (const auto& out : outputs) {
        const double n_inf = ArrSips_ninf(T[n], *out.p);
        const double z = terms[out.term_index].z;
        (*out.qi_sat)[n] = n_inf * z / (1 + denominator);
      }
    }
  }

} // end adsim namespace


// Register the derived type with the base type
CEREAL_REGISTER_TYPE(adsim::HenryIsotherm)
CEREAL_REGISTER_TYPE(adsim::InertIsotherm)
CEREAL_REGISTER_TYPE(adsim::LangmuirIsotherm)
CEREAL_REGISTER_TYPE(adsim::AntiLangmuirIsotherm)
CEREAL_REGISTER_TYPE(adsim::FreundlichIsotherm)
CEREAL_REGISTER_TYPE(adsim::SipsIsotherm)
CEREAL_REGISTER_TYPE(adsim::Sips2Isotherm)
CEREAL_REGISTER_TYPE(adsim::ExtendedSips2Isotherm)
CEREAL_REGISTER_TYPE(adsim::DualSiteSipsIsotherm)
CEREAL_REGISTER_TYPE(adsim::RedlichPetersonIsotherm)
CEREAL_REGISTER_TYPE(adsim::TothIsotherm)
CEREAL_REGISTER_TYPE(adsim::Toth2Isotherm)
CEREAL_REGISTER_TYPE(adsim::ArrheniusSipsIsotherm)
CEREAL_REGISTER_TYPE(adsim::ExtendedLangmuirFreundlichIsotherm)
CEREAL_REGISTER_TYPE(adsim::ExtendedLangmuirIsotherm)
CEREAL_REGISTER_TYPE(adsim::DualSiteLangmuirIsotherm)
CEREAL_REGISTER_TYPE(adsim::ExtendedDualSiteLangmuirIsotherm)
CEREAL_REGISTER_TYPE(adsim::LangmuirIsotherm_Concentration)
CEREAL_REGISTER_TYPE(adsim::ExtendedLangmuirIsotherm_Concentration)
CEREAL_REGISTER_TYPE(adsim::DualSiteLangmuirIsotherm_Concentration)
CEREAL_REGISTER_TYPE(adsim::ExtendedDualSiteLangmuirIsotherm_Concentration)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IIsotherm, adsim::InertIsotherm)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IIsotherm, adsim::HenryIsotherm)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IIsotherm, adsim::LangmuirIsotherm)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IIsotherm, adsim::AntiLangmuirIsotherm)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IIsotherm, adsim::FreundlichIsotherm)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IIsotherm, adsim::SipsIsotherm)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IIsotherm, adsim::Sips2Isotherm)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IIsotherm, adsim::ExtendedSips2Isotherm)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IIsotherm, adsim::DualSiteSipsIsotherm)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IIsotherm, adsim::RedlichPetersonIsotherm)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IIsotherm, adsim::TothIsotherm)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IIsotherm, adsim::Toth2Isotherm)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IIsotherm, adsim::ArrheniusSipsIsotherm)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IIsotherm, adsim::ExtendedLangmuirFreundlichIsotherm)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IIsotherm, adsim::ExtendedLangmuirIsotherm)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IIsotherm, adsim::DualSiteLangmuirIsotherm)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IIsotherm, adsim::ExtendedDualSiteLangmuirIsotherm)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IIsotherm, adsim::LangmuirIsotherm_Concentration)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IIsotherm, adsim::ExtendedLangmuirIsotherm_Concentration)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IIsotherm, adsim::DualSiteLangmuirIsotherm_Concentration)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IIsotherm, adsim::ExtendedDualSiteLangmuirIsotherm_Concentration)
CEREAL_CLASS_VERSION(adsim::InertIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::HenryIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::LangmuirIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::AntiLangmuirIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::FreundlichIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::SipsIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::Sips2Isotherm, 2);
CEREAL_CLASS_VERSION(adsim::ExtendedSips2Isotherm, 2);
CEREAL_CLASS_VERSION(adsim::DualSiteSipsIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::RedlichPetersonIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::TothIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::Toth2Isotherm, 2);
CEREAL_CLASS_VERSION(adsim::ArrheniusSipsIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::ExtendedLangmuirFreundlichIsotherm, 1);
CEREAL_CLASS_VERSION(adsim::ExtendedLangmuirIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::DualSiteLangmuirIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::ExtendedDualSiteLangmuirIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::LangmuirIsotherm_Concentration, 2);
CEREAL_CLASS_VERSION(adsim::ExtendedLangmuirIsotherm_Concentration, 2);
CEREAL_CLASS_VERSION(adsim::DualSiteLangmuirIsotherm_Concentration, 2);
CEREAL_CLASS_VERSION(adsim::ExtendedDualSiteLangmuirIsotherm_Concentration, 2);
