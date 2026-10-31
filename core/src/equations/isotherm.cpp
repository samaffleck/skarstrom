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

	void IsothermFunction::Henry(double& qi, double T, double P, double yi, const HenryIsothermParameters& ips) {
		qi = ips.k * P * yi * std::exp(ips.Q / (constants::kIdealGasConstant * T));
	}

	void IsothermFunction::Langmuir(double& qi, double T, double P, double yi, const LangmuirIsothermParameters& ips) {
		const double e1 = std::exp(ips.q1 / (constants::kIdealGasConstant * T));
		qi = ips.m1 * ips.b1 * e1 * P * yi / (1 + ips.b1 * e1 * P * yi);
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
    double yi_pos = softplus(yi, kSoftplusTau);
    double b_i = ips.B1 * std::exp(ips.B2 / (R * T));
    double s_i = ips.C1 + ips.C2 * ((T / ips.T_ref) - 1);
    double z = safe_pow(b_i * P * yi_pos, s_i);
		double q_max = ips.A1 * std::exp(ips.A2 * ((T / ips.T_ref) - 1));
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
		double yi_pos = softplus(yi, kSoftplusTau);
		double z = safe_pow(k_i * P * yi_pos, c_i);
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
		double yi_pos = softplus(yi, kSoftplusTau);
    double qs = ips.q_sat * std::exp(ips.x * (1 - T / ips.T_0));
    double b = ips.b0 * std::exp(ips.dH0 / (R * ips.T_0) * (ips.T_0 / T - 1));
    double w = ips.w0 + ips.alpha * (1 - ips.T_0 / T);
    double bp = b * P * yi_pos;
    double bpv = safe_pow(bp, w);
		double denom = safe_pow(1 + bpv, (1 / w));
		qi = qs * bp / denom;
	}

	void IsothermFunction::DualSiteLangmuir(double& qi, double T, double P, double yi, const DualSiteLangmuirIsothermParameters& ips) {
		const double e1 = std::exp(ips.q1 / (constants::kIdealGasConstant * T));
		const double e2 = std::exp(ips.q2 / (constants::kIdealGasConstant * T));
		qi = (ips.m1 * ips.b1 * e1 * P * yi / (1 + ips.b1 * e1 * P * yi)) +
			(ips.m2 * ips.b2 * e2 * P * yi / (1 + ips.b2 * e2 * P * yi));
	}

	void IsothermFunction::Langmuir_Concentration(double& qi, double T, double P, double yi, const LangmuirIsothermParametersConcentration& ips) {
		double ci = P * yi / (constants::kIdealGasConstant * T);
		const double e1 = std::exp(ips.q1 / (constants::kIdealGasConstant * T));
		qi = ips.m1 * ips.b1 * e1 * ci / (1 + ips.b1 * e1 * ci);
	}

	void IsothermFunction::DualSiteLangmuir_Concentration(double& qi, double T, double P, double yi, const DualSiteLangmuirIsothermParametersConcentration& ips) {
		double ci = P * yi / (constants::kIdealGasConstant * T);
		const double e1 = std::exp(ips.q1 / (constants::kIdealGasConstant * T));
		const double e2 = std::exp(ips.q2 / (constants::kIdealGasConstant * T));
		qi = (ips.m1 * ips.b1 * e1 * ci / (1 + ips.b1 * e1 * ci)) +
			(ips.m2 * ips.b2 * e2 * ci / (1 + ips.b2 * e2 * ci));
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

    for (int n = start_index; n <= end_index; ++n) {
      IsothermFunction::Henry(qi_sat[n], T[n], P[n], yi[n], params);
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

    for (int n = start_index; n <= end_index; ++n) {
      IsothermFunction::Langmuir(qi_sat[n], T[n], P[n], yi[n], params);
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

    for (int n = start_index; n <= end_index; ++n) {
      IsothermFunction::Sips2(qi_sat[n], T[n], P[n], yi[n], params);
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

    for (int n = start_index; n <= end_index; ++n) {
      IsothermFunction::Toth2(qi_sat[n], T[n], P[n], yi[n], params);
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

    for (int n = start_index; n <= end_index; ++n) {
      IsothermFunction::DualSiteLangmuir(qi_sat[n], T[n], P[n], yi[n], params);
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

    for (int n = start_index; n <= end_index; ++n) {
      // Langmuir_Concentration calculates ci internally from yi, P, T
      IsothermFunction::Langmuir_Concentration(qi_sat[n], T[n], P[n], yi[n], params);
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

    for (int n = start_index; n <= end_index; ++n) {
      IsothermFunction::DualSiteLangmuir_Concentration(qi_sat[n], T[n], P[n], yi[n], params);
    }
	}

  void ExtendedLangmuirIsotherm::Update(IsothermUpdateParams& state) {
    auto& component = state.component;
    auto& qi_sat = state.qi_sat_map.at(component).var;
    auto& T = state.T;
    auto& P = state.P;
    auto& params = params_[component];
    int start_index = state.start_index;
    int end_index = state.end_index;
    const double R = constants::kIdealGasConstant;

    for (int n = start_index; n <= end_index; ++n) {
      // Calculate denominator over all components
      double denominator = 0.0;
      for (const auto& [comp, comp_params] : params_) {
        auto& yi = state.yi_map.at(comp).var;
        denominator += comp_params.b1 * std::exp(comp_params.q1 / (R * T[n])) * P[n] * yi[n];
      }

      // Calculate equilibrium adsorbed concentration for current component
      auto& yi = state.yi_map.at(component).var;
      double qsat = params.m1 * params.b1 * std::exp(params.q1 / (R * T[n])) * P[n] * yi[n] / (1 + denominator);
      qi_sat[n] = qsat;
    }
	}

	void ExtendedSips2Isotherm::Update(IsothermUpdateParams& state) {
    const double& R = constants::kIdealGasConstant;
    auto& component = state.component;
    auto& qi_sat = state.qi_sat_map.at(component).var;
    auto& T = state.T;
    auto& P = state.P;
    auto& params = params_[component];
    int start_index = state.start_index;
    int end_index = state.end_index;

    // Regularized competitive Sips loading.
    // z_raw = (b_i*P*mf_pos)^s uses the solver-stable softplus (tau=1e-6).
    // Multiplied by a smooth step mf_pos²/(mf_pos²+δ²) that → 0 as mf → 0 and → 1
    // for mf >> δ. This makes z(mf=0) effectively zero while suppressing the Jacobian
    // spike at mf=0. δ=1e-3 (not 1e-4) keeps the peak Jacobian gradient at mf=δ to
    // ~500 per unit vs ~5000 with δ=1e-4, preventing solver stalls in multi-cycle PSA
    // when CO2 mole fractions pass through the near-zero transition region.
    constexpr double kDelta = 1e-3;
    auto get_z = [&](double temp, double press, double mf, const Sips2IsothermParameters& params) -> double {
      const double b_i = params.B1 * std::exp(params.B2 / (R * temp));
      const double s_i = params.C1 + params.C2 * ((temp / params.T_ref) - 1);
      const double mf_pos = softplus(mf, kSoftplusTau);
      const double z_raw = safe_pow(b_i * press * mf_pos, s_i);
      return z_raw * (mf_pos * mf_pos) / (mf_pos * mf_pos + kDelta * kDelta);
    };

    for (int n = start_index; n <= end_index; ++n) {
      // Calculate denominator over all components
      double denominator = 0.0;
      for (const auto& [comp, comp_params] : params_) {
        auto& yi = state.yi_map.at(comp).var;
        denominator += get_z(T[n], P[n], yi[n], comp_params);
      }

      // Calculate equilibrium adsorbed concentration for current component
      auto& yi = state.yi_map.at(component).var;
      double q_max = params.A1 * std::exp(params.A2 * ((T[n] / params.T_ref) - 1));
      double z = get_z(T[n], P[n], yi[n], params);
      double qsat = q_max * z / (1 + denominator);
      qi_sat[n] = qsat;
    }
	}

  void ExtendedLangmuirIsotherm::UpdateAll(IsothermUpdateParams& state, const std::vector<std::string>& components) {
    auto& T = state.T;
    auto& P = state.P;
    int start_index = state.start_index;
    int end_index = state.end_index;
    const double R = constants::kIdealGasConstant;

    // Materialise params for every requested component first (matches the
    // operator[] insertion behaviour of the per-component Update).
    for (const auto& comp : components) params_[comp];

    // Hoist map lookups out of the cell loop. The denominator terms iterate
    // params_ in map order, exactly as the per-component Update does, so the
    // shared sum is bit-identical.
    struct Term { const LangmuirIsothermParameters* p; const std::vector<double>* yi; double e = 0.0; };
    std::vector<Term> terms;
    terms.reserve(params_.size());
    for (const auto& [comp, comp_params] : params_) {
      terms.push_back({ &comp_params, &state.yi_map.at(comp).var });
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
      // Calculate denominator over all components, caching the exponential so
      // each output component reuses it with the original operation order.
      double denominator = 0.0;
      for (auto& t : terms) {
        t.e = std::exp(t.p->q1 / (R * T[n]));
        denominator += t.p->b1 * t.e * P[n] * (*t.yi)[n];
      }

      for (const auto& out : outputs) {
        const double e = terms[out.term_index].e;
        (*out.qi_sat)[n] = out.p->m1 * out.p->b1 * e * P[n] * (*out.yi)[n] / (1 + denominator);
      }
    }
  }

	void ExtendedDualSiteLangmuirIsotherm::Update(IsothermUpdateParams& state) {
    auto& component = state.component;
    auto& qi_sat = state.qi_sat_map.at(component).var;
    auto& T = state.T;
    auto& P = state.P;
    auto& params = params_[component];
    int start_index = state.start_index;
    int end_index = state.end_index;
    const double R = constants::kIdealGasConstant;

    for (int n = start_index; n <= end_index; ++n) {
      // Calculate denominators over all components
      double denominator1 = 0.0;
      double denominator2 = 0.0;
      for (const auto& [comp, comp_params] : params_) {
        auto& yi = state.yi_map.at(comp).var;
        denominator1 += comp_params.b1 * std::exp(comp_params.q1 / (R * T[n])) * P[n] * yi[n];
        denominator2 += comp_params.b2 * std::exp(comp_params.q2 / (R * T[n])) * P[n] * yi[n];
      }

      // Calculate equilibrium adsorbed concentration for current component
      auto& yi = state.yi_map.at(component).var;
      double qsat =
        params.m1 * params.b1 * std::exp(params.q1 / (R * T[n])) * P[n] * yi[n] / (1 + denominator1) +
        params.m2 * params.b2 * std::exp(params.q2 / (R * T[n])) * P[n] * yi[n] / (1 + denominator2);
      qi_sat[n] = qsat;
    }
	}

  void ExtendedDualSiteLangmuirIsotherm::UpdateAll(IsothermUpdateParams& state, const std::vector<std::string>& components) {
    auto& T = state.T;
    auto& P = state.P;
    int start_index = state.start_index;
    int end_index = state.end_index;
    const double R = constants::kIdealGasConstant;

    // Materialise params for every requested component first (matches the
    // operator[] insertion behaviour of the per-component Update).
    for (const auto& comp : components) params_[comp];

    // Hoist map lookups out of the cell loop. The denominator terms iterate
    // params_ in map order, exactly as the per-component Update does, so the
    // shared sums are bit-identical.
    struct Term { const DualSiteLangmuirIsothermParameters* p; const std::vector<double>* yi; double e1 = 0.0; double e2 = 0.0; };
    std::vector<Term> terms;
    terms.reserve(params_.size());
    for (const auto& [comp, comp_params] : params_) {
      terms.push_back({ &comp_params, &state.yi_map.at(comp).var });
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
      // Calculate denominators over all components, caching the exponentials
      // so each output component reuses them with the original operation order.
      double denominator1 = 0.0;
      double denominator2 = 0.0;
      for (auto& t : terms) {
        t.e1 = std::exp(t.p->q1 / (R * T[n]));
        t.e2 = std::exp(t.p->q2 / (R * T[n]));
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

  void ExtendedSips2Isotherm::UpdateAll(IsothermUpdateParams& state, const std::vector<std::string>& components) {
    const double& R = constants::kIdealGasConstant;
    auto& T = state.T;
    auto& P = state.P;
    int start_index = state.start_index;
    int end_index = state.end_index;

    // Materialise params for every requested component first (matches the
    // operator[] insertion behaviour of the per-component Update).
    for (const auto& comp : components) params_[comp];

    // Same regularised competitive Sips loading as the per-component Update —
    // see the comment there for the smooth-step rationale.
    constexpr double kDelta = 1e-3;
    auto get_z = [&](double temp, double press, double mf, const Sips2IsothermParameters& params) -> double {
      const double b_i = params.B1 * std::exp(params.B2 / (R * temp));
      const double s_i = params.C1 + params.C2 * ((temp / params.T_ref) - 1);
      const double mf_pos = softplus(mf, kSoftplusTau);
      const double z_raw = safe_pow(b_i * press * mf_pos, s_i);
      return z_raw * (mf_pos * mf_pos) / (mf_pos * mf_pos + kDelta * kDelta);
    };

    // Hoist map lookups out of the cell loop. The denominator terms iterate
    // params_ in map order, exactly as the per-component Update does, so the
    // shared sum is bit-identical; each output component reuses its own z.
    struct Term { const Sips2IsothermParameters* p; const std::vector<double>* yi; double z = 0.0; };
    std::vector<Term> terms;
    terms.reserve(params_.size());
    for (const auto& [comp, comp_params] : params_) {
      terms.push_back({ &comp_params, &state.yi_map.at(comp).var });
    }

    struct Output { const Sips2IsothermParameters* p; std::vector<double>* qi_sat; size_t term_index; };
    std::vector<Output> outputs;
    outputs.reserve(components.size());
    for (const auto& comp : components) {
      const auto* p = &params_.at(comp);
      size_t idx = 0;
      while (idx < terms.size() && terms[idx].p != p) ++idx;
      outputs.push_back({ p, &state.qi_sat_map.at(comp).var, idx });
    }

    for (int n = start_index; n <= end_index; ++n) {
      // Calculate denominator over all components
      double denominator = 0.0;
      for (auto& t : terms) {
        t.z = get_z(T[n], P[n], (*t.yi)[n], *t.p);
        denominator += t.z;
      }

      for (const auto& out : outputs) {
        double q_max = out.p->A1 * std::exp(out.p->A2 * ((T[n] / out.p->T_ref) - 1));
        double z = terms[out.term_index].z;
        (*out.qi_sat)[n] = q_max * z / (1 + denominator);
      }
    }
  }

  // --- Arrhenius Sips (Schell/Casas 2012) -----------------------------------
  // Shared helpers: evaluate the temperature-dependent Sips parameters and the
  // regularised competitive term z_i = (k_i·P·y_i)^{c_i}. The regularisation is
  // identical to ExtendedSips2Isotherm (softplus + smooth step, kDelta = 1e-3)
  // so the Jacobian stays well conditioned as a mole fraction passes through 0.
  namespace {
    constexpr double kArrSipsDelta = 1e-3;

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
      const double mf_pos = softplus(mf, kSoftplusTau);
      const double z_raw = safe_pow(k_i * P * mf_pos, c_i);
      return z_raw * (mf_pos * mf_pos) / (mf_pos * mf_pos + kArrSipsDelta * kArrSipsDelta);
    }
  } // namespace

  void ArrheniusSipsIsotherm::Update(IsothermUpdateParams& state) {
    auto& component = state.component;
    auto& qi_sat = state.qi_sat_map.at(component).var;
    auto& T = state.T;
    auto& P = state.P;
    auto& params = params_[component];
    int start_index = state.start_index;
    int end_index = state.end_index;

    for (int n = start_index; n <= end_index; ++n) {
      double denominator = 0.0;
      for (const auto& [comp, comp_params] : params_) {
        auto& yi = state.yi_map.at(comp).var;
        denominator += ArrSips_z(T[n], P[n], yi[n], comp_params);
      }
      auto& yi = state.yi_map.at(component).var;
      const double n_inf = ArrSips_ninf(T[n], params);
      const double z = ArrSips_z(T[n], P[n], yi[n], params);
      qi_sat[n] = n_inf * z / (1 + denominator);
    }
  }

  void ArrheniusSipsIsotherm::UpdateAll(IsothermUpdateParams& state,
                                        const std::vector<std::string>& components) {
    auto& T = state.T;
    auto& P = state.P;
    int start_index = state.start_index;
    int end_index = state.end_index;

    for (const auto& comp : components) params_[comp];

    // Hoist map lookups; iterate params_ in map order for a shared denominator
    // that is bit-identical to the per-component Update.
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

	void ExtendedLangmuirIsotherm_Concentration::Update(IsothermUpdateParams& state) {
    auto& component = state.component;
    auto& qi_sat = state.qi_sat_map.at(component).var;
    auto& T = state.T;
    auto& P = state.P;
    auto& params = params_[component];
    int start_index = state.start_index;
    int end_index = state.end_index;
    const double R = constants::kIdealGasConstant;

    for (int n = start_index; n <= end_index; ++n) {
      // Calculate denominator over all components (using concentration)
      double denominator = 0.0;
      for (const auto& [comp, comp_params] : params_) {
        auto& yi = state.yi_map.at(comp).var;
        double ci = P[n] * yi[n] / (R * T[n]);
        denominator += comp_params.b1 * std::exp(comp_params.q1 / (R * T[n])) * ci;
      }

      // Calculate equilibrium adsorbed concentration for current component
      auto& yi = state.yi_map.at(component).var;
      double ci = P[n] * yi[n] / (R * T[n]);
      double qsat = params.m1 * params.b1 * std::exp(params.q1 / (R * T[n])) * ci / (1 + denominator);
      qi_sat[n] = qsat;
    }
	}

  void ExtendedDualSiteLangmuirIsotherm_Concentration::Update(IsothermUpdateParams& state) {
    auto& component = state.component;
    auto& qi_sat = state.qi_sat_map.at(component).var;
    auto& T = state.T;
    auto& P = state.P;
    auto& params = params_[component];
    int start_index = state.start_index;
    int end_index = state.end_index;
    const double R = constants::kIdealGasConstant;

    for (int n = start_index; n <= end_index; ++n) {
      // Calculate denominators over all components (using concentration)
      double denominator1 = 0.0;
      double denominator2 = 0.0;
      for (const auto& [comp, comp_params] : params_) {
        auto& yi = state.yi_map.at(comp).var;
        double ci = P[n] * yi[n] / (R * T[n]);
        denominator1 += comp_params.b1 * std::exp(comp_params.q1 / (R * T[n])) * ci;
        denominator2 += comp_params.b2 * std::exp(comp_params.q2 / (R * T[n])) * ci;
      }

      // Calculate equilibrium adsorbed concentration for current component
      auto& yi = state.yi_map.at(component).var;
      double ci = P[n] * yi[n] / (R * T[n]);
      double qsat =
        params.m1 * params.b1 * std::exp(params.q1 / (R * T[n])) * ci / (1 + denominator1) +
        params.m2 * params.b2 * std::exp(params.q2 / (R * T[n])) * ci / (1 + denominator2);
      qi_sat[n] = qsat;
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
CEREAL_CLASS_VERSION(adsim::ExtendedLangmuirIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::DualSiteLangmuirIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::ExtendedDualSiteLangmuirIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::LangmuirIsotherm_Concentration, 2);
CEREAL_CLASS_VERSION(adsim::ExtendedLangmuirIsotherm_Concentration, 2);
CEREAL_CLASS_VERSION(adsim::DualSiteLangmuirIsotherm_Concentration, 2);
CEREAL_CLASS_VERSION(adsim::ExtendedDualSiteLangmuirIsotherm_Concentration, 2);
