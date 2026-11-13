/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/equations/mass_transfer_model.h"
#include "adsim/equations/equation.h"
#include "adsim/equations/isotherm.h"
#include "adsim/objects/reactor.h"
#include "adsim/objects/porous_media.h"
#include "adsim/objects/model_kinds.h"

#include <stdexcept>


namespace adsim {

  void MassTransferCoefficientConstant::UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) {
    auto& mtc = state.mass_transfer_coefficient.at(ads_id).at(component).var;
    auto& T = state.fluid_temperature.var;
    auto& P = state.pressure.var;
    auto& y = state.mole_fraction.at(component).var;

    if (constant_mtc_.IsRunConstant()) {
      // State-independent: one evaluation for the layer (same value every cell).
      EvalContext ctx{registry, T[layer.domain_.start_index], P[layer.domain_.start_index], 0.0, y[layer.domain_.start_index], 0.0};
      const double k = constant_mtc_.EvaluateAt(ctx);
      for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) mtc[n] = k;
      return;
    }
    for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
      EvalContext ctx{registry, T[n], P[n], 0.0, y[n], 0.0};
      mtc[n] = constant_mtc_.EvaluateAt(ctx);
    }
  }

  void MassTransferCoefficientArrhenius::UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) {
    auto& mtc = state.mass_transfer_coefficient.at(ads_id).at(component).var;
    auto& T = state.fluid_temperature.var;
    auto& P = state.pressure.var;
    auto& y = state.mole_fraction.at(component).var;
    double R = constants::kIdealGasConstant;

    for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
      EvalContext ctx{registry, T[n], P[n], 0.0, y[n], 0.0};
      const double k0 = pre_exp_factor_.EvaluateAt(ctx);
      const double Ea = activation_energy_.EvaluateAt(ctx);
      mtc[n] = k0 * std::exp(-Ea / (R * T[n]));
    }
  }

  void MassTransferCoefficientMacroporeKnudsen::UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) {
    auto& mtc = state.mass_transfer_coefficient.at(ads_id).at(component).var;
    auto& T = state.fluid_temperature.var;
    auto& Dmi = state.molecular_diffusion.at(component).var;
    double R = constants::kIdealGasConstant;

    double Dki = 0;
    double Dki_eff = 0;
    double Dmi_eff = 0;
    double Di_eff = 0;

    double ep = macropore_void_fraction_;
    // Glueckauf (1955) LDF: k = c D_eff / l^2, c = 60 on d_p for spheres
    // (packed), 3 on the half-thickness for a slab (structured layer).
    const auto geo = layer.GetDiffusionGeometry();
    const double c_ldf = geo.c;
    double dp = geo.length;
    double tortuosity = tortuosity_;
    double Mr = fluid.molar_mass_.at(component);
    const double Dm_const = molecular_diffusion_coefficient_;  // > 0: constant D_m in place of the state's

    int N = mtc.size();
    for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
      // Newton iterates can probe T < 0; skip just this cell rather than
      // abandoning the sweep and leaving the remaining cells stale.
      if (T[n] < 0)
        continue;

      Dki = (2.0 / 3.0) * pore_radius_ * sqrt(8.0 * R * T[n] / (constants::kPI * Mr));
      Dki_eff = Dki * ep / tortuosity;
      Dmi_eff = (Dm_const > 0.0 ? Dm_const : Dmi[n]) * ep / tortuosity;
      Di_eff = 1.0 / ((1.0 / Dki_eff) + (1.0 / Dmi_eff));
      mtc[n] = c_ldf * Di_eff / (dp * dp);
    }
  }

  void MassTransferCoefficientMacropore::UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) {
    auto& mtc = state.mass_transfer_coefficient.at(ads_id).at(component).var;
    auto& T = state.fluid_temperature.var;
    auto& Dmi = state.molecular_diffusion.at(component).var;
    double ep = macropore_void_fraction_;
    const auto geo = layer.GetDiffusionGeometry();  // Glueckauf (1955): sphere 60 / d_p^2, slab 3 / l^2
    const double c_ldf = geo.c;
    double dp = geo.length;
    double tortuosity = tortuosity_;
    const double Dm_const = molecular_diffusion_coefficient_;  // > 0: constant D_m in place of the state's

    int N = mtc.size();
    for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
      // Newton iterates can probe T < 0; skip just this cell rather than
      // abandoning the sweep and leaving the remaining cells stale.
      if (T[n] < 0)
        continue;

      // The plain macropore coefficient, 15 eps_p D_m / (tau r_p^2); the
      // capacity-factor variant is MassTransferCoefficientMacroporeCapacity.
      double Dmi_eff = (Dm_const > 0.0 ? Dm_const : Dmi[n]) * ep / tortuosity;
      mtc[n] = c_ldf * Dmi_eff / (dp * dp);
    }
  }

  void MassTransferCoefficientMacroporeCapacity::UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) {
    auto& mtc = state.mass_transfer_coefficient.at(ads_id).at(component).var;
    const auto& T = state.fluid_temperature.var;
    const auto& P = state.pressure.var;
    const auto& Dmi = state.molecular_diffusion.at(component).var;
    auto& qi_sat = state.eq_ads_concentration.at(ads_id).at(component).var;
    auto& yi = state.mole_fraction.at(component).var;
    const double R = constants::kIdealGasConstant;

    const double ep = macropore_void_fraction_;
    const auto geo = layer.GetDiffusionGeometry();  // Glueckauf (1955): sphere 60 / d_p^2, slab 3 / l^2
    const double c_ldf = geo.c;
    const double dp = geo.length;
    const double tortuosity = tortuosity_;
    const double c_floor = std::max(concentration_floor_.Get(), 0.0);
    const double rho_p = solid_density_;
    const double Dm_const = molecular_diffusion_coefficient_;  // > 0: constant D_m in place of the state's
    const int start = layer.domain_.start_index, end = layer.domain_.end_index;

    // q* at max(c_i, c_floor) from the isotherm itself, evaluated here so it is
    // consistent with the current y, P and T whatever the residual's update
    // order (the state's own q* may be the previous evaluation's). The mole
    // fraction of this component is floored in place, the isotherm run for
    // this component only, the result kept, and both vectors restored.
    const bool evaluate = isotherm_ && rho_p > 0.0;
    if (evaluate) {
      y_scratch_.assign(yi.begin(), yi.end());
      q_scratch_.assign(qi_sat.begin(), qi_sat.end());
      for (int n = start; n <= end; ++n) {
        if (!(T[n] > 0.0) || !(P[n] > 0.0)) continue;
        const double y_floor = c_floor * R * T[n] / P[n];
        if (yi[n] < y_floor) yi[n] = y_floor;
      }
      if (source_ && source_->mixing_rule_ == MixingRule::IAST) {
        // IAST couples the components: re-solve the whole adsorbent at the
        // floored state, keep this component's loading, restore every other one.
        auto& q_map = state.eq_ads_concentration.at(ads_id);
        const auto& components = fluid.components_;
        q_all_scratch_.resize(components.size());
        for (size_t i = 0; i < components.size(); ++i) {
          const auto& q = q_map.at(components[i]).var;
          q_all_scratch_[i].assign(q.begin(), q.end());
        }
        source_->UpdateIsothermAll(state, ads_id, components, start, end, registry);
        std::swap(q_scratch_, qi_sat);
        for (size_t i = 0; i < components.size(); ++i) {
          auto& q = q_map.at(components[i]).var;
          if (&q == &qi_sat) continue;  // restored by the swap above
          std::copy(q_all_scratch_[i].begin(), q_all_scratch_[i].end(), q.begin());
        }
      } else {
        // The adsorbent's loadings, for the kinds that read another
        // component's loading (null when the state has none: their dry form).
        const auto q_it = state.ads_concentration.find(ads_id);
        IsothermUpdateParams params(state.mole_fraction, state.eq_ads_concentration.at(ads_id),
                                    state.solid_temperature.var, state.pressure.var, component, start, end, registry,
                                    q_it == state.ads_concentration.end() ? nullptr : &q_it->second);
        isotherm_->Update(params);
        std::swap(q_scratch_, qi_sat);   // q_scratch_ = q* at the floored state; qi_sat = the original, restored
      }
      std::copy(y_scratch_.begin(), y_scratch_.end(), yi.begin());
    }

    for (int n = start; n <= end; ++n) {
      if (T[n] < 0) continue;  // Newton iterates can probe T < 0; skip the cell
      const double Dmi_eff = (Dm_const > 0.0 ? Dm_const : Dmi[n]) * ep / tortuosity;
      const double k_macro = c_ldf * Dmi_eff / (dp * dp);
      double ratio = 1.0;
      if (evaluate && T[n] > 0.0 && P[n] > 0.0) {
        // Gas-phase concentration from the state itself (y P / R T), floored
        // like the isotherm evaluation above.
        const double c = std::max(y_scratch_[n], 0.0) * P[n] / (R * T[n]);
        const double c_eff = std::max(c, c_floor);
        const double denom = rho_p * q_scratch_[n];
        ratio = denom > 0.0 ? std::min(c_eff / denom, 1.0) : 1.0;
      }
      mtc[n] = k_macro * ratio;
    }
  }

  void MassTransferCoefficient_Haynes::UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) {
    auto& mtc = state.mass_transfer_coefficient.at(ads_id).at(component).var;
    auto& T = state.fluid_temperature.var;
    auto& Dmi = state.molecular_diffusion.at(component).var;
    auto& Sh = state.sherwood_number.var;
    double R = constants::kIdealGasConstant;

    double ep = macropore_void_fraction_;
    // 1/k = (V/S) / k_f + tau l^2 / (c' eps_p D): film plus pore diffusion in
    // series. Packed spheres: V/S = r_p / 3, l = r_p, c' = 15 (Glueckauf
    // 1955), k_f = Sh D / d_p with the layer's film correlation Sh
    // (film_correlation.h; default Wakao-Funazkri). A slab (structured
    // layer): V/S = l, the half-thickness, c' = 3, k_f = Sh D / D_h with the
    // channel Sh.
    // A packed bed of infinite cylinders: V/S = r_f / 2, l = r_f, c' = 8,
    // k_f = Sh D / d_p (the packed correlation on d_p).
    if (layer.IsStructured() || layer.particle_shape_ == ParticleShape::Cylinder) {
      const auto geo = layer.GetDiffusionGeometry();
      const double l = geo.length;
      const double c = geo.c;
      const double dh = layer.IsStructured() ? layer.GetHydraulicDiameter() : layer.particle_diameter_.Get();
      for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
        const double kf = Sh[n] * Dmi[n] / dh;
        const double film = geo.volume_to_surface / kf;
        const double molecular = tortuosity_ * l * l / (c * ep * Dmi[n]);
        mtc[n] = 1.0 / (film + molecular);
      }
      return;
    }
    double dp = layer.particle_diameter_;
    double rp = dp / 2.0;

    int N = mtc.size();
    for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
      double kf = Sh[n] * Dmi[n] / dp;
      double film = rp / (3.0 * kf);
      double molecular = tortuosity_ * rp * rp / (15.0 * ep * Dmi[n]);
      mtc[n] = 1.0 / (film + molecular);
    }
  }

  void LegacySeriesResistanceV1::EvaluateLoading(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, int start, int end, const UserVariableRegistry* registry, std::vector<double>& out) {
    auto& qi_sat = state.eq_ads_concentration.at(ads_id).at(component).var;
    out.assign(qi_sat.begin(), qi_sat.end());
    if (source_ && source_->mixing_rule_ == MixingRule::IAST) {
      // IAST couples the components: solve the whole adsorbent, keep this
      // component's loading, restore every other one.
      auto& q_map = state.eq_ads_concentration.at(ads_id);
      const auto& components = fluid.components_;
      q_all_scratch_.resize(components.size());
      for (size_t i = 0; i < components.size(); ++i) {
        const auto& q = q_map.at(components[i]).var;
        q_all_scratch_[i].assign(q.begin(), q.end());
      }
      source_->UpdateIsothermAll(state, ads_id, components, start, end, registry);
      std::swap(out, qi_sat);  // out = the new q*; qi_sat = the original, restored
      for (size_t i = 0; i < components.size(); ++i) {
        auto& q = q_map.at(components[i]).var;
        if (&q == &qi_sat) continue;
        std::copy(q_all_scratch_[i].begin(), q_all_scratch_[i].end(), q.begin());
      }
      return;
    }
    const auto q_it = state.ads_concentration.find(ads_id);
    IsothermUpdateParams params(state.mole_fraction, state.eq_ads_concentration.at(ads_id),
                                state.solid_temperature.var, state.pressure.var, component, start, end, registry,
                                q_it == state.ads_concentration.end() ? nullptr : &q_it->second);
    isotherm_->Update(params);
    std::swap(out, qi_sat);  // out = the new q*; qi_sat = the original, restored
  }

  void LegacySeriesResistanceV1::CapacityFactor(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, int start, int end, const UserVariableRegistry* registry, std::vector<double>& k_out) {
    auto& yi = state.mole_fraction.at(component).var;
    const auto& T = state.fluid_temperature.var;
    const auto& P = state.pressure.var;
    const double R = constants::kIdealGasConstant;
    k_out.assign(yi.size(), 1.0);

    const double fixed = capacity_factor_.Get();
    if (fixed > 0.0) {
      std::fill(k_out.begin() + start, k_out.begin() + end + 1, fixed);
      return;
    }
    const double rho_p = solid_density_;
    if (!isotherm_ || !(rho_p > 0.0)) return;  // no isotherm context: K = 1

    // Forward difference of q* in c at c_eff = max(c, c_floor): this
    // component's mole fraction floored, then raised by the relative step h,
    // every other component at its actual value. Both vectors restored.
    constexpr double h = 1e-6;
    const double c_floor = std::max(concentration_floor_.Get(), 0.0);
    auto valid = [&](int n) { return T[n] > 0.0 && P[n] > 0.0; };
    auto floored = [&](int n) {
      const double y_floor = c_floor * R * T[n] / P[n];
      return yi[n] < y_floor ? y_floor : yi[n];
    };
    y_scratch_.assign(yi.begin(), yi.end());
    for (int n = start; n <= end; ++n)
      if (valid(n)) yi[n] = floored(n);
    EvaluateLoading(state, ads_id, fluid, component, start, end, registry, q_base_);
    for (int n = start; n <= end; ++n)
      if (valid(n)) yi[n] = yi[n] * (1.0 + h);
    EvaluateLoading(state, ads_id, fluid, component, start, end, registry, q_pert_);
    std::copy(y_scratch_.begin(), y_scratch_.end(), yi.begin());

    for (int n = start; n <= end; ++n) {
      if (!valid(n)) continue;
      const double y_eff = floored(n);
      const double y_pert = y_eff * (1.0 + h);
      const double to_c = P[n] / (R * T[n]);
      const double c_eff = y_eff * to_c;
      const double dc = (y_pert - y_eff) * to_c;
      double K = dc > 0.0 ? rho_p * (q_pert_[n] - q_base_[n]) / dc : 0.0;    // slope
      if (!(K > 0.0) || !std::isfinite(K)) K = c_eff > 0.0 ? rho_p * q_base_[n] / c_eff : 0.0;  // secant
      if (!(K >= 1.0) || !std::isfinite(K)) K = 1.0;  // floor: never above the K = 1 coefficient
      k_out[n] = K;
    }
  }

  void LegacySeriesResistanceV1::UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) {
    auto& mtc = state.mass_transfer_coefficient.at(ads_id).at(component).var;
    const auto& T = state.fluid_temperature.var;
    const auto& Dmi = state.molecular_diffusion.at(component).var;
    const auto& Sh = state.sherwood_number.var;
    const double R = constants::kIdealGasConstant;
    const int start = layer.domain_.start_index, end = layer.domain_.end_index;

    CapacityFactor(state, ads_id, fluid, component, start, end, registry, capacity_);

    const double f_film = film_resistance_factor_;
    const double tau = tortuosity_;
    const double ep = macropore_void_fraction_;
    const double Dm_const = molecular_diffusion_coefficient_;  // > 0: constant D_m in the macropore term
    const double rc = crystal_radius_;
    const double Dc0 = micropore_diffusivity_;
    const double Ea = micropore_activation_energy_;
    const bool film_on = f_film > 0.0;
    const bool macro_on = tau > 0.0;
    const bool micro_on = rc > 0.0 && Dc0 > 0.0;

    // The terms are written as in MassTransferCoefficient_Haynes so that with
    // f_film = 1, K = 1 and no micropore term the result is the same double.
    const bool structured = layer.IsStructured();
    // A packed bed of infinite cylinders takes the geometry's form (V/S = r_f / 2,
    // l = r_f, c' = 8) with the packed k_f = Sh D / d_p; spheres keep their own
    // expressions below.
    const bool cylinder = !structured && layer.particle_shape_ == ParticleShape::Cylinder;
    const auto geo = layer.GetDiffusionGeometry();
    const double l = geo.length;
    const double c_geo = geo.c;
    const double dh = structured ? layer.GetHydraulicDiameter() : 0.0;
    const double dp = layer.particle_diameter_;
    const double rp = dp / 2.0;
    for (int n = start; n <= end; ++n) {
      if (micro_on && !(T[n] > 0.0)) continue;  // Newton iterates can probe T < 0; skip the cell
      const double K = capacity_[n];
      const double Dp = Dm_const > 0.0 ? Dm_const : Dmi[n];
      double film = 0.0, macro = 0.0, micro = 0.0;
      if (structured) {
        if (film_on) {
          const double kf = Sh[n] * Dmi[n] / dh;
          film = geo.volume_to_surface / kf;
        }
        if (macro_on) macro = tau * l * l / (3.0 * ep * Dp);
      } else if (cylinder) {
        if (film_on) {
          const double kf = Sh[n] * Dmi[n] / dp;
          film = geo.volume_to_surface / kf;
        }
        if (macro_on) macro = tau * l * l / (c_geo * ep * Dp);
      } else {
        if (film_on) {
          const double kf = Sh[n] * Dmi[n] / dp;
          film = rp / (3.0 * kf);
        }
        if (macro_on) macro = tau * rp * rp / (15.0 * ep * Dp);
      }
      if (micro_on) micro = rc * rc / (15.0 * Dc0 * std::exp(-Ea / (R * T[n])));
      mtc[n] = 1.0 / ((f_film * K) * film + K * macro + micro);
    }
  }

  // ---- the generic coefficient model ----------------------------------------------------

  namespace {

    using Generic = MassTransferCoefficientSeriesResistance;

    // The switch set with the terms named on; the macropore term holds
    // molecular diffusion (the former kinds' form).
    MassTransferResistances Only(bool film, bool macropore, bool knudsen, bool micropore, bool lumped) {
      MassTransferResistances r;
      r.film = film; r.macropore = macropore; r.molecular = macropore; r.knudsen = knudsen; r.micropore = micropore;
      r.lumped = lumped;
      return r;
    }

    // True when the parameter is the plain number v (no expression).
    bool IsNumber(const ExpressionParam& p, double v) { return !p.HasExpression() && p.Get() == v; }

    bool* SwitchOf(MassTransferResistances& r, const std::string& tag) {
      if (tag == "film") return &r.film;
      if (tag == "macropore") return &r.macropore;
      if (tag == "molecular") return &r.molecular;
      if (tag == "knudsen") return &r.knudsen;
      if (tag == "micropore") return &r.micropore;
      if (tag == "lumped") return &r.lumped;
      return nullptr;
    }

  }  // namespace

  bool MassTransferCoefficientSeriesResistance::ActivationEnergyIsZero() const {
    return IsNumber(activation_energy_, 0.0);
  }

  bool MassTransferCoefficientSeriesResistance::Enabled(const std::string& tag) const {
    MassTransferResistances r = Resistances();
    const bool* on = SwitchOf(r, tag);
    return on && *on;
  }

  bool& MassTransferCoefficientSeriesResistance::Switch(const std::string& tag) {
    bool* on = SwitchOf(resistances_, tag);
    if (!on) throw std::out_of_range("unknown mass-transfer resistance '" + tag + "'");
    return *on;
  }

  CapacityFactorMode MassTransferCoefficientSeriesResistance::CapacityMode() const {
    if (capacity_mode_ < 0 || capacity_mode_ > static_cast<int>(CapacityFactorMode::FIXED)) return CapacityFactorMode::NONE;
    return static_cast<CapacityFactorMode>(capacity_mode_);
  }

  void MassTransferCoefficientSeriesResistance::SetCapacityMode(CapacityFactorMode mode) {
    capacity_mode_ = static_cast<int>(mode);
    if (mode == CapacityFactorMode::FIXED && !(capacity_factor_.Get() > 0.0)) capacity_factor_ = 1.0;
  }

  bool MassTransferCoefficientSeriesResistance::SetResistance(const std::string& tag, bool on, std::string& err) {
    MassTransferResistances next = Resistances();
    bool* sw = SwitchOf(next, tag);
    if (!sw) {
      err = "unknown mass-transfer resistance '" + tag + "' (film, macropore, molecular, knudsen, micropore, lumped)";
      return false;
    }
    const bool pore_term = tag == "molecular" || tag == "knudsen";
    if (on && pore_term && !next.macropore) {
      err = std::string(tag == "knudsen" ? "Knudsen" : "molecular") +
            " diffusion refines the macropore term: switch the macropore resistance on first";
      return false;
    }
    if (!on && pore_term && !next.macropore) return true;  // already off with the macropore term
    *sw = on;
    if (tag == "macropore" && !on) next.molecular = next.knudsen = false;
    if (tag == "macropore" && on && !next.molecular && !next.knudsen) next.molecular = true;
    if (next.macropore && !next.molecular && !next.knudsen) {
      err = "the macropore term needs molecular or Knudsen diffusion: at least one of them must stay on";
      return false;
    }
    if (!(next.film || next.macropore || next.micropore || next.lumped)) {
      err = "at least one mass-transfer resistance must stay on";
      return false;
    }
    const bool first_enable = on && !Enabled(tag);
    resistances_ = next;
    if (first_enable) {
      // Out-of-range values of the term just switched on take its defaults.
      if (tag == "film" && !(film_resistance_factor_.Get() > 0.0)) film_resistance_factor_ = 1.0;
      if (tag == "macropore") {
        if (!(tortuosity_.Get() > 0.0)) tortuosity_ = 3.0;
        const double ep = macropore_void_fraction_.Get();
        if (!(ep > 0.0 && ep <= 1.0)) macropore_void_fraction_ = 0.35;
      }
      if (tag == "knudsen" && !(pore_radius_.Get() > 0.0)) pore_radius_ = 1e-9;
      if (tag == "micropore") {
        if (!(crystal_radius_.Get() > 0.0)) crystal_radius_ = 1e-6;
        if (!(micropore_diffusivity_.Get() > 0.0)) micropore_diffusivity_ = 1e-13;
      }
    }
    return true;
  }

  std::vector<std::string> MassTransferCoefficientSeriesResistance::CapacityParameters() const {
    if (!UsesCapacityFactor()) return {};
    switch (CapacityMode()) {
      case CapacityFactorMode::SLOPE:
      case CapacityFactorMode::SECANT: return {"concentration_floor"};
      case CapacityFactorMode::FIXED: return {"capacity_factor"};
      default: return {};
    }
  }

  const std::vector<std::string>& MassTransferCoefficientSeriesResistance::ResistanceParameters(const std::string& tag) {
    static const std::vector<std::string> none;
    static const std::vector<std::string> film{"film_resistance_factor"};
    static const std::vector<std::string> macro{"tortuosity", "macropore_void_fraction"};
    static const std::vector<std::string> molecular{"molecular_diffusion_coefficient"};
    static const std::vector<std::string> knudsen{"pore_radius"};
    static const std::vector<std::string> micro{"crystal_radius", "micropore_diffusivity", "micropore_activation_energy"};
    static const std::vector<std::string> lumped{"k", "activation_energy"};
    if (tag == "film") return film;
    if (tag == "macropore") return macro;
    if (tag == "molecular") return molecular;
    if (tag == "knudsen") return knudsen;
    if (tag == "micropore") return micro;
    if (tag == "lumped") return lumped;
    return none;
  }

  void MassTransferCoefficientSeriesResistance::AdoptVersion1Switches() {
    // Version 1: a term is off when its parameter is zero; K from the slope
    // unless capacity_factor > 0.
    const MassTransferResistances r = Only(film_resistance_factor_.Get() > 0.0, tortuosity_.Get() > 0.0, false,
                                           crystal_radius_.Get() > 0.0 && micropore_diffusivity_.Get() > 0.0, false);
    resistances_ = r;
    capacity_mode_ = static_cast<int>(capacity_factor_.Get() > 0.0 ? CapacityFactorMode::FIXED : CapacityFactorMode::SLOPE);
    // A switched-off term's zero sentinels become its defaults.
    if (!r.film && !(film_resistance_factor_.Get() > 0.0)) film_resistance_factor_ = 1.0;
    if (!r.macropore && !(tortuosity_.Get() > 0.0)) tortuosity_ = 3.0;
    if (!r.micropore) {
      if (!(crystal_radius_.Get() > 0.0)) crystal_radius_ = 1e-6;
      if (!(micropore_diffusivity_.Get() > 0.0)) micropore_diffusivity_ = 1e-13;
    }
  }

  void MassTransferCoefficientSeriesResistance::ApplyPreset(MassTransferCoefficientType preset) {
    if (preset == MassTransferCoefficientType::SERIES_RESISTANCE) return;  // custom: as it is
    capacity_mode_ = static_cast<int>(CapacityFactorMode::NONE);
    switch (preset) {
      case MassTransferCoefficientType::CONSTANT:
        resistances_ = Only(false, false, false, false, true);
        activation_energy_ = 0.0;
        break;
      case MassTransferCoefficientType::ARRHENIUS:
        resistances_ = Only(false, false, false, false, true);
        break;
      case MassTransferCoefficientType::MACROPORE:
        resistances_ = Only(false, true, false, false, false);
        break;
      case MassTransferCoefficientType::MACROPORE_KNUDSEN:
        resistances_ = Only(false, true, true, false, false);
        break;
      case MassTransferCoefficientType::HAYNES_1973:
        resistances_ = Only(true, true, false, false, false);
        film_resistance_factor_ = 1.0;
        molecular_diffusion_coefficient_ = 0.0;
        break;
      case MassTransferCoefficientType::MACROPORE_CAPACITY:
        resistances_ = Only(false, true, false, false, false);
        capacity_mode_ = static_cast<int>(CapacityFactorMode::SECANT);
        break;
      case MassTransferCoefficientType::SERIES_RESISTANCE:
        break;
    }
  }

  MassTransferCoefficientType MassTransferCoefficientSeriesResistance::PresetType() const {
    using T = MassTransferCoefficientType;
    const MassTransferResistances r = Resistances();
    const CapacityFactorMode mode = CapacityMode();
    if (r == Only(false, false, false, false, true)) return ActivationEnergyIsZero() ? T::CONSTANT : T::ARRHENIUS;
    if (mode == CapacityFactorMode::NONE) {
      if (r == Only(false, true, false, false, false)) return T::MACROPORE;
      if (r == Only(false, true, true, false, false)) return T::MACROPORE_KNUDSEN;
      if (r == Only(true, true, false, false, false) && IsNumber(film_resistance_factor_, 1.0) &&
          IsNumber(molecular_diffusion_coefficient_, 0.0))
        return T::HAYNES_1973;
    }
    if (mode == CapacityFactorMode::SECANT && r == Only(false, true, false, false, false)) return T::MACROPORE_CAPACITY;
    return T::SERIES_RESISTANCE;
  }

  std::string MassTransferCoefficientSeriesResistance::PresetTag() const {
    const auto t = PresetType();
    if (t == MassTransferCoefficientType::SERIES_RESISTANCE) return "custom";
    const KindDescriptor* k = KindAt(ModelFamily::RateCoefficient, static_cast<std::size_t>(t));
    return k ? k->tag : std::string("custom");
  }

  unsigned MassTransferCoefficientSeriesResistance::ReadsTransportProperties() const {
    const MassTransferResistances& r = resistances_;
    unsigned reads = kNoTransportProperties;
    if (MolecularDiffusion()) reads |= kMolecularDiffusion;
    if (r.film) reads |= kMolecularDiffusion | kSherwoodNumber;
    return reads;
  }

  void MassTransferCoefficientSeriesResistance::EvaluateExpressions(const EvalContext& ctx) {
    film_resistance_factor_.Evaluate(ctx); tortuosity_.Evaluate(ctx); macropore_void_fraction_.Evaluate(ctx);
    molecular_diffusion_coefficient_.Evaluate(ctx); crystal_radius_.Evaluate(ctx);
    micropore_diffusivity_.Evaluate(ctx); micropore_activation_energy_.Evaluate(ctx);
    capacity_factor_.Evaluate(ctx); concentration_floor_.Evaluate(ctx);
    k_.Evaluate(ctx); activation_energy_.Evaluate(ctx); pore_radius_.Evaluate(ctx);
  }

  void MassTransferCoefficientSeriesResistance::UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) {
    // The form is chosen by the switches and the mode only (never by a stored
    // preset name); each single-term form is the old class's arithmetic.
    const MassTransferResistances r = Resistances();
    const CapacityFactorMode mode = CapacityMode();
    if (r == Only(false, false, false, false, true)) {
      UpdateLumped(state, ads_id, component, layer, registry);
    } else if (mode == CapacityFactorMode::NONE && r == Only(false, true, false, false, false)) {
      UpdateMacropore(state, ads_id, fluid, component, layer, false);
    } else if (mode == CapacityFactorMode::NONE && r == Only(false, true, true, false, false)) {
      UpdateMacropore(state, ads_id, fluid, component, layer, true);
    } else if (mode == CapacityFactorMode::SECANT && r == Only(false, true, false, false, false)) {
      UpdateMacroporeSecant(state, ads_id, fluid, component, layer, registry);
    } else {
      UpdateSeries(state, ads_id, fluid, component, layer, registry);
    }
  }

  // The lumped term alone: the former constant (E_a = 0) and arrhenius kinds.
  void MassTransferCoefficientSeriesResistance::UpdateLumped(Equations& state, AdsorbentID ads_id, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) {
    auto& mtc = state.mass_transfer_coefficient.at(ads_id).at(component).var;
    auto& T = state.fluid_temperature.var;
    auto& P = state.pressure.var;
    auto& y = state.mole_fraction.at(component).var;

    if (ActivationEnergyIsZero()) {
      if (k_.IsRunConstant()) {
        // State-independent: one evaluation for the layer (same value every cell).
        EvalContext ctx{registry, T[layer.domain_.start_index], P[layer.domain_.start_index], 0.0, y[layer.domain_.start_index], 0.0};
        const double k = k_.EvaluateAt(ctx);
        for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) mtc[n] = k;
        return;
      }
      for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
        EvalContext ctx{registry, T[n], P[n], 0.0, y[n], 0.0};
        mtc[n] = k_.EvaluateAt(ctx);
      }
      return;
    }
    double R = constants::kIdealGasConstant;
    for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
      EvalContext ctx{registry, T[n], P[n], 0.0, y[n], 0.0};
      const double k0 = k_.EvaluateAt(ctx);
      const double Ea = activation_energy_.EvaluateAt(ctx);
      mtc[n] = k0 * std::exp(-Ea / (R * T[n]));
    }
  }

  // The macropore term alone, mode none: the former macropore and
  // macropore_knudsen kinds, k = c D_eff / l^2 (Glueckauf 1955: sphere 60 / d_p^2,
  // slab 3 / l^2).
  void MassTransferCoefficientSeriesResistance::UpdateMacropore(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, bool knudsen) {
    auto& mtc = state.mass_transfer_coefficient.at(ads_id).at(component).var;
    auto& T = state.fluid_temperature.var;
    auto& Dmi = state.molecular_diffusion.at(component).var;
    double ep = macropore_void_fraction_;
    const auto geo = layer.GetDiffusionGeometry();
    const double c_ldf = geo.c;
    double dp = geo.length;
    double tortuosity = tortuosity_;
    const double Dm_const = molecular_diffusion_coefficient_;  // > 0: constant D_m in place of the state's

    if (knudsen) {
      double R = constants::kIdealGasConstant;
      double Mr = fluid.molar_mass_.at(component);
      for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
        // Newton iterates can probe T < 0; skip just this cell rather than
        // abandoning the sweep and leaving the remaining cells stale.
        if (T[n] < 0)
          continue;
        double Dki = (2.0 / 3.0) * pore_radius_ * sqrt(8.0 * R * T[n] / (constants::kPI * Mr));
        double Dki_eff = Dki * ep / tortuosity;
        double Dmi_eff = (Dm_const > 0.0 ? Dm_const : Dmi[n]) * ep / tortuosity;
        double Di_eff = 1.0 / ((1.0 / Dki_eff) + (1.0 / Dmi_eff));
        mtc[n] = c_ldf * Di_eff / (dp * dp);
      }
      return;
    }
    for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
      if (T[n] < 0)
        continue;
      double Dmi_eff = (Dm_const > 0.0 ? Dm_const : Dmi[n]) * ep / tortuosity;
      mtc[n] = c_ldf * Dmi_eff / (dp * dp);
    }
  }

  // The macropore term alone, secant mode: the former macropore_capacity kind
  // (Haghpanah et al. 2013; Estupinan Perez et al. 2019, SI Eq. S4),
  //   k_i = (15 eps_p D_m,i / (tau r_p^2)) * min(c_i / (rho_p q_i*), 1),
  // with c and q* taken at max(c_i, c_floor) (q* re-evaluated there, every other
  // component at its actual value; the IAST mixture under the IAST rule).
  void MassTransferCoefficientSeriesResistance::UpdateMacroporeSecant(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) {
    auto& mtc = state.mass_transfer_coefficient.at(ads_id).at(component).var;
    const auto& T = state.fluid_temperature.var;
    const auto& P = state.pressure.var;
    const auto& Dmi = state.molecular_diffusion.at(component).var;
    auto& qi_sat = state.eq_ads_concentration.at(ads_id).at(component).var;
    auto& yi = state.mole_fraction.at(component).var;
    const double R = constants::kIdealGasConstant;

    const double ep = macropore_void_fraction_;
    const auto geo = layer.GetDiffusionGeometry();
    const double c_ldf = geo.c;
    const double dp = geo.length;
    const double tortuosity = tortuosity_;
    const double c_floor = std::max(concentration_floor_.Get(), 0.0);
    const double rho_p = solid_density_;
    const double Dm_const = molecular_diffusion_coefficient_;
    const int start = layer.domain_.start_index, end = layer.domain_.end_index;

    // q* at max(c_i, c_floor) from the isotherm itself, consistent with the
    // current y, P and T; the mole fraction is floored in place and restored,
    // and so is q*.
    std::vector<double>& q_floor = q_base_;
    const bool evaluate = isotherm_ && rho_p > 0.0;
    if (evaluate) {
      y_scratch_.assign(yi.begin(), yi.end());
      q_floor.assign(qi_sat.begin(), qi_sat.end());
      for (int n = start; n <= end; ++n) {
        if (!(T[n] > 0.0) || !(P[n] > 0.0)) continue;
        const double y_floor = c_floor * R * T[n] / P[n];
        if (yi[n] < y_floor) yi[n] = y_floor;
      }
      if (source_ && source_->mixing_rule_ == MixingRule::IAST) {
        auto& q_map = state.eq_ads_concentration.at(ads_id);
        const auto& components = fluid.components_;
        q_all_scratch_.resize(components.size());
        for (size_t i = 0; i < components.size(); ++i) {
          const auto& q = q_map.at(components[i]).var;
          q_all_scratch_[i].assign(q.begin(), q.end());
        }
        source_->UpdateIsothermAll(state, ads_id, components, start, end, registry);
        std::swap(q_floor, qi_sat);
        for (size_t i = 0; i < components.size(); ++i) {
          auto& q = q_map.at(components[i]).var;
          if (&q == &qi_sat) continue;  // restored by the swap above
          std::copy(q_all_scratch_[i].begin(), q_all_scratch_[i].end(), q.begin());
        }
      } else {
        const auto q_it = state.ads_concentration.find(ads_id);
        IsothermUpdateParams params(state.mole_fraction, state.eq_ads_concentration.at(ads_id),
                                    state.solid_temperature.var, state.pressure.var, component, start, end, registry,
                                    q_it == state.ads_concentration.end() ? nullptr : &q_it->second);
        isotherm_->Update(params);
        std::swap(q_floor, qi_sat);   // q_floor = q* at the floored state; qi_sat = the original, restored
      }
      std::copy(y_scratch_.begin(), y_scratch_.end(), yi.begin());
    }

    for (int n = start; n <= end; ++n) {
      if (T[n] < 0) continue;  // Newton iterates can probe T < 0; skip the cell
      const double Dmi_eff = (Dm_const > 0.0 ? Dm_const : Dmi[n]) * ep / tortuosity;
      const double k_macro = c_ldf * Dmi_eff / (dp * dp);
      double ratio = 1.0;
      if (evaluate && T[n] > 0.0 && P[n] > 0.0) {
        const double c = std::max(y_scratch_[n], 0.0) * P[n] / (R * T[n]);
        const double c_eff = std::max(c, c_floor);
        const double denom = rho_p * q_floor[n];
        ratio = denom > 0.0 ? std::min(c_eff / denom, 1.0) : 1.0;
      }
      mtc[n] = k_macro * ratio;
    }
  }

  void MassTransferCoefficientSeriesResistance::EvaluateLoading(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, int start, int end, const UserVariableRegistry* registry, std::vector<double>& out) {
    auto& qi_sat = state.eq_ads_concentration.at(ads_id).at(component).var;
    out.assign(qi_sat.begin(), qi_sat.end());
    if (source_ && source_->mixing_rule_ == MixingRule::IAST) {
      // IAST couples the components: solve the whole adsorbent, keep this
      // component's loading, restore every other one.
      auto& q_map = state.eq_ads_concentration.at(ads_id);
      const auto& components = fluid.components_;
      q_all_scratch_.resize(components.size());
      for (size_t i = 0; i < components.size(); ++i) {
        const auto& q = q_map.at(components[i]).var;
        q_all_scratch_[i].assign(q.begin(), q.end());
      }
      source_->UpdateIsothermAll(state, ads_id, components, start, end, registry);
      std::swap(out, qi_sat);  // out = the new q*; qi_sat = the original, restored
      for (size_t i = 0; i < components.size(); ++i) {
        auto& q = q_map.at(components[i]).var;
        if (&q == &qi_sat) continue;
        std::copy(q_all_scratch_[i].begin(), q_all_scratch_[i].end(), q.begin());
      }
      return;
    }
    const auto q_it = state.ads_concentration.find(ads_id);
    IsothermUpdateParams params(state.mole_fraction, state.eq_ads_concentration.at(ads_id),
                                state.solid_temperature.var, state.pressure.var, component, start, end, registry,
                                q_it == state.ads_concentration.end() ? nullptr : &q_it->second);
    isotherm_->Update(params);
    std::swap(out, qi_sat);  // out = the new q*; qi_sat = the original, restored
  }

  void MassTransferCoefficientSeriesResistance::CapacityFactor(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, int start, int end, const UserVariableRegistry* registry, std::vector<double>& k_out) {
    auto& yi = state.mole_fraction.at(component).var;
    const auto& T = state.fluid_temperature.var;
    const auto& P = state.pressure.var;
    const double R = constants::kIdealGasConstant;
    k_out.assign(yi.size(), 1.0);

    const CapacityFactorMode mode = CapacityMode();
    if (mode == CapacityFactorMode::NONE) return;
    if (mode == CapacityFactorMode::FIXED) {
      std::fill(k_out.begin() + start, k_out.begin() + end + 1, capacity_factor_.Get());
      return;
    }
    const double rho_p = solid_density_;
    if (!isotherm_ || !(rho_p > 0.0)) return;  // no isotherm context: K = 1

    // q* at c_eff = max(c, c_floor): this component's mole fraction floored
    // (every other component at its actual value); for the slope also at
    // c_eff (1 + h), a forward difference. Both vectors restored.
    constexpr double h = 1e-6;
    const bool slope = mode == CapacityFactorMode::SLOPE;
    const double c_floor = std::max(concentration_floor_.Get(), 0.0);
    auto valid = [&](int n) { return T[n] > 0.0 && P[n] > 0.0; };
    auto floored = [&](int n) {
      const double y_floor = c_floor * R * T[n] / P[n];
      return yi[n] < y_floor ? y_floor : yi[n];
    };
    y_scratch_.assign(yi.begin(), yi.end());
    for (int n = start; n <= end; ++n)
      if (valid(n)) yi[n] = floored(n);
    EvaluateLoading(state, ads_id, fluid, component, start, end, registry, q_base_);
    if (slope) {
      for (int n = start; n <= end; ++n)
        if (valid(n)) yi[n] = yi[n] * (1.0 + h);
      EvaluateLoading(state, ads_id, fluid, component, start, end, registry, q_pert_);
    }
    std::copy(y_scratch_.begin(), y_scratch_.end(), yi.begin());

    for (int n = start; n <= end; ++n) {
      if (!valid(n)) continue;
      const double y_eff = floored(n);
      const double to_c = P[n] / (R * T[n]);
      const double c_eff = y_eff * to_c;
      double K = 0.0;
      if (slope) {
        const double y_pert = y_eff * (1.0 + h);
        const double dc = (y_pert - y_eff) * to_c;
        K = dc > 0.0 ? rho_p * (q_pert_[n] - q_base_[n]) / dc : 0.0;    // slope
      }
      if (!(K > 0.0) || !std::isfinite(K)) K = c_eff > 0.0 ? rho_p * q_base_[n] / c_eff : 0.0;  // secant
      if (!(K >= 1.0) || !std::isfinite(K)) K = 1.0;  // floor: never above the K = 1 coefficient
      k_out[n] = K;
    }
  }

  // The series sum, every other combination (the former haynes and
  // series_resistance kinds and every new one).
  void MassTransferCoefficientSeriesResistance::UpdateSeries(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) {
    auto& mtc = state.mass_transfer_coefficient.at(ads_id).at(component).var;
    const auto& T = state.fluid_temperature.var;
    const auto& P = state.pressure.var;
    const auto& y = state.mole_fraction.at(component).var;
    const auto& Dmi = state.molecular_diffusion.at(component).var;
    const auto& Sh = state.sherwood_number.var;
    const double R = constants::kIdealGasConstant;
    const int start = layer.domain_.start_index, end = layer.domain_.end_index;
    const MassTransferResistances r = Resistances();

    if (r.film || r.macropore) CapacityFactor(state, ads_id, fluid, component, start, end, registry, capacity_);
    else capacity_.assign(mtc.size(), 1.0);

    const double f_film = film_resistance_factor_;
    const double tau = tortuosity_;
    const double ep = macropore_void_fraction_;
    const double Dm_const = molecular_diffusion_coefficient_;  // > 0: constant D_m in the macropore term
    const double rc = crystal_radius_;
    const double Dc0 = micropore_diffusivity_;
    const double Ea = micropore_activation_energy_;
    const bool molecular = MolecularDiffusion();
    const bool knudsen = KnudsenDiffusion();
    const double r_pore = pore_radius_;
    const double Mr = knudsen ? fluid.molar_mass_.at(component) : 0.0;
    const bool lumped_ea = r.lumped && !ActivationEnergyIsZero();
    // Terms that divide by T: Newton iterates can probe T <= 0; skip the cell.
    const bool needs_t = r.micropore || knudsen || lumped_ea;

    // The terms are written as in MassTransferCoefficient_Haynes so that with
    // f_film = 1, K = 1 and only film and macropore on the result is the same
    // double; a switched-off term adds an exact 0.0.
    const bool structured = layer.IsStructured();
    // A packed bed of infinite cylinders takes the geometry's form (V/S = r_f / 2,
    // l = r_f, c' = 8) with the packed k_f = Sh D / d_p; spheres keep their own
    // expressions below.
    const bool cylinder = !structured && layer.particle_shape_ == ParticleShape::Cylinder;
    const auto geo = layer.GetDiffusionGeometry();
    const double l = geo.length;
    const double c_geo = geo.c;
    const double dh = structured ? layer.GetHydraulicDiameter() : 0.0;
    const double dp = layer.particle_diameter_;
    const double rp = dp / 2.0;
    for (int n = start; n <= end; ++n) {
      if (needs_t && !(T[n] > 0.0)) continue;
      const double K = capacity_[n];
      double Dp = molecular ? (Dm_const > 0.0 ? Dm_const : Dmi[n]) : 0.0;
      if (knudsen) {
        // Knudsen alone: D = D_K; with molecular diffusion the Bosanquet
        // mean 1/D = 1/D_K + 1/D_m.
        const double Dk = (2.0 / 3.0) * r_pore * std::sqrt(8.0 * R * T[n] / (constants::kPI * Mr));
        Dp = molecular ? 1.0 / (1.0 / Dk + 1.0 / Dp) : Dk;
      }
      double film = 0.0, macro = 0.0, micro = 0.0, lumped = 0.0;
      if (structured) {
        if (r.film) {
          const double kf = Sh[n] * Dmi[n] / dh;
          film = geo.volume_to_surface / kf;
        }
        if (r.macropore) macro = tau * l * l / (3.0 * ep * Dp);
      } else if (cylinder) {
        if (r.film) {
          const double kf = Sh[n] * Dmi[n] / dp;
          film = geo.volume_to_surface / kf;
        }
        if (r.macropore) macro = tau * l * l / (c_geo * ep * Dp);
      } else {
        if (r.film) {
          const double kf = Sh[n] * Dmi[n] / dp;
          film = rp / (3.0 * kf);
        }
        if (r.macropore) macro = tau * rp * rp / (15.0 * ep * Dp);
      }
      if (r.micropore) micro = rc * rc / (15.0 * Dc0 * std::exp(-Ea / (R * T[n])));
      if (r.lumped) {
        EvalContext ctx{registry, T[n], P[n], 0.0, y[n], 0.0};
        double kL = k_.EvaluateAt(ctx);
        if (lumped_ea) kL = kL * std::exp(-activation_energy_.EvaluateAt(ctx) / (R * T[n]));
        lumped = 1.0 / kL;
      }
      mtc[n] = 1.0 / ((f_film * K) * film + K * macro + micro + lumped);
    }
  }

  // ---- conversion and factories ---------------------------------------------------------

  std::unique_ptr<IMassTransferCoefficientModel> ToResistances(std::unique_ptr<IMassTransferCoefficientModel> model) {
    if (!model || dynamic_cast<Generic*>(model.get())) return model;
    using T = MassTransferCoefficientType;
    auto g = std::make_unique<Generic>();
    if (auto* c = dynamic_cast<MassTransferCoefficientConstant*>(model.get())) {
      g->ApplyPreset(T::CONSTANT);
      g->k_ = c->constant_mtc_;
    } else if (auto* c = dynamic_cast<MassTransferCoefficientArrhenius*>(model.get())) {
      g->ApplyPreset(T::ARRHENIUS);
      g->k_ = c->pre_exp_factor_;
      g->activation_energy_ = c->activation_energy_;
    } else if (auto* c = dynamic_cast<MassTransferCoefficientMacropore*>(model.get())) {
      g->ApplyPreset(T::MACROPORE);
      g->tortuosity_ = c->tortuosity_;
      g->macropore_void_fraction_ = c->macropore_void_fraction_;
      g->molecular_diffusion_coefficient_ = c->molecular_diffusion_coefficient_;
    } else if (auto* c = dynamic_cast<MassTransferCoefficientMacroporeKnudsen*>(model.get())) {
      g->ApplyPreset(T::MACROPORE_KNUDSEN);
      g->pore_radius_ = c->pore_radius_;
      g->tortuosity_ = c->tortuosity_;
      g->macropore_void_fraction_ = c->macropore_void_fraction_;
      g->molecular_diffusion_coefficient_ = c->molecular_diffusion_coefficient_;
    } else if (auto* c = dynamic_cast<MassTransferCoefficient_Haynes*>(model.get())) {
      g->ApplyPreset(T::HAYNES_1973);
      g->tortuosity_ = c->tortuosity_;
      g->macropore_void_fraction_ = c->macropore_void_fraction_;
    } else if (auto* c = dynamic_cast<MassTransferCoefficientMacroporeCapacity*>(model.get())) {
      g->ApplyPreset(T::MACROPORE_CAPACITY);
      g->tortuosity_ = c->tortuosity_;
      g->macropore_void_fraction_ = c->macropore_void_fraction_;
      g->concentration_floor_ = c->concentration_floor_;
      g->molecular_diffusion_coefficient_ = c->molecular_diffusion_coefficient_;
    } else if (auto* c = dynamic_cast<LegacySeriesResistanceV1*>(model.get())) {
      g->film_resistance_factor_ = c->film_resistance_factor_;
      g->tortuosity_ = c->tortuosity_;
      g->macropore_void_fraction_ = c->macropore_void_fraction_;
      g->molecular_diffusion_coefficient_ = c->molecular_diffusion_coefficient_;
      g->crystal_radius_ = c->crystal_radius_;
      g->micropore_diffusivity_ = c->micropore_diffusivity_;
      g->micropore_activation_energy_ = c->micropore_activation_energy_;
      g->capacity_factor_ = c->capacity_factor_;
      g->concentration_floor_ = c->concentration_floor_;
      g->AdoptVersion1Switches();
    } else {
      return model;  // not a coefficient class this function knows: left as it is
    }
    return g;
  }

  std::unique_ptr<IMassTransferCoefficientModel> CreateLegacyMassTransferCoefficientModel(MassTransferCoefficientType type) {
    switch (type) {
    case MassTransferCoefficientType::CONSTANT: return std::make_unique<MassTransferCoefficientConstant>();
    case MassTransferCoefficientType::MACROPORE: return std::make_unique<MassTransferCoefficientMacropore>();
    case MassTransferCoefficientType::MACROPORE_KNUDSEN: return std::make_unique<MassTransferCoefficientMacroporeKnudsen>();
    case MassTransferCoefficientType::HAYNES_1973: return std::make_unique<MassTransferCoefficient_Haynes>();
    case MassTransferCoefficientType::ARRHENIUS: return std::make_unique<MassTransferCoefficientArrhenius>();
    case MassTransferCoefficientType::MACROPORE_CAPACITY: return std::make_unique<MassTransferCoefficientMacroporeCapacity>();
    case MassTransferCoefficientType::SERIES_RESISTANCE: return std::make_unique<LegacySeriesResistanceV1>();
    default: return nullptr;
    }
  }

  std::unique_ptr<IMassTransferCoefficientModel> CreateMassTransferCoefficientModel(MassTransferCoefficientType type) {
    // The tag of the generic model: its defaults (the lumped term only).
    if (type == MassTransferCoefficientType::SERIES_RESISTANCE) return std::make_unique<Generic>();
    return ToResistances(CreateLegacyMassTransferCoefficientModel(type));
  }

  namespace {

    void CheckSite(const MassTransferCoefficientSeriesResistance& g, bool ldf2, const std::string& where,
                   std::vector<std::string>& errors, std::vector<std::string>& warnings) {
      const MassTransferResistances r = g.Resistances();
      auto error = [&](const std::string& what) { errors.push_back("Error: " + where + ": " + what + "."); };
      auto positive = [&](const ExpressionParam& p, const char* name, const char* term) {
        if (!(p.Get() > 0.0))
          error(std::string(name) + " must be positive with the " + term + " resistance on (it is " +
                ExpressionParam::FormatScalar(p.Get()) + ")");
      };
      if (!(r.film || r.macropore || r.micropore || r.lumped)) error("no mass-transfer resistance is switched on");
      if (r.film) {
        const double f = g.film_resistance_factor_.Get();
        if (f < 0.0 || std::isnan(f)) error("the film resistance factor must not be negative");
        else if (f == 0.0)
          warnings.push_back(where + ": the film resistance is on with a film resistance factor of 0, which removes"
                             " it; switch the film resistance off instead.");
      }
      if (r.macropore) {
        positive(g.tortuosity_, "the tortuosity", "macropore");
        const double ep = g.macropore_void_fraction_.Get();
        if (!(ep > 0.0 && ep <= 1.0))
          error("the macropore void fraction must be in (0, 1] with the macropore resistance on (it is " +
                ExpressionParam::FormatScalar(ep) + ")");
        if (r.knudsen) positive(g.pore_radius_, "the pore radius", "Knudsen");
      }
      if (r.micropore) {
        positive(g.crystal_radius_, "the crystal radius", "micropore");
        positive(g.micropore_diffusivity_, "the micropore diffusivity", "micropore");
      }
      if (r.lumped && (g.k_.Get() < 0.0 || std::isnan(g.k_.Get())))
        error("the lumped mass transfer coefficient must not be negative");
      if (g.UsesCapacityFactor() && g.CapacityMode() == CapacityFactorMode::FIXED && !(g.capacity_factor_.Get() > 0.0))
        error("the fixed capacity factor must be positive");
      if (ldf2 && (r.film || r.macropore || r.micropore))
        warnings.push_back(where + ": a second-order LDF with a physical resistance on: the resistances give a"
                           " coefficient in 1/s, but the second-order rate constant is in kg/(mol s).");
    }

  }  // namespace

  void CheckMassTransferModel(const IMassTransferModel& model, const std::string& where,
                              std::vector<std::string>& errors, std::vector<std::string>& warnings) {
    const bool ldf2 = dynamic_cast<const LinearDrivingForceSecondOrder*>(&model) != nullptr;
    if (auto* g = dynamic_cast<const Generic*>(model.mtc_model_.get())) CheckSite(*g, ldf2, where, errors, warnings);
    const auto* dual = dynamic_cast<const DualRateLinearDrivingForce*>(&model);
    if (!dual || !dual->HasSlowSite()) return;
    if (auto* g = dynamic_cast<const Generic*>(dual->slow_mtc_model_.get()))
      CheckSite(*g, ldf2, where + " (slow site)", errors, warnings);
  }

  LinearDrivingForce::LinearDrivingForce() {
    mtc_model_ = CreateMassTransferCoefficientModel(MassTransferCoefficientType::CONSTANT);
  }

  LinearDrivingForceSecondOrder::LinearDrivingForceSecondOrder() {
    mtc_model_ = CreateMassTransferCoefficientModel(MassTransferCoefficientType::CONSTANT);
  }

  DualRateLinearDrivingForce::DualRateLinearDrivingForce() {
    mtc_model_ = CreateMassTransferCoefficientModel(MassTransferCoefficientType::CONSTANT);
    slow_mtc_model_ = CreateMassTransferCoefficientModel(MassTransferCoefficientType::CONSTANT);
  }

  bool DualRateLinearDrivingForce::SyncSlowSite() {
    if (!mtc_model_ || slow_mtc_model_) return false;
    if (auto* g = dynamic_cast<Generic*>(mtc_model_.get())) {
      slow_mtc_model_ = std::make_unique<Generic>(*g);
      return true;
    }
    // An old class on the fast site (only built by hand): a generic copy of it.
    const auto kind = mtc_model_->GetModelType();
    const KindDescriptor* to = KindAt(ModelFamily::RateCoefficient, static_cast<std::size_t>(kind));
    auto fresh = CreateMassTransferCoefficientModel(kind);
    if (!fresh || !to) return false;
    CarryOverParameters(*to, ModelFields(mtc_model_.get()), *to, ModelFields(fresh.get()));
    slow_mtc_model_ = std::move(fresh);
    return true;
  }

} // end adsim namespace


// Register the derived type with the base type
CEREAL_REGISTER_TYPE(adsim::MassTransferCoefficientConstant)
CEREAL_REGISTER_TYPE(adsim::MassTransferCoefficientMacroporeKnudsen)
CEREAL_REGISTER_TYPE(adsim::MassTransferCoefficientMacropore)
CEREAL_REGISTER_TYPE(adsim::MassTransferCoefficient_Haynes)
CEREAL_REGISTER_TYPE(adsim::MassTransferCoefficientMacroporeCapacity)
CEREAL_REGISTER_TYPE(adsim::MassTransferCoefficientArrhenius)
CEREAL_REGISTER_TYPE(adsim::MassTransferCoefficientSeriesResistance)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IMassTransferCoefficientModel, adsim::MassTransferCoefficientConstant)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IMassTransferCoefficientModel, adsim::MassTransferCoefficientMacroporeKnudsen)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IMassTransferCoefficientModel, adsim::MassTransferCoefficientMacropore)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IMassTransferCoefficientModel, adsim::MassTransferCoefficient_Haynes)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IMassTransferCoefficientModel, adsim::MassTransferCoefficientMacroporeCapacity)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IMassTransferCoefficientModel, adsim::MassTransferCoefficientArrhenius)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IMassTransferCoefficientModel, adsim::MassTransferCoefficientSeriesResistance)

CEREAL_REGISTER_TYPE(adsim::LinearDrivingForce)
CEREAL_REGISTER_TYPE(adsim::LinearDrivingForceSecondOrder)
CEREAL_REGISTER_TYPE(adsim::DualRateLinearDrivingForce)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IMassTransferModel, adsim::LinearDrivingForce)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IMassTransferModel, adsim::LinearDrivingForceSecondOrder)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IMassTransferModel, adsim::DualRateLinearDrivingForce)
