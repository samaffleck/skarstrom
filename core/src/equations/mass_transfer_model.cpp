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
    double dp = layer.particle_diameter_;
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
      mtc[n] = 60.0 * Di_eff / (dp * dp);
    }
  }

  void MassTransferCoefficientMacropore::UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) {
    auto& mtc = state.mass_transfer_coefficient.at(ads_id).at(component).var;
    auto& T = state.fluid_temperature.var;
    auto& Dmi = state.molecular_diffusion.at(component).var;
    double ep = macropore_void_fraction_;
    double dp = layer.particle_diameter_;
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
      mtc[n] = 60.0 * Dmi_eff / (dp * dp);
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
    const double dp = layer.particle_diameter_;
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
      IsothermUpdateParams params(state.mole_fraction, state.eq_ads_concentration.at(ads_id),
                                  state.solid_temperature.var, state.pressure.var, component, start, end, registry);
      isotherm_->Update(params);
      std::swap(q_scratch_, qi_sat);   // q_scratch_ = q* at the floored state; qi_sat = the original, restored
      std::copy(y_scratch_.begin(), y_scratch_.end(), yi.begin());
    }

    for (int n = start; n <= end; ++n) {
      if (T[n] < 0) continue;  // Newton iterates can probe T < 0; skip the cell
      const double Dmi_eff = (Dm_const > 0.0 ? Dm_const : Dmi[n]) * ep / tortuosity;
      const double k_macro = 60.0 * Dmi_eff / (dp * dp);
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

  LinearDrivingForce::LinearDrivingForce() {
    mtc_model_ = std::make_unique<MassTransferCoefficientConstant>();
  }

  LinearDrivingForceSecondOrder::LinearDrivingForceSecondOrder() {
    mtc_model_ = std::make_unique<MassTransferCoefficientConstant>();
  }

} // end adsim namespace


// Register the derived type with the base type
CEREAL_REGISTER_TYPE(adsim::MassTransferCoefficientConstant)
CEREAL_REGISTER_TYPE(adsim::MassTransferCoefficientMacroporeKnudsen)
CEREAL_REGISTER_TYPE(adsim::MassTransferCoefficientMacropore)
CEREAL_REGISTER_TYPE(adsim::MassTransferCoefficient_Haynes)
CEREAL_REGISTER_TYPE(adsim::MassTransferCoefficientMacroporeCapacity)
CEREAL_REGISTER_TYPE(adsim::MassTransferCoefficientArrhenius)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IMassTransferCoefficientModel, adsim::MassTransferCoefficientConstant)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IMassTransferCoefficientModel, adsim::MassTransferCoefficientMacroporeKnudsen)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IMassTransferCoefficientModel, adsim::MassTransferCoefficientMacropore)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IMassTransferCoefficientModel, adsim::MassTransferCoefficient_Haynes)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IMassTransferCoefficientModel, adsim::MassTransferCoefficientMacroporeCapacity)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IMassTransferCoefficientModel, adsim::MassTransferCoefficientArrhenius)
CEREAL_CLASS_VERSION(adsim::MassTransferCoefficientConstant, 3);
CEREAL_CLASS_VERSION(adsim::MassTransferCoefficientMacroporeKnudsen, 3);
CEREAL_CLASS_VERSION(adsim::MassTransferCoefficientMacropore, 3);
CEREAL_CLASS_VERSION(adsim::MassTransferCoefficient_Haynes, 2);
CEREAL_CLASS_VERSION(adsim::MassTransferCoefficientMacroporeCapacity, 2);
CEREAL_CLASS_VERSION(adsim::MassTransferCoefficientArrhenius, 3);

CEREAL_REGISTER_TYPE(adsim::LinearDrivingForce)
CEREAL_REGISTER_TYPE(adsim::LinearDrivingForceSecondOrder)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IMassTransferModel, adsim::LinearDrivingForce)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IMassTransferModel, adsim::LinearDrivingForceSecondOrder)
CEREAL_CLASS_VERSION(adsim::LinearDrivingForce, 1);
CEREAL_CLASS_VERSION(adsim::LinearDrivingForceSecondOrder, 1);
