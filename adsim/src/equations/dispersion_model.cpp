/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/objects/reactor.h"
#include "adsim/objects/porous_media.h"


namespace adsim {

  void DispersionModelConstant::Update(Equations& state, const MixedAdsorbentLayer& layer) {		
    auto& Dx = state.axial_dispersion_coefficient.var;
    for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
      Dx[n] = dispersion_coefficient_;
    }
	}

	void HeatDispersionModelConstant::Update(Equations& state, const MixedAdsorbentLayer& layer) {
		auto& Lam = state.axial_thermal_conductivity.var;
    auto start = layer.domain_.start_index;
    auto end = layer.domain_.end_index;

		for (int n = start; n <= end; ++n) {
			Lam[n] = dispersion_coefficient_;
		}
	}

	void DispersionModelRuthven::Update(Equations& state, const MixedAdsorbentLayer& layer) {
		auto& Dx = state.axial_dispersion_coefficient.var;
		auto& u = state.velocity.var;
    auto& eb = state.domain_vectors.eb;
    auto& dp = state.domain_vectors.dp;

    for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
      Dx[n] = 0.7 * Dm + (std::abs(u[n]) / eb[n]) * 0.5 * dp[n];
    }
	}

	// void DispersionModelWakaoAndFunazkri::Update(Equations& state, const std::string& component, const MixedAdsorbentLayer& layer) {
	// 	auto& Dxi = state.dispersion_coefficient.at(component).var;
	// 	auto& Dmi = state.molecular_diffusion.at(component).var;
	// 	auto& u = state.velocity.var;
	// 	auto& Re = state.reynolds_number.var;
	// 	auto& Sc = state.schmidt_number.var;
  //   auto& eb = state.domain_vectors.eb;
  //   auto& dp = state.domain_vectors.dp;
    
  //   for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
  //     Dxi[n] = (Dmi[n] / eb[n]) * (20.0 + 0.5 * Sc[n] * Re[n]);
  //   }
	// }

	// void DispersionModelWakaoAndFunazkri::Update(Equations& state, const MixedAdsorbentLayer& layer) {
	// 	auto& Dx = state.axial_dispersion_coefficient.var;
	// 	auto& u = state.velocity.var;
	// 	auto& Re = state.reynolds_number.var;
	// 	auto& Sc = state.schmidt_number.var;
  //   auto& eb = state.domain_vectors.eb;
  //   auto& dp = state.domain_vectors.dp;
    
  //   for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
  //     Dx[n] = 20.0 + 0.5 * Sc[n] * Re[n]; // NOT CORRECT, DONT USE
  //   }
	// }

	void HeatDispersionModelRuthven::Update(Equations& state, const MixedAdsorbentLayer& layer) {
		auto& Lam = state.axial_thermal_conductivity.var;
		auto& Cpg = state.heat_capacity.var;
		auto& rho = state.fluid_density.var;
		auto& u = state.velocity.var;

    auto& eb = state.domain_vectors.eb;
    auto& dp = state.domain_vectors.dp;

    auto start = layer.domain_.start_index;
    auto end = layer.domain_.end_index;

		for (int n = start; n <= end; ++n) {
      double Dx = (std::abs(u[n]) / eb[n]) * 0.5 * dp[n];
      Lam[n] = Dx * Cpg[n] * rho[n];
    }
	}

	void HeatDispersionModelYagiandKunii::Update(Equations& state, const MixedAdsorbentLayer& layer) {
		auto& Lam = state.axial_thermal_conductivity.var;
		auto& Cpg = state.heat_capacity.var;
		auto& kg = state.fluid_thermal_conductivity.var;
		auto& vis = state.viscosity.var;
		auto& re = state.reynolds_number.var;

    auto& eb = state.domain_vectors.eb;
    auto& ks = state.domain_vectors.ks;
    
    auto start = layer.domain_.start_index;
    auto end = layer.domain_.end_index;

		for (int n = start; n <= end; ++n) {
      double exponent = 0.28 - 0.757 * std::log10(eb[n]) - 0.057 * std::log10(ks[n] / kg[n]);
      double K_z0 = kg[n] * std::pow((ks[n] / kg[n]), exponent);
      double pr = Cpg[n] * vis[n] / kg[n]; // Prandlt number, TODO: move into state.

      Lam[n] = K_z0 + 0.75 * kg[n] * pr * re[n];
    }
  }

} // end adsim namespace


// Register the derived type with the base type
CEREAL_REGISTER_TYPE(adsim::DispersionModelConstant)
CEREAL_REGISTER_TYPE(adsim::HeatDispersionModelConstant)
// CEREAL_REGISTER_TYPE(adsim::DispersionModelWakaoAndFunazkri)
CEREAL_REGISTER_TYPE(adsim::DispersionModelRuthven)
CEREAL_REGISTER_TYPE(adsim::HeatDispersionModelRuthven)
CEREAL_REGISTER_TYPE(adsim::HeatDispersionModelYagiandKunii)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IDispersionModel, adsim::DispersionModelConstant)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IHeatDispersionModel, adsim::HeatDispersionModelConstant)
// CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IDispersionModel, adsim::DispersionModelWakaoAndFunazkri)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IDispersionModel, adsim::DispersionModelRuthven)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IHeatDispersionModel, adsim::HeatDispersionModelRuthven)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IHeatDispersionModel, adsim::HeatDispersionModelYagiandKunii)
CEREAL_CLASS_VERSION(adsim::DispersionModelConstant, 1);
CEREAL_CLASS_VERSION(adsim::HeatDispersionModelConstant, 1);
// CEREAL_CLASS_VERSION(adsim::DispersionModelWakaoAndFunazkri, 1);
CEREAL_CLASS_VERSION(adsim::DispersionModelRuthven, 1);
CEREAL_CLASS_VERSION(adsim::HeatDispersionModelRuthven, 1);
CEREAL_CLASS_VERSION(adsim::HeatDispersionModelYagiandKunii, 1);
