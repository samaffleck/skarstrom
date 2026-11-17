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

	namespace {
		// sum_i y_i f(D_m,i) over the layer's cells, into out: the logged
		// mixture value of a per-component dispersion coefficient. Summed in
		// the map order of the mole fractions (deterministic).
		template <class F>
		void MoleFractionWeightedMean(Equations& state, const MixedAdsorbentLayer& layer, F f, std::vector<double>& out) {
			const int start = layer.domain_.start_index;
			const int end = layer.domain_.end_index;
			for (int n = start; n <= end; ++n) out[n] = 0.0;
			for (const auto& [component, yi_eq] : state.mole_fraction) {
				const auto& yi = yi_eq.var;
				const auto& Dmi = state.molecular_diffusion.at(component).var;
				for (int n = start; n <= end; ++n) out[n] += yi[n] * f(n, Dmi[n]);
			}
		}
	}

	// The logged mixture value sum_i y_i D_ax,i of Ruthven's
	// D_ax,i = 0.7 D_m,i + 0.5 d_p |u| / eb.
	void DispersionModelRuthven::Update(Equations& state, const MixedAdsorbentLayer& layer) {
		auto& Dx = state.axial_dispersion_coefficient.var;
		const auto& u = state.velocity.var;
		const auto& eb = state.domain_vectors.eb;
		const auto& dp = state.domain_vectors.dp;
		MoleFractionWeightedMean(state, layer,
			[&](int n, double dmi) { return 0.7 * dmi + (std::abs(u[n]) / eb[n]) * 0.5 * dp[n]; }, Dx);
	}

	void DispersionModelRuthven::UpdateComponent(Equations& state, const MixedAdsorbentLayer& layer,
	                                             const std::string& component, std::vector<double>& Dxi) {
		const auto& u = state.velocity.var;
		const auto& eb = state.domain_vectors.eb;
		const auto& dp = state.domain_vectors.dp;
		const auto& Dmi = state.molecular_diffusion.at(component).var;
		for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
			Dxi[n] = 0.7 * Dmi[n] + (std::abs(u[n]) / eb[n]) * 0.5 * dp[n];
		}
	}

	// D_ax = D_m + u_ch^2 D_h^2 / (kappa D_m) in cell n; D_m alone where the
	// cell is not a channel.
	double DispersionModelTaylorAris::TaylorAris(int n, double dm, const std::vector<double>& u,
	                                             const std::vector<double>& eb, const std::vector<double>& dh,
	                                             const std::vector<double>& kappa) {
		if (!(kappa[n] > 0.0) || !(dm > 0.0)) return dm;  // not a channel: molecular diffusion only
		const double u_ch = std::abs(u[n]) / eb[n];
		return dm + u_ch * u_ch * dh[n] * dh[n] / (kappa[n] * dm);
	}

	void DispersionModelTaylorAris::Update(Equations& state, const MixedAdsorbentLayer& layer) {
		auto& Dx = state.axial_dispersion_coefficient.var;
		const auto& u = state.velocity.var;
		const auto& eb = state.domain_vectors.eb;   // the open frontal area of a structured layer
		const auto& dh = state.domain_vectors.hydraulic_diameter;
		const auto& kappa = state.domain_vectors.taylor_aris_kappa;
		MoleFractionWeightedMean(state, layer, [&](int n, double dmi) { return TaylorAris(n, dmi, u, eb, dh, kappa); }, Dx);
	}

	void DispersionModelTaylorAris::UpdateComponent(Equations& state, const MixedAdsorbentLayer& layer,
	                                                const std::string& component, std::vector<double>& Dxi) {
		const auto& u = state.velocity.var;
		const auto& eb = state.domain_vectors.eb;
		const auto& dh = state.domain_vectors.hydraulic_diameter;
		const auto& kappa = state.domain_vectors.taylor_aris_kappa;
		const auto& Dmi = state.molecular_diffusion.at(component).var;
		for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
			Dxi[n] = TaylorAris(n, Dmi[n], u, eb, dh, kappa);
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

	void HeatDispersionModelChannel::Update(Equations& state, const MixedAdsorbentLayer& layer) {
		auto& Lam = state.axial_thermal_conductivity.var;
		const auto& kg = state.fluid_thermal_conductivity.var;
		const auto& dv = state.domain_vectors;
		const auto& ofa = dv.eb;                  // the open frontal area of a structured layer
		const auto& lambda_s = dv.solid_conductivity;  // (1 - OFA) k_s; 0 in a packed layer
		// LTE lumps the solid's conduction into the one energy balance; LTNE
		// gives it to the solid energy balance (UpdateSolidTemperature).
		const bool lte = dv.local_thermal_equilibrium;
		for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
			Lam[n] = ofa[n] * kg[n] + (lte ? lambda_s[n] : 0.0);
		}
	}

} // end adsim namespace


// Register the derived type with the base type
CEREAL_REGISTER_TYPE(adsim::DispersionModelConstant)
CEREAL_REGISTER_TYPE(adsim::HeatDispersionModelConstant)
// CEREAL_REGISTER_TYPE(adsim::DispersionModelWakaoAndFunazkri)
CEREAL_REGISTER_TYPE(adsim::DispersionModelRuthven)
CEREAL_REGISTER_TYPE(adsim::DispersionModelTaylorAris)
CEREAL_REGISTER_TYPE(adsim::HeatDispersionModelRuthven)
CEREAL_REGISTER_TYPE(adsim::HeatDispersionModelYagiandKunii)
CEREAL_REGISTER_TYPE(adsim::HeatDispersionModelChannel)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IDispersionModel, adsim::DispersionModelConstant)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IHeatDispersionModel, adsim::HeatDispersionModelConstant)
// CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IDispersionModel, adsim::DispersionModelWakaoAndFunazkri)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IDispersionModel, adsim::DispersionModelRuthven)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IDispersionModel, adsim::DispersionModelTaylorAris)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IHeatDispersionModel, adsim::HeatDispersionModelRuthven)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IHeatDispersionModel, adsim::HeatDispersionModelYagiandKunii)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IHeatDispersionModel, adsim::HeatDispersionModelChannel)
// CEREAL_CLASS_VERSION(adsim::DispersionModelWakaoAndFunazkri, 1);
