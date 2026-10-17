/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/objects/reactor.h"
#include "adsim/equations/equation.h"


namespace adsim {

	void ConstantInternalHeatTransferModel::Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall) {
		auto& hi = state.internal_heat_transfer_coefficient.var;

    auto start = layer.domain_.start_index;
    auto end = layer.domain_.end_index;
		for (int n = start; n <= end; ++n) {
			hi[n] = heat_transfer_coefficient_;
		}
	}

	void InternalHeatTransferModel_Dixon::Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall) {
		auto& hi = state.internal_heat_transfer_coefficient.var;
		const auto& Re = state.reynolds_number.var;
		const auto& k_g = state.fluid_thermal_conductivity.var;
		const auto& u = state.velocity.var;
		const auto& rho_g = state.fluid_density.var;
		const auto& Cpg = state.heat_capacity.var;
    const auto& eb = state.domain_vectors.eb;
    const auto& dx = state.domain_vectors.cell_width;
    const auto& dp = state.domain_vectors.dp;
    const auto& k_s = state.domain_vectors.ks;
    const auto& D_c = state.domain_vectors.di;
    const int N = hi.size();

    auto start = layer.domain_.start_index;
    auto end = layer.domain_.end_index;
		for (int n = start; n <= end; ++n) {
      const double re = std::max(0.0, Re[n]); // Sometimes Re can go negative as density goes negative as mole fraction goes negative. 
      
      // heat transfer coefficient at the wall
      double h_g = (k_g[n] / dp[n]) * (0.0835 * std::pow(re, 0.91)); // Valid for 10 < Re < 1200
      double theta_w = 0.0024 * std::pow((D_c[n] / dp[n]), 1.58);
      double h_0 = (k_g[n] / dp[n]) * (2 * eb[n] + (1 - eb[n]) / ((k_g[n] / k_s[n]) * (1.0 / 3.0) + theta_w));
			double h_w = h_0 + h_g;
			
      // effective thermal conductivity through the bed
      double phi = 0.22 * eb[n] * eb[n];
      double k_e0 = k_g[n] * (eb[n] + ( (1 - eb[n]) / ( (k_g[n] / k_s[n]) * (2.0 / 3.0) + phi ) ));
      double d = dp[n] / D_c[n];
      double Pe_H = 8.65 * (1 + 19.4 * d * d);
      double vel = 0.5 * (u[n] + u[n - 1]);
      double k_eg = std::abs(vel) * rho_g[n] * Cpg[n] * dp[n] / Pe_H;
      double k_e = k_e0 + k_eg;

      // Dixon (1996)
      double Bi = h_w * D_c[n] / (2.0 * k_e);
      // assert(h_w > 0.0);
      double h_l_invr = (1.0 / h_w) + (D_c[n] / (6.0 * k_e) * ((Bi + 3.0) / (Bi + 4.0))); 
      // assert(h_l_invr > 0.0);
      hi[n] = 1.0 / h_l_invr;
		}
    
    // Ghost cells
    hi[0] = hi[1];
    hi[N - 1] = hi[N - 2];
	}

	void InternalHeatTransferModel_Leva::Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall) {
		auto& hi = state.internal_heat_transfer_coefficient.var;
		const auto& Re = state.reynolds_number.var;
		const auto& k_g = state.fluid_thermal_conductivity.var;
		const auto& dp = state.domain_vectors.dp;
		const auto& D_c = state.domain_vectors.di;
		const int N = hi.size();

		const double e1 = eta1, e2 = eta2, ec = exp_coeff;

		auto start = layer.domain_.start_index;
		auto end = layer.domain_.end_index;
		for (int n = start; n <= end; ++n) {
			const double re = std::max(0.0, Re[n]);
			// Leva (1947): Nu = h_L·D_t/k_g = eta1·Re^eta2·exp(−exp_coeff·dp/D_t)
			const double Nu = e1 * std::pow(re, e2) * std::exp(-ec * dp[n] / D_c[n]);
			hi[n] = Nu * k_g[n] / D_c[n];
		}

		// Ghost cells
		hi[0] = hi[1];
		hi[N - 1] = hi[N - 2];
	}

  void ConstantExternalHeatTransferModel::Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall) {
		auto& ho = state.external_heat_transfer_coefficient.var;

    auto start = layer.domain_.start_index;
    auto end = layer.domain_.end_index;
		for (int n = start; n <= end; ++n) {
      ho[n] = heat_transfer_coefficient_;
		}
	}

	void ExternalHeatTransferModel_Holman::Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall) {
		auto& ho = state.external_heat_transfer_coefficient.var;

		// TODO: Should depend on user input temperature
		const double k_air = 0.0259;         // W/(m-K)
		const double rho_air = 1.204;        // kg/m
		const double vis_air = 1.81e-5;      // Pa-s
		const double cp_air = 1005;          // J/(kg-K)
		const double Pr_air = cp_air * vis_air / k_air;
    const auto& de = state.domain_vectors.de;

    auto start = layer.domain_.start_index;
    auto end = layer.domain_.end_index;
    const double pr_exp = n;  // Prandtl exponent (member); avoid shadowing by the loop index
		for (int i = start; i <= end; ++i) {
      double Re_air = de[i] * external_velocity * rho_air / vis_air;
      double Nu_air = c * std::pow(Re_air, m) * std::pow(Pr_air, pr_exp);
      double htc = Nu_air * k_air / de[i];
			ho[i] = htc;
		}
	}

  void ConstantParticleHeatTransferModel::Update(Equations& state, const MixedAdsorbentLayer& layer) {
    auto& hp = state.particle_heat_transfer_coefficient.var;
    const auto& eb = state.domain_vectors.eb;
    const auto& dp = state.domain_vectors.dp;

    auto start = layer.domain_.start_index;
    auto end = layer.domain_.end_index;
    for (int n = start; n <= end; ++n) {
      double av = 6.0 * (1.0 - eb[n]) / dp[n]; // Specific surface area [1/m]
      hp[n] = heat_transfer_coefficient_ * av;   // Volumetric HTC [W/(m3.K)]
    }
  }

} // end adsim namespace


CEREAL_REGISTER_TYPE(adsim::ConstantParticleHeatTransferModel)
CEREAL_REGISTER_TYPE(adsim::ConstantInternalHeatTransferModel)
CEREAL_REGISTER_TYPE(adsim::InternalHeatTransferModel_Dixon)
CEREAL_REGISTER_TYPE(adsim::InternalHeatTransferModel_Leva)
CEREAL_REGISTER_TYPE(adsim::ConstantExternalHeatTransferModel)
CEREAL_REGISTER_TYPE(adsim::ExternalHeatTransferModel_Holman)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IInternalHeatTransferModel, adsim::ConstantInternalHeatTransferModel)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IInternalHeatTransferModel, adsim::InternalHeatTransferModel_Dixon)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IInternalHeatTransferModel, adsim::InternalHeatTransferModel_Leva)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IExternalHeatTransferModel, adsim::ConstantExternalHeatTransferModel)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IExternalHeatTransferModel, adsim::ExternalHeatTransferModel_Holman)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IParticleHeatTransferModel, adsim::ConstantParticleHeatTransferModel)
CEREAL_CLASS_VERSION(adsim::ConstantParticleHeatTransferModel, 1);
CEREAL_CLASS_VERSION(adsim::ConstantInternalHeatTransferModel, 1);
CEREAL_CLASS_VERSION(adsim::InternalHeatTransferModel_Dixon, 1);
CEREAL_CLASS_VERSION(adsim::InternalHeatTransferModel_Leva, 1);
CEREAL_CLASS_VERSION(adsim::ConstantExternalHeatTransferModel, 1);
CEREAL_CLASS_VERSION(adsim::ExternalHeatTransferModel_Holman, 1);
