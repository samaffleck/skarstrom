/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/ambient_air.h"
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

	double LumpedWallHeatTransferCoefficient(double h_w, double k_e, double D_c) {
		if (h_w <= 0.0) return 0.0;  // no wall coefficient: no heat to the wall
		// Dixon (1996)
		double Bi = h_w * D_c / (2.0 * k_e);
		double h_l_invr = (1.0 / h_w) + (D_c / (6.0 * k_e) * ((Bi + 3.0) / (Bi + 4.0)));
		return 1.0 / h_l_invr;
	}

	double DixonWallCoefficient(double k_g, double k_s, double eb, double dp, double D_c, double re) {
		double h_g = (k_g / dp) * (0.0835 * std::pow(re, 0.91)); // Valid for 10 < Re < 1200
		double theta_w = 0.0024 * std::pow((D_c / dp), 1.58);
		double h_0 = (k_g / dp) * (2 * eb + (1 - eb) / ((k_g / k_s) * (1.0 / 3.0) + theta_w));
		return h_0 + h_g;
	}

	double DixonRadialConductivity(double k_g, double k_s, double eb, double dp, double D_c, double u,
	                               double rho_g, double cp_g) {
		double phi = 0.22 * eb * eb;
		double k_e0 = k_g * (eb + ((1 - eb) / ((k_g / k_s) * (2.0 / 3.0) + phi)));
		double d = dp / D_c;
		double Pe_H = 8.65 * (1 + 19.4 * d * d);
		double k_eg = std::abs(u) * rho_g * cp_g * dp / Pe_H;
		return k_e0 + k_eg;
	}

	double LiFinlaysonWallNusselt(double re, bool cylinders) {
		return cylinders ? 0.16 * std::pow(re, 0.93) : 0.17 * std::pow(re, 0.79);
	}

	double MartinNillesWallNusselt(double re, double pr, double k0_over_kg, double dp_over_D) {
		return (1.3 + 5.0 * dp_over_D) * k0_over_kg + 0.19 * std::pow(re, 0.75) * std::cbrt(pr);
	}

	double ZehnerSchlunderRadialConductivity(double k_g, double k0_over_kg, double re, double pr, double dp_over_D) {
		const double a = 1.0 - 2.0 * dp_over_D;
		const double K_r = 7.0 * (2.0 - a * a);
		return k_g * (k0_over_kg + re * pr / K_r);
	}

	void InternalHeatTransferModel_Dixon::Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall) {
		auto& hi = state.internal_heat_transfer_coefficient.var;
		const auto& Re = state.reynolds_number.var;
		const auto& k_g = state.fluid_thermal_conductivity.var;
		const auto& u = state.velocity.var;
		const auto& rho_g = state.fluid_density.var;
		const auto& Cpg = state.heat_capacity.var;
    const auto& eb = state.domain_vectors.eb;
    const auto& dp = state.domain_vectors.dp;
    const auto& k_s = state.domain_vectors.ks;
    const auto& D_c = state.domain_vectors.di;
    const int N = hi.size();

    auto start = layer.domain_.start_index;
    auto end = layer.domain_.end_index;
		for (int n = start; n <= end; ++n) {
      const double re = std::max(0.0, Re[n]); // Sometimes Re can go negative as density goes negative as mole fraction goes negative.
      const double h_w = DixonWallCoefficient(k_g[n], k_s[n], eb[n], dp[n], D_c[n], re);
      const double vel = 0.5 * (u[n] + u[n - 1]);
      const double k_e = DixonRadialConductivity(k_g[n], k_s[n], eb[n], dp[n], D_c[n], vel, rho_g[n], Cpg[n]);
      hi[n] = LumpedWallHeatTransferCoefficient(h_w, k_e, D_c[n]);
		}

    // Ghost cells
    hi[0] = hi[1];
    hi[N - 1] = hi[N - 2];
	}

	void InternalHeatTransferModel_LiFinlayson::Update(Equations& state, const MixedAdsorbentLayer& layer,
	                                                   const BedWall& wall) {
		auto& hi = state.internal_heat_transfer_coefficient.var;
		const auto& Re = state.reynolds_number.var;
		const auto& k_g = state.fluid_thermal_conductivity.var;
		const auto& u = state.velocity.var;
		const auto& rho_g = state.fluid_density.var;
		const auto& Cpg = state.heat_capacity.var;
		const auto& eb = state.domain_vectors.eb;
		const auto& dp = state.domain_vectors.dp;
		const auto& k_s = state.domain_vectors.ks;
		const auto& D_c = state.domain_vectors.di;
		const int N = hi.size();
		const bool cylinders = layer.particle_shape_ == ParticleShape::Cylinder;

		for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
			const double re = std::max(0.0, Re[n]);
			const double h_w = LiFinlaysonWallNusselt(re, cylinders) * k_g[n] / dp[n];
			const double vel = 0.5 * (u[n] + u[n - 1]);
			const double k_e = DixonRadialConductivity(k_g[n], k_s[n], eb[n], dp[n], D_c[n], vel, rho_g[n], Cpg[n]);
			hi[n] = LumpedWallHeatTransferCoefficient(h_w, k_e, D_c[n]);
		}

		// Ghost cells
		hi[0] = hi[1];
		hi[N - 1] = hi[N - 2];
	}

	void InternalHeatTransferModel_MartinNilles::Update(Equations& state, const MixedAdsorbentLayer& layer,
	                                                    const BedWall& wall) {
		auto& hi = state.internal_heat_transfer_coefficient.var;
		const auto& Re = state.reynolds_number.var;
		const auto& k_g = state.fluid_thermal_conductivity.var;
		const auto& Cpg = state.heat_capacity.var;
		const auto& vis = state.viscosity.var;
		const auto& eb = state.domain_vectors.eb;
		const auto& dp = state.domain_vectors.dp;
		const auto& k_s = state.domain_vectors.ks;
		const auto& D_c = state.domain_vectors.di;
		const int N = hi.size();

		for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
			const double re = std::max(0.0, Re[n]);
			const double pr = Cpg[n] * vis[n] / k_g[n];
			const double k0 = ZehnerSchlunderConductivityRatio(eb[n], k_s[n] / k_g[n]);
			const double d = dp[n] / D_c[n];
			const double h_w = MartinNillesWallNusselt(re, pr, k0, d) * k_g[n] / dp[n];
			const double k_e = ZehnerSchlunderRadialConductivity(k_g[n], k0, re, pr, d);
			hi[n] = LumpedWallHeatTransferCoefficient(h_w, k_e, D_c[n]);
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

  void ConstantExternalHeatTransferModel::Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall,
                                                 double ambient_temperature) {
    (void)wall;
    (void)ambient_temperature;
    auto& ho = state.external_heat_transfer_coefficient.var;

    auto start = layer.domain_.start_index;
    auto end = layer.domain_.end_index;
    for (int n = start; n <= end; ++n) {
      ho[n] = heat_transfer_coefficient_;
    }
  }

  namespace {
    // Film temperature of a cell's outer wall surface [K]. A non-physical wall
    // temperature (a trial Newton iterate) falls back to the ambient one so
    // the air properties stay defined.
    double FilmTemperature(double T_wall, double T_inf) {
      return T_wall > 0.0 ? 0.5 * (T_wall + T_inf) : T_inf;
    }
  }  // namespace

  void ExternalHeatTransferModel_Holman::Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall,
                                                double ambient_temperature) {
    (void)wall;
    auto& ho = state.external_heat_transfer_coefficient.var;
    const auto& Tw = state.wall_temperature.var;
    const auto& de = state.domain_vectors.de;

    auto start = layer.domain_.start_index;
    auto end = layer.domain_.end_index;
    const double pr_exp = n;  // Prandtl exponent (member); avoid shadowing by the loop index
    for (int i = start; i <= end; ++i) {
      // Air properties at the film temperature (ambient_air.h).
      const double T_film = FilmTemperature(Tw[i], ambient_temperature);
      const double k_air = ambient_air::ThermalConductivity(T_film);
      const double rho_air = ambient_air::Density(T_film);
      const double vis_air = ambient_air::Viscosity(T_film);
      const double Pr_air = ambient_air::Prandtl(T_film);
      const double Re_air = de[i] * external_velocity * rho_air / vis_air;
      const double Nu_air = c * std::pow(Re_air, m) * std::pow(Pr_air, pr_exp);
      ho[i] = Nu_air * k_air / de[i];
    }
  }

  void ExternalHeatTransferModel_ChurchillChu::Update(Equations& state, const MixedAdsorbentLayer& layer,
                                                      const BedWall& wall, double ambient_temperature) {
    (void)wall;
    auto& ho = state.external_heat_transfer_coefficient.var;
    const auto& Tw = state.wall_temperature.var;
    const auto& de = state.domain_vectors.de;
    const auto& dx = state.domain_vectors.cell_width;

    // The wall runs along every layer, so a vertical column is one plate of
    // the bed's height (all layers; the lumped void regions have no length).
    double height = 0.0;
    for (int i : state.domain_vectors.bed_cell_indices) height += dx[i];

    for (int i = layer.domain_.start_index; i <= layer.domain_.end_index; ++i) {
      const double T_wall = Tw[i] > 0.0 ? Tw[i] : ambient_temperature;
      ho[i] = ambient_air::ChurchillChuCoefficient(vertical_, T_wall, ambient_temperature, height, de[i]);
    }
  }

  void ExternalHeatTransferModel_YangCylinder::Update(Equations& state, const MixedAdsorbentLayer& layer,
                                                      const BedWall& wall, double ambient_temperature) {
    (void)wall;
    auto& ho = state.external_heat_transfer_coefficient.var;
    const auto& Tw = state.wall_temperature.var;
    const auto& de = state.domain_vectors.de;
    const auto& dx = state.domain_vectors.cell_width;

    // One cylinder of the bed's height (all layers), as churchill_chu_vertical.
    double height = 0.0;
    for (int i : state.domain_vectors.bed_cell_indices) height += dx[i];

    for (int i = layer.domain_.start_index; i <= layer.domain_.end_index; ++i) {
      const double T_wall = Tw[i] > 0.0 ? Tw[i] : ambient_temperature;
      ho[i] = ambient_air::YangCoefficient(T_wall, ambient_temperature, height, de[i]);
    }
  }

  void ExternalHeatTransferModel_ChurchillBernstein::Update(Equations& state, const MixedAdsorbentLayer& layer,
                                                            const BedWall& wall, double ambient_temperature) {
    (void)wall;
    auto& ho = state.external_heat_transfer_coefficient.var;
    const auto& Tw = state.wall_temperature.var;
    const auto& de = state.domain_vectors.de;
    const double v = external_velocity;
    for (int i = layer.domain_.start_index; i <= layer.domain_.end_index; ++i) {
      const double T_wall = Tw[i] > 0.0 ? Tw[i] : ambient_temperature;
      ho[i] = ambient_air::ChurchillBernsteinCoefficient(T_wall, ambient_temperature, v, de[i]);
    }
  }

  void ExternalHeatTransferModel_MixedConvection::Update(Equations& state, const MixedAdsorbentLayer& layer,
                                                         const BedWall& wall, double ambient_temperature) {
    (void)wall;
    auto& ho = state.external_heat_transfer_coefficient.var;
    const auto& Tw = state.wall_temperature.var;
    const auto& de = state.domain_vectors.de;
    const auto& dx = state.domain_vectors.cell_width;

    // One cylinder of the bed's height (all layers), as yang_vertical_cylinder.
    double height = 0.0;
    for (int i : state.domain_vectors.bed_cell_indices) height += dx[i];

    const double v = external_velocity;
    const double n = exponent;
    for (int i = layer.domain_.start_index; i <= layer.domain_.end_index; ++i) {
      const double T_wall = Tw[i] > 0.0 ? Tw[i] : ambient_temperature;
      ho[i] = ambient_air::MixedConvectionCoefficient(vertical_, T_wall, ambient_temperature, v, height, de[i], n);
    }
  }

  void ConstantParticleHeatTransferModel::Update(Equations& state, const MixedAdsorbentLayer& layer) {
    auto& hp = state.particle_heat_transfer_coefficient.var;
    // Specific surface area [1/m]: 6 (1 - eb) / dp for a packed bed, the
    // channel wall area per bed volume for a structured layer.
    const auto& av = state.domain_vectors.av;

    auto start = layer.domain_.start_index;
    auto end = layer.domain_.end_index;
    for (int n = start; n <= end; ++n) {
      hp[n] = heat_transfer_coefficient_ * av[n];   // Volumetric HTC [W/(m3.K)]
    }
  }

  void FilmParticleHeatTransferModel::Update(Equations& state, const MixedAdsorbentLayer& layer) {
    auto& hp = state.particle_heat_transfer_coefficient.var;
    const auto& dv = state.domain_vectors;
    const auto& Re = state.reynolds_number.var;
    const auto& k_g = state.fluid_thermal_conductivity.var;
    const auto& cp = state.heat_capacity.var;   // [J/(kg K)]
    const auto& mu = state.viscosity.var;
    const int start = layer.domain_.start_index, end = layer.domain_.end_index;
    for (int n = start; n <= end; ++n) {
      const double re = std::max(Re[n], 0.0);
      const double pr = k_g[n] > 0.0 ? std::max(cp[n] * mu[n] / k_g[n], 0.0) : 0.0;
      double h = 0.0;
      if (dv.channel_nu[n] > 0.0) {
        const double nu = ChannelFilmNumber(dv.channel_nu[n], re, pr, dv.hydraulic_diameter[n], dv.entrance_length[n]);
        h = nu * k_g[n] / dv.hydraulic_diameter[n];
      } else {
        const double nu = PackedFilmNumber(layer.film_correlation_, re, pr, layer.GetVoidFraction());
        h = nu * k_g[n] / dv.dp[n];
      }
      hp[n] = h * dv.av[n];
    }
  }

} // end adsim namespace


CEREAL_REGISTER_TYPE(adsim::ConstantParticleHeatTransferModel)
CEREAL_REGISTER_TYPE(adsim::FilmParticleHeatTransferModel)
CEREAL_REGISTER_TYPE(adsim::ConstantInternalHeatTransferModel)
CEREAL_REGISTER_TYPE(adsim::InternalHeatTransferModel_Dixon)
CEREAL_REGISTER_TYPE(adsim::InternalHeatTransferModel_Leva)
CEREAL_REGISTER_TYPE(adsim::InternalHeatTransferModel_LiFinlayson)
CEREAL_REGISTER_TYPE(adsim::InternalHeatTransferModel_MartinNilles)
CEREAL_REGISTER_TYPE(adsim::ConstantExternalHeatTransferModel)
CEREAL_REGISTER_TYPE(adsim::ExternalHeatTransferModel_Holman)
CEREAL_REGISTER_TYPE(adsim::ExternalHeatTransferModel_ChurchillChu)
CEREAL_REGISTER_TYPE(adsim::ExternalHeatTransferModel_YangCylinder)
CEREAL_REGISTER_TYPE(adsim::ExternalHeatTransferModel_ChurchillBernstein)
CEREAL_REGISTER_TYPE(adsim::ExternalHeatTransferModel_MixedConvection)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IInternalHeatTransferModel, adsim::ConstantInternalHeatTransferModel)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IInternalHeatTransferModel, adsim::InternalHeatTransferModel_Dixon)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IInternalHeatTransferModel, adsim::InternalHeatTransferModel_Leva)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IInternalHeatTransferModel, adsim::InternalHeatTransferModel_LiFinlayson)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IInternalHeatTransferModel, adsim::InternalHeatTransferModel_MartinNilles)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IExternalHeatTransferModel, adsim::ConstantExternalHeatTransferModel)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IExternalHeatTransferModel, adsim::ExternalHeatTransferModel_Holman)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IExternalHeatTransferModel, adsim::ExternalHeatTransferModel_ChurchillChu)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IExternalHeatTransferModel, adsim::ExternalHeatTransferModel_YangCylinder)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IExternalHeatTransferModel, adsim::ExternalHeatTransferModel_ChurchillBernstein)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IExternalHeatTransferModel, adsim::ExternalHeatTransferModel_MixedConvection)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IParticleHeatTransferModel, adsim::ConstantParticleHeatTransferModel)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IParticleHeatTransferModel, adsim::FilmParticleHeatTransferModel)
