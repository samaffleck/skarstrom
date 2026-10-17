/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/equations/binary_diffusion_model.h"
#include "adsim/objects/fluid.h"


namespace adsim {

	void BinaryDiffusionCoefficientModel::Initialise(int number_of_cells, const std::vector<std::string>& component_list) {
		m_DiffusionCoefficients.clear();

		// ComponentPair hashes and compares symmetrically, so only the upper
		// triangle is needed; every entry here is symmetric in (1, 2) anyway.
		//
		// A single-component fluid has no binary pair, so store the self pair
		// {i, i} instead: Chapman-Enskog with sigma_ii, eps_ii and M_AB = M_i is
		// the self-diffusion coefficient, which diffusion-based MTC models need
		// for a non-zero D_m. Multi-component fluids keep the off-diagonal pairs
		// only (no self pairs are added or evaluated for them).
		const bool self_pair_only = component_list.size() == 1;
		for (size_t i = 0; i < component_list.size(); ++i) {
			for (size_t j = self_pair_only ? i : i + 1; j < component_list.size(); ++j) {
				const auto& component1 = component_list[i];
				const auto& component2 = component_list[j];

				BinaryLenardJonesParamaters params{};
				params.omega.resize(number_of_cells, 0.0);
				params.Dm_ij.resize(number_of_cells, 0.0);

				double mr_1 = CoolProp::Props1SI(component1, "molar_mass") * 1e3;		// molecular weight in [g/mol]
				double mr_2 = CoolProp::Props1SI(component2, "molar_mass") * 1e3;		// molecular weight in [g/mol]
				params.mr_avr_invers_ij = (1 / mr_1) + (1 / mr_2);

				const ComponentData& cd_1 = Components::GetComponentData(component1);
				const ComponentData& cd_2 = Components::GetComponentData(component2);

				params.eps = std::sqrt(cd_1.eps * cd_2.eps);
				params.sigma = 0.5 * (cd_1.sigma + cd_2.sigma);

				m_DiffusionCoefficients[{component1, component2}] = std::move(params);
			}
		}
	}

	void BinaryDiffusionCoefficientModel::Update(const std::vector<double>& T, const std::vector<double>& P) {
		auto N = (int)T.size();

		// The numerator 2.66e-7 * T^1.5 is identical for every component pair —
		// compute it once per cell instead of once per pair per cell.
		top_scratch_.resize(N);
		for (int n = 0; n < N; ++n) {
			top_scratch_[n] = 2.66e-7 * std::pow(T[n], (3.0 / 2.0));
		}

		for (auto& [componentPair, params] : m_DiffusionCoefficients) {
			auto& omega_D = params.omega;
			auto& moldiff_ij = params.Dm_ij;

			double eps = params.eps;
			double sigma_ij = params.sigma;
			double mr_ij = params.mr_avr_invers_ij;

			// Pair-constant factors hoisted out of the cell loop.
			const double M_AB = 2 / mr_ij;
			const double sqrt_M_AB = std::sqrt(M_AB);

			for (int n = 0; n < N; ++n) {
        // Collision diameter correlation of Neufield et al. (1972)
        double Ts = T[n] / eps;
        double a = 1.06036 / (std::pow(Ts, 0.15610));
        double c = 0.193 / (std::exp(0.47635 * Ts));
        double e = 1.03587 / (std::exp(1.52996 * Ts));
        double g = 1.76474 / (std::exp(3.89411 * Ts));
        omega_D[n] = a + c + e + g;

        double P_bar = P[n] / 1.0e5;
        double top = top_scratch_[n];
        double bottom = P_bar * sqrt_M_AB * sigma_ij * sigma_ij * omega_D[n];
        moldiff_ij[n] = top / bottom;
			}
		}
	}

} // adsim namespace
