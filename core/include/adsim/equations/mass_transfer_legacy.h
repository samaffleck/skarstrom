/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// The rate-coefficient classes that each held one fixed combination of
// resistances (constant, arrhenius, macropore, macropore_knudsen, haynes,
// macropore_capacity) and the first version of the series-resistance class.
// They are load-only: decks that store them still load, and
// IMassTransferModel's loader converts each to the generic coefficient model
// (MassTransferCoefficientSeriesResistance) with the switches, capacity-factor
// mode and parameter values that reproduce it bit for bit (ToResistances).
// Their evaluation code is kept as the reference the generic model is tested
// against. Nothing creates them except the loader and the tests
// (CreateLegacyMassTransferCoefficientModel).

#include "adsim/equations/mass_transfer_model.h"

namespace adsim {

	class MassTransferCoefficientConstant : public IMassTransferCoefficientModel {
	public:
		void UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) override;
    MassTransferCoefficientType GetModelType() override { return MassTransferCoefficientType::CONSTANT; }
    unsigned ReadsTransportProperties() const override { return kNoTransportProperties; }
    void EvaluateExpressions(const EvalContext& ctx) override { constant_mtc_.Evaluate(ctx); }
    void CompileForRun(const UserVariableRegistry* registry) override { constant_mtc_.CompileForRun(registry); }

		ExpressionParam constant_mtc_{0.5};

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const /*version*/) {
			ar(cereal::virtual_base_class<IMassTransferCoefficientModel>(this),
			   cereal::make_nvp("Mass_Transfer_Coefficient", constant_mtc_));
		}
	};

	class MassTransferCoefficientArrhenius : public IMassTransferCoefficientModel {
	public:
    void UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) override;
		MassTransferCoefficientType GetModelType() override { return MassTransferCoefficientType::ARRHENIUS; }
		unsigned ReadsTransportProperties() const override { return kNoTransportProperties; }
		void EvaluateExpressions(const EvalContext& ctx) override {
			pre_exp_factor_.Evaluate(ctx); activation_energy_.Evaluate(ctx);
		}
		void CompileForRun(const UserVariableRegistry* registry) override {
			pre_exp_factor_.CompileForRun(registry); activation_energy_.CompileForRun(registry);
		}

		ExpressionParam pre_exp_factor_{120.0};
		ExpressionParam activation_energy_{1e4};

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const /*version*/) {
			ar(cereal::virtual_base_class<IMassTransferCoefficientModel>(this),
			   cereal::make_nvp("Pre_Exponential_Factor", pre_exp_factor_),
			   cereal::make_nvp("Activation_Energy", activation_energy_));
		}
	};

	class MassTransferCoefficientMacroporeKnudsen : public IMassTransferCoefficientModel {
	public:
    void UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) override;
		MassTransferCoefficientType GetModelType() override { return MassTransferCoefficientType::MACROPORE_KNUDSEN; } // Macropore limited
		unsigned ReadsTransportProperties() const override { return kMolecularDiffusion; }
		void EvaluateExpressions(const EvalContext& ctx) override {
			pore_radius_.Evaluate(ctx); tortuosity_.Evaluate(ctx); macropore_void_fraction_.Evaluate(ctx);
			molecular_diffusion_coefficient_.Evaluate(ctx);
		}

		ExpressionParam pore_radius_{1e-9};              // [m]
		ExpressionParam tortuosity_{3.0};                // particle tortuosity [-]
		ExpressionParam macropore_void_fraction_{0.2};   // [-]
		// Constant molecular diffusivity D_m [m2/s] used in place of the
		// state's (T, P)-dependent one when > 0; 0 (the default) keeps the
		// binary-diffusion model's value.
		ExpressionParam molecular_diffusion_coefficient_{0.0};

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IMassTransferCoefficientModel>(this),
				cereal::make_nvp("Pore_Radius", pore_radius_),
				cereal::make_nvp("Tortuosity", tortuosity_),
				cereal::make_nvp("Macropore_Void_Fraction", macropore_void_fraction_)
			);
			// v3 added the optional constant D_m; older archives keep 0 (state D_m).
			if (version >= 3) ar(cereal::make_nvp("Molecular_Diffusion_Coefficient", molecular_diffusion_coefficient_));
		}
	};

	class MassTransferCoefficientMacropore : public IMassTransferCoefficientModel {
	public:
    void UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) override;
		MassTransferCoefficientType GetModelType() override { return MassTransferCoefficientType::MACROPORE; } // Macropore limited
		unsigned ReadsTransportProperties() const override { return kMolecularDiffusion; }
		void EvaluateExpressions(const EvalContext& ctx) override {
			tortuosity_.Evaluate(ctx); macropore_void_fraction_.Evaluate(ctx);
			molecular_diffusion_coefficient_.Evaluate(ctx);
		}

		ExpressionParam tortuosity_{3.0};                // particle tortuosity [-]
		ExpressionParam macropore_void_fraction_{0.35};  // [-]
		// Constant molecular diffusivity D_m [m2/s] used in place of the
		// state's (T, P)-dependent one when > 0 (Haghpanah et al. 2013 and
		// Estupinan Perez et al. 2019 use a constant 1.6e-5 m2/s); 0 (the
		// default) keeps the binary-diffusion model's value.
		ExpressionParam molecular_diffusion_coefficient_{0.0};

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IMassTransferCoefficientModel>(this),
				cereal::make_nvp("Tortuosity", tortuosity_),
				cereal::make_nvp("Macropore_Void_Fraction", macropore_void_fraction_)
			);
			// v3 added the optional constant D_m; older archives keep 0 (state D_m).
			if (version >= 3) ar(cereal::make_nvp("Molecular_Diffusion_Coefficient", molecular_diffusion_coefficient_));
		}
	};

	// Macropore-diffusion LDF coefficient with the capacity factor of the
	// Haghpanah et al. (2013) / Estupinan Perez et al. (2019, SI Eq. S4) model,
	//   k_i = (15 eps_p D_m,i / (tau r_p^2)) * c_i / (rho_p q_i*),
	// the macropore coefficient scaled by the ratio of the gas-phase to the
	// adsorbed-phase capacity (mol per m3 of gas over mol per m3 of particle).
	// For a strongly adsorbed component it is orders of magnitude below the plain
	// macropore value. Evaluated per cell from the isotherm itself: as c -> 0 both
	// c and q* vanish and their ratio tends to 1 / (rho_p K_H), so the ratio is
	// taken at max(c_i, c_floor), with the isotherm re-evaluated at that
	// concentration (every other component at its actual value), which stays
	// finite and smooth where the state's own c / q* is 0 / 0. The ratio is
	// capped at 1 (the plain macropore coefficient). Under the adsorbent's IAST
	// mixing rule, q* is the IAST mixture loading at that state (all of the
	// adsorbent's components solved together), not the component's own
	// isotherm on its own.
	class MassTransferCoefficientMacroporeCapacity : public IMassTransferCoefficientModel {
	public:
    void UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) override;
		MassTransferCoefficientType GetModelType() override { return MassTransferCoefficientType::MACROPORE_CAPACITY; }
		unsigned ReadsTransportProperties() const override { return kMolecularDiffusion; }
		void SetEquilibriumContext(IIsotherm* isotherm, double solid_density, AdsorptionSource* source) override {
			isotherm_ = isotherm; solid_density_ = solid_density; source_ = source;
		}
		void EvaluateExpressions(const EvalContext& ctx) override {
			tortuosity_.Evaluate(ctx); macropore_void_fraction_.Evaluate(ctx); concentration_floor_.Evaluate(ctx);
			molecular_diffusion_coefficient_.Evaluate(ctx);
		}

		ExpressionParam tortuosity_{3.0};                // particle tortuosity [-]
		ExpressionParam macropore_void_fraction_{0.35};  // [-]
		ExpressionParam concentration_floor_{1e-3};      // [mol/m3] c below which the Henry-limit ratio is used
		// Constant molecular diffusivity D_m [m2/s] used in place of the
		// state's (T, P)-dependent one when > 0 (the papers above use a
		// constant 1.6e-5 m2/s); 0 (the default) keeps the binary-diffusion
		// model's value.
		ExpressionParam molecular_diffusion_coefficient_{0.0};

	private:
		IIsotherm* isotherm_ = nullptr;   // transient, set by AdsorptionSource
		AdsorptionSource* source_ = nullptr;  // transient: under IAST, q* comes from its UpdateIsothermAll
		double solid_density_ = 0.0;      // transient [kg/m3]
		std::vector<std::vector<double>> q_all_scratch_;  // transient: every component's q* before the IAST re-evaluation
		std::vector<double> y_scratch_;   // transient: the mole fraction before flooring
		std::vector<double> q_scratch_;   // transient: q* at the floored concentration

		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IMassTransferCoefficientModel>(this),
				cereal::make_nvp("Tortuosity", tortuosity_),
				cereal::make_nvp("Macropore_Void_Fraction", macropore_void_fraction_),
				cereal::make_nvp("Concentration_Floor", concentration_floor_)
			);
			// v2 added the optional constant D_m; older archives keep 0 (state D_m).
			if (version >= 2) ar(cereal::make_nvp("Molecular_Diffusion_Coefficient", molecular_diffusion_coefficient_));
		}
	};

	class MassTransferCoefficient_Haynes : public IMassTransferCoefficientModel {
	public:
    void UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) override;
		MassTransferCoefficientType GetModelType() override { return MassTransferCoefficientType::HAYNES_1973; } // Macropore limited
		unsigned ReadsTransportProperties() const override { return kMolecularDiffusion | kSherwoodNumber; }
		void EvaluateExpressions(const EvalContext& ctx) override {
			tortuosity_.Evaluate(ctx); macropore_void_fraction_.Evaluate(ctx);
		}

		ExpressionParam tortuosity_{3.0};                // particle tortuosity [-]
		ExpressionParam macropore_void_fraction_{0.2};   // [-]

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IMassTransferCoefficientModel>(this),
				cereal::make_nvp("Tortuosity", tortuosity_),
				cereal::make_nvp("Macropore_Void_Fraction", macropore_void_fraction_)
			);
		}
	};

	// The series-resistance coefficient as it was before the per-resistance
	// switches (cereal version 1 of MassTransferCoefficientSeriesResistance):
	// each term switched off by a zero parameter. Not archived: version-1
	// archives load into the generic class, which derives the switches from
	// these zeros. Kept as the reference the generic model is tested against.
	//
	// Series-resistance (Glueckauf) LDF coefficient: external film, macropore
	// and micropore resistances added in series, the first two scaled by the
	// local capacity factor K (Ruthven 1984; as used for direct-air-capture
	// breakthrough curves by Stampi-Bombelli et al. 2024, Ind. Eng. Chem. Res.
	// 63, 11637, doi:10.1021/acs.iecr.4c01392). For packed spheres
	//   1/k = f_film K r_p / (3 k_f) + K tau r_p^2 / (15 eps_p D_m)
	//         + r_c^2 / (15 D_c),        D_c = D_c0 exp(-E_a / (R T)),
	// with k_f = Sh D_m / d_p from the layer's film correlation (as haynes),
	// and for a structured layer (slab of half-thickness l) the haynes slab
	// form: (V/S) / k_f and tau l^2 / (3 eps_p D_m). K = rho_p dq*/dc, the
	// slope of the isotherm (q* per kg, c in mol per m3 of gas, rho_p the
	// particle density), by a forward difference at max(c, c_floor); where the
	// slope is not positive the secant rho_p q*/c is used, and K is floored at
	// 1 (the haynes coefficient), as macropore_capacity caps its ratio. A
	// capacity_factor > 0 fixes K instead. Each term is switched off by its
	// parameter: film_resistance_factor = 0, tortuosity = 0 (D_p = D_m / tau
	// infinite), crystal_radius = 0 or micropore_diffusivity = 0. With the
	// micropore term off and K = 1 the coefficient equals haynes exactly.
	class LegacySeriesResistanceV1 : public IMassTransferCoefficientModel {
	public:
    void UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) override;
		MassTransferCoefficientType GetModelType() override { return MassTransferCoefficientType::SERIES_RESISTANCE; }
		unsigned ReadsTransportProperties() const override { return kMolecularDiffusion | kSherwoodNumber; }
		void SetEquilibriumContext(IIsotherm* isotherm, double solid_density, AdsorptionSource* source) override {
			isotherm_ = isotherm; solid_density_ = solid_density; source_ = source;
		}
		void EvaluateExpressions(const EvalContext& ctx) override {
			film_resistance_factor_.Evaluate(ctx); tortuosity_.Evaluate(ctx); macropore_void_fraction_.Evaluate(ctx);
			molecular_diffusion_coefficient_.Evaluate(ctx); crystal_radius_.Evaluate(ctx);
			micropore_diffusivity_.Evaluate(ctx); micropore_activation_energy_.Evaluate(ctx);
			capacity_factor_.Evaluate(ctx); concentration_floor_.Evaluate(ctx);
		}

		// The capacity factor K at each cell of [start, end] (written into k_out),
		// from the isotherm as described above. Public for the unit tests.
		void CapacityFactor(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component,
		                    int start, int end, const UserVariableRegistry* registry, std::vector<double>& k_out);

		ExpressionParam film_resistance_factor_{1.0};       // [-] multiplies the film term; 0 drops it
		ExpressionParam tortuosity_{3.0};                   // particle tortuosity [-]; 0 drops the macropore term
		ExpressionParam macropore_void_fraction_{0.35};     // [-]
		// Constant molecular diffusivity D_m [m2/s] of the macropore term when
		// > 0; 0 (the default) keeps the state's binary-diffusion value. The
		// film term always uses the state's D_m, as the film correlation does.
		ExpressionParam molecular_diffusion_coefficient_{0.0};
		ExpressionParam crystal_radius_{0.0};               // r_c [m]; 0 drops the micropore term
		ExpressionParam micropore_diffusivity_{0.0};        // D_c0 [m2/s]; 0 drops the micropore term
		ExpressionParam micropore_activation_energy_{0.0};  // E_a [J/mol]
		ExpressionParam capacity_factor_{0.0};              // K [-] when > 0; 0: from the isotherm
		ExpressionParam concentration_floor_{1e-3};         // [mol/m3] c below which K is taken at c_floor

	private:
		// q* of this component at the state's current mole fractions into out;
		// the state's q* (every component) is left as it was.
		void EvaluateLoading(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component,
		                     int start, int end, const UserVariableRegistry* registry, std::vector<double>& out);

		IIsotherm* isotherm_ = nullptr;       // transient, set by AdsorptionSource
		AdsorptionSource* source_ = nullptr;  // transient: under IAST, q* comes from its UpdateIsothermAll
		double solid_density_ = 0.0;          // transient [kg/m3]
		std::vector<std::vector<double>> q_all_scratch_;  // transient
		std::vector<double> y_scratch_, q_base_, q_pert_, capacity_;  // transient

	};


	// The old class of a coefficient kind with its old defaults
	// (SERIES_RESISTANCE: LegacySeriesResistanceV1). For the loader's
	// conversion and the tests; decks never get these.
	std::unique_ptr<IMassTransferCoefficientModel> CreateLegacyMassTransferCoefficientModel(MassTransferCoefficientType type);

} // end adsim namespace

CEREAL_CLASS_VERSION(adsim::MassTransferCoefficientConstant, 3);
CEREAL_CLASS_VERSION(adsim::MassTransferCoefficientMacroporeKnudsen, 3);
CEREAL_CLASS_VERSION(adsim::MassTransferCoefficientMacropore, 3);
CEREAL_CLASS_VERSION(adsim::MassTransferCoefficient_Haynes, 2);
CEREAL_CLASS_VERSION(adsim::MassTransferCoefficientMacroporeCapacity, 2);
CEREAL_CLASS_VERSION(adsim::MassTransferCoefficientArrhenius, 3);
