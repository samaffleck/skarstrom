/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <string>
#include <memory>
#include <cmath>
#include <vector>

// Cereal includes
#include <cereal/types/memory.hpp>
#include <cereal/types/string.hpp>
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>

// adsim includes
#include "adsim/equations/transport_properties.h"
#include "adsim/math/expression_param.h"
#include "adsim/math/user_variables.h"


namespace adsim {

	// Forward declarations
	class PorousMedia;
	class Adsorbent;
  class MixedAdsorbentLayer;
  class Equations;
  class Fluid;
  class IIsotherm;
  class AdsorptionSource;
  using AdsorbentID = std::size_t; // mirrors porous_media.h; identical typedef redeclaration is legal
	
	enum class MassTransferModelType {
		LDF,
		LDF2
	};

	static const char* RateModelTypeList[] =  {
		"Linear Driving Force (1st Order)",
		"Linear Driving Force (2nd Order)"
	};

	enum class MassTransferCoefficientType {
		CONSTANT, 
		MACROPORE, 
		MACROPORE_KNUDSEN, 
		HAYNES_1973,
		ARRHENIUS,
		MACROPORE_CAPACITY
	};

	static const char* MassTransferCoefficientTypeList[] = {
		"Constant",		
		"Macropore Diffusion",
		"Macropore + Knudsen Diffusion",
		"Haynes (1973)",
		"Arrhenius",
		"Macropore Diffusion (capacity factor)"
	};

	class IMassTransferCoefficientModel {
	public:
		virtual ~IMassTransferCoefficientModel() = default;
		virtual MassTransferCoefficientType GetModelType() = 0;
		virtual unsigned ReadsTransportProperties() const { return kAllTransportProperties; }
		virtual void UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) = 0;
		// The adsorbent's isotherm for this component and its particle density
		// [kg/m3], handed over by AdsorptionSource before every update. Only the
		// capacity-factor model reads them (it evaluates the isotherm itself).
		// The component's isotherm, the adsorbent's solid density and the
		// adsorbent's AdsorptionSource (whose UpdateIsothermAll applies its
		// mixing rule), for the models that evaluate q* themselves.
		virtual void SetEquilibriumContext(IIsotherm* /*isotherm*/, double /*solid_density*/,
		                                   AdsorptionSource* /*source*/) {}
		virtual void EvaluateExpressions(const EvalContext& ctx) {}
		// Compiles the parameters evaluated per cell in the residual (see
		// ExpressionParam::CompileForRun). Called once at the start of a run.
		virtual void CompileForRun(const UserVariableRegistry* registry) {}

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar) {}
	};

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
		}

		ExpressionParam pore_radius_{1e-9};              // [m]
		ExpressionParam tortuosity_{3.0};                // particle tortuosity [-]
		ExpressionParam macropore_void_fraction_{0.2};   // [-]

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
		}
	};

	class MassTransferCoefficientMacropore : public IMassTransferCoefficientModel {
	public:
    void UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) override;
		MassTransferCoefficientType GetModelType() override { return MassTransferCoefficientType::MACROPORE; } // Macropore limited
		unsigned ReadsTransportProperties() const override { return kMolecularDiffusion; }
		void EvaluateExpressions(const EvalContext& ctx) override {
			tortuosity_.Evaluate(ctx); macropore_void_fraction_.Evaluate(ctx);
		}

		ExpressionParam tortuosity_{3.0};                // particle tortuosity [-]
		ExpressionParam macropore_void_fraction_{0.35};  // [-]

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
		}

		ExpressionParam tortuosity_{3.0};                // particle tortuosity [-]
		ExpressionParam macropore_void_fraction_{0.35};  // [-]
		ExpressionParam concentration_floor_{1e-3};      // [mol/m3] c below which the Henry-limit ratio is used

	private:
		IIsotherm* isotherm_ = nullptr;   // transient, set by AdsorptionSource
		AdsorptionSource* source_ = nullptr;  // transient: under IAST, q* comes from its UpdateIsothermAll
		double solid_density_ = 0.0;      // transient [kg/m3]
		std::vector<std::vector<double>> q_all_scratch_;  // transient: every component's q* before the IAST re-evaluation
		std::vector<double> y_scratch_;   // transient: the mole fraction before flooring
		std::vector<double> q_scratch_;   // transient: q* at the floored concentration

		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const /*version*/) {
			ar(
				cereal::virtual_base_class<IMassTransferCoefficientModel>(this),
				cereal::make_nvp("Tortuosity", tortuosity_),
				cereal::make_nvp("Macropore_Void_Fraction", macropore_void_fraction_),
				cereal::make_nvp("Concentration_Floor", concentration_floor_)
			);
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

	class IMassTransferModel {
	public:
		virtual ~IMassTransferModel() = default;
		virtual MassTransferModelType GetModelType() = 0;

		// Kinetic uptake rate dq/dt [mol/kg/s] for equilibrium loading q_star,
		// current loading q (both [mol/kg]) and coefficient k from mtc_model_.
		// This is the single source of truth for the loading residual
		// (UpdateAdsorbedConcentration); the energy and mass source terms all use
		// the resulting dq/dt, so they stay consistent automatically.
		virtual double Rate(double q_star, double q, double k) const = 0;

		void EvaluateExpressions(const EvalContext& ctx) {
			if (mtc_model_) mtc_model_->EvaluateExpressions(ctx);
		}
		void CompileForRun(const UserVariableRegistry* registry) {
			if (mtc_model_) mtc_model_->CompileForRun(registry);
		}

    std::unique_ptr<IMassTransferCoefficientModel> mtc_model_;

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar) {
      ar(
				cereal::make_nvp("Mass_Transfer_Coefficient_Model", mtc_model_)
      );
    }
	};

	class LinearDrivingForce : public IMassTransferModel {
	public:
		LinearDrivingForce();

    MassTransferModelType GetModelType() override { return MassTransferModelType::LDF; }

		// dq/dt = k (q* - q), k in [1/s].
		double Rate(double q_star, double q, double k) const override { return k * (q_star - q); }

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IMassTransferModel>(this)
			);
		}
	};

	class LinearDrivingForceSecondOrder : public IMassTransferModel {
	public:
		LinearDrivingForceSecondOrder();

		MassTransferModelType GetModelType() override { return MassTransferModelType::LDF2; }

		// dq/dt = k (q* - q) |q* - q|, k in [1/(s mol/kg)].
		// Sign-preserving form of the second-order driving force: the literal
		// k (q* - q)^2 is never negative and so could never desorb. The derivative
		// d(rate)/dq = -2 k |q* - q| is continuous (zero at q = q*), so plain abs()
		// needs no smoothing for the Newton iteration.
		double Rate(double q_star, double q, double k) const override {
			const double d = q_star - q;
			return k * d * std::abs(d);
		}

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IMassTransferModel>(this)
			);
		}
	};

	// Factory functions
	inline std::unique_ptr<IMassTransferModel> CreateMassTransferModel(MassTransferModelType type) {
		switch (type) {
		case MassTransferModelType::LDF: return std::make_unique<LinearDrivingForce>();
		case MassTransferModelType::LDF2: return std::make_unique<LinearDrivingForceSecondOrder>();
		default: return nullptr;
		}
	}

	inline std::unique_ptr<IMassTransferCoefficientModel> CreateMassTransferCoefficientModel(MassTransferCoefficientType type) {
		switch (type) {
		case MassTransferCoefficientType::CONSTANT: return std::make_unique<MassTransferCoefficientConstant>();
		case MassTransferCoefficientType::MACROPORE: return std::make_unique<MassTransferCoefficientMacropore>();
		case MassTransferCoefficientType::MACROPORE_KNUDSEN: return std::make_unique<MassTransferCoefficientMacroporeKnudsen>();
		case MassTransferCoefficientType::HAYNES_1973: return std::make_unique<MassTransferCoefficient_Haynes>();
		case MassTransferCoefficientType::ARRHENIUS: return std::make_unique<MassTransferCoefficientArrhenius>();
		case MassTransferCoefficientType::MACROPORE_CAPACITY: return std::make_unique<MassTransferCoefficientMacroporeCapacity>();
		default: return nullptr;
		}
	}

} // end adsim namespace
