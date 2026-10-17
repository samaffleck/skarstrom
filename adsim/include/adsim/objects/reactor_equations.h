/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <memory>

// adsim includes
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/equations/binary_diffusion_model.h"
#include "adsim/equations/wall_temperature_boundary_condition.h"
#include "adsim/math/expression_param.h"

// Cereal includes
#include <cereal/cereal.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>


namespace adsim {

  enum class ThermalEquilibriumMode {
    LOCAL,      // LTE: T_fluid = T_solid (single combined equation)
    NON_LOCAL   // LTNE: separate fluid and solid temperature equations
  };

  static const char* ThermalEquilibriumModeList[] = {
    "Local thermal equilibrium (LTE)",
    "Local thermal non-equilibrium (LTNE)"
  };

  // Whether the energy balance is solved at all. Orthogonal to
  // ThermalEquilibriumMode (LTE / LTNE), which only says HOW it is solved: an
  // isothermal run removes the fluid, solid and wall temperature unknowns from
  // the DAE and holds them at a fixed temperature, and switching back restores
  // the user's LTE / LTNE choice untouched. Values are stored in the deck.
  enum class EnergyBalanceMode {
    NON_ISOTHERMAL = 0, // solve the energy balance(s) (default)
    ISOTHERMAL = 1      // T_fluid = T_solid = T_wall = T_iso, no energy equations
  };

  static const char* EnergyBalanceModeList[] = {
    "Non-isothermal (solve energy balance)",
    "Isothermal (fixed temperature)"
  };

  class EnergySource {
  public:
    bool enabled_ = false;
    ExpressionParam power_{0.0}; // Total power [W]. Positive = heating, negative = cooling.

    void EvaluateExpressions(const EvalContext& ctx) { power_.Evaluate(ctx); }

  private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
        cereal::make_nvp("Power", power_)
			);
			if (version >= 1) {
				ar(cereal::make_nvp("Enabled", enabled_));
			}
		}
  };

	class ReactorEquations {
	public:
		ReactorEquations();

		BinaryDiffusionCoefficientModel binary_diffusion_model_{};
		std::unique_ptr<IWallTemperatureBoundaryCondition> wall_temperature_boundary_condition_;
		ThermalEquilibriumMode thermal_equilibrium_mode_ = ThermalEquilibriumMode::LOCAL;
		std::unique_ptr<IParticleHeatTransferModel> particle_heat_transfer_model_;

		// Isothermal switch. In ISOTHERMAL mode every temperature is held at
		// T_iso: the initial-condition (fluid) temperature when
		// isothermal_use_initial_temperature_ is set (the default), otherwise
		// isothermal_temperature_.
		EnergyBalanceMode energy_balance_mode_ = EnergyBalanceMode::NON_ISOTHERMAL;
		bool isothermal_use_initial_temperature_ = true;
		ExpressionParam isothermal_temperature_{298.15}; // [K]

		bool IsIsothermal() const { return energy_balance_mode_ == EnergyBalanceMode::ISOTHERMAL; }
		double GetIsothermalTemperature(double initial_temperature) const {
			return isothermal_use_initial_temperature_ ? initial_temperature : isothermal_temperature_.Get();
		}

		void EvaluateExpressions(const EvalContext& ctx) {
			if (particle_heat_transfer_model_) particle_heat_transfer_model_->EvaluateExpressions(ctx);
			if (wall_temperature_boundary_condition_) wall_temperature_boundary_condition_->EvaluateExpressions(ctx);
			isothermal_temperature_.Evaluate(ctx);
		}

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
        cereal::make_nvp("Wall_Temperature_Boundary_Condition", wall_temperature_boundary_condition_)
      );
			if (version >= 1) {
				ar(
					cereal::make_nvp("Thermal_Equilibrium_Mode", thermal_equilibrium_mode_),
					cereal::make_nvp("Particle_Heat_Transfer_Model", particle_heat_transfer_model_)
				);
			}
			if (version >= 2) {
				ar(
					cereal::make_nvp("Energy_Balance_Mode", energy_balance_mode_),
					cereal::make_nvp("Isothermal_Temperature_From_Initial_Conditions", isothermal_use_initial_temperature_),
					cereal::make_nvp("Isothermal_Temperature", isothermal_temperature_)
				);
			}
			// version < 2: the members keep their (non-isothermal) defaults.
		}
	};

	class LayerEquations {
	public:
		LayerEquations();

		std::unique_ptr<IDispersionModel> dispersion_model_;
		std::unique_ptr<IHeatDispersionModel> heat_dispersion_model_;
		std::unique_ptr<IInternalHeatTransferModel> internal_heat_transfer_model_;
		std::unique_ptr<IExternalHeatTransferModel> external_heat_transfer_model_;

		void EvaluateExpressions(const EvalContext& ctx) {
			if (dispersion_model_) dispersion_model_->EvaluateExpressions(ctx);
			if (heat_dispersion_model_) heat_dispersion_model_->EvaluateExpressions(ctx);
			if (internal_heat_transfer_model_) internal_heat_transfer_model_->EvaluateExpressions(ctx);
			if (external_heat_transfer_model_) external_heat_transfer_model_->EvaluateExpressions(ctx);
		}

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::make_nvp("Dispersion_Model", dispersion_model_),
				cereal::make_nvp("Heat_Dispersion_Model", heat_dispersion_model_),
				cereal::make_nvp("Internal_Heat_Transfer_Model", internal_heat_transfer_model_),
				cereal::make_nvp("External_Heat_Transfer_Model", external_heat_transfer_model_)
			);
			if (version >= 1) {
				// Backward compat: these fields moved to ReactorEquations, discard on load
				ThermalEquilibriumMode dummy_mode{};
				std::unique_ptr<IParticleHeatTransferModel> dummy_model;
				ar(
					cereal::make_nvp("Thermal_Equilibrium_Mode", dummy_mode),
					cereal::make_nvp("Particle_Heat_Transfer_Model", dummy_model)
				);
			}
		}
	};

} // end adsim namespace

CEREAL_CLASS_VERSION(adsim::EnergySource, 2);
CEREAL_CLASS_VERSION(adsim::ReactorEquations, 2);
CEREAL_CLASS_VERSION(adsim::LayerEquations, 1);
