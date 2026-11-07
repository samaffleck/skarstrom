/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/objects/fluid.h"
#include "adsim/objects/viscosity_coefficients.h"
#include "adsim/solver/initial_conditions.h"


namespace adsim {

	// Convenience getter: resolve the hard-coded coefficients and evaluate. Falls
	// back to a nominal dilute-gas value for an unknown component (one not in
	// kViscosityCoefficients).
	double Fluid::GetComponentViscosity(const std::string& component, double T) {
		const auto* k = ViscosityCoeffs(component);
		return k ? EvalViscosity(*k, T) : 1.0e-5;
	}

	// Convenience getter: hard-coded conductivity coefficients + eval. Falls back
	// to a nominal air-like value for an unknown component.
	double Fluid::GetComponentThermalConductivity(const std::string& component, double T) {
		const auto* k = ThermalConductivityCoeffs(component);
		return k ? EvalThermalConductivity(*k, T) : 0.024;
	}

	double GetGasHeatCapacity(const std::string& component, double P, double T) {
		try {
			const auto phase = static_cast<CoolProp::phases>(
				static_cast<int>(CoolProp::PropsSI("Phase", "P", P, "T", T, component)));
			if (phase == CoolProp::iphase_liquid || phase == CoolProp::iphase_twophase)
				return CoolProp::PropsSI("CP0MASS", "P", P, "T", T, component); // ideal-gas [J/(kg-K)]
		}
		catch (const CoolProp::CoolPropBaseError& e) {
			// No phase (or C_p0) for this state: fall back to the state value below.
			std::cerr << "CoolProp Error: " << e.what() << std::endl;
		}
		return GetHeatCapacity(component, P, T);
	}

	void Fluid::AddComponent(const std::string& component_name) {
		assert(IsComponentValid(component_name));

		if (IsComponentInList(component_name)) 
			return;
		components_.emplace_back(component_name);
	}

	void Fluid::RemoveComponent(const std::string& component_name) {
		assert(IsComponentValid(component_name));

    if (!IsComponentInList(component_name)) 
			return;
		components_.erase(std::find(components_.begin(), components_.end(), component_name));
	}

	void Fluid::Initialise(InitialConditions& ic) {
		molar_mass_.clear();
		heat_capacity_.clear();
		cp_polynomial_.clear();
		for (const std::string& component : components_) {
			molar_mass_[component] = GetMolarMass(component);
			auto cp_it = gas_heat_capacity_.try_emplace(component, kDefaultGasHeatCapacity).first;
			if (gas_heat_capacity_mode_ == GasHeatCapacityMode::CONSTANT) {
				heat_capacity_[component] = cp_it->second.Get() / molar_mass_[component];  // [J/(mol-K)] -> [J/(kg-K)]
			} else if (const NasaPolynomial* p = IdealGasCoeffs(component)) {
				// The solver evaluates C_p0(T) per cell; the cached value is for display.
				cp_polynomial_[component] = p;
				heat_capacity_[component] = EvalIdealGasHeatCapacity(*p, ic.temperature) / molar_mass_[component];
			} else {
				// No polynomial: CoolProp at the initial conditions, held for the run.
				heat_capacity_[component] = GetGasHeatCapacity(component, ic.pressure, ic.temperature);
			}
			// Viscosity mu_i(T) uses hard-coded polynomial coefficients
			// (viscosity_coefficients.h) — no CoolProp fitting needed here.

			if (adsorbed_heat_capacity_.find(component) == adsorbed_heat_capacity_.end()) {
				adsorbed_heat_capacity_[component] = 0.0;
      }
		}
	}

	bool Fluid::PerformChecks(std::string& error_message) const {
		if (components_.empty()) {
			error_message += "Component list is empty. You need to add at least one component to the system. ";
			return false;
		}
		if (gas_heat_capacity_mode_ == GasHeatCapacityMode::CONSTANT) {
			for (const auto& c : components_) {
				auto it = gas_heat_capacity_.find(c);
				if (it != gas_heat_capacity_.end() && !(it->second.Get() > 0.0)) {
					error_message += "The constant gas heat capacity of " + c + " must be positive. ";
					return false;
				}
			}
		}
		if (viscosity_mode_ == ViscosityMode::CONSTANT && !(viscosity_.Get() > 0.0)) {
			error_message += "The constant gas viscosity must be positive. ";
			return false;
		}

		return true;
	}

	bool Fluid::IsComponentValid(const std::string& component_name) const {
		for (const auto& component : Components::componentNames) {
			if (component == component_name)
				return true;
		}

		return false;
	}

	bool Fluid::IsComponentInList(const std::string& component_name) const {
		for (const auto& component : components_) {
			if (component == component_name)
				return true;
		}

		return false;
	}

} // end adsim namespace
