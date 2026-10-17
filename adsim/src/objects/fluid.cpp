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
		for (const std::string& component : components_) {
			molar_mass_[component] = GetMolarMass(component);
			heat_capacity_[component] = GetHeatCapacity(component, ic.pressure, ic.temperature);
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


CEREAL_CLASS_VERSION(adsim::Fluid, 3);
