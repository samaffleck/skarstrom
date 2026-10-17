/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/solver/initial_conditions.h"

namespace adsim {

  bool InitialConditions::PerformChecks(const Fluid& fluid, std::string& error_message) const {
    double total_mole_fraction = 0;
    for (const auto& [component, yi] : mole_fraction) {
      // only count components that are in the fluid list
      for (const auto& comp : fluid.components_) {
        if (component == comp)
          total_mole_fraction += yi;
      }
    }

    if (std::abs(1 - total_mole_fraction) > 1e-6) {
      error_message += "Sum of mole fractions in the initial conditions list does not sum to 1. ";
      
      return false;
    }

    return true;
  }

} // end adsim namepsace
