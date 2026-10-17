/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/objects/reactor_equations.h"


namespace adsim {

  ReactorEquations::ReactorEquations() {
    wall_temperature_boundary_condition_ = std::make_unique<WallTemperatureBoundaryConditionZeroFlux>();
    particle_heat_transfer_model_ = std::make_unique<ConstantParticleHeatTransferModel>();
  }

  LayerEquations::LayerEquations() {
    dispersion_model_ = std::make_unique<DispersionModelConstant>();
    heat_dispersion_model_ = std::make_unique<HeatDispersionModelConstant>();
    internal_heat_transfer_model_ = std::make_unique<ConstantInternalHeatTransferModel>();
    external_heat_transfer_model_ = std::make_unique<ConstantExternalHeatTransferModel>();
  }

} // end adsim namespace
