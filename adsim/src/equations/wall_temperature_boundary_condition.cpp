/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/equations/wall_temperature_boundary_condition.h"
#include "adsim/equations/equation.h"


namespace adsim {

	void WallTemperatureBoundaryConditionFixed::Update(Equations& state) {
    auto& Tw = state.wall_temperature.var;
    auto& res_Tw = state.wall_temperature.res_var;
    int N = (int)Tw.size() - 1;

    res_Tw[0] = Tw[0] + Tw[1] - 2 * temperature_;
    res_Tw[N] = Tw[N] + Tw[N - 1] - 2 * temperature_;
	}

	void WallTemperatureBoundaryConditionZeroFlux::Update(Equations& state) {
    auto& Tw = state.wall_temperature.var;
    auto& res_Tw = state.wall_temperature.res_var;
    int N = (int)Tw.size() - 1;

    res_Tw[0] = Tw[0] - Tw[1];
    res_Tw[N] = Tw[N] - Tw[N - 1];
	}

} // end adsim namespace


// Register the derived type with the base type
CEREAL_REGISTER_TYPE(adsim::WallTemperatureBoundaryConditionFixed)
CEREAL_REGISTER_TYPE(adsim::WallTemperatureBoundaryConditionZeroFlux)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IWallTemperatureBoundaryCondition, adsim::WallTemperatureBoundaryConditionFixed)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IWallTemperatureBoundaryCondition, adsim::WallTemperatureBoundaryConditionZeroFlux)
CEREAL_CLASS_VERSION(adsim::WallTemperatureBoundaryConditionFixed, 1);
CEREAL_CLASS_VERSION(adsim::WallTemperatureBoundaryConditionZeroFlux, 1);
