/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Which advisory warnings a run reports. None of them stops a run; a switched-
// off warning is neither logged nor listed in the operating-limit report.
// Saved with the deck (IDASettings, "Warnings"); registry paths
// settings.warn_<name>.

#include <cstdint>

#include <cereal/access.hpp>
#include <cereal/cereal.hpp>

namespace adsim {

  struct WarningSettings {
    // Start of a run (deck_warnings.h).
    bool pressure_step_change = true;  // a constant pressure boundary far from the bed pressure
    bool wall_thickness = true;        // the wall is thinner than its pressure design requires
    // End of a run, over the last logged cycle (operating_limits.h).
    bool design_pressure = true;       // pressure above the wall's design pressure
    bool vacuum = true;                // below atmospheric without vacuum service
    bool design_temperature = true;    // wall temperature above the design temperature
    bool fluidisation = true;          // approach to fluidisation above 1

  private:
    friend class ::cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/) {
      ar(
        cereal::make_nvp("Pressure_Step_Change", pressure_step_change),
        cereal::make_nvp("Wall_Thickness", wall_thickness),
        cereal::make_nvp("Design_Pressure", design_pressure),
        cereal::make_nvp("Vacuum", vacuum),
        cereal::make_nvp("Design_Temperature", design_temperature),
        cereal::make_nvp("Fluidisation", fluidisation)
      );
    }
  };

} // namespace adsim

CEREAL_CLASS_VERSION(adsim::WarningSettings, 0);
