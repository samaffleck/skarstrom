/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// The gas state of a dynamic tank as the bed's boundary conditions see it
// during a run, and the reference from a bed endpoint to a tank. Kept apart
// from tank.h so the boundary conditions can name a tank without including
// the tank's streams and serialisation.

#include <cstddef>
#include <string>
#include <unordered_map>

#include <cereal/access.hpp>
#include <cereal/cereal.hpp>

namespace adsim {

  // Live state of one tank, refreshed from the DAE unknowns on every residual
  // evaluation (Equations::UnpackBanded). Pressure from the ideal gas law,
  // P = N R T / V, and mole fractions n_i / N.
  struct TankLiveState {
    double pressure = 0.0;      // [Pa]
    double temperature = 0.0;   // [K]
    double total_moles = 0.0;   // [mol]
    std::unordered_map<std::string, double> mole_fraction{};

    double MoleFraction(const std::string& component) const {
      const auto it = mole_fraction.find(component);
      return it == mole_fraction.end() ? 0.0 : it->second;
    }
  };

  // A bed endpoint's port to a tank (IBoundaryStream::tank_port_,
  // PressureBoundaryCondition::tank_port_). tank_id_ = 0: no tank. The
  // endpoint's own direction and bed end say which way and where gas flows:
  // an INLET draws gas from the tank into the bed, an OUTLET sends the bed's
  // gas into the tank. Mutually exclusive with a connection and a buffer-tank
  // source (Cycle::PerformChecks).
  struct TankPort {
    size_t tank_id_ = 0;
    // Transient (never serialised): the tank's live state during a run, bound
    // by Equations::InitialiseTanks. nullptr outside a run.
    const TankLiveState* live_ = nullptr;

    bool IsSet() const noexcept { return tank_id_ != 0; }
    bool IsBound() const noexcept { return live_ != nullptr; }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/) {
      ar(cereal::make_nvp("Tank_ID", tank_id_));
    }
  };

} // end adsim namespace

CEREAL_CLASS_VERSION(adsim::TankPort, 1);
