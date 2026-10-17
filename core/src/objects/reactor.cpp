/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/objects/reactor.h"


namespace adsim {

  bool Bed::PerformChecks(const AdsorbentDirectory& dir, std::string& errorMsg) const {
    if (layers_.layers_.empty()) {
      errorMsg += "Error: Adsorbent layers are empty. You need to add at least one layer to the column.\n";
      return false;
    }

    // Check that each layer has at least one adsorbent selected
    for (size_t i = 0; i < layers_.layers_.size(); ++i) {
      const auto& layer = layers_.layers_[i];
      if (layer.composition_.empty()) {
        errorMsg += "Error: Layer \"" + layer.name_ + "\" has no adsorbent selected. Each layer must have at least one adsorbent.\n";
        return false;
      }
      // Every composition entry must name an adsorbent that still exists in the
      // library; otherwise the directory lookups during initialisation throw.
      for (const auto& [ads_id, frac] : layer.composition_) {
        if (!dir.Contains(ads_id)) {
          errorMsg += "Error: Layer \"" + layer.name_ + "\" references adsorbent id " + std::to_string(ads_id)
                    + ", which does not exist in the adsorbent library. Remove it from the layer composition or select another adsorbent.\n";
          return false;
        }
      }
    }

    // Check that LTNE mode has a valid particle heat transfer model
    if (equations_.thermal_equilibrium_mode_ == ThermalEquilibriumMode::NON_LOCAL) {
      if (!equations_.particle_heat_transfer_model_) {
        errorMsg += "Error: Local thermal non-equilibrium (LTNE) mode requires a particle heat transfer model.\n";
        return false;
      }
    }

    // A user-set isothermal temperature must be a physical (positive) absolute
    // temperature; the initial-condition temperature is checked with the ICs.
    if (equations_.IsIsothermal() && !equations_.isothermal_use_initial_temperature_
        && !(equations_.isothermal_temperature_.Get() > 0.0)) {
      errorMsg += "Error: The isothermal temperature must be greater than 0 K.\n";
      return false;
    }

    return true;
  }

  void Bed::Initialise(const AdsorbentDirectory& dir, Fluid& fluid) {
    auto& layers = layers_.layers_;

    int start_index = 1 + (bottom_void_.enabled_ ? 1 : 0);
    for (auto& layer : layers) {
      int end_index = start_index + layer.domain_.n_cells - 1;
      layer.Initialise(dir, fluid, start_index, end_index);
      start_index = end_index + 1;
    }

    equations_.binary_diffusion_model_.Initialise(GetNumberOfCells(), fluid.components_);
  }

  int Bed::GetNumberOfCells() const {
    auto& layers = layers_.layers_;
    int number_of_cells = 2;
    for (const auto& layer : layers) {
      number_of_cells += layer.domain_.n_cells;
    }
    if (bottom_void_.enabled_) ++number_of_cells;
    if (top_void_.enabled_) ++number_of_cells;
    return number_of_cells;
  }

  void Bed::RemoveComponent(const std::string& component, const AdsorbentDirectory& dir) {
    // Visit every adsorbent in the library, not only those currently in a layer:
    // one assigned to a layer later must not carry params for the removed component.
    for (auto& [ads_id, adsorbent] : dir.GetAll()) {
      const_cast<Adsorbent&>(adsorbent).adsorption_source_.RemoveComponent(component);
    }
  }

} // end adsim namespace


CEREAL_CLASS_VERSION(adsim::Bed, 4);
