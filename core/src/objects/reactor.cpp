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

    // Each layer's composition (CheckLayerComposition, porous_media.h): every
    // error of every layer is reported; warnings are left to the run-start
    // deck warnings (LayerCompositionWarnings).
    bool compositions_ok = true;
    for (const auto& layer : layers_.layers_) {
      for (const auto& issue : CheckLayerComposition(layer, dir)) {
        if (issue.severity != LayerCompositionIssue::Severity::Error) continue;
        errorMsg += issue.message + "\n";
        compositions_ok = false;
      }
    }
    if (!compositions_ok) return false;

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
