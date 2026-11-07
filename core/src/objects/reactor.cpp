/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/objects/reactor.h"
#include "adsim/equations/dispersion_model.h"


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

    // Structured layers: a geometry that leaves an open channel (porous_media.h).
    bool geometry_ok = true;
    for (const auto& layer : layers_.layers_) {
      for (const auto& e : CheckLayerGeometry(layer)) {
        errorMsg += e + "\n";
        geometry_ok = false;
      }
    }
    if (!geometry_ok) return false;

    // Closures that belong to one kind of layer. On a structured layer:
    // Ruthven's packed-bed dispersion (Taylor-Aris replaces it), Ruthven's and
    // Yagi-Kunii's packed-bed heat dispersion (channel conduction replaces
    // them) and Dixon's and Leva's packed-bed wall coefficients (a constant
    // h_i replaces them). On a packed bed: Taylor-Aris dispersion and channel
    // conduction (they need channels).
    for (const auto& layer : layers_.layers_) {
      const std::string kind = LayerKindList[static_cast<int>(layer.kind_)];
      const auto* heat = layer.equations_.heat_dispersion_model_.get();
      const auto* wall = layer.equations_.internal_heat_transfer_model_.get();
      if (layer.IsStructured()) {
        const char* packed_closure = nullptr;
        const char* replacement = nullptr;
        if (dynamic_cast<const HeatDispersionModelRuthven*>(heat)) {
          packed_closure = "Ruthven's packed-bed heat dispersion";
          replacement = "use channel conduction (channel) or a constant axial conductivity";
        } else if (dynamic_cast<const HeatDispersionModelYagiandKunii*>(heat)) {
          packed_closure = "Yagi and Kunii's packed-bed heat dispersion";
          replacement = "use channel conduction (channel) or a constant axial conductivity";
        } else if (dynamic_cast<const InternalHeatTransferModel_Dixon*>(wall)) {
          packed_closure = "Dixon's packed-bed wall heat transfer coefficient";
          replacement = "use a constant coefficient";
        } else if (dynamic_cast<const InternalHeatTransferModel_Leva*>(wall)) {
          packed_closure = "Leva's packed-bed wall heat transfer coefficient";
          replacement = "use a constant coefficient";
        }
        if (packed_closure) {
          errorMsg += "Error: Layer \"" + layer.name_ + "\" is a " + kind + " but uses " + packed_closure +
                      ", which reads the particle diameter; " + replacement + ".\n";
          return false;
        }
      } else if (dynamic_cast<const HeatDispersionModelChannel*>(heat)) {
        errorMsg += "Error: Layer \"" + layer.name_ + "\" is a packed bed but uses channel heat conduction, which is"
                    " for monoliths and laminates; use Ruthven or Yagi-Kunii heat dispersion or a constant axial"
                    " conductivity.\n";
        return false;
      }
    }
    for (const auto& layer : layers_.layers_) {
      const auto* disp = layer.equations_.dispersion_model_.get();
      if (layer.IsStructured() && dynamic_cast<const DispersionModelRuthven*>(disp)) {
        errorMsg += "Error: Layer \"" + layer.name_ + "\" is a " + LayerKindList[static_cast<int>(layer.kind_)] +
                    " but uses Ruthven's packed-bed dispersion; use Taylor-Aris dispersion (taylor_aris) or a constant coefficient.\n";
        return false;
      }
      if (!layer.IsStructured() && dynamic_cast<const DispersionModelTaylorAris*>(disp)) {
        errorMsg += "Error: Layer \"" + layer.name_ + "\" is a packed bed but uses Taylor-Aris dispersion, which is for"
                    " channels; use Ruthven dispersion or a constant coefficient.\n";
        return false;
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
