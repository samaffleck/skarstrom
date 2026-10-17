/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <string>

// SUNDIALS includes
#include <ida/ida.h>


namespace adsim {

  class UserData;
  class Equations;
  class Bed;
  using Component = std::string; // Forward declaration - Component is std::string from fluid.h
  using AdsorbentID = std::size_t; // mirrors porous_media.h; identical typedef redeclaration is legal

  void UpdateConcentration(UserData* udata);
  void UpdateVelocity(UserData* udata);

  // Momentum coefficients of the Ergun face `face` (between cells face and
  // face + 1), such that the face velocity u solves
  //   P[face] - P[face + 1] = A u + B |u| u.
  // Returns false, leaving A and B untouched, for a Darcy (Kozeny-Carman) face
  // (no inertial resistance on either side), whose residual keeps its legacy
  // closed form. Interior faces take the two half-cells in series; a face
  // with a void cell on one side takes the bed-side cell over its full width
  // plus, on a void-bed (not ghost-void) face, the void region's resistance
  // coefficient, read from `bed`. `bed` may be null only for the two
  // boundary faces (0 and n_cells - 2), which never carry that resistance.
  bool ErgunFaceCoefficients(const Equations& state, const Bed* bed, int face, double& A, double& B);
  void UpdateTemperature(UserData* udata);
  void UpdateFluidTemperature(UserData* udata);
  void UpdateSolidTemperature(UserData* udata);
  void UpdateWallTemperature(UserData* udata);
  void UpdateParticleHeatTransferCoefficient(UserData* udata);
  void UpdatePressure(UserData* udata);
  void UpdateMolarConcentration(UserData* udata, Component component);
  void UpdateAdsorbedConcentration(UserData* udata, AdsorbentID ads_id, Component component);
  void UpdateEquilibriumAdsorbedConcentration(UserData* udata, AdsorbentID ads_id, Component component);
  void UpdateMoleFraction(UserData* udata, Component component);
  void UpdateMoleFractionFaceUpwind(UserData* udata, Component component);
  void UpdateMoleFractionFaceWENO(UserData* udata, Component component);
  void UpdateMoleFractionFaceVanLeer(UserData* udata, Component component);
  // Dispatches to the scheme selected in udata->settings.advection_scheme_.
  void UpdateMoleFractionFace(UserData* udata, Component component);
  void UpdateMassTransferCoefficient(UserData* udata, AdsorbentID ads_id, Component component);
  void UpdateMolecularDiffusionCoefficient(UserData* udata, Component component);
  void UpdateTotalMoleFraction(UserData* udata);
  void UpdateFluidDensity(UserData* udata);
  void UpdateViscosity(UserData* udata);
  void UpdateHeatCapacity(UserData* udata);
  void UpdateSherwoodNumber(UserData* udata);
  void UpdateReynoldsNumber(UserData* udata);
  void UpdateSchmidtNumber(UserData* udata);
  void UpdateFluidThermalConductivity(UserData* udata);
  void UpdateAxialThermalConductivity(UserData* udata);
  void UpdateInternalHeatTransferCoefficient(UserData* udata);
  void UpdateExternalHeatTransferCoefficient(UserData* udata);
  void UpdateMolarMass(UserData* udata);
  void UpdateVolumetricFlowRate(UserData* udata);
  void UpdateMolarFlowRate(UserData* udata);
  void UpdateMassSource(UserData* udata);
  void UpdateEnergySource(UserData* udata);
  void UpdateApproachToFluidisation(UserData* udata);
  void UpdatePressureGradient(UserData* udata);
  void UpdateComponentMassSource(UserData* udata, Component component);
  void UpdateComponentEnergySource(UserData* udata, Component component);
  void UpdateAxialDispersionCoefficient(UserData* udata);

} // end adsim namespace
