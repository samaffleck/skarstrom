/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Building a deck from scratch.
//
// These functions create the STRUCTURE of a deck — components, adsorbents
// with their isotherm and rate models, layers with their transport models,
// the column, steps and their boundaries, streams, pumps, connections, energy
// sources and the KPI configuration — choosing each model by a kind name
// (lowercase snake_case, listed below). Every object starts with the model's
// default values; set the values afterwards through the parameter registry
// (adsim/parameters/parameter_registry.h), which lists every one of them by
// path. The Python package's builder (skarstrom.Simulation.new()) is a thin
// layer over these functions.
//
// Kinds (in enum order):
//   isotherm             IsothermKindString names (isotherm_kinds.h), e.g.
//                        "extended_dual_site_langmuir_c", "sips", "inert"
//   rate model           "ldf", "ldf2"
//   rate coefficient     "constant", "macropore", "macropore_knudsen", "haynes", "arrhenius"
//   momentum             "kozeny_carman", "ergun", "blake_kozeny"
//   dispersion           "constant", "ruthven"
//   heat dispersion      "constant", "ruthven", "yagi_kunii"
//   internal heat transfer "constant", "dixon", "leva"
//   external heat transfer "constant", "holman"
//   particle heat transfer "constant"
//   thermal equilibrium  "lte", "ltne"
//   bed shape            "cylinder", "square", "rectangle"
//   wall material        "custom", "ss316l", "ss304", "ss304l", "carbon_steel_a106",
//                        "carbon_steel_a53", "hastelloy_c276", "inconel_625",
//                        "aluminium_6061", "titanium_gr2"
//   wall boundary        "zero_flux", "fixed"
//   profile              "constant", "linear", "exponential", "log_linear", "double_exponential"
//   flow                 "velocity", "volumetric_flow", "molar_flow", "mass_flow"
//   pump                 "none", "compressor", "vacuum_pump"
//   pump energy          "isothermal", "adiabatic"
//   pump efficiency      "constant", "vacuum_pump", "polynomial"
//   location             "bottom" (feed end, z = 0), "top"
//   direction            "inlet", "outlet"
//
// Every function throws std::invalid_argument for an unknown kind (the message
// lists the valid ones) or a missing object (adsorbent, layer, step or
// endpoint ID). Steps are addressed by Step::id_, endpoints (flow streams and
// pressure boundaries) by their endpoint_id_; both are assigned on creation,
// in creation order, as Cycle::EnsureIds does. The functions do not evaluate
// expressions: refresh the parameter registry after a batch of edits.

#include <string>
#include <utility>
#include <vector>

#include "adsim/objects/porous_media.h"
#include "adsim/objects/user_data.h"

namespace adsim::build {

// ---- kind names -------------------------------------------------------------
std::vector<std::string> IsothermKinds();
const std::vector<std::string>& RateModelKinds();
const std::vector<std::string>& RateCoefficientKinds();
const std::vector<std::string>& MomentumKinds();
const std::vector<std::string>& DispersionKinds();
const std::vector<std::string>& HeatDispersionKinds();
const std::vector<std::string>& InternalHeatTransferKinds();
const std::vector<std::string>& ExternalHeatTransferKinds();
const std::vector<std::string>& ParticleHeatTransferKinds();
const std::vector<std::string>& ThermalEquilibriumKinds();
const std::vector<std::string>& BedShapeKinds();
const std::vector<std::string>& WallMaterialKinds();
const std::vector<std::string>& WallBoundaryKinds();
const std::vector<std::string>& ProfileKinds();
const std::vector<std::string>& FlowKinds();
const std::vector<std::string>& PumpKinds();
const std::vector<std::string>& PumpEnergyKinds();
const std::vector<std::string>& PumpEfficiencyKinds();
const std::vector<std::string>& PumpPowerKinds();
// Unit categories of a user variable (UnitCategory, snake_case).
const std::vector<std::string>& VariableCategoryKinds();
std::string VariableCategoryKind(UnitCategory category);

// ---- components and adsorbents ------------------------------------------------
// Add a fluid component (a CoolProp name) with zero mole fraction everywhere
// and default physics on every adsorbent (AddComponentToDeck).
void AddComponent(UserData& data, const std::string& name);

// A new adsorbent with default physics for every component: an inert isotherm,
// an LDF rate model with a constant coefficient, zero heat of adsorption.
AdsorbentID AddAdsorbent(UserData& data, const std::string& name);

// One isotherm of `kind` for `components` (a competitive kind couples them).
// The components leave their previous isotherm; an isotherm left with no
// component is removed.
void SetIsotherm(UserData& data, AdsorbentID adsorbent, const std::string& kind,
                 const std::vector<std::string>& components);

// The rate model of one component: `model` with a `coefficient` model.
void SetRateModel(UserData& data, AdsorbentID adsorbent, const std::string& component,
                  const std::string& model, const std::string& coefficient);

// ---- column ---------------------------------------------------------------------
void SetBedShape(UserData& data, const std::string& kind);
// Sets the material and its density, heat capacity and conductivity.
void SetWallMaterial(UserData& data, const std::string& kind);
void SetWallBoundary(UserData& data, const std::string& kind);
void SetThermalEquilibrium(UserData& data, const std::string& kind);
void SetParticleHeatTransfer(UserData& data, const std::string& kind);

// A new layer at the product end of the bed (no adsorbent until
// SetLayerComposition). Returns its LayerID.
size_t AddLayer(UserData& data, const std::string& name);
// The layer's adsorbents and their mass fractions (replaces the composition).
void SetLayerComposition(UserData& data, size_t layer,
                         const std::vector<std::pair<AdsorbentID, double>>& composition);
void SetMomentum(UserData& data, size_t layer, const std::string& kind);
void SetDispersion(UserData& data, size_t layer, const std::string& kind);
void SetHeatDispersion(UserData& data, size_t layer, const std::string& kind);
void SetInternalHeatTransfer(UserData& data, size_t layer, const std::string& kind);
void SetExternalHeatTransfer(UserData& data, size_t layer, const std::string& kind);

// ---- cycle ------------------------------------------------------------------------
// A new step at the end of the cycle with both ends closed. Returns its ID.
size_t AddStep(UserData& data, const std::string& name);
// Close one end of a step.
void SetClosed(UserData& data, size_t step, const std::string& location);
// Make one end a pressure boundary following a `profile`; returns its endpoint ID.
size_t SetPressureBoundary(UserData& data, size_t step, const std::string& location,
                           const std::string& direction, const std::string& profile);
// Make one end a flow boundary with no streams yet (AddStream / AddValveStream).
void SetFlowBoundary(UserData& data, size_t step, const std::string& location, const std::string& direction);
// A stream with a set `flow` following a `profile`; returns its endpoint ID.
size_t AddStream(UserData& data, size_t step, const std::string& location, const std::string& name,
                 const std::string& flow, const std::string& profile);
// A stream through a valve against a pressure following a `profile`.
size_t AddValveStream(UserData& data, size_t step, const std::string& location, const std::string& name,
                      const std::string& profile);
// The inlet temperature of a stream or pressure boundary: "constant" (no
// ramp), "linear" or "exponential".
void SetTemperatureProfile(UserData& data, size_t endpoint, const std::string& profile);
// A flow stream driven by a pump of given shaft power: a compressor on an
// inlet, a vacuum pump on an outlet. `power` is "constant" (a fixed power) or
// "polynomial" (power a fifth-order polynomial of the boundary pressure);
// `energy` and `efficiency` as SetPump. Returns its endpoint ID.
size_t AddPumpStream(UserData& data, size_t step, const std::string& location, const std::string& name,
                     const std::string& power, const std::string& energy, const std::string& efficiency);
// Feed an inlet (a stream or a pressure boundary) from a buffer tank: it
// draws a percentage (endpoint[N].buffer_tank.percentage, 100 % by default)
// of the moles endpoint `source` produced in its step's last run, spread over
// this step, at their average temperature and composition. Until the source
// has run, the inlet uses its own values. `source` 0 removes the tank. An
// endpoint with a connection cannot also have a buffer tank.
void SetBufferTank(UserData& data, size_t endpoint, size_t source);
// The pump of a flow stream or pressure boundary ("none" removes it).
void SetPump(UserData& data, size_t endpoint, const std::string& kind, const std::string& energy,
             const std::string& efficiency);
// Connect two endpoints (one supplies the other, e.g. light-product
// pressurisation or a purge drawn from another step's product).
void Connect(UserData& data, size_t endpoint_a, size_t endpoint_b);
// An energy source (heater or cooler) acting on a layer during a step.
void SetEnergySource(UserData& data, size_t step, size_t layer, bool enabled);
// Track a component's KPIs with these feed and product endpoints. Recycle
// endpoints carry product back into the bed (a purge or pressurisation fed
// from the product): recovery and productivity count product minus recycle.
void SetKPI(UserData& data, const std::string& component, const std::vector<size_t>& feed_endpoints,
            const std::vector<size_t>& product_endpoints, const std::vector<size_t>& recycle_endpoints = {});

// ---- user variables -----------------------------------------------------------------
// A named value other values' expressions can use ("2 * t_feed"). `formula`
// (may be empty) defines it from other variables; otherwise its value is set
// through the registry (variable[name]). Replaces a variable of the same name.
void AddVariable(UserData& data, const std::string& name, const std::string& category,
                 const std::string& formula, const std::string& description);

}  // namespace adsim::build
