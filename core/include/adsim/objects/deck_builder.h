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
//   rate model           "ldf", "ldf2", "dual_rate_ldf"
//   rate coefficient     "constant", "macropore", "macropore_knudsen", "haynes", "arrhenius",
//                        "macropore_capacity", "series_resistance"
//   reaction rate        "power_law", "reversible_power_law", "lhhw", "expression"
//                        (ReactionRateKindString, equations/reaction_kinds.h)
//   reaction basis       "partial_pressure", "concentration"
//   effectiveness        "constant", "thiele"
//   momentum             "kozeny_carman", "ergun", "blake_kozeny"
//   film correlation     "wakao_funazkri", "gunn", "ranz_marshall"
//   dispersion           "constant", "ruthven"
//   heat dispersion      "constant", "ruthven", "yagi_kunii"
//   internal heat transfer "constant", "dixon", "leva"
//   external heat transfer "constant", "holman"
//   particle heat transfer "constant", "film"
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

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "adsim/objects/porous_media.h"
#include "adsim/objects/user_data.h"

namespace adsim::build {

// ---- kind names -------------------------------------------------------------
std::vector<std::string> IsothermKinds();
const std::vector<std::string>& RateModelKinds();
// The coefficient kinds offered: "series_resistance" (the generic resistance
// model). The former kinds are its presets (RateCoefficientPresets), still
// accepted by SetRateModel.
const std::vector<std::string>& RateCoefficientKinds();
const std::vector<std::string>& RateCoefficientPresets();
const std::vector<std::string>& MixingRuleKinds();
const std::vector<std::string>& MomentumKinds();
const std::vector<std::string>& FilmCorrelationKinds(); // wakao_funazkri, gunn, ranz_marshall
const std::vector<std::string>& LayerKinds();          // packed, monolith, laminate
const std::vector<std::string>& ChannelShapeKinds();   // square, circular
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
const std::vector<std::string>& ReservoirKinds();  // source, sink
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

// How the adsorbent combines its components' isotherms into mixture loadings:
// "none" (each isotherm on its own; the extended kinds couple components
// themselves) or "iast" (Ideal Adsorbed Solution Theory over the pure
// isotherms, iast.h).
void SetMixingRule(UserData& data, AdsorbentID adsorbent, const std::string& kind);

// The rate model of one component: `model` with a `coefficient` model, the
// kind "series_resistance" (the resistance model at its defaults: the lumped
// term only) or a preset tag (RateCoefficientPresets), which sets the switches.
void SetRateModel(UserData& data, AdsorbentID adsorbent, const std::string& component,
                  const std::string& model, const std::string& coefficient);

// ---- reactions (objects/reaction.h) -------------------------------------------------
// A reaction on a solid `solid` (an adsorbent ID) with stoichiometry `nu`
// (component -> coefficient, < 0 reactant, > 0 product; exact zeros are
// dropped) and a rate law of `kind`. Every component must be in the deck.
// Orders default to -nu for reactants and 0 for products. Returns its ID.
// Kind, stoichiometry and the LHHW adsorption terms are structure (these
// functions); every value goes through the registry
// (adsorbent[ID].reaction[ID].*).
ReactionID AddReaction(UserData& data, AdsorbentID solid, const std::string& name,
                       const std::map<std::string, double>& nu, const std::string& kind = "power_law");
// A new rate law of `kind`; parameters the new layout shares keep their values.
void SetReactionKind(UserData& data, AdsorbentID solid, ReactionID reaction, const std::string& kind);
void SetReactionStoichiometry(UserData& data, AdsorbentID solid, ReactionID reaction,
                              const std::map<std::string, double>& nu);
void SetReactionBasis(UserData& data, AdsorbentID solid, ReactionID reaction, const std::string& basis);
// The effectiveness model: "constant" (eta through the registry) or "thiele"
// with its key reactant.
void SetEffectiveness(UserData& data, AdsorbentID solid, ReactionID reaction, const std::string& mode,
                      const std::string& key_component = "");
// A new LHHW denominator term (lhhw only); returns its index j (from 1).
size_t AddAdsorptionTerm(UserData& data, AdsorbentID solid, ReactionID reaction);
void RemoveAdsorptionTerm(UserData& data, AdsorbentID solid, ReactionID reaction, size_t j);
void RemoveReaction(UserData& data, AdsorbentID solid, ReactionID reaction);
std::vector<std::string> ReactionRateKinds();
const std::vector<std::string>& ReactionBasisKinds();
const std::vector<std::string>& EffectivenessKinds();

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
// The packed layer's gas-particle film correlation (FilmCorrelationKinds;
// film_correlation.h): the Nu of the film particle heat transfer and the Sh
// of the Haynes rate coefficient. Kept, unused, on a structured layer.
void SetFilmCorrelation(UserData& data, size_t layer, const std::string& kind);
// What the layer is packed with (LayerKinds: "packed", "monolith",
// "laminate"). The composition, domain, name and every other field are kept,
// including the structured geometry and the packed-bed fields, so switching
// back restores them. A structured kind uses laminar channel flow instead of
// the momentum model. Switching a packed layer to a structured kind installs
// the channel closures: Taylor-Aris dispersion on the layer (keeping a
// Ruthven model's molecular diffusion coefficient), channel heat conduction,
// a constant wall coefficient h_i (replacing Dixon or Leva; a constant model
// keeps its value) and the bed's film gas-solid heat transfer. Switching back
// restores the packed defaults: constant dispersion and heat-dispersion
// models instead of Taylor-Aris and channel conduction, and the constant
// gas-solid model once no structured layer remains.
void SetLayerKind(UserData& data, size_t layer, const std::string& kind);
// The monolith channel cross-section (ChannelShapeKinds: "square",
// "circular"); kept, unused, for other kinds.
void SetChannelShape(UserData& data, size_t layer, const std::string& shape);
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
// A named source ("source": inlets attach to it) or sink ("sink": outlets
// attach to it); see adsim/boundary_conditions/reservoir.h. Returns its ID.
// Its values are registry entries: reservoir[ID].temperature, .pressure,
// .y[Component].
size_t AddReservoir(UserData& data, const std::string& kind, const std::string& name);
// Attach an endpoint to a reservoir, so that its temperature, mole fractions
// and pressure target are the reservoir's; reservoir 0 detaches it.
void AttachReservoir(UserData& data, size_t endpoint, size_t reservoir);
// An energy source (heater or cooler) acting on a layer during a step.
void SetEnergySource(UserData& data, size_t step, size_t layer, bool enabled);
// Track a component's KPIs with these feed and product endpoints. Recycle
// endpoints carry product back into the bed (a purge or pressurisation fed
// from the product): recovery and productivity count product minus recycle.
void SetKPI(UserData& data, const std::string& component, const std::vector<size_t>& feed_endpoints,
            const std::vector<size_t>& product_endpoints, const std::vector<size_t>& recycle_endpoints = {});
// Reaction KPIs (KPIRecord::conversion, selectivity, reaction_yield,
// space_time_yield, peak_temperature): the key reactant, the products and the
// space-time-yield basis ("catalyst_mass" or "bed_volume"). An empty key
// reactant clears them.
void SetReactionKPI(UserData& data, const std::string& key_reactant, const std::vector<std::string>& products,
                    const std::string& basis = "catalyst_mass");

// ---- user variables -----------------------------------------------------------------
// A named value other values' expressions can use ("2 * t_feed"). `formula`
// (may be empty) defines it from other variables; otherwise its value is set
// through the registry (variable[name]). Replaces a variable of the same name.
void AddVariable(UserData& data, const std::string& name, const std::string& category,
                 const std::string& formula, const std::string& description);

}  // namespace adsim::build
