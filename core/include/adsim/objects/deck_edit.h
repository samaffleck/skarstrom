/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Structural deck edits that keep a UserData consistent.
//
// The model classes expose low-level operations (Fluid::AddComponent,
// AdsorbentDirectory::Add, ...) that each touch one object. A front end (or a
// script) editing a deck also needs the small cascades that keep the rest of
// the deck in step: a new component needs a mole fraction in every map and
// default physics on every adsorbent; a removed adsorbent must not leave an
// empty layer; disconnecting an endpoint must also disconnect its partner.
// These helpers are those cascades. They never evaluate expressions or run the
// deck checks: call UserData::EvaluateAllExpressions and PerformChecks after a
// batch of edits.
//
// Together with the parameter registry (values, by path) they are the whole
// editing API of a deck: a front end that changes structure only through
// these functions never builds model objects itself. Functions that can be
// refused return false with a message in `err` and leave the deck unchanged.
// Objects are addressed by their stable IDs: steps by Step::id_, endpoints
// (flow streams and pressure boundaries) by endpoint_id_, layers by LayerID.

#include <cstddef>
#include <string>
#include <vector>

#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_value.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/connection.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/equations/mass_transfer_model.h"
#include "adsim/objects/compressor.h"
#include "adsim/objects/model_kinds.h"
#include "adsim/objects/porous_media.h"
#include "adsim/objects/user_data.h"

namespace adsim {

// ---- components ---------------------------------------------------------------

// True if `name` is a fluid component the engine knows (Components::componentNames).
bool IsAvailableComponent(const std::string& name);

// Add (add=true) or erase `component` in every mole-fraction map of the deck:
// the initial conditions and each step's pressure boundaries and flow streams.
// An added entry is 0.
void ReconcileMoleFractions(UserData& data, const std::string& component, bool add);

// Add a component to the fluid, seed a 0 mole fraction everywhere, a 0
// adsorbed-phase heat capacity and default physics on every adsorbent (Inert
// isotherm, LDF rate model with a constant coefficient, 0 heat of adsorption;
// existing components are untouched).
// False + err for an unknown component name. Adding an existing component is a
// no-op that returns true.
bool AddComponentToDeck(UserData& data, const std::string& name, std::string& err);

// Remove a component from the fluid, from every adsorbent, from every
// mole-fraction map and from the adsorbed-phase heat capacities. A deck may
// end up with no components (it then fails the deck checks until one is added).
// Refused (false + err, nothing changed) when a reaction of any solid uses the
// component (a nonzero coefficient, order or LHHW adsorption exponent; see
// Reaction::UsesComponent), e.g. "Water is used by reactions SMR, WGS (on Ni
// catalyst); remove them first". Exact-zero entries of the component are
// dropped from the reactions. The reaction KPIs drop it too: a removed key
// reactant clears them, a removed product leaves their product list.
bool RemoveComponentFromDeck(UserData& data, const std::string& name, std::string& err);

// As above for a deck whose reactions do not use the component (asserted;
// in a release build a refused removal leaves the deck unchanged).
void RemoveComponentFromDeck(UserData& data, const std::string& name);

// The refusal message of RemoveComponentFromDeck, or "" when no reaction
// uses `name`.
std::string ComponentInUseByReactions(const UserData& data, const std::string& name);

// ---- layers and adsorbents -----------------------------------------------------

// Give an empty layer the lowest-id adsorbent of the directory at mass
// fraction 1.0 (no-op when the directory is empty or the layer already holds
// a solid).
void SeedLayerComposition(UserData& data, MixedAdsorbentLayer& layer);

// Blend edits. Each leaves the layer's mass fractions summing to 1, so a
// front end that edits blends only through them never builds a deck the
// composition check (CheckLayerComposition) refuses. The other solids are
// rescaled in proportion to their fractions (equal shares when they are all
// 0). False + err, with nothing changed, for an unknown solid or a fraction
// outside the stated range.
//
// A layer whose fractions include a formula (for example "x" and "1 - x" over
// a deck variable) is refused by all three with nothing changed: rescaling
// cannot keep a formula's meaning, and rescaling only the numbers would not
// keep the sum at 1 for every value of the variable. Edit such fractions
// directly (the registry path layer[ID].mass_fraction[ID], or
// MixedAdsorbentLayer::composition_) or make them numbers first.

// Add solid `id` to the layer at `mass_fraction` (0 < w <= 1), scaling the
// others by (1 - w). The first solid of an empty layer gets 1. A solid the
// layer already holds is set as by SetLayerMassFraction.
bool AddSolidToLayer(UserData& data, MixedAdsorbentLayer& layer, AdsorbentID id, double mass_fraction,
                     std::string& err);

// Set the mass fraction of a solid the layer holds (0 <= w <= 1) and rescale
// the others to 1 - w. A layer whose only solid this is accepts only 1.
bool SetLayerMassFraction(UserData& data, MixedAdsorbentLayer& layer, AdsorbentID id, double mass_fraction,
                          std::string& err);

// Take a solid out of the layer and rescale the rest to sum to 1. Refused for
// the layer's only solid (a layer needs one: add another first).
bool RemoveSolidFromLayer(UserData& data, MixedAdsorbentLayer& layer, AdsorbentID id, std::string& err);

// SeedLayerComposition for every layer whose composition is empty.
void SeedEmptyLayerCompositions(UserData& data);

// Add a new adsorbent with default physics for every component of the fluid,
// and assign it to any layer that has no adsorbent yet. Returns its id.
AdsorbentID AddAdsorbentToDeck(UserData& data, const std::string& name = "New Adsorbent");

// Deserialise a single-adsorbent JSON blob (the format of
// data/adsorbents/<key>.json, written by SaveAdsorbentToJSON) into `out`.
// False + err on a parse failure; never throws.
bool LoadAdsorbentFromString(const std::string& json, Adsorbent& out, std::string& err);

// Add an adsorbent from a single-adsorbent JSON blob. Components of the deck
// that the blob does not cover get default physics; parameters for components
// the deck's fluid lacks are kept but unused. Its reactions are kept too; one
// that uses a component the deck lacks is inactive and reported by the deck
// checks (ReactionChecks) once a layer uses the solid. On success `id` is the new
// adsorbent; on failure nothing is added and err says why.
bool AddAdsorbentFromJson(UserData& data, const std::string& label, const std::string& json,
                          AdsorbentID& id, std::string& err);

// Add an entry of the curated adsorbent library (data/adsorbents/library.json
// in `library_dir`) by its key, or a material's recommended entry by the
// material id (e.g. "zeolite-13x") or the alias of a removed entry. The
// adsorbent takes the name stored in the entry's file.
bool AddLibraryAdsorbent(UserData& data, const std::string& library_dir, const std::string& key,
                         AdsorbentID& id, std::string& err);

// Remove an adsorbent from the directory and every layer, then re-seed any
// layer left empty with another adsorbent (if one remains). Its reactions go
// with it.
// Refused (false + err, nothing changed) while the solid is part of a blend,
// i.e. a layer whose composition also holds another solid, e.g. "Zeolite 13X
// is part of the blend of layer \"Bottom\" (mass fraction 0.7); remove it
// from the layer's composition first". Taking one solid out of a blend
// changes what the layer is (and its fractions would no longer sum to 1), so
// the user decides the new composition; the engine never renormalises. A
// layer that holds only this solid has no blend to change: it is re-seeded.
bool RemoveAdsorbentFromDeck(UserData& data, AdsorbentID id, std::string& err);

// As above for a solid that no blend uses (asserted; in a release build a
// refused removal leaves the deck unchanged).
void RemoveAdsorbentFromDeck(UserData& data, AdsorbentID id);

// The refusal message of RemoveAdsorbentFromDeck, or "" when no blend uses
// the solid (or the id is unknown).
std::string AdsorbentInUseByMixtures(const UserData& data, AdsorbentID id);

// ---- connections ----------------------------------------------------------------

// Visit every connectable endpoint (each flow stream and each pressure
// boundary) of every step.
template <class Fn>
void ForEachEndpoint(UserData& data, Fn fn) {
  for (auto& step : data.cycle.step_list_) {
    for (auto* bc : {step.top_boundary_condition_.get(), step.bottom_boundary_condition_.get()}) {
      if (auto* pbc = dynamic_cast<PressureBoundaryCondition*>(bc)) {
        fn(MakeEndpoint(step, *pbc));
      } else if (auto* fbc = dynamic_cast<FlowBoundaryCondition*>(bc)) {
        for (auto& s : fbc->streams_)
          if (s) fn(MakeEndpoint(step, *s));
      }
    }
  }
}

// Disconnect `conn` together with its partner (every endpoint sharing the same
// non-zero pair id). No-op if not connected.
void DisconnectPair(UserData& data, Connection& conn);

// Disconnect every connection carried by a boundary (a pressure boundary, or
// each stream of a flow boundary) together with its partners. Call before
// replacing or erasing a boundary so no partner keeps a dangling pair.
void DisconnectBoundary(UserData& data, IBoundaryCondition* bc);

// Connect two endpoints (one supplies the other: a light-product
// pressurisation, a purge from another step's product) under a pair id no
// other connection of the deck uses. Refused when either is already connected
// or has a buffer tank, when both are pressure boundaries, or when they are
// the same endpoint.
bool ConnectEndpoints(UserData& data, size_t endpoint_a, size_t endpoint_b, std::string& err);

// ---- identity and repair ------------------------------------------------------------

// Stable IDs by first appearance: every step and endpoint whose ID is 0 or
// repeated gets the next free one (Cycle::EnsureIds), likewise every layer
// (and the layer counter is kept ahead of them), and the connection pair ids
// are renumbered 1..N in the order their first endpoint appears (steps in
// order, bottom before top, streams in order). Two decks with the same
// structure therefore carry the same IDs. Idempotent.
void EnsureIds(UserData& data);

// Re-link every connection from its pair id: a pair id held by exactly two
// endpoints (not both pressure boundaries) is reconnected to point at each
// other's stream data; any other connected endpoint (its partner is gone, or
// three share the id) is disconnected. Returns how many endpoints were
// disconnected. Call after editing step_list_ by hand.
size_t RepairConnections(UserData& data);

// Deactivate every buffer tank whose source endpoint no longer exists (a
// removed step or stream). Returns how many were deactivated. The step and
// stream removals here call it.
size_t PruneBufferTanks(UserData& data);

// ---- KPI configuration ------------------------------------------------------------

// One KPI entry per fluid component, in fluid order (new ones untracked);
// entries of removed components are dropped, and the reaction KPIs forget
// them (KPIConfig::SyncToComponents).
void SyncKpiConfigToComponents(UserData& data);
// Drop KPI endpoint references that no longer name an endpoint of the cycle
// (KPIConfig::Prune).
void PruneKpiConfig(UserData& data);

// ---- filling gaps --------------------------------------------------------------------

// Everything run-time initialisation (UserData::Initialise: Fluid::Initialise,
// Bed::Initialise, AdsorptionSource::Initialise) would otherwise fill in, done
// on the deck itself so an editor sees every entry: a mole fraction per
// component in the initial conditions and every endpoint (missing entries 0,
// entries of removed components dropped), an adsorbed-phase and a constant
// gas heat capacity per component, default physics for every component on
// every adsorbent, a composition for every empty layer, a KPI entry per
// component, stable IDs and repaired connections. A deck edited only through
// the cascades below already has all of this. Idempotent.
void FillDeckGaps(UserData& data);

// ---- layers ----------------------------------------------------------------------------

// A unique name from `base`: `base` itself if no other name in `taken` equals
// it, otherwise "<base> 2", "<base> 3", ...
std::string UniqueName(const std::vector<std::string>& taken, const std::string& base);

// Append a layer (MixedAdsorbentLayer defaults) at the product end. An empty
// `name` gives "Layer N" (N = the new layer count); a taken name is made
// unique (UniqueName). Its composition is seeded (SeedLayerComposition).
// Returns its LayerID.
LayerID AddLayer(UserData& data, const std::string& name = "");

// Remove a layer and everything that refers to it: the energy sources of
// every step on that layer. A deck may end up with no layers (it then fails
// the deck checks until one is added).
bool RemoveLayer(UserData& data, LayerID id, std::string& err);

// Move a layer one place towards the product end (offset +1) or the feed end
// (-1). Refused at the end of the bed or for another offset.
bool MoveLayer(UserData& data, LayerID id, int offset, std::string& err);

// ---- steps -------------------------------------------------------------------------------

// Append a step. An empty `name` gives "Step N"; a taken name is made unique.
// Its ends are the engine's default step: a flow inlet at the bottom with one
// velocity stream "Stream 1", a constant-pressure outlet at the top. Every
// mole-fraction map of the new step has an entry per component summing to 1
// (the first component 1, the rest 0; SeedMoleFractions). Returns its step ID.
size_t AddStep(UserData& data, const std::string& name = "");

// Remove a step: its connections are disconnected together with their
// partners, its KPI references dropped, and buffer tanks on other steps that
// drew from it deactivated.
bool RemoveStep(UserData& data, size_t step_id, std::string& err);

// Copy a step (UserData::DuplicateStep: inserted after it, named
// "<name> (copy)", fully disconnected, no KPI roles); `new_step_id` is the
// copy's ID.
bool DuplicateStep(UserData& data, size_t step_id, size_t& new_step_id, std::string& err);

// Move a step one place later (+1) or earlier (-1) in the cycle; KPI
// references follow it.
bool MoveStep(UserData& data, size_t step_id, int offset, std::string& err);

// ---- streams --------------------------------------------------------------------------------

// Give `y` an entry for every fluid component (and drop others). When the
// entries then sum to 0 (a new map) the first component is set to 1, so the
// map sums to 1.
void SeedMoleFractions(const UserData& data, ExpressionMap& y);

// Add a stream of `kind` (a set flow, a valve or a pump) to the flow boundary
// on `side` of a step. Its name is "Stream N" made unique on that boundary,
// it takes the boundary's direction, a pump that matches it (a compressor on
// an inlet, a vacuum pump on an outlet) and seeded mole fractions.
// `endpoint_id` is the new stream's endpoint ID. Refused when that end is not
// a flow boundary.
bool AddStream(UserData& data, size_t step_id, BoundaryLocation side, BoundaryStreamType kind, size_t& endpoint_id,
               std::string& err);

// Remove stream `index` of the flow boundary on `side` of a step, with its
// connection (and the partner's), the KPI references to it and the buffer
// tanks that drew from it.
bool RemoveStream(UserData& data, size_t step_id, BoundaryLocation side, size_t index, std::string& err);

// ---- kind switches -------------------------------------------------------------------------
//
// Each replaces a model by one of another kind and keeps every parameter the
// two kinds share: the same name and unit in their descriptors
// (model_kinds.h, CarryOverParameters); the others take the new kind's
// defaults. Switching to the kind already in place changes nothing.

// The isotherm of one component (the other components of its old isotherm
// keep it). Parameters with the same name and unit carry over, so switching to
// the paired kind (PairedIsothermKind) keeps every one.
bool SetIsothermKind(UserData& data, AdsorbentID adsorbent, const std::string& component, IsothermType kind,
                     std::string& err);

// The rate model of one component and its coefficient model (an LDF model is
// created when the component has none). The coefficient kind is a preset of
// the generic resistance model (mass_transfer_model.h): it sets the switches,
// the capacity-factor mode and the parameters the preset shows, which take the
// kind's defaults unless the old preset shared them (same name and unit); the
// parameters it does not show keep their stored values. The rate constant
// keeps its value only while the rate model keeps its order (its unit is 1/s
// for LDF, kg/(mol s) for LDF2). A dual-rate model's slow site takes the same
// preset: a new dual-rate model starts with it a copy of the fast site's, an
// existing one keeps the slow-site parameters the old and new presets share.
bool SetMassTransferKind(UserData& data, AdsorbentID adsorbent, const std::string& component,
                         MassTransferModelType rate_kind, MassTransferCoefficientType coefficient_kind,
                         std::string& err);

// The site of a dual-rate model a resistance edit applies to (single-site
// models have the fast site only).
enum class MassTransferSite { Fast, Slow };

// Switches one resistance ("film", "macropore", "knudsen", "micropore",
// "lumped") of a component's coefficient model on or off, on the fast or the
// slow site. Refused (false, err set): switching the last of film, macropore,
// micropore and lumped off, and knudsen without macropore; switching macropore
// off also switches knudsen off. Parameter values are kept, so switching a
// term off and on again restores them. Each dual-rate site has its own
// switches.
bool SetMassTransferResistance(UserData& data, AdsorbentID adsorbent, const std::string& component,
                               MassTransferSite site, const std::string& tag, bool on, std::string& err);

// The capacity-factor mode of a component's coefficient model (fast or slow
// site); the fixed mode with capacity_factor <= 0 sets it to 1.
bool SetCapacityFactorMode(UserData& data, AdsorbentID adsorbent, const std::string& component, MassTransferSite site,
                           CapacityFactorMode mode, std::string& err);

// The transport models of a layer and the bed's gas-solid heat transfer.
bool SetDispersionKind(UserData& data, LayerID layer, AxialDispersionType kind, std::string& err);
bool SetHeatDispersionKind(UserData& data, LayerID layer, AxialHeatDispersionType kind, std::string& err);
bool SetInternalHeatTransferKind(UserData& data, LayerID layer, InternalHeatTransferType kind, std::string& err);
bool SetExternalHeatTransferKind(UserData& data, LayerID layer, ExternalHeatTransferType kind, std::string& err);
bool SetParticleHeatTransferKind(UserData& data, ParticleHeatTransferType kind, std::string& err);

// The boundary on `side` of a step (a pressure, flow or closed end). The old
// end's connections are disconnected (with their partners), KPI references
// and buffer tanks that pointed at its endpoints dropped. The direction is
// kept (a closed end becomes an outlet). A new pressure end gets a constant
// pressure profile; a new flow end one set-flow stream "Stream 1". The inlet
// temperature and composition carry over between a pressure end and the
// first stream of a flow end.
bool SetBoundaryConditionKind(UserData& data, size_t step_id, BoundaryLocation side, BoundaryType kind,
                              std::string& err);

// The device of stream `index` on `side` of a step: a set flow, a valve or a
// pump. The stream keeps its endpoint ID, name, inlet temperature (and its
// profile), composition and buffer tank, and its pump moves between a set-flow
// stream and a pump stream; its connection is dropped (with the partner's).
bool SetStreamKind(UserData& data, size_t step_id, BoundaryLocation side, size_t index, BoundaryStreamType kind,
                   std::string& err);

// The flow quantity of a set-flow stream. Its profile keeps its kind and the
// parameters that are not a flow (ramp time, factor, rates, weight); the flow
// values take the new kind's defaults, as their unit changed.
bool SetFlowKind(UserData& data, size_t endpoint_id, FlowType kind, std::string& err);

// Which profile of an endpoint SetProfileKind changes.
enum class ProfileTarget {
  Pressure,     // a pressure boundary's pressure, or a valve stream's back pressure
  Flow,         // a set-flow stream's flow
  Temperature,  // the inlet temperature of a stream or pressure boundary
};

// The time profile of one value of an endpoint. Parameters carry over by name
// and unit; a constant's value and a ramp's end value carry over to each
// other (both are the value the profile settles at), and whether a ramp lasts
// the step time is kept. A constant inlet temperature has no profile object:
// the endpoint's temperature is the constant, and a new ramp starts and ends
// there over the step time (CreateInletTemperatureProfile).
bool SetProfileKind(UserData& data, size_t endpoint_id, ProfileTarget target, BoundaryValueType kind,
                    std::string& err);

// The pump of a set-flow stream or pressure boundary: "none" (kind index 0;
// switches the pump off and keeps it) or a PumpType (index type + 1). The
// polytropic index, energy model and efficiency model carry over. A pump
// stream always has a pump (refused for none); a valve stream has none.
bool SetPumpKind(UserData& data, size_t endpoint_id, std::size_t pump_kind, std::string& err);
bool SetPumpEfficiencyKind(UserData& data, size_t endpoint_id, PumpEfficiencyType kind, std::string& err);
// The shaft-power model of a pump stream.
bool SetPumpPowerKind(UserData& data, size_t endpoint_id, PumpPowerType kind, std::string& err);

}  // namespace adsim
