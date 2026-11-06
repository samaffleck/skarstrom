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

#include <string>

#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/connection.h"
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
// dropped from the reactions.
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

}  // namespace adsim
