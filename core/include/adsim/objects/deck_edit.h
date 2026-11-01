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
void RemoveComponentFromDeck(UserData& data, const std::string& name);

// ---- layers and adsorbents -----------------------------------------------------

// Give a layer the lowest-id adsorbent of the directory at mass fraction 1.0
// (no-op when the directory is empty).
void SeedLayerComposition(UserData& data, MixedAdsorbentLayer& layer);

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
// the deck's fluid lacks are kept but unused. On success `id` is the new
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
// layer left empty with another adsorbent (if one remains).
void RemoveAdsorbentFromDeck(UserData& data, AdsorbentID id);

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
