/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// The deck view: a read-only JSON description of a deck for display and
// editing in a front end (the web workspace, the Python package, the CLI).
//
// The view gives the deck's STRUCTURE -- components, adsorbents (per-component
// isotherm and rate-model kinds, reactions), the column and its wall, layers,
// steps with their bottom and top ends (closed, pressure or flow; streams,
// profiles, pumps, buffer tanks), heating media and energy sources,
// connections, the KPI configuration, user variables -- and attaches to each
// node the EDITABLE VALUES that belong to it. Every value is an entry of the
// parameter registry (adsim/parameters/parameter_registry.h), given with its
// registry path and the registry's value, text, unit, bounds, kind, choices
// and label; the view keeps no second table of names or units.
//
// Each node claims the registry entries under its path prefix
// ("endpoint[5].pump." for a pump, "endpoint[5]." for the rest of a stream,
// "layer[1].", "adsorbent[2].isotherm[CarbonDioxide].", ...). Every registry
// entry appears exactly once; an entry no node claims (a registry family this
// view does not know yet) is listed under "other", so nothing is hidden.
//
// Nothing in the view is writable in place: change a value by its path
// (Parameter::SetValue / SetText through ParameterRegistry, or the Simulation
// facade), change structure through objects/deck_builder.h and
// objects/deck_edit.h, then take a new view.
//
// Model kinds are given by name (the builder's kind names, deck_builder.h),
// enumerated values with their choice names. Ids are the objects' stable ids
// (step id, endpoint id, layer id, adsorbent id, reaction id), the same ids
// the registry paths use.
//
// The view's shape is described by DeckViewSchema(), a JSON Schema (draft 4,
// "definitions" + "$ref") whose enumerations are the kind names; front ends
// can generate their types from it. kDeckViewVersion changes when the shape
// changes incompatibly.

#include <string>

namespace adsim {

  class UserData;

  inline constexpr int kDeckViewVersion = 1;

  struct DeckViewOptions {
    bool pretty = false;  // indented JSON (for reading); compact otherwise
  };

  // The view of `data` as JSON text. Values are the cached ones: call
  // UserData::EvaluateAllExpressions first if expressions changed since. The
  // deck is not modified, except that missing step and endpoint ids are
  // assigned (Cycle::EnsureIds) and buffer-tank sources resolved, as the
  // registry and the builder do.
  std::string DeckView(UserData& data, const DeckViewOptions& options = {});

  // JSON Schema (draft 4) of DeckView's output.
  std::string DeckViewSchema(const DeckViewOptions& options = {});

}  // namespace adsim
