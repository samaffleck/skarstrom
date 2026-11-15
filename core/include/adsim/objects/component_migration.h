/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Load-time rename of components that earlier engines had and this one
// replaces with another (Components::renamedComponents in fluid.h, e.g. Air ->
// Nitrogen). Every deck, adsorbent and cycle load runs the JSON text through
// MigrateRenamedComponents before cereal reads it, so an old file loads with
// the new component in every place that names the old one.

#include <string>

namespace adsim {

// Rewrite `json` (a deck, adsorbent or cycle archive) so that each renamed
// component carries its new name. Returns true if the text changed.
//
// The rename covers:
//   * every component map ({"key": ..., "value": ...} arrays): the key;
//   * every array of names (Fluid "Components", KPI "Products", ...);
//   * the name fields "Component", "Key_Component", "Product_Component" and
//     "Key_Reactant";
//   * parameter paths ("Path"): "[Air]" becomes "[Nitrogen]".
// Stream and endpoint names are not component names and are left alone.
//
// A deck whose fluid has a component list without the old name is not
// touched. When the deck already has the new component, the two are merged:
// the old component's mole fractions are added to the new one's (in the maps
// whose name contains "Mole_Fraction" or "Composition"), and its other
// entries (isotherm, kinetics, heats, heat capacities) are dropped, so the
// new component keeps its own physics.
//
// Text that is not valid JSON is returned unchanged (false): cereal then
// reports the parse error.
bool MigrateRenamedComponents(std::string& json);

}  // namespace adsim
