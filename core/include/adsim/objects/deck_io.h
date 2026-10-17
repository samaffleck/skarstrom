/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Deck I/O from JSON text, accepting both on-disk root shapes.
//
// Decks exist with two cereal root names:
//   * "simulation" - written by adsim::SaveToString (front ends that save
//     decks through a string, e.g. a browser download);
//   * "value0"     - cereal's default root, written by adsim::SaveToJson
//     (the bundled examples/decks/*.json and anything saved to a file).
// adsim::LoadFromString only accepts "simulation". LoadDeckFlexible accepts
// either, so a front end can hand over whatever bytes the user opened without
// parsing them. adsim::Simulation::LoadJson applies the same rule.

#include <string>

#include "adsim/objects/user_data.h"

namespace adsim {

// Load a deck from JSON text into `data`. Tries the "simulation" root first,
// then the default "value0" root, and runs the same post-load steps as
// LoadFromString (legacy KPI configuration upgrade, endpoint id resolution,
// load-generation bump). Throws std::invalid_argument on empty input and
// cereal::Exception / std::exception on malformed decks.
void LoadDeckFlexible(UserData& data, const std::string& json);

}  // namespace adsim
