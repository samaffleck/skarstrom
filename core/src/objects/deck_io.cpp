/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "adsim/objects/deck_io.h"

#include <sstream>
#include <stdexcept>

#include <cereal/archives/json.hpp>

namespace adsim {

void LoadDeckFlexible(UserData& data, const std::string& json) {
  if (json.empty()) {
    throw std::invalid_argument("LoadDeckFlexible: input JSON string is empty");
  }

  try {
    // Preferred: string-format root ("simulation"). The NVP lookup happens at
    // the very first archive read, so if it is absent this throws before `data`
    // is mutated, keeping the fallback safe.
    std::istringstream iss(json);
    cereal::JSONInputArchive ar(iss);
    ar(cereal::make_nvp("simulation", data));
  } catch (const cereal::Exception&) {
    // Fall back to the file-format default root ("value0"). Fresh stream +
    // archive so we re-read from the top.
    std::istringstream iss(json);
    cereal::JSONInputArchive ar(iss);
    ar(data);
  }

  // Mirror adsim::LoadFromString's post-load steps.
  data.kpi_config.FinaliseLegacyLoad(data.fluid.components_);
  data.kpi_config.ResolveIds(data.cycle.step_list_);
  ++data.load_generation_;
}

}  // namespace adsim
