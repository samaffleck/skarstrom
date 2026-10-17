/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

#include <string>
#include <unordered_map>

#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>


namespace adsim {

  class StreamData;

  // Scalar aggregates summarising one endpoint's most-recently-completed run.
  // Updated at the end of every IntegrateStep so that downstream BufferTankSource
  // consumers can read time-averaged inlet conditions for the next time the
  // consumer step fires (this cycle or later).
  struct CycleEndpointSummary {
    bool   is_valid = false;
    double total_moles = 0.0;          // [mol] signed integral of molar flow over the run
    double avg_temperature = 0.0;      // [K] molar-flow-magnitude-weighted average
    std::unordered_map<std::string, double> avg_mole_fraction; // [-] molar-flow-weighted
    int    source_cycle_index = -1;

    void Clear();

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::make_nvp("Is_Valid", is_valid),
        cereal::make_nvp("Total_Moles", total_moles),
        cereal::make_nvp("Average_Temperature", avg_temperature),
        cereal::make_nvp("Average_Mole_Fraction", avg_mole_fraction),
        cereal::make_nvp("Source_Cycle_Index", source_cycle_index)
      );
    }
  };

  // Compute the summary from a StreamData time series. The components in the
  // resulting avg_mole_fraction map are taken from the streamed data itself.
  CycleEndpointSummary ComputeSummary(const StreamData& stream_data, int cycle_index);

} // end adsim namespace

CEREAL_CLASS_VERSION(adsim::CycleEndpointSummary, 1);
