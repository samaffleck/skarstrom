/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>

// Eigen includes
#include <Eigen/Core>
#include "adsim/equations/temperature_cache.h"


namespace adsim {

  struct ComponentPair {
    std::string component1{};
    std::string component2{};

    bool operator==(const ComponentPair& other) const {
      return (component1 == other.component1 && component2 == other.component2) ||
        (component1 == other.component2 && component2 == other.component1);
    }
  };

  struct ComponentPairHash {
    std::size_t operator()(const ComponentPair& pair) const {
      auto hash1 = std::hash<std::string>()(pair.component1);
      auto hash2 = std::hash<std::string>()(pair.component2);
      // XOR the hashes of the two components
      return hash1 ^ hash2;
    }
  };

  struct BinaryLenardJonesParamaters {
    std::vector<double> Dm_ij{};
    std::vector<double> omega{};
    double eps{};
    double sigma{};
    double mr_avr_invers_ij{}; // Mr_i^-1 + Mr_j^-1 
    TemperatureCache omega_cache{}; // per-cell memo of the collision integral omega(T / eps)
  };

  class BinaryDiffusionCoefficientModel {
  public:
    void Initialise(int numberOfCells, const std::vector<std::string>& componentList);
    void Update(const std::vector<double>& T, const std::vector<double>& P);

    std::unordered_map<ComponentPair, BinaryLenardJonesParamaters, ComponentPairHash> m_DiffusionCoefficients{};

  private:
    std::vector<double> top_scratch_{}; // per-cell numerator scratch, reused across Update calls
    TemperatureCache top_cache_{};      // per-cell memo of 2.66e-7 T^1.5
  };

} // adsim namespace
