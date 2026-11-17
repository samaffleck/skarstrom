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
#include <cmath>
#include <cstddef>

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

  // Scale of the smooth floor on the partner mole fractions in Wilke's equation
  // (dimensionless). It is 100 times the absolute tolerance of the solver on a
  // mole fraction (1e-8), which is also the difference-quotient increment of a
  // mole-fraction column near zero. Below about 1e-8 a mole fraction is solver
  // noise, and Wilke's equation weights the partners by their ratios; this floor
  // keeps those noise ratios out of D_m,i, and the Jacobian's increment
  // resolves the transition (1 % of its width).
  inline constexpr double kWilkeMoleFractionFloor = 1e-6;

  // Smooth floor of a mole fraction: 0.5 (y + sqrt(y^2 + d^2)). It is
  // infinitely differentiable and > 0 for every y. It equals
  // y + d^2 / (4 y) + O(d^4 / y^3) for y >> d, d / 2 at y = 0, and about
  // d^2 / (4 |y|) for y << -d (a solver overshoot). The y < 0 branch is the
  // same function written without the cancellation of y + sqrt(y^2 + d^2).
  inline double SmoothMoleFractionFloor(double y, double d) {
    const double r = std::sqrt(y * y + d * d);
    return y >= 0.0 ? 0.5 * (y + r) : 0.5 * d * d / (r - y);
  }

  // Mixture-averaged molecular diffusivity of component i by Wilke's equation
  // (Wilke 1950; Blanc's law, Blanc 1908, in the trace limit),
  //   D_m,i = sum_j y_j / sum_j (y_j / D_ij),   j != i,
  // with each partner mole fraction y_j replaced by its smooth floor
  // (SmoothMoleFractionFloor, scale kWilkeMoleFractionFloor). sum_j y_j is
  // 1 - y_i when the mole fractions sum to 1.
  //
  // Without the floor the formula is 0 / 0 when every partner is absent (a
  // cell of pure i in a deck of three or more components), and the old form
  // returned 0, so every diffusion-based mass-transfer coefficient of i was 0.
  // With the floor:
  //  - D_m,i is a weighted harmonic mean of the D_ij with weights > 0, so it is
  //    finite and lies between the smallest and the largest D_ij;
  //  - at the pure limit it is the unweighted harmonic mean of the D_ij;
  //  - it is infinitely differentiable in every y_j and does not depend on the
  //    order of the partners;
  //  - where every partner's fraction is >= 1e-3 it is Wilke's equation to a
  //    relative 2e-7 (2e-5 at >= 1e-4, 2e-3 at >= 1e-5): only a nearly pure
  //    cell sees the floor;
  //  - the solver's difference quotient (increment about 1e-8 near y = 0)
  //    gives its derivative to within 1 % across the transition.
  // A Newton probe with P < 0 makes every D_ij < 0; the magnitude of the
  // denominator keeps D_m,i > 0 there, as the previous form did.
  inline double WilkeMixtureDiffusivity(const double* y, const double* D_ij, std::size_t n_partners,
                                        double floor = kWilkeMoleFractionFloor) {
    double num = 0.0, den = 0.0;
    for (std::size_t j = 0; j < n_partners; ++j) {
      const double yj = SmoothMoleFractionFloor(y[j], floor);
      num += yj;
      den += yj / D_ij[j];
    }
    return num / std::abs(den);
  }

  class BinaryDiffusionCoefficientModel {
  public:
    void Initialise(int numberOfCells, const std::vector<std::string>& componentList);
    void Update(const std::vector<double>& T, const std::vector<double>& P);

    std::unordered_map<ComponentPair, BinaryLenardJonesParamaters, ComponentPairHash> m_DiffusionCoefficients{};

  private:
    std::vector<double> top_scratch_{}; // per-cell numerator scratch, reused across Update calls
    TemperatureCache top_cache_{};      // per-cell memo of C T^1.5 (the Chapman-Enskog numerator)
  };

} // adsim namespace
