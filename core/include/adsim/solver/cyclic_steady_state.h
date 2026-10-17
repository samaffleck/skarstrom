/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <vector>
#include <string>

// Cereal includes
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>

#include "adsim/math/expression_param.h"


namespace adsim {

  class StateLogger;
  class EquationLog;
  class Equations;

  enum class CSSMethod {
    NONE,                     // No CSS checking
    MOLAR_PROFILE,            // Compare end-of-cycle molar concentration profiles
    MASS_BALANCE_TOTAL,       // |total_moles_in - total_moles_out| / total_moles_in
    MASS_BALANCE_COMPONENT    // Per-component mass balance check
  };

  static const char* CSSMethodList[] = {
    "None",
    "Molar Profile",
    "Mass Balance (Total)",
    "Mass Balance (Per-Component)"
  };

  // Formula used to accelerate the approach to cyclic steady state. NONE disables
  // acceleration; add further entries here (and a case in Accelerate) to offer
  // additional extrapolation schemes.
  enum class AccelerationMethod {
    NONE,     // No acceleration (plain successive substitution)
    AITKEN    // Aitken delta^2 / Graves-Morris vector extrapolation
  };

  static const char* AccelerationMethodList[] = {
    "None",
    "Aitken (Graves-Morris)"
  };

  struct CSSResult {
    bool converged = false;
    double error = 0.0;
    double tolerance = 0.0;
    std::string method_name;
  };

  class CyclicSteadyState {
  public:
    CSSMethod method_ = CSSMethod::NONE;

    // Profile method settings
    ExpressionParam molar_absolute_tolerance{1e-4};
    ExpressionParam molar_relative_tolerance{1e-4};

    // Mass balance method settings
    ExpressionParam mass_balance_tolerance_{0.005};  // 0.5%

    // Minimum number of cycles to run before ANY method may declare convergence.
    // Applies to every CSS method (a floor on early stopping). A large value hides
    // the benefit of CSS acceleration, which is why it is surfaced in the GUI.
    int min_cycles_ = 5;

    void EvaluateExpressions(const EvalContext& ctx) {
      molar_absolute_tolerance.Evaluate(ctx);
      molar_relative_tolerance.Evaluate(ctx);
      mass_balance_tolerance_.Evaluate(ctx);
    }

    bool IsAtCSS(const EquationLog& logger, int cycle_index);
    CSSResult CheckCSS(const EquationLog& logger, int cycle_index);

  private:
    CSSResult CheckProfileResult(const EquationLog& logger, int cycle_index);
    CSSResult CheckMassBalanceTotalResult(const EquationLog& logger, int cycle_index);
    CSSResult CheckMassBalanceComponentResult(const EquationLog& logger, int cycle_index);

    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::make_nvp("CSS_Method", method_),
        cereal::make_nvp("molar_absolute_tolerance", molar_absolute_tolerance),
        cereal::make_nvp("molar_relative_tolerance", molar_relative_tolerance),
        cereal::make_nvp("mass_balance_tolerance", mass_balance_tolerance_),
        cereal::make_nvp("min_cycles", min_cycles_)
      );
    }
  };

  // Accelerates the approach to cyclic steady state by Aitken Delta^2
  // extrapolation of the end-of-cycle state toward the fixed point of the
  // cycle-to-cycle map.
  //
  // The extrapolation is "restarted": three successive NATURAL end-of-cycle
  // states are collected, the slowly-converging fields are extrapolated with a
  // single shared coefficient, the solver is re-seeded, and the history is cleared
  // so the next extrapolation uses a fresh triple. Restarting is what makes it
  // work -- it spaces the jumps ~3 cycles apart so the bed relaxes between them,
  // and it keeps the three samples on the un-perturbed cycle trajectory (mixing a
  // jumped state into the triple would corrupt the geometric-mode fit).
  //
  // The moved fields are the slow ones (adsorbed loadings, solid and wall
  // temperature, gas composition); the total gas concentration and fluid
  // temperature are deliberately left alone (they set the algebraic pressure and
  // are fast, so re-seeding them makes the state non-integrable for little gain).
  // Moving all the slow fields together, along the self-consistent natural step
  // direction, keeps the re-seed balanced so the fields do not fight and snap
  // back. Each re-seed is validated in a throwaway integrator before it is
  // committed (see CoupledSolver), so acceleration never destabilises a run.
  class AccelerateCyclicSteadyState {
  public:
    AccelerationMethod method_ = AccelerationMethod::NONE;

    bool Accelerate(Equations& state, int cycle_index);

  private:
    // Rolling history of natural end-of-cycle snapshots. Each snapshot holds one
    // profile per equation in Equations::GetEquations() order. Transient: rebuilt
    // every run, never serialized.
    std::vector<std::vector<std::vector<double>>> history_;
    int last_cycle_index_ = -1;

    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("Method", method_));
    }
  };

} // end adsim namespace
