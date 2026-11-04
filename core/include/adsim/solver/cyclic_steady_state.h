/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <functional>
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
    // True when the method's stopping rule holds at this cycle: its single-cycle
    // test passed on `required_cycles` consecutive cycles ending here, and the
    // minimum number of cycles has run.
    bool converged = false;
    // The method's own measure at this cycle (the relative balance error for the
    // mass-balance methods, the profile change for Molar Profile) and its tolerance.
    double error = 0.0;
    double tolerance = 0.0;
    std::string method_name;
    // The relative cycle-to-cycle change of the end-of-cycle bed state
    // (EndOfCycleStateChange), checked alongside the balance by the
    // mass-balance methods. state_tolerance is 0 when the state test is off.
    double state_change = 0.0;
    double state_tolerance = 0.0;
    // Consecutive cycles, ending at this one, on which the single-cycle test
    // passed (capped at required_cycles), and how many the rule needs.
    int passing_cycles = 0;
    int required_cycles = 1;
  };

  // Cyclic-steady-state detection: the stopping rule applied after every cycle.
  //
  // Each method has a single-cycle test; the run stops at the first cycle at
  // which that test has passed on consecutive_cycles_ cycles in a row (and at
  // least min_cycles_ cycles have run).
  //
  //  * Molar Profile: the total-concentration profile at the end of the first
  //    step changed by at most abs + rel * |profile| since the previous cycle.
  //  * Mass Balance (Total / Per-Component): the relative balance error of the
  //    cycle is within mass_balance_tolerance_, AND the end-of-cycle bed state
  //    changed by at most state_tolerance_ since the previous cycle
  //    (EndOfCycleStateChange). The balance alone is not a safe stopping test:
  //    on the way to CSS it can swing through zero (accumulation changing sign),
  //    so a single cycle can pass by accident while purity and recovery are
  //    still moving by a point or more. The state change is a norm, so it does
  //    not cross zero, and it measures the drift directly. state_tolerance_ = 0
  //    switches the state test off (the balance-only test of class version < 2).
  //
  // The check is a pure function of the equation log, so it needs the last
  // CyclesToLog() cycles to be logged; the solver's logs keep at least that many.
  class CyclicSteadyState {
  public:
    CSSMethod method_ = CSSMethod::NONE;

    // Profile method settings
    ExpressionParam molar_absolute_tolerance{1e-4};
    ExpressionParam molar_relative_tolerance{1e-4};

    // Mass balance method settings
    ExpressionParam mass_balance_tolerance_{0.005};  // 0.5%

    // Largest relative cycle-to-cycle change of the end-of-cycle bed state that
    // the mass-balance methods accept (EndOfCycleStateChange). 0 disables it.
    ExpressionParam state_tolerance_{kDefaultStateTolerance};

    // Minimum number of cycles to run before ANY method may declare convergence.
    // Applies to every CSS method (a floor on early stopping). A large value hides
    // the benefit of CSS acceleration, which is why it is surfaced in the GUI.
    int min_cycles_ = 5;

    // Number of consecutive cycles on which the method's single-cycle test must
    // pass before CSS is declared (values below 1 count as 1).
    int consecutive_cycles_ = kDefaultConsecutiveCycles;

    static constexpr double kDefaultStateTolerance = 2e-5;
    static constexpr int kDefaultConsecutiveCycles = 3;

    void EvaluateExpressions(const EvalContext& ctx) {
      molar_absolute_tolerance.Evaluate(ctx);
      molar_relative_tolerance.Evaluate(ctx);
      mass_balance_tolerance_.Evaluate(ctx);
      state_tolerance_.Evaluate(ctx);
    }

    bool IsAtCSS(const EquationLog& logger, int cycle_index);
    // The stopping rule at cycle_index (0-based): see the class comment.
    CSSResult CheckCSS(const EquationLog& logger, int cycle_index);
    // The method's single-cycle test at cycle_index, without the
    // consecutive-cycles and minimum-cycles rules (converged = the test passed).
    CSSResult CheckCycle(const EquationLog& logger, int cycle_index);

    // Number of most recent cycles the check reads from the log (0 for None).
    int CyclesToLog() const;

    // The rule's two building blocks, on plain numbers (CheckCSS applies them to
    // the log; tests apply them to synthetic histories).
    //
    // Single-cycle test of the mass-balance methods: the balance error within its
    // tolerance and, when state_tolerance > 0, the state change within its own
    // (a negative state_change, meaning "not available", fails).
    static bool MassBalanceCyclePasses(double balance_error, double balance_tolerance,
                                       double state_change, double state_tolerance);
    // Number of consecutive cycles ending at cycle_index (counting back to cycle
    // 0 at most) for which passes(c) is true, capped at `required`.
    static int TrailingPasses(int cycle_index, int required, const std::function<bool(int)>& passes);

    // Relative change of the bed state at the end of the cycle (the end of its
    // last step) since the previous cycle: the largest, over the slow fields, of
    // RMS(x_c - x_{c-1}) / scale. The fields and scales are the adsorbed loading
    // of each component (scaled by the RMS total loading, so a trace component's
    // round-off cannot dominate), the gas mole fractions (absolute) and the
    // fluid, solid and wall temperatures (scaled by their RMS). Returns a
    // negative value when either cycle is not in the log.
    // Includes the tanks' states (TankStateChange): a large tank approaches
    // CSS slowly, and the bed alone could pass while the tank still drifts.
    static double EndOfCycleStateChange(const EquationLog& logger, int cycle_index);
    // The tanks' part of EndOfCycleStateChange: the largest, over the tanks,
    // of the relative change of the total moles, the absolute change of the
    // mole fractions and the relative change of the temperature (and wall
    // temperature) at the end of the cycle since the previous cycle. 0 when
    // there are no tanks, negative when either cycle is not in the log.
    static double TankStateChange(const EquationLog& logger, int cycle_index);
    // The bed's part of EndOfCycleStateChange (no tanks).
    static double BedStateChange(const EquationLog& logger, int cycle_index);
    // Largest relative balance error of a tank over the cycle,
    // |in - out| / in over every port and external stream (per component
    // when `component` is non-empty); 0 when no tank had any inflow.
    static double TankBalanceError(const EquationLog& logger, int cycle_index, const std::string& component = {});

  private:
    CSSResult CheckProfileResult(const EquationLog& logger, int cycle_index);
    CSSResult CheckMassBalanceTotalResult(const EquationLog& logger, int cycle_index);
    CSSResult CheckMassBalanceComponentResult(const EquationLog& logger, int cycle_index);
    void ApplyStateTest(const EquationLog& logger, int cycle_index, CSSResult& result);

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
      // Version 2 added the state test and the consecutive-cycles rule. Decks
      // saved before that load with the defaults, i.e. the robust rule.
      if (version >= 2) {
        ar(
          cereal::make_nvp("state_tolerance", state_tolerance_),
          cereal::make_nvp("consecutive_cycles", consecutive_cycles_)
        );
      } else if constexpr (Archive::is_loading::value) {
        state_tolerance_ = ExpressionParam{kDefaultStateTolerance};
        consecutive_cycles_ = kDefaultConsecutiveCycles;
      }
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

// Here rather than in the .cpp, so that every archive that serialises the CSS
// settings sees the version. 2: state_tolerance and consecutive_cycles.
CEREAL_CLASS_VERSION(adsim::CyclicSteadyState, 2);
