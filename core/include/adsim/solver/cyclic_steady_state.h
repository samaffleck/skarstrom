/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <functional>
#include <limits>
#include <map>
#include <memory>
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
  // The values are stored in decks by number: append new methods, never
  // reorder.
  enum class AccelerationMethod {
    NONE,     // No acceleration (plain successive substitution)
    AITKEN,   // Aitken delta^2 / Graves-Morris vector extrapolation
    ANDERSON, // Anderson mixing over the last few cycles (see AccelerateCyclicSteadyState)
    NEWTON_KRYLOV  // Newton-Krylov on G(x) - x = 0 with finite-difference J v (see AccelerateCyclicSteadyState)
  };

  static const char* AccelerationMethodList[] = {
    "None",
    "Aitken (Graves-Morris)",
    "Anderson",
    "Newton-Krylov"
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
    // Diagnostic only (it never changes the stopping rule): the estimated
    // remaining distance to the cyclic steady state, from the decay of the
    // end-of-cycle state change (see CyclicSteadyState::EstimateDistance).
    // distance is in the norm of state_change; distance_ratio is distance /
    // state_tolerance (NaN when the state test is off). When not estimated,
    // the three numbers are NaN and distance_note says why.
    bool distance_estimated = false;
    double convergence_rate = std::numeric_limits<double>::quiet_NaN();
    double distance = std::numeric_limits<double>::quiet_NaN();
    double distance_ratio = std::numeric_limits<double>::quiet_NaN();
    std::string distance_note;
  };

  // The CSS check of each completed cycle, recorded when the cycle completed.
  // The equation log evicts old cycles (settings.max_cycles_to_keep) but this
  // record is never evicted: one CSSResult per cycle (about 200 bytes with its two
  // short strings), so a 2000-cycle run adds well under 1 MB.
  class CSSHistory {
  public:
    void Clear() { records_.clear(); }
    void Record(int cycle, CSSResult result) { records_[cycle] = std::move(result); }
    bool Has(int cycle) const { return records_.find(cycle) != records_.end(); }
    const CSSResult& Cycle(int cycle) const { return records_.at(cycle); }
    const std::map<int, CSSResult>& Records() const { return records_; }

  private:
    std::map<int, CSSResult> records_;
  };

  // The remaining-distance estimate on its own (CyclicSteadyState::EstimateDistance).
  struct CSSDistanceEstimate {
    bool estimated = false;
    double rate = std::numeric_limits<double>::quiet_NaN();      // rho, in [0, kMaxRate]
    double distance = std::numeric_limits<double>::quiet_NaN();  // last change * rho / (1 - rho)
    double ratio = std::numeric_limits<double>::quiet_NaN();     // distance / tolerance
    int changes_used = 0;   // length of the window rho was fitted on
    std::string note;       // why it was not estimated (empty when estimated)
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
  //    cycle, |in + G - out| / in with G the moles the reactions made in the
  //    bed over the cycle (EquationLog::ReactionGeneration; 0 without
  //    reactions), is within mass_balance_tolerance_, AND the end-of-cycle bed state
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

    // Predicate on a cycle index: true when that cycle started from a state
    // moved by CSS acceleration (a committed jump), so its change from the
    // previous cycle is not a step of the natural cycle-to-cycle map.
    using JumpPredicate = std::function<bool(int)>;

    bool IsAtCSS(const EquationLog& logger, int cycle_index);
    // The stopping rule at cycle_index (0-based): see the class comment. Also
    // fills the remaining-distance diagnostic (EstimateDistance), which does
    // not affect `converged`; started_from_jump marks the cycles to leave out
    // of its window (none when empty).
    CSSResult CheckCSS(const EquationLog& logger, int cycle_index,
                       const JumpPredicate& started_from_jump = {});
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

    // Remaining-distance diagnostic. With a geometric approach at rate rho, a
    // run whose last cycle-to-cycle change is D is still about
    // D * rho / (1 - rho) from the cyclic steady state: the stopping rule
    // bounds the last step, not the distance still to go, so slow decks
    // (rho near 1) stop further from their CSS. The estimate:
    //  * changes: the end-of-cycle state changes (EndOfCycleStateChange) of
    //    consecutive cycles, oldest first, ending at the cycle checked;
    //  * window: the longest trailing run of strictly decreasing, positive,
    //    finite changes, at most kDistanceMaxChanges long. Fewer than
    //    kDistanceMinChanges changes: not estimated ("history too short" when
    //    there were not enough changes, "not a monotone decay" otherwise);
    //  * rho: the geometric mean of the window's ratios,
    //    (D_last / D_first)^(1 / (n - 1)), clipped to [0, kDistanceMaxRate];
    //  * distance = D_last * rho / (1 - rho), ratio = distance / tolerance
    //    (NaN unless tolerance > 0).
    // The log-reading overload collects the changes back from cycle_index
    // while the cycles are logged, and stops before a change whose cycle
    // started from an accelerated state (a jump is not a step of the natural
    // map, and the cycle after it is not either). It uses state_tolerance_ as
    // the tolerance. Near the integrator's noise floor (about 1e-6 at the
    // default reltol) the changes stop decaying, so the window breaks up and
    // the estimate is usually "not estimated" there.
    static constexpr int kDistanceMinChanges = 3;
    static constexpr int kDistanceMaxChanges = 5;
    static constexpr double kDistanceMaxRate = 0.999;
    // At the CSS stop, the solver logs a warning when distance / tolerance
    // exceeds this (a stop more than this many tolerances from the CSS).
    static constexpr double kDistanceWarningRatio = 100.0;
    static CSSDistanceEstimate EstimateDistance(const std::vector<double>& changes, double tolerance);
    CSSDistanceEstimate EstimateDistance(const EquationLog& logger, int cycle_index,
                                         const JumpPredicate& started_from_jump = {}) const;

    // Relative change of the bed state at the end of the cycle (the end of its
    // last step) since the previous cycle: the largest, over the slow fields, of
    // RMS(x_c - x_{c-1}) / scale. The fields and scales are the adsorbed loading
    // of each component (scaled by the RMS total loading, so a trace component's
    // round-off cannot dominate, floored at the solver's absolute tolerance on
    // loadings, 1e-8 mol/kg, so an inert bed's round-off loadings do not either),
    // the gas mole fractions (absolute) and the
    // fluid, solid and wall temperatures (scaled by their RMS). Returns a
    // negative value when either cycle is not in the log.
    static double EndOfCycleStateChange(const EquationLog& logger, int cycle_index);

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
  //
  // The jump is limited to kMaxStepMultiple natural steps: when the fitted
  // coefficient asks for more (a slow mode, rate above about 0.95), the jump is
  // shortened to that limit rather than dropped. A shortened jump along a
  // contracting mode undershoots the fixed point, so it still closes part of
  // the gap. A window that is not a geometric decay is skipped (see
  // Accelerate).
  //
  // Anderson (AccelerationMethod::ANDERSON) moves the same fields along a
  // combination of up to kAndersonDepth earlier natural cycle-to-cycle steps,
  // so it can remove several slow modes at once where Aitken removes one.
  // It is not restarted: it keeps pairs (x, G(x)) of a cycle's start and end
  // state, where G is the cycle-to-cycle map. A pair is kept only when x is
  // known exactly: when the previous cycle ended without a jump, x is that
  // cycle's natural end state. The cycle that follows a jump gives no pair (the
  // caller may have rejected the jump, so its start state is not known), but
  // the older pairs remain valid samples of G and are kept. Jumps are thus at
  // most every other cycle. With residuals f = G(x) - x, the newest pair n and
  // the older pairs j, the jump is
  //     x* = G(x_n) - sum_j gamma_j (G(x_n) - G(x_j)),
  //     gamma = argmin || f_n - sum_j gamma_j (f_n - f_j) ||,
  // in the same field-scaled norm as Aitken. The jump is skipped when the
  // fitted combination does not reduce the residual, and shortened to
  // kMaxStepMultiple natural steps. It goes through the same projection and
  // the caller's re-seed probe as an Aitken jump.
  //
  // Newton-Krylov (AccelerationMethod::NEWTON_KRYLOV) solves F(x) = G(x) - x = 0
  // by Newton's method with GMRES, forming each product J v by a finite
  // difference of the cycle map, one plain cycle per Krylov vector
  // (SolveFixedPointNewtonKrylov in newton_krylov.h). The extra cycles run in
  // a scratch copy of the simulation (CycleMapEvaluator, provided by
  // Simulation for its runs), started from the state the live run holds at
  // the end of the cycle with the unknowns replaced; the live run and its log
  // are never touched.
  //  * Unknowns: every differential field (the slow fields above AND the
  //    total gas concentration and fluid temperature), at its differential,
  //    active entries (ghost and void cells are left out); the algebraic
  //    fields follow from the consistent-IC solve at the cycle start. The
  //    concentration and fluid temperature are not fast enough to hold fixed:
  //    on Perez 2019 A a solve over the slow fields alone converged to 6e-7
  //    in the scratch map, yet the next live cycle moved by 2e-3 again; with
  //    them, the cycle after the jump moved by 3e-7. The endpoint histories
  //    (connected streams, buffer tanks) are held at the live values during a
  //    solve and updated by the logged cycles between solves.
  //  * Scaling: as the state test (CyclicSteadyState::EndOfCycleStateChange):
  //    loadings by the RMS total loading (floored at 1e-8 mol/kg), mole
  //    fractions absolute, temperatures and the concentration by their RMS.
  //    The residual norm is the state test's: the largest, over the fields,
  //    of the RMS scaled residual.
  //  * The difference increment is scaled to the integrator's noise:
  //    fd_relative_step = sqrt(reltol), with reltol the scratch cycles' own.
  //    They run at settings.reltol, until a solve fails to converge; later
  //    solves run them at newton_reltol_factor_ times it. A slow mode of
  //    multiplier mu puts an eigenvalue 1 - mu in I - J, and the difference
  //    J v carries an error of about sqrt(reltol): at the deck's reltol 1e-6
  //    that swamps the slow mode of Nguyen 2023 X (mu = 0.998) and the Newton
  //    steps stop reducing |F|, while on Perez 2019 A and E the first solve
  //    converges and a tighter tolerance would only make the cycles dearer.
  //  * A solve starts from the live end state x0 (one scratch cycle for
  //    G(x0)), runs up to kNewtonMaxIterations damped Newton steps (Glover
  //    rule, one halving when a step does not reduce |F|), and stops at
  //    |F| <= newton_tolerance_ (by default 10 times the scratch cycles'
  //    reltol, so 10 settings.reltol until tightened) or when a step
  //    fails to reduce |F|. The jump is then G of the best iterate (one free
  //    plain cycle beyond it), through the same projection (clipping, mole
  //    fractions renormalised) and re-seed probe as the other methods.
  //  * Solves start at cycle kNewtonFirstCycle. After a solve that did not
  //    converge, Anderson runs for kNewtonCooldownCycles cycles before the
  //    next attempt. After one that did, the run continues with plain cycles
  //    and solves again only when the natural cycle-to-cycle change exceeds
  //    kNewtonRetryFactor times the residual that solve reached (or times
  //    10 settings.reltol, the live cycles' noise, when that is larger).
  //  * Without an evaluator (an engine-level CoupledSolver run, not one
  //    started through Simulation) it behaves as Anderson.
  struct AccelerationNote {
    std::string skipped;  // set when a full three-cycle window is rejected: the reason
    std::string fit;      // set for every full window: the fitted coefficient and rates
  };

  // One unknown of the Newton-Krylov solve: an entry of one equation's
  // profile (Equations::GetEquations() index and cell).
  struct NewtonUnknown {
    size_t equation = 0;
    size_t cell = 0;
  };

  // The cycle map G for the Newton-Krylov acceleration: plain cycles run
  // outside the live run (see AccelerateCyclicSteadyState). Simulation
  // provides one for its runs.
  class CycleMapEvaluator {
  public:
    virtual ~CycleMapEvaluator() = default;
    // Takes the live state at the cycle boundary (the end of the cycle just
    // completed) as the start state of the evaluations that follow; unknowns
    // are the entries Evaluate sets and reads. false when it cannot.
    virtual bool Begin(const std::vector<NewtonUnknown>& unknowns) = 0;
    // gx = G(x): one plain cycle from the Begin state with the unknowns set
    // to x, and the unknowns at its end. false when the cycle fails or a stop
    // was requested.
    virtual bool Evaluate(const std::vector<double>& x, std::vector<double>& gx) = 0;
    // The integrator's relative tolerance (settings.reltol of the run).
    virtual double RelativeTolerance() const = 0;
    // The relative tolerance of the evaluations that follow (<= 0: the
    // run's own).
    virtual void SetEvaluationTolerance(double reltol) = 0;
    // A message for the run log (level medium).
    virtual void Log(const std::string& message) = 0;
  };

  // The Newton-Krylov solve made after one cycle.
  struct NewtonKrylovRecord {
    int cycle = 0;                       // the cycle after which it ran (0-based)
    int newton_iterations = 0;           // accepted Newton steps
    int cycles = 0;                      // scratch cycles run (cycle equivalents)
    std::vector<double> residuals;       // |F| at each accepted iterate, the start first
    std::vector<int> krylov_vectors;     // per Newton step
    double floquet_multiplier = std::numeric_limits<double>::quiet_NaN();
    bool converged = false;              // |F| <= the Newton tolerance
    bool jumped = false;                 // a jump was proposed (subject to the re-seed probe)
    std::string outcome;
  };

  class AccelerateCyclicSteadyState {
  public:
    AccelerationMethod method_ = AccelerationMethod::NONE;

    // Newton-Krylov: stop a solve when the residual (state-test norm) is at
    // most this. 0 (the default) means kNewtonDefaultToleranceFactor times the
    // run's settings.reltol: the state norm's noise floor is about reltol
    // (5e-7 to 1e-6 measured at 1e-6), where Newton steps stop paying, so a
    // deeper CSS is asked for by tightening settings.reltol (down to about
    // 1e-8, where the fixed absolute tolerances take over). Not serialised
    // (an override for tests and C++ callers).
    double newton_tolerance_ = 0.0;
    // Newton-Krylov: after a solve that did not converge, the scratch cycles
    // of later solves run at settings.reltol times this (1: never tighter).
    // Not serialised.
    double newton_reltol_factor_ = kNewtonTightReltolFactor;

    // Largest jump, in multiples of the last natural cycle-to-cycle step.
    static constexpr double kMaxStepMultiple = 20.0;
    // Anderson: the most earlier pairs combined in one jump (the depth m).
    static constexpr int kAndersonDepth = 4;
    // Anderson: skip a jump whose fit leaves more than this fraction of the
    // residual. Such a fit gains little over the natural cycle, and a jump
    // late in a run delays the stopping rule (it needs consecutive quiet
    // cycles). 0.5 was chosen on Perez 2019 A and E and Nguyen 2023 X: 0.7 and
    // 0.9 stopped Perez 2019 A and E as late or later, and 0.9 stopped Nguyen
    // 2023 X at cycle 234 instead of 190 (see docs/numerics.md).
    static constexpr double kAndersonMaxResidualRatio = 0.5;

    // Newton-Krylov (see the comment above AccelerationNote).
    static constexpr int kNewtonFirstCycle = 3;
    static constexpr int kNewtonMaxIterations = 8;
    static constexpr int kNewtonMaxKrylov = 30;
    static constexpr int kNewtonCooldownCycles = 5;
    static constexpr double kNewtonRetryFactor = 10.0;
    static constexpr double kNewtonTightReltolFactor = 0.01;
    static constexpr double kNewtonDefaultToleranceFactor = 10.0;

    // Called after every cycle with the natural end-of-cycle state. Returns
    // true when it has moved the extrapolable fields of state (the caller then
    // validates and commits the re-seed). note, when given, says why a full
    // window was skipped and what was fitted.
    bool Accelerate(Equations& state, int cycle_index, AccelerationNote* note = nullptr);

    // Newton-Krylov: the cycle map for the next run (nullptr: none, then the
    // method behaves as Anderson). Not serialised.
    void SetCycleMap(std::shared_ptr<CycleMapEvaluator> evaluator) { cycle_map_ = std::move(evaluator); }
    // Forgets everything collected in a run (the histories, the Newton
    // records); Simulation calls it before each run.
    void ResetRun();
    // The Newton-Krylov solves of the last run, by cycle.
    const std::map<int, NewtonKrylovRecord>& NewtonRecords() const { return newton_records_; }

  private:
    bool AccelerateAitken(Equations& state, int cycle_index, AccelerationNote* note);
    bool AccelerateAnderson(Equations& state, int cycle_index, AccelerationNote* note);
    // Anderson's bookkeeping and jump without the new-run check;
    // allow_jump = false only records the pair.
    bool AndersonStep(Equations& state, int cycle_index, AccelerationNote* note, bool allow_jump);
    bool AccelerateNewtonKrylov(Equations& state, int cycle_index, AccelerationNote* note);
    void ResetHistories();

    // Rolling history of natural end-of-cycle snapshots. Each snapshot holds one
    // profile per equation in Equations::GetEquations() order. Transient: rebuilt
    // every run, never serialized.
    std::vector<std::vector<std::vector<double>>> history_;
    int last_cycle_index_ = -1;

    // Anderson (transient, like history_): the slow fields flattened into one
    // vector, in GetEquations() order.
    struct AndersonPair {
      std::vector<double> x;  // start-of-cycle state
      std::vector<double> g;  // end-of-cycle state, G(x)
    };
    std::vector<AndersonPair> pairs_;   // oldest first, at most kAndersonDepth + 1
    std::vector<double> previous_end_;  // the previous cycle's end state
    bool previous_jumped_ = false;      // the previous call returned true

    // Newton-Krylov (transient).
    std::shared_ptr<CycleMapEvaluator> cycle_map_;
    std::map<int, NewtonKrylovRecord> newton_records_;
    std::vector<double> newton_previous_end_;  // the previous cycle's end (the unknowns)
    bool newton_previous_jumped_ = false;
    int newton_cooldown_ = 0;                  // Anderson cycles left before the next attempt
    double newton_settled_residual_ = -1.0;    // residual of the last converged solve (< 0: none)
    bool newton_tight_ = false;                // a solve did not converge: later ones run tighter

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
CEREAL_CLASS_VERSION(adsim::AccelerateCyclicSteadyState, 0);
