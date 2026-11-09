#pragma once

#include <functional>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <sundials/sundials_context.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

namespace adsim {

  class Equations;
  class UserData;

  // How the Newton-iteration Jacobian is built. A solver implementation detail
  // (not serialised); switchable for benchmarking and verification.
  enum class JacobianMode {
    IdaInternal,  // IDA's own banded difference-quotient Jacobian
    BandDQ,       // our replica of IDA's banded DQ Jacobian (same column groups)
    Coloured,     // DQ over a learned sparsity pattern with coloured column groups
    ColouredCheck,  // Coloured, plus a band DQ comparison at every evaluation (diagnostic)
  };

  // Difference-quotient Jacobian for the banded IDA system, installed with
  // IDASetJacFn. Each column is perturbed with the increment IDA's
  // idaLsBandDQJac uses, so every entry written is bit-identical to IDA's,
  // except the pressure columns: their increment is sqrt(eps) |P| floored at
  // the pressure's absolute tolerance, not IDA's reltol |P| + abstol, which
  // is coarser than the pressure differences the flow depends on (see
  // Increment in jacobian.cpp), and the reaction columns below.
  //
  // Reaction columns. Rate laws clip each activity at 0 and floor an
  // activity raised to a negative power at kMinActivityBar /
  // kMinActivityConcentration (reaction_kinds.h), so a rate is flat in a
  // mole fraction below zero (and below the floor). A reaction participant
  // can sit far below the mole fractions' absolute tolerance (1e-8): a
  // product that the bed does not hold yet, a reactant ahead of its front.
  // IDA's increment, floored at the error weight reltol |y| + abstol,
  // then spans many times the mole fraction, and its sign follows h y':
  //  * a falling mole fraction is perturbed downwards, into y < 0, where the
  //    rate is clipped: the secant of even a linear term is y / increment of
  //    the tangent (1e-4 at y = 1e-12);
  //  * a power below 1 (H2 at -2.5 in steam-reforming kinetics, a square
  //    root in a reverse term) has an unbounded derivative at zero: the
  //    secant across [y, y + 1e-8] is far flatter than the tangent (25x for
  //    y^-2.5 at y = 1e-9).
  // Both starve the modified Newton iteration of the reaction's coupling.
  // On a reaction that relaxes much faster than the flow (a fast catalyst
  // starting from a product-free bed) it fails to converge, and IDA keeps
  // the step at the reaction's time scale. So the mole-fraction column of a
  // component that a reaction in the column's cell involves (its
  // stoichiometry or IReactionRate::MinExponent) is perturbed away from
  // zero (upwards) whenever IDA's increment would cross it, and, when the
  // smallest power is below 1, by the relative increment
  //   kReactionRelIncrement max(y, y_floor)   (y_floor: the activity floor's
  //   mole fraction for a negative power, ActivityFloorMoleFraction, else 0)
  // upwards, whenever that is smaller than IDA's. 1e-4 gives the tangent of
  // y^-n to (n + 1) 5e-5 and stays far above the round-off of the other
  // couplings in the column. Every other column, and every column of a deck
  // without reactions, keeps IDA's increment bit for bit.
  //
  // In Coloured mode the coupling pattern is learned from band DQ evaluations:
  // at a few randomised "generic" states near the current one (so couplings
  // that vanish only at special states such as u = 0 or y = 0 still show up)
  // and at the first few real evaluations of every cycle step in the first
  // cycle (each step has its own boundary conditions). Columns whose row
  // sets do not overlap are then perturbed together, which needs far fewer
  // residual calls per Jacobian than IDA's fixed band stride. Interior columns
  // share a translation-invariant pattern keyed by (variable, cell offset);
  // columns near either end of the bed (ghost/void/boundary cells) keep their
  // own exact row sets because boundary conditions overwrite those residuals,
  // and also take the interior pattern (a coupling one boundary column did not
  // show while learning, e.g. of the WENO / van Leer faces next to a ghost
  // cell, has been seen all along the interior).
  // Every kAuditInterval-th Jacobian is a full band DQ evaluation that checks
  // the pattern; a coupling it finds outside the pattern is added and the
  // columns are recoloured, so a missed dependency is self-correcting. The
  // Jacobian requested after a Newton convergence failure is also a band
  // evaluation: couplings that only exist in rare regimes (e.g. a boundary
  // whose flow reverses at equilibrium) are then present exactly when Newton
  // needs them, so robustness never falls below IDA's own DQ Jacobian.
  class BandedJacobian {
  public:
    static constexpr int kLearnPerStep = 10;    // band evaluations per step in cycle 0
    static constexpr int kAuditInterval = 50;   // every Nth evaluation is a band audit
    static constexpr int kBoundaryCells = 3;    // cells at each end with exact patterns
    static constexpr int kProbeStates = 3;      // generic states probed per step in cycle 0
    static constexpr double kReactionRelIncrement = 1e-4;  // reaction columns (see the class comment)

    BandedJacobian() = default;
    ~BandedJacobian();
    BandedJacobian(const BandedJacobian&) = delete;
    BandedJacobian& operator=(const BandedJacobian&) = delete;

    // `ida_mem` is queried for the current step size and error weights.
    // (mu, ml) are the matrix half-bandwidths; (scan_mu, scan_ml) >= (mu, ml)
    // is the reach band evaluations scan, so couplings outside the matrix
    // band are counted (NumOutOfBand) rather than silently dropped.
    void Initialise(void* ida_mem, const Equations& eqs, SUNContext ctx,
                    sunindextype mu, sunindextype ml, JacobianMode mode,
                    bool record_pattern, sunindextype scan_mu = -1, sunindextype scan_ml = -1);

    // Widens the matrix band to (mu, ml) <= the scan reach after the band was
    // sized too narrow (NumOutOfBand() > 0): the learned pattern is dropped
    // and re-learned, the evaluation counters are kept. The caller swaps the
    // SUNMatrix / linear solver; this only re-sizes the pattern tables.
    void ResizeBand(sunindextype mu, sunindextype ml);

    // Measures the coupling reach of the residual by band difference quotients
    // (stride 2 reach + 1) at the state (y, yp) and at kProbeStates jittered
    // copies of it; returns false if no residual evaluation succeeded.
    // `abstol` sets the increment floor. Leaves udata->state at (y, yp).
    static bool ProbeBandwidth(const Equations& eqs, UserData* udata, SUNContext ctx,
                               N_Vector y, N_Vector yp, N_Vector abstol, sunindextype reach,
                               sunindextype& mu, sunindextype& ml);

    // Called at the start of every cycle step; re-opens pattern learning for
    // the first steps of cycle 0 (each step has its own boundary conditions).
    void OnStepStart(int cycle_index);

    // Testing only (IDASettings::jacobian_learn_once_for_testing_): learn the
    // pattern from the first band evaluation alone, so the audits must
    // complete it. Call after Initialise.
    void SetLearnOnceForTesting(bool on);

    int Evaluate(sunrealtype tt, sunrealtype c_j, N_Vector yy, N_Vector yp, N_Vector rr,
                 SUNMatrix Jac, UserData* udata, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

    // Testing only: a read-only view of a Jacobian as IDA receives it, passed
    // to the observer installed with SetEvaluationObserverForTesting after
    // every successful evaluation. `increments` holds the DQ increment of
    // every column (the same for a band and a coloured evaluation), so an
    // independent difference quotient can reproduce each entry. The observer
    // may evaluate residuals through `udata` (that moves udata->state, which
    // IDA's next residual call overwrites, as ColouredCheck's extra
    // evaluations do) but must not modify the vectors or the matrix.
    struct EvaluationView {
      enum class Kind { Band, Coloured };  // Band: BandDQ, learning or audit evaluation
      Kind kind = Kind::Band;
      sunrealtype tt = 0.0;
      sunrealtype c_j = 0.0;
      N_Vector yy = nullptr;
      N_Vector yp = nullptr;
      SUNMatrix jac = nullptr;
      sunindextype mu = 0;  // matrix half-bandwidths
      sunindextype ml = 0;
      const std::vector<sunrealtype>* increments = nullptr;
      UserData* udata = nullptr;
    };
    using EvaluationObserver = std::function<void(const EvaluationView&)>;
    // Process-wide (the integrator owns its BandedJacobian, tests cannot reach
    // it before the run); pass an empty function to remove. Not thread-safe:
    // set it only around a single-threaded run.
    static void SetEvaluationObserverForTesting(EvaluationObserver observer);

    // Testing only: false gives the reaction columns IDA's increment, as
    // before they were treated specially (to reproduce the time-step
    // collapse). Process-wide, default true; read when a run's Jacobian is
    // first evaluated; not thread-safe.
    static void SetReactionIncrementForTesting(bool on);

    // The state layout, for diagnostics: the variable (equation index) and
    // cell of state entry j, the short name of a variable ("C", "P", "y0",
    // ...), and the component of a mole-fraction variable (empty otherwise).
    sunindextype Size() const { return n_; }
    int VariableOf(sunindextype j) const { return var_[j]; }
    int CellOf(sunindextype j) const { return cell_[j]; }
    const std::string& VariableNameOf(int e) const { return var_names_[e]; }
    const std::string& ComponentOf(int e) const;
    // Mole-fraction columns with a reaction treatment (0 before the first
    // Jacobian evaluation).
    long NumReactionColumns() const { return num_reaction_columns_; }

    // Residual calls spent on Jacobians (not counting ColouredCheck's comparison).
    long NumResidualEvals() const { return num_res_evals_; }
    long NumAudits() const { return num_audits_; }
    long NumFailureAudits() const { return num_failure_audits_; }
    // Pattern entries discovered by audits after learning (0 = pattern held).
    long NumMisses() const { return num_misses_; }
    // Nonzero couplings band evaluations found outside the matrix band (only
    // possible when the band was sized below the scan reach).
    long NumOutOfBand() const { return num_out_of_band_; }
    int NumColours() const { return static_cast<int>(colour_groups_.size()); }
    // Half-bandwidths of the banded system (as passed to Initialise).
    sunindextype UpperBandwidth() const { return mu_; }
    sunindextype LowerBandwidth() const { return ml_; }

    // Observed nonzero coupling over all band evaluations, keyed by
    // (row variable, column variable, row cell - column cell) -> count.
    // Only populated when Initialise was called with record_pattern = true.
    using PatternKey = std::tuple<std::string, std::string, int>;
    const std::map<PatternKey, long>& PatternCounts() const { return pattern_counts_; }

    // ColouredCheck mode: entries where the coloured and band matrices differed.
    // The comparison only observes: the pattern, its audits and every counter
    // above evolve exactly as in Coloured mode, so a coupling the pattern
    // misses is reported here until an audit adds it (NumMisses).
    long NumCheckMismatches() const { return num_check_mismatches_; }
    double MaxCheckRelativeError() const { return max_check_rel_; }
    const std::vector<std::string>& CheckReport() const { return check_report_; }

  private:
    // Evaluate without the testing observer.
    int EvaluateMatrix(sunrealtype tt, sunrealtype c_j, N_Vector yy, N_Vector yp, N_Vector rr,
                       SUNMatrix Jac, UserData* udata, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    void InitialiseBand();  // pattern tables and column groups for (mu_, ml_, scan_mu_, scan_ml_)
    bool IsBoundaryColumn(sunindextype j) const;
    // Records the nonzeros of a band evaluation; returns true if any was new.
    bool Learn(sunindextype j, sunindextype i);
    void Recolour();
    // Band DQ at kProbeStates randomised states around (yy, yp); learns only.
    void ProbeGenericStates(sunrealtype tt, sunrealtype c_j, N_Vector yy, N_Vector yp,
                            UserData* udata, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    // What a group evaluation is for. Band: the full band (BandDQ mode, or a
    // learning / audit evaluation, which learns the nonzeros it finds outside
    // the pattern). Coloured: only the pattern rows of each column are read.
    // Check: the full band for the ColouredCheck comparison, purely observing
    // (no learning, no counters), so the diagnostic does not change the run.
    enum class Pass { Band, Coloured, Check };
    int EvaluateGroups(const std::vector<std::vector<sunindextype>>& groups, Pass pass,
                       sunrealtype tt, sunrealtype c_j, N_Vector yy, N_Vector yp, N_Vector rr,
                       SUNMatrix Jac, UserData* udata, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

    // Fills reaction_column_ / reaction_floor_ from the reactions of each
    // cell (once per Initialise, at the first evaluation).
    void MarkReactionColumns(UserData& udata);
    // The DQ increment of column j: IDA's (Increment in jacobian.cpp), or the
    // reaction increment of a reaction column.
    sunrealtype ColumnIncrement(sunindextype j, sunrealtype yj, sunrealtype ypj, sunrealtype ewtj,
                                sunrealtype hh, sunrealtype srur) const;

    void* ida_mem_ = nullptr;
    N_Vector ewt_ = nullptr;
    JacobianMode mode_ = JacobianMode::BandDQ;
    sunindextype n_ = 0;
    sunindextype mu_ = 0;
    sunindextype ml_ = 0;
    sunindextype scan_mu_ = 0;
    sunindextype scan_ml_ = 0;
    long num_out_of_band_ = 0;

    // Banded layout: variable (equation index) and cell of each state entry,
    // and the reverse map index_of_[cell * n_vars_ + var] (-1 if inactive).
    int n_vars_ = 0;
    int n_cells_ = 0;
    int velocity_var_ = -1;  // equation index of the face velocity (-1 if absent)
    int pressure_var_ = -1;  // equation index of the pressure (-1 if absent)
    double pressure_abstol_ = 0.0;  // the pressure's absolute tolerance (its DQ increment floor)
    std::vector<int> var_;
    std::vector<int> cell_;
    std::vector<sunindextype> index_of_;
    std::vector<std::string> var_names_;
    std::vector<std::string> var_components_;  // component of a mole-fraction variable, else empty
    // Per column (empty: not marked yet): its reaction treatment, and the
    // activity floor's mole fraction of a Singular column with a negative
    // power (0 otherwise).
    enum class ReactionColumn : char { None = 0, Participant = 1, Singular = 2 };
    std::vector<ReactionColumn> reaction_column_;
    std::vector<double> reaction_floor_;
    long num_reaction_columns_ = 0;

    // Learned pattern. Interior: interior_[(col_var * n_vars_ + row_var) * n_offsets_ + (d + max_offset_)].
    // Boundary columns: boundary_[j][i - j + mu_] over the band.
    int max_offset_ = 0;
    int n_offsets_ = 0;
    std::vector<char> interior_;
    std::vector<std::vector<char>> boundary_;

    // Rows read back per column (CSR) and the column groups for each mode.
    std::vector<sunindextype> row_ptr_;
    std::vector<sunindextype> rows_;
    std::vector<std::vector<sunindextype>> band_groups_;
    std::vector<std::vector<sunindextype>> colour_groups_;

    int learn_budget_ = 0;
    bool learn_once_ = false;
    bool probe_pending_ = false;
    std::vector<N_Vector> probe_vecs_;  // y, yp, r at a probe state
    bool learning_done_ = false;
    int evals_since_audit_ = 0;
    long new_entries_ = 0;  // pattern entries learned by the current band evaluation
    long num_res_evals_ = 0;
    long num_audits_ = 0;
    long num_failure_audits_ = 0;
    long last_conv_fails_ = 0;  // IDA convergence-failure count at the previous evaluation
    long num_misses_ = 0;

    bool record_pattern_ = false;
    std::map<PatternKey, long> pattern_counts_;

    bool last_pass_coloured_ = false;  // the latest Evaluate was a coloured pass

    SUNMatrix check_matrix_ = nullptr;
    long num_check_mismatches_ = 0;
    double max_check_rel_ = 0.0;
    std::vector<std::string> check_report_;
  };

  // IDALsJacFn trampoline; user_data is the UserData whose jacobian_ is set.
  int BandedJacobianFn(sunrealtype tt, sunrealtype c_j, N_Vector yy, N_Vector yp,
                       N_Vector rr, SUNMatrix Jac, void* user_data,
                       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

} // end adsim namespace
