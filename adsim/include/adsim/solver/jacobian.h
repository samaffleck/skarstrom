#pragma once

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
  // IDASetJacFn. Each column is perturbed with exactly the increment IDA's
  // idaLsBandDQJac uses, so every entry written is bit-identical to IDA's.
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
  // own exact row sets because boundary conditions overwrite those residuals.
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

    BandedJacobian() = default;
    ~BandedJacobian();
    BandedJacobian(const BandedJacobian&) = delete;
    BandedJacobian& operator=(const BandedJacobian&) = delete;

    // `ida_mem` is queried for the current step size and error weights.
    void Initialise(void* ida_mem, const Equations& eqs, SUNContext ctx,
                    sunindextype mu, sunindextype ml, JacobianMode mode,
                    bool record_pattern);

    // Called at the start of every cycle step; re-opens pattern learning for
    // the first steps of cycle 0 (each step has its own boundary conditions).
    void OnStepStart(int cycle_index);

    int Evaluate(sunrealtype tt, sunrealtype c_j, N_Vector yy, N_Vector yp, N_Vector rr,
                 SUNMatrix Jac, UserData* udata, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

    long NumResidualEvals() const { return num_res_evals_; }
    long NumAudits() const { return num_audits_; }
    long NumFailureAudits() const { return num_failure_audits_; }
    // Pattern entries discovered by audits after learning (0 = pattern held).
    long NumMisses() const { return num_misses_; }
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
    long NumCheckMismatches() const { return num_check_mismatches_; }
    double MaxCheckRelativeError() const { return max_check_rel_; }
    const std::vector<std::string>& CheckReport() const { return check_report_; }

  private:
    bool IsBoundaryColumn(sunindextype j) const;
    // Records the nonzeros of a band evaluation; returns true if any was new.
    bool Learn(sunindextype j, sunindextype i);
    void Recolour();
    // Band DQ at kProbeStates randomised states around (yy, yp); learns only.
    void ProbeGenericStates(sunrealtype tt, sunrealtype c_j, N_Vector yy, N_Vector yp,
                            UserData* udata, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    int EvaluateGroups(const std::vector<std::vector<sunindextype>>& groups, bool band,
                       sunrealtype tt, sunrealtype c_j, N_Vector yy, N_Vector yp, N_Vector rr,
                       SUNMatrix Jac, UserData* udata, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

    void* ida_mem_ = nullptr;
    N_Vector ewt_ = nullptr;
    JacobianMode mode_ = JacobianMode::BandDQ;
    sunindextype n_ = 0;
    sunindextype mu_ = 0;
    sunindextype ml_ = 0;

    // Banded layout: variable (equation index) and cell of each state entry,
    // and the reverse map index_of_[cell * n_vars_ + var] (-1 if inactive).
    int n_vars_ = 0;
    int n_cells_ = 0;
    std::vector<int> var_;
    std::vector<int> cell_;
    std::vector<sunindextype> index_of_;
    std::vector<std::string> var_names_;

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
