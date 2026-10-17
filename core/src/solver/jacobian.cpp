#include "../pch.h"
#include "adsim/solver/jacobian.h"
#include "adsim/objects/user_data.h"

#include <ida/ida.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_band.h>

#include <cmath>
#include <random>

namespace adsim {

  namespace {

    // Short, deck-independent name for equation-type index `e` (registration
    // order in Equations::Initialise: C, P, Tf, Ts, Tw, u, y[k]..., q[a][k]...;
    // an isothermal run registers no Tf / Ts / Tw, so the fixed block is
    // matched by identity and the per-component block starts after it).
    std::string VariableName(const Equations& eqs, size_t e) {
      const auto& equations = eqs.GetEquations();
      const std::pair<const Equation*, const char*> fixed[] = {
        { &eqs.concentration, "C" }, { &eqs.pressure, "P" },
        { &eqs.fluid_temperature, "Tf" }, { &eqs.solid_temperature, "Ts" },
        { &eqs.wall_temperature, "Tw" }, { &eqs.velocity, "u" } };
      size_t n_fixed = 0;
      for (const auto* eq : equations)
        for (const auto& [ptr, name] : fixed)
          if (eq == ptr) ++n_fixed;
      if (e < equations.size())
        for (const auto& [ptr, name] : fixed)
          if (equations[e] == ptr) return name;
      const size_t n_components = eqs.mole_fraction.size();
      const size_t k = e - n_fixed;
      if (k < n_components) return "y" + std::to_string(k);
      const size_t q = k - n_components;
      return "q" + std::to_string(q / n_components) + "_" + std::to_string(q % n_components);
    }

    // IDA's DQ increment for column j (idaLsBandDQJac, no constraints).
    inline sunrealtype Increment(sunrealtype yj, sunrealtype ypj, sunrealtype ewtj,
                                 sunrealtype hh, sunrealtype srur) {
      sunrealtype inc = SUNMAX(srur * SUNMAX(SUNRabs(yj), SUNRabs(hh * ypj)), SUN_RCONST(1.0) / ewtj);
      if (hh * ypj < SUN_RCONST(0.0)) { inc = -inc; }
      return (yj + inc) - yj;
    }

  } // namespace

  BandedJacobian::~BandedJacobian() {
    if (ewt_) N_VDestroy(ewt_);
    for (N_Vector v : probe_vecs_) N_VDestroy(v);
    if (check_matrix_) SUNMatDestroy(check_matrix_);
  }

  void BandedJacobian::Initialise(void* ida_mem, const Equations& eqs, SUNContext ctx,
                                  sunindextype mu, sunindextype ml, JacobianMode mode,
                                  bool record_pattern) {
    ida_mem_ = ida_mem;
    mode_ = mode;
    mu_ = mu;
    ml_ = ml;
    record_pattern_ = record_pattern;
    num_res_evals_ = 0;
    num_audits_ = 0;
    num_failure_audits_ = 0;
    last_conv_fails_ = 0;
    num_misses_ = 0;
    pattern_counts_.clear();

    // Banded layout: same loop order as Equations::PackBanded.
    const auto& equations = eqs.GetEquations();
    const auto* void_mask = &eqs.domain_vectors.is_void_cell;
    size_t max_size = 0;
    for (const auto* eq : equations) max_size = std::max(max_size, eq->Size());

    n_vars_ = static_cast<int>(equations.size());
    n_cells_ = static_cast<int>(max_size);
    var_.clear();
    cell_.clear();
    index_of_.assign(static_cast<size_t>(n_cells_) * n_vars_, -1);
    for (int cell = 0; cell < n_cells_; ++cell) {
      for (int e = 0; e < n_vars_; ++e) {
        if (!equations[e]->IsActiveAt(cell, void_mask)) continue;
        index_of_[static_cast<size_t>(cell) * n_vars_ + e] = static_cast<sunindextype>(var_.size());
        var_.push_back(e);
        cell_.push_back(cell);
      }
    }
    n_ = static_cast<sunindextype>(var_.size());

    var_names_.clear();
    for (int e = 0; e < n_vars_; ++e)
      var_names_.push_back(VariableName(eqs, e));

    if (ewt_) N_VDestroy(ewt_);
    ewt_ = N_VNew_Serial(n_, ctx);
    for (N_Vector v : probe_vecs_) N_VDestroy(v);
    probe_vecs_.clear();
    if (mode_ != JacobianMode::BandDQ)
      for (int k = 0; k < 3; ++k) probe_vecs_.push_back(N_VNew_Serial(n_, ctx));
    if (check_matrix_) SUNMatDestroy(check_matrix_);
    check_matrix_ = nullptr;
    num_check_mismatches_ = 0;
    max_check_rel_ = 0.0;
    check_report_.clear();
    if (mode_ == JacobianMode::ColouredCheck) {
      // Same storage bandwidth IDA's matrix gets from SUNBandMatrix(n, mu, ml).
      check_matrix_ = SUNBandMatrix(n_, mu_, ml_, ctx);
    }

    // IDA's column groups: stride mu + ml + 1.
    const sunindextype width = mu_ + ml_ + 1;
    band_groups_.assign(SUNMIN(width, n_), {});
    for (sunindextype g = 0; g < static_cast<sunindextype>(band_groups_.size()); ++g)
      for (sunindextype j = g; j < n_; j += width) band_groups_[g].push_back(j);

    // Largest cell offset reachable inside the band, for the interior table.
    max_offset_ = 0;
    for (sunindextype j = 0; j < n_; ++j) {
      max_offset_ = std::max(max_offset_, cell_[SUNMIN(j + ml_, n_ - 1)] - cell_[j]);
      max_offset_ = std::max(max_offset_, cell_[j] - cell_[SUNMAX(0, j - mu_)]);
    }
    n_offsets_ = 2 * max_offset_ + 1;
    interior_.assign(static_cast<size_t>(n_vars_) * n_vars_ * n_offsets_, 0);
    boundary_.assign(n_, {});
    for (sunindextype j = 0; j < n_; ++j)
      if (IsBoundaryColumn(j)) boundary_[j].assign(mu_ + ml_ + 1, 0);

    row_ptr_.clear();
    rows_.clear();
    colour_groups_.clear();
    learn_budget_ = kLearnPerStep;  // covers IDACalcIC and the first step
    probe_pending_ = mode_ != JacobianMode::BandDQ;
    evals_since_audit_ = 0;
  }

  void BandedJacobian::OnStepStart(int cycle_index) {
    if (mode_ != JacobianMode::BandDQ && cycle_index == 0) {
      learn_budget_ = kLearnPerStep;
      probe_pending_ = true;
    }
  }

  bool BandedJacobian::IsBoundaryColumn(sunindextype j) const {
    return cell_[j] < kBoundaryCells || cell_[j] >= n_cells_ - kBoundaryCells;
  }

  bool BandedJacobian::Learn(sunindextype j, sunindextype i) {
    char* slot = nullptr;
    if (IsBoundaryColumn(j)) {
      slot = &boundary_[j][i - j + mu_];
    } else {
      const int d = cell_[i] - cell_[j];
      slot = &interior_[(static_cast<size_t>(var_[j]) * n_vars_ + var_[i]) * n_offsets_ + d + max_offset_];
    }
    if (*slot) return false;
    *slot = 1;
    return true;
  }

  void BandedJacobian::Recolour() {
    // Rows of each column under the current pattern, clipped to the band.
    row_ptr_.assign(n_ + 1, 0);
    rows_.clear();
    for (sunindextype j = 0; j < n_; ++j) {
      const sunindextype lo = SUNMAX(0, j - mu_);
      const sunindextype hi = SUNMIN(j + ml_, n_ - 1);
      if (IsBoundaryColumn(j)) {
        for (sunindextype i = lo; i <= hi; ++i)
          if (boundary_[j][i - j + mu_]) rows_.push_back(i);
      } else {
        const size_t base = static_cast<size_t>(var_[j]) * n_vars_ * n_offsets_;
        for (int e = 0; e < n_vars_; ++e) {
          for (int d = -max_offset_; d <= max_offset_; ++d) {
            if (!interior_[base + static_cast<size_t>(e) * n_offsets_ + d + max_offset_]) continue;
            const int c = cell_[j] + d;
            if (c < 0 || c >= n_cells_) continue;
            const sunindextype i = index_of_[static_cast<size_t>(c) * n_vars_ + e];
            if (i >= lo && i <= hi) rows_.push_back(i);
          }
        }
      }
      row_ptr_[j + 1] = static_cast<sunindextype>(rows_.size());
    }

    // Greedy colouring: columns sharing a row must be perturbed separately.
    // Columns with no learned rows are never perturbed.
    std::vector<int> colour(n_, -1);
    std::vector<std::vector<sunindextype>> row_cols(n_);
    std::vector<sunindextype> stamp;  // stamp[c] == j  <=> colour c is taken for column j
    colour_groups_.clear();
    for (sunindextype j = 0; j < n_; ++j) {
      if (row_ptr_[j] == row_ptr_[j + 1]) continue;
      for (sunindextype k = row_ptr_[j]; k < row_ptr_[j + 1]; ++k)
        for (sunindextype other : row_cols[rows_[k]]) stamp[colour[other]] = j;
      int c = 0;
      while (c < static_cast<int>(stamp.size()) && stamp[c] == j) ++c;
      if (c == static_cast<int>(stamp.size())) {
        stamp.push_back(-1);
        colour_groups_.emplace_back();
      }
      colour[j] = c;
      colour_groups_[c].push_back(j);
      for (sunindextype k = row_ptr_[j]; k < row_ptr_[j + 1]; ++k) row_cols[rows_[k]].push_back(j);
    }
  }

  int BandedJacobian::Evaluate(sunrealtype tt, sunrealtype c_j, N_Vector yy, N_Vector yp,
                               N_Vector rr, SUNMatrix Jac, UserData* udata,
                               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    if (mode_ == JacobianMode::BandDQ) {
      return EvaluateGroups(band_groups_, true, tt, c_j, yy, yp, rr, Jac, udata, tmp1, tmp2, tmp3);
    }

    if (probe_pending_) {
      probe_pending_ = false;
      ProbeGenericStates(tt, c_j, yy, yp, udata, tmp1, tmp2, tmp3);
    }

    // A Newton convergence failure since the previous Jacobian means the
    // iteration matrix was not good enough: rebuild this one from the full band.
    long conv_fails = 0;
    IDAGetNumNonlinSolvConvFails(ida_mem_, &conv_fails);
    const bool after_failure = conv_fails > last_conv_fails_;
    last_conv_fails_ = conv_fails;

    const bool learning = learn_budget_ > 0 || colour_groups_.empty();
    const bool audit = !learning && (after_failure || ++evals_since_audit_ >= kAuditInterval);
    if (audit && after_failure) ++num_failure_audits_;
    if (!learning && !audit) {
      int retval = EvaluateGroups(colour_groups_, false, tt, c_j, yy, yp, rr, Jac, udata, tmp1, tmp2, tmp3);
      if (retval == 0 && mode_ == JacobianMode::ColouredCheck) {
        SUNMatZero(check_matrix_);
        const long saved_new = new_entries_;
        retval = EvaluateGroups(band_groups_, true, tt, c_j, yy, yp, rr, check_matrix_, udata, tmp1, tmp2, tmp3);
        new_entries_ = saved_new;
        for (sunindextype j = 0; j < n_; ++j) {
          const sunindextype i1 = SUNMAX(0, j - mu_);
          const sunindextype i2 = SUNMIN(j + ml_, n_ - 1);
          for (sunindextype i = i1; i <= i2; ++i) {
            const sunrealtype a = SM_ELEMENT_B(Jac, i, j);
            const sunrealtype b = SM_ELEMENT_B(check_matrix_, i, j);
            if (a == b) continue;
            ++num_check_mismatches_;
            sunrealtype col_max = 0.0;
            for (sunindextype k = i1; k <= i2; ++k)
              col_max = SUNMAX(col_max, SUNRabs(SM_ELEMENT_B(check_matrix_, k, j)));
            const double rel = col_max > 0.0 ? SUNRabs(a - b) / col_max : 0.0;
            max_check_rel_ = SUNMAX(max_check_rel_, rel);
            if (rel > 1e-6) {
              char buf[256];
              std::snprintf(buf, sizeof(buf), "row %s@%d col %s@%d coloured=%.3e band=%.3e rel=%.2e",
                            var_names_[var_[i]].c_str(), cell_[i], var_names_[var_[j]].c_str(), cell_[j],
                            a, b, rel);
              if (check_report_.size() < 60) check_report_.push_back(buf);
            }
          }
        }
      }
      return retval;
    }

    // Band evaluation: learn the pattern from it, or audit the pattern.
    if (learn_budget_ > 0) --learn_budget_;
    if (audit) {
      ++num_audits_;
      evals_since_audit_ = 0;
    }
    const long misses_before = num_misses_;
    new_entries_ = 0;
    const int retval = EvaluateGroups(band_groups_, true, tt, c_j, yy, yp, rr, Jac, udata, tmp1, tmp2, tmp3);
    if (audit) num_misses_ = misses_before + new_entries_;
    if (new_entries_ > 0 || colour_groups_.empty()) Recolour();
    return retval;
  }

  void BandedJacobian::ProbeGenericStates(sunrealtype tt, sunrealtype c_j, N_Vector yy, N_Vector yp,
                                         UserData* udata, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    constexpr int kVelocity = 5;            // equation index of the face velocity
    constexpr sunrealtype kRelJitter = 1e-2;
    constexpr sunrealtype kVelocityJitter = 0.05;  // [m/s], applied with random sign

    sunrealtype hh = 0.0;
    IDAGetCurrentStep(ida_mem_, &hh);
    IDAGetErrWeights(ida_mem_, ewt_);
    const sunrealtype* ewt_data = N_VGetArrayPointer(ewt_);
    const sunrealtype* y_data = N_VGetArrayPointer(yy);
    const sunrealtype* yp_data = N_VGetArrayPointer(yp);
    sunrealtype* ys = N_VGetArrayPointer(probe_vecs_[0]);
    sunrealtype* yps = N_VGetArrayPointer(probe_vecs_[1]);
    sunrealtype* rs = N_VGetArrayPointer(probe_vecs_[2]);
    sunrealtype* rtemp = N_VGetArrayPointer(tmp1);
    sunrealtype* ytemp = N_VGetArrayPointer(tmp2);
    sunrealtype* yptemp = N_VGetArrayPointer(tmp3);
    const sunrealtype srur = SUNRsqrt(SUN_UNIT_ROUNDOFF);

    std::mt19937_64 rng(0x5eed + num_res_evals_);
    std::uniform_real_distribution<sunrealtype> jitter(-1.0, 1.0);

    for (int s = 0; s < kProbeStates; ++s) {
      for (sunindextype j = 0; j < n_; ++j) {
        const sunrealtype abs_scale = SUN_RCONST(1.0) / ewt_data[j];
        if (var_[j] == kVelocity) {
          ys[j] = y_data[j] + (SUNRabs(y_data[j]) + kVelocityJitter) * jitter(rng);
        } else {
          ys[j] = y_data[j] + (kRelJitter * SUNRabs(y_data[j]) + abs_scale) * jitter(rng);
        }
        yps[j] = yp_data[j] + (kRelJitter * SUNRabs(yp_data[j]) + abs_scale) * jitter(rng);
      }

      // Residual evaluation at an arbitrary state may fail; such a probe is
      // simply skipped (the learning windows and audits still apply).
      try {
        if (UpdateResiduals(tt, probe_vecs_[0], probe_vecs_[1], probe_vecs_[2], udata) != 0) continue;
        ++num_res_evals_;
        std::copy(ys, ys + n_, ytemp);
        std::copy(yps, yps + n_, yptemp);
        for (const auto& group : band_groups_) {
          for (sunindextype j : group) {
            const sunrealtype inc = Increment(ys[j], yps[j], ewt_data[j], hh, srur);
            ytemp[j] += inc;
            yptemp[j] += c_j * inc;
          }
          if (UpdateResiduals(tt, tmp2, tmp3, tmp1, udata) != 0) break;
          ++num_res_evals_;
          for (sunindextype j : group) {
            ytemp[j] = ys[j];
            yptemp[j] = yps[j];
            const sunindextype i1 = SUNMAX(0, j - mu_);
            const sunindextype i2 = SUNMIN(j + ml_, n_ - 1);
            for (sunindextype i = i1; i <= i2; ++i) {
              if (std::isnan(rs[i])) continue;
              const sunrealtype diff = rtemp[i] - rs[i];
              if (diff != SUN_RCONST(0.0)) Learn(j, i);  // NaN != 0 counts as a coupling
            }
          }
        }
      } catch (...) {
        continue;
      }
    }
    Recolour();
  }

  int BandedJacobian::EvaluateGroups(const std::vector<std::vector<sunindextype>>& groups, bool band,
                                     sunrealtype tt, sunrealtype c_j, N_Vector yy, N_Vector yp, N_Vector rr,
                                     SUNMatrix Jac, UserData* udata,
                                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    sunrealtype hh = 0.0;
    IDAGetCurrentStep(ida_mem_, &hh);
    IDAGetErrWeights(ida_mem_, ewt_);

    N_Vector rtemp = tmp1;
    N_Vector ytemp = tmp2;
    N_Vector yptemp = tmp3;

    const sunrealtype* ewt_data = N_VGetArrayPointer(ewt_);
    const sunrealtype* r_data = N_VGetArrayPointer(rr);
    const sunrealtype* y_data = N_VGetArrayPointer(yy);
    const sunrealtype* yp_data = N_VGetArrayPointer(yp);
    sunrealtype* rtemp_data = N_VGetArrayPointer(rtemp);
    sunrealtype* ytemp_data = N_VGetArrayPointer(ytemp);
    sunrealtype* yptemp_data = N_VGetArrayPointer(yptemp);

    N_VScale(SUN_RCONST(1.0), yy, ytemp);
    N_VScale(SUN_RCONST(1.0), yp, yptemp);

    const sunrealtype srur = SUNRsqrt(SUN_UNIT_ROUNDOFF);
    const bool learn = band && mode_ != JacobianMode::BandDQ;
    int retval = 0;

    for (const auto& group : groups) {
      for (sunindextype j : group) {
        const sunrealtype inc = Increment(y_data[j], yp_data[j], ewt_data[j], hh, srur);
        ytemp_data[j] += inc;
        yptemp_data[j] += c_j * inc;
      }

      retval = UpdateResiduals(tt, ytemp, yptemp, rtemp, udata);
      ++num_res_evals_;
      if (retval != 0) break;

      for (sunindextype j : group) {
        ytemp_data[j] = y_data[j];
        yptemp_data[j] = yp_data[j];
        const sunrealtype inc = Increment(y_data[j], yp_data[j], ewt_data[j], hh, srur);
        const sunrealtype inc_inv = SUN_RCONST(1.0) / inc;
        sunrealtype* col_j = SUNBandMatrix_Column(Jac, j);

        if (band) {
          const sunindextype i1 = SUNMAX(0, j - mu_);
          const sunindextype i2 = SUNMIN(j + ml_, n_ - 1);
          for (sunindextype i = i1; i <= i2; ++i) {
            const sunrealtype value = inc_inv * (rtemp_data[i] - r_data[i]);
            SM_COLUMN_ELEMENT_B(col_j, i, j) = value;
            if (value == SUN_RCONST(0.0)) continue;
            if (learn && Learn(j, i)) ++new_entries_;
            if (record_pattern_)
              ++pattern_counts_[{ var_names_[var_[i]], var_names_[var_[j]], cell_[i] - cell_[j] }];
          }
        } else {
          for (sunindextype k = row_ptr_[j]; k < row_ptr_[j + 1]; ++k) {
            const sunindextype i = rows_[k];
            SM_COLUMN_ELEMENT_B(col_j, i, j) = inc_inv * (rtemp_data[i] - r_data[i]);
          }
        }
      }
    }
    return retval;
  }

  int BandedJacobianFn(sunrealtype tt, sunrealtype c_j, N_Vector yy, N_Vector yp,
                       N_Vector rr, SUNMatrix Jac, void* user_data,
                       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    auto* udata = static_cast<UserData*>(user_data);
    return udata->jacobian_->Evaluate(tt, c_j, yy, yp, rr, Jac, udata, tmp1, tmp2, tmp3);
  }

} // end adsim namespace
