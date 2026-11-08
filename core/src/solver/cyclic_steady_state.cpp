/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/solver/cyclic_steady_state.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

#include <Eigen/Dense>

#include "adsim/solver/cycle.h"
#include "adsim/solver/log.h"
#include "adsim/solver/kpis.h"
#include "adsim/objects/reactor.h"


namespace adsim {

  bool CyclicSteadyState::IsAtCSS(const EquationLog& logger, int cycle_index) {
    return CheckCSS(logger, cycle_index).converged;
  }

  int CyclicSteadyState::CyclesToLog() const {
    if (method_ == CSSMethod::NONE) return 0;
    const int k = std::max(1, consecutive_cycles_);
    // Every single-cycle test except the balance-only one differences the cycle
    // against the previous one, so the window of k tests spans k + 1 cycles.
    const bool balance_only = (method_ == CSSMethod::MASS_BALANCE_TOTAL ||
                               method_ == CSSMethod::MASS_BALANCE_COMPONENT) &&
                              !(state_tolerance_ > 0.0);
    return balance_only ? k : k + 1;
  }

  CSSResult CyclicSteadyState::CheckCycle(const EquationLog& logger, int cycle_index) {
    switch (method_) {
      case CSSMethod::NONE:
        return {false, 0.0, 0.0, "None"};
      case CSSMethod::MOLAR_PROFILE:
        return CheckProfileResult(logger, cycle_index);
      case CSSMethod::MASS_BALANCE_TOTAL:
        return CheckMassBalanceTotalResult(logger, cycle_index);
      case CSSMethod::MASS_BALANCE_COMPONENT:
        return CheckMassBalanceComponentResult(logger, cycle_index);
    }
    return {false, 0.0, 0.0, "Unknown"};
  }

  bool CyclicSteadyState::MassBalanceCyclePasses(double balance_error, double balance_tolerance,
                                                 double state_change, double state_tolerance) {
    if (!(balance_error <= balance_tolerance)) return false;
    if (!(state_tolerance > 0.0)) return true;  // state test off: the balance decides alone
    return state_change >= 0.0 && state_change <= state_tolerance;
  }

  int CyclicSteadyState::TrailingPasses(int cycle_index, int required,
                                        const std::function<bool(int)>& passes) {
    int count = 0;
    for (int c = cycle_index; c >= 0 && count < required; --c) {
      if (!passes(c)) break;
      ++count;
    }
    return count;
  }

  CSSResult CyclicSteadyState::CheckCSS(const EquationLog& logger, int cycle_index,
                                        const JumpPredicate& started_from_jump) {
    CSSResult result = CheckCycle(logger, cycle_index);
    if (method_ == CSSMethod::NONE) return result;

    // The diagnostic reads the log only; it does not enter the rule below.
    const CSSDistanceEstimate d = EstimateDistance(logger, cycle_index, started_from_jump);
    result.distance_estimated = d.estimated;
    result.convergence_rate = d.rate;
    result.distance = d.distance;
    result.distance_ratio = d.ratio;
    result.distance_note = d.note;

    const int required = std::max(1, consecutive_cycles_);
    result.required_cycles = required;

    // Count the consecutive passes ending at this cycle. The single-cycle test is
    // a pure function of the log, so the earlier cycles of the window are simply
    // re-evaluated; a cycle no longer in the log fails.
    const bool passed_here = result.converged;
    result.passing_cycles = TrailingPasses(cycle_index, required, [&](int c) {
      return c == cycle_index ? passed_here : CheckCycle(logger, c).converged;
    });

    // No method may declare convergence before min_cycles_ cycles have run.
    result.converged = result.passing_cycles >= required && cycle_index >= min_cycles_ - 1;
    return result;
  }

  CSSDistanceEstimate CyclicSteadyState::EstimateDistance(const std::vector<double>& changes,
                                                          double tolerance) {
    CSSDistanceEstimate out;
    const int n = static_cast<int>(changes.size());
    auto usable = [](double v) { return std::isfinite(v) && v > 0.0; };
    // Longest trailing run of strictly decreasing positive changes, capped.
    int len = 0;
    if (n > 0 && usable(changes[n - 1])) {
      len = 1;
      for (int i = n - 2; i >= 0 && len < kDistanceMaxChanges; --i) {
        if (!usable(changes[i]) || !(changes[i] > changes[i + 1])) break;
        ++len;
      }
    }
    out.changes_used = len;
    if (len < kDistanceMinChanges) {
      out.note = n < kDistanceMinChanges ? "history too short" : "not a monotone decay";
      return out;
    }
    const double first = changes[n - len];
    const double last = changes[n - 1];
    double rho = std::pow(last / first, 1.0 / static_cast<double>(len - 1));
    rho = std::clamp(rho, 0.0, kDistanceMaxRate);
    out.estimated = true;
    out.rate = rho;
    out.distance = last * rho / (1.0 - rho);
    if (tolerance > 0.0) out.ratio = out.distance / tolerance;
    return out;
  }

  CSSDistanceEstimate CyclicSteadyState::EstimateDistance(const EquationLog& logger, int cycle_index,
                                                          const JumpPredicate& started_from_jump) const {
    // Collect the changes back from cycle_index, newest first, while they are
    // logged and on the natural trajectory.
    std::vector<double> newest_first;
    bool hit_jump = false;
    for (int c = cycle_index; c >= 1 && static_cast<int>(newest_first.size()) < kDistanceMaxChanges; --c) {
      if (started_from_jump && started_from_jump(c)) { hit_jump = true; break; }
      const double change = EndOfCycleStateChange(logger, c);
      if (change < 0.0) break;  // not in the log
      newest_first.push_back(change);
    }
    std::vector<double> changes(newest_first.rbegin(), newest_first.rend());
    const double tolerance = state_tolerance_ > 0.0 ? static_cast<double>(state_tolerance_) : 0.0;
    CSSDistanceEstimate out = EstimateDistance(changes, tolerance);
    if (!out.estimated && hit_jump && static_cast<int>(changes.size()) < kDistanceMinChanges)
      out.note = "too few cycles since the last acceleration jump";
    return out;
  }

  double CyclicSteadyState::EndOfCycleStateChange(const EquationLog& logger, int cycle_index) {
    if (cycle_index < 1) return -1.0;
    const int last_step = static_cast<int>(logger.sys_data.step_list_.size()) - 1;
    if (last_step < 0) return -1.0;
    const Log& curr = logger.GetLog(cycle_index, last_step);
    const Log& prev = logger.GetLog(cycle_index - 1, last_step);
    if (curr.step_data.times.empty() || prev.step_data.times.empty()) return -1.0;

    // Last logged row (the end of the step) of one variable, or nullptr.
    auto last_row = [](const VariableData* v) -> const std::vector<double>* {
      return (v && !v->values.empty()) ? &v->values.back() : nullptr;
    };
    auto bulk = [&](const Log& log, VariableType type) -> const VariableData* {
      auto it = log.variables.find(type);
      return it == log.variables.end() ? nullptr : &it->second;
    };
    auto component = [&](const Log& log, ComponentVariableType type,
                         const std::string& name) -> const VariableData* {
      auto c = log.component_variables.find(name);
      if (c == log.component_variables.end()) return nullptr;
      auto it = c->second.find(type);
      return it == c->second.end() ? nullptr : &it->second;
    };
    // RMS over the finite cells of a row, and of the difference of two rows.
    auto rms = [](const std::vector<double>& x) {
      double s = 0.0; size_t n = 0;
      for (double v : x) if (std::isfinite(v)) { s += v * v; ++n; }
      return n ? std::sqrt(s / static_cast<double>(n)) : 0.0;
    };
    auto rms_diff = [](const std::vector<double>& a, const std::vector<double>& b) {
      double s = 0.0; size_t n = 0;
      const size_t m = std::min(a.size(), b.size());
      for (size_t i = 0; i < m; ++i) {
        const double d = a[i] - b[i];
        if (std::isfinite(d)) { s += d * d; ++n; }
      }
      return n ? std::sqrt(s / static_cast<double>(n)) : 0.0;
    };

    double change = 0.0;
    auto add = [&](const std::vector<double>* a, const std::vector<double>* b, double scale) {
      if (!a || !b || a->size() != b->size()) return;
      if (!(scale > 1e-300)) return;  // an all-zero field (e.g. no loading) cannot drift
      change = std::max(change, rms_diff(*a, *b) / scale);
    };

    for (VariableType type : {VariableType::TEMPERATURE, VariableType::SOLID_TEMPERATURE,
                              VariableType::WALL_TEMPERATURE}) {
      const auto* a = last_row(bulk(curr, type));
      const auto* b = last_row(bulk(prev, type));
      if (a) add(a, b, rms(*a));
    }

    const auto& components = logger.sys_data.components_;
    std::vector<double> total_loading;
    for (const auto& name : components) {
      const auto* q = last_row(component(curr, ComponentVariableType::ADSORBED_CONCENTRATION, name));
      if (!q) continue;
      if (total_loading.empty()) total_loading.assign(q->size(), 0.0);
      for (size_t i = 0; i < q->size() && i < total_loading.size(); ++i) total_loading[i] += (*q)[i];
    }
    // The scale is floored at the integrator's absolute tolerance on the
    // loading unknowns (1e-8 mol/kg, set in Equations::Initialise): a loading
    // below it is not resolved, and an inert packing (no isotherm) holds only
    // round-off (1e-26 to 1e-22 mol/kg, of either sign), which divided by its own RMS
    // read as an O(1) change and kept the state test from ever passing.
    constexpr double kLoadingAbsTol = 1e-8;  // [mol/kg]
    const double loading_scale = std::max(rms(total_loading), kLoadingAbsTol);
    for (const auto& name : components) {
      add(last_row(component(curr, ComponentVariableType::ADSORBED_CONCENTRATION, name)),
          last_row(component(prev, ComponentVariableType::ADSORBED_CONCENTRATION, name)), loading_scale);
      add(last_row(component(curr, ComponentVariableType::MOLE_FRACTION, name)),
          last_row(component(prev, ComponentVariableType::MOLE_FRACTION, name)), 1.0);
      // The slow-site loading of a dual-rate LDF component: a state of its
      // own, which can still drift while the total loading looks settled.
      // Not logged (no-op) for decks without one.
      add(last_row(component(curr, ComponentVariableType::SLOW_ADSORBED_CONCENTRATION, name)),
          last_row(component(prev, ComponentVariableType::SLOW_ADSORBED_CONCENTRATION, name)), loading_scale);
    }
    return change;
  }

  void CyclicSteadyState::ApplyStateTest(const EquationLog& logger, int cycle_index, CSSResult& result) {
    if (state_tolerance_ > 0.0) {
      result.state_tolerance = state_tolerance_;
      result.state_change = EndOfCycleStateChange(logger, cycle_index);
    }
    result.converged = MassBalanceCyclePasses(result.error, result.tolerance,
                                              result.state_change, result.state_tolerance);
    if (result.state_change < 0.0)  // no previous cycle to compare with
      result.state_change = std::numeric_limits<double>::quiet_NaN();
  }

  CSSResult CyclicSteadyState::CheckProfileResult(const EquationLog& logger, int cycle_index) {
    CSSResult result;
    result.method_name = "Molar Profile";

    // This method differences consecutive cycles. (The minimum-cycles floor is
    // applied by CheckCSS.)
    result.tolerance = molar_absolute_tolerance;
    if (cycle_index < 1)
      return result;

    auto curr_molar_conc_data = logger.GetVariable(cycle_index, 0, VariableType::MOLAR_CONCENTRATION);
    auto prev_molar_conc_data = logger.GetVariable(cycle_index - 1, 0, VariableType::MOLAR_CONCENTRATION);

    const auto& curr_molar_conc_timeseries = curr_molar_conc_data.values;
    const auto& prev_molar_conc_timeseries = prev_molar_conc_data.values;

    if (curr_molar_conc_timeseries.empty() || prev_molar_conc_timeseries.empty())
      return result;

    const auto& curr_molar_conc_profile = curr_molar_conc_timeseries.back();
    const auto& prev_molar_conc_profile = prev_molar_conc_timeseries.back();

    if (curr_molar_conc_profile.size() != prev_molar_conc_profile.size())
      return result;

    Eigen::VectorXd curr_profile = Eigen::Map<const Eigen::VectorXd>(curr_molar_conc_profile.data(), curr_molar_conc_profile.size());
    Eigen::VectorXd prev_profile = Eigen::Map<const Eigen::VectorXd>(prev_molar_conc_profile.data(), prev_molar_conc_profile.size());

    result.error = (curr_profile - prev_profile).norm();
    result.tolerance = molar_absolute_tolerance + molar_relative_tolerance * curr_profile.norm();
    result.converged = (result.error <= result.tolerance);

    return result;
  }

  CSSResult CyclicSteadyState::CheckMassBalanceTotalResult(const EquationLog& logger, int cycle_index) {
    CSSResult result;
    result.method_name = "Mass Balance (Total)";
    result.tolerance = mass_balance_tolerance_;

    const auto& step_indices = logger.sys_data.step_index_set_;

    double total_moles_in = 0.0;
    double total_moles_out = 0.0;

    for (int step_index : step_indices) {
      double bottom = GetTotalMolesIn(logger, cycle_index, step_index, BoundaryLocation::BOTTOM);
      double top = GetTotalMolesIn(logger, cycle_index, step_index, BoundaryLocation::TOP);

      if (bottom > 0.0) total_moles_in += bottom;
      else total_moles_out += std::abs(bottom);

      if (top > 0.0) total_moles_in += top;
      else total_moles_out += std::abs(top);
    }

    if (total_moles_in < 1e-30) return result;

    result.error = std::abs(total_moles_in - total_moles_out) / total_moles_in;
    result.converged = (result.error <= result.tolerance);
    ApplyStateTest(logger, cycle_index, result);
    return result;
  }

  CSSResult CyclicSteadyState::CheckMassBalanceComponentResult(const EquationLog& logger, int cycle_index) {
    CSSResult result;
    result.method_name = "Mass Balance (Component)";
    result.tolerance = mass_balance_tolerance_;

    const auto& components = logger.sys_data.components_;
    const auto& step_indices = logger.sys_data.step_index_set_;

    double max_error = 0.0;

    for (const auto& component : components) {
      double total_moles_in = 0.0;
      double total_moles_out = 0.0;

      for (int step_index : step_indices) {
        double bottom = GetTotalComponentMolesIn(logger, cycle_index, step_index, BoundaryLocation::BOTTOM, component);
        double top = GetTotalComponentMolesIn(logger, cycle_index, step_index, BoundaryLocation::TOP, component);

        if (bottom > 0.0) total_moles_in += bottom;
        else total_moles_out += std::abs(bottom);

        if (top > 0.0) total_moles_in += top;
        else total_moles_out += std::abs(top);
      }

      // A component that is never fed (trace species, inert, product-only) has
      // no meaningful relative error; skip it rather than abandoning the whole
      // check, which would leave result.converged false for every cycle.
      if (total_moles_in < 1e-30) continue;

      double comp_error = std::abs(total_moles_in - total_moles_out) / total_moles_in;
      max_error = std::max(max_error, comp_error);
    }

    result.error = max_error;
    result.converged = (result.error <= result.tolerance);
    ApplyStateTest(logger, cycle_index, result);
    return result;
  }

  namespace {
    // Which fields an acceleration is allowed to move. We extrapolate the SLOW
    // differential fields -- adsorbed loadings, solid and wall temperature, gas
    // composition -- but deliberately leave the total gas concentration and the
    // fluid temperature untouched, and skip the algebraic fields (pressure,
    // velocity) entirely.
    //
    // The gas concentration and fluid temperature set the pressure through the
    // equation of state (p = c R T), and the pressure and velocity are solved
    // *algebraically*. Re-seeding via IDACalcIC(IDA_YA_YDP_INIT) recomputes the
    // algebraic fields from the (fixed) differential ones; if the extrapolation
    // has perturbed c or the fluid temperature, that consistency solve cannot
    // reconcile them with the momentum balance and the run aborts. Those two
    // fields are also FAST -- they re-equilibrate to the re-seeded slow state
    // within the first step of the next cycle -- so leaving them out costs almost
    // nothing while keeping the seed integrable. The slow fields carry the
    // cycle-to-cycle memory that actually governs the approach to CSS.
    bool IsExtrapolable(const Equations& state, const Equation* eq) {
      return eq->GetEquationType() == EquationType::Differential
          && eq != &state.concentration
          && eq != &state.fluid_temperature;
    }

    // Mole fractions must remain a valid composition after a jump: renormalise
    // them to sum to one in every cell, or the species balance seeds an
    // impossible gas mixture. (Negative values are clipped where the jump is
    // applied.)
    void RenormaliseMoleFractions(Equations& state) {
      if (state.mole_fraction.empty()) return;
      const size_t n_cells = state.mole_fraction.begin()->second.var.size();
      for (size_t c = 0; c < n_cells; ++c) {
        double sum = 0.0;
        for (auto& [component, eq] : state.mole_fraction) sum += eq.var[c];
        if (sum > 1e-30) {
          for (auto& [component, eq] : state.mole_fraction) eq.var[c] /= sum;
        }
      }
    }

    std::string FormatNumber(double x, int digits = 3) {
      std::ostringstream ss;
      ss << std::setprecision(digits) << x;
      return ss.str();
    }
  }

  bool AccelerateCyclicSteadyState::Accelerate(Equations& state, int cycle_index, AccelerationNote* note) {
    switch (method_) {
      case AccelerationMethod::AITKEN: return AccelerateAitken(state, cycle_index, note);
      case AccelerationMethod::ANDERSON: return AccelerateAnderson(state, cycle_index, note);
      case AccelerationMethod::NONE: break;
    }
    return false;
  }

  bool AccelerateCyclicSteadyState::AccelerateAitken(Equations& state, int cycle_index, AccelerationNote* note) {
    // A cycle index that does not advance means a new run has started; drop the
    // history collected for the previous one.
    if (cycle_index <= last_cycle_index_) history_.clear();
    last_cycle_index_ = cycle_index;

    // Skip the initial transient: the extrapolation assumes near-geometric decay
    // of the slowest error mode, which only holds once the start-up front has
    // washed through.
    if (cycle_index < 2) return false;

    const std::vector<Equation*>& equations = state.GetEquations();

    // The fields the extrapolation moves (see IsExtrapolable).
    auto extrapolable = [&](const Equation* eq) { return IsExtrapolable(state, eq); };

    // Snapshot the current (natural) end-of-cycle state -- one profile per
    // extrapolable equation, in GetEquations() order (an empty placeholder keeps
    // the indexing aligned for the skipped fields).
    std::vector<std::vector<double>> snapshot;
    snapshot.reserve(equations.size());
    for (const Equation* eq : equations) {
      if (extrapolable(eq)) {
        snapshot.push_back(eq->var);
      } else {
        snapshot.emplace_back();  // placeholder to keep indices aligned
      }
    }
    history_.push_back(std::move(snapshot));

    // Need three successive natural end-of-cycle states before extrapolating.
    if (history_.size() < 3) return false;

    const auto& s0 = history_[0];
    const auto& s1 = history_[1];
    const auto& s2 = history_[2];

    // Aitken Delta^2 extrapolation of the extrapolable differential fields toward
    // the fixed point of the cycle-to-cycle map (the scalar-weighted vector form of
    // the Graves-Morris vector-epsilon method).
    // REF: https://psecommunity.org/wp-content/plugins/wpor/includes/file/2408/LAPSE-2024.1530-1v2.pdf
    //
    // All moved fields share a SINGLE coefficient mu, i.e. each is advanced by the
    // same multiple of its own natural cycle-to-cycle step r2 = x2 - x1. This is
    // essential: the natural step is a self-consistent direction (it is an actual
    // segment of the cycle trajectory), so scaling it keeps the moved fields in
    // mutual balance -- the loading stays consistent with the composition and solid
    // temperature it was coupled to. Giving each field its own coefficient (or
    // moving the loading alone) breaks that balance and seeds a stiff, often
    // non-integrable state.
    //
    // mu is computed from magnitude-scaled inner products so that no single field
    // (temperature ~ 300 vs a mole fraction ~ 0.2) dominates the fit:
    //     mu = sum_f w_f (r2_f . d_f) / sum_f w_f (d_f . d_f),   w_f = 1 / scale_f^2
    // with scale_f the RMS magnitude of field f. Applied per field as
    // x*_f = x2_f - mu r2_f, this is exact for a single coupled geometric mode.
    double num = 0.0, den = 0.0;
    double scaled_r1_sq = 0.0, scaled_r2_sq = 0.0, scaled_r1r2 = 0.0;
    for (size_t i = 0; i < equations.size(); ++i) {
      if (!extrapolable(equations[i])) continue;
      const auto& x0 = s0[i]; const auto& x1 = s1[i]; const auto& x2 = s2[i];
      if (x0.size() != x1.size() || x1.size() != x2.size() || x2.empty()) continue;

      double mag_sq = 0.0, r1_sq = 0.0, r2_sq = 0.0, r1r2 = 0.0, r2d = 0.0, dd = 0.0;
      for (size_t c = 0; c < x2.size(); ++c) {
        const double r1 = x1[c] - x0[c];
        const double r2 = x2[c] - x1[c];
        const double d  = r2 - r1;
        mag_sq += x2[c] * x2[c];
        r1_sq  += r1 * r1;
        r2_sq  += r2 * r2;
        r1r2   += r1 * r2;
        r2d    += r2 * d;
        dd     += d * d;
      }
      const double scale_sq = std::max(mag_sq / static_cast<double>(x2.size()), 1e-30);
      const double w = 1.0 / scale_sq;
      num += w * r2d;
      den += w * dd;
      scaled_r1_sq += w * r1_sq;
      scaled_r2_sq += w * r2_sq;
      scaled_r1r2  += w * r1r2;
    }

    // Every full window is consumed whether or not it is used, so the next
    // window is three fresh cycles. (s0, s1, s2 refer into history_: clear it
    // only on the way out.)
    auto skip = [&](const std::string& reason) {
      history_.clear();
      if (note) note->skipped = reason;
      return false;
    };
    auto fmt = [](double x, int digits = 3) { return FormatNumber(x, digits); };

    if (scaled_r2_sq < 1e-300) return skip("the state no longer changes");
    if (den < 1e-300) return skip("no curvature to fit");

    // The scaled step must be shrinking, or the sequence is not contracting
    // toward a fixed point and there is nothing to extrapolate to.
    const double step_ratio = std::sqrt(scaled_r2_sq / scaled_r1_sq);
    if (!(step_ratio < 1.0)) {
      return skip("steps not shrinking: step ratio " + fmt(step_ratio));
    }

    double mu = num / den;

    // The decay rate a single geometric mode would have, from the projection of
    // the second step on the first: rate = (r1 . r2) / (r1 . r1). It lies below
    // the step ratio, so below 1. For one mode of rate rho, r2 = rho r1 and
    // mu = rho / (rho - 1): negative for monotone decay (0 < rho < 1), when the
    // jump continues along the steps, and positive only for an oscillating mode
    // (rho < 0). A coefficient of the other sign means the three samples are not
    // one geometric mode -- typically a front that drifts by the same amount
    // every cycle, where d = r2 - r1 is mostly noise -- and a jump would move
    // the state BACK along its own trajectory (undoing natural progress). Skip.
    const double rate = scaled_r1r2 / scaled_r1_sq;
    if (note) {
      note->fit = "mu = " + fmt(mu) + ", step ratio " + fmt(step_ratio, 4)
                + ", single-mode rate " + fmt(rate, 4);
    }
    //
    // A positive mu belongs to an oscillating mode, rho < 0, for which
    // mu = rho / (rho - 1) lies in (0, 1/2) while |rho| < 1. A positive mu of 1
    // or more has no contracting single mode behind it, so it is skipped too,
    // never limited: the limit below applies to negative mu only. (With the
    // weighted inner product used throughout, mu >= 1/2 is the same as
    // |r2| >= |r1|, so after the step-ratio test this guard only catches
    // rounding; it states the rule rather than relying on that identity.)
    if ((rate > 0.0) != (mu < 0.0) || mu >= 1.0) {
      return skip("steps do not follow one geometric mode: mu = " + fmt(mu)
                  + ", single-mode rate " + fmt(rate));
    }

    // Limit the jump. |mu| grows without bound as the rate approaches 1
    // (|mu| = 5 at a rate of 0.83, 20 at 0.95, 28 at 0.965, 140 at 0.993).
    // Such slow modes are where acceleration matters most, so the jump is
    // SHORTENED to kMaxStepMultiple natural steps rather than dropped (an
    // earlier version skipped every window with |mu| > 5, which stopped
    // acceleration altogether once the run settled on its slow mode). Along a
    // contracting mode of rate rho, a jump of m < |mu| steps multiplies the
    // remaining error by 1 - m (1 - rho) / rho, which is between 0 and 1: it
    // undershoots, and never jumps past the fixed point. The limit bounds what
    // the jump does to a faster mode of rate rho_f still in the samples
    // (multiplied by 1 - m (1 - rho_f) / rho_f). 20 was chosen on measured
    // decks: up to 20 the KPIs at the CSS stop stayed within the plain run's
    // stopping error; at 50 one deck (Perez 2019 F, 10 cells) stopped 1.6x
    // further from its CSS than the plain run.
    if (!std::isfinite(mu)) return skip("no curvature to fit");
    if (mu < -kMaxStepMultiple) {
      if (note) note->fit += " (jump limited to mu = " + fmt(-kMaxStepMultiple) + ")";
      mu = -kMaxStepMultiple;
    }

    // Apply x*_f = x2_f - mu r2_f per field, then PROJECT back onto the physical
    // manifold: a jump can still nudge a loading negative or a mole fraction out
    // of [0, 1], and IDA cannot integrate such a seed. All extrapolated fields
    // here (temperatures, mole fractions, loadings) are non-negative.
    for (size_t i = 0; i < equations.size(); ++i) {
      if (!extrapolable(equations[i])) continue;
      const auto& x1 = s1[i]; const auto& x2 = s2[i];
      std::vector<double>& var = equations[i]->var;
      for (size_t c = 0; c < var.size(); ++c) {
        const double r2 = x2[c] - x1[c];
        var[c] = std::max(x2[c] - mu * r2, 0.0);
      }
    }

    RenormaliseMoleFractions(state);

    // Restart: the three samples are consumed, so the next triple is collected
    // fresh from the (post-jump) trajectory. This keeps the geometric-mode fit
    // clean and spaces the jumps ~3 cycles apart, letting the bed relax between
    // them.
    history_.clear();

    return true;
  }

  bool AccelerateCyclicSteadyState::AccelerateAnderson(Equations& state, int cycle_index, AccelerationNote* note) {
    // A cycle index that does not advance means a new run has started.
    if (cycle_index <= last_cycle_index_) {
      pairs_.clear();
      previous_end_.clear();
      previous_jumped_ = false;
    }
    last_cycle_index_ = cycle_index;

    const std::vector<Equation*>& equations = state.GetEquations();

    // The natural end state of this cycle, the slow fields flattened in
    // GetEquations() order, and each field's range in the flat vector.
    struct Range { size_t eq, begin, end; };
    std::vector<Range> ranges;
    std::vector<double> end_state;
    for (size_t i = 0; i < equations.size(); ++i) {
      if (!IsExtrapolable(state, equations[i])) continue;
      const auto& var = equations[i]->var;
      ranges.push_back({i, end_state.size(), end_state.size() + var.size()});
      end_state.insert(end_state.end(), var.begin(), var.end());
    }

    // A pair (x, G(x)) needs this cycle's start state x exactly: the previous
    // cycle's natural end state, when no jump came between. Skip the start-up
    // transient (cycles 0 and 1) as Aitken does: the cycle map is far from
    // linear there.
    const bool layout_changed = !previous_end_.empty() && previous_end_.size() != end_state.size();
    if (layout_changed) pairs_.clear();
    const bool new_pair = cycle_index >= 2 && !previous_jumped_ && !layout_changed && !previous_end_.empty();
    if (new_pair) {
      pairs_.push_back({previous_end_, end_state});
      if (pairs_.size() > static_cast<size_t>(kAndersonDepth) + 1) pairs_.erase(pairs_.begin());
    }
    previous_end_ = end_state;
    previous_jumped_ = false;

    // A jump starts from the newest pair's end state, which must be the state
    // the solver holds now, and needs at least one earlier pair.
    if (!new_pair || pairs_.size() < 2) return false;

    auto skip = [&](const std::string& reason) {
      if (note) note->skipped = reason;
      return false;
    };

    // Field-scaled norm, as in Aitken: each field weighted by 1 / RMS magnitude,
    // so a temperature (~300) does not outweigh a mole fraction (~0.2).
    std::vector<double> weight(end_state.size());
    for (const Range& r : ranges) {
      double mag_sq = 0.0;
      for (size_t k = r.begin; k < r.end; ++k) mag_sq += end_state[k] * end_state[k];
      const size_t n = r.end - r.begin;
      const double scale = std::sqrt(std::max(n > 0 ? mag_sq / static_cast<double>(n) : 0.0, 1e-30));
      for (size_t k = r.begin; k < r.end; ++k) weight[k] = 1.0 / scale;
    }

    const AndersonPair& newest = pairs_.back();
    const Eigen::Index n_rows = static_cast<Eigen::Index>(end_state.size());
    Eigen::VectorXd f_n(n_rows);
    for (Eigen::Index k = 0; k < n_rows; ++k) f_n[k] = weight[k] * (newest.g[k] - newest.x[k]);
    const double residual = f_n.norm();
    // Below 1e-12 RMS (relative to each field's magnitude) the steps are
    // rounding: there is nothing left to fit.
    if (!(residual > 1e-12 * std::sqrt(static_cast<double>(n_rows))))
      return skip("the state no longer changes");

    // Columns: differences of the newest residual (and end state) with each
    // earlier pair's, newest first, so a rank-deficient fit drops the oldest.
    const int n_cols = static_cast<int>(pairs_.size()) - 1;
    Eigen::MatrixXd dF(n_rows, n_cols);
    for (int j = 0; j < n_cols; ++j) {
      const AndersonPair& older = pairs_[pairs_.size() - 2 - j];
      for (Eigen::Index k = 0; k < n_rows; ++k)
        dF(k, j) = weight[k] * ((newest.g[k] - newest.x[k]) - (older.g[k] - older.x[k]));
    }

    // Least squares by column-pivoted QR. Columns that are (nearly) linear
    // combinations of others carry no new information: drop the oldest until
    // the fit has full rank.
    int used = n_cols;
    Eigen::VectorXd gamma;
    while (used > 0) {
      Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(dF.leftCols(used));
      qr.setThreshold(1e-10);
      if (qr.rank() == used) {
        gamma = qr.solve(f_n);
        break;
      }
      --used;
    }
    if (used == 0 || !gamma.allFinite()) return skip("no independent steps to combine");

    const double fit_residual = (f_n - dF.leftCols(used) * gamma).norm();
    const double ratio = fit_residual / residual;

    // The jump from the newest end state: -sum_j gamma_j (g_n - g_j), and its
    // size in natural steps (the newest step g_n - x_n, in the same norm).
    std::vector<double> jump(end_state.size(), 0.0);
    for (int j = 0; j < used; ++j) {
      const AndersonPair& older = pairs_[pairs_.size() - 2 - j];
      for (size_t k = 0; k < jump.size(); ++k) jump[k] -= gamma[j] * (newest.g[k] - older.g[k]);
    }
    double jump_norm = 0.0;
    for (size_t k = 0; k < jump.size(); ++k) jump_norm += weight[k] * weight[k] * jump[k] * jump[k];
    jump_norm = std::sqrt(jump_norm);
    double multiple = jump_norm / residual;

    if (note) {
      note->fit = "Anderson depth " + std::to_string(used) + ", residual ratio " + FormatNumber(ratio, 4)
                + ", jump " + FormatNumber(multiple) + " natural steps";
    }
    if (!(ratio < kAndersonMaxResidualRatio)) {
      return skip("the fit does not reduce the residual: ratio " + FormatNumber(ratio));
    }
    if (!(jump_norm > 0.0) || !std::isfinite(jump_norm)) return skip("no jump to make");

    // Shorten, never drop, a jump longer than kMaxStepMultiple natural steps
    // (as for Aitken: it then closes part of the gap along the same direction).
    double factor = 1.0;
    if (multiple > kMaxStepMultiple) {
      factor = kMaxStepMultiple / multiple;
      if (note) note->fit += " (limited to " + FormatNumber(kMaxStepMultiple) + ")";
    }

    // Apply from the newest end state (the state the solver holds), clip at
    // zero and renormalise the mole fractions, as for Aitken.
    for (const Range& r : ranges) {
      std::vector<double>& var = equations[r.eq]->var;
      for (size_t k = r.begin; k < r.end; ++k)
        var[k - r.begin] = std::max(newest.g[k] + factor * jump[k], 0.0);
    }
    RenormaliseMoleFractions(state);

    previous_jumped_ = true;
    return true;
  }

} // end adsim namespace
