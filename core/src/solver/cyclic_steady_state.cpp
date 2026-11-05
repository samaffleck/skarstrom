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

  CSSResult CyclicSteadyState::CheckCSS(const EquationLog& logger, int cycle_index) {
    CSSResult result = CheckCycle(logger, cycle_index);
    if (method_ == CSSMethod::NONE) return result;

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

  bool AccelerateCyclicSteadyState::Accelerate(Equations& state, int cycle_index, AccelerationNote* note) {
    if (method_ == AccelerationMethod::NONE) return false;
    // Only the Aitken/Graves-Morris scheme is implemented; it is the body below.
    // Future methods would branch here on method_.

    // A cycle index that does not advance means a new run has started; drop the
    // history collected for the previous one.
    if (cycle_index <= last_cycle_index_) history_.clear();
    last_cycle_index_ = cycle_index;

    // Skip the initial transient: the extrapolation assumes near-geometric decay
    // of the slowest error mode, which only holds once the start-up front has
    // washed through.
    if (cycle_index < 2) return false;

    const std::vector<Equation*>& equations = state.GetEquations();

    // Which fields the extrapolation is allowed to move. We extrapolate the SLOW
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
    auto extrapolable = [&](const Equation* eq) {
      return eq->GetEquationType() == EquationType::Differential
          && eq != &state.concentration
          && eq != &state.fluid_temperature;
    };

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
    auto fmt = [](double x, int digits = 3) {
      std::ostringstream ss;
      ss << std::setprecision(digits) << x;
      return ss.str();
    };

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
    if ((rate > 0.0) != (mu < 0.0)) {
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
    if (std::abs(mu) > kMaxStepMultiple) {
      const double limited = mu < 0.0 ? -kMaxStepMultiple : kMaxStepMultiple;
      if (note) note->fit += " (jump limited to mu = " + fmt(limited) + ")";
      mu = limited;
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

    // Mole fractions must remain a valid composition: renormalise to sum to one in
    // every cell, or the species balance seeds an impossible gas mixture.
    if (!state.mole_fraction.empty()) {
      const size_t n_cells = state.mole_fraction.begin()->second.var.size();
      for (size_t c = 0; c < n_cells; ++c) {
        double sum = 0.0;
        for (auto& [component, eq] : state.mole_fraction) sum += eq.var[c];
        if (sum > 1e-30) {
          for (auto& [component, eq] : state.mole_fraction) eq.var[c] /= sum;
        }
      }
    }

    // Restart: the three samples are consumed, so the next triple is collected
    // fresh from the (post-jump) trajectory. This keeps the geometric-mode fit
    // clean and spaces the jumps ~3 cycles apart, letting the bed relax between
    // them.
    history_.clear();

    return true;
  }

} // end adsim namespace

CEREAL_CLASS_VERSION(adsim::AccelerateCyclicSteadyState, 0);
