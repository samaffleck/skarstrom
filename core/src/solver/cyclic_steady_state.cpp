/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/solver/cyclic_steady_state.h"
#include "adsim/solver/cycle.h"
#include "adsim/solver/log.h"
#include "adsim/solver/kpis.h"
#include "adsim/objects/reactor.h"


namespace adsim {

  bool CyclicSteadyState::IsAtCSS(const EquationLog& logger, int cycle_index) {
    return CheckCSS(logger, cycle_index).converged;
  }

  CSSResult CyclicSteadyState::CheckCSS(const EquationLog& logger, int cycle_index) {
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

  CSSResult CyclicSteadyState::CheckProfileResult(const EquationLog& logger, int cycle_index) {
    CSSResult result;
    result.method_name = "Molar Profile";

    // Honour the minimum-cycles floor (shared with the mass-balance methods) so a
    // run is never declared converged before min_cycles_ cycles. cycle_index >= 1
    // is also required because this method differences consecutive cycles.
    if (cycle_index < 1 || cycle_index < min_cycles_ - 1)
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

    if (cycle_index < min_cycles_ - 1) return result;

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
    return result;
  }

  CSSResult CyclicSteadyState::CheckMassBalanceComponentResult(const EquationLog& logger, int cycle_index) {
    CSSResult result;
    result.method_name = "Mass Balance (Component)";
    result.tolerance = mass_balance_tolerance_;

    if (cycle_index < min_cycles_ - 1) return result;

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
    return result;
  }

  bool AccelerateCyclicSteadyState::Accelerate(Equations& state, int cycle_index) {
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
    double scaled_r1_sq = 0.0, scaled_r2_sq = 0.0;
    for (size_t i = 0; i < equations.size(); ++i) {
      if (!extrapolable(equations[i])) continue;
      const auto& x0 = s0[i]; const auto& x1 = s1[i]; const auto& x2 = s2[i];
      if (x0.size() != x1.size() || x1.size() != x2.size() || x2.empty()) continue;

      double mag_sq = 0.0, r1_sq = 0.0, r2_sq = 0.0, r2d = 0.0, dd = 0.0;
      for (size_t c = 0; c < x2.size(); ++c) {
        const double r1 = x1[c] - x0[c];
        const double r2 = x2[c] - x1[c];
        const double d  = r2 - r1;
        mag_sq += x2[c] * x2[c];
        r1_sq  += r1 * r1;
        r2_sq  += r2 * r2;
        r2d    += r2 * d;
        dd     += d * d;
      }
      const double scale_sq = std::max(mag_sq / static_cast<double>(x2.size()), 1e-30);
      const double w = 1.0 / scale_sq;
      num += w * r2d;
      den += w * dd;
      scaled_r1_sq += w * r1_sq;
      scaled_r2_sq += w * r2_sq;
    }

    // No curvature to fit, or the sequence is not contracting toward a fixed point
    // (the scaled step must be shrinking). Either way, do not extrapolate.
    if (den < 1e-300 || scaled_r2_sq >= scaled_r1_sq || scaled_r2_sq < 1e-300) {
      history_.clear();
      return false;
    }

    const double mu = num / den;

    // Safeguard against near-singular curvature: for a single geometric mode of
    // rate rho the coefficient is mu = rho / (rho - 1), so |mu| grows without
    // bound as rho -> 1. Cap |mu| at 5 (rho ~ 0.83); a larger implied jump means
    // the three points are poorly described by one mode, so skip.
    if (!(std::abs(mu) <= 5.0)) {
      history_.clear();
      return false;
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

CEREAL_CLASS_VERSION(adsim::CyclicSteadyState, 1);
CEREAL_CLASS_VERSION(adsim::AccelerateCyclicSteadyState, 0);
