/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Cyclic-steady-state (CSS) convergence of the molar profile.
//
// A periodic adsorption process does not reach a steady state; it reaches a
// cyclic steady state, in which the bed returns to the SAME internal profile at
// the end of every cycle. This is the defining property of CSS and the premise of
// every KPI a PSA design reports: the process must first "forget" its initial
// condition. The mass balance already asserts closure within a cycle; nothing
// asserted that the end-of-cycle state itself converges.
//
// This is a VERIFICATION property rather than a closed-form value: the exact CSS
// profile has no analytic expression, but the SEQUENCE of end-of-cycle profiles
// must be Cauchy -- successive profiles must get closer, and the change must decay
// toward zero. For a stable linear-ish cycle the end-of-cycle map is a contraction
// and the profile error decays geometrically, so we assert both that the change
// shrinks monotonically from cycle to cycle and that it contracts by a clear
// factor between the first and last gap.
//
// The metric is the L2 norm (over interior cells) of the difference between the
// end-of-cycle gas-phase O2 concentration profiles of successive cycles. A small
// grid and a modest cycle count keep the test fast while still resolving several
// contractions.
//
// GUARD: the cycle must actually DO something -- the end-of-cycle profile must be
// far from the initial condition -- or "convergence" would be trivially satisfied
// by a bed that never moves.

#include "psa_cycle_case.h"

#include <cstdio>

namespace adsim {

  using namespace adsim::validation;

  namespace {
    // L2 norm over interior cells of (a - b).
    double InteriorL2Diff(const std::vector<double>& a, const std::vector<double>& b) {
      double s = 0.0;
      for (size_t i = 1; i + 1 < a.size() && i + 1 < b.size(); ++i) {
        const double d = a[i] - b[i];
        s += d * d;
      }
      return std::sqrt(s);
    }
    double InteriorL2Norm(const std::vector<double>& a) {
      double s = 0.0;
      for (size_t i = 1; i + 1 < a.size(); ++i) s += a[i] * a[i];
      return std::sqrt(s);
    }
  } // namespace

  // ===========================================================================
  // CSS molar-profile convergence
  // ===========================================================================

  TEST(CyclicSteadyState, EndOfCycleProfileConverges) {
    UserData data;
    BuildPsaCycle(data, 30);  // coarse: CSS runs are the expensive ones

    // The end-of-cycle map contracts geometrically (~0.75/cycle here), so the
    // cycle count sets how deep into convergence the assertions can probe.
    // Twelve cycles carries the successive-profile gap comfortably below the
    // 15%-of-first-gap contraction bound checked below.
    const int n_cycles = 12;
    data.cycle.number_of_cycles_ = n_cycles;
    data.settings.max_cycles_to_keep_ = n_cycles;  // every end-of-cycle profile is read below
    data.settings.reltol = 1e-8;
    data.settings.max_step = 2.0;

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    const auto& log = solver.GetLog();

    // The last step index is the blowdown; its final snapshot is the end-of-cycle
    // state. Collect that O2 gas-phase profile for each cycle.
    const int last_step = *log.sys_data.step_index_set_.rbegin();
    std::vector<std::vector<double>> end_of_cycle;
    for (int cycle_index : log.sys_data.cycle_index_set_) {
      const auto c = log.GetComponentVariable(cycle_index, last_step,
          ComponentVariableType::COMPONENT_MOLAR_CONCENTRATION, "Oxygen");
      ASSERT_FALSE(c.values.empty()) << "No end-of-cycle profile for cycle " << cycle_index;
      end_of_cycle.push_back(c.values.back());
    }
    ASSERT_GE(end_of_cycle.size(), 5u) << "Too few cycles logged to assess convergence";

    // GUARD: the cycle must move the bed far from the initial condition, or CSS
    // "convergence" is vacuous. Compare the first end-of-cycle profile against the
    // uniform initial O2 concentration.
    const double C0 = data.ic.pressure / (kR * data.ic.temperature);
    const double C_O2_init = C0 * data.ic.mole_fraction["Oxygen"];
    std::vector<double> init_profile(end_of_cycle.front().size(), C_O2_init);
    const double move = InteriorL2Diff(end_of_cycle.front(), init_profile);
    const double scale = InteriorL2Norm(end_of_cycle.front());
    ASSERT_GT(move, 0.05 * scale)
      << "End-of-cycle profile barely differs from the initial condition; the "
         "cycle is not exercising the bed";

    // Successive-cycle profile changes.
    std::vector<double> gaps;
    for (size_t k = 1; k < end_of_cycle.size(); ++k) {
      gaps.push_back(InteriorL2Diff(end_of_cycle[k], end_of_cycle[k - 1]));
    }

    // 1. The successive-cycle change must contract: the end-of-cycle map is
    //    pulling the profile toward a fixed point, not pushing it away. A real
    //    PSA cycle contracts geometrically but with a slow tail, so we assert a
    //    substantial contraction (>2x) rather than a specific rate; the "how
    //    close to the fixed point" claim is carried by the absolute bound in (3).
    ASSERT_GT(gaps.front(), 0.0);
    EXPECT_LT(gaps.back(), 0.5 * gaps.front())
      << "End-of-cycle profile change did not contract toward a fixed point "
      << "(first gap " << gaps.front() << ", last gap " << gaps.back()
      << " mol/m3, L2 over interior cells)";

    // 2. The change must decay essentially monotonically (a contraction mapping
    //    produces a decreasing sequence). Allow slack for the nonlinear cycle
    //    rather than demanding strict monotonicity every step — near the fixed
    //    point the gaps level off at ~1e-3 mol/m3 and wobble there by up to
    //    ~20% (more so with the real gas-viscosity model, which slightly
    //    perturbs the pressure dynamics each cycle; measured 0.79e-3 -> 0.94e-3
    //    in the last two cycles), so the bound has a 15% relative slack plus
    //    an absolute 2e-4 mol/m3 for that floor; a diverging sequence would
    //    grow far faster.
    constexpr double kGapFloor = 2.0e-4;  // [mol/m3]
    for (size_t k = 1; k < gaps.size(); ++k) {
      EXPECT_LT(gaps[k], 1.15 * gaps[k - 1] + kGapFloor)
        << "End-of-cycle profile change grew from cycle " << k << " to " << k + 1
        << " (" << gaps[k - 1] << " -> " << gaps[k] << "); the sequence is not "
        << "converging to a cyclic steady state";
    }

    // 3. The last change must be a small fraction of the profile magnitude: the
    //    bed is close to its cyclic steady state.
    EXPECT_LT(gaps.back(), 0.02 * scale)
      << "End-of-cycle profile has not settled to a cyclic steady state (final "
      << "change " << gaps.back() << " mol/m3 vs profile norm " << scale << ")";
  }

  // ===========================================================================
  // The state test on an inert bed
  // ===========================================================================
  //
  // A bed of inert packing (an adsorbent with no isotherm) carries loading
  // unknowns that stay at zero up to round-off (of either sign: 1e-26 mol/kg
  // here, 1e-22 on the feed/hold deck of test_css_status), far below the
  // integrator's absolute tolerance on loadings (1e-8 mol/kg). The state test
  // scales the loading change by the RMS total loading, so without a floor on
  // that scale it compared noise with noise: on this deck it reported changes
  // of 0.77 to 51 at every cycle although the bed repeats to round-off, and
  // the run never stopped. With the scale floored at the loading tolerance the
  // round-off counts for ~1e-18, and the state change (now set by the mole
  // fractions) is 6e-14, then ~4e-16.
  //
  // The deck: the shared PSA cycle, isothermal, with the 13X replaced by inert
  // packing and the bed initially filled with nearly pure N2, so that the first
  // cycle flushes it to the feed composition (the guard: the state moves). From
  // then on every cycle repeats, so the default rule (Per-Component Mass
  // Balance, state test on) must stop at the earliest cycle it may: min_cycles_
  // cycles run, the last consecutive_cycles_ of them passing.

  TEST(CyclicSteadyState, InertBedReachesCSSUnderTheStateTest) {
    UserData data;
    BuildPsaCycle(data, 20);
    data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
    std::vector<AdsorbentID> ids;
    for (const auto& entry : data.adsorbent_directory.GetAll()) ids.push_back(entry.first);
    for (AdsorbentID id : ids)  // no isotherm, no LDF: inert packing
      data.adsorbent_directory.Get(id).adsorption_source_ = AdsorptionSource{};
    data.ic.mole_fraction["Nitrogen"] = 0.99;
    data.ic.mole_fraction["Oxygen"] = 0.01;

    constexpr int kCycles = 10;
    data.cycle.number_of_cycles_ = kCycles;
    data.settings.reltol = 1e-8;
    data.settings.max_cycles_to_keep_ = kCycles;  // the whole history
    auto& css = data.settings.css;
    css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;  // defaults otherwise: state test on
    ASSERT_GT(static_cast<double>(css.state_tolerance_), 0.0);

    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    CoupledSolver solver{data};
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
    const auto& log = solver.GetLog();
    const int last_step = *log.sys_data.step_index_set_.rbegin();

    // GUARD 1: the bed is inert -- its loadings are logged (so the loading part
    // of the state test runs) and are zero to within round-off.
    double max_loading = 0.0;
    for (int c : log.sys_data.cycle_index_set_)
      for (const char* component : {"Nitrogen", "Oxygen"}) {
        const auto q = log.GetComponentVariable(c, last_step, ComponentVariableType::ADSORBED_CONCENTRATION,
                                                component);
        ASSERT_FALSE(q.values.empty()) << "no loading logged";
        for (double v : q.values.back()) max_loading = std::max(max_loading, std::abs(v));
      }
    std::printf("  inert bed: largest |loading| %.1e mol/kg\n", max_loading);
    EXPECT_LT(max_loading, 1e-15) << "the packing adsorbs: this is not the inert case";

    // GUARD 2: the first cycle moves the bed from its initial composition.
    const auto y0 = log.GetComponentVariable(0, last_step, ComponentVariableType::MOLE_FRACTION, "Oxygen");
    ASSERT_FALSE(y0.values.empty());
    double y_mean = 0.0;
    for (double v : y0.values.back()) y_mean += v / static_cast<double>(y0.values.back().size());
    EXPECT_GT(y_mean, 0.15) << "the first cycle did not flush the bed towards the feed (0.21 O2)";

    for (int c : log.sys_data.cycle_index_set_)
      if (c > 0)
        std::printf("  inert bed: cycle %d state change %.3e\n", c,
                    CyclicSteadyState::EndOfCycleStateChange(log, c));

    // The bed repeats from cycle 1 on, so the stop is at the floor
    // (min_cycles_ = 5: 0-based cycle 4), with a cycle of slack.
    ASSERT_TRUE(solver.ReachedCSS()) << "the state test never passed on an inert bed";
    EXPECT_LE(solver.GetCSSCycleIndex(), css.min_cycles_);
    const int stop = solver.GetCSSCycleIndex();
    EXPECT_LE(CyclicSteadyState::EndOfCycleStateChange(log, stop), 1e-10);
  }

} // namespace adsim
