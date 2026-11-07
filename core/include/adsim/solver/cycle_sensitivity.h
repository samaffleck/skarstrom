/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Prototype (internal, not a stable API): the cycle map G and its Jacobian
// applied to a few directions, from IDAS forward sensitivities.
//
// x is the differential part of the DAE state at the start of a cycle (after
// the consistent-IC solve), G(x) the differential state at the end of one
// cycle. Integrating the cycle with IDAS forward sensitivities whose initial
// values are the directions v_k gives the Jacobian-vector products
// Phi_x v_k = dG/dx v_k along with G(x), for the price of the sensitivity
// corrector solves (the step sequence and the Jacobian factorisations are the
// base run's). This is what a Newton-Krylov cyclic-steady-state solver needs
// per Krylov vector.
//
// Scope: a single bed without connected streams or buffer tanks (whose
// histories carry state between steps that the sensitivities would also have
// to follow), no logging, no CSS check, no acceleration. The step switches
// follow the solver: at every switch IDA restarts from the end-of-step point
// with IDACalcIC, and here the sensitivities restart with it (IDASensReInit
// at the end-of-step sensitivities, then the same IDACalcIC makes the
// algebraic sensitivities consistent). Boundary values a step derives from
// the bed state when it starts (a pressure ramp's start pressure, say) are
// held fixed by the sensitivities, so they are not differentiated.
//
// See docs/cycle-sensitivity.md for the measured cost and accuracy.

#include <string>
#include <vector>

namespace adsim {

  class UserData;

  struct CycleSensitivityOptions {
    enum class Corrector { Staggered, Simultaneous };
    enum class DifferenceQuotient { Forward, Centred };

    // false: integrate the plain cycle only (the baseline, and the runs a
    // finite-difference check uses); directions are then ignored.
    bool sensitivities = true;
    Corrector corrector = Corrector::Staggered;
    // Difference quotient of the engine's residual for the sensitivity
    // residual F_y s + F_y' s': one extra residual per direction and
    // corrector iteration (forward) or two (centred).
    DifferenceQuotient dq = DifferenceQuotient::Centred;
    // Use IDAS's internal quotient instead (with dummy parameters). Its
    // increment, about reltol * |y|, is too small for this residual: the
    // first consistent-IC solve fails (kept to show it).
    bool idas_internal_dq = false;
    // Added to x before the cycle (differential components, length
    // n_differential); empty for none. For finite differences of G.
    std::vector<double> shift;
    // Override settings.reltol for this run (<= 0: the deck's own).
    double reltol = 0.0;
    // Plain cycles run before the measured one (not counted in its cost), so
    // that it starts at a cycle boundary of a running cycle. From the deck's
    // bed at rest the velocities sit exactly on the upwind switch at the
    // first step's start, where the residual is not differentiable and the
    // sensitivity corrector stalls.
    int warmup_cycles = 0;
    // Relative tolerance of the sensitivity corrector (<= 0: the run's
    // reltol); the absolute tolerances scale with it.
    double sens_reltol = 0.0;
  };

  struct CycleSensitivityResult {
    bool ok = false;
    std::string message;
    // Banded indices of the differential components (the coordinates of x,
    // G and the directions) and their descriptions ("<equation> @ cell i").
    std::vector<size_t> differential_index;
    std::vector<std::string> description;
    std::vector<double> abstol;          // the solver's absolute tolerances
    std::vector<double> x;               // start of the cycle (after the shift and IDACalcIC)
    std::vector<double> g;               // end of the cycle
    std::vector<std::vector<double>> jv; // Phi_x v_k, one per direction

    // Cost of the measured cycle (its consistent-IC solves included).
    double wall_s = 0.0;
    long residual_calls = 0;      // every residual evaluation (IDA, Jacobian, sensitivity DQ)
    long steps = 0;               // accepted internal steps
    long nonlin_iters = 0;
    long jac_evals = 0;
    long err_test_fails = 0;
    long conv_fails = 0;
    long sens_res_evals = 0;      // residual evaluations of the sensitivity DQ
    long sens_nonlin_iters = 0;   // sensitivity corrector iterations
    long sens_conv_fails = 0;
    int reinits = 0;              // step-switch restarts carried through
  };

  // One cycle of the deck in `data` from its initial conditions (plus
  // options.shift), with forward sensitivities along `directions` (each of
  // length n_differential; may be empty). `data` is re-initialised; the deck
  // itself is not changed.
  CycleSensitivityResult CycleJacobianVector(UserData& data,
                                             const std::vector<std::vector<double>>& directions,
                                             const CycleSensitivityOptions& options = {});

} // namespace adsim
