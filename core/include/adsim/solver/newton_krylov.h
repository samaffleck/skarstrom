/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Newton-Krylov solution of a fixed point x = G(x), for the cyclic steady
// state (AccelerationMethod::NEWTON_KRYLOV, see cyclic_steady_state.h).
//
// Newton's method on F(x) = G(x) - x: each iteration solves
//     (I - J) d = F(x),   J = dG/dx,
// by GMRES, and moves to x + alpha d. GMRES needs only products J v, which are
// formed by a forward difference of the map,
//     J v ~ (G(x + h v) - G(x)) / h,
// so each Krylov vector costs one evaluation of G (one plain cycle). Nothing
// else about G is needed: it may be any piecewise-smooth map, such as a cycle
// of the integrator with its step switches.
//
// Everything is done in scaled coordinates u = w .* x (w the weights given),
// so that a temperature does not outweigh a mole fraction.
//
//  * Difference increment: h is chosen so that the RMS size of the scaled
//    perturbation h v is fd_relative_step * (1 + RMS(u)). With G evaluated to
//    a relative accuracy of about reltol, fd_relative_step = sqrt(reltol)
//    balances the truncation error of the difference (of order h) against the
//    noise of G divided by h.
//  * Inexact Newton: GMRES stops when its residual is below eta * |F|, with
//    the forcing term eta of Eisenstat and Walker (1996, choice 2):
//    eta_k = 0.9 (|F_k| / |F_{k-1}|)^2, safeguarded by 0.9 eta_{k-1}^2 when
//    that is above 0.1, and kept within [eta_min, eta_max] and above
//    0.5 tolerance / |F_k| (no point solving the linear system more
//    accurately than the stopping test needs).
//  * Damping (Glover, Papathanasiou and Pini): when some entries would go
//    negative, alpha = min over those entries of x_i / (2 (x_i - x_i^new)),
//    i.e. half the step to zero; otherwise alpha = 1.
//  * A step is accepted when it reduces the residual norm. Otherwise the step
//    is halved once; if that does not reduce it either, the solve ends there
//    (the caller falls back to another method).
//
// The residual norm used for the stopping test and the reported history is
// the caller's (`norm`, on the scaled residual); GMRES and the forcing term
// use the Euclidean norm of the scaled vectors.

#include <functional>
#include <limits>
#include <string>
#include <vector>

#include <Eigen/Dense>

namespace adsim {

  struct GmresResult {
    Eigen::VectorXd x;             // the approximate solution (zero start)
    int iterations = 0;            // operator applications (Krylov vectors)
    double residual = 0.0;         // |b - A x| (Euclidean), from the Givens recurrence
    bool converged = false;        // residual <= tolerance
    bool failed = false;           // an operator application failed (x is the best before it)
    // The (iterations x iterations) upper Hessenberg matrix of the Arnoldi
    // process: its eigenvalues (Ritz values) approximate those of A.
    Eigen::MatrixXd hessenberg;
  };

  // GMRES without restart, from x = 0, for A x = b. apply(v, Av) computes A v
  // and returns false on failure. Stops when the residual is <= tolerance
  // (absolute), after max_iterations applications, or on a breakdown (the
  // Krylov space is invariant: x is then exact).
  GmresResult Gmres(const std::function<bool(const Eigen::VectorXd&, Eigen::VectorXd&)>& apply,
                    const Eigen::VectorXd& b, double tolerance, int max_iterations);

  // The dominant Floquet multiplier estimate from the Hessenberg matrix of
  // GMRES on A = I - J: the largest |1 - theta| over its eigenvalues theta.
  // NaN for an empty matrix.
  double DominantMultiplierFromHessenberg(const Eigen::MatrixXd& hessenberg);

  // The damping factor of a step d from x, for entries that must stay
  // non-negative: min over {i : x_i + d_i < 0} of x_i / (2 (-d_i)), else 1.
  // Entries with x_i <= 0 already are ignored (they cannot be halved towards
  // zero; the caller clips them).
  double GloverDamping(const std::vector<double>& x, const std::vector<double>& d);

  struct FixedPointNewtonOptions {
    // Stop when norm(F) <= tolerance.
    double tolerance = 1e-8;
    // Relative size of the difference increment (see the file comment).
    double fd_relative_step = 1e-3;
    int max_newton_iterations = 8;
    int max_krylov = 30;
    double eta_initial = 0.1;
    double eta_max = 0.5;
    double eta_min = 1e-3;
    // Damp steps that would make entries negative (all entries are
    // non-negative quantities).
    bool nonnegative = true;
  };

  struct FixedPointNewtonResult {
    std::vector<double> x;    // the iterate with the smallest residual
    std::vector<double> gx;   // G(x) of that iterate
    double initial_residual = std::numeric_limits<double>::quiet_NaN();
    double residual = std::numeric_limits<double>::quiet_NaN();  // norm(F) at x
    std::vector<double> residual_history;  // norm(F) at each accepted iterate, the start first
    std::vector<int> krylov_vectors;       // per Newton step
    std::vector<double> damping;           // alpha of each accepted step
    int map_evaluations = 0;               // evaluations of G (cycle equivalents)
    int newton_iterations = 0;             // accepted steps
    // From the last GMRES solve (DominantMultiplierFromHessenberg).
    double dominant_multiplier = std::numeric_limits<double>::quiet_NaN();
    bool converged = false;
    std::string outcome;  // why the solve ended
  };

  // Solves x = G(x) from x0 by Newton-Krylov (see the file comment).
  //  * map(x, gx): gx = G(x); false when G cannot be evaluated at x.
  //  * weights: the scaling w (length of x0, positive).
  //  * norm(f): the residual norm on a scaled residual f = w .* (G(x) - x).
  //  * project(x): applied to every new iterate (clipping, renormalising);
  //    may be empty.
  //  * g0, when not null, is G(x0) already evaluated (saves one evaluation).
  FixedPointNewtonResult SolveFixedPointNewtonKrylov(
      const std::function<bool(const std::vector<double>&, std::vector<double>&)>& map,
      const std::vector<double>& x0, const std::vector<double>& weights,
      const std::function<double(const std::vector<double>&)>& norm,
      const std::function<void(std::vector<double>&)>& project,
      const FixedPointNewtonOptions& options, const std::vector<double>* g0 = nullptr);

}  // namespace adsim
