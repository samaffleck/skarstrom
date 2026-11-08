/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/solver/newton_krylov.h"

#include <algorithm>
#include <cmath>
#include <complex>

#include <Eigen/Eigenvalues>

namespace adsim {

  GmresResult Gmres(const std::function<bool(const Eigen::VectorXd&, Eigen::VectorXd&)>& apply,
                    const Eigen::VectorXd& b, double tolerance, int max_iterations) {
    GmresResult result;
    const Eigen::Index n = b.size();
    result.x = Eigen::VectorXd::Zero(n);
    const double beta = b.norm();
    result.residual = beta;
    if (n == 0 || !(beta > tolerance) || !(beta > 0.0)) {
      result.converged = !(beta > tolerance);
      return result;
    }
    const int m = std::max(1, max_iterations);

    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n, m + 1);  // Krylov basis
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(m + 1, m);  // Arnoldi (kept for the Ritz values)
    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(m + 1, m);  // H after the Givens rotations
    Eigen::VectorXd c = Eigen::VectorXd::Zero(m), s = Eigen::VectorXd::Zero(m);
    Eigen::VectorXd g = Eigen::VectorXd::Zero(m + 1);
    g[0] = beta;
    V.col(0) = b / beta;

    int k = 0;  // completed columns
    Eigen::VectorXd w;
    for (int j = 0; j < m; ++j) {
      if (!apply(V.col(j), w) || w.size() != n || !w.allFinite()) {
        result.failed = true;
        break;
      }
      // Modified Gram-Schmidt, twice: each application is a whole cycle, so
      // the orthogonalisation costs nothing by comparison and the second pass
      // keeps the basis orthogonal when the cycle map's columns are close.
      for (int pass = 0; pass < 2; ++pass) {
        for (int i = 0; i <= j; ++i) {
          const double hij = w.dot(V.col(i));
          H(i, j) += hij;
          w -= hij * V.col(i);
        }
      }
      H(j + 1, j) = w.norm();
      k = j + 1;

      // QR of H by Givens rotations: the earlier rotations on the new column,
      // then a new rotation to zero its subdiagonal entry.
      for (int i = 0; i <= j + 1; ++i) R(i, j) = H(i, j);
      for (int i = 0; i < j; ++i) {
        const double t = c[i] * R(i, j) + s[i] * R(i + 1, j);
        R(i + 1, j) = -s[i] * R(i, j) + c[i] * R(i + 1, j);
        R(i, j) = t;
      }
      const double r = std::hypot(R(j, j), R(j + 1, j));
      if (r > 0.0) {
        c[j] = R(j, j) / r;
        s[j] = R(j + 1, j) / r;
      } else {
        c[j] = 1.0;
        s[j] = 0.0;
      }
      R(j, j) = r;
      R(j + 1, j) = 0.0;
      g[j + 1] = -s[j] * g[j];
      g[j] = c[j] * g[j];
      result.residual = std::abs(g[j + 1]);

      // Breakdown: the Krylov space is invariant, the solution is exact.
      const bool breakdown = !(H(j + 1, j) > 1e-14 * beta);
      if (result.residual <= tolerance || breakdown) {
        result.converged = true;
        break;
      }
      V.col(j + 1) = w / H(j + 1, j);
    }
    result.iterations = k;
    if (k > 0) {
      // Back substitution on the triangular R (a zero pivot means a singular
      // operator on the Krylov space: leave that coordinate at zero).
      Eigen::VectorXd y = Eigen::VectorXd::Zero(k);
      for (int i = k - 1; i >= 0; --i) {
        double sum = g[i];
        for (int l = i + 1; l < k; ++l) sum -= R(i, l) * y[l];
        y[i] = (std::abs(R(i, i)) > 0.0) ? sum / R(i, i) : 0.0;
      }
      result.x = V.leftCols(k) * y;
      result.hessenberg = H.topLeftCorner(k, k);
    }
    // A failed application leaves the residual of the last completed column.
    return result;
  }

  double DominantMultiplierFromHessenberg(const Eigen::MatrixXd& hessenberg) {
    if (hessenberg.rows() == 0 || hessenberg.rows() != hessenberg.cols())
      return std::numeric_limits<double>::quiet_NaN();
    Eigen::EigenSolver<Eigen::MatrixXd> es(hessenberg, false);
    if (es.info() != Eigen::Success) return std::numeric_limits<double>::quiet_NaN();
    double dominant = 0.0;
    for (Eigen::Index i = 0; i < es.eigenvalues().size(); ++i)
      dominant = std::max(dominant, std::abs(std::complex<double>(1.0, 0.0) - es.eigenvalues()[i]));
    return dominant;
  }

  double GloverDamping(const std::vector<double>& x, const std::vector<double>& d) {
    double alpha = 1.0;
    const size_t n = std::min(x.size(), d.size());
    for (size_t i = 0; i < n; ++i) {
      if (x[i] > 0.0 && x[i] + d[i] < 0.0) alpha = std::min(alpha, 0.5 * x[i] / (-d[i]));
    }
    return alpha;
  }

  FixedPointNewtonResult SolveFixedPointNewtonKrylov(
      const std::function<bool(const std::vector<double>&, std::vector<double>&)>& map,
      const std::vector<double>& x0, const std::vector<double>& weights,
      const std::function<double(const std::vector<double>&)>& norm,
      const std::function<void(std::vector<double>&)>& project,
      const FixedPointNewtonOptions& options, const std::vector<double>* g0) {
    FixedPointNewtonResult result;
    const size_t n = x0.size();
    result.x = x0;
    if (n == 0 || weights.size() != n) {
      result.outcome = "no unknowns";
      return result;
    }

    std::vector<double> x = x0, gx;
    if (g0 && g0->size() == n) {
      gx = *g0;
    } else {
      ++result.map_evaluations;
      if (!map(x, gx) || gx.size() != n) {
        result.outcome = "the cycle from the start state failed";
        return result;
      }
    }

    auto scaled_residual = [&](const std::vector<double>& xv, const std::vector<double>& gv) {
      Eigen::VectorXd f(static_cast<Eigen::Index>(n));
      for (size_t i = 0; i < n; ++i) f[static_cast<Eigen::Index>(i)] = weights[i] * (gv[i] - xv[i]);
      return f;
    };
    auto as_vector = [](const Eigen::VectorXd& v) { return std::vector<double>(v.data(), v.data() + v.size()); };

    Eigen::VectorXd f = scaled_residual(x, gx);
    double r = norm(as_vector(f));
    result.initial_residual = r;
    result.residual = r;
    result.residual_history.push_back(r);
    result.gx = gx;

    double eta_prev = options.eta_initial;
    double f_euclid_prev = 0.0;
    const double sqrt_n = std::sqrt(static_cast<double>(n));
    bool stopped = false;

    for (int k = 0; k < options.max_newton_iterations; ++k) {
      if (r <= options.tolerance) break;
      const double f_euclid = f.norm();

      // Forcing term (Eisenstat-Walker, choice 2, with its safeguard).
      double eta = options.eta_initial;
      if (k > 0 && f_euclid_prev > 0.0) {
        const double ratio = f_euclid / f_euclid_prev;
        eta = 0.9 * ratio * ratio;
        const double safeguard = 0.9 * eta_prev * eta_prev;
        if (safeguard > 0.1) eta = std::max(eta, safeguard);
      }
      eta = std::clamp(eta, options.eta_min, options.eta_max);
      if (r > 0.0) eta = std::max(eta, std::min(options.eta_max, 0.5 * options.tolerance / r));

      // Difference increment: the scaled perturbation has RMS
      // fd_relative_step * (1 + RMS(u)); GMRES vectors have unit norm.
      double u_sq = 0.0;
      for (size_t i = 0; i < n; ++i) u_sq += (weights[i] * x[i]) * (weights[i] * x[i]);
      const double h = options.fd_relative_step * (1.0 + std::sqrt(u_sq) / sqrt_n) * sqrt_n;

      std::vector<double> xp(n), gp;
      auto apply = [&](const Eigen::VectorXd& v, Eigen::VectorXd& av) {
        for (size_t i = 0; i < n; ++i) xp[i] = x[i] + h * v[static_cast<Eigen::Index>(i)] / weights[i];
        ++result.map_evaluations;
        if (!map(xp, gp) || gp.size() != n) return false;
        av.resize(static_cast<Eigen::Index>(n));
        for (size_t i = 0; i < n; ++i) {
          const Eigen::Index e = static_cast<Eigen::Index>(i);
          av[e] = v[e] - weights[i] * (gp[i] - gx[i]) / h;
        }
        return true;
      };
      const GmresResult gm = Gmres(apply, f, eta * f_euclid, options.max_krylov);
      result.krylov_vectors.push_back(gm.iterations);
      if (gm.hessenberg.size() > 0) result.dominant_multiplier = DominantMultiplierFromHessenberg(gm.hessenberg);
      if (gm.iterations == 0) {
        result.outcome = "a Krylov cycle failed";
        stopped = true;
        break;
      }

      std::vector<double> d(n);
      for (size_t i = 0; i < n; ++i) d[i] = gm.x[static_cast<Eigen::Index>(i)] / weights[i];
      double alpha = options.nonnegative ? GloverDamping(x, d) : 1.0;

      bool accepted = false;
      std::vector<double> xn(n), gn;
      Eigen::VectorXd fn;
      double rn = 0.0;
      for (int trial = 0; trial < 2 && !accepted; ++trial, alpha *= 0.5) {
        for (size_t i = 0; i < n; ++i) xn[i] = x[i] + alpha * d[i];
        if (project) project(xn);
        ++result.map_evaluations;
        if (!map(xn, gn) || gn.size() != n) continue;
        fn = scaled_residual(xn, gn);
        rn = norm(as_vector(fn));
        if (rn < r) {
          accepted = true;
          result.damping.push_back(alpha);
        }
      }
      if (!accepted) {
        result.outcome = "the Newton step did not reduce the residual";
        stopped = true;
        break;
      }
      x = xn;
      gx = gn;
      f = fn;
      f_euclid_prev = f_euclid;
      eta_prev = eta;
      r = rn;
      result.residual_history.push_back(r);
      ++result.newton_iterations;
      if (r < result.residual) {
        result.residual = r;
        result.x = x;
        result.gx = gx;
      }
    }

    result.converged = result.residual <= options.tolerance;
    if (result.converged) result.outcome = "converged";
    else if (!stopped) result.outcome = "iteration limit";
    return result;
  }

}  // namespace adsim
