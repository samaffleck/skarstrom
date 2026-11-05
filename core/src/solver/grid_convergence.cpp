/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "adsim/solver/grid_convergence.h"

#include <cmath>

namespace adsim {

  namespace {

    constexpr double kNaN = GridConvergenceResult::kNaN;

    // log_r(x), with log2 for r = 2 so the common case is exact to the bit.
    double LogBase(double x, double r) { return r == 2.0 ? std::log2(x) : std::log(x) / std::log(r); }

    bool ValidRatio(double r) { return std::isfinite(r) && r > 1.0; }

  } // namespace

  std::string_view ToString(GridConvergenceStatus status) {
    switch (status) {
      case GridConvergenceStatus::Monotone: return "monotone";
      case GridConvergenceStatus::GridIndependent: return "grid independent";
      case GridConvergenceStatus::Oscillatory: return "oscillatory";
      case GridConvergenceStatus::Divergent: return "divergent";
      case GridConvergenceStatus::TooFewLevels: return "too few levels";
      case GridConvergenceStatus::Invalid: return "invalid";
    }
    return "invalid";
  }

  double ObservedOrder(double f1, double f2, double f3, double r) {
    const double d_coarse = f2 - f1, d_fine = f3 - f2;
    if (!ValidRatio(r) || !std::isfinite(d_coarse) || !std::isfinite(d_fine) || d_fine == 0.0) return kNaN;
    const double q = d_coarse / d_fine;
    if (!(q > 0.0)) return kNaN;
    return LogBase(q, r);
  }

  double RichardsonExtrapolation(double f2, double f3, double r, double order) {
    if (!ValidRatio(r) || !(order > 0.0)) return kNaN;
    return f3 + (f3 - f2) / (std::pow(r, order) - 1.0);
  }

  double GridConvergenceIndex(double f2, double f3, double r, double order, double safety_factor) {
    if (!ValidRatio(r) || !(order > 0.0) || f3 == 0.0) return kNaN;
    return safety_factor * std::abs((f3 - f2) / f3) / (std::pow(r, order) - 1.0);
  }

  double CellsForTargetError(double finest, double extrapolated, double finest_cells, double order,
                             double target_error) {
    if (!(order > 0.0) || !(finest_cells > 0.0) || !(target_error > 0.0)) return kNaN;
    const double error = std::abs(finest - extrapolated);
    if (!std::isfinite(error)) return kNaN;
    return finest_cells * std::pow(error / target_error, 1.0 / order);
  }

  GridConvergenceResult AnalyseGridConvergence(const std::vector<double>& values, double refinement_ratio,
                                               const GridConvergenceOptions& options, double finest_cells) {
    GridConvergenceResult r;
    r.refinement_ratio = refinement_ratio;
    const size_t n = values.size();
    if (n > 0) r.finest = values.back();
    for (double v : values) {
      if (!std::isfinite(v)) return r;  // Invalid
    }
    if (!ValidRatio(refinement_ratio)) return r;  // Invalid
    for (size_t i = 0; i + 1 < n; ++i) r.differences.push_back(values[i + 1] - values[i]);
    for (size_t i = 0; i + 2 < n; ++i)
      r.orders.push_back(ObservedOrder(values[i], values[i + 1], values[i + 2], refinement_ratio));
    if (n < 3) {
      r.status = GridConvergenceStatus::TooFewLevels;
      return r;
    }

    const auto& d = r.differences;
    r.monotone = true;
    r.shrinking = true;
    double sign = 0.0;
    for (size_t i = 0; i < d.size(); ++i) {
      if (d[i] != 0.0) {
        const double s = d[i] > 0.0 ? 1.0 : -1.0;
        if (sign != 0.0 && s != sign) r.monotone = false;
        sign = s;
      }
      if (i > 0 && !(std::abs(d[i]) < std::abs(d[i - 1]) || (d[i] == 0.0 && d[i - 1] == 0.0)))
        r.shrinking = false;
    }

    const double f2 = values[n - 2], f3 = values[n - 1];
    const double d_fine = d.back();
    if (!r.monotone) {
      r.status = GridConvergenceStatus::Oscillatory;
      if (d_fine != 0.0) {
        r.order_assumed = true;
        r.gci = GridConvergenceIndex(f2, f3, refinement_ratio, options.assumed_order, options.safety_factor);
      }
      else {
        r.gci = 0.0;
      }
    }
    else if (d_fine == 0.0) {
      // Includes an unchanging sequence. Only exact equality of the two finest
      // values counts: the differences are not compared with any noise level.
      r.status = GridConvergenceStatus::GridIndependent;
      r.gci = 0.0;
      r.extrapolated = f3;
    }
    else if (!r.shrinking) {
      r.status = GridConvergenceStatus::Divergent;
    }
    else {
      r.status = GridConvergenceStatus::Monotone;
      r.order = r.orders.back();
      r.gci = GridConvergenceIndex(f2, f3, refinement_ratio, r.order, options.safety_factor);
      r.extrapolated = RichardsonExtrapolation(f2, f3, refinement_ratio, r.order);
      r.cells_for_target = CellsForTargetError(f3, r.extrapolated, finest_cells, r.order, options.target_error);
    }
    r.error_band = r.gci * std::abs(f3);
    return r;
  }

  GridConvergenceResult AnalyseGridConvergence(const std::vector<int>& cells, const std::vector<double>& values,
                                               const GridConvergenceOptions& options, double ratio_tolerance) {
    if (cells.size() != values.size()) return GridConvergenceResult{};  // Invalid
    if (cells.size() < 2)
      return AnalyseGridConvergence(values, 2.0, options, cells.empty() ? kNaN : double(cells.back()));
    for (int c : cells) {
      if (c <= 0) return GridConvergenceResult{};
    }
    // The ratio of every pair; exact when they are all equal (cell counts
    // doubled), otherwise their geometric mean, each within ratio_tolerance of it.
    std::vector<double> ratios;
    for (size_t i = 0; i + 1 < cells.size(); ++i) ratios.push_back(double(cells[i + 1]) / double(cells[i]));
    double ratio = ratios.front();
    bool all_equal = true;
    for (double q : ratios) all_equal = all_equal && q == ratio;
    if (!all_equal) {
      ratio = std::pow(double(cells.back()) / double(cells.front()), 1.0 / double(ratios.size()));
      for (double q : ratios) {
        if (!(std::abs(q / ratio - 1.0) <= ratio_tolerance)) {
          GridConvergenceResult r;
          r.refinement_ratio = kNaN;
          r.finest = values.back();
          return r;  // Invalid: not a constant ratio
        }
      }
    }
    return AnalyseGridConvergence(values, ratio, options, double(cells.back()));
  }

} // namespace adsim
