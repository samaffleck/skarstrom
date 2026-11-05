/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Grid convergence analysis of a scalar result (a cycle KPI at cyclic steady
// state, say) computed on three or more grids, each finer than the last by the
// same refinement ratio r = h_coarse / h_fine (= N_fine / N_coarse for a
// uniform refinement of every layer). The standard verification analysis
// (Roache, "Verification and Validation in Computational Science and
// Engineering", 1998; ASME V&V 20-2009):
//
//   values f_0 ... f_{n-1}, coarse to fine; differences d_i = f_{i+1} - f_i.
//
//   * Monotone convergence: the nonzero differences keep one sign
//     (otherwise: oscillatory) and their magnitudes shrink from each grid to
//     the next (otherwise: divergent). Either failure means the grids are not
//     in the asymptotic range, and no order is reported.
//   * Observed order from the three finest grids (Richardson):
//         p = ln(d_{n-3} / d_{n-2}) / ln r.
//   * Richardson extrapolation, the estimate of the grid-converged value:
//         f* = f_{n-1} + d_{n-2} / (r^p - 1).
//   * Roache's grid convergence index of the finest value, relative to it,
//     with safety factor Fs (1.25 for a study of three or more grids):
//         GCI = Fs |d_{n-2} / f_{n-1}| / (r^p - 1).
//     GCI |f_{n-1}| is a band around f_{n-1} meant to contain f*.
//   * Cells for a target error: with |f_N - f*| ~ |f_{n-1} - f*| (N_{n-1}/N)^p,
//     the error falls to `target` at
//         N = N_{n-1} (|f_{n-1} - f*| / target)^(1/p).
//
// The observed order of every consecutive triple is also returned
// (`orders`); with four or more grids, an order that settles from one triple
// to the next is further evidence of the asymptotic range.
//
// The analysis says nothing about the noise in the values. A difference
// smaller than the noise of the values it is taken from (the integrator and
// CSS tolerances of the runs) measures noise, not the grid: check the finest
// difference against it before trusting an order (test_grid_convergence.cpp
// does). Values should come from runs that reached cyclic steady state.

#include <limits>
#include <string_view>
#include <vector>

namespace adsim {

  enum class GridConvergenceStatus {
    Monotone,         // in the asymptotic range: order, GCI and extrapolation defined
    GridIndependent,  // finest difference exactly zero: f* = finest value, GCI 0, no order
    Oscillatory,      // a sign change in the differences: no order (GCI with the assumed order)
    Divergent,        // a difference did not shrink: no order, GCI or extrapolation
    TooFewLevels,     // fewer than three grids
    Invalid,          // a non-finite value, or a refinement ratio not above 1 / not constant
  };

  // "monotone", "grid independent", "oscillatory", "divergent",
  // "too few levels", "invalid".
  std::string_view ToString(GridConvergenceStatus status);

  struct GridConvergenceOptions {
    double safety_factor = 1.25;  // Fs of the GCI (Roache: 1.25 for three or more grids)
    // Order used for the GCI of an oscillatory sequence, whose observed order
    // is undefined (1: the formal order of the bed-end faces of every scheme).
    double assumed_order = 1.0;
    // Absolute error, in the units of the values, for cells_for_target
    // (1e-3: 0.1 percentage point of a purity or recovery).
    double target_error = 1e-3;
  };

  struct GridConvergenceResult {
    static constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    GridConvergenceStatus status = GridConvergenceStatus::Invalid;
    double refinement_ratio = kNaN;
    std::vector<double> differences;  // d_i = f_{i+1} - f_i, coarse to fine
    std::vector<double> orders;       // observed order of each consecutive triple (NaN if undefined)
    bool monotone = false;            // no sign change among the nonzero differences
    bool shrinking = false;           // |d_i| < |d_{i-1}| for every i (or both zero)
    double finest = kNaN;             // f_{n-1}
    double order = kNaN;              // observed order of the three finest grids (Monotone only)
    bool order_assumed = false;       // gci uses GridConvergenceOptions::assumed_order (Oscillatory)
    double gci = kNaN;                // relative to |finest|
    double error_band = kNaN;         // gci |finest|: absolute half-width of the band around finest
    double extrapolated = kNaN;       // f* (Monotone, GridIndependent)
    // Cells at which the error falls to options.target_error (Monotone, and
    // only when the finest cell count is known). May be below or above the
    // finest grid; not rounded.
    double cells_for_target = kNaN;

    // True when order, GCI and extrapolation are all defined from the data.
    bool InAsymptoticRange() const { return status == GridConvergenceStatus::Monotone; }
  };

  // Building blocks (r = refinement ratio, f1, f2, f3 coarse to fine).
  // NaN where undefined (differences of opposite sign, zero, non-finite).
  double ObservedOrder(double f1, double f2, double f3, double r);
  double RichardsonExtrapolation(double f2, double f3, double r, double order);
  // Relative to |f3|; NaN when f3 is zero.
  double GridConvergenceIndex(double f2, double f3, double r, double order, double safety_factor = 1.25);
  // Cells at which |f_N - f*| falls to target_error, from the finest grid's error.
  double CellsForTargetError(double finest, double extrapolated, double finest_cells, double order,
                             double target_error);

  // The analysis of values on grids refined by a constant ratio r > 1,
  // ordered coarse to fine. finest_cells (optional) gives cells_for_target.
  GridConvergenceResult AnalyseGridConvergence(const std::vector<double>& values, double refinement_ratio,
                                               const GridConvergenceOptions& options = {},
                                               double finest_cells = GridConvergenceResult::kNaN);

  // The same with the cell count of each grid, coarse to fine: the ratio is
  // that of successive cell counts, which must be the same for every pair
  // within ratio_tolerance (relative; cell counts are rounded when a ratio
  // other than 2 is applied to each layer), otherwise the status is Invalid.
  GridConvergenceResult AnalyseGridConvergence(const std::vector<int>& cells, const std::vector<double>& values,
                                               const GridConvergenceOptions& options = {},
                                               double ratio_tolerance = 0.05);

} // namespace adsim
