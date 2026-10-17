/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

#include <cmath>

namespace adsim {

  namespace constants {
    constexpr double kIdealGasConstant = 8.31446261815324;
    constexpr double kGravity = 9.80665;
    constexpr double kPI = 3.14159265358979;

    // Robustness floor for pressure/temperature used as a denominator in the
    // governing-equation residuals. The implicit BDF solver can transiently
    // overshoot to a zero or negative P/T while iterating; dividing by that raw
    // value spawns Inf/NaN that then poisons the whole residual vector. Clamping
    // the *denominator* to this floor keeps the residual finite so the solver
    // can recover (error-test-fail -> shrink step) instead of failing silently.
    //
    // The floor is chosen far below any physical value (real columns operate at
    // ~1e3-1e7 Pa and ~10-1e3 K), so it never perturbs a physical (positive)
    // state: the clamp is a no-op unless the state is already unphysical.
    constexpr double kMinDenomPressure = 1e-6;    // [Pa]
    constexpr double kMinDenomTemperature = 1e-6; // [K]

    // Robustness floor for the wall energy-balance denominator (d_o^2 - d_i^2)
    // [m^2]. For a vanishingly thin wall (d_o -> d_i) this difference tends to
    // zero and the per-wall-volume flux terms would blow up to Inf/NaN. The floor
    // is far below any physical annulus cross-section (a real column wall spans
    // millimetres, giving d_o^2 - d_i^2 ~ 1e-4 m^2 or larger), so the clamp is a
    // no-op for physical geometry and only engages for a degenerate wall.
    constexpr double kMinDenomWallArea = 1e-12; // [m^2]

    // Width of the smoothed upwind switch [m/s]. A hard upwind selector
    // ((u > 0) ? left : right) makes the convective residual non-differentiable
    // at u = 0; where the flow stalls (e.g. the end of a pressurisation step as
    // the bed reaches the boundary pressure) Newton then chatters across the
    // kink and IDA collapses the step size. The face value instead blends the
    // two upwind candidates over |u| ~ this velocity. Equal to the velocity
    // absolute tolerance and far below any physical superficial velocity, so for
    // |u| > 20x this value the blend is exactly the hard upwind choice.
    constexpr double kUpwindSmoothingVelocity = 1e-6; // [m/s]
  } // end constants namespace

  namespace detail {
    // Return |x| clamped up to a small positive floor, preserving the value for
    // any physical (positive, >= floor) input. Used only to guard denominators.
    inline double DenomFloor(double x, double floor) {
      return (x > floor) ? x : floor;
    }

    // Weight of the left (u > 0 upstream) value in a smoothed upwind face value:
    // face = w * left + (1 - w) * right. Exactly 1 or 0 once |u| exceeds
    // 20 * kUpwindSmoothingVelocity (tanh(20) rounds to 1), so it reproduces the
    // hard upwind choice bit-for-bit away from stagnant flow.
    inline double UpwindWeight(double u) {
      const double x = u / constants::kUpwindSmoothingVelocity;
      if (x >= 20.0) return 1.0;
      if (x <= -20.0) return 0.0;
      return 0.5 * (1.0 + std::tanh(x));
    }

    // Weight in [0, 1] with which an inlet boundary's stream properties
    // (composition, temperature) replace the interior values as the flow
    // through it rises above `scale`: 0 at zero flow (a stagnant inlet must not
    // inject the stream's composition by dispersion), exactly 1 once
    // |flow| >= 5 * scale, smooth in between so the residual stays
    // differentiable when the flow through the boundary decays to zero.
    inline double InflowWeight(double flow, double scale) {
      if (!(scale > 0.0)) return flow != 0.0 ? 1.0 : 0.0;
      const double x = flow / scale;
      const double x2 = x * x;
      if (x2 >= 25.0) return 1.0;
      return std::tanh(x2);
    }
  } // end detail namespace

} // end adsim namespace
