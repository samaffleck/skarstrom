/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

namespace adsim {

  // Per-cell memo of a function of the temperature and up to two parameters.
  //
  // The residual is evaluated at many states that differ only in a few cells'
  // temperatures: a coloured-Jacobian column group perturbs the temperature of
  // about one cell in three, and most groups perturb no temperature at all, so a
  // transcendental that depends on T alone (an isotherm's exp(q / (R T)), the
  // Neufield collision integral) is recomputed with exactly the same arguments
  // most of the time. Get(n, T, k1, k2, f) returns the value stored for cell n
  // when it was computed for the same T and parameter bits, otherwise calls f()
  // and stores the result. Either way the value is the double f() returns, so
  // the results are bit-identical to the uncached code.
  //
  // Keying on the parameters as well as T keeps the memo correct across
  // parameter edits between runs (optimiser, parametric study, isotherm fitter,
  // the GUI) without an explicit reset. Keys are compared as bit patterns, so
  // +0 and -0 are different keys and a NaN only matches the same NaN.
  class TemperatureCache {
  public:
    template <class F>
    double Get(std::size_t n, double T, double k1, double k2, F&& compute) {
      if (n >= slots_.size()) slots_.resize(n + 1);
      Slot& s = slots_[n];
      if (SameBits(s.T, T) && SameBits(s.k1, k1) && SameBits(s.k2, k2)) return s.value;
      const double v = compute();
      s.T = T;
      s.k1 = k1;
      s.k2 = k2;
      s.value = v;
      return v;
    }

    template <class F>
    double Get(std::size_t n, double T, double k1, F&& compute) {
      return Get(n, T, k1, 0.0, static_cast<F&&>(compute));
    }

    void Clear() { slots_.clear(); }
    std::size_t Size() const { return slots_.size(); }

  private:
    struct Slot {
      double T = std::numeric_limits<double>::quiet_NaN();
      double k1 = std::numeric_limits<double>::quiet_NaN();
      double k2 = std::numeric_limits<double>::quiet_NaN();
      double value = std::numeric_limits<double>::quiet_NaN();
    };

    static bool SameBits(double a, double b) {
      std::uint64_t ua, ub;
      std::memcpy(&ua, &a, sizeof ua);
      std::memcpy(&ub, &b, sizeof ub);
      return ua == ub;
    }

    std::vector<Slot> slots_;
  };

} // end adsim namespace
