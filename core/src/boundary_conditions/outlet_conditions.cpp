/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/boundary_conditions/outlet_conditions.h"


namespace adsim {

  void OutletConditions::Clear() {
    u_.clear();
    T_.clear();
    P_.clear();
    t_.clear();
    yi_.clear();
  }

  void OutletConditions::LogOutletConditions(double u,
    double T, double P, double t,
    const std::unordered_map<std::string, double>& yi) {
    u_.emplace_back(u);
    T_.emplace_back(T);
    P_.emplace_back(P);
    t_.emplace_back(t);

    for (const auto& [component, value] : yi) {
      yi_[component].emplace_back(value);
    }
  }

  double OutletConditions::GetP(double time) const {
    return InterpolateTimeSeries(t_, P_, time);
  }

  double OutletConditions::GetT(double time) const {
    return InterpolateTimeSeries(t_, T_, time);
  }

  double OutletConditions::GetU(double time) const {
    return InterpolateTimeSeries(t_, u_, time);
  }

  double OutletConditions::GetYi(double time, const std::string& component) const {
    return InterpolateTimeSeries(t_, yi_.at(component), time);
  }

  double InterpolateTimeSeries(const std::vector<double>& times, const std::vector<double>& values, double time) {
    if (times.size() != values.size()) {
      throw std::runtime_error("Mismatched time and value vector sizes.");
    }

    if (times.empty()) {
      throw std::runtime_error("No data available for interpolation.");
    }

    if (time <= times.front()) return values.front();
    if (time >= times.back()) return values.back();

    // Find interval [i, i+1] such that times[i] <= time <= times[i+1]
    auto upper = std::upper_bound(times.begin(), times.end(), time);
    size_t i = std::distance(times.begin(), upper) - 1;

    double t0 = times[i];
    double t1 = times[i + 1];
    double v0 = values[i];
    double v1 = values[i + 1];

    double factor = (time - t0) / (t1 - t0);
    return v0 + factor * (v1 - v0);
  }

} // end adsim namespace
