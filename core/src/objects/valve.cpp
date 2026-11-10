/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/objects/valve.h"
#include "adsim/math/constants.h"


namespace adsim {

  // double Valve::GetMolarFlowRate(double upstream_pressure, double downstream_pressure, 
  //   double molar_mass, double upstream_temperature, double heat_capacity) const {
  //   double Cpg = heat_capacity * molar_mass; // Convert from J/(kg-K) to J/(mol-K)
  //   double Cvg = Cpg - constants::kIdealGasConstant;
  //   double Fk = (Cpg / Cvg) / 1.4;

  //   double x_ = std::max((upstream_pressure - downstream_pressure) / upstream_pressure, 0.0);
  //   double y = 1 - (x_ / (3 * Fk * xt_));

  //   double y_star = std::max(0.66666667, y);
  //   double x_star = std::min(Fk * xt_, x_);

  //   double molar_flow = 8.30834615847399E-06 * cv_ * upstream_pressure * y_star * std::sqrt(x_star / (molar_mass * upstream_temperature));

  //   // return std::max(molar_flow, 1e-8);
  //   return molar_flow;
  // }

  double Valve::GetMolarFlowRate(double upstream_pressure, double downstream_pressure, 
    double molar_mass, double upstream_temperature, double heat_capacity) const 
  {
    double molar_flow = cv_ * (upstream_pressure - downstream_pressure) / 1e5;
    return molar_flow;
  }

  double Valve::OpeningFraction(double time) const {
    return OpeningFraction(time, opening_start_, opening_ramp_time_);
  }

  double Valve::OpeningFraction(double time, double start, double ramp_time) {
    if (!(ramp_time > 0.0) || start >= 1.0) return 1.0;  // no ramp: open at once (defaults)
    if (time >= ramp_time) return 1.0;
    if (time <= 0.0) return start;
    const double x = time / ramp_time;
    return start + (1.0 - start) * x * x * (3.0 - 2.0 * x);
  }

  double Valve::CheckValvePressureDrop(double dp) {
    if (dp <= 0.0) return 0.0;  // closed (a NaN falls through and stays NaN)
    if (dp >= kCheckValveBand) return dp;
    const double x = dp / kCheckValveBand;
    return kCheckValveBand * x * x * (2.0 - x);
  }

} // end adsim namespace
