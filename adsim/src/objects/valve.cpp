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

} // end adsim namespace

CEREAL_CLASS_VERSION(adsim::Valve, 1);
