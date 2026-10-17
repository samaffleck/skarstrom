/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/boundary_conditions/boundary_value.h"


namespace adsim{

  double ConstantBoundaryValue::GetValue(double time) {
    return constant_value_.Get();
  }

  double LinearBoundaryValue::GetValue(double time) {
    const double sv = runtime_start_value_.value_or(start_value_.Get());
    const double ev = end_value_.Get();
    const double ramp = runtime_ramp_time_.value_or(ramp_time_.Get());
    // A zero (or negative) ramp is an instantaneous step to the end value.
    if (!(ramp > 0.0)) return ev;
    // Hold the end value once the ramp is complete instead of extrapolating.
    const double tau = std::clamp(time, 0.0, ramp);
    return sv + (ev - sv) * (tau / ramp);
  }

  double ExponentialBoundaryValue::GetValue(double time) {
    const double sv = runtime_start_value_.value_or(start_value_.Get());
    const double ev = end_value_.Get();
    // step_time_ is the time to reach 99% of the change (rate = ln(100) / step_time_).
    // The curve saturates on its own, so no upper clamp is needed (clamping would
    // freeze it at 99% and introduce a jump); a zero time constant is a step.
    const double ts = runtime_step_time_.value_or(step_time_.Get());
    if (!(ts > 0.0)) return ev;
    const double t = std::max(time, 0.0);
    return sv + (ev - sv) * (1.0 - std::exp(-(std::log(100.0) / ts) * t));
  }

} // end adsim namespace

CEREAL_REGISTER_TYPE(adsim::ConstantBoundaryValue)
CEREAL_REGISTER_TYPE(adsim::LinearBoundaryValue)
CEREAL_REGISTER_TYPE(adsim::ExponentialBoundaryValue)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IBoundaryValue, adsim::ConstantBoundaryValue)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IBoundaryValue, adsim::LinearBoundaryValue)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IBoundaryValue, adsim::ExponentialBoundaryValue)
CEREAL_CLASS_VERSION(adsim::ConstantBoundaryValue, 2);
CEREAL_CLASS_VERSION(adsim::LinearBoundaryValue, 3);
CEREAL_CLASS_VERSION(adsim::ExponentialBoundaryValue, 3);
