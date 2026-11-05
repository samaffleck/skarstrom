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
    // After the ramp time step_time_ the change left is 1/F of the initial one
    // (rate = ln(F) / step_time_; F = 100 covers 99 %). The curve saturates on
    // its own, so no upper clamp is needed (clamping would freeze it short of
    // the end and introduce a jump); a zero time constant is a step. F <= 1 is
    // rejected by CheckProfileFactor; here it holds the start value.
    const double ts = runtime_step_time_.value_or(step_time_.Get());
    if (!(ts > 0.0)) return ev;
    const double t = std::max(time, 0.0);
    const double f = factor_.Get();
    const double rate = (f > 1.0) ? std::log(f) / ts : 0.0;
    return sv + (ev - sv) * (1.0 - std::exp(-rate * t));
  }

  double DoubleExponentialBoundaryValue::GetValue(double time) {
    const double sv = runtime_start_value_.value_or(start_value_.Get());
    const double ev = end_value_.Get();
    const double t = std::max(time, 0.0);
    const double w = std::min(std::max(weight_.Get(), 0.0), 1.0);
    const double r1 = std::max(rate_1_.Get(), 0.0);
    const double r2 = std::max(rate_2_.Get(), 0.0);
    return ev + (sv - ev) * (w * std::exp(-r1 * t) + (1.0 - w) * std::exp(-r2 * t));
  }

  double LogLinearBoundaryValue::GetValue(double time) {
    const double sv = runtime_start_value_.value_or(start_value_.Get());
    const double ev = end_value_.Get();
    const double ramp = runtime_ramp_time_.value_or(ramp_time_.Get());
    if (!(ramp > 0.0)) return ev;
    const double f = std::clamp(time, 0.0, ramp) / ramp;
    // The ratio ev / sv needs both values positive; otherwise ramp linearly.
    if (!(sv > 0.0 && ev > 0.0)) return sv + (ev - sv) * f;
    // Pin the ends exactly (pow would leave rounding at f = 1).
    if (f >= 1.0) return ev;
    return sv * std::pow(ev / sv, f);
  }

} // end adsim namespace

CEREAL_REGISTER_TYPE(adsim::ConstantBoundaryValue)
CEREAL_REGISTER_TYPE(adsim::LinearBoundaryValue)
CEREAL_REGISTER_TYPE(adsim::ExponentialBoundaryValue)
CEREAL_REGISTER_TYPE(adsim::LogLinearBoundaryValue)
CEREAL_REGISTER_TYPE(adsim::DoubleExponentialBoundaryValue)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IBoundaryValue, adsim::ConstantBoundaryValue)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IBoundaryValue, adsim::LinearBoundaryValue)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IBoundaryValue, adsim::ExponentialBoundaryValue)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IBoundaryValue, adsim::LogLinearBoundaryValue)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IBoundaryValue, adsim::DoubleExponentialBoundaryValue)
CEREAL_CLASS_VERSION(adsim::ConstantBoundaryValue, 2);
CEREAL_CLASS_VERSION(adsim::LinearBoundaryValue, 4);
CEREAL_CLASS_VERSION(adsim::ExponentialBoundaryValue, 4);
CEREAL_CLASS_VERSION(adsim::LogLinearBoundaryValue, 2);
CEREAL_CLASS_VERSION(adsim::DoubleExponentialBoundaryValue, 1);
