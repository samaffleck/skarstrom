/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/boundary_conditions/flow_boundary_value.h"
#include "adsim/math/constants.h"
#include "adsim/equations/equation.h"


namespace adsim {

  double GetArea(const Equations& state, BoundaryLocation loc) {
    auto& A = state.domain_vectors.area;
    double area = A[0];
    if (loc == BoundaryLocation::TOP) area = A[A.size() - 1];
    return area;
  }

  double GetMolarMassFace(const Equations& state, BoundaryLocation loc) {
    auto& molar_mass = state.molar_mass.var;
    double mr = molar_mass[0];
    if (loc == BoundaryLocation::TOP) mr = molar_mass[molar_mass.size() - 1];
    return mr;
  }

  double GetMolarFace(const Equations& state, BoundaryLocation loc) {
    auto& T = state.fluid_temperature.var;
    auto& P = state.pressure.var;
    auto& C = state.concentration.var;
    auto& u = state.velocity.var;
    auto N = (int)P.size() - 1;

    int i = 0; // Ghost cell index
    int ui = 0; // Ghost cell index
    int d = 1; // Interior direction
    if (loc == BoundaryLocation::TOP) {
      i = N; d = -1; ui = N - 1;
    }

    double T_boundary = 0.5 * (T[i] + T[i + d]);
    double P_boundary = 0.5 * (P[i] + P[i + d]);
    double C_boundary = P_boundary / (constants::kIdealGasConstant * T_boundary);
    return C_boundary;
  }

  double VelocityBoundaryValue::GetVelocity(const Equations& state, double time, BoundaryLocation loc) {
    return profile_->GetValue(time);
  }

  double VelocityBoundaryValue::GetVolumeFlow(const Equations& state, double time, BoundaryLocation loc) {
    double area = GetArea(state, loc);
    return profile_->GetValue(time) * area;
  }

  double VelocityBoundaryValue::GetMolarFlow(const Equations& state, double time, BoundaryLocation loc) {
    double area = GetArea(state, loc);
    double Cf = GetMolarFace(state, loc);
    return profile_->GetValue(time) * area * Cf;
  }

  double VolumetricFlowBoundaryValue::GetVelocity(const Equations& state, double time, BoundaryLocation loc) {
    double area = GetArea(state, loc);
    return profile_->GetValue(time) / area;
  }

  double VolumetricFlowBoundaryValue::GetVolumeFlow(const Equations& state, double time, BoundaryLocation loc) {
    return profile_->GetValue(time);
  }

  double VolumetricFlowBoundaryValue::GetMolarFlow(const Equations& state, double time, BoundaryLocation loc) {
    double Cf = GetMolarFace(state, loc);
    return profile_->GetValue(time) * Cf;
  }

  double MolarFlowBoundaryValue::GetVelocity(const Equations& state, double time, BoundaryLocation loc) {
    double Cf = GetMolarFace(state, loc);
    double area = GetArea(state, loc);
    return profile_->GetValue(time) / (Cf * area);
  }

  double MolarFlowBoundaryValue::GetVolumeFlow(const Equations& state, double time, BoundaryLocation loc) {
    double Cf = GetMolarFace(state, loc);
    return profile_->GetValue(time) / Cf;
  }

  double MolarFlowBoundaryValue::GetMolarFlow(const Equations& state, double time, BoundaryLocation loc) {
    return profile_->GetValue(time);
  }

  double MassFlowBoundaryValue::GetVelocity(const Equations& state, double time, BoundaryLocation loc) {
    double Mrf = GetMolarMassFace(state, loc); // [kg/mol]
    double molar_flow = profile_->GetValue(time) / Mrf;
    double Cf = GetMolarFace(state, loc);
    double area = GetArea(state, loc);
    return molar_flow / (Cf * area);
  }

  double MassFlowBoundaryValue::GetVolumeFlow(const Equations& state, double time, BoundaryLocation loc) {
    double Mrf = GetMolarMassFace(state, loc); // [kg/mol]
    double molar_flow = profile_->GetValue(time) / Mrf;
    double Cf = GetMolarFace(state, loc);
    return molar_flow / Cf;
  }

  double MassFlowBoundaryValue::GetMolarFlow(const Equations& state, double time, BoundaryLocation loc) {
    double Mrf = GetMolarMassFace(state, loc); // [kg/mol]
    return profile_->GetValue(time) / Mrf;
  }

} // end adsim namespace


CEREAL_REGISTER_TYPE(adsim::VelocityBoundaryValue)
CEREAL_REGISTER_TYPE(adsim::VolumetricFlowBoundaryValue)
CEREAL_REGISTER_TYPE(adsim::MolarFlowBoundaryValue)
CEREAL_REGISTER_TYPE(adsim::MassFlowBoundaryValue)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IFlowBoundaryValue, adsim::VelocityBoundaryValue)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IFlowBoundaryValue, adsim::VolumetricFlowBoundaryValue)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IFlowBoundaryValue, adsim::MolarFlowBoundaryValue)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IFlowBoundaryValue, adsim::MassFlowBoundaryValue)
CEREAL_CLASS_VERSION(adsim::VelocityBoundaryValue, 2);
CEREAL_CLASS_VERSION(adsim::VolumetricFlowBoundaryValue, 2);
CEREAL_CLASS_VERSION(adsim::MolarFlowBoundaryValue, 2);
CEREAL_CLASS_VERSION(adsim::MassFlowBoundaryValue, 2);
