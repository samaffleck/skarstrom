/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/solver/cycle.h"
#include "adsim/math/constants.h"
#include "adsim/objects/reactor.h"
#include "adsim/objects/valve.h"


namespace adsim {

  // Boundary-face momentum residual of an Ergun layer under a boundary
  // condition that fixes the face's pressure drop g (a pressure condition:
  // g = 2 d (P_bc - P_first) from the ghost-cell relation P_ghost + P_first =
  // 2 P_bc; a closed end: g = 0). It replaces UpdateVelocity's implicit
  // residual u - (g - B|u|u)/A on that face by the explicit root
  //   u - ErgunVelocity(g, A, B),
  // written with g from the boundary value rather than the ghost pressure.
  // Both have the same solution, but at a step switch, when the condition
  // changes and the face velocity must jump, they behave very differently in
  // IDA's corrector, a modified Newton iteration whose Jacobian is evaluated
  // at the predicted (pre-switch) state:
  //  * the implicit residual has slope (A + 2 B |u|) / A in u, which changes
  //    by an order of magnitude between the old and the new velocity of a
  //    high-Reynolds face (~33 -> 1, i.e. 1 + 2 x 16, when the 15 bar biogas feed face of
  //    BiogasPSA_8Step_Skarstrom_Ergun is closed),
  //    so the iteration contracts too slowly (or overshoots) and IDA declares
  //    a convergence failure at every step size down to hmin;
  //  * the explicit root written in the ghost pressure has a square-root
  //    characteristic in g (slope 1/A at g = 0, ~1/sqrt(B|g|) beyond), so a
  //    jump of the ghost pressure is mis-linearised the same way.
  // Written in u and the boundary value, the residual is linear in u with
  // slope 1 and depends only on the first-cell pressure, which is continuous
  // across a switch (it follows the differential concentration and
  // temperature): Newton resolves the jump in one iteration. A flow
  // condition instead fixes the face velocity through its molar-flow
  // residual, which is where the implicit form (linear in the pressures) is
  // the well-conditioned one, so flow conditions keep it. Darcy (Kozeny-
  // Carman) faces are left untouched.
  static void ImposeBoundaryFaceVelocity(Equations& state, BoundaryLocation location, double g) {
    auto& u = state.velocity.var;
    auto& res_u = state.velocity.res_var;
    if (u.empty() || res_u.size() != u.size()) return;
    const int face = (location == BoundaryLocation::BOTTOM) ? 0 : static_cast<int>(u.size()) - 1;
    double A = 0.0;
    double B = 0.0;
    if (!ErgunFaceCoefficients(state, nullptr, face, A, B)) return;
    res_u[face] = u[face] - ErgunVelocity(g, A, B);
  }

  // Pressure Boundary Conditions
  PressureBoundaryCondition::PressureBoundaryCondition() {
    pressure_ = std::make_unique<ConstantBoundaryValue>(101325.0);
    pump_ = CreatePump(PumpType::COMPRESSOR);
  }

  void PressureBoundaryCondition::Initialise(Equations& state, double step_time) {
    auto& P = state.pressure.var;

    int i = 0;
    int d = 1;
    int N = P.size() - 1;
    if (location_ == BoundaryLocation::TOP) {
      i = N; d = -1;
    }

    double Pstart = 0.5 * (P[i] + P[i + d]);
    pressure_->InitialiseStartValue(Pstart);
    pressure_->InitialiseStepTime(step_time);
    step_time_cached_ = step_time;
    isothermal_temperature_ = state.IsothermalTemperature();
  }

  void PressureBoundaryCondition::Apply(Equations& state, double time) {
    auto& molar_concentration_map = state.molar_concentration;
    
    auto& temperature = state.fluid_temperature;
    auto& T = temperature.var;
    auto& res_T = temperature.res_var;
    
    auto& pressure = state.pressure;
    auto& P = pressure.var;
    auto& res_P = pressure.res_var;

    auto& mole_frac = state.mole_fraction;

    auto& concentration = state.concentration;
    auto& C = concentration.var;
    auto& res_C = concentration.res_var;

    const double R = constants::kIdealGasConstant;

    int N = (int)T.size() - 1;  // Last valid index (includes ghost cells: indices 0 and N)

    int i = 0;
    int d = 1;
    if (location_ == BoundaryLocation::TOP) {
     i = N;
     d = -1; 
    }

    // Calculate outlet pressure for pump calculations (if needed)
    double outlet_pressure = 0.5 * (P[i] + P[i + d]);

    // Get inlet temperature (for INLET boundaries)
    double inlet_temperature = 0.0;
    if (direction_ == BoundaryDirection::INLET) {
      inlet_temperature = InletTemperature(time); // Same source as Log()

      if (pump_->use_pump_ && !isothermal_temperature_) { // no compression heating in an isothermal run
        // Increases the inlet temperature due to compression
        double outlet_temperature = pump_->CalculateOutletTemperature(outlet_pressure, inlet_temperature);
        inlet_temperature = outlet_temperature; // The outlet temp from the compressor is used as the inlet temp to the bed.
      }
    }

    // Concentration boundary (based on specified pressure)
    double Tf = 0.0;
    if (direction_ == BoundaryDirection::INLET) {
      Tf = inlet_temperature; // Use specified inlet temperature for INLET
    } else {
      Tf = 0.5 * (T[i] + T[i + d]); // Use average temperature for OUTLET
    }
    
    // double Cf = pressure_->GetValue(time) / (R * Tf);
    // res_C[i] = C[i] + C[i + d] - 2 * Cf;

    double Pf = pressure_->GetValue(time);
    res_P[i] = P[i] + P[i + d] - 2 * Pf;

    // Ergun layer: this condition fixes the boundary-face pressure drop, so
    // the face velocity is imposed as the explicit Ergun root of it (see
    // ImposeBoundaryFaceVelocity). Darcy faces keep their residual.
    ImposeBoundaryFaceVelocity(state, location_, 2.0 * d * (Pf - P[i + d]));

    // Temperature boundary condition
    if (direction_ == BoundaryDirection::INLET) {
      res_T[i] = T[i] + T[i + d] - 2 * inlet_temperature;
    } else { // OUTLET
      res_T[i] = T[i] - T[i + d];
    }

    // Component molar concentration boundary conditions
    for (auto& [comp, molar_conc] : molar_concentration_map) {
      auto& Yi = mole_frac.at(comp).var;
      auto& res_Yi = mole_frac.at(comp).res_var;
      auto& Ci = molar_conc.var;
      // auto& res_Ci = molar_conc.res_var;
      
      if (direction_ == BoundaryDirection::INLET) {
        // Inlet mole fraction: buffer tank, connected stream or user value
        // (same source as Log()).
        const double Y_inlet = InletMoleFraction(time, comp);
        
        // double Cf_inlet = 0.5 * (C[i] + C[i + d]);
        // double Ci_inlet = Cf_inlet * Y_inlet;
        // res_Ci[i] = Ci[i] + Ci[i + d] - 2 * Ci_inlet;

        // Enforce inlet composition directly on the ghost cell.
        // This matches the upwind face reconstruction used by the species equation.
        res_Yi[i] = Yi[i] - Y_inlet;
      } else { // OUTLET
        // Zero flux boundary condition
        // res_Ci[i] = Ci[i] - Yi[i + d] * C[i];

        res_Yi[i] = Yi[i] - Yi[i + d];
      }
    }
  }

  bool PressureBoundaryCondition::PerformChecks(std::string& error_msg, const Fluid& fluid) {
    bool connection_checks_passed = connection_.PerformChecks(error_msg);
    if (!connection_checks_passed)
      return false;

    if (!buffer_tank_.PerformChecks(error_msg))
      return false;

    if (pump_ && !pump_->PerformChecks(error_msg)) {
      error_msg += "(pump of the pressure boundary) ";
      return false;
    }

    // Key the inlet composition exactly by the fluid's components (back-fill
    // missing ones with 0, drop removed ones). Without this a map that sums to
    // 1 but lacks a fluid component passes the check below and Apply() would
    // meet a missing key inside the IDA residual.
    ReconcileMoleFractionMap(mole_fraction_, fluid.components_);

    if (direction_ == BoundaryDirection::OUTLET)
      return true;

    double total_mole_fraction = 0;
    for (const auto& [comp, yi] : mole_fraction_) {
      (void)comp;
      total_mole_fraction += yi;
    }

    if (std::abs(1 - total_mole_fraction) > 1e-10) {
      error_msg += "Sum of mole fractions in boundary condition does not equal 1. ";
      return false;
    }

    return true;
  }

  void PressureBoundaryCondition::Log(const Equations& state, double time) {
    auto& molar_flowrate = state.molar_flow_rate.var;
    auto& component_mole_fractions = state.mole_fraction;
    auto& T = state.fluid_temperature.var;
    auto& P = state.pressure.var;
    
    double molar_flow = 0.0;
    double temperature = 0.0;
    double suction_temperature = 0.0; // Pump suction-side (pre-compression) temperature
    double pressure = 0.0;
    double enthalpy_flow = 0.0;
    std::unordered_map<std::string, double> component_molar_flow{};
    std::unordered_map<std::string, double> mole_fraction{};

    if (location_ == BoundaryLocation::BOTTOM) {
      molar_flow = molar_flowrate[0];
      pressure = 0.5 * (P[0] + P[1]); 
    } else {
      int N = T.size() - 1;
      molar_flow = molar_flowrate[molar_flowrate.size() - 1];
      pressure = 0.5 * (P[N] + P[N - 1]);
    }

    int i = 0;
    int d = 1;
    int N = static_cast<int>(T.size()) - 1;
    if (location_ == BoundaryLocation::TOP) {
      i = N;
      d = -1;
    }

    // Log the stream-side values:
    // - Inlet: the specified/connected inlet stream values (the values used by Apply()).
    // - Outlet: boundary-adjacent bed values.
    if (direction_ == BoundaryDirection::INLET) {
      temperature = InletTemperature(time);

      // The source stream is the pump's suction side; keep its temperature for
      // the power calculation and log the post-compression temperature.
      suction_temperature = temperature;
      if (pump_->use_pump_ && !isothermal_temperature_) { // no compression heating in an isothermal run
        temperature = pump_->CalculateOutletTemperature(pressure, temperature);
      }

      for (const auto& [component, yi] : component_mole_fractions) {
        (void)yi;
        double y = InletMoleFraction(time, component);
        if (!std::isfinite(y)) y = 0.0;
        if (y < 0.0 && std::abs(y) < 1e-12) y = 0.0;
        if (y < 0.0) y = 0.0;
        mole_fraction[component] = y;
      }
    } else {
      temperature = 0.5 * (T[i] + T[i + d]);
      suction_temperature = temperature; // The bed is the pump's suction side.
      for (const auto& [component, yi] : component_mole_fractions) {
        mole_fraction[component] = 0.5 * (yi.var[i] + yi.var[i + d]);
      }
    }

    double sum_y = 0.0;
    for (const auto& [component, y] : mole_fraction) {
      (void)component;
      sum_y += y;
    }
    if (sum_y > 0.0) {
      for (auto& [component, y] : mole_fraction) {
        (void)component;
        y /= sum_y;
      }
    }

    for (const auto& [component, y] : mole_fraction) {
      component_molar_flow[component] = molar_flow * y;
    }

    enthalpy_flow = molar_flow * temperature;

    double power = 0.0;
    if (pump_->use_pump_) {
      power = pump_->CalculatePower(pressure, molar_flow, suction_temperature);
    }

    stream_data_.Log(molar_flow, enthalpy_flow, temperature, pressure, power, time, component_molar_flow, mole_fraction);
  }

  void PressureBoundaryCondition::ClearLog() {
    stream_data_.Clear();
  }

  // Flow Boundary Conditions
  void FlowBoundaryCondition::Initialise(Equations& state, double step_time) {
    for (const auto& stream : streams_) {
      stream->isothermal_temperature_ = state.IsothermalTemperature();
      stream->Initialise(state, step_time, location_, direction_);
    }
  }

  void FlowBoundaryCondition::Apply(Equations& state, double time) {
    auto& molar_concentration_map = state.molar_concentration;
    
    auto& temperature = state.fluid_temperature;
    auto& T = temperature.var;
    auto& res_T = temperature.res_var;
    
    auto& pressure = state.pressure;
    auto& P = pressure.var;
    auto& res_P = pressure.res_var;

    auto& mole_frac = state.mole_fraction;

    auto& concentration = state.concentration;
    auto& C = concentration.var;
    auto& res_C = concentration.res_var;

    const auto& velocity = state.velocity;
    const auto& u = velocity.var;

    const auto& viscosity = state.viscosity;
    const auto& vis = viscosity.var;

    const auto& dx = state.domain_vectors.cell_width;
    const auto& K = state.domain_vectors.permeability;
    const auto& eb = state.domain_vectors.eb;
    const auto& A = state.domain_vectors.area;

    int N = (int)T.size() - 1;  // Last valid index (includes ghost cells: indices 0 and N)
    int N_u = (int)u.size() - 1;  // Last valid velocity face index

    int i = 0; // Ghost cell index
    int ui = 0; // Face index at boundary
    int d = 1; // Interior direction
    if (location_ == BoundaryLocation::TOP) {
      i = N;  
      ui = N - 1;  
      d = -1;
    }
    
    // Pressure boundary condition
    // Match interior equation form: F = (T/P) * (P_face * u / T_face)
    // At boundary, use face-averaged values
    double P_face = 0.5 * (P[i] + P[i + d]);
    double T_face = 0.5 * (T[i] + T[i + d]);
    double u_face = u[ui];
    
    // Use first interior cell's properties for material parameters
    int interior_idx = i + d;
    double P_interior = P[interior_idx];
    double T_interior = T[interior_idx];
    
    double Cf = P_face / (R * T_face);
    double molar_flow_face = 0.0;
    const double beta_i = state.domain_vectors.inertial_coefficient[i];
    if (beta_i == 0.0) {
      // Darcy (Kozeny-Carman) layer: legacy closed form, unchanged.
      molar_flow_face = (Cf * K[i] * A[i] / (vis[i] * dx[i])) * d * (P[i] - P[i + d]);
    } else {
      // Ergun layer: the flow through the boundary face is carried by the
      // face's own velocity unknown, whose (implicit) Ergun momentum residual
      // UpdateVelocity sets between the ghost and the first bed cell. Using
      // it, rather than re-solving the quadratic here, keeps one momentum law
      // per face and the residual free of the square-root characteristic.
      molar_flow_face = Cf * A[i] * u_face;
    }
    double total_molar_flow = GetTotalMolarFlow(state, time);
    res_P[i] = molar_flow_face - total_molar_flow;

    // The inlet stream's temperature/composition (flow-weighted over the
    // streams) applies while gas flows in; at zero flow the ghost cell takes the
    // interior value instead (no dispersive injection through a stagnant inlet).
    // Blend smoothly between the two: a hard switch at a tiny flow threshold
    // makes the residual jump whenever the boundary flow decays through zero
    // (e.g. a valve at pressure equalisation) and Newton cannot converge.
    const double inflow_scale = Cf * A[i] * constants::kUpwindSmoothingVelocity;  // [mol/s]
    const double w_in = detail::InflowWeight(total_molar_flow, inflow_scale);
    auto inlet_value = [&](double stream_flow_weighted, double interior) {
      if (w_in == 0.0) return interior;
      const double stream_value = stream_flow_weighted / total_molar_flow;
      if (!std::isfinite(stream_value)) return interior;
      if (w_in == 1.0) return stream_value;
      return w_in * stream_value + (1.0 - w_in) * interior;
    };

    // Temperature boundary condition
    if (direction_ == BoundaryDirection::INLET) {
      double total_enthalpy_flow = GetTotalEnthalpyFlow(state, time);
      double Tin = inlet_value(total_enthalpy_flow, T[interior_idx]);
      res_T[i] = T[i] + T[interior_idx] - 2 * Tin;
    } else { // OUTLET
      // Zero flux boundary condition for outlet (zero gradient)
      res_T[i] = T[i] - T[interior_idx];
    }

    // Mole fraction boundary conditions
    for (auto& [comp, molar_conc] : molar_concentration_map) {
      auto& Yi = mole_frac.at(comp).var;
      auto& res_Yi = mole_frac.at(comp).res_var;
      double total_component_molar_flow = GetTotalComponentMolarFlow(state, comp, time);
      
      // Use first interior cell's properties
      int interior_idx = i + d;
      double P_interior = P[interior_idx];
      double T_interior = T[interior_idx];
      double P_face = 0.5 * (P[i] + P[interior_idx]);
      double T_face = 0.5 * (T[i] + T[interior_idx]);
      
      if (direction_ == BoundaryDirection::INLET) {
        // Calculate inlet mole fraction from specified component and total molar flows
        double Y_inlet = inlet_value(total_component_molar_flow, Yi[interior_idx]);

        // Enforce inlet composition directly on the ghost cell.
        // This matches the upwind face reconstruction used by the species equation.
        res_Yi[i] = Yi[i] - Y_inlet;
      } else { // OUTLET
        // Zero flux boundary condition (zero gradient)
        // This ensures both advective and dispersive fluxes are zero
        res_Yi[i] = Yi[i] - Yi[interior_idx];
      }
    }
  }

  void FlowBoundaryCondition::AddStream(const std::string& stream_name) {
    auto it = std::find_if(streams_.begin(), streams_.end(),
      [&stream_name](const auto& stream) {
      return stream->name_ == stream_name;
    });

    if (it == streams_.end()) { // wasnt found in list 
      auto new_stream = CreateStream(BoundaryStreamType::NONE); // No valve is default
      new_stream->name_ = stream_name;
      streams_.push_back(std::move(new_stream));
    }
  }

  void FlowBoundaryCondition::RemoveStream(const std::string& stream_name) {
    auto it = std::find_if(streams_.begin(), streams_.end(),
      [&stream_name](const auto& stream) {
      return stream->name_ == stream_name;
    });

    if (it != streams_.end()) 
      streams_.erase(it);
  }

  double FlowBoundaryCondition::GetTotalMolarFlow(const Equations& state, double time) const {
    double total = 0.0;

    for (auto& stream : streams_) {
      total += stream->GetMolarFlow(state, time, location_);
    }

    return total;
  }

  double FlowBoundaryCondition::GetTotalComponentMolarFlow(const Equations& state, const std::string& component, double time) const {
    double total = 0.0;

    for (auto& stream : streams_) {
      total += stream->GetComponentMolarFlow(state, time, location_, component);
    }

    return total;
  }

  double FlowBoundaryCondition::GetTotalEnthalpyFlow(const Equations& state, double time) const {
    double total = 0.0;

    for (auto& stream : streams_) {
      total += stream->GetEnthalpyFlow(state, time, location_);
    }

    return total;
  }

  bool FlowBoundaryCondition::PerformChecks(std::string& error_msg, const Fluid& fluid) {
    bool checksSuccess = true;

    if (streams_.empty()) {
      error_msg += "Stream list is empty. Add at least one stream. ";
      checksSuccess = false;
    }

    for (const auto& stream : streams_) {
      stream->SetDirection(direction_);
      if (stream->PerformChecks(error_msg, fluid) == false) {
        checksSuccess = false;
      }
    }

    return checksSuccess;
  }

  void FlowBoundaryCondition::Log(const Equations& state, double time) {
    for (const auto& stream : streams_) {
      stream->Log(state, time, location_);
    }
  }

  void FlowBoundaryCondition::ClearLog() {
    for (const auto& stream : streams_) {
      stream->stream_data_.Clear();
    }
  }

  // Closed Boundary Conditions
  void ClosedBoundaryCondition::Apply(Equations& state, double time) {
    auto& molar_concentration_map = state.molar_concentration;
    
    auto& mole_frac = state.mole_fraction;
    
    auto& temperature = state.fluid_temperature;
    auto& T = temperature.var;
    auto& res_T = temperature.res_var;
    
    auto& pressure = state.pressure;
    auto& P = pressure.var;
    auto& res_P = pressure.res_var;

    auto N = (int)T.size() - 1;

    int i = 0;
    int d = 1;
    if (location_ == BoundaryLocation::TOP) {
     i = N;
     d = -1; 
    }
    
    for (auto& [comp, molar_concentration] : molar_concentration_map) {
      auto& Yi = mole_frac.at(comp).var;
      auto& res_Yi = mole_frac.at(comp).res_var;
      // auto& Ci = molar_concentration.var;
      // auto& res_Ci = molar_concentration.res_var;
      // res_Ci[i] = Ci[i] - Yi[i + d] * C[i];

      res_Yi[i] = Yi[i] - Yi[i + d];
    }
    
    // res_C[i] = C[i] - C[i + d];
    res_P[i] = P[i] - P[i + d];
    res_T[i] = T[i] - T[i + d];

    // Ergun layer: no pressure drop across a closed end, so its boundary-face
    // velocity is imposed as zero (see ImposeBoundaryFaceVelocity).
    ImposeBoundaryFaceVelocity(state, location_, 0.0);
  }

  bool ClosedBoundaryCondition::PerformChecks(std::string& error_msg, const Fluid& fluid) {
    return true;
  }

  void ClosedBoundaryCondition::Log(const Equations& state, double time) {}

  void ClosedBoundaryCondition::ClearLog() {}

} // end adsim namespace


// Register the derived type with the base type
CEREAL_REGISTER_TYPE(adsim::PressureBoundaryCondition)
CEREAL_REGISTER_TYPE(adsim::FlowBoundaryCondition)
CEREAL_REGISTER_TYPE(adsim::ClosedBoundaryCondition)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IBoundaryCondition, adsim::PressureBoundaryCondition)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IBoundaryCondition, adsim::FlowBoundaryCondition)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IBoundaryCondition, adsim::ClosedBoundaryCondition)
CEREAL_CLASS_VERSION(adsim::PressureBoundaryCondition, 5);
CEREAL_CLASS_VERSION(adsim::FlowBoundaryCondition, 1);
CEREAL_CLASS_VERSION(adsim::ClosedBoundaryCondition, 1);
