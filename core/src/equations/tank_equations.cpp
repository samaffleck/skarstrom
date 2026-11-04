/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Residuals of the dynamic tanks (see adsim/objects/tank.h and docs/tanks.md).
//
// Species moles:  dn_i/dt = sum_p F_p,i,
// energy:         d/dt( sum_i n_i c_v,i T ) = sum_p sum_i F_p,i c_p,i T_p + Q,
//
// with F_p,i the molar flow of component i into the tank through p (a bed
// port or an external stream), T_p the temperature of the gas that crosses p
// and c_v,i = c_p,i - R (ideal gas, as the bed). Written in the internal
// energy, the balance is exactly conservative; it reduces to
//   N c_v dT/dt = sum_in F c_p (T_in - T) + R T dN/dt + Q.
//
// A live bed port takes exactly the flow the bed's boundary condition passes:
// a flow-boundary stream's own GetComponentMolarFlow (the function the bed's
// face flux is set from), a pressure boundary's face flux (face velocity,
// face concentration and the face mole fraction the species balance uses).
// So whatever leaves the bed enters the tank, including at flow reversal.

#include "../pch.h"
#include "adsim/equations/equation.h"
#include "adsim/objects/user_data.h"
#include "adsim/objects/tank.h"
#include "adsim/boundary_conditions/step.h"

#include <algorithm>

namespace adsim {

  void Equations::UpdateTanks(UserData* udata, double time) {
    if (tanks_.empty()) return;
    (void)udata;
    const auto& P = pressure.var;
    const auto& T = fluid_temperature.var;
    const auto& u = velocity.var;
    const auto& area = domain_vectors.area;
    const int N = static_cast<int>(T.size()) - 1;

    static thread_local std::vector<double> net;  // net flow into the tank per component

    for (auto& tk_ptr : tanks_) {
      TankUnknowns& tk = *tk_ptr;
      const Tank& tank = *tk.tank;
      const size_t nc = tk.n.size();
      const double T_tank = tk.live.temperature;
      net.assign(nc, 0.0);
      std::fill(tk.flow_in.begin(), tk.flow_in.end(), 0.0);
      std::fill(tk.flow_out.begin(), tk.flow_out.end(), 0.0);
      double H = 0.0;  // enthalpy flow into the tank [W]

      // One component's flow F into the tank, carrying gas at Tp.
      auto add = [&](size_t c, double F, double Tp) {
        net[c] += F;
        if (F >= 0.0) tk.flow_in[c] += F;
        else tk.flow_out[c] -= F;
        H += F * tk.cp[c] * Tp;
      };

      // Live ports of the current step.
      if (current_step_ && N >= 1) {
        for (IBoundaryCondition* bc : {current_step_->bottom_boundary_condition_.get(),
                                       current_step_->top_boundary_condition_.get()}) {
          if (!bc) continue;
          const BoundaryLocation loc = bc->location_;
          const int i = (loc == BoundaryLocation::TOP) ? N : 0;  // ghost cell
          const int d = (loc == BoundaryLocation::TOP) ? -1 : 1;
          // Axial flow (positive upwards) leaves the bed into the tank at the
          // top and enters the bed from it at the bottom.
          const double sign = (loc == BoundaryLocation::TOP) ? 1.0 : -1.0;
          const double T_end = 0.5 * (T[i] + T[i + d]);

          if (auto* fbc = dynamic_cast<FlowBoundaryCondition*>(bc)) {
            for (auto& stream : fbc->streams_) {
              if (!stream || stream->tank_port_.tank_id_ != tk.tank_id) continue;
              // Gas crossing: the tank's for an inlet, the bed end's for an outlet.
              const double Tp = (bc->direction_ == BoundaryDirection::INLET) ? stream->InletTemperature(time) : T_end;
              for (size_t c = 0; c < nc; ++c)
                add(c, sign * stream->GetComponentMolarFlow(*this, time, loc, tk.components[c]), Tp);
            }
          } else if (auto* pbc = dynamic_cast<PressureBoundaryCondition*>(bc)) {
            if (pbc->tank_port_.tank_id_ != tk.tank_id || u.empty()) continue;
            const int ui = (loc == BoundaryLocation::TOP) ? static_cast<int>(u.size()) - 1 : 0;
            const double Pf = 0.5 * (P[i] + P[i + d]);
            const double Tf = detail::DenomFloor(T_end, constants::kMinDenomTemperature);
            const double axial = area[i] * u[ui] * Pf / (R * Tf);
            for (size_t c = 0; c < nc; ++c)
              add(c, sign * axial * mole_fraction_face.at(tk.components[c]).var[ui], Tf);
          }
        }
      }

      // External streams.
      for (const auto& stream : tank.streams_) {
        if (!stream) continue;
        const double F = stream->FlowIn(tk.live);
        if (F > 0.0 && stream->CanFlowIn()) {
          const double Tp = tk.T_unknown ? stream->temperature_.Get() : T_tank;
          for (size_t c = 0; c < nc; ++c) {
            const auto it = stream->mole_fraction_.find(tk.components[c]);
            add(c, F * (it == stream->mole_fraction_.end() ? 0.0 : it->second.Get()), Tp);
          }
        } else {
          for (size_t c = 0; c < nc; ++c) add(c, F * tk.live.MoleFraction(tk.components[c]), T_tank);
        }
      }

      // Heat to the gas.
      double Q = 0.0;
      if (tk.T_unknown) {
        if (tank.thermal_ == TankThermalMode::HEAT_LOSS) {
          Q = tank.ua_ * (tank.ambient_temperature_ - tk.T);
        } else if (tank.thermal_ == TankThermalMode::WALL) {
          Q = tank.h_in_ * tank.wall_area_ * (tk.Tw - tk.T);
        }
      }
      tk.heat = Q;

      for (size_t c = 0; c < nc; ++c) tk.res_n[c] = tk.dn[c] - net[c];

      if (tk.T_unknown) {
        double dU = 0.0;  // d/dt (sum_i n_i c_v,i T)
        for (size_t c = 0; c < nc; ++c) {
          const double cv = tk.cp[c] - R;
          dU += cv * (tk.dn[c] * tk.T + tk.n[c] * tk.dT);
        }
        tk.res_T = dU - H - Q;
      }
      if (tk.Tw_unknown) {
        const double A = tank.wall_area_;
        tk.res_Tw = tank.wall_heat_capacity_ * tk.dTw
                  - tank.h_in_ * A * (tk.T - tk.Tw)
                  + tank.h_out_ * A * (tk.Tw - tank.ambient_temperature_);
      }
    }
  }

} // end adsim namespace
