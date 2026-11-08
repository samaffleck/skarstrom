/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/solver/kpi_quadrature.h"

#include <algorithm>
#include <cmath>

#include "adsim/math/constants.h"
#include "adsim/objects/user_data.h"
#include "adsim/boundary_conditions/step.h"

namespace adsim {

  KpiQuadratureLayout KpiQuadratureLayout::For(const UserData& user_data) {
    KpiQuadratureLayout layout;
    layout.n_components = static_cast<int>(user_data.fluid.components_.size());
    for (const auto& step : user_data.cycle.step_list_) {
      // LoggedEndpoints is non-const only because it hands out mutable
      // pointers; counting does not change the step.
      const int n = static_cast<int>(const_cast<Step&>(step).LoggedEndpoints().size());
      layout.max_endpoints = std::max(layout.max_endpoints, n);
    }
    return layout;
  }

  void EvaluateKpiQuadratureIntegrands(UserData& user_data, double t_relative,
                                       const KpiQuadratureLayout& layout, double* out) {
    const int size = layout.Size();
    std::fill(out, out + size, 0.0);
    auto& state = user_data.state;
    const auto& components = user_data.fluid.components_;
    const int nc = layout.n_components;

    // Species flux through the two bed faces as the species balance sees it:
    // the face molar flow times the face mole fraction (the convective flux
    // alone; no dispersion crosses a bed face).
    const auto& n_face = state.molar_flow_rate.var;
    if (!n_face.empty()) {
      const size_t top = n_face.size() - 1;
      for (int c = 0; c < nc; ++c) {
        const auto it = state.mole_fraction_face.find(components[c]);
        if (it == state.mole_fraction_face.end() || it->second.var.size() != n_face.size()) continue;
        const auto& y_face = it->second.var;
        out[layout.Face(0, c)] = n_face[0] * y_face[0];
        out[layout.Face(1, c)] = -n_face[top] * y_face[top];
      }
    }

    // Endpoints of the current step, as the stream log samples them.
    // CurrentStep is the step the residual applies (Equations only hands out
    // a const pointer; sampling evaluates its flow laws, as logging does).
    Step* step = const_cast<Step*>(state.CurrentStep());
    if (!step) return;
    static thread_local std::vector<EndpointSample> samples;
    samples.clear();
    step->SampleEndpoints(state, t_relative, samples);
    static thread_local std::vector<Step::LoggedEndpoint> endpoints;
    endpoints = step->LoggedEndpoints();
    const int n_ep = std::min<int>(static_cast<int>(samples.size()), layout.max_endpoints);
    for (int k = 0; k < n_ep; ++k) {
      const EndpointSample& s = samples[k];
      const double sign = (endpoints[k].location == BoundaryLocation::BOTTOM) ? 1.0 : -1.0;
      double* b = out + layout.Block(k);
      const double entering = sign * s.molar_flow;
      b[KpiQuadratureLayout::In()] = std::max(entering, 0.0);
      b[KpiQuadratureLayout::Out()] = std::max(-entering, 0.0);
      for (int c = 0; c < nc; ++c) {
        const auto it = s.component_molar_flow.find(components[c]);
        const double e = (it == s.component_molar_flow.end()) ? 0.0 : sign * it->second;
        b[layout.ComponentIn(c)] = std::max(e, 0.0);
        b[layout.ComponentOut(c)] = std::max(-e, 0.0);
      }
      b[layout.Energy()] = s.power;
      b[layout.AbsMoles()] = std::abs(s.molar_flow);
      b[layout.TemperatureMoles()] = s.temperature * std::abs(s.molar_flow);
    }
  }

  std::vector<double> ColumnInventory(const UserData& user_data) {
    const auto& state = user_data.state;
    const auto& components = user_data.fluid.components_;
    const auto& dv = state.domain_vectors;
    const auto& P = state.pressure.var;
    const auto& T = state.fluid_temperature.var;
    const double R = constants::kIdealGasConstant;
    std::vector<double> inventory(components.size(), 0.0);
    const int N = static_cast<int>(P.size());
    for (size_t c = 0; c < components.size(); ++c) {
      const auto y_it = state.mole_fraction.find(components[c]);
      if (y_it == state.mole_fraction.end()) continue;
      const auto& y = y_it->second.var;
      double total = 0.0;
      for (int n = 1; n < N - 1; ++n) {  // interior cells; 0 and N-1 are ghosts
        double amount = dv.et[n] * y[n] * P[n] / (R * T[n]);
        if (n < static_cast<int>(dv.cell_adsorbents.size())) {
          for (const auto& ca : dv.cell_adsorbents[n]) {
            const auto ads = state.ads_concentration.find(ca.id);
            if (ads == state.ads_concentration.end()) continue;
            const auto q = ads->second.find(components[c]);
            if (q == ads->second.end()) continue;
            amount += ca.bulk_density * q->second.var[n];
          }
        }
        total += dv.volume[n] * amount;
      }
      inventory[c] = total;
    }
    return inventory;
  }

  namespace {
    double Sum(const std::vector<double>& v) {
      double s = 0.0;
      for (double x : v) s += x;
      return s;
    }
  }

  double StepBalance::TotalInventoryChange() const { return Sum(inventory_end) - Sum(inventory_start); }
  double StepBalance::TotalEndpointNetIn() const { return Sum(endpoint_in) - Sum(endpoint_out); }
  double StepBalance::TotalFaceNetIn() const { return Sum(face_in_bottom) + Sum(face_in_top); }
  double StepBalance::TotalThroughput() const { return Sum(endpoint_in) + Sum(endpoint_out); }

  double StepBalance::RelativeEndpointDiscrepancy() const {
    const double scale = std::max(std::abs(TotalInventoryChange()), TotalThroughput());
    return scale > 0.0 ? EndpointDiscrepancy() / scale : 0.0;
  }

  void DistributeKpiQuadratures(Step& step, const KpiQuadratureLayout& layout,
                                const std::vector<std::string>& components,
                                const std::vector<double>& q, bool valid, StepBalance& balance) {
    const int nc = layout.n_components;
    balance.components = components;
    balance.face_in_bottom.assign(nc, 0.0);
    balance.face_in_top.assign(nc, 0.0);
    balance.endpoint_in.assign(nc, 0.0);
    balance.endpoint_out.assign(nc, 0.0);
    balance.open_ends = step.NumberOfOpenEnds();
    valid = valid && static_cast<int>(q.size()) == layout.Size();
    const auto endpoints = step.LoggedEndpoints();
    for (int k = 0; k < static_cast<int>(endpoints.size()); ++k) {
      EndpointIntegrals ints;
      if (valid && k < layout.max_endpoints) {
        const double* b = q.data() + layout.Block(k);
        ints.valid = true;
        ints.moles_in = b[KpiQuadratureLayout::In()];
        ints.moles_out = b[KpiQuadratureLayout::Out()];
        for (int c = 0; c < nc; ++c) {
          ints.component_in[components[c]] = b[layout.ComponentIn(c)];
          ints.component_out[components[c]] = b[layout.ComponentOut(c)];
          balance.endpoint_in[c] += b[layout.ComponentIn(c)];
          balance.endpoint_out[c] += b[layout.ComponentOut(c)];
        }
        ints.energy = b[layout.Energy()];
        ints.abs_axial_moles = b[layout.AbsMoles()];
        ints.temperature_weighted_moles = b[layout.TemperatureMoles()];
      }
      endpoints[k].data->integrals_ = std::move(ints);
    }
    if (valid) {
      for (int c = 0; c < nc; ++c) {
        balance.face_in_bottom[c] = q[layout.Face(0, c)];
        balance.face_in_top[c] = q[layout.Face(1, c)];
      }
    }
    balance.valid = valid && balance.inventory_start.size() == static_cast<size_t>(nc)
                    && balance.inventory_end.size() == static_cast<size_t>(nc);
  }

}  // namespace adsim
