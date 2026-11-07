/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// KPI quadratures: the time integrals behind the KPIs integrated by IDAS as
// quadrature states alongside the solution, instead of the trapezoid over the
// logged boundary samples.
//
// The integrands are what the stream log records at each logged instant (each
// endpoint's EndpointSample) plus the species flux through the two bed faces,
// evaluated at IDA's corrected solution of every step, so the integrals carry
// the integrator's own accuracy rather than the log's sampling. Quadrature
// error control is off by default (IDASettings::quadrature_error_control_):
// the quadratures then do not enter IDA's error test or step-size choice, and
// the solution, the logs and the reported (trapezoid) KPIs are exactly those of
// a run without them.
//
// Each cycle step starts the integrals from zero. The vector layout is the same
// for every step: the face block first, then one endpoint block per logged
// endpoint of the current step (Step::LoggedEndpoints order; unused blocks
// stay zero). The layout is a plain vector of integrals so that quadrature
// sensitivities (IDAQuadSensInit) could be added over the same states later.

#include <string>
#include <vector>

namespace adsim {

  class UserData;
  class Step;

  struct KpiQuadratureLayout {
    int n_components = 0;
    int max_endpoints = 0;  // most logged endpoints of any step of the cycle

    // Face block: species flux entering the bed through the bottom (end 0) and
    // top (end 1) face, per component.
    int Face(int end, int component) const { return end * n_components + component; }
    // Endpoint block k: moles in, moles out, component in (n), component out
    // (n), pump energy, integral of |molar flow|, integral of T |molar flow|.
    int BlockSize() const { return 5 + 2 * n_components; }
    int Block(int k) const { return 2 * n_components + k * BlockSize(); }
    static int In() { return 0; }
    static int Out() { return 1; }
    int ComponentIn(int c) const { return 2 + c; }
    int ComponentOut(int c) const { return 2 + n_components + c; }
    int Energy() const { return 2 + 2 * n_components; }
    int AbsMoles() const { return 3 + 2 * n_components; }
    int TemperatureMoles() const { return 4 + 2 * n_components; }
    int Size() const { return 2 * n_components + max_endpoints * BlockSize(); }

    // The layout for the deck's components and cycle.
    static KpiQuadratureLayout For(const UserData& user_data);
  };

  // Writes the integrands at the current state into out[0 .. layout.Size()).
  // The state must be unpacked with its independent and diagnostic equations
  // updated (as for logging); t_relative is the time since the step start.
  void EvaluateKpiQuadratureIntegrands(UserData& user_data, double t_relative,
                                       const KpiQuadratureLayout& layout, double* out);

  // Gas plus adsorbed amount of each component in the column [mol], in
  // fluid.components_ order, from the current state: the sum over the bed cells
  // of V (et y P / (R T) + sum_k rho_b,k q_k), the quantity whose time
  // derivative the species balance holds.
  std::vector<double> ColumnInventory(const UserData& user_data);

  // Mole balance of one cycle step from the quadratures and the column
  // inventory at the step's start and end (the per-step accumulation form of
  // Glover, Papathanasiou & Pini 2026): what entered through the endpoints (or
  // through the faces) must equal the change of the inventory. Per component,
  // in `components` order.
  struct StepBalance {
    bool valid = false;
    int open_ends = 0;  // bed ends not closed (Step::NumberOfOpenEnds)
    std::vector<std::string> components;
    std::vector<double> inventory_start;   // [mol]
    std::vector<double> inventory_end;     // [mol]
    std::vector<double> face_in_bottom;    // [mol] species flux into the bed through the bottom face
    std::vector<double> face_in_top;       // [mol] ... through the top face
    std::vector<double> endpoint_in;       // [mol] entering through the step's endpoints (gross)
    std::vector<double> endpoint_out;      // [mol] leaving through the step's endpoints (gross)

    double InventoryChange(size_t c) const { return inventory_end[c] - inventory_start[c]; }
    double EndpointNetIn(size_t c) const { return endpoint_in[c] - endpoint_out[c]; }
    double FaceNetIn(size_t c) const { return face_in_bottom[c] + face_in_top[c]; }
    // Totals over the components [mol].
    double TotalInventoryChange() const;
    double TotalEndpointNetIn() const;
    double TotalFaceNetIn() const;
    double TotalThroughput() const;  // gross moles in + out through the endpoints
    // Endpoint (face) moles in minus the inventory change [mol]: zero for an
    // exactly conservative integration. With one open end this is the
    // inventory-difference cross-check of that end's outflow.
    double EndpointDiscrepancy() const { return TotalEndpointNetIn() - TotalInventoryChange(); }
    double FaceDiscrepancy() const { return TotalFaceNetIn() - TotalInventoryChange(); }
    // EndpointDiscrepancy relative to the larger of |inventory change| and the
    // endpoint throughput (0 when both are 0).
    double RelativeEndpointDiscrepancy() const;
  };

  // Writes the step's quadrature vector `q` into the integrals of each logged
  // endpoint of `step` (StreamData::integrals_) and the face and endpoint
  // totals into `balance`. `valid` false marks the integrals invalid.
  void DistributeKpiQuadratures(Step& step, const KpiQuadratureLayout& layout,
                                const std::vector<std::string>& components,
                                const std::vector<double>& q, bool valid, StepBalance& balance);

}  // namespace adsim
