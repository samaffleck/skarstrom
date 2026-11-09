/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// A minimal single-adsorbent O2/N2 pressure-swing cycle (pressurise, feed,
// blowdown) on a short 13X bed, built from scratch like the breakthrough case
// so every value is pinned in one visible place. Shared by the tests that need
// a cycle that genuinely separates and reverses flow: cyclic-steady-state
// convergence, CSS acceleration and the bounds invariants.
//
// The caller sets the number of cycles and the solver tolerances. A single
// adsorbent is deliberate: CSS acceleration only re-seeds the adsorbed loading
// of a single-adsorbent bed. AddO2ProductKpis turns on the KPIs of the O2
// product (with a feed compressor, so the specific energy is defined) for the
// grid- and CSS-convergence tests.

#include "validation_helpers.h"

namespace adsim::validation {

  // n_cells: the grid (coarse grids keep the many-cycle runs cheap);
  // mtc: the linear-driving-force coefficient [1/s] of both components.
  inline void BuildPsaCycle(UserData& data, int n_cells = 30, double mtc = 0.5) {
    data.bed.wall_.thickness_ = 0.0175;
    data.bed.wall_.heat_capacity_ = 537;
    data.bed.wall_.density_ = 7800;
    data.bed.wall_.thermal_conductivity_ = 16.2;

    data.fluid.AddComponent("Nitrogen");
    data.fluid.AddComponent("Oxygen");

    auto ads = data.adsorbent_directory.Add("13X");
    auto lid = data.bed.layers_.Add("13X");
    auto& layer = data.bed.layers_.Get(lid);
    layer.momentum_model_ = MomentumModel::KarmenCozeny;  // fixture predates Ergun: keep its Darcy / Kozeny-Carman pressure drop
    layer.domain_.length = 1.0;
    layer.domain_.n_cells = n_cells;
    layer.eb_ = 0.4;
    layer.particle_diameter_ = 0.002;
    layer.composition_[ads] = 1.0;
    std::get<CylinderBed>(data.bed.shape_).diameter = 0.1;

    auto& adsorbent = data.adsorbent_directory.Get(ads);
    adsorbent.ep_ = 0.35;
    adsorbent.solid_density_ = 1130.0;
    adsorbent.heat_capacity_ = 1070.0;

    // Both components adsorb, with different affinities, so the cycle genuinely
    // separates the feed.
    auto isotherm = std::make_unique<LangmuirIsotherm>();
    isotherm->params_["Oxygen"].m1 = 3.0;
    isotherm->params_["Oxygen"].b1 = 3.0e-11;
    isotherm->params_["Oxygen"].q1 = 20000.0;
    isotherm->params_["Nitrogen"].m1 = 3.0;
    isotherm->params_["Nitrogen"].b1 = 5.0e-11;
    isotherm->params_["Nitrogen"].q1 = 24000.0;
    adsorbent.adsorption_source_.AddIsotherm({"Oxygen", "Nitrogen"}, std::move(isotherm));
    for (const auto& component : {"Oxygen", "Nitrogen"}) {
      auto ldf = std::make_unique<LinearDrivingForce>();
      auto mtc_model = std::make_unique<MassTransferCoefficientConstant>();
      mtc_model->constant_mtc_ = mtc;
      ldf->mtc_model_ = ToResistances(std::move(mtc_model));  // generic model, as the example decks store it
      adsorbent.adsorption_source_.mass_transfer_models_[component] = std::move(ldf);
      adsorbent.adsorption_source_.heat_of_adsorption_[component] =
        (std::string(component) == "Nitrogen") ? 24000.0 : 20000.0;
    }

    auto ihtc = std::make_unique<ConstantInternalHeatTransferModel>();
    ihtc->heat_transfer_coefficient_ = 20.0;
    layer.equations_.internal_heat_transfer_model_ = std::move(ihtc);
    auto ehtc = std::make_unique<ConstantExternalHeatTransferModel>();
    ehtc->heat_transfer_coefficient_ = 10.0;
    layer.equations_.external_heat_transfer_model_ = std::move(ehtc);

    data.bed.ambient_temperature_ = 298.0;
    data.ic.pressure = 1.0e5;
    data.ic.temperature = 298.0;
    data.ic.wall_temperature = 298.0;
    data.ic.velocity = 0.0;
    data.ic.mole_fraction["Nitrogen"] = 0.79;
    data.ic.mole_fraction["Oxygen"] = 0.21;

    auto make_flow_inlet = [&](double molar_flow, BoundaryLocation location) {
      auto stream = std::make_unique<BoundaryStream>();
      auto flow = std::make_unique<MolarFlowBoundaryValue>();
      static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = molar_flow;
      stream->boundary_flow_ = std::move(flow);
      stream->mole_fraction["Nitrogen"] = 0.79;
      stream->mole_fraction["Oxygen"] = 0.21;
      stream->temperature_ = 298.0;
      auto bc = std::make_unique<FlowBoundaryCondition>();
      bc->direction_ = BoundaryDirection::INLET;
      bc->location_ = location;
      bc->streams_.push_back(std::move(stream));
      return bc;
    };
    auto make_pressure_outlet = [&](double pressure, BoundaryLocation location) {
      auto bc = std::make_unique<PressureBoundaryCondition>();
      bc->direction_ = BoundaryDirection::OUTLET;
      bc->location_ = location;
      auto value = std::make_unique<ConstantBoundaryValue>();
      value->constant_value_ = pressure;
      bc->pressure_ = std::move(value);
      bc->mole_fraction_["Nitrogen"] = 0.79;
      bc->mole_fraction_["Oxygen"] = 0.21;
      return bc;
    };
    auto make_closed = [&](BoundaryLocation location) {
      auto bc = std::make_unique<ClosedBoundaryCondition>();
      bc->location_ = location;
      return bc;
    };

    // Pressurise (closed top) / feed (open top) / blow down (reverse flow out
    // the bottom) -- a minimal Skarstrom-like cycle.
    data.cycle.AddStep("Pressurise");
    Step& pressurise = data.cycle.GetStep("Pressurise");
    pressurise.step_time_ = 30.0;
    pressurise.bottom_boundary_condition_ = make_flow_inlet(0.02, BoundaryLocation::BOTTOM);
    pressurise.top_boundary_condition_ = make_closed(BoundaryLocation::TOP);

    data.cycle.AddStep("Feed");
    Step& feed = data.cycle.GetStep("Feed");
    feed.step_time_ = 120.0;
    feed.bottom_boundary_condition_ = make_flow_inlet(0.02, BoundaryLocation::BOTTOM);
    feed.top_boundary_condition_ = make_pressure_outlet(3.0e5, BoundaryLocation::TOP);

    data.cycle.AddStep("Blowdown");
    Step& blowdown = data.cycle.GetStep("Blowdown");
    blowdown.step_time_ = 120.0;
    blowdown.bottom_boundary_condition_ = make_pressure_outlet(1.0e5, BoundaryLocation::BOTTOM);
    blowdown.top_boundary_condition_ = make_closed(BoundaryLocation::TOP);
  }

  // The O2 product KPIs of the cycle above: feed = the bottom inlet of
  // Pressurise and Feed, product = the top outlet of Feed. A compressor on the
  // feed stream (suction 1 bar) makes the specific energy defined.
  inline void AddO2ProductKpis(UserData& data) {
    for (const char* step : {"Pressurise", "Feed"}) {
      auto* bc = dynamic_cast<FlowBoundaryCondition*>(
          data.cycle.GetStep(step).bottom_boundary_condition_.get());
      for (auto& stream : bc->streams_) {
        auto* s = dynamic_cast<BoundaryStream*>(stream.get());
        s->pump_ = CreatePump(PumpType::COMPRESSOR);
        s->pump_->use_pump_ = true;
      }
    }
    data.kpi_config.SyncToComponents(data.fluid.components_);
    ComponentKPIConfig* kpi = data.kpi_config.FindByComponent("Oxygen");
    kpi->tracked = true;
    for (int s : {0, 1})
      for (const auto& ep : data.cycle.step_list_[s].GetEndpoints())
        if (ep.location == BoundaryLocation::BOTTOM) kpi->feed_endpoints.push_back({s, ep.name, ep.location});
    for (const auto& ep : data.cycle.step_list_[1].GetEndpoints())
      if (ep.location == BoundaryLocation::TOP) kpi->product_endpoints.push_back({1, ep.name, ep.location});
  }

  struct O2Kpis {
    double purity = 0.0;           // [mol/mol]
    double recovery = 0.0;         // [mol/mol]
    double specific_energy = 0.0;  // [J/kg O2]
    bool defined = false;
  };

  // The O2 product KPIs of one logged cycle (AddO2ProductKpis).
  inline O2Kpis O2KpisAt(const UserData& data, const EquationLog& log, const BoundaryLog& b_log,
                         int cycle) {
    const ComponentKPIConfig* cfg = data.kpi_config.FindByComponent("Oxygen");
    KPIParams p(log, b_log);
    p.cycle_index_ = cycle;
    p.product_comp_ = "Oxygen";
    p.feed_endpoints = cfg->feed_endpoints;
    p.product_endpoints = cfg->product_endpoints;
    bool ok_p = false, ok_r = false, ok_e = false;
    O2Kpis k;
    k.purity = CalculatePurity(p, &ok_p);
    k.recovery = CalculateRecovery(p, &ok_r);
    k.specific_energy = CalculateSpecificEnergyConsumption(p, &ok_e);
    k.defined = ok_p && ok_r && ok_e;
    return k;
  }

} // namespace adsim::validation
