/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// A single, deliberately plain isothermal breakthrough case, built from scratch
// rather than on the shared UserDataTest fixture so that each test can vary
// exactly one thing and everything else is pinned in one visible place.
//
// Physical dispersion, wall heat transfer and axial conduction are all switched
// off by default: what is left is convection, the isotherm (optionally), and the
// discretisation itself.

#include "validation_helpers.h"

namespace adsim::validation {

  struct BreakthroughCase {
    int n_cells = 80;
    double length = 1.0;          // [m]
    double diameter = 0.1;        // [m]
    double molar_flow = 0.005;    // [mol/s]
    double run_time = 400.0;      // [s]
    double reltol = 1e-9;
    double max_step = 0.5;
    double y_feed = 0.05;
    double feed_temperature = 298.0;
    double bed_temperature = 298.0;
    bool feed_from_top = false;
    bool adsorbing = false;
  };

  inline constexpr double kCasePressure = 1.0e5;      // [Pa]
  inline constexpr double kCaseVoidFraction = 0.4;    // [-]
  inline constexpr double kCaseSolidDensity = 1000.0; // [kg/m3]
  inline constexpr double kCaseHenryK = 4.0e-10;      // [mol/(kg-Pa)]
  inline constexpr double kCaseHenryQ = 20000.0;      // [J/mol]

  inline void Build(UserData& data, const BreakthroughCase& c) {
    data.bed.wall_.thickness_ = 0.0175;
    data.bed.wall_.heat_capacity_ = 537;
    data.bed.wall_.density_ = 7800;
    data.bed.wall_.thermal_conductivity_ = 16.2;

    data.fluid.AddComponent("Nitrogen");
    data.fluid.AddComponent("Oxygen");

    auto ads = data.adsorbent_directory.Add("Packing");
    auto lid = data.bed.layers_.Add("Packing");
    auto& layer = data.bed.layers_.Get(lid);
    layer.momentum_model_ = MomentumModel::KarmenCozeny;  // fixture predates Ergun: keep its Darcy / Kozeny-Carman pressure drop
    layer.domain_.length = c.length;
    layer.domain_.n_cells = c.n_cells;
    layer.eb_ = kCaseVoidFraction;
    layer.particle_diameter_ = 0.002;
    layer.composition_[ads] = 1.0;
    std::get<CylinderBed>(data.bed.shape_).diameter = c.diameter;

    auto& adsorbent = data.adsorbent_directory.Get(ads);
    adsorbent.ep_ = 0.0;
    adsorbent.solid_density_ = kCaseSolidDensity;
    adsorbent.heat_capacity_ = 900.0;
    adsorbent.thermal_conductivity_ = 0.0;

    if (c.adsorbing) {
      auto isotherm = std::make_unique<HenryIsotherm>();
      isotherm->params_["Oxygen"].k = kCaseHenryK;
      isotherm->params_["Oxygen"].Q = kCaseHenryQ;
      adsorbent.adsorption_source_.AddIsotherm({"Oxygen"}, std::move(isotherm));
      auto ldf = std::make_unique<LinearDrivingForce>();
      auto mtc = std::make_unique<MassTransferCoefficientConstant>();
      mtc->constant_mtc_ = 0.5;
      ldf->mtc_model_ = std::move(mtc);
      adsorbent.adsorption_source_.mass_transfer_models_["Oxygen"] = std::move(ldf);
      adsorbent.adsorption_source_.heat_of_adsorption_["Oxygen"] = 0.0;
    }

    // Isothermal and dispersion-free: the only smearing left is numerical.
    auto ihtc = std::make_unique<ConstantInternalHeatTransferModel>();
    ihtc->heat_transfer_coefficient_ = 0.0;
    layer.equations_.internal_heat_transfer_model_ = std::move(ihtc);
    auto ehtc = std::make_unique<ConstantExternalHeatTransferModel>();
    ehtc->heat_transfer_coefficient_ = 0.0;
    layer.equations_.external_heat_transfer_model_ = std::move(ehtc);
    auto disp = std::make_unique<DispersionModelConstant>();
    disp->dispersion_coefficient_ = 0.0;
    layer.equations_.dispersion_model_ = std::move(disp);
    auto heat_disp = std::make_unique<HeatDispersionModelConstant>();
    heat_disp->dispersion_coefficient_ = 0.0;
    layer.equations_.heat_dispersion_model_ = std::move(heat_disp);

    data.ic.pressure = kCasePressure;
    data.ic.temperature = c.bed_temperature;
    data.ic.wall_temperature = c.bed_temperature;
    data.ic.velocity = 0.0;
    data.ic.mole_fraction["Nitrogen"] = 1.0;
    data.ic.mole_fraction["Oxygen"] = 0.0;

    data.cycle.AddStep("Feed");
    Step& step = data.cycle.GetStep("Feed");
    step.step_time_ = c.run_time;

    auto inlet = std::make_unique<BoundaryStream>();
    auto flow = std::make_unique<MolarFlowBoundaryValue>();
    static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = c.molar_flow;
    inlet->boundary_flow_ = std::move(flow);
    inlet->mole_fraction["Nitrogen"] = 1.0 - c.y_feed;
    inlet->mole_fraction["Oxygen"] = c.y_feed;
    inlet->temperature_ = c.feed_temperature;

    const BoundaryLocation inlet_location =
      c.feed_from_top ? BoundaryLocation::TOP : BoundaryLocation::BOTTOM;
    const BoundaryLocation outlet_location =
      c.feed_from_top ? BoundaryLocation::BOTTOM : BoundaryLocation::TOP;

    auto inlet_bc = std::make_unique<FlowBoundaryCondition>();
    inlet_bc->direction_ = BoundaryDirection::INLET;
    inlet_bc->location_ = inlet_location;
    inlet_bc->streams_.push_back(std::move(inlet));

    auto outlet_bc = std::make_unique<PressureBoundaryCondition>();
    outlet_bc->direction_ = BoundaryDirection::OUTLET;
    outlet_bc->location_ = outlet_location;
    auto outlet_p = std::make_unique<ConstantBoundaryValue>();
    outlet_p->constant_value_ = kCasePressure;
    outlet_bc->pressure_ = std::move(outlet_p);
    outlet_bc->mole_fraction_["Nitrogen"] = 1.0 - c.y_feed;
    outlet_bc->mole_fraction_["Oxygen"] = c.y_feed;

    if (c.feed_from_top) {
      step.top_boundary_condition_ = std::move(inlet_bc);
      step.bottom_boundary_condition_ = std::move(outlet_bc);
    } else {
      step.bottom_boundary_condition_ = std::move(inlet_bc);
      step.top_boundary_condition_ = std::move(outlet_bc);
    }

    data.cycle.number_of_cycles_ = 1;
    data.settings.reltol = c.reltol;
    data.settings.max_step = c.max_step;
  }

  // Outlet mole fraction series, taken from whichever end is downstream.
  inline std::pair<std::vector<double>, std::vector<double>> OutletSeries(
      const EquationLog& log, bool feed_from_top) {
    const auto y = log.GetComponentFaceVariable(0, 0, ComponentVariableType::MOLE_FRACTION,
                                                "Oxygen");
    std::vector<double> series;
    series.reserve(y.values.size());
    for (const auto& snapshot : y.values) {
      if (snapshot.empty()) { series.push_back(0.0); continue; }
      series.push_back(feed_from_top ? snapshot.front() : snapshot.back());
    }
    return {y.times, series};
  }

} // namespace adsim::validation
