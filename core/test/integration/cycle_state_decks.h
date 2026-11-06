// cycle_state_decks.h
//
// Small, fast decks shared by the cycle-state tests (round trip, grid remap,
// warm start). Both have a buffer-tank purge fed from the feed vent, with the
// consumer step before its source, so a restart must carry the previous
// cycle's summary.
#pragma once

#include <map>
#include <memory>
#include <string>

#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/boundary_value.h"
#include "adsim/boundary_conditions/buffer_tank_source.h"
#include "adsim/boundary_conditions/flow_boundary_value.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/equations/isotherm.h"
#include "adsim/objects/user_data.h"
#include "adsim/simulation.h"

namespace adsim::test_decks {

// ---- a buffer-tank deck (consumer before source) ---------------------------
// Purge feeds the top of the bed from a tank holding 60 % of the moles Feed
// vented through its top in its last run; Feed runs after Purge, so from cycle
// 2 on Purge replays the previous cycle's summary.

inline std::unique_ptr<ConstantBoundaryValue> Constant(double v) {
  auto bv = std::make_unique<ConstantBoundaryValue>();
  bv->constant_value_ = v;
  return bv;
}

inline std::unique_ptr<VelocityBoundaryValue> Velocity(double v) {
  auto vel = std::make_unique<VelocityBoundaryValue>();
  static_cast<ConstantBoundaryValue*>(vel->profile_.get())->constant_value_ = v;
  return vel;
}

inline std::unique_ptr<PressureBoundaryCondition> PressureOutlet(BoundaryLocation loc) {
  auto pbc = std::make_unique<PressureBoundaryCondition>();
  pbc->direction_ = BoundaryDirection::OUTLET;
  pbc->location_ = loc;
  pbc->pressure_ = Constant(101325.0);
  pbc->mole_fraction_["Oxygen"] = 0.21;
  pbc->mole_fraction_["Nitrogen"] = 0.79;
  return pbc;
}

inline std::unique_ptr<FlowBoundaryCondition> FlowInlet(BoundaryLocation loc, std::unique_ptr<BoundaryStream> stream) {
  auto fbc = std::make_unique<FlowBoundaryCondition>();
  fbc->direction_ = BoundaryDirection::INLET;
  fbc->location_ = loc;
  fbc->streams_.push_back(std::move(stream));
  return fbc;
}

inline UserData BufferTankDeck() {
  UserData data;
  data.bed.wall_.thickness_ = 0.0175;
  data.bed.wall_.heat_capacity_ = 537;
  data.bed.wall_.density_ = 7800;
  data.bed.wall_.thermal_conductivity_ = 16.2;

  data.fluid.AddComponent("Oxygen");
  data.fluid.AddComponent("Nitrogen");

  auto ads_id = data.adsorbent_directory.Add("13X");
  auto layer_id = data.bed.layers_.Add("13X");
  auto& layer = data.bed.layers_.Get(layer_id);
  layer.momentum_model_ = MomentumModel::KarmenCozeny;
  layer.domain_.length = 1.0;
  layer.domain_.n_cells = 10;
  layer.particle_diameter_ = 0.0011;
  layer.eb_ = 0.37;
  layer.composition_[ads_id] = 1.0;

  auto& ads = data.adsorbent_directory.Get(ads_id);
  ads.ep_ = 0.37;
  auto langmuir = std::make_unique<LangmuirIsotherm>();
  langmuir->params_["Oxygen"].m1 = 0.5;
  langmuir->params_["Oxygen"].b1 = 3e-11;
  langmuir->params_["Oxygen"].q1 = 30000;
  langmuir->params_["Nitrogen"].m1 = 1.5;
  langmuir->params_["Nitrogen"].b1 = 3e-11;
  langmuir->params_["Nitrogen"].q1 = 35000;
  ads.adsorption_source_.AddIsotherm({"Oxygen", "Nitrogen"}, std::move(langmuir));
  for (const char* comp : {"Oxygen", "Nitrogen"}) {
    auto ldf = std::make_unique<LinearDrivingForce>();
    auto mtc = std::make_unique<MassTransferCoefficientConstant>();
    mtc->constant_mtc_ = 0.1;
    ldf->mtc_model_ = std::move(mtc);
    ads.adsorption_source_.mass_transfer_models_[comp] = std::move(ldf);
  }
  auto ihtc = std::make_unique<ConstantInternalHeatTransferModel>();
  ihtc->heat_transfer_coefficient_ = 0.0;
  layer.equations_.internal_heat_transfer_model_ = std::move(ihtc);
  auto ehtc = std::make_unique<ConstantExternalHeatTransferModel>();
  ehtc->heat_transfer_coefficient_ = 0.0;
  layer.equations_.external_heat_transfer_model_ = std::move(ehtc);

  data.ic.pressure = 101325.0;
  data.ic.temperature = 298.0;
  data.ic.wall_temperature = 298.0;
  data.ic.velocity = 0.0;
  data.ic.mole_fraction["Oxygen"] = 0.21;
  data.ic.mole_fraction["Nitrogen"] = 0.79;

  data.cycle.AddStep("Purge");
  {
    Step& purge = data.cycle.GetStep("Purge");
    purge.step_time_ = 20.0;
    auto in = std::make_unique<BoundaryStream>();
    in->name_ = "Purge In";
    in->boundary_flow_ = Velocity(0.01);  // first-cycle assumption
    in->mole_fraction["Oxygen"] = 1.0;
    in->mole_fraction["Nitrogen"] = 0.0;
    in->temperature_ = 298.0;
    in->buffer_tank_.Activate("Feed", kPressureBoundaryEndpointName);
    in->buffer_tank_.percentage_ = 60.0;
    purge.top_boundary_condition_ = FlowInlet(BoundaryLocation::TOP, std::move(in));
    purge.bottom_boundary_condition_ = PressureOutlet(BoundaryLocation::BOTTOM);
  }
  data.cycle.AddStep("Feed");
  {
    Step& feed = data.cycle.GetStep("Feed");
    feed.step_time_ = 20.0;
    auto air = std::make_unique<BoundaryStream>();
    air->name_ = "Feed In";
    air->boundary_flow_ = Velocity(0.1);
    air->mole_fraction["Oxygen"] = 0.21;
    air->mole_fraction["Nitrogen"] = 0.79;
    air->temperature_ = 298.0;
    feed.bottom_boundary_condition_ = FlowInlet(BoundaryLocation::BOTTOM, std::move(air));
    feed.top_boundary_condition_ = PressureOutlet(BoundaryLocation::TOP);
  }
  data.cycle.EnsureIds();
  data.cycle.number_of_cycles_ = 3;
  data.settings.max_step = 1.0;
  return data;
}

inline Simulation BufferTankSimulation() {
  Simulation sim;
  sim.SetData(BufferTankDeck());
  return sim;
}

inline void AddLangmuirAdsorbent(UserData& data, AdsorbentID id, double scale) {
  auto& ads = data.adsorbent_directory.Get(id);
  ads.ep_ = 0.37;
  auto langmuir = std::make_unique<LangmuirIsotherm>();
  langmuir->params_["Oxygen"].m1 = 0.5 * scale;
  langmuir->params_["Oxygen"].b1 = 3e-11;
  langmuir->params_["Oxygen"].q1 = 30000;
  langmuir->params_["Nitrogen"].m1 = 1.5 * scale;
  langmuir->params_["Nitrogen"].b1 = 3e-11;
  langmuir->params_["Nitrogen"].q1 = 35000;
  ads.adsorption_source_.AddIsotherm({"Oxygen", "Nitrogen"}, std::move(langmuir));
  for (const char* comp : {"Oxygen", "Nitrogen"}) {
    auto ldf = std::make_unique<LinearDrivingForce>();
    auto mtc = std::make_unique<MassTransferCoefficientConstant>();
    mtc->constant_mtc_ = 0.1;
    ldf->mtc_model_ = std::move(mtc);
    ads.adsorption_source_.mass_transfer_models_[comp] = std::move(ldf);
  }
}

// Two layers (cells1, cells2): "13X" (adsorbent A) below a mixed layer of A
// and B (B only there), an optional top void, and a buffer-tank purge fed from
// the feed vent (consumer before source).
inline UserData TwoLayerDeck(int cells1, int cells2, bool top_void = true) {
  UserData data;
  data.bed.wall_.thickness_ = 0.0175;
  data.bed.wall_.heat_capacity_ = 537;
  data.bed.wall_.density_ = 7800;
  data.bed.wall_.thermal_conductivity_ = 16.2;
  data.fluid.AddComponent("Oxygen");
  data.fluid.AddComponent("Nitrogen");

  const auto a = data.adsorbent_directory.Add("13X");
  const auto b = data.adsorbent_directory.Add("5A");
  AddLangmuirAdsorbent(data, a, 1.0);
  AddLangmuirAdsorbent(data, b, 0.6);

  auto make_layer = [&](const std::string& name, double length, int cells, std::map<AdsorbentID, double> comp) {
    auto id = data.bed.layers_.Add(name);
    auto& layer = data.bed.layers_.Get(id);
    layer.momentum_model_ = MomentumModel::KarmenCozeny;
    layer.domain_.length = length;
    layer.domain_.n_cells = cells;
    layer.particle_diameter_ = 0.0011;
    layer.eb_ = 0.37;
    layer.composition_ = std::move(comp);
    auto ihtc = std::make_unique<ConstantInternalHeatTransferModel>();
    ihtc->heat_transfer_coefficient_ = 0.0;
    layer.equations_.internal_heat_transfer_model_ = std::move(ihtc);
    auto ehtc = std::make_unique<ConstantExternalHeatTransferModel>();
    ehtc->heat_transfer_coefficient_ = 0.0;
    layer.equations_.external_heat_transfer_model_ = std::move(ehtc);
  };
  make_layer("Bottom", 0.6, cells1, {{a, 1.0}});
  make_layer("Top", 0.4, cells2, {{a, 0.5}, {b, 0.5}});
  if (top_void) {
    data.bed.top_void_.enabled_ = true;
    data.bed.top_void_.volume_ = 2e-4;
  }

  data.ic.pressure = 101325.0;
  data.ic.temperature = 298.0;
  data.ic.wall_temperature = 298.0;
  data.ic.velocity = 0.0;
  data.ic.mole_fraction["Oxygen"] = 0.21;
  data.ic.mole_fraction["Nitrogen"] = 0.79;

  data.cycle.AddStep("Purge");
  {
    Step& purge = data.cycle.GetStep("Purge");
    purge.step_time_ = 20.0;
    auto in = std::make_unique<BoundaryStream>();
    in->name_ = "Purge In";
    in->boundary_flow_ = Velocity(0.01);
    in->mole_fraction["Oxygen"] = 1.0;
    in->mole_fraction["Nitrogen"] = 0.0;
    in->temperature_ = 298.0;
    in->buffer_tank_.Activate("Feed", kPressureBoundaryEndpointName);
    in->buffer_tank_.percentage_ = 60.0;
    purge.top_boundary_condition_ = FlowInlet(BoundaryLocation::TOP, std::move(in));
    purge.bottom_boundary_condition_ = PressureOutlet(BoundaryLocation::BOTTOM);
  }
  data.cycle.AddStep("Feed");
  {
    Step& feed = data.cycle.GetStep("Feed");
    feed.step_time_ = 20.0;
    auto air = std::make_unique<BoundaryStream>();
    air->name_ = "Feed In";
    air->boundary_flow_ = Velocity(0.1);
    air->mole_fraction["Oxygen"] = 0.21;
    air->mole_fraction["Nitrogen"] = 0.79;
    air->temperature_ = 298.0;
    feed.bottom_boundary_condition_ = FlowInlet(BoundaryLocation::BOTTOM, std::move(air));
    feed.top_boundary_condition_ = PressureOutlet(BoundaryLocation::TOP);
  }
  data.cycle.EnsureIds();
  data.cycle.number_of_cycles_ = 2;
  data.settings.max_step = 1.0;
  return data;
}

inline Simulation TwoLayerSimulation(int cells1, int cells2, bool top_void = true) {
  Simulation sim;
  sim.SetData(TwoLayerDeck(cells1, cells2, top_void));
  return sim;
}

}  // namespace adsim::test_decks
