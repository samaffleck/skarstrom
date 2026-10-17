/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "adsim/objects/deck_builder.h"

#include <algorithm>
#include <cctype>
#include <stdexcept>

#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/boundary_value.h"
#include "adsim/boundary_conditions/connection.h"
#include "adsim/boundary_conditions/flow_boundary_value.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/equations/isotherm.h"
#include "adsim/equations/isotherm_kinds.h"
#include "adsim/equations/mass_transfer_model.h"
#include "adsim/equations/wall_temperature_boundary_condition.h"
#include "adsim/objects/compressor.h"
#include "adsim/objects/deck_edit.h"
#include "adsim/objects/wall.h"

namespace adsim::build {

namespace {

// Index of `kind` in `kinds`, or std::invalid_argument naming the valid kinds.
size_t KindIndex(const std::vector<std::string>& kinds, const std::string& kind, const char* family) {
  auto it = std::find(kinds.begin(), kinds.end(), kind);
  if (it != kinds.end()) return static_cast<size_t>(it - kinds.begin());
  std::string valid;
  for (const auto& k : kinds) valid += (valid.empty() ? "" : ", ") + k;
  throw std::invalid_argument("unknown " + std::string(family) + " kind '" + kind + "'; one of: " + valid);
}

template <class E>
E KindEnum(const std::vector<std::string>& kinds, const std::string& kind, const char* family) {
  return static_cast<E>(KindIndex(kinds, kind, family));
}

BoundaryLocation Location(const std::string& s) {
  if (s == "bottom") return BoundaryLocation::BOTTOM;
  if (s == "top") return BoundaryLocation::TOP;
  throw std::invalid_argument("unknown location '" + s + "'; one of: bottom, top");
}

BoundaryDirection Direction(const std::string& s) {
  if (s == "inlet") return BoundaryDirection::INLET;
  if (s == "outlet") return BoundaryDirection::OUTLET;
  throw std::invalid_argument("unknown direction '" + s + "'; one of: inlet, outlet");
}

Adsorbent& GetAdsorbent(UserData& data, AdsorbentID id) {
  if (!data.adsorbent_directory.Contains(id))
    throw std::invalid_argument("no adsorbent with ID " + std::to_string(id));
  return data.adsorbent_directory.Get(id);
}

MixedAdsorbentLayer& GetLayer(UserData& data, size_t id) {
  for (auto& layer : data.bed.layers_.layers_)
    if (layer.id_ == id) return layer;
  throw std::invalid_argument("no layer with ID " + std::to_string(id));
}

Step& GetStep(UserData& data, size_t id) {
  for (auto& step : data.cycle.step_list_)
    if (step.id_ == id) return step;
  throw std::invalid_argument("no step with ID " + std::to_string(id));
}

void CheckComponent(const UserData& data, const std::string& component) {
  const auto& c = data.fluid.components_;
  if (std::find(c.begin(), c.end(), component) == c.end())
    throw std::invalid_argument("'" + component + "' is not a component of this deck");
}

std::unique_ptr<IBoundaryValue> MakeProfile(const std::string& kind) {
  switch (KindIndex(ProfileKinds(), kind, "profile")) {
    case 0: return std::make_unique<ConstantBoundaryValue>();
    case 1: return std::make_unique<LinearBoundaryValue>();
    default: return std::make_unique<ExponentialBoundaryValue>();
  }
}

std::unique_ptr<IFlowBoundaryValue> MakeFlow(const std::string& kind, std::unique_ptr<IBoundaryValue> profile) {
  auto set = [&](auto flow) -> std::unique_ptr<IFlowBoundaryValue> {
    flow->profile_ = std::move(profile);
    return flow;
  };
  switch (KindIndex(FlowKinds(), kind, "flow")) {
    case 0: return set(std::make_unique<VelocityBoundaryValue>());
    case 1: return set(std::make_unique<VolumetricFlowBoundaryValue>());
    case 2: return set(std::make_unique<MolarFlowBoundaryValue>());
    default: return set(std::make_unique<MassFlowBoundaryValue>());
  }
}

// Replace one end of a step, dropping any connection the old boundary held.
void ReplaceEnd(UserData& data, Step& step, BoundaryLocation location, std::shared_ptr<IBoundaryCondition> bc) {
  bc->location_ = location;
  auto& slot = location == BoundaryLocation::TOP ? step.top_boundary_condition_ : step.bottom_boundary_condition_;
  if (slot) DisconnectBoundary(data, slot.get());
  slot = std::move(bc);
}

std::shared_ptr<IBoundaryCondition>& EndOf(Step& step, BoundaryLocation location) {
  return location == BoundaryLocation::TOP ? step.top_boundary_condition_ : step.bottom_boundary_condition_;
}

void SeedMoleFractions(const UserData& data, std::unordered_map<std::string, double>& y) {
  for (const auto& c : data.fluid.components_) y.emplace(c, 0.0);
}

FlowBoundaryCondition& FlowEnd(UserData& data, size_t step_id, const std::string& location) {
  Step& step = GetStep(data, step_id);
  auto* fbc = dynamic_cast<FlowBoundaryCondition*>(EndOf(step, Location(location)).get());
  if (!fbc)
    throw std::invalid_argument("the " + location + " of step '" + step.step_name_ +
                                "' is not a flow boundary; make it one first (SetFlowBoundary)");
  return *fbc;
}

size_t AddStreamTo(UserData& data, FlowBoundaryCondition& fbc, std::unique_ptr<IBoundaryStream> stream,
                   const std::string& name) {
  stream->name_ = name;
  stream->SetDirection(fbc.direction_);
  SeedMoleFractions(data, stream->mole_fraction);
  IBoundaryStream* raw = stream.get();
  fbc.streams_.push_back(std::move(stream));
  data.cycle.EnsureIds();
  return raw->endpoint_id_;
}

Cycle::EndpointHandle FindEndpoint(UserData& data, size_t endpoint) {
  auto h = data.cycle.FindEndpoint(endpoint);
  if (!h) throw std::invalid_argument("no stream or pressure boundary with endpoint ID " + std::to_string(endpoint));
  return *h;
}

}  // namespace

// ---- kind names -------------------------------------------------------------------

std::vector<std::string> IsothermKinds() {
  std::vector<std::string> out;
  for (auto t : AllIsothermTypes())
    if (IsSupportedIsothermKind(t)) out.emplace_back(IsothermKindString(t));
  return out;
}

const std::vector<std::string>& RateModelKinds() {
  static const std::vector<std::string> k{"ldf", "ldf2"};
  return k;
}
const std::vector<std::string>& RateCoefficientKinds() {
  static const std::vector<std::string> k{"constant", "macropore", "macropore_knudsen", "haynes", "arrhenius"};
  return k;
}
const std::vector<std::string>& MomentumKinds() {
  static const std::vector<std::string> k{"kozeny_carman", "ergun"};
  return k;
}
const std::vector<std::string>& DispersionKinds() {
  static const std::vector<std::string> k{"constant", "ruthven"};
  return k;
}
const std::vector<std::string>& HeatDispersionKinds() {
  static const std::vector<std::string> k{"constant", "ruthven", "yagi_kunii"};
  return k;
}
const std::vector<std::string>& InternalHeatTransferKinds() {
  static const std::vector<std::string> k{"constant", "dixon", "leva"};
  return k;
}
const std::vector<std::string>& ExternalHeatTransferKinds() {
  static const std::vector<std::string> k{"constant", "holman"};
  return k;
}
const std::vector<std::string>& ParticleHeatTransferKinds() {
  static const std::vector<std::string> k{"constant"};
  return k;
}
const std::vector<std::string>& ThermalEquilibriumKinds() {
  static const std::vector<std::string> k{"lte", "ltne"};
  return k;
}
const std::vector<std::string>& BedShapeKinds() {
  static const std::vector<std::string> k{"cylinder", "square", "rectangle"};
  return k;
}
const std::vector<std::string>& WallMaterialKinds() {
  static const std::vector<std::string> k{"custom",           "ss316l",         "ss304",       "ss304l",
                                          "carbon_steel_a106", "carbon_steel_a53", "hastelloy_c276",
                                          "inconel_625",       "aluminium_6061",  "titanium_gr2"};
  return k;
}
const std::vector<std::string>& WallBoundaryKinds() {
  static const std::vector<std::string> k{"zero_flux", "fixed"};
  return k;
}
const std::vector<std::string>& ProfileKinds() {
  static const std::vector<std::string> k{"constant", "linear", "exponential"};
  return k;
}
const std::vector<std::string>& FlowKinds() {
  static const std::vector<std::string> k{"velocity", "volumetric_flow", "molar_flow", "mass_flow"};
  return k;
}
const std::vector<std::string>& PumpKinds() {
  static const std::vector<std::string> k{"none", "compressor", "vacuum_pump"};
  return k;
}
const std::vector<std::string>& PumpEnergyKinds() {
  static const std::vector<std::string> k{"isothermal", "adiabatic"};
  return k;
}
const std::vector<std::string>& PumpEfficiencyKinds() {
  static const std::vector<std::string> k{"constant", "vacuum_pump", "polynomial"};
  return k;
}

namespace {
// UnitCategory values are not contiguous (Custom = 8 sits in the middle).
const std::vector<std::pair<std::string, UnitCategory>>& VariableCategories() {
  using U = UnitCategory;
  static const std::vector<std::pair<std::string, UnitCategory>> k{
      {"dimensionless", U::Dimensionless},
      {"time", U::Time},
      {"length", U::Length},
      {"temperature", U::Temperature},
      {"pressure", U::Pressure},
      {"energy", U::Energy},
      {"specific_molar_energy", U::SpecificMolarEnergy},
      {"inverse_time", U::InverseTime},
      {"velocity", U::Velocity},
      {"volumetric_flow_rate", U::VolumetricFlowRate},
      {"molar_flow_rate", U::MolarFlowRate},
      {"mass_flow_rate", U::MassFlowRate},
      {"mole_fraction", U::MoleFraction},
      {"viscosity", U::Viscosity},
      {"density", U::Density},
      {"molar_concentration", U::MolarConcentration},
      {"heat_capacity", U::HeatCapacity},
      {"molar_mass", U::MolarMass},
      {"thermal_conductivity", U::ThermalConductivity},
      {"heat_transfer_coefficient", U::HeatTransferCoefficient},
      {"pressure_gradient", U::PressureGradient},
      {"adsorbed_concentration", U::AdsorbedConcentration},
      {"mass", U::Mass},
      {"volume", U::Volume},
      {"molar_amount", U::MolarAmount},
      {"area", U::Area},
      {"diffusivity", U::Diffusivity},
      {"specific_energy", U::SpecificEnergy},
      {"custom", U::Custom},
  };
  return k;
}
}  // namespace

const std::vector<std::string>& VariableCategoryKinds() {
  static const std::vector<std::string> k = [] {
    std::vector<std::string> out;
    for (const auto& [name, category] : VariableCategories()) out.push_back(name);
    return out;
  }();
  return k;
}

std::string VariableCategoryKind(UnitCategory category) {
  for (const auto& [name, c] : VariableCategories())
    if (c == category) return name;
  return "custom";
}

// ---- components and adsorbents --------------------------------------------------------

void AddComponent(UserData& data, const std::string& name) {
  std::string err;
  if (!AddComponentToDeck(data, name, err)) throw std::invalid_argument(err);
  // 0 is what Fluid::Initialise assumes for a component without a value.
  data.fluid.adsorbed_heat_capacity_.emplace(name, 0.0);
  data.kpi_config.SyncToComponents(data.fluid.components_);
}

AdsorbentID AddAdsorbent(UserData& data, const std::string& name) {
  // No layer is seeded here: a layer gets its adsorbents from SetLayerComposition.
  const AdsorbentID id = data.adsorbent_directory.Add(name);
  data.adsorbent_directory.Get(id).adsorption_source_.Initialise(data.fluid);
  return id;
}

void SetIsotherm(UserData& data, AdsorbentID adsorbent, const std::string& kind,
                 const std::vector<std::string>& components) {
  IsothermType type{};
  if (!KindFromString(kind, type) || !IsSupportedIsothermKind(type))
    KindIndex(IsothermKinds(), kind, "isotherm");  // throws with the list
  if (components.empty()) throw std::invalid_argument("SetIsotherm needs at least one component");
  for (const auto& c : components) CheckComponent(data, c);
  GetAdsorbent(data, adsorbent).adsorption_source_.AssignIsotherm(components, type);
}

void SetRateModel(UserData& data, AdsorbentID adsorbent, const std::string& component, const std::string& model,
                  const std::string& coefficient) {
  CheckComponent(data, component);
  auto m = CreateMassTransferModel(KindEnum<MassTransferModelType>(RateModelKinds(), model, "rate model"));
  m->mtc_model_ = CreateMassTransferCoefficientModel(
      KindEnum<MassTransferCoefficientType>(RateCoefficientKinds(), coefficient, "rate coefficient"));
  GetAdsorbent(data, adsorbent).adsorption_source_.mass_transfer_models_[component] = std::move(m);
}

// ---- column -------------------------------------------------------------------------------

void SetBedShape(UserData& data, const std::string& kind) {
  switch (KindIndex(BedShapeKinds(), kind, "bed shape")) {
    case 0: data.bed.shape_ = CylinderBed{}; break;
    case 1: data.bed.shape_ = SquareBed{}; break;
    default: data.bed.shape_ = RectangularBed{}; break;
  }
}

void SetWallMaterial(UserData& data, const std::string& kind) {
  data.bed.wall_.ApplyMaterial(KindEnum<WallMaterial>(WallMaterialKinds(), kind, "wall material"));
}

void SetWallBoundary(UserData& data, const std::string& kind) {
  data.bed.equations_.wall_temperature_boundary_condition_ = CreateWallTemperatureBoundaryCondition(
      KindEnum<WallTemperatureBoundaryConditionType>(WallBoundaryKinds(), kind, "wall boundary"));
}

void SetThermalEquilibrium(UserData& data, const std::string& kind) {
  data.bed.equations_.thermal_equilibrium_mode_ =
      KindEnum<ThermalEquilibriumMode>(ThermalEquilibriumKinds(), kind, "thermal equilibrium");
}

void SetParticleHeatTransfer(UserData& data, const std::string& kind) {
  data.bed.equations_.particle_heat_transfer_model_ = CreateParticleHeatTransferModel(
      KindEnum<ParticleHeatTransferType>(ParticleHeatTransferKinds(), kind, "particle heat transfer"));
}

size_t AddLayer(UserData& data, const std::string& name) {
  const LayerID id = data.bed.layers_.Add(name);
  data.bed.layers_.Get(id).composition_.clear();
  return id;
}

void SetLayerComposition(UserData& data, size_t layer,
                         const std::vector<std::pair<AdsorbentID, double>>& composition) {
  auto& l = GetLayer(data, layer);
  for (const auto& [id, fraction] : composition) {
    (void)fraction;
    GetAdsorbent(data, id);  // throws for an unknown ID
  }
  l.composition_.clear();
  for (const auto& [id, fraction] : composition) l.composition_[id] = fraction;
}

void SetMomentum(UserData& data, size_t layer, const std::string& kind) {
  GetLayer(data, layer).momentum_model_ = KindEnum<MomentumModel>(MomentumKinds(), kind, "momentum");
}

void SetDispersion(UserData& data, size_t layer, const std::string& kind) {
  GetLayer(data, layer).equations_.dispersion_model_ =
      CreateAxialDispersionModel(KindEnum<AxialDispersionType>(DispersionKinds(), kind, "dispersion"));
}

void SetHeatDispersion(UserData& data, size_t layer, const std::string& kind) {
  GetLayer(data, layer).equations_.heat_dispersion_model_ = CreateAxialHeatDispersionModel(
      KindEnum<AxialHeatDispersionType>(HeatDispersionKinds(), kind, "heat dispersion"));
}

void SetInternalHeatTransfer(UserData& data, size_t layer, const std::string& kind) {
  GetLayer(data, layer).equations_.internal_heat_transfer_model_ = CreateInternalHeatTransferModel(
      KindEnum<InternalHeatTransferType>(InternalHeatTransferKinds(), kind, "internal heat transfer"));
}

void SetExternalHeatTransfer(UserData& data, size_t layer, const std::string& kind) {
  GetLayer(data, layer).equations_.external_heat_transfer_model_ = CreateExternalHeatTransferModel(
      KindEnum<ExternalHeatTransferType>(ExternalHeatTransferKinds(), kind, "external heat transfer"));
}

// ---- cycle ----------------------------------------------------------------------------------

size_t AddStep(UserData& data, const std::string& name) {
  for (const auto& s : data.cycle.step_list_)
    if (s.step_name_ == name) throw std::invalid_argument("the cycle already has a step named '" + name + "'");
  Step step;  // comes with a default pressure top and a one-stream bottom: replace both
  step.step_name_ = name;
  for (auto location : {BoundaryLocation::BOTTOM, BoundaryLocation::TOP}) {
    auto closed = std::make_shared<ClosedBoundaryCondition>();
    closed->location_ = location;
    EndOf(step, location) = std::move(closed);
  }
  data.cycle.step_list_.push_back(std::move(step));
  data.cycle.EnsureIds();
  return data.cycle.step_list_.back().id_;
}

void SetClosed(UserData& data, size_t step, const std::string& location) {
  Step& s = GetStep(data, step);
  ReplaceEnd(data, s, Location(location), std::make_shared<ClosedBoundaryCondition>());
}

size_t SetPressureBoundary(UserData& data, size_t step, const std::string& location, const std::string& direction,
                           const std::string& profile) {
  Step& s = GetStep(data, step);
  auto bc = std::make_shared<PressureBoundaryCondition>();
  bc->direction_ = Direction(direction);
  bc->pressure_ = MakeProfile(profile);
  SeedMoleFractions(data, bc->mole_fraction_);
  PressureBoundaryCondition* raw = bc.get();
  ReplaceEnd(data, s, Location(location), std::move(bc));
  data.cycle.EnsureIds();
  return raw->endpoint_id_;
}

void SetFlowBoundary(UserData& data, size_t step, const std::string& location, const std::string& direction) {
  Step& s = GetStep(data, step);
  auto bc = std::make_shared<FlowBoundaryCondition>();
  bc->direction_ = Direction(direction);
  ReplaceEnd(data, s, Location(location), std::move(bc));
}

size_t AddStream(UserData& data, size_t step, const std::string& location, const std::string& name,
                 const std::string& flow, const std::string& profile) {
  auto& fbc = FlowEnd(data, step, location);
  auto stream = std::make_unique<BoundaryStream>();
  stream->boundary_flow_ = MakeFlow(flow, MakeProfile(profile));
  return AddStreamTo(data, fbc, std::move(stream), name);
}

size_t AddValveStream(UserData& data, size_t step, const std::string& location, const std::string& name,
                      const std::string& profile) {
  auto& fbc = FlowEnd(data, step, location);
  auto stream = std::make_unique<ValveBoundaryStream>();
  stream->boundary_value_ = MakeProfile(profile);
  return AddStreamTo(data, fbc, std::move(stream), name);
}

void SetTemperatureProfile(UserData& data, size_t endpoint, const std::string& profile) {
  auto h = FindEndpoint(data, endpoint);
  std::unique_ptr<IBoundaryValue> value;
  if (KindIndex(ProfileKinds(), profile, "profile") != 0) value = MakeProfile(profile);
  if (h.stream) h.stream->temperature_profile_ = std::move(value);
  else h.pressure->temperature_profile_ = std::move(value);
}

void SetPump(UserData& data, size_t endpoint, const std::string& kind, const std::string& energy,
             const std::string& efficiency) {
  auto h = FindEndpoint(data, endpoint);
  std::unique_ptr<IPump>* slot = nullptr;
  if (h.pressure) {
    slot = &h.pressure->pump_;
  } else if (auto* s = dynamic_cast<BoundaryStream*>(h.stream)) {
    slot = &s->pump_;
  } else {
    throw std::invalid_argument("endpoint " + std::to_string(endpoint) + " ('" + h.log_name +
                                "') is a valve stream, which has no pump");
  }
  const size_t k = KindIndex(PumpKinds(), kind, "pump");
  if (k == 0) {
    if (*slot) (*slot)->use_pump_ = false;
    return;
  }
  auto pump = CreatePump(static_cast<PumpType>(k - 1));
  pump->use_pump_ = true;
  pump->energy_eq_ = CreateCompressorEnergyModel(KindEnum<CompressorEnergyType>(PumpEnergyKinds(), energy, "pump energy"));
  pump->efficiency_model_ =
      CreatePumpEfficiencyModel(KindEnum<PumpEfficiencyType>(PumpEfficiencyKinds(), efficiency, "pump efficiency"));
  *slot = std::move(pump);
}

void Connect(UserData& data, size_t endpoint_a, size_t endpoint_b) {
  auto a = FindEndpoint(data, endpoint_a);
  auto b = FindEndpoint(data, endpoint_b);
  auto make = [](Cycle::EndpointHandle& h) {
    return h.stream ? MakeEndpoint(*h.step, *h.stream) : MakeEndpoint(*h.step, *h.pressure);
  };
  ConnectionEndpoint ea = make(a), eb = make(b);
  if (!MakeMutualConnection(ea, eb))
    throw std::invalid_argument("cannot connect endpoints " + std::to_string(endpoint_a) + " and " +
                                std::to_string(endpoint_b) +
                                ": one is already connected, or both are pressure boundaries");
}

void SetEnergySource(UserData& data, size_t step, size_t layer, bool enabled) {
  GetLayer(data, layer);
  GetStep(data, step).energy_sources_[layer].enabled_ = enabled;
}

void SetKPI(UserData& data, const std::string& component, const std::vector<size_t>& feed_endpoints,
            const std::vector<size_t>& product_endpoints, const std::vector<size_t>& recycle_endpoints) {
  CheckComponent(data, component);
  data.kpi_config.SyncToComponents(data.fluid.components_);
  ComponentKPIConfig* kpi = data.kpi_config.FindByComponent(component);
  if (!kpi) throw std::invalid_argument("no KPI entry for '" + component + "'");
  auto refs = [&](const std::vector<size_t>& ids) {
    std::vector<KPIEndpointRef> out;
    for (size_t id : ids) {
      KPIEndpointRef ref;
      ref.endpoint_id = id;
      if (!ref.Resolve(data.cycle.step_list_))
        throw std::invalid_argument("no stream or pressure boundary with endpoint ID " + std::to_string(id));
      out.push_back(ref);
    }
    return out;
  };
  kpi->tracked = true;
  kpi->feed_endpoints = refs(feed_endpoints);
  kpi->product_endpoints = refs(product_endpoints);
  kpi->recycle_endpoints = refs(recycle_endpoints);
}

// ---- user variables -----------------------------------------------------------------

void AddVariable(UserData& data, const std::string& name, const std::string& category,
                 const std::string& formula, const std::string& description) {
  const bool identifier = !name.empty() && (std::isalpha(static_cast<unsigned char>(name[0])) || name[0] == '_') &&
                          std::all_of(name.begin(), name.end(), [](char c) {
                            return std::isalnum(static_cast<unsigned char>(c)) || c == '_';
                          });
  if (!identifier) throw std::invalid_argument("a variable name is letters, digits and _ ('" + name + "')");
  if (UserVariableRegistry::IsReservedName(name))
    throw std::invalid_argument("'" + name + "' is a built-in name (T, P, y, R, PI...); choose another");
  UserVariable var;
  var.name = name;
  var.unit_category = VariableCategories()[KindIndex(VariableCategoryKinds(), category, "variable category")].second;
  var.value_expr = formula;
  var.description = description;
  if (const UserVariable* old = data.variable_registry.Find(name)) var.value = old->value;
  data.variable_registry.Add(std::move(var));
}

}  // namespace adsim::build
