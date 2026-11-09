/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "adsim/objects/deck_edit.h"

#include <algorithm>
#include <cassert>
#include <filesystem>
#include <map>
#include <optional>
#include <fstream>
#include <sstream>
#include <unordered_map>

#include <cereal/archives/json.hpp>

#include "adsim/equations/isotherm_kinds.h"
#include "adsim/objects/adsorbent_library.h"

namespace adsim {

bool IsAvailableComponent(const std::string& name) {
  for (const auto& c : Components::componentNames)
    if (c == name) return true;
  return false;
}

void ReconcileMoleFractions(UserData& data, const std::string& component, bool add) {
  auto touch = [&](ExpressionMap& mf) {
    if (add) {
      if (!mf.count(component)) mf[component] = 0.0;
    } else {
      mf.erase(component);
    }
  };
  touch(data.ic.mole_fraction);
  for (auto& step : data.cycle.step_list_) {
    for (IBoundaryCondition* bc : {step.top_boundary_condition_.get(), step.bottom_boundary_condition_.get()}) {
      if (!bc) continue;
      if (auto* p = dynamic_cast<PressureBoundaryCondition*>(bc)) {
        touch(p->mole_fraction_);
      } else if (auto* f = dynamic_cast<FlowBoundaryCondition*>(bc)) {
        for (auto& s : f->streams_)
          if (s) touch(s->mole_fraction);
      }
    }
  }
}

bool AddComponentToDeck(UserData& data, const std::string& name, std::string& err) {
  if (!IsAvailableComponent(name)) {
    err = "Unknown component: " + name;
    return false;
  }
  data.fluid.AddComponent(name);  // dedupes internally
  ReconcileMoleFractions(data, name, true);
  // 0 is what Fluid::Initialise assumes for a component without a value; seeding
  // it here makes the entry editable before the first run (a front end lists
  // the map's entries, so a component without one had no field).
  data.fluid.adsorbed_heat_capacity_.emplace(name, 0.0);
  data.fluid.gas_heat_capacity_.emplace(name, Fluid::kDefaultGasHeatCapacity);
  // Initialise is idempotent: existing components keep their models while the
  // new one gets a full, editable default set.
  for (auto& [id, ads_const] : data.adsorbent_directory.GetAll()) {
    auto& ads = const_cast<Adsorbent&>(ads_const);
    ads.adsorption_source_.Initialise(data.fluid);
  }
  return true;
}

std::string ComponentInUseByReactions(const UserData& data, const std::string& name) {
  std::vector<AdsorbentID> ids;
  for (const auto& [id, ads] : data.adsorbent_directory.GetAll()) ids.push_back(id);
  std::sort(ids.begin(), ids.end());
  std::string where;
  size_t count = 0;
  for (AdsorbentID id : ids) {
    const Adsorbent& ads = data.adsorbent_directory.Get(id);
    const auto names = ads.reaction_source_.ReactionsUsing(name);
    if (names.empty()) continue;
    if (!where.empty()) where += "; ";
    for (size_t i = 0; i < names.size(); ++i) where += (i ? ", " : "") + names[i];
    where += " (on " + ads.name_ + ")";
    count += names.size();
  }
  if (count == 0) return "";
  return name + " is used by " + (count == 1 ? "reaction " : "reactions ") + where + "; remove " +
         (count == 1 ? "it" : "them") + " first";
}

bool RemoveComponentFromDeck(UserData& data, const std::string& name, std::string& err) {
  const std::string in_use = ComponentInUseByReactions(data, name);
  if (!in_use.empty()) {
    err = in_use;
    return false;
  }
  for (auto& [id, ads_const] : data.adsorbent_directory.GetAll())
    const_cast<Adsorbent&>(ads_const).reaction_source_.RemoveComponent(name);
  data.fluid.RemoveComponent(name);
  data.bed.RemoveComponent(name, data.adsorbent_directory);  // strips it from every adsorbent
  ReconcileMoleFractions(data, name, false);
  data.fluid.adsorbed_heat_capacity_.erase(name);
  data.fluid.gas_heat_capacity_.erase(name);
  data.kpi_config.reaction.RemoveComponent(name);
  return true;
}

void RemoveComponentFromDeck(UserData& data, const std::string& name) {
  std::string err;
  [[maybe_unused]] const bool removed = RemoveComponentFromDeck(data, name, err);
  assert(removed && "RemoveComponentFromDeck: a reaction uses the component; use the checked form");
}

void SeedLayerComposition(UserData& data, MixedAdsorbentLayer& layer) {
  const auto& dir = data.adsorbent_directory.GetAll();
  if (dir.empty() || !layer.composition_.empty()) return;
  AdsorbentID lowest = dir.begin()->first;
  for (const auto& [id, a] : dir) lowest = std::min(lowest, id);
  layer.AddAdsorbent(lowest, 1.0);
}

namespace {

// Scale every fraction except `keep`'s so they sum to `target` (in proportion;
// equal shares when they are all 0).
void RescaleOthers(MixedAdsorbentLayer& layer, AdsorbentID keep, double target) {
  double others = 0.0;
  size_t n = 0;
  for (const auto& [id, w] : layer.composition_)
    if (id != keep) { others += w.Get(); ++n; }
  if (n == 0) return;
  for (auto& [id, w] : layer.composition_) {
    if (id == keep) continue;
    w = others > 0.0 ? w.Get() * (target / others) : target / static_cast<double>(n);
  }
}

bool ValidFraction(double w, bool allow_zero, std::string& err) {
  if (w <= 1.0 && (allow_zero ? w >= 0.0 : w > 0.0)) return true;
  std::ostringstream ss;
  ss << "a mass fraction must be " << (allow_zero ? "between 0 and 1" : "above 0 and at most 1") << ", not " << w;
  err = ss.str();
  return false;
}

// The refusal of the rescaling edits for a layer whose fractions include a
// formula (deck_edit.h).
bool RefuseFormulaBlend(const MixedAdsorbentLayer& layer, std::string& err) {
  if (!layer.HasCompositionExpression()) return false;
  std::string terms;
  for (const auto& [id, w] : layer.composition_) terms += (terms.empty() ? "" : ", ") + w.GetExpr();
  err = "the mass fractions of layer \"" + layer.name_ + "\" include expressions (" + terms +
        "), which cannot be rescaled; set each fraction directly, or make them numbers first";
  return true;
}

}  // namespace

bool AddSolidToLayer(UserData& data, MixedAdsorbentLayer& layer, AdsorbentID id, double mass_fraction,
                     std::string& err) {
  if (!data.adsorbent_directory.Contains(id)) {
    err = "no adsorbent with id " + std::to_string(id);
    return false;
  }
  if (layer.composition_.count(id)) return SetLayerMassFraction(data, layer, id, mass_fraction, err);
  if (layer.composition_.empty()) {
    layer.composition_[id] = 1.0;
    return true;
  }
  if (RefuseFormulaBlend(layer, err)) return false;
  if (!ValidFraction(mass_fraction, false, err)) return false;
  layer.composition_[id] = mass_fraction;
  RescaleOthers(layer, id, 1.0 - mass_fraction);
  return true;
}

bool SetLayerMassFraction(UserData& data, MixedAdsorbentLayer& layer, AdsorbentID id, double mass_fraction,
                          std::string& err) {
  (void)data;
  auto it = layer.composition_.find(id);
  if (it == layer.composition_.end()) {
    err = "layer \"" + layer.name_ + "\" does not hold adsorbent id " + std::to_string(id);
    return false;
  }
  if (RefuseFormulaBlend(layer, err)) return false;
  if (!ValidFraction(mass_fraction, true, err)) return false;
  if (layer.composition_.size() == 1 && mass_fraction != 1.0) {
    err = "the only solid of layer \"" + layer.name_ + "\" has mass fraction 1; add another solid first";
    return false;
  }
  it->second = mass_fraction;
  RescaleOthers(layer, id, 1.0 - mass_fraction);
  return true;
}

bool RemoveSolidFromLayer(UserData& data, MixedAdsorbentLayer& layer, AdsorbentID id, std::string& err) {
  (void)data;
  if (!layer.composition_.count(id)) {
    err = "layer \"" + layer.name_ + "\" does not hold adsorbent id " + std::to_string(id);
    return false;
  }
  if (layer.composition_.size() == 1) {
    err = "it is the only solid of layer \"" + layer.name_ + "\"; add another solid first";
    return false;
  }
  if (RefuseFormulaBlend(layer, err)) return false;
  layer.composition_.erase(id);
  RescaleOthers(layer, id, 1.0);
  return true;
}

void SeedEmptyLayerCompositions(UserData& data) {
  for (auto& layer : data.bed.layers_.layers_)
    if (layer.composition_.empty()) SeedLayerComposition(data, layer);
}

AdsorbentID AddAdsorbentToDeck(UserData& data, const std::string& name) {
  const AdsorbentID id = data.adsorbent_directory.Add(name);
  data.adsorbent_directory.Get(id).adsorption_source_.Initialise(data.fluid);
  SeedEmptyLayerCompositions(data);
  return id;
}

bool LoadAdsorbentFromString(const std::string& json, Adsorbent& out, std::string& err) {
  try {
    std::istringstream iss(json);
    cereal::JSONInputArchive archive(iss);
    archive(out);
    return true;
  } catch (const std::exception& e) {
    err = std::string("failed to parse adsorbent: ") + e.what();
    return false;
  }
}

bool AddAdsorbentFromJson(UserData& data, const std::string& label, const std::string& json,
                          AdsorbentID& id, std::string& err) {
  id = data.adsorbent_directory.Add(label);
  auto& ads = data.adsorbent_directory.Get(id);
  if (!LoadAdsorbentFromString(json, ads, err)) {
    data.adsorbent_directory.Remove(id);
    return false;
  }
  ads.adsorption_source_.Initialise(data.fluid);
  ads.reaction_source_.Initialise(data.fluid);
  SeedEmptyLayerCompositions(data);
  return true;
}

bool AddLibraryAdsorbent(UserData& data, const std::string& library_dir, const std::string& key,
                         AdsorbentID& id, std::string& err) {
  namespace fs = std::filesystem;
  try {
    const AdsorbentLibrary lib = LoadAdsorbentLibraryFile((fs::path(library_dir) / "library.json").string());
    const AdsorbentLibraryEntry* e = lib.Resolve(key);
    if (!e) {
      err = "no adsorbent or material '" + key + "' in the library";
      return false;
    }
    std::ifstream in(fs::path(library_dir) / e->file, std::ios::binary);
    if (!in) {
      err = "cannot open " + (fs::path(library_dir) / e->file).string();
      return false;
    }
    std::ostringstream ss;
    ss << in.rdbuf();
    return AddAdsorbentFromJson(data, e->label, ss.str(), id, err);
  } catch (const std::exception& ex) {
    err = ex.what();
    return false;
  }
}

std::string AdsorbentInUseByMixtures(const UserData& data, AdsorbentID id) {
  if (!data.adsorbent_directory.Contains(id)) return "";
  std::string where;
  size_t count = 0;
  for (const auto& layer : data.bed.layers_.layers_) {
    auto it = layer.composition_.find(id);
    if (it == layer.composition_.end() || layer.composition_.size() < 2) continue;
    std::ostringstream frac;
    if (it->second.HasExpression()) frac << it->second.GetExpr() << " = ";
    frac << it->second.Get();
    where += (count ? "; " : "") + ("\"" + layer.name_ + "\" (mass fraction " + frac.str() + ")");
    ++count;
  }
  if (count == 0) return "";
  return data.adsorbent_directory.Get(id).name_ + " is part of the blend of " +
         (count == 1 ? "layer " : "layers ") + where + "; remove it from the " +
         (count == 1 ? "layer's composition" : "layers' compositions") + " first";
}

bool RemoveAdsorbentFromDeck(UserData& data, AdsorbentID id, std::string& err) {
  const std::string in_use = AdsorbentInUseByMixtures(data, id);
  if (!in_use.empty()) {
    err = in_use;
    return false;
  }
  data.RemoveAdsorbent(id);
  SeedEmptyLayerCompositions(data);
  return true;
}

void RemoveAdsorbentFromDeck(UserData& data, AdsorbentID id) {
  std::string err;
  [[maybe_unused]] const bool removed = RemoveAdsorbentFromDeck(data, id, err);
  assert(removed && "RemoveAdsorbentFromDeck: the solid is part of a blend; use the checked form");
}

void DisconnectPair(UserData& data, Connection& conn) {
  if (!conn.IsConnected()) return;
  const size_t pid = conn.GetID();
  ForEachEndpoint(data, [&](ConnectionEndpoint e) {
    if (e.conn && e.conn->IsConnected() && e.conn->GetID() == pid) e.conn->Disconnect();
  });
}

void DisconnectBoundary(UserData& data, IBoundaryCondition* bc) {
  if (!bc) return;
  if (auto* pbc = dynamic_cast<PressureBoundaryCondition*>(bc)) {
    DisconnectPair(data, pbc->connection_);
  } else if (auto* fbc = dynamic_cast<FlowBoundaryCondition*>(bc)) {
    for (auto& s : fbc->streams_)
      if (s) DisconnectPair(data, s->connection_);
  }
}


// ---- connections, identity and repair ------------------------------------------------

namespace {

Connection& ConnectionOf(const Cycle::EndpointHandle& h) {
  return h.stream ? h.stream->connection_ : h.pressure->connection_;
}

BufferTankSource& TankOf(const Cycle::EndpointHandle& h) {
  return h.stream ? h.stream->buffer_tank_ : h.pressure->buffer_tank_;
}

// Point `h`'s connection at `partner`'s stream data under pair id `id`.
void LinkTo(const Cycle::EndpointHandle& h, const Cycle::EndpointHandle& partner, size_t id) {
  if (partner.stream) ConnectionOf(h).MakeConnection(id, *partner.stream);
  else ConnectionOf(h).MakeConnection(id, *partner.pressure);
}

// Re-link every connection from its pair id (RepairConnections); with
// `renumber`, the surviving pairs get ids 1..N in order of appearance.
size_t RelinkConnections(UserData& data, bool renumber) {
  data.cycle.EnsureIds();
  std::vector<Cycle::EndpointHandle> eps = data.cycle.Endpoints();
  std::vector<size_t> order;                         // pair ids by first appearance
  std::map<size_t, std::vector<size_t>> groups;      // pair id -> endpoint indices
  for (size_t i = 0; i < eps.size(); ++i) {
    const Connection& c = ConnectionOf(eps[i]);
    if (!c.IsConnected()) continue;
    auto& g = groups[c.GetID()];
    if (g.empty()) order.push_back(c.GetID());
    g.push_back(i);
  }
  size_t disconnected = 0, next = 0;
  for (size_t id : order) {
    const auto& g = groups[id];
    const bool pair = g.size() == 2 && (eps[g[0]].stream || eps[g[1]].stream);
    if (!pair) {
      for (size_t i : g) ConnectionOf(eps[i]).Disconnect();
      disconnected += g.size();
      continue;
    }
    const size_t new_id = renumber ? ++next : id;
    LinkTo(eps[g[0]], eps[g[1]], new_id);
    LinkTo(eps[g[1]], eps[g[0]], new_id);
  }
  return disconnected;
}

// Steps, endpoints and layers get IDs (no pair-id renumbering).
void EnsureObjectIds(UserData& data) {
  data.cycle.EnsureIds();
  auto& layers = data.bed.layers_;
  LayerID next = 0;
  for (const auto& l : layers.layers_) next = std::max(next, l.id_);
  std::vector<LayerID> seen;
  for (auto& l : layers.layers_) {
    if (l.id_ == 0 || std::find(seen.begin(), seen.end(), l.id_) != seen.end()) l.id_ = ++next;
    seen.push_back(l.id_);
  }
  layers.EnsureNextIdAtLeast(next);
}

std::shared_ptr<IBoundaryCondition>& EndOf(Step& step, BoundaryLocation side) {
  return side == BoundaryLocation::TOP ? step.top_boundary_condition_ : step.bottom_boundary_condition_;
}

const char* SideName(BoundaryLocation side) { return side == BoundaryLocation::TOP ? "top" : "bottom"; }

Step* StepOrError(UserData& data, size_t step_id, std::string& err) {
  Step* s = data.cycle.FindStepById(step_id);
  if (!s) err = "no step with ID " + std::to_string(step_id);
  return s;
}

MixedAdsorbentLayer* LayerOrError(UserData& data, LayerID id, std::string& err) {
  for (auto& l : data.bed.layers_.layers_)
    if (l.id_ == id) return &l;
  err = "no layer with ID " + std::to_string(id);
  return nullptr;
}

Adsorbent* AdsorbentOrError(UserData& data, AdsorbentID id, std::string& err) {
  if (data.adsorbent_directory.Contains(id)) return &data.adsorbent_directory.Get(id);
  err = "no adsorbent with ID " + std::to_string(id);
  return nullptr;
}

bool HasComponent(const UserData& data, const std::string& c, std::string& err) {
  const auto& comps = data.fluid.components_;
  if (std::find(comps.begin(), comps.end(), c) != comps.end()) return true;
  err = "'" + c + "' is not a component of this deck";
  return false;
}

FlowBoundaryCondition* FlowEndOrError(UserData& data, size_t step_id, BoundaryLocation side, std::string& err) {
  Step* step = StepOrError(data, step_id, err);
  if (!step) return nullptr;
  auto* fbc = dynamic_cast<FlowBoundaryCondition*>(EndOf(*step, side).get());
  if (!fbc) err = std::string("the ") + SideName(side) + " of step '" + step->step_name_ + "' is not a flow boundary";
  return fbc;
}

std::optional<Cycle::EndpointHandle> EndpointOrError(UserData& data, size_t endpoint_id, std::string& err) {
  data.cycle.EnsureIds();
  auto h = data.cycle.FindEndpoint(endpoint_id);
  if (!h) err = "no stream or pressure boundary with endpoint ID " + std::to_string(endpoint_id);
  return h;
}

// After an endpoint, a stream or a whole end was removed or replaced: new
// IDs, no KPI reference or buffer tank left pointing at nothing.
void AfterEndpointsRemoved(UserData& data) {
  data.cycle.EnsureIds();
  data.kpi_config.Prune(data.cycle.step_list_);
  PruneBufferTanks(data);
}

// The value `y` carries for an endpoint end: the map and the inlet
// temperature (and its profile) of a pressure boundary or a stream.
struct InletState {
  ExpressionMap y;
  ExpressionParam temperature{293.0};
  std::unique_ptr<IBoundaryValue> temperature_profile;
  bool set = false;
};

}  // namespace

bool ConnectEndpoints(UserData& data, size_t endpoint_a, size_t endpoint_b, std::string& err) {
  auto a = EndpointOrError(data, endpoint_a, err);
  if (!a) return false;
  auto b = EndpointOrError(data, endpoint_b, err);
  if (!b) return false;
  if (endpoint_a == endpoint_b) {
    err = "an endpoint cannot be connected to itself";
    return false;
  }
  if (a->pressure && b->pressure) {
    err = "two pressure boundaries cannot be connected; one end must be a stream";
    return false;
  }
  for (const auto* h : {&*a, &*b}) {
    if (ConnectionOf(*h).IsConnected()) {
      err = "endpoint " + std::to_string(h->id) + " ('" + h->log_name + "') is already connected";
      return false;
    }
    if (TankOf(*h).IsActive()) {
      err = "endpoint " + std::to_string(h->id) + " ('" + h->log_name + "') draws from a buffer tank";
      return false;
    }
  }
  size_t id = 0;
  for (const auto& h : data.cycle.Endpoints()) id = std::max(id, ConnectionOf(h).GetID());
  ++id;
  LinkTo(*a, *b, id);
  LinkTo(*b, *a, id);
  return true;
}

void EnsureIds(UserData& data) {
  EnsureObjectIds(data);
  RelinkConnections(data, true);
}

size_t RepairConnections(UserData& data) { return RelinkConnections(data, false); }

size_t PruneBufferTanks(UserData& data) {
  data.cycle.ResolveBufferTankReferences();  // also EnsureIds; stores IDs for name-only sources
  const auto eps = data.cycle.Endpoints();
  size_t n = 0;
  for (const auto& h : eps) {
    BufferTankSource& tank = TankOf(h);
    if (!tank.IsActive()) continue;
    const bool exists = tank.source_endpoint_id_ != 0 &&
                        std::any_of(eps.begin(), eps.end(), [&](const auto& e) { return e.id == tank.source_endpoint_id_; });
    if (exists) continue;
    tank.Deactivate();
    ++n;
  }
  return n;
}

void SyncKpiConfigToComponents(UserData& data) { data.kpi_config.SyncToComponents(data.fluid.components_); }

void PruneKpiConfig(UserData& data) {
  data.cycle.EnsureIds();
  data.kpi_config.Prune(data.cycle.step_list_);
}

void FillDeckGaps(UserData& data) {
  const auto& comps = data.fluid.components_;
  ReconcileMoleFractionMap(data.ic.mole_fraction, comps);
  for (const auto& h : data.cycle.Endpoints())
    ReconcileMoleFractionMap(h.stream ? h.stream->mole_fraction : h.pressure->mole_fraction_, comps);
  for (const auto& c : comps) {
    data.fluid.adsorbed_heat_capacity_.emplace(c, 0.0);
    data.fluid.gas_heat_capacity_.emplace(c, Fluid::kDefaultGasHeatCapacity);
  }
  for (auto& [id, ads_const] : data.adsorbent_directory.GetAll()) {
    auto& ads = const_cast<Adsorbent&>(ads_const);
    ads.adsorption_source_.Initialise(data.fluid);
    ads.reaction_source_.Initialise(data.fluid);
  }
  SeedEmptyLayerCompositions(data);
  SyncKpiConfigToComponents(data);
  EnsureObjectIds(data);
  RepairConnections(data);
}

// ---- layers ---------------------------------------------------------------------------------

std::string UniqueName(const std::vector<std::string>& taken, const std::string& base) {
  auto used = [&](const std::string& n) { return std::find(taken.begin(), taken.end(), n) != taken.end(); };
  if (!used(base)) return base;
  for (int n = 2;; ++n) {
    const std::string candidate = base + " " + std::to_string(n);
    if (!used(candidate)) return candidate;
  }
}

LayerID AddLayer(UserData& data, const std::string& name) {
  EnsureObjectIds(data);
  auto& layers = data.bed.layers_;
  std::vector<std::string> taken;
  for (const auto& l : layers.layers_) taken.push_back(l.name_);
  std::string base = name;
  if (base.empty()) {
    // "Layer N", the next free number from the new layer count.
    for (size_t n = layers.layers_.size() + 1;; ++n) {
      base = "Layer " + std::to_string(n);
      if (std::find(taken.begin(), taken.end(), base) == taken.end()) break;
    }
  }
  const LayerID id = layers.Add(UniqueName(taken, base));
  SeedLayerComposition(data, layers.Get(id));
  return id;
}

bool RemoveLayer(UserData& data, LayerID id, std::string& err) {
  if (!LayerOrError(data, id, err)) return false;
  data.bed.layers_.Remove(id);
  for (auto& step : data.cycle.step_list_) step.energy_sources_.erase(id);
  return true;
}

bool MoveLayer(UserData& data, LayerID id, int offset, std::string& err) {
  if (!LayerOrError(data, id, err)) return false;
  if (offset != 1 && offset != -1) {
    err = "a layer moves one place at a time (offset +1 or -1)";
    return false;
  }
  auto& v = data.bed.layers_.layers_;
  const auto it = std::find_if(v.begin(), v.end(), [id](const auto& l) { return l.id_ == id; });
  const auto i = static_cast<std::ptrdiff_t>(it - v.begin());
  const auto j = i + offset;
  if (j < 0 || j >= static_cast<std::ptrdiff_t>(v.size())) {
    err = "layer \"" + it->name_ + "\" is already at the " + (offset > 0 ? "product" : "feed") + " end of the bed";
    return false;
  }
  std::swap(v[static_cast<size_t>(i)], v[static_cast<size_t>(j)]);
  return true;
}

// ---- steps ------------------------------------------------------------------------------------

void SeedMoleFractions(const UserData& data, ExpressionMap& y) {
  const auto& comps = data.fluid.components_;
  ReconcileMoleFractionMap(y, comps);
  double sum = 0.0;
  for (const auto& [c, v] : y) sum += v.Get();
  if (sum == 0.0 && !comps.empty()) y[comps.front()] = 1.0;
}

namespace {

void SeedEndpointsOf(const UserData& data, Step& step) {
  for (IBoundaryCondition* bc : {step.bottom_boundary_condition_.get(), step.top_boundary_condition_.get()}) {
    if (auto* p = dynamic_cast<PressureBoundaryCondition*>(bc)) {
      SeedMoleFractions(data, p->mole_fraction_);
    } else if (auto* f = dynamic_cast<FlowBoundaryCondition*>(bc)) {
      for (auto& s : f->streams_)
        if (s) SeedMoleFractions(data, s->mole_fraction);
    }
  }
}

}  // namespace

size_t AddStep(UserData& data, const std::string& name) {
  auto& steps = data.cycle.step_list_;
  std::vector<std::string> taken;
  for (const auto& s : steps) taken.push_back(s.step_name_);
  std::string base = name;
  if (base.empty()) {
    for (size_t n = steps.size() + 1;; ++n) {
      base = "Step " + std::to_string(n);
      if (std::find(taken.begin(), taken.end(), base) == taken.end()) break;
    }
  }
  Step step;  // a flow inlet with "Stream 1" at the bottom, a constant-pressure outlet at the top
  step.step_name_ = UniqueName(taken, base);
  SeedEndpointsOf(data, step);
  steps.push_back(std::move(step));
  data.cycle.EnsureIds();
  return steps.back().id_;
}

bool RemoveStep(UserData& data, size_t step_id, std::string& err) {
  const int index = data.cycle.StepIndexById(step_id);
  if (index < 0) {
    err = "no step with ID " + std::to_string(step_id);
    return false;
  }
  data.RemoveStep(static_cast<size_t>(index));
  AfterEndpointsRemoved(data);
  return true;
}

bool DuplicateStep(UserData& data, size_t step_id, size_t& new_step_id, std::string& err) {
  const int index = data.cycle.StepIndexById(step_id);
  if (index < 0) {
    err = "no step with ID " + std::to_string(step_id);
    return false;
  }
  const size_t at = data.DuplicateStep(static_cast<size_t>(index));
  data.cycle.EnsureIds();
  new_step_id = data.cycle.step_list_[at].id_;
  return true;
}

bool MoveStep(UserData& data, size_t step_id, int offset, std::string& err) {
  const int index = data.cycle.StepIndexById(step_id);
  if (index < 0) {
    err = "no step with ID " + std::to_string(step_id);
    return false;
  }
  if (offset != 1 && offset != -1) {
    err = "a step moves one place at a time (offset +1 or -1)";
    return false;
  }
  const int to = index + offset;
  if (to < 0 || to >= static_cast<int>(data.cycle.step_list_.size())) {
    err = "step '" + data.cycle.step_list_[static_cast<size_t>(index)].step_name_ + "' is already the " +
          (offset > 0 ? "last" : "first") + " step";
    return false;
  }
  data.SwapSteps(static_cast<size_t>(index), static_cast<size_t>(to));
  return true;
}

// ---- streams ------------------------------------------------------------------------------------

namespace {

std::unique_ptr<IBoundaryStream> MakeStream(BoundaryStreamType kind, BoundaryDirection direction) {
  std::unique_ptr<IBoundaryStream> s = CreateStream(kind);
  if (!s) return s;
  s->SetDirection(direction);
  if (auto* b = dynamic_cast<BoundaryStream*>(s.get())) EnsurePumpMatchesDirection(b->pump_, direction);
  if (auto* p = dynamic_cast<PumpBoundaryStream*>(s.get())) {
    if (EnsurePumpMatchesDirection(p->pump_, direction)) p->pump_->use_pump_ = true;
  }
  return s;
}

}  // namespace

bool AddStream(UserData& data, size_t step_id, BoundaryLocation side, BoundaryStreamType kind, size_t& endpoint_id,
               std::string& err) {
  FlowBoundaryCondition* fbc = FlowEndOrError(data, step_id, side, err);
  if (!fbc) return false;
  auto stream = MakeStream(kind, fbc->direction_);
  if (!stream) {
    err = "unknown stream kind " + std::to_string(static_cast<int>(kind));
    return false;
  }
  std::vector<std::string> taken;
  for (const auto& s : fbc->streams_)
    if (s) taken.push_back(s->name_);
  std::string base;
  for (size_t n = fbc->streams_.size() + 1;; ++n) {
    base = "Stream " + std::to_string(n);
    if (std::find(taken.begin(), taken.end(), base) == taken.end()) break;
  }
  stream->name_ = base;
  SeedMoleFractions(data, stream->mole_fraction);
  IBoundaryStream* raw = stream.get();
  fbc->streams_.push_back(std::move(stream));
  data.cycle.EnsureIds();
  endpoint_id = raw->endpoint_id_;
  return true;
}

bool RemoveStream(UserData& data, size_t step_id, BoundaryLocation side, size_t index, std::string& err) {
  FlowBoundaryCondition* fbc = FlowEndOrError(data, step_id, side, err);
  if (!fbc) return false;
  if (index >= fbc->streams_.size()) {
    err = std::string("the ") + SideName(side) + " of step ID " + std::to_string(step_id) + " has no stream " +
          std::to_string(index);
    return false;
  }
  data.cycle.RemoveStream(*fbc, index);
  AfterEndpointsRemoved(data);
  return true;
}

// ---- kind switches -------------------------------------------------------------------------------

namespace {

// Replace a model held in a unique_ptr by `fresh`, carrying over the shared
// parameters (both bound through ModelFields).
template <class Ptr, class New>
void SwitchModel(Ptr& slot, New fresh, ModelFamily family, std::size_t old_kind, std::size_t new_kind) {
  const KindDescriptor* from = KindAt(family, old_kind);
  const KindDescriptor* to = KindAt(family, new_kind);
  if (slot && fresh && from && to) CarryOverParameters(*from, ModelFields(slot.get()), *to, ModelFields(fresh.get()));
  slot = std::move(fresh);
}

bool BadKind(std::string& err, const char* family, int value) {
  err = std::string("unknown ") + family + " kind " + std::to_string(value);
  return false;
}

// The rate-coefficient descriptor as the given rate model reads it (the LDF2
// rate constant has another unit).
KindDescriptor CoefficientKind(MassTransferCoefficientType kind, MassTransferModelType rate) {
  KindDescriptor d = *KindAt(ModelFamily::RateCoefficient, static_cast<std::size_t>(kind));
  if (rate == MassTransferModelType::LDF2)
    for (auto& p : d.parameters) p = Ldf2RateConstant(p);
  return d;
}

}  // namespace

bool SetIsothermKind(UserData& data, AdsorbentID adsorbent, const std::string& component, IsothermType kind,
                     std::string& err) {
  Adsorbent* ads = AdsorbentOrError(data, adsorbent, err);
  if (!ads || !HasComponent(data, component, err)) return false;
  if (static_cast<std::size_t>(kind) >= IsothermKindDescriptors().size() || !IsSupportedIsothermKind(kind))
    return BadKind(err, "isotherm", static_cast<int>(kind));
  auto& src = ads->adsorption_source_;
  auto it = src.IsothermMap().find(component);
  if (it != src.IsothermMap().end() && it->second == kind) return true;

  // The component's current parameters, kept by value: the old isotherm drops
  // the component when it switches.
  std::vector<std::pair<std::string, ExpressionParam>> saved;
  IsothermType old_kind = IsothermType::INERT;
  if (it != src.IsothermMap().end()) {
    old_kind = it->second;
    for (const auto& [name, p] : IsothermFields(src.GetIsotherm(component).get(), component))
      saved.emplace_back(name, *p);
  }
  src.AssignIsotherm({component}, kind);

  ModelFieldList from;
  for (auto& [name, p] : saved) from.emplace_back(name, &p);
  CarryOverParameters(IsothermKind(old_kind), from, IsothermKind(kind),
                      IsothermFields(src.GetIsotherm(component).get(), component));
  return true;
}

namespace {

// Applies a coefficient preset to one site's coefficient model, as a kind
// switch did before the generic model: the parameters the new preset shows
// (under its kind's names) take its defaults, except those the old preset
// shared with it (same name and unit, the rate model's unit for k); the
// parameters it does not show keep their stored values. `fresh` is the new
// preset's model with its defaults. series_resistance (custom) keeps a generic
// site as it is (switches, mode and values), except that the lumped k takes
// its default when the rate model changes order (its unit changes).
void ApplyCoefficientPreset(std::unique_ptr<IMassTransferCoefficientModel>& site,
                            std::unique_ptr<IMassTransferCoefficientModel> fresh, MassTransferModelType old_rate,
                            MassTransferCoefficientType coefficient_kind, MassTransferModelType rate_kind) {
  using T = MassTransferCoefficientType;
  const bool ldf2_old = old_rate == MassTransferModelType::LDF2, ldf2_new = rate_kind == MassTransferModelType::LDF2;
  auto* g_new = dynamic_cast<MassTransferCoefficientSeriesResistance*>(fresh.get());
  if (coefficient_kind == T::SERIES_RESISTANCE && site) {
    if (auto* g_old = dynamic_cast<MassTransferCoefficientSeriesResistance*>(site.get())) {
      if (ldf2_old != ldf2_new && g_new) g_old->k_ = g_new->k_;
      return;
    }
  }
  if (site && g_new) {
    if (auto* g_old = dynamic_cast<MassTransferCoefficientSeriesResistance*>(site.get())) {
      // Hidden parameters: the stored values. Shown ones and those the preset
      // fixes: the preset's.
      std::vector<const ExpressionParam*> fixed;
      for (const auto& [name, p] : PresetFields(fresh.get(), coefficient_kind)) fixed.push_back(p);
      if (coefficient_kind == T::CONSTANT) fixed.push_back(&g_new->activation_energy_);
      if (coefficient_kind == T::HAYNES_1973) {
        fixed.push_back(&g_new->film_resistance_factor_);
        fixed.push_back(&g_new->molecular_diffusion_coefficient_);
      }
      const ModelFieldList from_all = ResistanceFields(g_old);
      const ModelFieldList to_all = ResistanceFields(g_new);
      for (std::size_t i = 0; i < to_all.size() && i < from_all.size(); ++i)
        if (std::find(fixed.begin(), fixed.end(), to_all[i].second) == fixed.end()) *to_all[i].second = *from_all[i].second;
    }
  }
  if (site) {
    const T old_preset = CoefficientPreset(site.get());
    const KindDescriptor from = CoefficientKind(old_preset, old_rate);
    const KindDescriptor to = CoefficientKind(coefficient_kind, rate_kind);
    CarryOverParameters(from, PresetFields(site.get(), old_preset), to, PresetFields(fresh.get(), coefficient_kind));
  }
  site = std::move(fresh);
}

MassTransferCoefficientSeriesResistance* ResistanceSite(UserData& data, AdsorbentID adsorbent,
                                                        const std::string& component, MassTransferSite site,
                                                        std::string& err) {
  Adsorbent* ads = AdsorbentOrError(data, adsorbent, err);
  if (!ads || !HasComponent(data, component, err)) return nullptr;
  auto& models = ads->adsorption_source_.mass_transfer_models_;
  auto it = models.find(component);
  if (it == models.end() || !it->second) {
    err = "component '" + component + "' has no mass-transfer model";
    return nullptr;
  }
  IMassTransferModel& rate = *it->second;
  std::unique_ptr<IMassTransferCoefficientModel>* slot = &rate.mtc_model_;
  if (site == MassTransferSite::Slow) {
    auto* dual = dynamic_cast<DualRateLinearDrivingForce*>(&rate);
    if (!dual) {
      err = "the mass-transfer model of '" + component + "' has no slow site (only dual_rate_ldf has one)";
      return nullptr;
    }
    slot = &dual->slow_mtc_model_;
  }
  if (!*slot) *slot = CreateMassTransferCoefficientModel(MassTransferCoefficientType::CONSTANT);
  *slot = ToResistances(std::move(*slot));
  auto* g = dynamic_cast<MassTransferCoefficientSeriesResistance*>(slot->get());
  if (!g) err = "the coefficient model of '" + component + "' is not a resistance model";
  return g;
}

}  // namespace

bool SetMassTransferKind(UserData& data, AdsorbentID adsorbent, const std::string& component,
                         MassTransferModelType rate_kind, MassTransferCoefficientType coefficient_kind,
                         std::string& err) {
  Adsorbent* ads = AdsorbentOrError(data, adsorbent, err);
  if (!ads || !HasComponent(data, component, err)) return false;
  if (!KindAt(ModelFamily::RateModel, static_cast<std::size_t>(rate_kind)))
    return BadKind(err, "rate model", static_cast<int>(rate_kind));
  if (!KindAt(ModelFamily::RateCoefficient, static_cast<std::size_t>(coefficient_kind)))
    return BadKind(err, "rate coefficient", static_cast<int>(coefficient_kind));
  auto& slot = ads->adsorption_source_.mass_transfer_models_[component];
  if (!slot) slot = std::make_unique<LinearDrivingForce>();
  const MassTransferModelType old_rate = slot->GetModelType();
  std::unique_ptr<IMassTransferCoefficientModel> old_coeff = std::move(slot->mtc_model_);
  std::unique_ptr<IMassTransferCoefficientModel> old_slow;
  if (auto* dual = dynamic_cast<DualRateLinearDrivingForce*>(slot.get())) old_slow = std::move(dual->slow_mtc_model_);
  if (old_rate != rate_kind) slot = CreateMassTransferModel(rate_kind);
  auto* dual = dynamic_cast<DualRateLinearDrivingForce*>(slot.get());

  if (old_coeff && CoefficientPreset(old_coeff.get()) == coefficient_kind && old_rate == rate_kind) {
    slot->mtc_model_ = std::move(old_coeff);
    if (dual) dual->slow_mtc_model_ = std::move(old_slow);
    if (dual) dual->SyncSlowSite();
    return true;
  }
  // The coefficient kind is a preset of the generic model (mass_transfer_model.h).
  slot->mtc_model_ = std::move(old_coeff);
  ApplyCoefficientPreset(slot->mtc_model_, CreateMassTransferCoefficientModel(coefficient_kind), old_rate,
                         coefficient_kind, rate_kind);
  // A dual-rate model's slow site takes the preset too: a new dual-rate model
  // starts with the slow site a copy of the fast one (so it behaves as the LDF
  // until split or the slow coefficient changes); an existing one keeps the
  // slow-site values the old and new presets share, the rest from the fast
  // site. After that the two sites' switches are independent.
  if (dual) {
    dual->slow_mtc_model_.reset();
    if (old_rate == rate_kind && old_slow) {
      dual->SyncSlowSite();  // a copy of the new fast site
      auto fresh = std::move(dual->slow_mtc_model_);
      dual->slow_mtc_model_ = std::move(old_slow);
      ApplyCoefficientPreset(dual->slow_mtc_model_, std::move(fresh), old_rate, coefficient_kind, rate_kind);
    } else {
      dual->SyncSlowSite();
    }
  }
  return true;
}

bool SetMassTransferResistance(UserData& data, AdsorbentID adsorbent, const std::string& component,
                               MassTransferSite site, const std::string& tag, bool on, std::string& err) {
  MassTransferCoefficientSeriesResistance* g = ResistanceSite(data, adsorbent, component, site, err);
  return g && g->SetResistance(tag, on, err);
}

bool SetCapacityFactorMode(UserData& data, AdsorbentID adsorbent, const std::string& component, MassTransferSite site,
                           CapacityFactorMode mode, std::string& err) {
  MassTransferCoefficientSeriesResistance* g = ResistanceSite(data, adsorbent, component, site, err);
  if (!g) return false;
  if (static_cast<int>(mode) < 0 || mode > CapacityFactorMode::FIXED)
    return BadKind(err, "capacity-factor mode", static_cast<int>(mode));
  g->SetCapacityMode(mode);
  return true;
}

bool SetDispersionKind(UserData& data, LayerID layer, AxialDispersionType kind, std::string& err) {
  MixedAdsorbentLayer* l = LayerOrError(data, layer, err);
  if (!l) return false;
  auto fresh = CreateAxialDispersionModel(kind);
  if (!fresh) return BadKind(err, "dispersion", static_cast<int>(kind));
  auto& slot = l->equations_.dispersion_model_;
  if (slot && slot->GetModelType() == kind) return true;
  SwitchModel(slot, std::move(fresh), ModelFamily::Dispersion,
              slot ? static_cast<std::size_t>(slot->GetModelType()) : 0, static_cast<std::size_t>(kind));
  return true;
}

bool SetHeatDispersionKind(UserData& data, LayerID layer, AxialHeatDispersionType kind, std::string& err) {
  MixedAdsorbentLayer* l = LayerOrError(data, layer, err);
  if (!l) return false;
  auto fresh = CreateAxialHeatDispersionModel(kind);
  if (!fresh) return BadKind(err, "heat dispersion", static_cast<int>(kind));
  auto& slot = l->equations_.heat_dispersion_model_;
  if (slot && slot->GetModelType() == kind) return true;
  SwitchModel(slot, std::move(fresh), ModelFamily::HeatDispersion,
              slot ? static_cast<std::size_t>(slot->GetModelType()) : 0, static_cast<std::size_t>(kind));
  return true;
}

bool SetInternalHeatTransferKind(UserData& data, LayerID layer, InternalHeatTransferType kind, std::string& err) {
  MixedAdsorbentLayer* l = LayerOrError(data, layer, err);
  if (!l) return false;
  auto fresh = CreateInternalHeatTransferModel(kind);
  if (!fresh) return BadKind(err, "internal heat transfer", static_cast<int>(kind));
  auto& slot = l->equations_.internal_heat_transfer_model_;
  if (slot && slot->GetModelType() == kind) return true;
  SwitchModel(slot, std::move(fresh), ModelFamily::InternalHeatTransfer,
              slot ? static_cast<std::size_t>(slot->GetModelType()) : 0, static_cast<std::size_t>(kind));
  return true;
}

bool SetExternalHeatTransferKind(UserData& data, LayerID layer, ExternalHeatTransferType kind, std::string& err) {
  MixedAdsorbentLayer* l = LayerOrError(data, layer, err);
  if (!l) return false;
  auto fresh = CreateExternalHeatTransferModel(kind);
  if (!fresh) return BadKind(err, "external heat transfer", static_cast<int>(kind));
  auto& slot = l->equations_.external_heat_transfer_model_;
  if (slot && slot->GetModelType() == kind) return true;
  SwitchModel(slot, std::move(fresh), ModelFamily::ExternalHeatTransfer,
              slot ? static_cast<std::size_t>(slot->GetModelType()) : 0, static_cast<std::size_t>(kind));
  return true;
}

bool SetParticleHeatTransferKind(UserData& data, ParticleHeatTransferType kind, std::string& err) {
  auto fresh = CreateParticleHeatTransferModel(kind);
  if (!fresh) return BadKind(err, "particle heat transfer", static_cast<int>(kind));
  auto& slot = data.bed.equations_.particle_heat_transfer_model_;
  if (slot && slot->GetModelType() == kind) return true;
  SwitchModel(slot, std::move(fresh), ModelFamily::ParticleHeatTransfer,
              slot ? static_cast<std::size_t>(slot->GetModelType()) : 0, static_cast<std::size_t>(kind));
  return true;
}

namespace {

// The inlet state of an end (a pressure boundary, or the first stream of a
// flow boundary), moved out of it.
InletState TakeInlet(IBoundaryCondition* bc) {
  InletState s;
  if (auto* p = dynamic_cast<PressureBoundaryCondition*>(bc)) {
    s.y = p->mole_fraction_;
    s.temperature = p->temperature_;
    s.temperature_profile = std::move(p->temperature_profile_);
    s.set = true;
  } else if (auto* f = dynamic_cast<FlowBoundaryCondition*>(bc)) {
    if (!f->streams_.empty() && f->streams_.front()) {
      auto& st = *f->streams_.front();
      s.y = st.mole_fraction;
      s.temperature = st.temperature_;
      s.temperature_profile = std::move(st.temperature_profile_);
      s.set = true;
    }
  }
  return s;
}

}  // namespace

bool SetBoundaryConditionKind(UserData& data, size_t step_id, BoundaryLocation side, BoundaryType kind,
                              std::string& err) {
  Step* step = StepOrError(data, step_id, err);
  if (!step) return false;
  std::shared_ptr<IBoundaryCondition> fresh = CreateBoundaryCondition(kind);
  if (!fresh) return BadKind(err, "boundary", static_cast<int>(kind));
  auto& slot = EndOf(*step, side);
  if (slot && slot->GetType() == kind) return true;

  const BoundaryDirection direction = slot ? slot->direction_ : BoundaryDirection::OUTLET;
  if (slot) DisconnectBoundary(data, slot.get());
  InletState inlet = TakeInlet(slot.get());
  fresh->location_ = side;
  if (kind != BoundaryType::CLOSED) fresh->direction_ = direction;
  if (auto* p = dynamic_cast<PressureBoundaryCondition*>(fresh.get())) {
    if (inlet.set) {
      p->mole_fraction_ = std::move(inlet.y);
      p->temperature_ = inlet.temperature;
      p->temperature_profile_ = std::move(inlet.temperature_profile);
    }
    SeedMoleFractions(data, p->mole_fraction_);
  } else if (auto* f = dynamic_cast<FlowBoundaryCondition*>(fresh.get())) {
    auto stream = MakeStream(BoundaryStreamType::NONE, f->direction_);
    stream->name_ = "Stream 1";
    if (inlet.set) {
      stream->mole_fraction = std::move(inlet.y);
      stream->temperature_ = inlet.temperature;
      stream->temperature_profile_ = std::move(inlet.temperature_profile);
    }
    SeedMoleFractions(data, stream->mole_fraction);
    f->streams_.push_back(std::move(stream));
  }
  slot = std::move(fresh);
  AfterEndpointsRemoved(data);
  return true;
}

bool SetStreamKind(UserData& data, size_t step_id, BoundaryLocation side, size_t index, BoundaryStreamType kind,
                   std::string& err) {
  FlowBoundaryCondition* fbc = FlowEndOrError(data, step_id, side, err);
  if (!fbc) return false;
  if (index >= fbc->streams_.size() || !fbc->streams_[index]) {
    err = std::string("the ") + SideName(side) + " of step ID " + std::to_string(step_id) + " has no stream " +
          std::to_string(index);
    return false;
  }
  auto& slot = fbc->streams_[index];
  if (slot->GetType() == kind) return true;
  auto fresh = MakeStream(kind, fbc->direction_);
  if (!fresh) return BadKind(err, "stream", static_cast<int>(kind));
  fresh->name_ = slot->name_;
  fresh->temperature_ = slot->temperature_;
  fresh->temperature_profile_ = std::move(slot->temperature_profile_);
  fresh->mole_fraction = slot->mole_fraction;
  fresh->buffer_tank_ = slot->buffer_tank_;
  // The pump moves between a set-flow stream and a pump stream.
  std::unique_ptr<IPump>* old_pump = nullptr;
  if (auto* b = dynamic_cast<BoundaryStream*>(slot.get())) old_pump = &b->pump_;
  if (auto* p = dynamic_cast<PumpBoundaryStream*>(slot.get())) old_pump = &p->pump_;
  if (old_pump && *old_pump) {
    if (auto* b = dynamic_cast<BoundaryStream*>(fresh.get())) b->pump_ = std::move(*old_pump);
    if (auto* p = dynamic_cast<PumpBoundaryStream*>(fresh.get())) {
      p->pump_ = std::move(*old_pump);
      p->pump_->use_pump_ = true;  // a pump stream is driven by its pump
    }
  }
  data.cycle.ReplaceStream(slot, std::move(fresh));  // keeps the endpoint ID, drops the connection
  AfterEndpointsRemoved(data);
  return true;
}

namespace {

// Carry a profile's parameters to a fresh profile: by name and unit, a
// constant's value <-> a ramp's end, and use_step_time between ramps.
// `values_too` false keeps only the parameters that are not the profiled
// quantity (the flow kind changed, so its values mean something else).
void CarryProfile(IBoundaryValue* from, IBoundaryValue* to, bool values_too) {
  if (!from || !to) return;
  const KindDescriptor* fk = KindAt(ModelFamily::Profile, static_cast<std::size_t>(from->GetType()));
  const KindDescriptor* tk = KindAt(ModelFamily::Profile, static_cast<std::size_t>(to->GetType()));
  if (!fk || !tk) return;
  ModelFieldList f = ModelFields(from), t = ModelFields(to);
  auto rename = [](ModelFieldList& l, const char* a, const char* b) {
    for (auto& [n, p] : l)
      if (n == a) n = b;
  };
  // A constant's value is the value a ramp settles at: its end.
  if (from->GetType() == BoundaryValueType::CONSTANT) rename(f, "value", "end");
  if (to->GetType() == BoundaryValueType::CONSTANT) rename(t, "value", "end");
  KindDescriptor fd = *fk, td = *tk;
  for (auto* d : {&fd, &td})
    for (auto& p : d->parameters) {
      if (p.name == "value") p.name = "end";
      if (!values_too && p.unit == kProfileQuantityUnit) p.name = "#" + p.name;  // never matches
    }
  CarryOverParameters(fd, f, td, t);
  auto use_step_time = [](IBoundaryValue* v) -> bool* {
    if (auto* x = dynamic_cast<LinearBoundaryValue*>(v)) return &x->use_step_time_;
    if (auto* x = dynamic_cast<ExponentialBoundaryValue*>(v)) return &x->use_step_time_;
    if (auto* x = dynamic_cast<LogLinearBoundaryValue*>(v)) return &x->use_step_time_;
    return nullptr;
  };
  if (bool* a = use_step_time(from))
    if (bool* b = use_step_time(to)) *b = *a;
}

}  // namespace

bool SetFlowKind(UserData& data, size_t endpoint_id, FlowType kind, std::string& err) {
  auto h = EndpointOrError(data, endpoint_id, err);
  if (!h) return false;
  auto* s = dynamic_cast<BoundaryStream*>(h->stream);
  if (!s) {
    err = "endpoint " + std::to_string(endpoint_id) + " ('" + h->log_name + "') is not a set-flow stream";
    return false;
  }
  auto fresh = CreateFlowBoundaryValue(kind);
  if (!fresh) return BadKind(err, "flow", static_cast<int>(kind));
  if (s->boundary_flow_ && s->boundary_flow_->GetType() == kind) return true;
  IBoundaryValue* old_profile = s->boundary_flow_ ? s->boundary_flow_->profile_.get() : nullptr;
  if (old_profile && old_profile->GetType() != BoundaryValueType::CONSTANT) {
    auto ramp = CreateBoundaryValue(old_profile->GetType());
    CarryProfile(old_profile, ramp.get(), false);
    fresh->profile_ = std::move(ramp);
  }
  s->boundary_flow_ = std::move(fresh);
  return true;
}

bool SetProfileKind(UserData& data, size_t endpoint_id, ProfileTarget target, BoundaryValueType kind,
                    std::string& err) {
  auto h = EndpointOrError(data, endpoint_id, err);
  if (!h) return false;
  if (!KindAt(ModelFamily::Profile, static_cast<std::size_t>(kind))) return BadKind(err, "profile", static_cast<int>(kind));
  auto refuse = [&](const char* what) {
    err = "endpoint " + std::to_string(endpoint_id) + " ('" + h->log_name + "') has no " + what;
    return false;
  };
  if (target == ProfileTarget::Temperature) {
    ExpressionParam& constant = h->stream ? h->stream->temperature_ : h->pressure->temperature_;
    auto& slot = h->stream ? h->stream->temperature_profile_ : h->pressure->temperature_profile_;
    if (InletTemperatureProfileType(slot.get()) == kind) return true;
    if (kind == BoundaryValueType::CONSTANT) {
      slot.reset();  // the constant is the endpoint's temperature
      return true;
    }
    const Step& step = *h->step;
    auto fresh = CreateInletTemperatureProfile(kind, constant.Get(), step.step_time_.Get());
    if (!fresh) fresh = CreateBoundaryValue(kind);
    if (slot) CarryProfile(slot.get(), fresh.get(), true);
    slot = std::move(fresh);
    return true;
  }
  std::unique_ptr<IBoundaryValue>* slot = nullptr;
  if (target == ProfileTarget::Pressure) {
    if (h->pressure) slot = &h->pressure->pressure_;
    else if (auto* v = dynamic_cast<ValveBoundaryStream*>(h->stream)) slot = &v->boundary_value_;
    else return refuse("pressure profile (it is not a pressure boundary or a valve stream)");
  } else {
    auto* s = dynamic_cast<BoundaryStream*>(h->stream);
    if (!s || !s->boundary_flow_) return refuse("flow profile (it is not a set-flow stream)");
    slot = &s->boundary_flow_->profile_;
  }
  if (*slot && (*slot)->GetType() == kind) return true;
  auto fresh = CreateBoundaryValue(kind);
  CarryProfile(slot->get(), fresh.get(), true);
  *slot = std::move(fresh);
  return true;
}

namespace {

// The pump slot of an endpoint, or nullptr with err.
std::unique_ptr<IPump>* PumpSlot(const Cycle::EndpointHandle& h, size_t endpoint_id, std::string& err) {
  if (h.pressure) return &h.pressure->pump_;
  if (auto* b = dynamic_cast<BoundaryStream*>(h.stream)) return &b->pump_;
  if (auto* p = dynamic_cast<PumpBoundaryStream*>(h.stream)) return &p->pump_;
  err = "endpoint " + std::to_string(endpoint_id) + " ('" + h.log_name + "') is a valve stream, which has no pump";
  return nullptr;
}

}  // namespace

bool SetPumpKind(UserData& data, size_t endpoint_id, std::size_t pump_kind, std::string& err) {
  auto h = EndpointOrError(data, endpoint_id, err);
  if (!h) return false;
  if (!KindAt(ModelFamily::Pump, pump_kind)) return BadKind(err, "pump", static_cast<int>(pump_kind));
  std::unique_ptr<IPump>* slot = PumpSlot(*h, endpoint_id, err);
  if (!slot) return false;
  const bool pump_stream = dynamic_cast<PumpBoundaryStream*>(h->stream) != nullptr;
  if (pump_kind == 0) {
    if (pump_stream) {
      err = "endpoint " + std::to_string(endpoint_id) + " ('" + h->log_name +
            "') is a pump-driven stream, which always has a pump";
      return false;
    }
    if (*slot) (*slot)->use_pump_ = false;
    return true;
  }
  const auto type = static_cast<PumpType>(pump_kind - 1);
  if (*slot && (*slot)->GetType() == type) {
    (*slot)->use_pump_ = true;
    return true;
  }
  auto fresh = CreatePump(type);
  if (*slot) {
    IPump& old = **slot;
    fresh->poly_index_ = old.poly_index_;
    fresh->energy_eq_ = std::move(old.energy_eq_);
    fresh->efficiency_model_ = std::move(old.efficiency_model_);
    CarryOverParameters(*KindAt(ModelFamily::Pump, static_cast<std::size_t>(old.GetType()) + 1), ModelFields(&old),
                        *KindAt(ModelFamily::Pump, pump_kind), ModelFields(fresh.get()));
  }
  if (!fresh->energy_eq_) fresh->energy_eq_ = CreateCompressorEnergyModel(CompressorEnergyType::ISOTHERMAL);
  if (!fresh->efficiency_model_) fresh->efficiency_model_ = CreatePumpEfficiencyModel(PumpEfficiencyType::CONSTANT);
  fresh->use_pump_ = true;
  *slot = std::move(fresh);
  return true;
}

bool SetPumpEfficiencyKind(UserData& data, size_t endpoint_id, PumpEfficiencyType kind, std::string& err) {
  auto h = EndpointOrError(data, endpoint_id, err);
  if (!h) return false;
  std::unique_ptr<IPump>* slot = PumpSlot(*h, endpoint_id, err);
  if (!slot) return false;
  if (!*slot) {
    err = "endpoint " + std::to_string(endpoint_id) + " ('" + h->log_name + "') has no pump; set one first";
    return false;
  }
  auto fresh = CreatePumpEfficiencyModel(kind);
  if (!fresh) return BadKind(err, "pump efficiency", static_cast<int>(kind));
  auto& eff = (*slot)->efficiency_model_;
  if (eff && eff->GetType() == kind) return true;
  SwitchModel(eff, std::move(fresh), ModelFamily::PumpEfficiency,
              eff ? static_cast<std::size_t>(eff->GetType()) : 0, static_cast<std::size_t>(kind));
  return true;
}

bool SetPumpPowerKind(UserData& data, size_t endpoint_id, PumpPowerType kind, std::string& err) {
  auto h = EndpointOrError(data, endpoint_id, err);
  if (!h) return false;
  auto* p = dynamic_cast<PumpBoundaryStream*>(h->stream);
  if (!p) {
    err = "endpoint " + std::to_string(endpoint_id) + " ('" + h->log_name + "') is not a pump-driven stream";
    return false;
  }
  auto fresh = CreatePumpPowerModel(kind);
  if (!fresh) return BadKind(err, "pump power", static_cast<int>(kind));
  if (p->power_ && p->power_->GetType() == kind) return true;
  SwitchModel(p->power_, std::move(fresh), ModelFamily::PumpPower,
              p->power_ ? static_cast<std::size_t>(p->power_->GetType()) : 0, static_cast<std::size_t>(kind));
  return true;
}

}  // namespace adsim
