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
#include <fstream>
#include <sstream>
#include <unordered_map>

#include <cereal/archives/json.hpp>

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
    if (id != keep) { others += w; ++n; }
  if (n == 0) return;
  for (auto& [id, w] : layer.composition_) {
    if (id == keep) continue;
    w = others > 0.0 ? w * (target / others) : target / static_cast<double>(n);
  }
}

bool ValidFraction(double w, bool allow_zero, std::string& err) {
  if (w <= 1.0 && (allow_zero ? w >= 0.0 : w > 0.0)) return true;
  std::ostringstream ss;
  ss << "a mass fraction must be " << (allow_zero ? "between 0 and 1" : "above 0 and at most 1") << ", not " << w;
  err = ss.str();
  return false;
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
    frac << it->second;
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

}  // namespace adsim
