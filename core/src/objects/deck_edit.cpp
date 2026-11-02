/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "adsim/objects/deck_edit.h"

#include <algorithm>
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
  // Initialise is idempotent: existing components keep their models while the
  // new one gets a full, editable default set.
  for (auto& [id, ads_const] : data.adsorbent_directory.GetAll()) {
    auto& ads = const_cast<Adsorbent&>(ads_const);
    ads.adsorption_source_.Initialise(data.fluid);
  }
  return true;
}

void RemoveComponentFromDeck(UserData& data, const std::string& name) {
  data.fluid.RemoveComponent(name);
  data.bed.RemoveComponent(name, data.adsorbent_directory);  // strips it from every adsorbent
  ReconcileMoleFractions(data, name, false);
  data.fluid.adsorbed_heat_capacity_.erase(name);
}

void SeedLayerComposition(UserData& data, MixedAdsorbentLayer& layer) {
  const auto& dir = data.adsorbent_directory.GetAll();
  if (dir.empty()) return;
  AdsorbentID lowest = dir.begin()->first;
  for (const auto& [id, a] : dir) lowest = std::min(lowest, id);
  layer.AddAdsorbent(lowest, 1.0);
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

void RemoveAdsorbentFromDeck(UserData& data, AdsorbentID id) {
  data.RemoveAdsorbent(id);
  SeedEmptyLayerCompositions(data);
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
