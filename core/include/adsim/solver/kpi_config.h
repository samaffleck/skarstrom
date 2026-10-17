#pragma once

// STL includes
#include <algorithm>
#include <optional>
#include <string>
#include <vector>

// adsim includes
#include "adsim/boundary_conditions/flow_boundary_value.h"
#include "adsim/boundary_conditions/step.h"

// Cereal includes
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/cereal.hpp>


namespace adsim {

  // A single KPI throughput target: one endpoint (a flow stream, or a pressure
  // boundary named kPressureBoundaryEndpointName) on a given cycle step. The
  // location records which bed end the endpoint sits on so the KPI maths can
  // convert the stream's velocity-signed molar flow into "entering the bed"
  // (positive) / "leaving the bed" (negative). Refreshed by Prune().
  //
  // Identity: endpoint_id (the stream's / pressure boundary's endpoint_id_)
  // is the reference. step_index, endpoint_name and location are caches that
  // Resolve() refreshes from it, so renaming or reordering steps and streams
  // does not break the ref. A ref without IDs (an old deck, the web document,
  // SetRole) is resolved once by (step_index, endpoint_name) and gets its IDs
  // then.
  struct KPIEndpointRef {
    int step_index = 0;
    std::string endpoint_name{};
    BoundaryLocation location = BoundaryLocation::BOTTOM;
    size_t step_id = 0;
    size_t endpoint_id = 0;

    bool Matches(int step, const std::string& name) const {
      return step_index == step && endpoint_name == name;
    }

    // Re-binds the ref to the current cycle: by endpoint_id when set,
    // otherwise by (step_index, endpoint_name). Returns false (ref unchanged)
    // when it does not resolve.
    bool Resolve(const std::vector<Step>& steps) {
      if (endpoint_id != 0) {
        for (size_t i = 0; i < steps.size(); ++i) {
          for (const auto& ep : steps[i].GetEndpoints()) {
            if (ep.id != endpoint_id) continue;
            step_index = static_cast<int>(i);
            step_id = steps[i].id_;
            endpoint_name = ep.name;
            location = ep.location;
            return true;
          }
        }
        return false;
      }
      if (step_index < 0 || step_index >= static_cast<int>(steps.size())) return false;
      const Step& s = steps[static_cast<size_t>(step_index)];
      // A name shared by both ends of the step: the cached location decides.
      const std::vector<Step::Endpoint> eps = s.GetEndpoints();
      const Step::Endpoint* hit = nullptr;
      for (const auto& ep : eps) {
        if (ep.name != endpoint_name) continue;
        if (!hit || (ep.location == location && hit->location != location)) hit = &ep;
      }
      if (!hit) return false;
      location = hit->location;
      step_id = s.id_;
      endpoint_id = hit->id;
      return true;
    }

  private:
    friend class ::cereal::access;
    template<class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::make_nvp("Step", step_index),
        cereal::make_nvp("Endpoint", endpoint_name),
        cereal::make_nvp("Location", location)
      );
      if (version >= 1) {
        ar(
          cereal::make_nvp("Step_ID", step_id),
          cereal::make_nvp("Endpoint_ID", endpoint_id)
        );
      }
      // version 0: IDs stay 0 and are resolved by name on load.
    }
  };

  // Per-component KPI configuration: which endpoints supply this component's
  // feed, product and recycle throughput. A recycle endpoint carries product
  // back into the bed (a purge or light-product pressurisation fed from the
  // product): it does not change the product's purity, but its moles are
  // subtracted from the product for recovery, productivity and production
  // rate (net product).
  struct ComponentKPIConfig {
    std::string product_component{};
    bool tracked = false;                                    // included in KPI tabs/results when true

    std::vector<KPIEndpointRef> feed_endpoints{};
    std::vector<KPIEndpointRef> product_endpoints{};
    std::vector<KPIEndpointRef> recycle_endpoints{};

    static bool Contains(const std::vector<KPIEndpointRef>& v, int step, const std::string& name) {
      return std::any_of(v.begin(), v.end(),
                         [&](const KPIEndpointRef& r) { return r.Matches(step, name); });
    }

    bool IsFeed(int step, const std::string& name) const { return Contains(feed_endpoints, step, name); }
    bool IsProduct(int step, const std::string& name) const { return Contains(product_endpoints, step, name); }
    bool IsRecycle(int step, const std::string& name) const { return Contains(recycle_endpoints, step, name); }

    bool HasAnyFeed() const { return !feed_endpoints.empty(); }
    bool HasAnyProduct() const { return !product_endpoints.empty(); }

    // Add/remove an endpoint from a role list (idempotent). Updates the cached
    // location when toggling on an already-present ref.
    static void SetRole(std::vector<KPIEndpointRef>& v, int step, const std::string& name,
                        BoundaryLocation loc, bool on) {
      auto it = std::find_if(v.begin(), v.end(),
                             [&](const KPIEndpointRef& r) { return r.Matches(step, name); });
      if (on) {
        if (it == v.end()) v.push_back({step, name, loc});
        else               it->location = loc;
      } else if (it != v.end()) {
        v.erase(it);
      }
    }

    void SetFeed(int step, const std::string& name, BoundaryLocation loc, bool on) {
      SetRole(feed_endpoints, step, name, loc, on);
    }
    void SetProduct(int step, const std::string& name, BoundaryLocation loc, bool on) {
      SetRole(product_endpoints, step, name, loc, on);
    }
    void SetRecycle(int step, const std::string& name, BoundaryLocation loc, bool on) {
      SetRole(recycle_endpoints, step, name, loc, on);
    }

    // ---- Structural cycle edits ------------------------------------------
    // Refs are keyed by step index, so every edit that moves steps must remap
    // them, or a ref silently resolves to whichever step now sits at the old
    // index when it has a same-named endpoint ("Stream 1"). Callers should use
    // the UserData wrappers (DuplicateStep / RemoveStep / SwapSteps /
    // RenameStep / SwapStepBoundaries / RenameStream), which pair the cycle
    // edit with the matching remap here.

    // Apply `map(old_index) -> new_index` to every ref; a negative result
    // drops the ref.
    template<class Fn>
    void RemapStepIndices(Fn map) {
      auto remap = [&](std::vector<KPIEndpointRef>& v) {
        std::vector<KPIEndpointRef> kept;
        kept.reserve(v.size());
        for (auto& r : v) {
          const int idx = map(r.step_index);
          if (idx < 0) continue;
          r.step_index = idx;
          kept.push_back(std::move(r));
        }
        v = std::move(kept);
      };
      remap(feed_endpoints);
      remap(product_endpoints);
      remap(recycle_endpoints);
    }

    // A step was inserted at `index`: refs at or after it shift +1.
    void OnStepInserted(int index) {
      RemapStepIndices([index](int i) { return i >= index ? i + 1 : i; });
    }

    // The step at `index` was removed: its refs are dropped, later ones shift -1.
    void OnStepRemoved(int index) {
      RemapStepIndices([index](int i) { return i == index ? -1 : (i > index ? i - 1 : i); });
    }

    // Steps `a` and `b` exchanged positions: their refs follow them.
    void OnStepsSwapped(int a, int b) {
      RemapStepIndices([a, b](int i) { return i == a ? b : (i == b ? a : i); });
    }

    // The step at `index` was renamed. Pressure-boundary endpoint names embed
    // the step name ("<step> Pressure Boundary Top|Bottom"), so rewrite those
    // refs to the new name; flow-stream refs are unaffected.
    void OnStepRenamed(int index, const std::string& old_name, const std::string& new_name) {
      auto rename = [&](std::vector<KPIEndpointRef>& v) {
        for (auto& r : v) {
          if (r.step_index != index) continue;
          for (BoundaryLocation loc : { BoundaryLocation::TOP, BoundaryLocation::BOTTOM }) {
            if (r.endpoint_name == Step::PressureEndpointName(old_name, loc)) {
              r.endpoint_name = Step::PressureEndpointName(new_name, loc);
              break;
            }
          }
        }
      };
      rename(feed_endpoints);
      rename(product_endpoints);
      rename(recycle_endpoints);
    }

    // The top and bottom boundary conditions of step `index` (named
    // `step_name`) were exchanged. Every endpoint on that step moves to the
    // other bed end: a pressure-boundary ref's side-keyed name flips
    // Top<->Bottom, and every ref on the step flips its location (the KPI
    // maths uses it to sign the flow as entering/leaving the bed). Flow-stream
    // names are unchanged. Mirrors swapKpiEndpointRef (web).
    void OnStepBoundariesSwapped(int index, const std::string& step_name) {
      const std::string top = Step::PressureEndpointName(step_name, BoundaryLocation::TOP);
      const std::string bottom = Step::PressureEndpointName(step_name, BoundaryLocation::BOTTOM);
      auto swap = [&](std::vector<KPIEndpointRef>& v) {
        for (auto& r : v) {
          if (r.step_index != index) continue;
          if (r.endpoint_name == top)         r.endpoint_name = bottom;
          else if (r.endpoint_name == bottom) r.endpoint_name = top;
          r.location = (r.location == BoundaryLocation::TOP) ? BoundaryLocation::BOTTOM
                                                             : BoundaryLocation::TOP;
        }
      };
      swap(feed_endpoints);
      swap(product_endpoints);
      swap(recycle_endpoints);
    }

    // A flow stream on step `index` was renamed from `old_name` to `new_name`.
    // Flow-stream refs are keyed by stream name, so rewrite them. When
    // `location` is given, only refs at that bed end are rewritten: the caller
    // passes it when another endpoint on the same step still carries
    // `old_name`, whose refs must stay put.
    void OnStreamRenamed(int index, const std::string& old_name, const std::string& new_name,
                         std::optional<BoundaryLocation> location = std::nullopt) {
      auto rename = [&](std::vector<KPIEndpointRef>& v) {
        for (auto& r : v) {
          if (r.step_index != index || r.endpoint_name != old_name) continue;
          if (location && r.location != *location) continue;
          r.endpoint_name = new_name;
        }
      };
      rename(feed_endpoints);
      rename(product_endpoints);
      rename(recycle_endpoints);
    }

    // Re-bind every ref to the current cycle (KPIEndpointRef::Resolve: by ID,
    // or by step index + name for refs without IDs) and drop the refs whose
    // endpoint no longer exists. Replaces the old fixed-size Resize().
    void Prune(const std::vector<Step>& steps) {
      auto prune = [&](std::vector<KPIEndpointRef>& v) {
        v.erase(std::remove_if(v.begin(), v.end(), [&](KPIEndpointRef& r) { return !r.Resolve(steps); }),
                v.end());
      };
      prune(feed_endpoints);
      prune(product_endpoints);
      prune(recycle_endpoints);
    }

    // As Prune, but keeps refs that do not resolve (used on load and before
    // computing KPIs, where nothing should be silently discarded).
    void ResolveIds(const std::vector<Step>& steps) {
      for (auto& r : feed_endpoints) r.Resolve(steps);
      for (auto& r : product_endpoints) r.Resolve(steps);
      for (auto& r : recycle_endpoints) r.Resolve(steps);
    }

  private:
    friend class ::cereal::access;

    template<class Archive>
    void save(Archive& ar, std::uint32_t const /*version*/) const {
      ar(
        cereal::make_nvp("Product_Component", product_component),
        cereal::make_nvp("Tracked", tracked),
        cereal::make_nvp("Feed_Endpoints", feed_endpoints),
        cereal::make_nvp("Product_Endpoints", product_endpoints),
        cereal::make_nvp("Recycle_Endpoints", recycle_endpoints)
      );
    }

    template<class Archive>
    void load(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::make_nvp("Product_Component", product_component),
        cereal::make_nvp("Tracked", tracked)
      );
      if (version >= 1) {
        ar(
          cereal::make_nvp("Feed_Endpoints", feed_endpoints),
          cereal::make_nvp("Product_Endpoints", product_endpoints)
        );
        if (version >= 2) ar(cereal::make_nvp("Recycle_Endpoints", recycle_endpoints));
        return;
      }

      // Legacy v0: per-step feed/product/recycle checkbox + Top/Bottom vectors.
      // Clean break — consume the old fields (so the archive stays in sync) and
      // discard them; the user re-specifies endpoints. See CLAUDE.md decision.
      std::vector<bool> discard_bool;
      std::vector<int>  discard_int;
      ar(
        cereal::make_nvp("Is_Feed_Step", discard_bool),
        cereal::make_nvp("Is_Product_Step", discard_bool),
        cereal::make_nvp("Is_Recycled_Product_Step", discard_bool),
        cereal::make_nvp("Feed_Location", discard_int),
        cereal::make_nvp("Product_Location", discard_int),
        cereal::make_nvp("Recycled_Product_Location", discard_int)
      );
    }
  };

  // KPI configuration container: one ComponentKPIConfig per fluid component the
  // user may track. Synchronised with the fluid component list on demand.
  struct KPIConfig {
    std::vector<ComponentKPIConfig> per_component;

    // Ensure there is exactly one entry per component name (in the given order),
    // creating missing entries (tracked=false by default) and removing entries
    // whose component is no longer present.
    void SyncToComponents(const std::vector<std::string>& components) {
      std::vector<ComponentKPIConfig> rebuilt;
      rebuilt.reserve(components.size());
      for (const auto& name : components) {
        auto it = std::find_if(
          per_component.begin(), per_component.end(),
          [&](const ComponentKPIConfig& c) { return c.product_component == name; }
        );
        if (it != per_component.end()) {
          rebuilt.push_back(std::move(*it));
        } else {
          ComponentKPIConfig fresh;
          fresh.product_component = name;
          rebuilt.push_back(std::move(fresh));
        }
      }
      per_component = std::move(rebuilt);
    }

    // Drop endpoint refs that no longer resolve against the current cycle and
    // refresh survivors' locations. Replaces the old Resize(n_steps).
    void Prune(const std::vector<Step>& steps) {
      for (auto& c : per_component) c.Prune(steps);
    }
    void ResolveIds(const std::vector<Step>& steps) {
      for (auto& c : per_component) c.ResolveIds(steps);
    }

    // Structural cycle edits: see ComponentKPIConfig::OnStep*.
    void OnStepInserted(int index) { for (auto& c : per_component) c.OnStepInserted(index); }
    void OnStepRemoved(int index)  { for (auto& c : per_component) c.OnStepRemoved(index); }
    void OnStepsSwapped(int a, int b) { for (auto& c : per_component) c.OnStepsSwapped(a, b); }
    void OnStepRenamed(int index, const std::string& old_name, const std::string& new_name) {
      for (auto& c : per_component) c.OnStepRenamed(index, old_name, new_name);
    }
    void OnStepBoundariesSwapped(int index, const std::string& step_name) {
      for (auto& c : per_component) c.OnStepBoundariesSwapped(index, step_name);
    }
    void OnStreamRenamed(int index, const std::string& old_name, const std::string& new_name,
                         std::optional<BoundaryLocation> location = std::nullopt) {
      for (auto& c : per_component) c.OnStreamRenamed(index, old_name, new_name, location);
    }

    ComponentKPIConfig* FindByComponent(const std::string& name) {
      for (auto& c : per_component)
        if (c.product_component == name) return &c;
      return nullptr;
    }

    const ComponentKPIConfig* FindByComponent(const std::string& name) const {
      for (const auto& c : per_component)
        if (c.product_component == name) return &c;
      return nullptr;
    }

    // First tracked entry, or nullptr if none tracked. Used as the default
    // target when callers don't specify a component (legacy paths).
    const ComponentKPIConfig* FirstTracked() const {
      for (const auto& c : per_component)
        if (c.tracked) return &c;
      return nullptr;
    }

    ComponentKPIConfig* FirstTracked() {
      for (auto& c : per_component)
        if (c.tracked) return &c;
      return nullptr;
    }

  private:
    friend class ::cereal::access;

    template<class Archive>
    void save(Archive& ar, std::uint32_t const /*version*/) const {
      ar(cereal::make_nvp("Per_Component", per_component));
    }

    template<class Archive>
    void load(Archive& ar, std::uint32_t const version) {
      if (version >= 2) {
        ar(cereal::make_nvp("Per_Component", per_component));
        return;
      }

      // Legacy v0/v1: a single flat config. Clean break — consume and discard.
      int discard_index = 0;
      std::vector<bool> discard_bool;
      std::vector<int>  discard_int;
      ar(
        cereal::make_nvp("Product_Component_Index", discard_index),
        cereal::make_nvp("Is_Feed_Step", discard_bool),
        cereal::make_nvp("Is_Product_Step", discard_bool),
        cereal::make_nvp("Feed_Location", discard_int),
        cereal::make_nvp("Product_Location", discard_int)
      );
      if (version >= 1) {
        ar(
          cereal::make_nvp("Is_Recycled_Product_Step", discard_bool),
          cereal::make_nvp("Recycled_Product_Location", discard_int)
        );
      }
      per_component.clear();
      legacy_pending_index_ = -1;
    }

  public:
    // After loading a project + the Fluid component list, call this to reconcile
    // the config against the components. Legacy decks come back with no KPI
    // selections (clean break); the user re-specifies. Safe to call repeatedly.
    void FinaliseLegacyLoad(const std::vector<std::string>& components) {
      legacy_pending_index_ = -1;
      SyncToComponents(components);
    }

  private:
    int legacy_pending_index_ = -1;
  };

} // end adsim namespace

CEREAL_CLASS_VERSION(adsim::KPIEndpointRef, 1);
CEREAL_CLASS_VERSION(adsim::ComponentKPIConfig, 2);
CEREAL_CLASS_VERSION(adsim::KPIConfig, 3);
