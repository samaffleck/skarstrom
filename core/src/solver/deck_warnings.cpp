/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/solver/deck_warnings.h"
#include "adsim/solver/cycle.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_value.h"
#include "adsim/objects/reactor.h"
#include "adsim/objects/wall_design.h"
#include "adsim/objects/porous_media.h"
#include "adsim/objects/fluid.h"
#include "adsim/equations/isotherm_kinds.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/equations/mass_transfer_model.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <iomanip>
#include <optional>
#include <set>
#include <sstream>

namespace adsim {

  namespace {

    const PressureBoundaryCondition* PressureBC(const Step& step, BoundaryLocation end) {
      const auto& bc = (end == BoundaryLocation::TOP) ? step.top_boundary_condition_
                                                      : step.bottom_boundary_condition_;
      return dynamic_cast<const PressureBoundaryCondition*>(bc.get());
    }

    // The value a pressure boundary holds when its step ends (a ramp's end
    // value, a constant's value); nullopt for a kind this check does not know.
    std::optional<double> ValueAtStepEnd(const PressureBoundaryCondition& bc) {
      const IBoundaryValue* v = bc.pressure_.get();
      if (!v) return std::nullopt;
      if (const auto* c = dynamic_cast<const ConstantBoundaryValue*>(v)) return c->constant_value_.Get();
      if (const auto* l = dynamic_cast<const LinearBoundaryValue*>(v)) return l->end_value_.Get();
      if (const auto* e = dynamic_cast<const ExponentialBoundaryValue*>(v)) return e->end_value_.Get();
      if (const auto* g = dynamic_cast<const LogLinearBoundaryValue*>(v)) return g->end_value_.Get();
      return std::nullopt;
    }

    std::optional<double> ConstantValue(const PressureBoundaryCondition& bc) {
      const auto* c = dynamic_cast<const ConstantBoundaryValue*>(bc.pressure_.get());
      if (!c) return std::nullopt;
      return c->constant_value_.Get();
    }

    BoundaryLocation Other(BoundaryLocation end) {
      return end == BoundaryLocation::TOP ? BoundaryLocation::BOTTOM : BoundaryLocation::TOP;
    }

    // The bed pressure step k is expected to start at, from the deck alone:
    // the end value of the most recent pressure boundary before it, looking
    // back cyclically. Also returns the step that set it.
    struct PriorPressure {
      double value;
      const Step* step;
    };
    std::optional<PriorPressure> PriorBoundaryPressure(const std::vector<Step>& steps, size_t k,
                                                       BoundaryLocation end) {
      const size_t n = steps.size();
      for (size_t back = 1; back < n; ++back) {
        const Step& prev = steps[(k + n - back) % n];
        for (BoundaryLocation at : { end, Other(end) }) {
          if (const auto* bc = PressureBC(prev, at)) {
            if (auto p = ValueAtStepEnd(*bc)) return PriorPressure{ *p, &prev };
          }
        }
      }
      return std::nullopt;
    }

    std::string Bar(double pascal) {
      std::ostringstream ss;
      ss << std::setprecision(3) << pascal / 1e5 << " bar";
      return ss.str();
    }

    bool IsStepChange(double boundary, double reference, double relative_threshold) {
      if (!std::isfinite(boundary) || !std::isfinite(reference)) return false;
      const double scale = std::max(std::abs(reference), 1.0);
      return std::abs(boundary - reference) / scale > relative_threshold;
    }

    std::string Warning(const Step& step, BoundaryLocation end, double boundary, const std::string& from) {
      std::ostringstream ss;
      ss << "Warning: step \"" << step.step_name_ << "\" ("
         << (end == BoundaryLocation::TOP ? "top" : "bottom")
         << ") starts with a step change in boundary pressure: its constant "
         << Bar(boundary) << " against " << from
         << ". The solver has to resolve an instantaneous pressure change there; "
         << "a linear or exponential ramp to " << Bar(boundary) << " avoids it.";
      return ss.str();
    }

  } // namespace

  std::vector<std::string> BoundaryPressureStepChangeWarnings(const Cycle& cycle, double initial_pressure,
                                                              double relative_threshold) {
    std::vector<std::string> warnings;
    const auto& steps = cycle.step_list_;
    for (size_t k = 0; k < steps.size(); ++k) {
      for (BoundaryLocation end : { BoundaryLocation::BOTTOM, BoundaryLocation::TOP }) {
        const auto* bc = PressureBC(steps[k], end);
        if (!bc) continue;
        const auto boundary = ConstantValue(*bc);
        if (!boundary) continue;  // a ramp starts from the bed pressure

        // The first step of the first cycle starts from the initial conditions.
        if (k == 0 && IsStepChange(*boundary, initial_pressure, relative_threshold)) {
          warnings.push_back(Warning(steps[k], end, *boundary,
                                     "the initial pressure " + Bar(initial_pressure)));
        }
        const auto prior = PriorBoundaryPressure(steps, k, end);
        if (!prior || prior->step == &steps[k]) continue;
        if (IsStepChange(*boundary, prior->value, relative_threshold)) {
          warnings.push_back(Warning(steps[k], end, *boundary,
                                     "about " + Bar(prior->value) + " at the end of step \""
                                     + prior->step->step_name_ + "\""));
        }
      }
    }
    return warnings;
  }

  std::vector<std::string> WallDesignWarnings(const Bed& bed) {
    std::vector<std::string> warnings;
    const auto& wall = bed.wall_;
    if (!wall.DesignInUse()) return warnings;
    const WallDesignResult design = DesignBedWall(bed);
    if (!design.valid) {
      warnings.push_back("Warning: the wall's minimum thickness cannot be checked: " + design.message + ".");
      return warnings;
    }
    const double t = wall.thickness_.Get();
    if (t < design.required_thickness * (1.0 - 1e-9)) {
      std::ostringstream ss;
      ss << std::setprecision(3) << "Warning: the wall is " << t * 1e3 << " mm thick, thinner than the "
         << design.required_thickness * 1e3 << " mm " << kWallDesignCodeList[static_cast<int>(wall.design_code_)]
         << " requires for its design" << (design.external_governs ? " (full vacuum governs)" : "")
         << ". Use Calculate in the wall design, or a thicker wall.";
      warnings.push_back(ss.str());
    }
    return warnings;
  }


  std::vector<std::string> MixingRuleWarnings(const AdsorbentDirectory& adsorbents, const Fluid& fluid) {
    std::vector<std::string> warnings;
    std::vector<AdsorbentID> ids;
    for (const auto& [id, ads] : adsorbents.GetAll()) ids.push_back(id);
    std::sort(ids.begin(), ids.end());
    for (AdsorbentID id : ids) {
      const Adsorbent& ads = adsorbents.Get(id);
      const auto& source = ads.adsorption_source_;
      if (source.mixing_rule_ != MixingRule::IAST) continue;
      std::map<IsothermType, std::vector<std::string>> extended;
      for (const auto& component : fluid.components_) {
        auto it = source.IsothermMap().find(component);
        if (it == source.IsothermMap().end()) continue;
        if (IsCompetitiveIsothermKind(it->second)) extended[it->second].push_back(component);
      }
      for (const auto& [type, components] : extended) {
        std::ostringstream ss;
        ss << "Warning: adsorbent '" << ads.name_ << "' uses the IAST mixing rule with the competitive isotherm '"
           << isothermModelNames[static_cast<int>(type)] << "' (" << IsothermKindString(type) << ") for ";
        for (size_t i = 0; i < components.size(); ++i) ss << (i ? ", " : "") << components[i];
        ss << ". IAST takes every isotherm in its pure-component limit, so that kind's own competition term"
              " is not applied";
        IsothermType pure{};
        if (PairedIsothermKind(type, pure))
          ss << "; use '" << isothermModelNames[static_cast<int>(pure)] << "' (" << IsothermKindString(pure)
             << "), which takes the same parameters.";
        else
          ss << "; give these components a pure-component isotherm.";
        warnings.push_back(ss.str());
      }
    }
    return warnings;
  }

  std::vector<std::string> LayerCompositionWarnings(const AdsorbentDirectory& adsorbents, const Bed& bed) {
    std::vector<std::string> warnings;
    for (const auto& layer : bed.layers_.layers_)
      for (const auto& issue : CheckLayerComposition(layer, adsorbents))
        if (issue.severity == LayerCompositionIssue::Severity::Warning) warnings.push_back(issue.message);
    return warnings;
  }

  std::vector<std::string> StructuredLayerWarnings(const AdsorbentDirectory& adsorbents, const Bed& bed) {
    std::vector<std::string> warnings;
    bool any_structured = false;
    for (const auto& layer : bed.layers_.layers_) {
      if (!layer.IsStructured()) continue;
      any_structured = true;
      const std::string what = std::string(LayerKindList[static_cast<int>(layer.kind_)]);
      std::vector<std::string> uses;
      const auto& eqs = layer.equations_;
      if (dynamic_cast<const DispersionModelRuthven*>(eqs.dispersion_model_.get()))
        uses.push_back("Ruthven axial dispersion");
      if (dynamic_cast<const HeatDispersionModelRuthven*>(eqs.heat_dispersion_model_.get()))
        uses.push_back("Ruthven heat dispersion");
      if (dynamic_cast<const HeatDispersionModelYagiandKunii*>(eqs.heat_dispersion_model_.get()))
        uses.push_back("Yagi-Kunii heat dispersion");
      if (dynamic_cast<const InternalHeatTransferModel_Dixon*>(eqs.internal_heat_transfer_model_.get()))
        uses.push_back("Dixon wall heat transfer");
      if (dynamic_cast<const InternalHeatTransferModel_Leva*>(eqs.internal_heat_transfer_model_.get()))
        uses.push_back("Leva wall heat transfer");
      if (bed.equations_.thermal_equilibrium_mode_ == ThermalEquilibriumMode::NON_LOCAL &&
          dynamic_cast<const ConstantParticleHeatTransferModel*>(bed.equations_.particle_heat_transfer_model_.get()))
        uses.push_back("the gas-solid heat transfer area 6 (1 - eb) / dp");
      for (const auto& [id, frac] : layer.composition_) {
        if (!adsorbents.Contains(id)) continue;
        const Adsorbent& ads = adsorbents.Get(id);
        std::vector<std::string> comps;
        std::vector<std::string> names;
        for (const auto& [c, model] : ads.adsorption_source_.mass_transfer_models_) names.push_back(c);
        std::sort(names.begin(), names.end());
        for (const auto& c : names) {
          const auto& model = ads.adsorption_source_.mass_transfer_models_.at(c);
          const auto* k = model ? model->mtc_model_.get() : nullptr;
          if (dynamic_cast<const MassTransferCoefficientMacropore*>(k) ||
              dynamic_cast<const MassTransferCoefficientMacroporeKnudsen*>(k) ||
              dynamic_cast<const MassTransferCoefficientMacroporeCapacity*>(k) ||
              dynamic_cast<const MassTransferCoefficient_Haynes*>(k))
            comps.push_back(c);
        }
        if (comps.empty()) continue;
        std::string list;
        for (size_t i = 0; i < comps.size(); ++i) list += (i ? ", " : "") + comps[i];
        uses.push_back("the particle-radius rate coefficient of '" + ads.name_ + "' for " + list);
      }
      if (uses.empty()) continue;
      std::ostringstream ss;
      ss << "Warning: layer '" << layer.name_ << "' is a " << what << " but uses packed-bed closures that read the"
         << " particle diameter (" << std::setprecision(4) << layer.particle_diameter_.Get() << " m): ";
      for (size_t i = 0; i < uses.size(); ++i) ss << (i ? "; " : "") << uses[i];
      ss << ". Closures for channels are not available yet; set the particle diameter to the length these"
            " correlations should use, or choose constant coefficients.";
      warnings.push_back(ss.str());
    }
    if (any_structured && bed.equations_.thermal_equilibrium_mode_ != ThermalEquilibriumMode::NON_LOCAL)
      warnings.push_back("Warning: the bed has a structured layer (monolith or laminate) and uses local thermal"
                         " equilibrium (LTE), which takes heat exchange between the gas and the thin channel walls"
                         " as instantaneous; local thermal non-equilibrium (LTNE) is recommended.");
    return warnings;
  }

  std::vector<std::string> ReactionWarnings(const AdsorbentDirectory& adsorbents, const Bed& bed,
                                            bool order_notes) {
    std::vector<std::string> warnings;
    std::set<AdsorbentID> in_layer;
    for (const auto& layer : bed.layers_.layers_)
      for (const auto& [id, frac] : layer.composition_) in_layer.insert(id);
    std::vector<AdsorbentID> ids;
    for (const auto& [id, ads] : adsorbents.GetAll()) ids.push_back(id);
    std::sort(ids.begin(), ids.end());
    const bool isothermal = bed.equations_.IsIsothermal();
    for (AdsorbentID id : ids) {
      const Adsorbent& ads = adsorbents.Get(id);
      const auto& reactions = ads.reaction_source_.reactions_;
      if (reactions.empty()) continue;
      if (!in_layer.count(id)) {
        warnings.push_back("Warning: solid '" + ads.name_ + "' has " + std::to_string(reactions.size()) +
                           (reactions.size() == 1 ? " reaction" : " reactions") +
                           " but no layer uses it, so they take no part in the run.");
        continue;
      }
      for (const auto& r : reactions) {
        const std::string who = "reaction '" + r.name_ + "' (on " + ads.name_ + ")";
        if (isothermal && r.heat_of_reaction_.Get() != 0.0)
          warnings.push_back("Warning: " + who + " has a nonzero heat of reaction, but the run is isothermal,"
                             " so the heat is not applied.");
        const ElementBalance eb = r.CheckElements();
        if (!eb.checked) {
          std::string list;
          for (size_t i = 0; i < eb.unknown.size(); ++i) list += (i ? ", " : "") + eb.unknown[i];
          warnings.push_back("Note: the element balance of " + who + " is not checked: no element formula for " +
                             list + ".");
        }
        if (!order_notes || !r.rate_) continue;
        std::vector<std::string> differing;
        for (auto& [leaf, param] : ReactionRateFields(r.rate_.get())) {
          if (leaf.rfind("order[", 0) != 0) continue;
          const std::string c = leaf.substr(6, leaf.size() - 7);
          auto it = r.nu_.find(c);
          const double nu = it == r.nu_.end() ? 0.0 : it->second;
          const double elementary = nu < 0.0 ? -nu : 0.0;
          if (param->Get() != elementary) differing.push_back(c);
        }
        if (!differing.empty()) {
          std::string list;
          for (size_t i = 0; i < differing.size(); ++i) list += (i ? ", " : "") + differing[i];
          warnings.push_back("Note: " + who + " is not elementary: the orders of " + list +
                             " differ from -nu (reactants) or 0 (products).");
        }
      }
    }
    return warnings;
  }

} // namespace adsim
