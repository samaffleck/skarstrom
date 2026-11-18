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
#include "adsim/equations/ambient_air.h"
#include "adsim/equations/wall_temperature_boundary_condition.h"
#include "adsim/equations/mass_transfer_model.h"
#include "adsim/solver/initial_conditions.h"
#include "adsim/solver/run_health.h"

#include <algorithm>
#include <cmath>
#include <limits>
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

  std::vector<std::string> HeatOfAdsorptionWarnings(const AdsorbentDirectory& adsorbents, const Bed& bed,
                                                    const Fluid& fluid, double T_ref, double P_ref) {
    std::vector<std::string> warnings;
    if (bed.equations_.IsIsothermal()) return warnings;
    std::set<AdsorbentID> in_layer;
    for (const auto& layer : bed.layers_.layers_)
      for (const auto& [id, frac] : layer.composition_) in_layer.insert(id);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    bool any = false;
    for (AdsorbentID id : in_layer) {
      if (!adsorbents.Contains(id)) continue;
      const Adsorbent& ads = adsorbents.Get(id);
      const auto& source = ads.adsorption_source_;
      for (const auto& component : fluid.components_) {
        if (source.HeatOfAdsorptionModeOf(component) != HeatOfAdsorptionMode::Isotherm) continue;
        IIsotherm* iso = source.FindIsotherm(component);
        if (!iso || iso->GetIsothermType() == IsothermType::INERT) continue;
        any = true;
        const std::string who = "the " + component + " heat of adsorption of adsorbent '" + ads.name_ + "'";
        const double at_ref = IsostericHeatAtState(iso, component, T_ref, P_ref, 1.0, 0.0, nan);
        const double dilute = IsostericHeatAtState(iso, component, T_ref, 1.0, 1.0, 0.0, nan);
        if (!std::isfinite(at_ref) || !std::isfinite(dilute)) {
          std::ostringstream ss;
          ss << "Warning: " << who << " comes from the isotherm, but its isosteric heat cannot be evaluated at "
             << T_ref << " K (pure component at " << P_ref << " Pa, or at 1 Pa); the stored constant is used"
             << " wherever that happens.";
          warnings.push_back(ss.str());
        } else if (std::abs(at_ref) < 1.0 && std::abs(dilute) < 1.0) {
          warnings.push_back("Warning: " + who + " comes from the isotherm, but the isotherm ('" +
                             std::string(IsothermKindString(iso->GetIsothermType())) +
                             "' with these parameters) does not depend on temperature, so its isosteric heat is 0"
                             " and adsorption releases no heat. Use the constant mode, or an isotherm with a"
                             " temperature dependence.");
        }
      }
    }
    if (any && fluid.heat_of_adsorption_convention_ != HeatOfAdsorptionConvention::ISOSTERIC)
      warnings.push_back("Warning: some heats of adsorption come from the isotherm, and the fluid's heat of"
                         " adsorption convention is not isosteric. The convention applies only to constant"
                         " heats: a heat from the isotherm is the isosteric heat, and the energy balance uses it"
                         " as it is.");
    return warnings;
  }

  std::vector<std::string> LayerCompositionWarnings(const AdsorbentDirectory& adsorbents, const Bed& bed) {
    std::vector<std::string> warnings;
    for (const auto& layer : bed.layers_.layers_)
      for (const auto& issue : CheckLayerComposition(layer, adsorbents))
        if (issue.severity == LayerCompositionIssue::Severity::Warning) warnings.push_back(issue.message);
    return warnings;
  }

  std::vector<std::string> MassTransferWarnings(const AdsorbentDirectory& adsorbents, const Bed& bed) {
    std::set<AdsorbentID> in_layer;
    for (const auto& layer : bed.layers_.layers_)
      for (const auto& [id, frac] : layer.composition_) in_layer.insert(id);
    std::vector<std::string> errors, warnings;
    for (AdsorbentID id : in_layer) {
      if (!adsorbents.Contains(id)) continue;
      const Adsorbent& ads = adsorbents.Get(id);
      for (const auto& [component, model] : ads.adsorption_source_.mass_transfer_models_)
        if (model) CheckMassTransferModel(*model, "the " + component + " kinetics of " + ads.name_, errors, warnings);
    }
    for (auto& w : warnings) w = "Warning: " + w;
    return warnings;
  }

  std::vector<std::string> FilmCorrelationWarnings(const Bed& bed) {
    std::vector<std::string> warnings;
    for (const auto& layer : bed.layers_.layers_) {
      if (layer.IsStructured() || layer.film_correlation_ != FilmCorrelation::RanzMarshall) continue;
      warnings.push_back("Warning: layer '" + layer.name_ + "' uses the Ranz-Marshall film correlation"
                         " (Nu, Sh = 2 + 0.6 Re^(1/2) Pr^(1/3) or Sc^(1/3)), which is for a single sphere, not a"
                         " packed bed: in a bed it under-predicts the particle heat and mass transfer coefficients at"
                         " moderate Reynolds number. Wakao-Funazkri or Gunn are packed-bed correlations.");
    }
    return warnings;
  }

  std::vector<std::string> NaturalConvectionWarnings(const Bed& bed) {
    std::vector<std::string> warnings;
    if (bed.equations_.IsIsothermal()) return warnings;
    bool vertical = false;
    for (const auto& layer : bed.layers_.layers_) {
      const auto* cc = dynamic_cast<const ExternalHeatTransferModel_ChurchillChu*>(
          layer.equations_.external_heat_transfer_model_.get());
      if (cc && cc->Vertical()) vertical = true;
    }
    if (!vertical) return warnings;
    double height = 0.0;
    for (const auto& layer : bed.layers_.layers_) height += layer.domain_.length.Get();
    const double d_o = bed.GetHydraulicDiameter() + 2.0 * bed.wall_.thickness_;
    if (!(height > 0.0) || !(d_o > 0.0)) return warnings;
    const double T_inf = bed.ambient_temperature_;
    const double Gr = ambient_air::Grashof(T_inf + kNaturalConvectionCheckDeltaT, T_inf, height);
    const double d_min = 35.0 * height / std::pow(Gr, 0.25);
    if (d_o >= d_min) return warnings;
    std::ostringstream os;
    os << std::setprecision(3) << "Warning: the column (outer diameter " << d_o << " m, bed height " << height
       << " m) is too slender for the vertical-plate Churchill-Chu natural convection correlation: a plate needs"
          " an outer diameter of at least "
       << d_min << " m (D/L >= 35/Gr_L^(1/4) at a wall "
       << kNaturalConvectionCheckDeltaT << " K above ambient). The curvature of a slender cylinder raises the"
          " heat loss above the plate value, so the model under-predicts the wall-to-ambient coefficient. The"
          " yang_vertical_cylinder kind holds at any diameter-to-height ratio.";
    warnings.push_back(os.str());
    return warnings;
  }

  std::vector<std::string> StructuredLayerWarnings(const AdsorbentDirectory& adsorbents, const Bed& bed,
                                                   const Cycle* cycle) {
    (void)adsorbents;
    std::vector<std::string> warnings;
    // The packed-bed heat closures (Ruthven and Yagi-Kunii heat dispersion,
    // Dixon and Leva wall coefficients) are deck-check errors on a structured
    // layer (Bed::PerformChecks), so only the thermal mode is left to note.
    bool any_structured = false;
    for (const auto& layer : bed.layers_.layers_) any_structured = any_structured || layer.IsStructured();
    if (any_structured && bed.equations_.thermal_equilibrium_mode_ != ThermalEquilibriumMode::NON_LOCAL)
      warnings.push_back("Warning: the bed has a structured layer (monolith or laminate) and uses local thermal"
                         " equilibrium (LTE), which takes heat exchange between the gas and the thin channel walls"
                         " as instantaneous; local thermal non-equilibrium (LTNE) is recommended.");

    // A structured layer cut off from a wall that does exchange heat.
    if (any_structured && !bed.equations_.IsIsothermal()) {
      bool wall_exchanges = false;
      for (const auto& layer : bed.layers_.layers_) {
        const auto* ho = layer.equations_.external_heat_transfer_model_.get();
        const auto* c = dynamic_cast<const ConstantExternalHeatTransferModel*>(ho);
        if (ho && (!c || c->heat_transfer_coefficient_.Get() != 0.0)) wall_exchanges = true;
      }
      if (dynamic_cast<const WallTemperatureBoundaryConditionFixed*>(
              bed.equations_.wall_temperature_boundary_condition_.get()))
        wall_exchanges = true;
      if (cycle)
        for (const auto& step : cycle->step_list_)
          if (step.heating_medium_.enabled_) wall_exchanges = true;
      if (wall_exchanges) {
        for (const auto& layer : bed.layers_.layers_) {
          if (!layer.IsStructured()) continue;
          const auto* hi = dynamic_cast<const ConstantInternalHeatTransferModel*>(
              layer.equations_.internal_heat_transfer_model_.get());
          if (!hi || hi->heat_transfer_coefficient_.Get() != 0.0) continue;
          warnings.push_back("Warning: layer '" + layer.name_ + "': the bed-to-wall heat transfer coefficient is 0, so"
                             " the column wall exchanges no heat with this layer; set layer[" +
                             std::to_string(layer.id_) + "].internal_htc.");
        }
      }
    }
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

  namespace {

    // The largest relative saturation of a gas (y over `components`) at
    // pressure P and temperature T; NaN when no component can condense.
    struct SaturationAt { double S = std::numeric_limits<double>::quiet_NaN(); std::string component; };
    SaturationAt LargestSaturation(const std::vector<std::string>& components, const ExpressionMap& y, double P,
                                   double T) {
      SaturationAt out;
      for (const auto& c : components) {
        const auto it = y.find(c);
        if (it == y.end()) continue;
        const double S = RelativeSaturation(c, std::max(it->second.Get(), 0.0) * P, T);
        if (std::isfinite(S) && (std::isnan(out.S) || S > out.S)) out = {S, c};
      }
      return out;
    }

    std::string Fixed(double v, int digits) {
      std::ostringstream ss;
      ss << std::fixed << std::setprecision(digits) << v;
      return ss.str();
    }

    std::string SaturationText(const SaturationAt& s, double P, double T) {
      return s.component + " at " + Fixed(s.S, 2) + " times its saturation pressure (" + Fixed(P / 1e5, 3)
             + " bar, " + Fixed(T, 2) + " K)";
    }

    const char* kSupersaturationTail =
        " Condensation is not modelled: the gas keeps the supersaturated mole fraction and no heat of condensation "
        "is released (see SUPERSATURATION in docs/run-health.md).";

  } // namespace

  std::vector<std::string> SupersaturationWarnings(const Cycle& cycle, const Bed& bed, const Fluid& fluid,
                                                   const InitialConditions& ic) {
    std::vector<std::string> warnings;
    const auto& comps = fluid.components_;
    const bool isothermal = bed.equations_.IsIsothermal();
    const double T_iso = bed.equations_.GetIsothermalTemperature(ic.temperature.Get());
    const double P0 = ic.pressure.Get();

    // The initial gas, at the gas and at the wall temperature.
    {
      const double Tg = isothermal ? T_iso : ic.temperature.Get();
      const double Tw = isothermal ? T_iso : ic.wall_temperature.Get();
      const SaturationAt gas = LargestSaturation(comps, ic.mole_fraction, P0, Tg);
      const SaturationAt wall = LargestSaturation(comps, ic.mole_fraction, P0, Tw);
      if (gas.S > 1.0 || wall.S > 1.0) {
        const bool at_wall = !(gas.S > 1.0);  // the gas measure first
        warnings.push_back("Warning: the initial gas is supersaturated" + std::string(at_wall ? " at the initial wall temperature" : "")
                           + ": " + (at_wall ? SaturationText(wall, P0, Tw) : SaturationText(gas, P0, Tg)) + "."
                           + kSupersaturationTail);
      }
    }

    // Inlets that carry their own gas.
    for (const auto& step : cycle.step_list_) {
      const double t_end = step.step_time_.Get();
      for (const auto end : {BoundaryLocation::BOTTOM, BoundaryLocation::TOP}) {
        const auto& bc = end == BoundaryLocation::TOP ? step.top_boundary_condition_ : step.bottom_boundary_condition_;
        if (!bc || bc->direction_ != BoundaryDirection::INLET) continue;
        const std::string end_name = end == BoundaryLocation::TOP ? "top" : "bottom";
        auto check = [&](const std::string& what, const ExpressionMap& y, double P, double T) {
          const SaturationAt s = LargestSaturation(comps, y, P, T);
          if (s.S > 1.0)
            warnings.push_back("Warning: " + what + " (step \"" + step.step_name_ + "\", " + end_name
                               + ") is supersaturated: " + SaturationText(s, P, T) + "." + kSupersaturationTail);
        };
        if (const auto* pbc = dynamic_cast<const PressureBoundaryCondition*>(bc.get())) {
          if (pbc->connection_.IsConnected() || pbc->buffer_tank_.IsActive()) continue;
          const auto P = ValueAtStepEnd(*pbc);
          if (!P) continue;
          const double T = isothermal ? T_iso
                                      : std::max(pbc->AssumedInletTemperature(0.0), pbc->AssumedInletTemperature(t_end));
          check("the inlet gas of the pressure boundary", pbc->mole_fraction_, *P, T);
        } else if (const auto* fbc = dynamic_cast<const FlowBoundaryCondition*>(bc.get())) {
          // The bed pressure the stream enters at: the other end's pressure
          // boundary, else the initial pressure.
          double P = P0;
          if (const auto* other = PressureBC(step, Other(end)))
            if (const auto v = ValueAtStepEnd(*other)) P = *v;
          for (const auto& stream : fbc->streams_) {
            if (!stream || stream->connection_.IsConnected() || stream->buffer_tank_.IsActive()) continue;
            const double T = isothermal ? T_iso
                                        : std::max(stream->AssumedInletTemperature(0.0),
                                                   stream->AssumedInletTemperature(t_end));
            check("inlet stream \"" + stream->name_ + "\"", stream->mole_fraction, P, T);
          }
        }
      }
    }
    return warnings;
  }

  std::vector<std::string> ViscousMacroporeWarnings(const AdsorbentDirectory& adsorbents, const Bed& bed, const Fluid& fluid,
                                                    const InitialConditions& ic, const Cycle& cycle) {
    std::vector<std::string> warnings;
    const bool isothermal = bed.equations_.IsIsothermal();
    const double T = isothermal ? bed.equations_.GetIsothermalTemperature(ic.temperature.Get()) : ic.temperature.Get();
    if (!(T > 0.0)) return warnings;

    // The highest pressure of the deck: the initial pressure and every
    // pressure boundary value (a ramp's end value, a constant's value).
    double P_max = ic.pressure.Get();
    for (const auto& step : cycle.step_list_)
      for (const auto end : {BoundaryLocation::BOTTOM, BoundaryLocation::TOP})
        if (const auto* pbc = PressureBC(step, end))
          if (const auto v = ValueAtStepEnd(*pbc)) P_max = std::max(P_max, *v);
    if (!(P_max > 0.0)) return warnings;

    // The gas viscosity of the initial state: the constant value, or the mole
    // fraction average of the pure-component polynomials (as the solver does).
    double mu = fluid.ConstantViscosity();
    if (!(mu > 0.0)) {
      double sum_y = 0.0, sum_mu = 0.0;
      for (const auto& c : fluid.components_) {
        const auto it = ic.mole_fraction.find(c);
        if (it == ic.mole_fraction.end()) continue;
        const double y = std::max(it->second.Get(), 0.0);
        sum_y += y;
        sum_mu += y * Fluid::GetComponentViscosity(c, T);
      }
      mu = sum_y > 0.0 ? sum_mu / sum_y : 0.0;
    }
    if (!(mu > 0.0)) return warnings;

    constexpr double kPi = 3.14159265358979323846;
    constexpr double kR = 8.31446261815324;
    std::set<AdsorbentID> in_layer;
    for (const auto& layer : bed.layers_.layers_)
      for (const auto& [id, frac] : layer.composition_) in_layer.insert(id);
    for (AdsorbentID id : in_layer) {
      if (!adsorbents.Contains(id)) continue;
      const Adsorbent& ads = adsorbents.Get(id);
      for (const auto& [component, model] : ads.adsorption_source_.mass_transfer_models_) {
        if (!model) continue;
        const double M = GetMolarMass(component);
        if (!(M > 0.0)) continue;
        const double v_bar = std::sqrt(8.0 * kR * T / (kPi * M));
        auto check = [&](const IMassTransferCoefficientModel* mtc, const std::string& site) {
          const auto* g = dynamic_cast<const MassTransferCoefficientSeriesResistance*>(mtc);
          if (!g) return;
          const auto res = g->Resistances();
          if (!res.macropore || !res.knudsen) return;
          const double r = g->pore_radius_.Get();
          if (!(r > 0.0)) return;
          const double ratio = 3.0 * r * P_max / (16.0 * mu * v_bar);
          if (!(ratio > kViscousPoreFlowRatio)) return;
          std::ostringstream os;
          os << std::setprecision(3) << "Warning: the macropore resistance of " << component << " in " << ads.name_ << site
             << " ignores viscous (Poiseuille) flow in the pores, which may matter here: the ratio of viscous to"
                " Knudsen pore transport 3 r P / (16 mu v) is "
             << ratio << " (pore radius " << r << " m, highest pressure " << P_max / 1.0e5 << " bar, viscosity " << mu
             << " Pa s, mean molecular speed " << v_bar << " m/s), above " << kViscousPoreFlowRatio
             << ". A viscous term adds to the uptake only for a pressure gradient across the particle, so the model"
                " can under-predict the uptake of a pressure step; a larger pore radius or a higher pressure raises the"
                " ratio. See docs/kinetic-models.md.";
          warnings.push_back(os.str());
        };
        check(model->mtc_model_.get(), "");
        if (const auto* dual = dynamic_cast<const DualRateLinearDrivingForce*>(model.get()))
          if (dual->HasSlowSite()) check(dual->slow_mtc_model_.get(), " (slow site)");
      }
    }
    return warnings;
  }

} // namespace adsim
