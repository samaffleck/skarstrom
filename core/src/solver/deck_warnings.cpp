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

#include <cmath>
#include <iomanip>
#include <optional>
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

} // namespace adsim
