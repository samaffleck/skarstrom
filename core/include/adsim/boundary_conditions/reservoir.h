/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

#include <cstddef>
#include <string>

#include <cereal/types/string.hpp>
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>

#include "adsim/math/expression_map.h"
#include "adsim/math/expression_param.h"

namespace adsim {

  // A reservoir is a named gas state outside the bed: a SOURCE that inlets
  // draw from (a feed header, a purge gas supply) or a SINK that outlets
  // discharge to (ambient, a vacuum header, a product line). It holds the
  // state only: temperature, pressure and composition. How much gas passes,
  // and through what (a set flow, a valve, a pump, a pressure ramp), stays on
  // the endpoint that is attached to it, so two feed steps can share one
  // source and still run at different flows.
  //
  // An endpoint (a flow stream or a pressure boundary) names its reservoir by
  // ID (reservoir_id_, 0 = none). While attached, the endpoint's own
  // stream-side values are the reservoir's (Cycle::ApplyReservoirs copies
  // them before every expression evaluation and run):
  //   * its temperature (constant; a temperature ramp is dropped),
  //   * its mole fractions (an inlet's composition, an outlet's backflow
  //     composition),
  //   * its pressure target, when it has one: the pressure of a pressure
  //     boundary or the back pressure of a valve stream (the constant value,
  //     or the end value of a ramp; the ramp's shape and time stay the
  //     endpoint's). A set-flow or pump stream has no pressure value, so the
  //     reservoir's pressure is not used there.
  // The registry lists those values once, on the reservoir
  // (reservoir[ID].temperature, .pressure, .y[Component]), and not on the
  // attached endpoints. A sink lists its pressure only: an outlet is
  // zero-gradient in temperature and composition whatever the direction of
  // the flow, so the engine does not use a sink's temperature and mole
  // fractions.
  //
  // A connection to another step or a buffer tank still takes precedence for
  // the gas an inlet receives (SelectInletTemperature); the reservoir then
  // gives the first-cycle assumption, as the endpoint's own values do.
  enum class ReservoirKind {
    SOURCE = 0,  // inlets attach to it
    SINK = 1     // outlets attach to it
  };

  static const char* ReservoirKindList[] = {
    "Source",  // 0
    "Sink"     // 1
  };

  class Reservoir {
  public:
    // Stable identity within the cycle (0 = not yet assigned, see
    // Cycle::EnsureIds). Endpoints and registry paths ("reservoir[<id>]")
    // refer to it, so a rename does not break them.
    size_t id_ = 0;
    std::string name_{};
    ReservoirKind kind_ = ReservoirKind::SOURCE;
    ExpressionParam temperature_{293.0};  // [K]
    ExpressionParam pressure_{1e5};       // [Pa]
    ExpressionMap mole_fraction_{};       // [-] per component

    void EvaluateExpressions(const EvalContext& ctx) {
      temperature_.Evaluate(ctx);
      pressure_.Evaluate(ctx);
      EvaluateExpressionMap(mole_fraction_, ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/) {
      ar(
        cereal::make_nvp("ID", id_),
        cereal::make_nvp("Name", name_),
        cereal::make_nvp("Kind", kind_),
        cereal::make_nvp("Temperature", temperature_),
        cereal::make_nvp("Pressure", pressure_)
      );
      SerializeExpressionMap(ar, "Mole_Fraction", mole_fraction_, false);
    }
  };

}  // namespace adsim

CEREAL_CLASS_VERSION(adsim::Reservoir, 1);
