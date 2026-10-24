/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <string>
#include <unordered_map>

// Cereal includes
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/string.hpp>
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>

// adsim includes
#include "adsim/objects/fluid.h"
#include "adsim/math/expression_map.h"


namespace adsim {

  class InitialConditions {
  public:
    bool PerformChecks(const Fluid& fluid, std::string& error_message) const;

	  ExpressionParam temperature{293.0};       // [K], fluid and solid temperature
	  ExpressionParam wall_temperature{293.0};  // [K]
	  ExpressionParam pressure{101325.0};       // [Pa]
	  ExpressionParam velocity{0.0};            // [m/s]
	  ExpressionMap mole_fraction{}; // [mol/mol]
	  std::unordered_map<std::string, double> adsorbed_concentration{}; // [mol/kg]

    void EvaluateExpressions(const EvalContext& ctx) {
      temperature.Evaluate(ctx);
      wall_temperature.Evaluate(ctx);
      pressure.Evaluate(ctx);
      velocity.Evaluate(ctx);
      EvaluateExpressionMap(mole_fraction, ctx);
    }

  private:
    friend class ::cereal::access;
    template<class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::make_nvp("Initial_Temperature", temperature),
        cereal::make_nvp("Initial_Wall_Temperature", wall_temperature),
        cereal::make_nvp("Initial_Pressure", pressure),
        cereal::make_nvp("Initial_Velocity", velocity)
      );
      // version < 3: the mole fractions were plain numbers.
      SerializeExpressionMap(ar, "Initial_Composition", mole_fraction, version < 3);
      ar(cereal::make_nvp("Adsorbed_Concentration", adsorbed_concentration));
    }
  };

} // end adsim namespace

CEREAL_CLASS_VERSION(adsim::InitialConditions, 3);
