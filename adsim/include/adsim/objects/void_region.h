/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

#include "adsim/math/expression_param.h"

#include <cereal/access.hpp>
#include <cereal/cereal.hpp>
#include <cstdint>


namespace adsim {

  // Lumped 0D well-mixed void region (e.g. manifold/dead volume between the
  // inlet/outlet valve and the packed bed at one end of a column). When
  // enabled, contributes one cell to the bed's global cell array.
  //
  // The cell has fluid (gas) only — no adsorbent, so solid temperature and
  // adsorbed concentrations are not defined at the void cell. Treated as
  // adiabatic: zero convective HTC between the void fluid and the wall.
  //
  // Pressure inside the void is uniform — Darcy's law does not apply to the
  // empty volume. The face velocity into/out of the void is governed by the
  // adjacent bed cell's Darcy resistance. If the user wants to model real
  // manifold pressure drop (a screen, a fitting, packed support), they can
  // set `resistance_coefficient_` > 0 — that value is added in series with the
  // bed-side resistance at the void↔bed face.
  //
  // Units of resistance_coefficient_: Pa·s/m  (equivalent to µ·dx/K for an
  // additional in-series resistor). To target a known pressure drop ΔP at a
  // known superficial velocity u: resistance_coefficient_ = ΔP / u.
  struct VoidRegion {
    bool enabled_ = false;
    ExpressionParam volume_{1e-4};                  // [m^3]
    ExpressionParam resistance_coefficient_{0.0};   // [Pa·s/m] — optional drop

    void EvaluateExpressions(const EvalContext& ctx) {
      volume_.Evaluate(ctx);
      resistance_coefficient_.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::make_nvp("Enabled", enabled_),
        cereal::make_nvp("Volume", volume_)
      );
      if (version >= 2) {
        ar(cereal::make_nvp("Resistance_Coefficient", resistance_coefficient_));
      } else {
        // Legacy v1 saved files carry an unused Permeability_Override field;
        // accept and discard it so old cycles continue to load.
        ExpressionParam legacy_perm{0.0};
        ar(cereal::make_nvp("Permeability_Override", legacy_perm));
      }
    }
  };

} // end adsim namespace

CEREAL_CLASS_VERSION(adsim::VoidRegion, 2);
