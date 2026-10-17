/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

namespace adsim {

  // Bulk transport properties the closure models (mass/heat transfer
  // coefficients, dispersion) read from the state. Each model reports the ones
  // it reads via ReadsTransportProperties(); Equations::Initialise computes a
  // property no selected model reads only when the state is logged, not on
  // every residual evaluation. The interfaces default to all of them, so a
  // model that does not override it is always given up-to-date values.
  enum TransportProperty : unsigned {
    kNoTransportProperties  = 0,
    kMolecularDiffusion     = 1u << 0,
    kSchmidtNumber          = 1u << 1,
    kSherwoodNumber         = 1u << 2,
    kReynoldsNumber         = 1u << 3,
    kGasThermalConductivity = 1u << 4,
    kAllTransportProperties = (1u << 5) - 1,
  };

} // namespace adsim
