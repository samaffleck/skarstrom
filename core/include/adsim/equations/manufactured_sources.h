/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// VERIFICATION HOOK -- not a modelling feature.
//
// Method of manufactured solutions (MMS): a verification test picks smooth
// fields P(z,t), y_i(z,t), q_i(z,t), u(z,t) (and T(z,t)), substitutes them
// into the continuous governing equations, and calls whatever does not
// vanish a source term. Adding exactly those sources to the discrete
// residuals makes the chosen fields the exact solution of the continuous
// problem, so the difference between the discrete solution and them is the
// discretisation error alone, and its observed order can be measured on a
// sequence of grids with every term of the model switched on at once.
//
// These sources are the only way such a test can reach the residual; they
// are deliberately NOT part of a deck: UserData holds them unserialised, and
// nothing in the parameter registry, the Python package or the web front end
// can set them. Empty (the default) means no source anywhere, and each
// residual then does exactly what it did without the hook (bit-identical
// results; the check per residual function is one branch on an empty
// std::function / map). See core/test/validation/test_manufactured_solutions.cpp.

// STL includes
#include <functional>
#include <string>
#include <unordered_map>


namespace adsim {

  struct ManufacturedSources {
    // A source as a function of the axial position z [m] (measured from the
    // bottom bed-end face, the coordinate of DomainVectors::positions) and the
    // solver time t [s] (the integrator's time, as passed to UpdateResiduals).
    using Field = std::function<double(double z, double t)>;

    // Each source is a term on the RIGHT-hand side of the continuous equation
    // named, in that equation's own units, evaluated at the cell centre
    // (balances) or the face (momentum) and added to the residual in the same
    // scaling as the neighbouring physical terms (see adsorbent_equations.cpp).

    // Overall mole balance, et dC/dt + d(u C)/dz + rho_b sum_i dq_i/dt = total
    // [mol/(m3 s)] (UpdatePressure).
    Field total;
    // Component balance of each gas component, et d(C y_i)/dt + d(u C y_i)/dz
    // - d(eb Dax C dy_i/dz)/dz + rho_b dq_i/dt = species[i] [mol/(m3 s)]
    // (UpdateMoleFraction). Keyed by component name.
    std::unordered_map<std::string, Field> species;
    // Rate law of each adsorbed component, dq_i/dt = rate(q_i*, q_i) +
    // loading[i] [mol/(kg s)] (UpdateAdsorbedConcentration), applied to the
    // loading of that component on every adsorbent. Keyed by component name.
    std::unordered_map<std::string, Field> loading;
    // Momentum, as a pressure-gradient (body-force) source [Pa/m]:
    //   -dP/dz = (mu/K) u + beta rho |u| u - momentum
    // (UpdateVelocity), at each face whose residual UpdateVelocity sets: every
    // Darcy face, and the interior faces of an Ergun layer. The bed-end faces
    // of an Ergun layer under a pressure, closed or flow condition impose the
    // Ergun law itself (boundary_conditions.cpp) and do not see this source,
    // so a manufactured solution for an Ergun layer must satisfy the momentum
    // law exactly at the two bed ends. Faces next to a void cell are not
    // covered either.
    Field momentum;
    // Bed energy balance in local thermal equilibrium [W/m3], the term added
    // to the right-hand side next to the heater's power (UpdateTemperature).
    // Not applied in LTNE (UpdateFluidTemperature) or to the wall.
    Field energy;

    bool Any() const {
      return static_cast<bool>(total) || !species.empty() || !loading.empty()
          || static_cast<bool>(momentum) || static_cast<bool>(energy);
    }
  };

} // end adsim namespace
