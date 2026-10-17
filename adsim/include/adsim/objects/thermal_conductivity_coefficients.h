/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

#include <array>
#include <string>
#include <unordered_map>

// ---------------------------------------------------------------------------
// Hard-coded gas-phase THERMAL CONDUCTIVITY coefficients (GENERATED ONCE, offline).
//
// Each entry maps a component to the coefficients [a0..a4] of a degree-4
// polynomial k_i(T) = a0 + a1*T + a2*T^2 + a3*T^3 + a4*T^4  [W/(m.K)], valid over
// T in [150, 650] K (Fluid::kViscTmin/kViscTmax). Evaluating a polynomial at
// runtime costs ~1 ns; calling CoolProp per cell per residual costs ~20 us and
// would add minutes to a solve, so these are baked in rather than refit.
//
// Provenance (regenerate with .context/gen_cond_coeffs.cpp against CoolProp):
//   - For components CoolProp has a conductivity correlation for, the polynomial
//     is a least-squares fit of CoolProp's dilute-gas conductivity (sampled at
//     low pressure to stay in the vapour phase across the whole window).
//   - For the 10 components CoolProp has NO conductivity model for (returns inf),
//     the polynomial fits the Eucken relation k = mu*(cp_mass + 1.25*R/M), reusing
//     the pure-gas viscosity (CoolProp or Chapman-Enskog) and CoolProp Cp.
//   - Max relative fit error across the window is <=1.8% for every component.
//
// Gas conductivity is weakly pressure-dependent over PSA/TSA swings, so the fit
// is temperature-only. Unknown components fall back to a nominal value (see
// Fluid::GetComponentThermalConductivity).
// ---------------------------------------------------------------------------

namespace adsim {

  // C++17 inline variable: one shared instance across all translation units.
  inline const std::unordered_map<std::string, std::array<double, 5>> kThermalConductivityCoefficients = {
    {"Air", {{-1.103689343e-03, 1.150872000e-04, -1.017278775e-07, 8.892823023e-11, -3.489346002e-14}}}, // k300=26.39 mW/m.K
    {"Ammonia", {{3.364189134e-02, -1.526062879e-04, 3.784640531e-07, 2.804148391e-10, -5.422946728e-13}}}, // k300=25.10 mW/m.K
    {"Argon", {{-7.125512209e-04, 7.752521344e-05, -6.538144016e-08, 4.867387913e-11, -1.701578264e-14}}}, // k300=17.84 mW/m.K
    {"Benzene", {{1.298818600e-02, -1.453267222e-04, 6.602809907e-07, -7.751658646e-10, 3.470083372e-13}}}, // k300=10.70 mW/m.K
    {"Butane", {{2.171007845e-03, -3.984701616e-06, 1.842877141e-07, -3.543541947e-11, 1.860946374e-14}}}, // k300=16.76 mW/m.K
    {"CarbonDioxide", {{9.416515217e-04, 1.366165474e-05, 1.900929996e-07, -2.284300458e-10, 9.819116439e-14}}}, // k300=16.78 mW/m.K
    {"CarbonMonoxide", {{-1.376844404e-03, 1.273088543e-04, -1.909338667e-07, 2.292556706e-10, -1.062151930e-13}}}, // k300=24.96 mW/m.K (Eucken)
    {"Cyclohexane", {{7.296756009e-03, -6.238471551e-05, 2.925499196e-07, -1.984080354e-10, 3.591056476e-14}}}, // k300=9.84 mW/m.K (Eucken)
    {"Ethane", {{1.284237111e-03, 9.409162784e-06, 1.599536977e-07, 1.661591530e-10, -2.083007400e-13}}}, // k300=21.30 mW/m.K
    {"Ethanol", {{1.183679899e-02, -1.147692558e-04, 5.990535289e-07, -6.806635415e-10, 3.033441519e-13}}}, // k300=15.40 mW/m.K
    {"Ethylene", {{5.728143734e-03, -2.688768639e-05, 3.246114713e-07, -3.030534938e-10, 9.855228207e-14}}}, // k300=19.49 mW/m.K (Eucken)
    {"Fluorine", {{3.873825544e-04, 8.816101660e-05, 1.156439994e-09, -4.807312824e-11, 2.769617365e-14}}}, // k300=25.87 mW/m.K (Eucken)
    {"Helium", {{2.277101233e-02, 5.654886831e-04, -5.498894004e-07, 5.592509200e-10, -2.506192097e-13}}}, // k300=156.00 mW/m.K
    {"Hexane", {{5.356463592e-03, -2.977925585e-05, 2.043950646e-07, -9.003284588e-12, -7.953449793e-14}}}, // k300=13.93 mW/m.K
    {"Hydrogen", {{-2.815427557e-02, 1.071286869e-03, -1.651895725e-06, 1.779587551e-09, -7.351566094e-13}}}, // k300=186.66 mW/m.K
    {"HydrogenSulfide", {{4.135948337e-03, 2.030482065e-06, 1.898541429e-07, -2.618401637e-10, 1.357468661e-13}}}, // k300=15.86 mW/m.K (Eucken)
    {"Isobutane", {{-1.916044286e-03, 2.117536815e-05, 1.457239015e-07, -1.832318819e-11, 6.720988992e-15}}}, // k300=17.11 mW/m.K
    {"Isopentane", {{1.780281893e-03, -2.285042130e-05, 2.226711585e-07, 9.485397529e-12, -7.823995723e-14}}}, // k300=14.59 mW/m.K
    {"Krypton", {{-8.198950855e-04, 4.071440023e-05, -2.404924216e-08, 4.336560210e-12, 4.740714021e-15}}}, // k300=9.39 mW/m.K (Eucken)
    {"Methane", {{-1.342352025e-03, 1.309864327e-04, -1.853001437e-07, 6.135535740e-10, -4.171245103e-13}}}, // k300=34.46 mW/m.K
    {"Methanol", {{7.149426974e-03, -4.860312275e-05, 3.068597586e-07, -1.962327540e-10, 5.129961260e-14}}}, // k300=15.30 mW/m.K
    {"Neon", {{5.177034927e-03, 2.034306502e-04, -2.839832274e-07, 3.308219126e-10, -1.599750404e-13}}}, // k300=48.28 mW/m.K (Eucken)
    {"Nitrogen", {{-1.140235134e-03, 1.150931704e-04, -1.074526406e-07, 9.459915977e-11, -3.721974272e-14}}}, // k300=25.97 mW/m.K
    {"NitrousOxide", {{-1.073241639e-03, 3.959371542e-05, 1.047639242e-07, -1.580727671e-10, 7.935528197e-14}}}, // k300=16.61 mW/m.K (Eucken)
    {"Octane", {{1.514313723e-02, -1.412034281e-04, 5.769547248e-07, -5.727057956e-10, 2.249106792e-13}}}, // k300=11.07 mW/m.K
    {"Oxygen", {{-1.337308378e-03, 1.110772666e-04, -7.704117649e-08, 5.946003230e-11, -2.125278802e-14}}}, // k300=26.49 mW/m.K
    {"Pentane", {{6.083391917e-03, -5.726343072e-05, 3.451353097e-07, -2.059411437e-10, 3.085914190e-14}}}, // k300=14.66 mW/m.K
    {"Propane", {{-6.332823724e-04, 1.426336110e-05, 1.779276446e-07, -5.130067726e-11, 2.878203637e-14}}}, // k300=18.51 mW/m.K
    {"Propylene", {{5.618233006e-03, -4.293905616e-05, 3.428177340e-07, -2.533491288e-10, 7.086168179e-14}}}, // k300=17.32 mW/m.K
    {"SulfurDioxide", {{1.656943128e-03, 6.152581683e-06, 1.071122546e-07, -1.237629746e-10, 4.839014155e-14}}}, // k300=10.19 mW/m.K (Eucken)
    {"SulfurHexafluoride", {{-2.893111203e-03, 9.567804820e-06, 2.379813880e-07, -3.598758976e-10, 1.800158661e-13}}}, // k300=13.14 mW/m.K
    {"Toluene", {{2.437868018e-03, -2.370139831e-05, 1.914645592e-07, -5.045649611e-11, -3.937756297e-14}}}, // k300=10.88 mW/m.K
    {"Water", {{1.992543704e-02, -1.304602429e-04, 6.227676609e-07, -8.012497117e-10, 4.136292535e-13}}}, // k300=18.55 mW/m.K
    {"Xenon", {{-2.204179156e-04, 1.914318761e-05, 5.478124120e-09, -2.466753816e-11, 1.641775460e-14}}}, // k300=5.48 mW/m.K (Eucken)
  };

} // namespace adsim
