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
// Hard-coded gas-phase viscosity coefficients (GENERATED ONCE, offline).
//
// Each entry maps a component to the coefficients [a0..a4] of a degree-4
// polynomial mu_i(T) = a0 + a1*T + a2*T^2 + a3*T^3 + a4*T^4  [Pa-s], valid over
// T in [150, 650] K (see Fluid::kViscTmin/kViscTmax). Evaluating a polynomial
// at runtime costs ~1 ns; calling CoolProp per cell per residual costs ~20 us
// and would add minutes to a solve, so these are baked in rather than refit.
//
// Provenance (regenerate with .context/gen_visc_coeffs.cpp against CoolProp):
//   - For components CoolProp has a viscosity correlation for, the polynomial is
//     a least-squares fit of CoolProp's dilute-gas viscosity (sampled at low
//     pressure to stay in the vapour phase across the whole window).
//   - For the 8 components CoolProp has NO viscosity model for (CarbonMonoxide,
//     Ethylene, Fluorine, Krypton, Neon, NitrousOxide, SulfurDioxide, Xenon), the
//     polynomial fits Chapman-Enskog kinetic theory using the Lennard-Jones
//     eps/sigma from Components::componentData.
//   - Max relative fit error across the window is <=1% for every component.
//
// Gas viscosity is essentially pressure-independent over PSA/TSA swings (<~1%
// across 0.1-10 bar), so the fit is temperature-only. Unknown components fall
// back to a nominal dilute-gas value (see Fluid::GetComponentViscosity).
// ---------------------------------------------------------------------------

namespace adsim {

  // C++17 inline variable: one shared instance across all translation units.
  inline const std::unordered_map<std::string, std::array<double, 5>> kViscosityCoefficients = {
    {"Air", {{-3.975668308e-07, 8.355550932e-08, -8.946761945e-11, 8.105793652e-14, -3.291047444e-17}}}, // mu300=18.54 uPa.s
    {"Ammonia", {{5.662838731e-06, -2.243235441e-08, 1.897691980e-10, -2.536813623e-13, 1.230607119e-16}}}, // mu300=10.16
    {"Argon", {{-8.761247680e-07, 9.881118067e-08, -8.412117788e-11, 6.410174481e-14, -2.296012007e-17}}}, // mu300=22.74
    {"Benzene", {{1.254467318e-06, 1.638907665e-08, 2.301358987e-11, -2.659481940e-14, 1.084363444e-17}}}, // mu300=7.61
    {"Butane", {{1.156191447e-07, 2.297894162e-08, 1.056158257e-11, -2.219275149e-14, 1.244823019e-17}}}, // mu300=7.46
    {"CarbonDioxide", {{1.308709777e-07, 4.820694113e-08, 1.759231875e-11, -5.203688805e-14, 2.874773551e-17}}}, // mu300=15.00
    {"CarbonMonoxide", {{-6.092967425e-07, 8.487113566e-08, -1.080963739e-10, 1.092869458e-13, -4.764574635e-17}}}, // mu300=17.69 (Chapman-Enskog)
    {"Cyclohexane", {{4.473278869e-06, -9.547866920e-09, 8.972370936e-11, -1.123234660e-13, 4.951796980e-17}}}, // mu300=7.05
    {"Ethane", {{-3.428385448e-08, 3.184711463e-08, 6.571797844e-12, -3.240126156e-14, 2.075877142e-17}}}, // mu300=9.40
    {"Ethanol", {{-1.059293726e-06, 3.562551730e-08, -1.129025515e-11, 1.013812599e-14, -6.916014690e-18}}}, // mu300=8.83
    {"Ethylene", {{-3.009241108e-07, 3.471829047e-08, 1.382501935e-11, -5.314037293e-14, 3.534840835e-17}}}, // mu300=10.21 (Chapman-Enskog)
    {"Fluorine", {{-1.729438454e-06, 1.165001496e-07, -1.461342425e-10, 1.488435836e-13, -6.621532263e-17}}}, // mu300=23.55 (Chapman-Enskog)
    {"Helium", {{3.258489119e-06, 7.019132557e-08, -6.643826261e-11, 6.845952015e-14, -3.105348082e-17}}}, // mu300=19.93
    {"Hexane", {{-3.648210448e-07, 2.513173431e-08, -1.674071278e-11, 2.929517370e-14, -1.982747468e-17}}}, // mu300=6.30
    {"Hydrogen", {{1.063502253e-06, 3.469806338e-08, -3.932579965e-11, 4.337354026e-14, -2.025121124e-17}}}, // mu300=8.94
    {"HydrogenSulfide", {{2.570644577e-06, 7.734227358e-09, 1.359101931e-10, -2.202066452e-13, 1.187091121e-16}}}, // mu300=12.14
    {"Isobutane", {{-1.114114530e-07, 2.646426604e-08, -1.489597263e-12, -6.509845136e-15, 4.182420795e-18}}}, // mu300=7.55
    {"Isopentane", {{7.793850349e-07, 1.309729644e-08, 4.315222167e-11, -6.834421373e-14, 3.426868483e-17}}}, // mu300=7.02
    {"Krypton", {{-2.217551756e-06, 1.095455071e-07, -6.503966870e-11, 1.225257133e-14, 1.241648758e-17}}}, // mu300=25.22 (Chapman-Enskog)
    {"Methane", {{-1.113330524e-06, 5.500995567e-08, -6.266647772e-11, 6.466274551e-14, -2.931415673e-17}}}, // mu300=11.26
    {"Methanol", {{8.255187980e-07, 2.844654012e-08, -3.736111742e-12, 3.362927905e-14, -3.085984234e-17}}}, // mu300=9.68
    {"Neon", {{3.350023905e-06, 1.316628940e-07, -1.838039767e-10, 2.141223765e-13, -1.035434873e-16}}}, // mu300=31.25 (Chapman-Enskog)
    {"Nitrogen", {{-2.489909434e-07, 8.016074213e-08, -8.622093890e-11, 7.807414792e-14, -3.158186321e-17}}}, // mu300=17.89
    {"NitrousOxide", {{-2.805803269e-07, 4.846025334e-08, 2.851470960e-11, -8.857142835e-14, 5.687294341e-17}}}, // mu300=14.89 (Chapman-Enskog)
    {"Octane", {{-1.079308585e-06, 3.187610322e-08, -4.949795456e-11, 7.746855436e-14, -4.370967753e-17}}}, // mu300=5.77
    {"Oxygen", {{-6.424959437e-07, 9.219506868e-08, -9.156634332e-11, 7.888483229e-14, -3.102898027e-17}}}, // mu300=20.65
    {"Pentane", {{-3.193560173e-07, 2.589662993e-08, -1.395470524e-11, 2.753700815e-14, -2.045304206e-17}}}, // mu300=6.77
    {"Propane", {{9.477457721e-07, 1.723091423e-08, 4.052645690e-11, -6.898194921e-14, 3.552897634e-17}}}, // mu300=8.19
    {"Propylene", {{1.040721774e-06, 1.361304903e-08, 6.932026272e-11, -1.227180457e-13, 6.835768150e-17}}}, // mu300=8.60
    {"SulfurDioxide", {{1.081053654e-06, 2.763383941e-08, 6.942644824e-11, -1.162432677e-13, 5.962990370e-17}}}, // mu300=12.96 (Chapman-Enskog)
    {"SulfurHexafluoride", {{8.115148060e-07, 4.938162982e-08, 3.676082818e-12, -2.907000788e-14, 1.630163748e-17}}}, // mu300=15.30
    {"Toluene", {{3.460682436e-07, 2.562857369e-08, -2.334484737e-11, 5.045965298e-14, -3.465751228e-17}}}, // mu300=7.02
    {"Water", {{6.685582872e-06, -2.114915253e-08, 1.421511602e-10, -1.395610837e-13, 4.941363610e-17}}}, // mu300=9.77
    {"Xenon", {{-4.791230607e-07, 7.569624049e-08, 4.190674025e-11, -1.342713884e-13, 8.671853481e-17}}}, // mu300=23.08 (Chapman-Enskog)
  };

} // namespace adsim
