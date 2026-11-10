/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

#include <array>
#include <string>
#include <string_view>
#include <unordered_map>

// ---------------------------------------------------------------------------
// Ideal-gas heat capacities C_p0(T): NASA 7-coefficient polynomials.
//
//   C_p0 / R = a0 + a1 T + a2 T^2 + a3 T^3 + a4 T^4
//   H0 / (R T) = a0 + a1 T / 2 + a2 T^2 / 3 + a3 T^3 / 4 + a4 T^4 / 5 + a5 / T
//
// with ONE coefficient set per species, for t_min <= T <= t_max = 1000 K: the
// low-temperature set of TM-4513. Adsorption processes stay well below
// 1000 K, so the report's second set (1000-6000 K) is not carried. Source:
//
//   B. J. McBride, S. Gordon and M. A. Reno, "Coefficients for Calculating
//   Thermodynamic and Transport Properties of Individual Species", NASA
//   Technical Memorandum 4513, October 1993 (NTRS 19940013151), species
//   tables of the gas-phase database.
//
// Each entry names its TM-4513 species and the data source code TM-4513 gives
// for it (the report's reference list resolves the codes: TPIS = Gurvich et
// al., Thermodynamic Properties of Individual Substances; J = JANAF tables;
// L = evaluated by the NASA Lewis group; X = other, see the report). The
// coefficients are copied unchanged from the electronic form of the database
// distributed with Cantera (data/nasa_gas.yaml, which states TM-4513 as its
// source), and were cross-checked digit by digit against the printed TM-4513
// tables: a0-a5, the ranges and the species agree; the entropy constant a6
// of F2, He, n-pentane and i-pentane differed in the electronic copy and is
// given here as printed (a6 enters neither C_p0 nor H0). Argon, helium,
// krypton and neon have one 200-6000 K range in TM-4513 (a constant 5/2 R),
// cut at 1000 K like the others. Air, absent from TM-4513, is the ideal-gas mixture of the N2, O2
// and Ar entries at the composition of Lemmon, Jacobsen, Penoncello and Friend,
// J. Phys. Chem. Ref. Data 29 (2000) 331, the air model of CoolProp; a
// mole-fraction mixture of ideal gases mixes the coefficients linearly, so the
// mixed polynomial is exact.
//
// Agreement with CoolProp's ideal-gas C_p0 (the reference equations of state)
// is checked in test/unit/test_ideal_gas_heat_capacity.cpp (against values
// recorded from CoolProp) and tabulated in docs/equation-of-state.md. Every
// component of Components::componentNames has an entry.
// ---------------------------------------------------------------------------

namespace adsim {

  struct NasaPolynomial {
    double t_min = 0.0;  // [K] lower limit of the fit
    double t_max = 0.0;  // [K] upper limit of the fit (1000 K)
    std::array<double, 7> a{};  // a0..a6 for t_min <= T <= t_max
    std::string_view species{};    // TM-4513 species label (Air, absent from TM-4513: its mixture)
    std::string_view source{};     // TM-4513 data-source code, or the note for a mixed entry
  };

  // Reference of the table; Air's `source` note gives its mixing reference.
  inline constexpr std::string_view kIdealGasHeatCapacityReference =
      "B. J. McBride, S. Gordon and M. A. Reno, Coefficients for Calculating Thermodynamic and Transport "
      "Properties of Individual Species, NASA Technical Memorandum 4513 (1993)";

  // C++17 inline variable: one shared instance across all translation units.
  inline const std::unordered_map<std::string, NasaPolynomial> kIdealGasHeatCapacity = {
    // N2 0.7812 + O2 0.2096 + Ar 0.0092 (mole fractions of Lemmon et al. 2000), mixed from the N2, O2 and Ar entries
    {"Air", {200, 1000,
      {3.57422417779, -0.000724719440884, 1.67105133902e-06, -1.26738307824e-10, -4.20678743564e-13, -1047.75789011, 3.1251330622}, "N2 0.7812 + O2 0.2096 + Ar 0.0092", "mixed from N2/O2/Ar; Lemmon 2000"}},
    // TM-4513 species NH3, data source TPIS89
    {"Ammonia", {200, 1000,
      {4.30177808, -0.0047712733, 2.19341619e-05, -2.29856489e-08, 8.28992268e-12, -6748.06394, -0.690644393}, "NH3", "TPIS89"}},
    // TM-4513 species Ar, data source L 6/88
    {"Argon", {200, 1000,
      {2.5, 0, 0, 0, 0, -745.375, 4.37967491}, "Ar", "L 6/88"}},
    // TM-4513 species C6H6, data source L 1/91
    {"Benzene", {200, 1000,
      {0.503469664, 0.0185142363, 7.37864409e-05, -1.18106127e-07, 5.07182527e-11, 8552.66293, 21.6481796}, "C6H6", "L 1/91"}},
    // TM-4513 species C4H10,n-butane, data source L 6/90
    {"Butane", {200, 1000,
      {6.14746806, 0.000155947389, 9.67913517e-05, -1.2548391e-07, 4.97816555e-11, -17599.4402, -1.09409879}, "C4H10,n-butane", "L 6/90"}},
    // TM-4513 species CO2, data source L 7/88
    {"CarbonDioxide", {200, 1000,
      {2.35677352, 0.00898459677, -7.12356269e-06, 2.45919022e-09, -1.43699548e-13, -48371.9697, 9.90105222}, "CO2", "L 7/88"}},
    // TM-4513 species CO, data source TPIS79
    {"CarbonMonoxide", {200, 1000,
      {3.57953347, -0.00061035368, 1.01681433e-06, 9.07005884e-10, -9.04424499e-13, -14344.086, 3.50840928}, "CO", "TPIS79"}},
    // TM-4513 species C6H12,cyclo-, data source L 6/90
    {"Cyclohexane", {200, 1000,
      {4.04348764, -0.00619527424, 0.000176621086, -2.22967809e-07, 8.6366739e-11, -16920.2872, 8.52566766}, "C6H12,cyclo-", "L 6/90"}},
    // TM-4513 species C2H6, data source L 8/88
    {"Ethane", {200, 1000,
      {4.29142492, -0.0055015427, 5.99438288e-05, -7.08466285e-08, 2.68685771e-11, -11522.2055, 2.66682316}, "C2H6", "L 8/88"}},
    // TM-4513 species C2H5OH, data source L 8/88
    {"Ethanol", {200, 1000,
      {4.85868178, -0.0037400674, 6.95550267e-05, -8.86541147e-08, 3.5168443e-11, -29996.1309, 4.80192294}, "C2H5OH", "L 8/88"}},
    // TM-4513 species C2H4, data source L 1/91
    {"Ethylene", {200, 1000,
      {3.95920148, -0.00757052247, 5.70990292e-05, -6.91588753e-08, 2.69884373e-11, 5089.77593, 4.09733096}, "C2H4", "L 1/91"}},
    // TM-4513 species F2, data source TPIS89
    {"Fluorine", {200, 1000,
      {3.20832415, 0.00125919179, 3.89747979e-06, -7.22184984e-09, 3.31837862e-12, -1034.25794, 5.61903603}, "F2", "TPIS89"}},
    // TM-4513 species He, data source L10/90
    {"Helium", {200, 1000,
      {2.5, 0, 0, 0, 0, -745.375, 0.928723974}, "He", "L10/90"}},
    // TM-4513 species H2, data source TPIS78
    {"Hydrogen", {200, 1000,
      {2.34433112, 0.00798052075, -1.9478151e-05, 2.01572094e-08, -7.37611761e-12, -917.935173, 0.683010238}, "H2", "TPIS78"}},
    // TM-4513 species H2S, data source J 6/77
    {"HydrogenSulfide", {300, 1000,
      {3.9323476, -0.00050260905, 4.5928473e-06, -3.1807214e-09, 6.6497561e-13, -3650.5359, 2.3157905}, "H2S", "J 6/77"}},
    // TM-4513 species C4H10,isobutane, data source L 6/90
    {"Isobutane", {200, 1000,
      {4.45479276, 0.00826057985, 8.29886664e-05, -1.14647642e-07, 4.64570101e-11, -18459.3931, 4.92743175}, "C4H10,isobutane", "L 6/90"}},
    // TM-4513 species C5H12,i-pentane, data source X10/85
    {"Isopentane", {298.15, 1000,
      {1.0832882, 0.044571076, 8.2389934e-06, -3.5258047e-08, 1.5785762e-11, -20807.535, 21.795145}, "C5H12,i-pentane", "X10/85"}},
    // TM-4513 species Kr, data source L10/90
    {"Krypton", {200, 1000,
      {2.5, 0, 0, 0, 0, -745.375, 5.49095651}, "Kr", "L10/90"}},
    // TM-4513 species CH4, data source L 8/88
    {"Methane", {200, 1000,
      {5.14987613, -0.0136709788, 4.91800599e-05, -4.84743026e-08, 1.66693956e-11, -10246.6476, -4.64130376}, "CH4", "L 8/88"}},
    // TM-4513 species CH3OH, data source L 8/88
    {"Methanol", {200, 1000,
      {5.71539582, -0.0152309129, 6.52441155e-05, -7.10806889e-08, 2.61352698e-11, -25642.7656, -1.50409823}, "CH3OH", "L 8/88"}},
    // TM-4513 species Ne, data source L10/90
    {"Neon", {200, 1000,
      {2.5, 0, 0, 0, 0, -745.375, 3.35532272}, "Ne", "L10/90"}},
    // TM-4513 species N2, data source TPIS78
    {"Nitrogen", {200, 1000,
      {3.53100528, -0.000123660987, -5.02999437e-07, 2.43530612e-09, -1.40881235e-12, -1046.97628, 2.96747468}, "N2", "TPIS78"}},
    // TM-4513 species N2O, data source L 7/88
    {"NitrousOxide", {200, 1000,
      {2.2571686, 0.0113046338, -1.3671035e-05, 9.68162098e-09, -2.93055583e-12, 8741.77146, 10.7579154}, "N2O", "L 7/88"}},
    // TM-4513 species C8H18,n-octane, data source X 4/85
    {"Octane", {200, 1000,
      {12.5244908, -0.0101018365, 0.000221991595, -2.8486242e-07, 1.12409624e-10, -29843.3034, -19.7108554}, "C8H18,n-octane", "X 4/85"}},
    // TM-4513 species O2, data source TPIS89
    {"Oxygen", {200, 1000,
      {3.78245636, -0.00299673415, 9.847302e-06, -9.68129508e-09, 3.24372836e-12, -1063.94356, 3.65767573}, "O2", "TPIS89"}},
    // TM-4513 species C5H12,n-pentane, data source X10/85
    {"Pentane", {298.15, 1000,
      {1.8983679, 0.041203037, 1.2312175e-05, -3.6589501e-08, 1.5042509e-11, -20091.5, 18.679072}, "C5H12,n-pentane", "X10/85"}},
    // TM-4513 species C3H8, data source L 6/90
    {"Propane", {200, 1000,
      {4.2110262, 0.00171599803, 7.06183472e-05, -9.19594116e-08, 3.64421372e-11, -14381.2106, 5.60930491}, "C3H8", "L 6/90"}},
    // TM-4513 species C3H6,propylene, data source L 7/90
    {"Propylene", {200, 1000,
      {3.83464524, 0.00329078405, 5.05228184e-05, -6.66251418e-08, 2.63707585e-11, 753.838295, 7.53410995}, "C3H6,propylene", "L 7/90"}},
    // TM-4513 species SO2, data source J 6/61
    {"SulfurDioxide", {300, 1000,
      {3.2665338, 0.0053237902, 6.8437552e-07, -5.2810047e-09, 2.5590454e-12, -36908.148, 9.66465108}, "SO2", "J 6/61"}},
    // TM-4513 species SF6, data source J 6/76
    {"SulfurHexafluoride", {300, 1000,
      {-3.8388088, 0.083221721, -0.00013181689, 9.9636154e-08, -2.9248767e-11, -148364.77, 37.1611426}, "SF6", "J 6/76"}},
    // TM-4513 species C7H8, data source L 1/93
    {"Toluene", {200, 1000,
      {1.611914, 0.0211188902, 8.53221453e-05, -1.32566876e-07, 5.59406109e-11, 4096.51976, 20.2973614}, "C7H8", "L 1/93"}},
    // TM-4513 species H2O, data source L 8/89
    {"Water", {200, 1000,
      {4.19864056, -0.0020364341, 6.52040211e-06, -5.48797062e-09, 1.77197817e-12, -30293.7267, -0.849032208}, "H2O", "L 8/89"}},
    // TM-4513 species Xe, data source L12/91
    {"Xenon", {200, 1000,
      {2.5, -8.9914133e-14, 2.5219686e-16, -2.92186662e-19, 1.18949218e-22, -745.375, 6.16441993}, "Xe", "L12/91"}},

  };

}  // namespace adsim
