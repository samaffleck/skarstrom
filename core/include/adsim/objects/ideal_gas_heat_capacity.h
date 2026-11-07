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
// with one coefficient set for t_min <= T < t_mid and another for
// t_mid <= T <= t_max (a single-range species repeats its set). Source:
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
// krypton and neon have one range in TM-4513: t_mid = t_max = 6000 K and both
// sets equal. Air, absent from TM-4513, is the ideal-gas mixture of the N2, O2
// and Ar entries at the composition of Lemmon, Jacobsen, Penoncello and Friend,
// J. Phys. Chem. Ref. Data 29 (2000) 331, the air model of CoolProp; a
// mole-fraction mixture of ideal gases mixes the coefficients linearly, so the
// mixed polynomial is exact.
//
// Agreement with CoolProp's ideal-gas C_p0 (the reference equations of state)
// is checked in test/unit/test_ideal_gas_heat_capacity.cpp and tabulated in
// docs/equation-of-state.md. Components of Components::componentNames without
// an entry (Hexane) keep the CoolProp value held at the initial conditions
// (Fluid::Initialise), and the run-start warnings say so.
// ---------------------------------------------------------------------------

namespace adsim {

  struct NasaPolynomial {
    double t_min = 0.0;  // [K] lower limit of the fit
    double t_mid = 0.0;  // [K] switch from `low` to `high`
    double t_max = 0.0;  // [K] upper limit of the fit
    std::array<double, 7> low{};   // a0..a6 for t_min <= T < t_mid
    std::array<double, 7> high{};  // a0..a6 for t_mid <= T <= t_max
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
    {"Air", {200, 1000, 6000,
      {3.57422417779, -0.000724719440884, 1.67105133902e-06, -1.26738307824e-10, -4.20678743564e-13, -1047.75789011, 3.1251330622},
      {3.09688996428, 0.0012288329389, -4.14428809065e-07, 6.57166487817e-11, -3.87171873546e-15, -983.514963074, 5.34327528746}, "N2 0.7812 + O2 0.2096 + Ar 0.0092", "mixed from N2/O2/Ar; Lemmon 2000"}},
    // TM-4513 species NH3, data source TPIS89
    {"Ammonia", {200, 1000, 6000,
      {4.30177808, -0.0047712733, 2.19341619e-05, -2.29856489e-08, 8.28992268e-12, -6748.06394, -0.690644393},
      {2.71709692, 0.00556856338, -1.76886396e-06, 2.6741726e-10, -1.52731419e-14, -6584.51989, 6.09289837}, "NH3", "TPIS89"}},
    // TM-4513 species Ar, data source L 6/88
    {"Argon", {200, 6000, 6000,
      {2.5, 0, 0, 0, 0, -745.375, 4.37967491},
      {2.5, 0, 0, 0, 0, -745.375, 4.37967491}, "Ar", "L 6/88"}},
    // TM-4513 species C6H6, data source L 1/91
    {"Benzene", {200, 1000, 6000,
      {0.503469664, 0.0185142363, 7.37864409e-05, -1.18106127e-07, 5.07182527e-11, 8552.66293, 21.6481796},
      {11.0771708, 0.0207067895, -7.516251e-06, 1.22209416e-09, -7.35312513e-14, 4309.88395, -40.011695}, "C6H6", "L 1/91"}},
    // TM-4513 species C4H10,n-butane, data source L 6/90
    {"Butane", {200, 1000, 6000,
      {6.14746806, 0.000155947389, 9.67913517e-05, -1.2548391e-07, 4.97816555e-11, -17599.4402, -1.09409879},
      {9.44535834, 0.0257858073, -9.23619122e-06, 1.48632755e-09, -8.87897158e-14, -20138.2165, -26.3470076}, "C4H10,n-butane", "L 6/90"}},
    // TM-4513 species CO2, data source L 7/88
    {"CarbonDioxide", {200, 1000, 6000,
      {2.35677352, 0.00898459677, -7.12356269e-06, 2.45919022e-09, -1.43699548e-13, -48371.9697, 9.90105222},
      {4.63659493, 0.00274131991, -9.95828531e-07, 1.60373011e-10, -9.16103468e-15, -49024.9341, -1.93534855}, "CO2", "L 7/88"}},
    // TM-4513 species CO, data source TPIS79
    {"CarbonMonoxide", {200, 1000, 6000,
      {3.57953347, -0.00061035368, 1.01681433e-06, 9.07005884e-10, -9.04424499e-13, -14344.086, 3.50840928},
      {3.04848583, 0.00135172818, -4.85794075e-07, 7.88536486e-11, -4.69807489e-15, -14266.1171, 6.0170979}, "CO", "TPIS79"}},
    // TM-4513 species C6H12,cyclo-, data source L 6/90
    {"Cyclohexane", {200, 1000, 6000,
      {4.04348764, -0.00619527424, 0.000176621086, -2.22967809e-07, 8.6366739e-11, -16920.2872, 8.52566766},
      {13.2147562, 0.035824241, -1.32110595e-05, 2.17202254e-09, -1.3173054e-13, -22809.1954, -55.3526464}, "C6H12,cyclo-", "L 6/90"}},
    // TM-4513 species C2H6, data source L 8/88
    {"Ethane", {200, 1000, 6000,
      {4.29142492, -0.0055015427, 5.99438288e-05, -7.08466285e-08, 2.68685771e-11, -11522.2055, 2.66682316},
      {4.04666674, 0.0153538766, -5.47039321e-06, 8.77826228e-10, -5.23167305e-14, -12447.3512, -0.968683607}, "C2H6", "L 8/88"}},
    // TM-4513 species C2H5OH, data source L 8/88
    {"Ethanol", {200, 1000, 6000,
      {4.85868178, -0.0037400674, 6.95550267e-05, -8.86541147e-08, 3.5168443e-11, -29996.1309, 4.80192294},
      {6.5628977, 0.0152034264, -5.38922247e-06, 8.62150224e-10, -5.12824683e-14, -31525.7984, -9.47557644}, "C2H5OH", "L 8/88"}},
    // TM-4513 species C2H4, data source L 1/91
    {"Ethylene", {200, 1000, 6000,
      {3.95920148, -0.00757052247, 5.70990292e-05, -6.91588753e-08, 2.69884373e-11, 5089.77593, 4.09733096},
      {3.99182761, 0.010483391, -3.71721385e-06, 5.94628514e-10, -3.53630526e-14, 4268.65819, -0.269052151}, "C2H4", "L 1/91"}},
    // TM-4513 species F2, data source TPIS89
    {"Fluorine", {200, 1000, 6000,
      {3.20832415, 0.00125919179, 3.89747979e-06, -7.22184984e-09, 3.31837862e-12, -1034.25794, 5.61903603},
      {3.86166219, 0.000788367679, -1.8198294e-07, -9.1743656e-12, 2.65193472e-15, -1232.38655, 2.04119869}, "F2", "TPIS89"}},
    // TM-4513 species He, data source L10/90
    {"Helium", {200, 6000, 6000,
      {2.5, 0, 0, 0, 0, -745.375, 0.928723974},
      {2.5, 0, 0, 0, 0, -745.375, 0.928723974}, "He", "L10/90"}},
    // TM-4513 species H2, data source TPIS78
    {"Hydrogen", {200, 1000, 6000,
      {2.34433112, 0.00798052075, -1.9478151e-05, 2.01572094e-08, -7.37611761e-12, -917.935173, 0.683010238},
      {2.93286579, 0.000826607967, -1.46402335e-07, 1.54100359e-11, -6.88804432e-16, -813.065597, -1.02432887}, "H2", "TPIS78"}},
    // TM-4513 species H2S, data source J 6/77
    {"HydrogenSulfide", {300, 1000, 5000,
      {3.9323476, -0.00050260905, 4.5928473e-06, -3.1807214e-09, 6.6497561e-13, -3650.5359, 2.3157905},
      {2.7452199, 0.0040434607, -1.538451e-06, 2.7520249e-10, -1.8592095e-14, -3419.9444, 8.0546745}, "H2S", "J 6/77"}},
    // TM-4513 species C4H10,isobutane, data source L 6/90
    {"Isobutane", {200, 1000, 6000,
      {4.45479276, 0.00826057985, 8.29886664e-05, -1.14647642e-07, 4.64570101e-11, -18459.3931, 4.92743175},
      {9.76991245, 0.025499721, -9.14142932e-06, 1.47328271e-09, -8.80800188e-14, -21405.2647, -30.0329101}, "C4H10,isobutane", "L 6/90"}},
    // TM-4513 species C5H12,i-pentane, data source X10/85
    {"Isopentane", {298.15, 1000, 5000,
      {1.0832882, 0.044571076, 8.2389934e-06, -3.5258047e-08, 1.5785762e-11, -20807.535, 21.795145},
      {12.327787, 0.030613087, -9.8415785e-06, 1.3919776e-09, -7.0337345e-14, -25037.492, -41.133504}, "C5H12,i-pentane", "X10/85"}},
    // TM-4513 species Kr, data source L10/90
    {"Krypton", {200, 6000, 6000,
      {2.5, 0, 0, 0, 0, -745.375, 5.49095651},
      {2.5, 0, 0, 0, 0, -745.375, 5.49095651}, "Kr", "L10/90"}},
    // TM-4513 species CH4, data source L 8/88
    {"Methane", {200, 1000, 6000,
      {5.14987613, -0.0136709788, 4.91800599e-05, -4.84743026e-08, 1.66693956e-11, -10246.6476, -4.64130376},
      {1.63552643, 0.0100842795, -3.36916254e-06, 5.34958667e-10, -3.15518833e-14, -10005.6455, 9.99313326}, "CH4", "L 8/88"}},
    // TM-4513 species CH3OH, data source L 8/88
    {"Methanol", {200, 1000, 6000,
      {5.71539582, -0.0152309129, 6.52441155e-05, -7.10806889e-08, 2.61352698e-11, -25642.7656, -1.50409823},
      {3.60134486, 0.0102430954, -3.59985517e-06, 5.72505986e-10, -3.3911764e-14, -25997.191, 4.70512253}, "CH3OH", "L 8/88"}},
    // TM-4513 species Ne, data source L10/90
    {"Neon", {200, 6000, 6000,
      {2.5, 0, 0, 0, 0, -745.375, 3.35532272},
      {2.5, 0, 0, 0, 0, -745.375, 3.35532272}, "Ne", "L10/90"}},
    // TM-4513 species N2, data source TPIS78
    {"Nitrogen", {200, 1000, 6000,
      {3.53100528, -0.000123660987, -5.02999437e-07, 2.43530612e-09, -1.40881235e-12, -1046.97628, 2.96747468},
      {2.95257626, 0.00139690057, -4.92631691e-07, 7.86010367e-11, -4.60755321e-15, -923.948645, 5.87189252}, "N2", "TPIS78"}},
    // TM-4513 species N2O, data source L 7/88
    {"NitrousOxide", {200, 1000, 6000,
      {2.2571686, 0.0113046338, -1.3671035e-05, 9.68162098e-09, -2.93055583e-12, 8741.77146, 10.7579154},
      {4.82318873, 0.00262685279, -9.58426058e-07, 1.59991296e-10, -9.77416939e-15, 8073.35662, -2.202366}, "N2O", "L 7/88"}},
    // TM-4513 species C8H18,n-octane, data source X 4/85
    {"Octane", {200, 1000, 6000,
      {12.5244908, -0.0101018365, 0.000221991595, -2.8486242e-07, 1.12409624e-10, -29843.3034, -19.7108554},
      {22.1755407, 0.0424426161, -1.49161103e-05, 2.40376673e-09, -1.44359037e-13, -36103.0944, -88.0854457}, "C8H18,n-octane", "X 4/85"}},
    // TM-4513 species O2, data source TPIS89
    {"Oxygen", {200, 1000, 6000,
      {3.78245636, -0.00299673415, 9.847302e-06, -9.68129508e-09, 3.24372836e-12, -1063.94356, 3.65767573},
      {3.66096083, 0.000656365523, -1.41149485e-07, 2.05797658e-11, -1.29913248e-15, -1215.97725, 3.41536184}, "O2", "TPIS89"}},
    // TM-4513 species C5H12,n-pentane, data source X10/85
    {"Pentane", {298.15, 1000, 5000,
      {1.8983679, 0.041203037, 1.2312175e-05, -3.6589501e-08, 1.5042509e-11, -20091.5, 18.679072},
      {13.546998, 0.028421786, -9.4174648e-06, 1.3893589e-09, -7.4212609e-14, -24577.68, -47.021185}, "C5H12,n-pentane", "X10/85"}},
    // TM-4513 species C3H8, data source L 6/90
    {"Propane", {200, 1000, 6000,
      {4.2110262, 0.00171599803, 7.06183472e-05, -9.19594116e-08, 3.64421372e-11, -14381.2106, 5.60930491},
      {6.66789363, 0.0206120214, -7.36553027e-06, 1.18440761e-09, -7.0695321e-14, -16274.8521, -13.1859503}, "C3H8", "L 6/90"}},
    // TM-4513 species C3H6,propylene, data source L 7/90
    {"Propylene", {200, 1000, 6000,
      {3.83464524, 0.00329078405, 5.05228184e-05, -6.66251418e-08, 2.63707585e-11, 753.838295, 7.53410995},
      {6.03870499, 0.0162963895, -5.82130624e-06, 9.35936483e-10, -5.58602903e-14, -776.595092, -8.43824322}, "C3H6,propylene", "L 7/90"}},
    // TM-4513 species SO2, data source J 6/61
    {"SulfurDioxide", {300, 1000, 5000,
      {3.2665338, 0.0053237902, 6.8437552e-07, -5.2810047e-09, 2.5590454e-12, -36908.148, 9.66465108},
      {5.2451364, 0.0019704204, -8.0375769e-07, 1.5149969e-10, -1.0558004e-14, -37558.227, -1.07404892}, "SO2", "J 6/61"}},
    // TM-4513 species SF6, data source J 6/76
    {"SulfurHexafluoride", {300, 1000, 5000,
      {-3.8388088, 0.083221721, -0.00013181689, 9.9636154e-08, -2.9248767e-11, -148364.77, 37.1611426},
      {15.16295, 0.0043842318, -1.9486337e-06, 3.8247196e-10, -2.760505e-14, -152268.01, -54.4157194}, "SF6", "J 6/76"}},
    // TM-4513 species C7H8, data source L 1/93
    {"Toluene", {200, 1000, 6000,
      {1.611914, 0.0211188902, 8.53221453e-05, -1.32566876e-07, 5.59406109e-11, 4096.51976, 20.2973614},
      {12.939475, 0.0266921558, -9.68420108e-06, 1.5739214e-09, -9.46670482e-14, -677.035769, -46.7255302}, "C7H8", "L 1/93"}},
    // TM-4513 species H2O, data source L 8/89
    {"Water", {200, 1000, 6000,
      {4.19864056, -0.0020364341, 6.52040211e-06, -5.48797062e-09, 1.77197817e-12, -30293.7267, -0.849032208},
      {2.67703787, 0.00297318329, -7.7376969e-07, 9.44336689e-11, -4.26900959e-15, -29885.8938, 6.88255571}, "H2O", "L 8/89"}},
    // TM-4513 species Xe, data source L12/91
    {"Xenon", {200, 1000, 6000,
      {2.5, -8.9914133e-14, 2.5219686e-16, -2.92186662e-19, 1.18949218e-22, -745.375, 6.16441993},
      {2.50005322, -1.05136544e-07, 6.75326897e-11, -1.70944909e-14, 1.47681049e-18, -745.394186, 6.16412898}, "Xe", "L12/91"}},

  };

}  // namespace adsim
