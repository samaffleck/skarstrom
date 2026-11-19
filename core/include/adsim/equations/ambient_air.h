/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Properties of the dry air around a column, and the Churchill-Chu
// natural-convection and Churchill-Bernstein forced-convection correlations
// (and Churchill's blend of the two) that the wall-to-ambient heat transfer
// models use. Pure functions, so the tests can check them directly.

#include <algorithm>
#include <cmath>

namespace adsim::ambient_air {

  // The surroundings are taken at one standard atmosphere.
  inline constexpr double kPressure = 101325.0;        // [Pa]
  inline constexpr double kMolarMass = 0.0289647;      // dry air [kg/mol]
  inline constexpr double kGasConstant = 8.314462618;  // [J/(mol K)]
  inline constexpr double kGravity = 9.80665;          // standard gravity [m/s2]
  // c_p of dry air at 300 K (Incropera et al., Table A.4: 1006 to 1014
  // J/(kg K) from 250 to 400 K, so within 1 % over that range).
  inline constexpr double kHeatCapacity = 1007.0;      // [J/(kg K)]

  // Dynamic viscosity [Pa s] by Sutherland's law, mu0 = 1.716e-5 Pa s at
  // T0 = 273 K, S = 111 K (White, Viscous Fluid Flow, 3rd ed., 2006,
  // Sec. 1-3). Within 2 % of Incropera et al.'s Table A.4 from 300 to 400 K.
  inline double Viscosity(double T) {
    constexpr double mu0 = 1.716e-5, T0 = 273.0, S = 111.0;
    return mu0 * std::pow(T / T0, 1.5) * (T0 + S) / (T + S);
  }

  // Thermal conductivity [W/(m K)] by Sutherland's law, k0 = 0.0241 W/(m K)
  // at T0 = 273 K, S = 194 K (White 2006, Sec. 1-3). Within 2 % of
  // Incropera et al.'s Table A.4 from 300 to 400 K.
  inline double ThermalConductivity(double T) {
    constexpr double k0 = 0.0241, T0 = 273.0, S = 194.0;
    return k0 * std::pow(T / T0, 1.5) * (T0 + S) / (T + S);
  }

  // Ideal-gas density [kg/m3] at the ambient pressure.
  inline double Density(double T) { return kPressure * kMolarMass / (kGasConstant * T); }

  inline double Prandtl(double T) { return kHeatCapacity * Viscosity(T) / ThermalConductivity(T); }

  // Rayleigh number on the length L [m] for a surface at T_wall in air at
  // T_inf, with the properties at the film temperature (T_wall + T_inf) / 2
  // and beta = 1 / T_film (ideal gas): Ra = g beta |dT| L^3 Pr / nu^2.
  inline double Rayleigh(double T_wall, double T_inf, double L) {
    const double T_film = 0.5 * (T_wall + T_inf);
    const double nu = Viscosity(T_film) / Density(T_film);
    return kGravity * std::abs(T_wall - T_inf) / T_film * L * L * L * Prandtl(T_film) / (nu * nu);
  }

  // Churchill and Chu (1975), average Nusselt number over a vertical plate of
  // height L, laminar and turbulent (all Ra):
  //   Nu_L = {0.825 + 0.387 Ra_L^(1/6) / [1 + (0.492/Pr)^(9/16)]^(8/27)}^2.
  // S.W. Churchill, H.H.S. Chu, Int. J. Heat Mass Transfer 18 (1975) 1323-1329.
  inline double ChurchillChuVerticalPlateNu(double Ra, double Pr) {
    const double f = std::pow(1.0 + std::pow(0.492 / Pr, 9.0 / 16.0), 8.0 / 27.0);
    const double s = 0.825 + 0.387 * std::pow(std::max(Ra, 0.0), 1.0 / 6.0) / f;
    return s * s;
  }

  // Churchill and Chu (1975), average Nusselt number over a long horizontal
  // cylinder of diameter D (Ra_D up to about 1e12):
  //   Nu_D = {0.60 + 0.387 Ra_D^(1/6) / [1 + (0.559/Pr)^(9/16)]^(8/27)}^2.
  // S.W. Churchill, H.H.S. Chu, Int. J. Heat Mass Transfer 18 (1975) 1049-1053.
  inline double ChurchillChuHorizontalCylinderNu(double Ra, double Pr) {
    const double f = std::pow(1.0 + std::pow(0.559 / Pr, 9.0 / 16.0), 8.0 / 27.0);
    const double s = 0.60 + 0.387 * std::pow(std::max(Ra, 0.0), 1.0 / 6.0) / f;
    return s * s;
  }

  // Yang (1985), average Nusselt number over a vertical cylinder of height H
  // and diameter d, at any H/d:
  //   Nu_H = || 0.36 H/d , 0.150 [Ra_H / ||1, 0.492/Pr||_{9/16}]^(1/3) ||_{1/2},
  // with ||a, b||_p = (a^p + b^p)^(1/p). This is Churchill and Chu's plate
  // formula (0.150 ~ 0.387^2) with the plate's conduction term 0.825^2
  // replaced by 0.36 H/d: at small Ra, h -> 0.36 k/d (the conduction limit
  // of a horizontal cylinder), so a slender column loses more heat than a
  // plate of its height. S.M. Yang, "General correlating equations for free
  // convection heat transfer from a vertical cylinder", Int. Symp. Heat
  // Transfer, Peking, 1985, pp. 153-159; as quoted by Popiel, Wojtkowiak
  // and Bober (2007), Exp. Therm. Fluid Sci. 32, 607, and Jaffer and Jaffer
  // (2026), arXiv:2512.06019. Original not checked.
  inline double YangVerticalCylinderNu(double Ra, double Pr, double H_over_d) {
    const double plate_denominator = std::pow(1.0 + std::pow(0.492 / Pr, 9.0 / 16.0), 16.0 / 9.0);
    const double conduction = 0.36 * H_over_d;
    const double convection = 0.150 * std::cbrt(std::max(Ra, 0.0) / plate_denominator);
    const double s = std::sqrt(conduction) + std::sqrt(convection);
    return s * s;
  }

  // h [W/(m2 K)] of a vertical column wall of height H and outer diameter d
  // at T_wall in still air at T_inf, by Yang (1985); Ra on H, film properties.
  inline double YangCoefficient(double T_wall, double T_inf, double height, double outer_diameter) {
    const double T_film = 0.5 * (T_wall + T_inf);
    const double Nu = YangVerticalCylinderNu(Rayleigh(T_wall, T_inf, height), Prandtl(T_film), height / outer_diameter);
    return Nu * ThermalConductivity(T_film) / height;
  }

  // Natural-convection coefficient h = Nu k / L_c [W/(m2 K)] of a column wall
  // at T_wall in still air at T_inf. A vertical column is a plate of height
  // `height`; a horizontal one is a cylinder of diameter `outer_diameter`.
  inline double ChurchillChuCoefficient(bool vertical, double T_wall, double T_inf, double height,
                                        double outer_diameter) {
    const double L = vertical ? height : outer_diameter;
    const double T_film = 0.5 * (T_wall + T_inf);
    const double Ra = Rayleigh(T_wall, T_inf, L);
    const double Pr = Prandtl(T_film);
    const double Nu = vertical ? ChurchillChuVerticalPlateNu(Ra, Pr) : ChurchillChuHorizontalCylinderNu(Ra, Pr);
    return Nu * ThermalConductivity(T_film) / L;
  }

  // Churchill and Bernstein (1977), average Nusselt number over a long
  // cylinder of diameter D in crossflow, for Re_D Pr > 0.2:
  //   Nu_D = 0.3 + 0.62 Re^(1/2) Pr^(1/3) / [1 + (0.4/Pr)^(2/3)]^(1/4)
  //              * [1 + (Re/282000)^(5/8)]^(4/5).
  // S.W. Churchill, M. Bernstein, "A correlating equation for forced
  // convection from gases and liquids to a circular cylinder in crossflow",
  // J. Heat Transfer 99 (1977) 300-306 (equation number to be confirmed).
  inline double ChurchillBernsteinNu(double Re, double Pr) {
    const double re = std::max(Re, 0.0);
    const double laminar = 0.62 * std::sqrt(re) * std::cbrt(Pr) / std::pow(1.0 + std::pow(0.4 / Pr, 2.0 / 3.0), 0.25);
    return 0.3 + laminar * std::pow(1.0 + std::pow(re / 282000.0, 5.0 / 8.0), 4.0 / 5.0);
  }

  // Reynolds number on the length L [m] of air at T [K] flowing at v [m/s].
  inline double Reynolds(double velocity, double L, double T) { return Density(T) * std::abs(velocity) * L / Viscosity(T); }

  // Forced-convection coefficient h [W/(m2 K)] of a column wall of outer
  // diameter d at T_wall in air at T_inf crossing it at v: Churchill and
  // Bernstein, with the air properties at the film temperature.
  inline double ChurchillBernsteinCoefficient(double T_wall, double T_inf, double velocity, double outer_diameter) {
    const double T_film = 0.5 * (T_wall + T_inf);
    const double Nu = ChurchillBernsteinNu(Reynolds(velocity, outer_diameter, T_film), Prandtl(T_film));
    return Nu * ThermalConductivity(T_film) / outer_diameter;
  }

  // Mixed (combined natural and forced) convection from a column wall of
  // outer diameter d and height H at T_wall, in air at T_inf crossing it at
  // v, by Churchill's (1977) blend of the two coefficients:
  //   h = (h_F^n + h_N^n)^(1/n),
  // always adding (the wind crosses the column, normal to the buoyant flow).
  // h_N is the natural-convection coefficient of the orientation: Yang's
  // vertical cylinder (Ra on H) or Churchill and Chu's horizontal cylinder
  // (Ra on d). h_F is Churchill and Bernstein's forced-convection
  // coefficient on d without its stagnant constant 0.3, which h_N already
  // holds as its conduction limit: h_F -> 0 as v -> 0, so v = 0 gives the
  // natural coefficient exactly. The blend is of h, not Nu: the two use
  // different lengths. S.W. Churchill, "A comprehensive correlating equation
  // for laminar, assisting, forced and free convection", AIChE J. 23 (1977)
  // 10-16; n = 3 for crossflow over a cylinder, as Incropera et al.,
  // Sec. 9.9 (to be confirmed).
  inline double MixedConvectionCoefficient(bool vertical, double T_wall, double T_inf, double velocity,
                                           double height, double outer_diameter, double exponent) {
    const double T_film = 0.5 * (T_wall + T_inf);
    const double h_F = ChurchillBernsteinCoefficient(T_wall, T_inf, velocity, outer_diameter) -
                       0.3 * ThermalConductivity(T_film) / outer_diameter;
    const double h_N = vertical ? YangCoefficient(T_wall, T_inf, height, outer_diameter)
                                : ChurchillChuCoefficient(false, T_wall, T_inf, height, outer_diameter);
    const double n = exponent;
    return std::pow(std::pow(std::max(h_F, 0.0), n) + std::pow(h_N, n), 1.0 / n);
  }

  // Grashof number on L [m] (film properties, as Rayleigh): Gr = Ra / Pr.
  inline double Grashof(double T_wall, double T_inf, double L) {
    return Rayleigh(T_wall, T_inf, L) / Prandtl(0.5 * (T_wall + T_inf));
  }

}  // namespace adsim::ambient_air
