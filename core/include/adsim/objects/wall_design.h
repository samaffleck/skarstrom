/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Minimum wall thickness of a cylindrical column from a design pressure, for
// preliminary sizing. A final design must be checked by a qualified engineer
// against the edition of the code in use.
//
// Internal pressure (the design pressure, gauge):
//   ASME BPVC VIII-1 UG-27(c)(1), circumferential stress:
//     t = P R / (S E - 0.6 P),   R = inside radius + corrosion allowance,
//   and Appendix 1-2(a)(1) for thick shells (P > 0.385 S E):
//     t = R (sqrt(Z) - 1),       Z = (S E + P) / (S E - P).
//   ASME B31.3 304.1.2, eq. (3b) (inside diameter d):
//     t = P (d + 2c) / (2 (S E W - P (1 - Y))),   W = 1, Y = 0.4,
//   valid while t < D/6 and P / (S E) <= 0.385 (304.1.2(b)).
//   The nominal thickness is t + c (the corrosion allowance).
//
// External pressure (vacuum service: full vacuum, 101 325 Pa outside):
//   The formula method of ASME Code Case 2286 / VIII-2 4.4.5.1, the chart-free
//   equivalent of UG-28 (B31.3 304.1.3 refers to the same rules):
//     Mx  = L / sqrt(Ro t),  Fhe = 1.6 Ch E t / Do  (elastic buckling stress)
//     Fic from Fhe / Sy (inelastic buckling stress), FS from Fic / Sy,
//     Pa  = 2 (Fic / FS) (t / Do).
//   t is the corroded thickness, Do the outside diameter and L the column
//   length (the packed length; heads and stiffening rings are not modelled).
//   Sy is the lower bound YieldStrengthLowerBound (wall.h).
//
// Allowable stresses are ASME II-D (VIII-1) values for both codes; B31.3's
// Table A-1 is equal for the austenitic steels and higher for carbon steel, so
// this is conservative for B31.3.

#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "adsim/objects/wall.h"

namespace adsim {

  class Bed;

  // Full vacuum: the external design pressure of a wall in vacuum service [Pa].
  inline constexpr double kFullVacuumExternalPressure = 101325.0;

  struct WallDesignInput {
    WallDesignCode code = WallDesignCode::ASME_VIII_1;
    WallMaterial material = WallMaterial::SS316L;
    double inner_diameter = 0.0;      // [m]
    double length = 0.0;              // [m] column length (external pressure only)
    double design_pressure = 0.0;     // [Pa gauge] internal
    double design_temperature = 293.0;// [K]
    double joint_efficiency = 1.0;    // E [-]
    double corrosion_allowance = 0.0; // c [m]
    bool vacuum_service = false;
  };

  struct WallDesignResult {
    bool valid = false;
    std::string message;              // why the design is not valid; empty when valid
    std::vector<std::string> notes;   // assumptions and limits of the result
    double allowable_stress = std::numeric_limits<double>::quiet_NaN();     // S [Pa]
    double internal_thickness = 0.0;  // [m] nominal (incl. c) for the design pressure; 0 if none
    double external_thickness = 0.0;  // [m] nominal (incl. c) for full vacuum; 0 if not in vacuum service
    double required_thickness = std::numeric_limits<double>::quiet_NaN();   // [m] the larger of the two
    bool external_governs = false;
    bool thick_shell = false;         // VIII-1 Appendix 1-2 was used
  };

  WallDesignResult DesignWall(const WallDesignInput& in);

  // Maximum allowable internal working pressure [Pa gauge] of a wall of
  // nominal thickness `thickness` (the inverse of the internal-pressure rule).
  // NaN when it cannot be evaluated (Custom material, thickness <= c, or the
  // code's validity limit is exceeded).
  double MaxAllowableInternalPressure(const WallDesignInput& in, double thickness);

  // Allowable external pressure [Pa] of a wall of nominal thickness
  // `thickness` (Code Case 2286). NaN when it cannot be evaluated.
  double AllowableExternalPressure(const WallDesignInput& in, double thickness);

  // The design input of a bed's wall: its design fields, the inside diameter
  // and the packed length. nullopt for a non-cylindrical column (the formulas
  // are for cylindrical shells). Expressions must already be evaluated.
  std::optional<WallDesignInput> WallDesignInputFor(const Bed& bed);

  // DesignWall(WallDesignInputFor(bed)); a not-valid result whose message says
  // why for a non-cylindrical column or a wall with no design in use.
  WallDesignResult DesignBedWall(const Bed& bed);

} // namespace adsim
