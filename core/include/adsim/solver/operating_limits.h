/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Operating limits of a finished run, from one logged cycle (the last by
// default): the extremes of pressure, wall temperature and approach to
// fluidisation in every step, and warnings where the cycle leaves the wall's
// pressure design (wall_design.h) or lifts the bed. The solver logs the
// warnings once at the end of a run that completes (LogLevel::LOW).
//
// - Pressure is taken over the column (packed cells and void cells), wall
//   temperature over the packed cells, and approach to fluidisation over the
//   faces between packed cells (the faces at a void or ghost cell compare the
//   pressure across a different region).
// - The design checks run only when the wall's design is in use
//   (BedWall::DesignInUse: a positive design pressure or vacuum service).
//   Gauge pressure is relative to kAtmosphericPressure.
// - Approach to fluidisation is the upward pressure gradient over the weight
//   of the packing less the gas, (-dP/dz) / ((rho_b - rho) g); above 1 an
//   unrestrained bed in upward flow lifts. It is checked whatever the design.
// - Values are those at the logged times; a peak between two log points is
//   not seen. A step's first logged point is the state the step before ended
//   in: an extreme a step reaches only there is in its summary but gives no
//   warning, since the step before reports it.

#include <limits>
#include <string>
#include <vector>

namespace adsim {

  class UserData;
  class EquationLog;

  // The reference for gauge pressures [Pa].
  inline constexpr double kAtmosphericPressure = 101325.0;

  struct StepOperatingSummary {
    int step = 0;                     // step index in the cycle
    std::string step_name;
    double max_pressure = std::numeric_limits<double>::quiet_NaN();          // [Pa] absolute
    double min_pressure = std::numeric_limits<double>::quiet_NaN();          // [Pa] absolute
    double max_wall_temperature = std::numeric_limits<double>::quiet_NaN();  // [K]
    double max_approach_to_fluidisation = std::numeric_limits<double>::quiet_NaN();  // [-]
    double fluidisation_position = std::numeric_limits<double>::quiet_NaN(); // [m] from the bottom of the bed
    double fluidisation_time = std::numeric_limits<double>::quiet_NaN();     // [s] into the step
  };

  enum class OperatingLimit {
    DesignPressure,     // pressure above the design pressure
    Vacuum,             // pressure below atmospheric, wall not in vacuum service
    DesignTemperature,  // wall temperature above the design temperature
    Fluidisation        // approach to fluidisation above 1
  };

  inline constexpr const char* kOperatingLimitList[] = {
    "design_pressure",
    "vacuum",
    "design_temperature",
    "fluidisation"
  };

  struct OperatingWarning {
    OperatingLimit limit = OperatingLimit::Fluidisation;
    int step = 0;                     // step index in the cycle
    double value = 0.0;               // the extreme reached (SI; gauge Pa for DesignPressure, absolute Pa for Vacuum)
    double limit_value = 0.0;         // the limit it passed (same units)
    std::string message;
  };

  struct OperatingLimitsReport {
    int cycle = -1;                   // 0-based cycle checked; -1 when nothing is logged
    bool design_checked = false;      // the wall's pressure design was in use
    std::vector<StepOperatingSummary> steps;
    std::vector<OperatingWarning> warnings;
  };

  // `cycle` < 0 checks the last logged cycle.
  OperatingLimitsReport CheckOperatingLimits(const UserData& data, const EquationLog& log, int cycle = -1);

} // namespace adsim
