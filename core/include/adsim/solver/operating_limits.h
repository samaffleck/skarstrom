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
//   pressure across a different region). Faces touching a structured layer
//   (monolith, laminate) are not checked: its solid is fixed and cannot lift.
// - The channel Reynolds number (the logged Reynolds number of a structured
//   layer's cells, rho |u| D_h / (OFA mu)) is checked against
//   kLaminarChannelReynolds, over which its laminar closures do not hold.
// - The design checks run only when the wall's design is in use
//   (BedWall::DesignInUse: a positive design pressure or vacuum service).
//   Gauge pressure is relative to kAtmosphericPressure. A pressure limit
//   counts as passed by more than 2 kPa (decks often take 1 bar as
//   atmospheric), a temperature limit by more than 0.01 K.
// - Outlet backflow: gas that enters the bed through an outlet endpoint with
//   no connection to another step (a valve outlet whose bed fell below the
//   downstream pressure, a pressure outlet above the bed). The gas that
//   flows back carries the bed-end composition (an outlet is zero-gradient),
//   so it is not the stream's gas; a valve's check valve stops it. Checked
//   for the last logged cycle only, from the endpoints' own logs, against
//   kOutletBackflowFraction of the cycle feed. Connected endpoints are not
//   checked: flow both ways through a connection is expected.
// - Connection pressure: a diagnostic. A connected endpoint replays its
//   partner's logged run and does not enforce the order of the pressures,
//   so gas can pass from the lower to the higher pressure. Reported for the
//   last logged cycle only, from the endpoints' own logs, when the gross
//   moles that enter the bed through one connected endpoint in one step
//   while its bed-end pressure is above the partner's logged pressure at
//   the same instant (by more than kConnectionPressureMargin) pass
//   kConnectionReverseFlowFraction of the cycle feed. Pump streams, which
//   raise the pressure by design, are not checked.
// - Approach to fluidisation is the upward pressure gradient over the weight
//   of the packing less the gas, (-dP/dz) / ((rho_b - rho) g); above 1 an
//   unrestrained bed in upward flow lifts. It is checked whatever the design.
// - A warning switched off in the deck's settings (IDASettings::warnings_,
//   warning_settings.h) is not reported; the step summaries are complete.
// - Values are those at the logged times; a peak between two log points is
//   not seen. A step's first logged point is the state the step before ended
//   in: an extreme a step reaches only there is in its summary but gives no
//   warning, since the step before reports it.

#include <limits>
#include <string>
#include <vector>

#include "adsim/math/constants.h"

namespace adsim {

  class UserData;
  class EquationLog;

  // The reference for gauge pressures [Pa].
  inline constexpr double kAtmosphericPressure = constants::kStandardAtmosphere;

  struct StepOperatingSummary {
    int step = 0;                     // step index in the cycle
    std::string step_name;
    double max_pressure = std::numeric_limits<double>::quiet_NaN();          // [Pa] absolute
    double min_pressure = std::numeric_limits<double>::quiet_NaN();          // [Pa] absolute
    double max_wall_temperature = std::numeric_limits<double>::quiet_NaN();  // [K]
    double max_approach_to_fluidisation = std::numeric_limits<double>::quiet_NaN();  // [-]
    double fluidisation_position = std::numeric_limits<double>::quiet_NaN(); // [m] from the bottom of the bed
    double fluidisation_time = std::numeric_limits<double>::quiet_NaN();     // [s] into the step
    // Largest channel Reynolds number rho |u| D_h / (OFA mu) over the
    // structured layers' cells; NaN without a structured layer.
    double max_channel_reynolds = std::numeric_limits<double>::quiet_NaN();  // [-]
  };

  enum class OperatingLimit {
    DesignPressure,     // pressure above the design pressure
    Vacuum,             // pressure below atmospheric, wall not in vacuum service
    DesignTemperature,  // wall temperature above the design temperature
    Fluidisation,       // approach to fluidisation above 1
    ChannelReynolds,    // channel Reynolds number above 2000 in a structured layer
    OutletBackflow,     // gas flows back into the bed through an unconnected outlet
    ConnectionPressure  // gas through a connection from the lower to the higher pressure
  };

  inline constexpr const char* kOperatingLimitList[] = {
    "design_pressure",
    "vacuum",
    "design_temperature",
    "fluidisation",
    "channel_reynolds",
    "outlet_backflow",
    "connection_pressure"
  };

  // Above this channel Reynolds number the laminar closures of a structured
  // layer (Hagen-Poiseuille, the fully developed Nu and Sh, Taylor-Aris) no
  // longer hold: flow in straight channels turns transitional near 2000.
  inline constexpr double kLaminarChannelReynolds = 2000.0;

  // Backflow through an outlet is reported when the gross moles that enter
  // the bed through one unconnected outlet endpoint in one step pass this
  // fraction of the cycle feed (the gross moles in through every unconnected
  // inlet endpoint of the cycle).
  inline constexpr double kOutletBackflowFraction = 1e-3;

  // Flow from the lower to the higher pressure through a connection is
  // reported when the gross moles that enter the bed through one connected
  // endpoint in one step, while it is the higher, pass this fraction of the
  // cycle feed. The same share as an outlet backflow: below it the gas moved
  // does not show in the KPIs at the precision they are reported to.
  inline constexpr double kConnectionReverseFlowFraction = 1e-3;

  // A connected bed counts as the higher only when its bed-end pressure is
  // above the partner's logged pressure by more than this, so the lag of
  // the logged partner and the solver tolerance near equal pressures do not
  // count as reverse flow.
  inline constexpr double kConnectionPressureMargin = 1000.0;  // [Pa]

  struct OperatingWarning {
    OperatingLimit limit = OperatingLimit::Fluidisation;
    int step = 0;                     // step index in the cycle
    double value = 0.0;               // the extreme reached (SI; gauge Pa for DesignPressure, absolute Pa for Vacuum;
                                      // a fraction of the cycle feed for OutletBackflow and ConnectionPressure)
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
