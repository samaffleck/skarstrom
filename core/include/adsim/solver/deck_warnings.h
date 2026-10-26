/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Advisory checks of a deck that do not stop a run. The solver logs them once,
// at the start of a run (LogLevel::LOW), after "Solver initialised".

#include <string>
#include <vector>

namespace adsim {

  class Bed;
  class Cycle;

  // One warning per step end whose CONSTANT pressure boundary value differs
  // from the bed pressure the step is expected to start at by more than
  // `relative_threshold` (0.25 = 25 %). A constant boundary value applies from
  // the first instant of the step, so the bed end sees an instantaneous
  // pressure change the integrator has to resolve with very small steps;
  // linear and exponential boundary values start from the bed pressure and
  // never trigger this. The expected start pressure is estimated from the
  // deck, not the solution: the end value of the most recent pressure boundary
  // before the step (the same bed end if it has one, else the other end),
  // walking back cyclically through steps with no pressure boundary; for the
  // first step the initial pressure is checked as well. Steps whose
  // predecessors have no pressure boundary at all give no warning. Boundary
  // values must already be evaluated (UserData::Initialise).
  std::vector<std::string> BoundaryPressureStepChangeWarnings(
      const Cycle& cycle, double initial_pressure, double relative_threshold = 0.25);

  // When the wall's pressure design is in use (BedWall::DesignInUse): one
  // warning if the wall is thinner than the minimum its design code requires
  // (wall_design.h), or if the requirement cannot be evaluated (a column that
  // is not cylindrical, a Custom material, a pressure no wall can hold).
  // Expressions must already be evaluated.
  std::vector<std::string> WallDesignWarnings(const Bed& bed);

} // namespace adsim
