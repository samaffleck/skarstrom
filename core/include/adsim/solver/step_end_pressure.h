/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// The pressure each step of a cycle ends at.
//
// A plant study often bounds a step's final pressure (an evacuation must not
// go below the vacuum pump's floor, a pressurisation must reach the feed
// pressure). The solvers record these values with the KPIs when each cycle
// completes (KPIHistory), so they cover every cycle of a run, not only the
// logged ones, and the optimiser and the sensitivity studies read them as
// outputs: "step[<step>].<name>" with <name> one of StepOutputNames().

#include <limits>
#include <string>
#include <string_view>
#include <vector>

namespace adsim {

  class EquationLog;
  struct BoundaryLog;
  class Step;

  // Pressures of one step at its last logged instant [Pa]. `bottom` and `top`
  // are the bed-end values that Simulation::Series("pressure", ..., Bottom or
  // Top) gives (the endpoint's logged face pressure where the step has an
  // endpoint at that end, else the face value); `bed_min`, `bed_mean` (cell
  // width weighted) and `bed_max` are over the packed cells (no ghost or void
  // cells). NaN when the step was not logged.
  struct StepEndPressure {
    int cycle = 0;
    int step = 0;              // index of the step in the cycle
    std::size_t step_id = 0;   // the step's deck ID (Step::id_)
    std::string step_name;
    double time = std::numeric_limits<double>::quiet_NaN();      // [s] the step's end, from the cycle start
    double bottom = std::numeric_limits<double>::quiet_NaN();    // [Pa]
    double top = std::numeric_limits<double>::quiet_NaN();       // [Pa]
    double bed_min = std::numeric_limits<double>::quiet_NaN();   // [Pa]
    double bed_mean = std::numeric_limits<double>::quiet_NaN();  // [Pa]
    double bed_max = std::numeric_limits<double>::quiet_NaN();   // [Pa]
  };

  // One record per step of `steps` for a logged cycle (empty when the cycle
  // is not logged).
  std::vector<StepEndPressure> ComputeStepEndPressures(const std::vector<Step>& steps, const EquationLog& eq_log,
                                                       const BoundaryLog& b_log, int cycle);

  // The step outputs a study or the optimiser can watch: end_pressure_bottom,
  // end_pressure_top, end_pressure_min, end_pressure_mean, end_pressure_max
  // (StepEndPressure::bottom ... bed_max).
  const std::vector<std::string>& StepOutputNames();

  // The value of a step output of a record; throws std::invalid_argument for
  // a name not in StepOutputNames().
  double StepOutputValue(const StepEndPressure& record, std::string_view name);

  // True when an output's subject is a step reference, "step[...]" (instead
  // of a component name).
  bool IsStepReference(std::string_view subject);

  // The index in `steps` of the step a reference "step[<ID or name>]" names:
  // the step with that deck ID when the text is a number and such a step
  // exists, otherwise the step with that name. Throws std::invalid_argument
  // when it names no step (or is not a step reference).
  int ResolveStepReference(const std::vector<Step>& steps, std::string_view subject);

}  // namespace adsim
