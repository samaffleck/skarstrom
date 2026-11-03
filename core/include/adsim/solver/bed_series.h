/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Bed variables read back from the equation log at the right axial positions.
//
// The log stores each variable as the solver holds it:
//
//   * cell variables (pressure, temperature, mole fraction, adsorbed
//     concentration, ...): one value per cell, N values per row, where cells 0
//     and N-1 are the ghost cells that carry the boundary conditions. A ghost
//     value is not a value at any point of the bed (an inlet temperature ghost
//     holds 2 T_feed - T_1), so it is never returned as one. Cell n (1..N-2)
//     sits at its centre, (z_{n-1} + z_n) / 2, with z the face positions.
//   * face variables (velocity, volumetric and molar flow rate, approach to
//     fluidisation): one value per face, N-1 values per row, at the face
//     positions z_0 = 0 (bottom of the column) to z_{N-2} = L (top).
//
// The bed ends (z = 0 and z = L) are the physical boundary values:
//
//   * face variables: the boundary face value.
//   * pressure, fluid temperature and mole fraction, where the step has a
//     logged endpoint (a flow stream or a pressure boundary) at that end: the
//     endpoint's logged value, so the series agrees with stream_series. At an
//     inlet that is the value the boundary imposes (the feed composition, the
//     feed temperature after any compressor, the set pressure); at an outlet
//     it is the bed-end face value. Several streams at one end are mixed by
//     their molar flows.
//   * otherwise (a closed end, other variables, or no boundary log): the
//     face value, the mean of the ghost and the end cell, which is how the
//     boundary conditions and fluxes see the face. The exceptions: a mole
//     fraction takes the ghost value, which is the face value itself (the
//     inlet composition, or a copy of the end cell at an outlet or closed
//     end), and adsorbent-only variables (adsorbed and equilibrium loading,
//     mass-transfer coefficient, solid temperature, particle heat-transfer
//     coefficient) take the value of the nearest packed cell.
//
// A position between a bed end and the centre of the end cell interpolates
// between the boundary value and that cell.

#include <string>
#include <vector>

#include "adsim/solver/log.h"

namespace adsim {

  class Step;

  // Spatial location for time-series extraction.
  enum class ParamSpatialLoc { Top = 0, Bottom = 1, Position = 2 };

  // One logged bed variable: a bulk VariableType, or a ComponentVariableType
  // of one component.
  struct BedVariable {
    bool is_component = false;
    VariableType bulk = VariableType::PRESSURE;
    ComponentVariableType comp = ComponentVariableType::MOLE_FRACTION;
    std::string component;

    static BedVariable Bulk(VariableType t) { BedVariable v; v.bulk = t; return v; }
    static BedVariable Component(ComponentVariableType t, std::string name) {
      BedVariable v;
      v.is_component = true;
      v.comp = t;
      v.component = std::move(name);
      return v;
    }
  };

  // Where the bed-end values of pressure, temperature and mole fraction come
  // from: the boundary-stream log of the run and the steps that ran (for the
  // endpoints at each bed end). Without them the bed ends are face values.
  struct BedBoundarySources {
    const BoundaryLog* boundary_log = nullptr;
    const std::vector<Step>* steps = nullptr;
  };

  // Centres of the bed cells 1..N-2 [m] (no ghost cells), from the face
  // positions of the domain.
  std::vector<double> CellCentrePositions(const DomainVectors& domain);

  // True when a logged row of `row_size` values lives on the faces (velocity,
  // flow rates), false when it lives on the cells (ghost cells included).
  bool IsFaceRow(const DomainVectors& domain, size_t row_size);

  // The log entry of `var` over one step of one cycle (empty when not logged).
  VariableLogEntry GetBedVariable(const EquationLog& eq_log, int cycle_index, int step_index,
                                  const BedVariable& var);

  // Time series of `var` over one cycle (the final logged cycle when
  // cycle_index < 0), concatenated across its steps, with times relative to
  // the start of that cycle (the logged times are unchanged otherwise). Each
  // row is evaluated at the bed end (Top: z = L, Bottom: z = 0, as above) or
  // at z_pos [m] (Position; clamped to [0, L]). Returns false, with empty
  // outputs, when there is no data.
  bool ExtractBedTimeSeries(const EquationLog& eq_log, const BedVariable& var, ParamSpatialLoc where,
                            double z_pos, std::vector<double>& out_times, std::vector<double>& out_values,
                            int cycle_index = -1, const BedBoundarySources& sources = {});

}  // namespace adsim
