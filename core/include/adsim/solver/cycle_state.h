/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// CycleState: the full simulation state at a cycle boundary, i.e. everything
// the next cycle depends on, so that a run can be restarted from it.
//
//   adsim::Simulation sim = adsim::Simulation::FromFile("deck.json");
//   sim.Run({.cycles = 20});
//   adsim::CycleState x = sim.FinalState();      // after cycle 20
//   sim.RunCycle(x);                             // exactly one cycle from x
//   adsim::CycleState gx = sim.FinalState();     // the cycle map G(x)
//
//   other.SetInitialState(x);                    // warm start: later runs of
//   other.Run();                                 // `other` start from x
//
// The state at a cycle boundary is more than the DAE vector:
//
//   * the DAE unknowns y and their time derivatives yp, in the solver's banded
//     (cell-major) order;
//   * the endpoint histories: every endpoint's logged time series of its last
//     run (StreamData). A connected inlet replays its partner's logged outlet,
//     and when the partner runs later in the cycle that log comes from the
//     previous cycle;
//   * the endpoint summaries (CycleEndpointSummary) that buffer-tank inlets
//     read from their source step.
//
// It also carries the metadata a compatibility check needs: components,
// layers with their cell counts, the step and endpoint structure and the step
// times, the number of cycles behind the state, and the DAE layout (fields).
//
// A CycleState lives for the session only: it is never written to the deck
// file (it is not part of UserData's serialisation).

#include <cstddef>
#include <string>
#include <vector>

#include "adsim/boundary_conditions/boundary_stream.h"         // StreamData
#include "adsim/boundary_conditions/cycle_endpoint_summary.h"  // CycleEndpointSummary

namespace adsim {

  class UserData;

  // One DAE field of the state (e.g. "Pressure", "CarbonDioxide Mole
  // Fraction", "CarbonDioxide Adsorbed Concentration" of one adsorbent) and
  // where its entries sit in CycleState::y.
  struct CycleStateField {
    std::string name;            // the equation name
    std::string component;       // mole fraction / adsorbed concentration; else empty
    size_t adsorbent = 0;        // adsorbed concentration: the adsorbent ID; else 0
    bool differential = false;   // the equation type (its ghost cells are algebraic)
    bool face = false;           // a face field: `cell` holds the face index (face c lies between cells c and c + 1)
    std::vector<size_t> index;   // entries of y, in cell order
    std::vector<size_t> cell;    // the cell (or face) of each entry
  };

  // The logged history and summary of one endpoint (flow stream or pressure
  // boundary), keyed by its stable endpoint ID.
  struct CycleStateEndpoint {
    size_t endpoint_id = 0;
    size_t step_id = 0;
    std::string name;            // display name (the BoundaryLog name)
    StreamData history;          // last run of the endpoint (step-relative time)
    CycleEndpointSummary summary;
  };

  struct CycleStateLayer {
    size_t id = 0;
    std::string name;
    int cells = 0;
  };

  struct CycleStateStep {
    size_t id = 0;
    std::string name;
    double time = 0.0;                  // [s] step duration when the state was taken
    std::vector<size_t> endpoint_ids;   // bottom boundary first, then top
  };

  class CycleState {
  public:
    // ---- metadata ---------------------------------------------------------
    std::vector<std::string> components;
    std::vector<CycleStateLayer> layers;
    std::vector<CycleStateStep> steps;
    // Void regions at the bottom / top of the column (one cell each).
    bool bottom_void = false;
    bool top_void = false;
    // Cycles behind the state: the run it came from had completed this many
    // cycles (including those of any warm start it began from).
    int cycles = 0;

    // ---- DAE state ----------------------------------------------------------
    std::vector<double> y;
    std::vector<double> yp;
    std::vector<CycleStateField> fields;
    // Entries of y that are differential unknowns (IDA's id = 1), ascending.
    // Ghost cells and the algebraic fields are left out.
    std::vector<size_t> differential_index;

    // ---- histories ----------------------------------------------------------
    std::vector<CycleStateEndpoint> endpoints;

    bool Empty() const { return y.empty(); }
    size_t Size() const { return y.size(); }

    // True when any endpoint carries logged history or a valid summary.
    bool HasHistory() const;

    // The differential entries of y (y[differential_index]) and the way back.
    // These are the unknowns the next cycle depends on: the algebraic ones
    // (pressure, velocity, ghost cells) are recomputed from them by the
    // consistent-initial-condition solve at the start of the run.
    std::vector<double> Differential() const;
    // Throws std::invalid_argument when the size does not match.
    void SetDifferential(const std::vector<double>& values);

    // The field with this equation name (and component / adsorbent when
    // given), or nullptr.
    const CycleStateField* FindField(const std::string& name, const std::string& component = {},
                                     size_t adsorbent = 0) const;
    // The values of one field (y[field.index]); throws std::out_of_range for
    // an unknown field.
    std::vector<double> FieldValues(const std::string& name, const std::string& component = {},
                                    size_t adsorbent = 0) const;
  };

  // Captures the state of `data` now: data.state must hold the end of a cycle
  // (as it does after a completed run) with its layout initialised.
  CycleState CaptureCycleState(UserData& data, int cycles);

  // Why `state` cannot start a run of `data`, or an empty string when it can.
  // Checks what the deck alone determines: the components, the layers and the
  // void regions, the step and endpoint structure, and the step times when the
  // state carries history that a connection or buffer tank reads. Different
  // cell counts are not a reason: such a state is remapped (RemapCycleState).
  std::string CycleStateIncompatibility(const CycleState& state, UserData& data);

  // True when the state's cells per layer differ from the deck's, i.e. the
  // state must be remapped (RemapCycleState) before it can start a run.
  bool CycleStateNeedsRemap(const CycleState& state, const UserData& data);

  // The state transferred onto the deck's grid (its cells per layer), for a
  // state that passes CycleStateIncompatibility. Initialises data.state (as a
  // run does) to obtain the deck's DAE layout; entries the state cannot supply
  // keep the deck's initial conditions.
  //
  // Within each layer every per-cell field is remapped by the overlap-weighted
  // piecewise-constant remap on the layer's uniform cells, which keeps each
  // field's layer average (sum of value x cell width) and is monotone. The gas
  // is remapped as component concentrations c_i = C y_i (C the total gas
  // concentration, or P / (R T)); the remapped composition is y_i = c_i / sum c
  // and the pressure P = (sum c) R T with the remapped fluid temperature. So
  // each component's gas and adsorbed inventory in each layer is kept to
  // rounding. Face fields (velocity: algebraic, only a guess for the
  // consistent-IC solve) are interpolated linearly in the layer coordinate.
  // Ghost and void cells (one each on both grids) are copied, and so are the
  // endpoint histories and summaries (they belong to steps, not cells).
  // Refinement by an integer ratio followed by the reverse coarsening gives
  // the original state back. Throws std::invalid_argument for an incompatible
  // state.
  CycleState RemapCycleState(const CycleState& state, UserData& data);

  // Why the DAE layout of `state` does not match `data.state` (which must be
  // initialised, as in the solver after UserData::Initialise), or empty.
  std::string CycleStateLayoutMismatch(const CycleState& state, const UserData& data);

  // Writes the DAE vector into data.state (layout must match).
  void ApplyCycleStateDAE(const CycleState& state, UserData& data);
  // Writes the endpoint histories and summaries into the deck's endpoints
  // (matched by endpoint ID; an endpoint the state does not know is cleared).
  void ApplyCycleStateHistories(const CycleState& state, UserData& data);

}  // namespace adsim
