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

  class AdsorbentDirectory;
  class Bed;
  class Cycle;
  class Fluid;

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

  // One warning per adsorbent whose mixing rule is IAST while one of its
  // components uses an extended (competitive) isotherm kind: IAST evaluates
  // every isotherm in its pure-component limit, so the kind's own competition
  // term is not applied and the plain kind with the same parameters is meant.
  std::vector<std::string> MixingRuleWarnings(const AdsorbentDirectory& adsorbents, const Fluid& fluid);

  // One warning per component whose gas heat capacity is held at the initial
  // conditions instead of following the temperature: in the ideal-gas mode, a
  // component without an ideal-gas polynomial (Fluid::
  // ComponentsWithHeldHeatCapacity; ideal_gas_heat_capacity.h). Requires
  // Fluid::Initialise. Empty in the constant mode.
  std::vector<std::string> GasHeatCapacityWarnings(const Fluid& fluid);

  // The warnings of CheckLayerComposition (porous_media.h) for every layer, in
  // layer order: one per solid at mass fraction 0 (it takes no part in the
  // run). Layer composition errors stop the run in the deck checks instead.
  std::vector<std::string> LayerCompositionWarnings(const AdsorbentDirectory& adsorbents, const Bed& bed);

  // One warning per packed layer whose film correlation is Ranz-Marshall
  // (film_correlation.h): a single-sphere correlation that under-predicts the
  // particle heat and mass transfer coefficients in a bed at moderate Re. It
  // is offered to reproduce papers that use it; not an error.
  std::vector<std::string> FilmCorrelationWarnings(const Bed& bed);

  // Structured layers (monolith, laminate; porous_media.h): one warning for
  // local thermal equilibrium (LTE) in a bed with a structured layer: the
  // thin walls exchange heat with the gas through the channel film, which
  // LTE takes as infinitely fast; LTNE is recommended.
  // (The packed-bed closures that read the particle diameter, Ruthven axial
  // and heat dispersion, Yagi-Kunii heat dispersion and Dixon or Leva wall
  // coefficients, are deck errors on a structured layer.)
  // - per structured layer whose bed-to-wall coefficient is a constant 0 while
  //   the column wall is not adiabatic (an outer-face coefficient h_o > 0 or
  //   a correlation, a fixed-temperature wall end, or, given `cycle`, a step
  //   with a heating medium): the wall then exchanges no heat with that
  //   layer. (Switching a layer to a structured kind can leave h_i = 0, the
  //   constant model's default.) Not for an isothermal run.
  std::vector<std::string> StructuredLayerWarnings(const AdsorbentDirectory& adsorbents, const Bed& bed,
                                                   const Cycle* cycle = nullptr);

  // Reactions (objects/reaction.h): one warning per solid with reactions that
  // no layer uses (its reactions take no part, and the deck checks skip
  // them); one per reaction with a nonzero heat of reaction in an isothermal
  // run (the heat is not applied); one per reaction whose element balance
  // could not be checked (a component without an element formula, e.g. Air).
  // With `order_notes`, also an informational note per reaction whose orders
  // differ from -nu (reactants) and 0 (products), i.e. a non-elementary rate.
  // Expressions must already be evaluated.
  std::vector<std::string> ReactionWarnings(const AdsorbentDirectory& adsorbents, const Bed& bed,
                                            bool order_notes = false);

} // namespace adsim
