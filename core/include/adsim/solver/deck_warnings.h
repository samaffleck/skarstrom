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
  class InitialConditions;

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

  // The warnings of CheckLayerComposition (porous_media.h) for every layer, in
  // layer order: one per solid at mass fraction 0 (it takes no part in the
  // run). Layer composition errors stop the run in the deck checks instead.
  std::vector<std::string> LayerCompositionWarnings(const AdsorbentDirectory& adsorbents, const Bed& bed);

  // The warnings of CheckMassTransferModel (mass_transfer_model.h) for every
  // solid a layer uses: the film resistance on with a film resistance factor
  // of 0, and a second-order LDF with a physical resistance on. Its errors
  // stop the run in the deck checks instead.
  std::vector<std::string> MassTransferWarnings(const AdsorbentDirectory& adsorbents, const Bed& bed);

  // One warning per packed layer whose film correlation is Ranz-Marshall
  // (film_correlation.h): a single-sphere correlation that under-predicts the
  // particle heat and mass transfer coefficients in a bed at moderate Re. It
  // is offered to reproduce papers that use it; not an error.
  std::vector<std::string> FilmCorrelationWarnings(const Bed& bed);

  // One warning when a layer uses the vertical-column Churchill-Chu natural
  // convection kind and the column is too slender to be a vertical plate:
  // D / L < 35 / Gr_L^(1/4) (Incropera et al., Fundamentals of Heat and Mass
  // Transfer, Sec. 9.6.1), with D the outer wall diameter, L the bed height
  // and Gr_L at a wall 10 K above the ambient temperature
  // (kNaturalConvectionCheckDeltaT). The curvature of a slender cylinder
  // raises h above the plate value, so the model then under-predicts the heat
  // loss; the warning names the yang_vertical_cylinder kind, which holds at
  // any D/L. Not for an isothermal run.
  inline constexpr double kNaturalConvectionCheckDeltaT = 10.0;  // [K]
  std::vector<std::string> NaturalConvectionWarnings(const Bed& bed);

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
  // could not be checked (a component without an element formula).
  // With `order_notes`, also an informational note per reaction whose orders
  // differ from -nu (reactants) and 0 (products), i.e. a non-elementary rate.
  // Expressions must already be evaluated.
  std::vector<std::string> ReactionWarnings(const AdsorbentDirectory& adsorbents, const Bed& bed,
                                            bool order_notes = false);

  // Heats of adsorption taken from the isotherm (HeatOfAdsorptionMode::Isotherm,
  // objects/porous_media.h), for the adsorbents that some layer uses; none for
  // an isothermal run:
  // - one warning per component whose isotherm gives no isosteric heat, i.e.
  //   |q_st| < 1 J/mol at the reference state (T_ref, P_ref, pure component)
  //   and at the dilute limit: its isotherm does not depend on temperature,
  //   so the energy balance would release no heat for it;
  // - one per component whose isosteric heat is not finite at the reference
  //   state (the stored constant is used where that happens);
  // - one if any component uses the isotherm and the fluid's heat-of-adsorption
  //   convention is not isosteric: the convention applies to constant values
  //   only, and the isosteric heat enters as it is.
  std::vector<std::string> HeatOfAdsorptionWarnings(const AdsorbentDirectory& adsorbents, const Bed& bed,
                                                    const Fluid& fluid, double T_ref, double P_ref);

  // Supersaturated gas at the start of a run (condensation is not modelled;
  // see SUPERSATURATION in run_health.h): the relative saturation
  // y_i P / p_sat,i(T) (RelativeSaturation) of every component that can
  // condense at T, above 1, gives
  // - one warning for the initial gas (initial.y at initial.pressure), at the
  //   initial temperature or, if larger, at the initial wall temperature;
  // - one per inlet stream and per inlet pressure boundary that carries its
  //   own gas (not a connected stream or a buffer-tank inlet, which carries
  //   another endpoint's gas), at its inlet temperature (the larger of the
  //   step's start and end values of a temperature ramp; the isothermal
  //   temperature in an isothermal run). The pressure is the boundary's own
  //   (the end value of a ramp); a flow inlet enters at the bed pressure,
  //   taken as the step's pressure boundary at the other end (its end
  //   value), else the initial pressure.
  // Expressions must already be evaluated.
  std::vector<std::string> SupersaturationWarnings(const Cycle& cycle, const Bed& bed, const Fluid& fluid,
                                                   const InitialConditions& ic);

} // namespace adsim
