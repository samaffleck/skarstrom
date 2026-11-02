/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Ideal Adsorbed Solution Theory (Myers & Prausnitz, AIChE J. 11 (1965) 121):
// the mixture loadings of an adsorbent from its components' PURE isotherms.
//
// At the cell temperature T, total pressure P and gas mole fractions y_i, the
// adsorbed phase is the ideal solution whose components all sit at the same
// reduced spreading pressure psi [mol/kg],
//
//   psi_i(P) = ∫_0^P q_i°(p) / p dp,          (pure isotherm q_i° at T)
//   P y_i    = x_i P_i°,   psi_i(P_i°) = psi for every i,   Σ x_i = 1,
//   1 / q_T  = Σ x_i / q_i°(P_i°),              q_i = x_i q_T.
//
// The unknown is the scalar psi (the fictitious pure pressures P_i° follow by
// inverting psi_i). psi is bracketed by max_i psi_i(y_i P) and max_i psi_i(Y P)
// (Y the adsorbing components' mole-fraction sum) and found by a bracketed
// Newton iteration; psi_i is inverted analytically where the kind allows and
// by Newton in ln P otherwise. IAST needs no mixture data and, unlike the
// extended Langmuir forms with unequal saturation capacities, is
// thermodynamically consistent; it assumes an ideal adsorbed phase.
//
// AdsorptionSource runs this per cell when its mixing rule is MixingRule::IAST
// (porous_media.h). The free functions below are also the public evaluation
// entry points (Python's skarstrom.evaluate_iast, front ends).

#include <string>
#include <vector>

#include "adsim/equations/isotherm.h"

namespace adsim {

  // One component's pure isotherm as the deck stores it: the kind and its raw
  // parameter values in IsothermFields order (isotherm_kinds.h). Resolved once
  // per update from the isotherm object, so the per-cell work does no map
  // lookups or casts. `kind` is INERT for a null isotherm or a component the
  // isotherm holds no parameters for.
  struct PureIsothermRef {
    IsothermType kind = IsothermType::INERT;
    double p[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  };
  PureIsothermRef ResolvePureIsotherm(IIsotherm* iso, const std::string& component);

  // The same isotherm with its temperature-dependent factors resolved at T [K],
  // so loading, spreading pressure and its inverse need no exp calls.
  struct PureIsotherm {
    IsothermType kind = IsothermType::INERT;
    double c[6] = {0, 0, 0, 0, 0, 0};  // meaning depends on kind (iast.cpp)
  };
  PureIsotherm PreparePureIsotherm(const PureIsothermRef& ref, double T);
  // Resolve and prepare in one call. Extended (competitive) kinds give their
  // pure-component limit.
  PureIsotherm PreparePureIsotherm(IIsotherm* iso, const std::string& component, double T);

  // Pure loading q°(P) [mol/kg] at partial pressure P [Pa] (yi = 1).
  double PureLoading(const PureIsotherm& iso, double P);

  // Reduced spreading pressure psi(P) = ∫_0^P q°(p)/p dp [mol/kg]. Closed form
  // for Henry, Langmuir and dual-site Langmuir (both bases), Freundlich, Sips,
  // dual-site Sips, Sips 2, Arrhenius Sips and anti-Langmuir; Gauss-Legendre
  // quadrature in ln p for Toth, Toth 2 and Redlich-Peterson. 0 for inert.
  double ReducedSpreadingPressure(const PureIsotherm& iso, double P);

  // The pressure P° [Pa] at which the pure isotherm reaches reduced spreading
  // pressure psi (the inverse of ReducedSpreadingPressure). `guess` (> 0) seeds
  // the Newton iteration of the kinds without a closed-form inverse.
  double InvertSpreadingPressure(const PureIsotherm& iso, double psi, double guess = 0.0);

  // lim_{P->0} q°/P [mol/(kg Pa)], the Henry constant; 0 for an inert and for
  // an isotherm with an infinite initial slope (Freundlich or Sips with n < 1).
  double DiluteSlope(const PureIsotherm& iso);

  // Same as the three above, straight from an isotherm object (convenience for
  // tests and tools; the solver prepares once per cell and reuses).
  double ReducedSpreadingPressure(IIsotherm* iso, const std::string& component, double T, double P);

  struct IastResult {
    bool converged = false;
    int iterations = 0;   // outer (psi) Newton iterations
    double psi = 0.0;     // reduced spreading pressure of the mixture [mol/kg]
  };

  // Scratch space reused across cells (no allocation inside the solve).
  struct IastScratch {
    std::vector<double> p0;    // fictitious pure pressures P_i°
    std::vector<double> q0;    // q_i°(P_i°)
    std::vector<double> yv;    // clamped mole fractions
  };

  // Solve IAST for `n` components at total pressure P [Pa] with gas mole
  // fractions y[i]: q[i] receives the mixture loading [mol/kg] of isos[i].
  // An inert isotherm gives q[i] = 0 and P <= 0 gives all zeros. A component
  // with y[i] <= 0 gets the infinite-dilution limit q_i = y_i P q_T / P_i°
  // (y_i P times the Henry constant when nothing else adsorbs), linear in y_i,
  // so the loading is smooth through y_i = 0 and negative below it, as the
  // extended isotherms are.
  // A single adsorbing component reduces to its pure isotherm at y P.
  // `psi_guess` > 0 warm-starts the iteration (the previous solve of the same
  // cell). On failure (non-finite inputs) q is NaN and converged is false.
  IastResult SolveIast(const PureIsotherm* isos, int n, double P, const double* y, double* q,
                       IastScratch& scratch, double psi_guess = 0.0);

  // Convenience: IAST loadings of several isotherm objects at (T, P, y).
  // isos[i] / components[i] name each pure isotherm; y and q have n entries.
  IastResult EvalIastMixture(const std::vector<IIsotherm*>& isos,
                             const std::vector<std::string>& components, double T, double P,
                             const std::vector<double>& y, std::vector<double>& q);

}  // namespace adsim
