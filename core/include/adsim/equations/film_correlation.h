/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Gas-particle film correlations of a packed layer. One correlation per layer
// (MixedAdsorbentLayer::film_correlation_) gives both the Nusselt number of
// the particle heat-transfer film model and the Sherwood number of the
// Haynes mass-transfer coefficient (its k_f), by the heat-mass analogy: the
// same function with Pr gives Nu and with Sc gives Sh. A structured layer
// does not use it (its channels have their own Nu = Sh; porous_media.h,
// ChannelFilmNumber). Sources and validity: docs/film-correlations.md.

#include <cmath>
#include <string>
#include <vector>

namespace adsim {

  // Stored values; never renumber.
  enum class FilmCorrelation {
    WakaoFunazkri = 0,  // default
    Gunn = 1,
    RanzMarshall = 2
  };

  // Kind names in enum order (deck builder, parameter registry).
  inline const std::vector<std::string>& FilmCorrelationKinds() {
    static const std::vector<std::string> k{"wakao_funazkri", "gunn", "ranz_marshall"};
    return k;
  }

  // Wakao and Funazkri (1978), mass, and Wakao, Kaguei and Funazkri (1979),
  // heat: N = 2 + 1.1 X^(1/3) Re^0.6, Re = rho u d_p / mu on the superficial
  // velocity. The cube root is evaluated as pow(X, 0.3333333), the form the
  // engine has always used for Sh, so decks keep their results bit for bit
  // (it differs from X^(1/3) by a relative 3.3e-8 ln X, below 1e-7 for gases).
  inline double WakaoFunazkriNumber(double re, double x) {
    return 2.0 + 1.1 * std::pow(x, 0.3333333) * std::pow(re, 0.6);
  }

  // Gunn (1978): N = (7 - 10 e + 5 e^2)(1 + 0.7 Re^0.2 X^(1/3))
  //                + (1.33 - 2.4 e + 1.2 e^2) Re^0.7 X^(1/3),
  // with e the bed void fraction and Re on the superficial velocity and d_p
  // (as the engine's Re). N = 2 at e = 1, Re = 0 (an isolated sphere).
  inline double GunnNumber(double re, double x, double eps) {
    const double x3 = std::cbrt(x);
    const double a = 7.0 - 10.0 * eps + 5.0 * eps * eps;
    const double b = 1.33 - 2.4 * eps + 1.2 * eps * eps;
    return a * (1.0 + 0.7 * std::pow(re, 0.2) * x3) + b * std::pow(re, 0.7) * x3;
  }

  // Ranz and Marshall (1952), a single sphere: N = 2 + 0.6 Re^(1/2) X^(1/3).
  // Not a packed-bed correlation (it under-predicts in a bed at moderate
  // Re); offered to reproduce papers that use it, with a deck warning
  // (FilmCorrelationWarnings, deck_warnings.h). Re on the superficial velocity.
  inline double RanzMarshallNumber(double re, double x) {
    return 2.0 + 0.6 * std::sqrt(re) * std::cbrt(x);
  }

  // Nu (x = Pr) or Sh (x = Sc) of a packed layer by `c`; `re` >= 0 is the
  // particle Reynolds number on the superficial velocity and d_p, `eps` the
  // bed void fraction (read by Gunn only).
  inline double PackedFilmNumber(FilmCorrelation c, double re, double x, double eps) {
    switch (c) {
    case FilmCorrelation::Gunn: return GunnNumber(re, x, eps);
    case FilmCorrelation::RanzMarshall: return RanzMarshallNumber(re, x);
    case FilmCorrelation::WakaoFunazkri:
    default: return WakaoFunazkriNumber(re, x);
    }
  }

} // namespace adsim
