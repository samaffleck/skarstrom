/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/equations/iast.h"
#include "adsim/equations/isotherm_kinds.h"
#include "adsim/math/constants.h"

#include <array>
#include <cmath>
#include <limits>

namespace adsim {

  namespace {

    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    constexpr double kInf = std::numeric_limits<double>::infinity();

    // PureIsotherm::c layout by kind (every pressure is the pure partial pressure P [Pa]):
    //   HENRY                                  q = K P                                c0 = K
    //   Langmuir family (P or C basis,         q = m B P / (1 + B P)                  c0 = m, c1 = B
    //     extended or not)                     B = b exp(q/RT) [/ (RT) on the C basis]
    //   dual-site Langmuir family              q = sum of two Langmuir terms          c0 = m1, c1 = B1, c2 = m2, c3 = B2
    //   FREUNDLICH                             q = (a P)^n                            c0 = a, c1 = n = 1/v
    //   SIPS, SIPS_2, EXTENDED_SIPS_2,         q = qs z / (1 + z), z = (b P)^n        c0 = qs, c1 = b, c2 = n
    //     ARRHENIUS_SIPS
    //   DUAL_SITE_SIPS                         q = sum of two Sips terms              c0 = qs1, c1 = b1, c2 = n1, c3 = qs2, c4 = b2, c5 = n2
    //   ANTI_LANGMUIR                          q = a P / (1 - b P)                    c0 = a, c1 = b
    //   REDLICH_PETERSON                       q = a P / (1 + (b P)^v)                c0 = a, c1 = b, c2 = v
    //   TOTH, TOTH_2                           q = qs b P / (1 + (b P)^v)^(1/v)       c0 = qs, c1 = b, c2 = v
    //   INERT                                  q = 0

    enum class Family { Inert, Henry, Langmuir, DualSiteLangmuir, Freundlich, Sips, DualSiteSips,
                        AntiLangmuir, RedlichPeterson, Toth };

    Family FamilyOf(IsothermType t) {
      using T = IsothermType;
      switch (t) {
        case T::INERT: return Family::Inert;
        case T::HENRY: return Family::Henry;
        case T::SINGLE_SITE_LANGMUIR_P:
        case T::SINGLE_SITE_LANGMUIR_C:
        case T::EXTENDED_SINGLE_SITE_LANGMUIR_P:
        case T::EXTENDED_SINGLE_SITE_LANGMUIR_C: return Family::Langmuir;
        case T::DUAL_SITE_LANGMUIR_P:
        case T::DUAL_SITE_LANGMUIR_C:
        case T::EXTENDED_DUAL_SITE_LANGMUIR_P:
        case T::EXTENDED_DUAL_SITE_LANGMUIR_C: return Family::DualSiteLangmuir;
        case T::FREUNDLICH: return Family::Freundlich;
        case T::SIPS:
        case T::SIPS_2:
        case T::EXTENDED_SIPS_2:
        case T::ARRHENIUS_SIPS: return Family::Sips;
        case T::DUAL_SITE_SIPS: return Family::DualSiteSips;
        case T::ANTI_LANGMUIR: return Family::AntiLangmuir;
        case T::REDLICH_PETERSON: return Family::RedlichPeterson;
        case T::TOTH:
        case T::TOTH_2: return Family::Toth;
      }
      return Family::Inert;
    }

    // ---- Gauss-Legendre quadrature for the kinds without a closed-form integral ----
    //
    // psi(P) = ∫_0^P q(p)/p dp = ∫_{-inf}^{ln P} q(e^u) du. In u the integrand is
    // smooth (the (b p)^v terms that have an infinite slope at p = 0 become
    // exponentials), so a Gauss rule on the window [ln P - L, ln P] converges
    // fast. Below the window the isotherm is in its Henry regime, q ≈ K p, whose
    // integral from -inf is q(p_min) exactly; at p_min = P e^{-L} = 1e-10 P the
    // tail is negligible anyway.

    constexpr int kGaussN = 48;
    constexpr double kLogWindow = 23.025850929940457;  // ln(1e10)

    struct GaussRule {
      std::array<double, kGaussN> x{}, w{};
      GaussRule() {
        // Nodes by Newton on P_n(x) from the Chebyshev-like initial guesses
        // (Numerical Recipes gauleg), symmetric about 0.
        const int n = kGaussN;
        const int m = (n + 1) / 2;
        for (int i = 0; i < m; ++i) {
          double z = std::cos(M_PI * (i + 0.75) / (n + 0.5));
          double pp = 0.0;
          for (int it = 0; it < 100; ++it) {
            double p1 = 1.0, p2 = 0.0;
            for (int j = 0; j < n; ++j) {
              const double p3 = p2;
              p2 = p1;
              p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1.0);
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            const double z1 = z;
            z = z1 - p1 / pp;
            if (std::abs(z - z1) < 1e-15) break;
          }
          x[i] = -z;
          x[n - 1 - i] = z;
          w[i] = w[n - 1 - i] = 2.0 / ((1.0 - z * z) * pp * pp);
        }
      }
    };

    const GaussRule& Gauss() {
      static const GaussRule rule;
      return rule;
    }

    double SpreadingPressureByQuadrature(const PureIsotherm& iso, double P) {
      const auto& g = Gauss();
      const double u_max = std::log(P);
      const double u_min = u_max - kLogWindow;
      const double half = 0.5 * (u_max - u_min);
      const double mid = 0.5 * (u_max + u_min);
      double sum = 0.0;
      for (int i = 0; i < kGaussN; ++i) sum += g.w[i] * PureLoading(iso, std::exp(mid + half * g.x[i]));
      return half * sum + PureLoading(iso, std::exp(u_min));
    }

    // Henry's constant lim_{p->0} q/p where the kind has one (the start of the
    // numeric inversion); a finite-pressure secant otherwise.
    double HenryConstant(const PureIsotherm& iso) {
      switch (FamilyOf(iso.kind)) {
        case Family::DualSiteLangmuir: return iso.c[0] * iso.c[1] + iso.c[2] * iso.c[3];
        case Family::RedlichPeterson: return iso.c[0];
        case Family::Toth: return iso.c[0] * iso.c[1];
        default: {
          const double p = 1e-3;
          return PureLoading(iso, p) / p;
        }
      }
    }

    // P° with psi(P°) = psi by Newton in s = ln P: f(s) = psi(e^s) - psi is
    // increasing and convex (f' = q(e^s) > 0, f'' = p q'(p) > 0), so from any
    // start the first step lands at or right of the root and the iteration then
    // descends monotonically. Steps are capped so a start far to the left (the
    // Henry estimate is a lower bound) cannot overflow.
    double InvertByNewton(const PureIsotherm& iso, double psi, double guess) {
      double s;
      if (guess > 0.0 && std::isfinite(guess)) {
        s = std::log(guess);
      } else {
        const double K = HenryConstant(iso);
        s = (K > 0.0 && std::isfinite(K)) ? std::log(std::max(psi / K, 1e-300)) : 0.0;
      }
      s = std::min(s, 700.0);
      double best_P = kNaN, best_err = kInf;
      for (int it = 0; it < 200; ++it) {
        const double P = std::exp(s);
        const double f = ReducedSpreadingPressure(iso, P) - psi;
        if (!std::isfinite(f)) break;
        const double err = std::abs(f) / psi;
        if (err < best_err) { best_err = err; best_P = P; }
        if (err <= 1e-14) return P;
        const double q = PureLoading(iso, P);
        if (!(q > 0.0)) break;
        double ds = -f / q;
        if (ds > 10.0) ds = 10.0;
        if (ds < -10.0) ds = -10.0;
        s = std::min(s + ds, 700.0);
        if (std::abs(ds) < 1e-15) break;
      }
      // Roundoff in the quadrature can keep |f| just above the target: accept
      // the best iterate when it is converged for every practical purpose.
      return (best_err <= 1e-9) ? best_P : kNaN;
    }

  }  // namespace

  PureIsothermRef ResolvePureIsotherm(IIsotherm* iso, const std::string& component) {
    PureIsothermRef ref;
    if (!iso) return ref;
    const IsothermType kind = iso->GetIsothermType();
    if (kind == IsothermType::INERT) return ref;
    const auto fields = IsothermFields(iso, component);
    if (fields.empty()) return ref;  // no parameters for this component: inert
    ref.kind = kind;
    for (size_t i = 0; i < fields.size() && i < 8; ++i) ref.p[i] = fields[i].second->Get();
    return ref;
  }

  PureIsotherm PreparePureIsotherm(const PureIsothermRef& ref, double T) {
    using K = IsothermType;
    const double R = constants::kIdealGasConstant;
    const double RT = R * T;
    const double* p = ref.p;
    PureIsotherm iso;
    iso.kind = ref.kind;
    double* c = iso.c;
    switch (ref.kind) {
      case K::INERT:
        break;
      case K::HENRY:  // k, Q
        c[0] = p[0] * std::exp(p[1] / RT);
        break;
      case K::SINGLE_SITE_LANGMUIR_P:
      case K::EXTENDED_SINGLE_SITE_LANGMUIR_P:  // m1, b1, q1
        c[0] = p[0];
        c[1] = p[1] * std::exp(p[2] / RT);
        break;
      case K::SINGLE_SITE_LANGMUIR_C:
      case K::EXTENDED_SINGLE_SITE_LANGMUIR_C:  // m1, b1, q1 on c = P/(RT)
        c[0] = p[0];
        c[1] = p[1] * std::exp(p[2] / RT) / RT;
        break;
      case K::DUAL_SITE_LANGMUIR_P:
      case K::EXTENDED_DUAL_SITE_LANGMUIR_P:  // m1, b1, q1, m2, b2, q2
        c[0] = p[0];
        c[1] = p[1] * std::exp(p[2] / RT);
        c[2] = p[3];
        c[3] = p[4] * std::exp(p[5] / RT);
        break;
      case K::DUAL_SITE_LANGMUIR_C:
      case K::EXTENDED_DUAL_SITE_LANGMUIR_C:
        c[0] = p[0];
        c[1] = p[1] * std::exp(p[2] / RT) / RT;
        c[2] = p[3];
        c[3] = p[4] * std::exp(p[5] / RT) / RT;
        break;
      case K::FREUNDLICH:  // a, v
        c[0] = p[0];
        c[1] = 1.0 / p[1];
        break;
      case K::SIPS:  // q_sat, b, v
      case K::TOTH:  // q_sat, b, v (c2 = v itself for Toth)
        c[0] = p[0];
        c[1] = p[1];
        c[2] = (ref.kind == K::SIPS) ? 1.0 / p[2] : p[2];
        break;
      case K::SIPS_2:
      case K::EXTENDED_SIPS_2: {  // A1, A2, B1, B2, C1, C2, T_ref
        const double T_ref = p[6];
        c[0] = p[0] * std::exp(p[1] * ((T / T_ref) - 1.0));
        c[1] = p[2] * std::exp(p[3] / RT);
        c[2] = p[4] + p[5] * ((T / T_ref) - 1.0);
        break;
      }
      case K::ARRHENIUS_SIPS:  // n_inf_a, n_inf_b, k_A, k_B, c_ref, c_amp, c_beta, c_Tref
        c[0] = p[0] * std::exp(p[1] / RT);
        c[1] = p[2] * std::exp(p[3] / RT);
        c[2] = p[4] + p[5] * (2.0 / M_PI) * std::atan(p[6] * (T - p[7]));
        break;
      case K::DUAL_SITE_SIPS:  // q_sat1, b1, v1, q_sat2, b2, v2
        c[0] = p[0]; c[1] = p[1]; c[2] = 1.0 / p[2];
        c[3] = p[3]; c[4] = p[4]; c[5] = 1.0 / p[5];
        break;
      case K::ANTI_LANGMUIR:  // a, b
        c[0] = p[0];
        c[1] = p[1];
        break;
      case K::REDLICH_PETERSON:  // a, b, v
        c[0] = p[0];
        c[1] = p[1];
        c[2] = p[2];
        break;
      case K::TOTH_2: {  // q_sat, x, T_0, b0, dH0, w0, alpha
        const double T_0 = p[2];
        c[0] = p[0] * std::exp(p[1] * (1.0 - T / T_0));
        c[1] = p[3] * std::exp(p[4] / (R * T_0) * (T_0 / T - 1.0));
        c[2] = p[5] + p[6] * (1.0 - T_0 / T);
        break;
      }
    }
    return iso;
  }

  PureIsotherm PreparePureIsotherm(IIsotherm* iso, const std::string& component, double T) {
    return PreparePureIsotherm(ResolvePureIsotherm(iso, component), T);
  }

  double PureLoading(const PureIsotherm& iso, double P) {
    const double* c = iso.c;
    switch (FamilyOf(iso.kind)) {
      case Family::Inert: return 0.0;
      case Family::Henry: return c[0] * P;
      case Family::Langmuir: {
        const double BP = c[1] * P;
        return c[0] * BP / (1.0 + BP);
      }
      case Family::DualSiteLangmuir: {
        const double B1P = c[1] * P, B2P = c[3] * P;
        return c[0] * B1P / (1.0 + B1P) + c[2] * B2P / (1.0 + B2P);
      }
      case Family::Freundlich: return std::pow(c[0] * P, c[1]);
      case Family::Sips: {
        const double z = std::pow(c[1] * P, c[2]);
        return c[0] * z / (1.0 + z);
      }
      case Family::DualSiteSips: {
        const double z1 = std::pow(c[1] * P, c[2]), z2 = std::pow(c[4] * P, c[5]);
        return c[0] * z1 / (1.0 + z1) + c[3] * z2 / (1.0 + z2);
      }
      case Family::AntiLangmuir: return c[0] * P / (1.0 - c[1] * P);
      case Family::RedlichPeterson: return c[0] * P / (1.0 + std::pow(c[1] * P, c[2]));
      case Family::Toth: {
        const double bp = c[1] * P;
        return c[0] * bp / std::pow(1.0 + std::pow(bp, c[2]), 1.0 / c[2]);
      }
    }
    return 0.0;
  }

  double ReducedSpreadingPressure(const PureIsotherm& iso, double P) {
    if (std::isnan(P)) return kNaN;
    if (!(P > 0.0)) return 0.0;
    const double* c = iso.c;
    switch (FamilyOf(iso.kind)) {
      case Family::Inert: return 0.0;
      case Family::Henry: return c[0] * P;
      case Family::Langmuir: return c[0] * std::log1p(c[1] * P);
      case Family::DualSiteLangmuir: return c[0] * std::log1p(c[1] * P) + c[2] * std::log1p(c[3] * P);
      case Family::Freundlich: return std::pow(c[0] * P, c[1]) / c[1];
      case Family::Sips: return c[0] / c[2] * std::log1p(std::pow(c[1] * P, c[2]));
      case Family::DualSiteSips:
        return c[0] / c[2] * std::log1p(std::pow(c[1] * P, c[2])) +
               c[3] / c[5] * std::log1p(std::pow(c[4] * P, c[5]));
      case Family::AntiLangmuir: {
        if (c[1] == 0.0) return c[0] * P;
        const double bP = c[1] * P;
        if (bP >= 1.0) return kInf;
        return -(c[0] / c[1]) * std::log1p(-bP);
      }
      case Family::RedlichPeterson:
      case Family::Toth:
        return SpreadingPressureByQuadrature(iso, P);
    }
    return 0.0;
  }

  double InvertSpreadingPressure(const PureIsotherm& iso, double psi, double guess) {
    if (std::isnan(psi)) return kNaN;
    if (!(psi > 0.0)) return 0.0;
    const double* c = iso.c;
    switch (FamilyOf(iso.kind)) {
      case Family::Inert: return kInf;  // never reached at any pressure
      case Family::Henry: return psi / c[0];
      case Family::Langmuir: return std::expm1(psi / c[0]) / c[1];
      case Family::Freundlich: return std::pow(c[1] * psi, 1.0 / c[1]) / c[0];
      case Family::Sips: return std::pow(std::expm1(c[2] * psi / c[0]), 1.0 / c[2]) / c[1];
      case Family::AntiLangmuir:
        return (c[1] == 0.0) ? psi / c[0] : -std::expm1(-c[1] * psi / c[0]) / c[1];
      case Family::DualSiteLangmuir:
      case Family::DualSiteSips:
      case Family::RedlichPeterson:
      case Family::Toth:
        return InvertByNewton(iso, psi, guess);
    }
    return kNaN;
  }

  double ReducedSpreadingPressure(IIsotherm* iso, const std::string& component, double T, double P) {
    return ReducedSpreadingPressure(PreparePureIsotherm(iso, component, T), P);
  }

  double DiluteSlope(const PureIsotherm& iso) {
    const double* c = iso.c;
    // lim q/p of scale * p^n as p -> 0: the scale for n = 1, 0 for n > 1, and
    // infinite for n < 1, reported as 0 (no linear extension below y = 0).
    auto power_law = [](double scale, double n) { return n == 1.0 ? scale : 0.0; };
    switch (FamilyOf(iso.kind)) {
      case Family::Inert: return 0.0;
      case Family::Henry: return c[0];
      case Family::Langmuir: return c[0] * c[1];
      case Family::DualSiteLangmuir: return c[0] * c[1] + c[2] * c[3];
      case Family::Freundlich: return power_law(c[0], c[1]);
      case Family::Sips: return power_law(c[0] * c[1], c[2]);
      case Family::DualSiteSips: return power_law(c[0] * c[1], c[2]) + power_law(c[3] * c[4], c[5]);
      case Family::AntiLangmuir: return c[0];
      case Family::RedlichPeterson: return c[0];
      case Family::Toth: return c[0] * c[1];
    }
    return 0.0;
  }

  IastResult SolveIast(const PureIsotherm* isos, int n, double P, const double* y, double* q,
                       IastScratch& scratch, double psi_guess) {
    IastResult result;
    if (n <= 0) { result.converged = true; return result; }
    if (scratch.p0.size() != static_cast<size_t>(n)) scratch.p0.assign(static_cast<size_t>(n), 0.0);
    scratch.q0.assign(static_cast<size_t>(n), 0.0);
    scratch.yv.assign(static_cast<size_t>(n), 0.0);
    double* p0 = scratch.p0.data();
    double* q0 = scratch.q0.data();
    double* yv = scratch.yv.data();

    bool bad = !std::isfinite(P);
    for (int i = 0; i < n; ++i) if (!std::isfinite(y[i])) bad = true;
    if (bad) {
      for (int i = 0; i < n; ++i) q[i] = kNaN;
      return result;
    }
    for (int i = 0; i < n; ++i) q[i] = 0.0;
    if (!(P > 0.0)) { result.converged = true; return result; }

    // Components that adsorb at all (not inert, a non-zero capacity: judged at
    // the total pressure so the choice never depends on the composition). Of
    // these, the ones in the gas (y > 0) set the spreading pressure; the rest
    // (y <= 0, which the integrator produces ahead of a front) get the
    // infinite-dilution limit of IAST, q_i = y_i P q_T / P_i°(psi), linear in
    // y_i. That is the first-order limit of the full solve as y_i -> 0+, so the
    // loading is continuous with a continuous slope through y_i = 0, as the
    // extended isotherms (linear in y_i) are.
    int count = 0, last = -1;
    double Y = 0.0;
    for (int i = 0; i < n; ++i) {
      // q°(P) > 0 exactly when psi(P) > 0, and costs one evaluation instead
      // of a quadrature for the Toth-type kinds.
      const bool adsorbs = isos[i].kind != IsothermType::INERT && PureLoading(isos[i], P) > 0.0;
      yv[i] = adsorbs ? y[i] : 0.0;
      if (adsorbs && y[i] > 0.0) { ++count; last = i; Y += y[i]; }
    }
    auto dilute = [&](double psi, double qT) {
      for (int i = 0; i < n; ++i) {
        if (yv[i] >= 0.0) continue;
        if (psi > 0.0) {
          const double Pi = InvertSpreadingPressure(isos[i], psi, p0[i]);
          if (Pi > 0.0 && std::isfinite(Pi)) { p0[i] = Pi; q[i] = yv[i] * P * qT / Pi; }
        } else {
          q[i] = yv[i] * P * DiluteSlope(isos[i]);
        }
      }
    };
    if (count == 0) {
      dilute(0.0, 0.0);
      result.converged = true;
      return result;
    }
    if (count == 1) {
      // Raoult's analogue with x = 1: the pure isotherm at the partial pressure.
      q[last] = PureLoading(isos[last], yv[last] * P);
      result.psi = ReducedSpreadingPressure(isos[last], yv[last] * P);
      result.converged = std::isfinite(q[last]);
      dilute(result.psi, q[last]);
      return result;
    }

    // Bracket: x_i <= 1 gives P_i° >= y_i P, so psi >= psi_i(y_i P) for every i;
    // some component has x_i >= y_i / Y (both sets sum to 1), so for it
    // P_i° <= Y P and psi <= psi_i(Y P). hi is infinite only for an
    // anti-Langmuir isotherm past its pole; the search then expands upwards.
    double lo = 0.0, hi = 0.0;
    for (int i = 0; i < n; ++i) {
      if (!(yv[i] > 0.0)) continue;
      lo = std::max(lo, ReducedSpreadingPressure(isos[i], yv[i] * P));
      hi = std::max(hi, ReducedSpreadingPressure(isos[i], Y * P));
    }
    if (!std::isfinite(lo) || std::isnan(hi)) {
      for (int i = 0; i < n; ++i) q[i] = kNaN;
      return result;
    }
    if (hi < lo) hi = lo;

    // S(psi) = sum_i x_i with x_i = y_i P / P_i°(psi). H = ln S is decreasing and
    // convex in psi (ln of a sum of exponentials of convex functions) and close
    // to linear, so Newton on H from the left of the root (S >= 1) converges
    // monotonically and fast; a start right of the root steps left at once.
    auto evaluate = [&](double psi, double& S, double& dS) {
      S = 0.0;
      dS = 0.0;
      for (int i = 0; i < n; ++i) {
        if (!(yv[i] > 0.0)) continue;
        const double Pi = InvertSpreadingPressure(isos[i], psi, p0[i]);
        const double qi = PureLoading(isos[i], Pi);
        if (!(Pi > 0.0) || !(qi > 0.0) || !std::isfinite(Pi)) { S = kNaN; return; }
        p0[i] = Pi;
        q0[i] = qi;
        const double x = yv[i] * P / Pi;
        S += x;
        dS -= x / qi;  // d/dpsi (y P / P°) = -y P / (P° q°)
      }
    };
    auto inside = [&](double lo_, double hi_) {
      return std::isfinite(hi_) ? 0.5 * (lo_ + hi_) : 2.0 * lo_ + 1.0;
    };

    double psi = (psi_guess > 0.0 && std::isfinite(psi_guess)) ? std::min(std::max(psi_guess, lo), hi) : lo;
    if (!(psi > 0.0)) psi = inside(lo, hi);
    double S = 0.0, dS = 0.0;
    for (int it = 0; it < 100; ++it) {
      ++result.iterations;
      evaluate(psi, S, dS);
      if (std::isnan(S)) break;
      const double H = std::log(S);
      if (std::abs(H) <= 1e-14) { result.converged = true; break; }
      if (H > 0.0) lo = psi; else hi = psi;
      double next = (dS < 0.0 && std::isfinite(dS)) ? psi - H * S / dS : kNaN;
      if (!(next >= lo && next <= hi) || next == psi) {
        if (std::isfinite(hi) && hi - lo <= 4e-16 * hi) {
          // The root is pinned to the last bit: S is 1 to within roundoff.
          result.converged = std::abs(H) < 1e-12;
          break;
        }
        next = inside(lo, hi);
      } else if (std::abs(next - psi) <= 1e-15 * psi && std::abs(H) < 1e-12) {
        psi = next;
        evaluate(psi, S, dS);
        result.converged = std::isfinite(S);
        ++result.iterations;
        break;
      }
      psi = next;
    }
    result.psi = psi;
    if (!result.converged) {
      for (int i = 0; i < n; ++i) q[i] = kNaN;
      return result;
    }

    // Adsorbed-phase composition (normalised so it sums to 1 exactly, which
    // removes the last-bit residual of the iteration) and the total loading.
    double inv_qT = 0.0;
    for (int i = 0; i < n; ++i) {
      if (!(yv[i] > 0.0)) continue;
      inv_qT += (yv[i] * P / p0[i]) / S / q0[i];
    }
    const double qT = 1.0 / inv_qT;
    for (int i = 0; i < n; ++i) {
      if (!(yv[i] > 0.0)) continue;
      q[i] = (yv[i] * P / p0[i]) / S * qT;
    }
    dilute(psi, qT);
    return result;
  }

  IastResult EvalIastMixture(const std::vector<IIsotherm*>& isos, const std::vector<std::string>& components,
                             double T, double P, const std::vector<double>& y, std::vector<double>& q) {
    const size_t n = isos.size();
    std::vector<PureIsotherm> pure(n);
    for (size_t i = 0; i < n; ++i) pure[i] = PreparePureIsotherm(isos[i], components[i], T);
    q.assign(n, 0.0);
    std::vector<double> yv(n, 0.0);
    for (size_t i = 0; i < n && i < y.size(); ++i) yv[i] = y[i];
    IastScratch scratch;
    return SolveIast(pure.data(), static_cast<int>(n), P, yv.data(), q.data(), scratch);
  }

}  // namespace adsim
