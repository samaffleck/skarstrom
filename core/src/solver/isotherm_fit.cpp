/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "adsim/solver/isotherm_fit.h"

#include <algorithm>
#include <cmath>

#include <Eigen/Dense>
#include <dlib/global_optimization.h>

#include "adsim/math/constants.h"

namespace adsim {
namespace isofit {

const std::vector<IsothermType>& FittableKinds() {
  using T = IsothermType;
  static const std::vector<T> kinds = {
      T::HENRY,
      T::SINGLE_SITE_LANGMUIR_P,
      T::DUAL_SITE_LANGMUIR_P,
      T::SINGLE_SITE_LANGMUIR_C,
      T::DUAL_SITE_LANGMUIR_C,
      T::EXTENDED_SINGLE_SITE_LANGMUIR_P,
      T::EXTENDED_DUAL_SITE_LANGMUIR_P,
      T::FREUNDLICH,
      T::SIPS,
      T::DUAL_SITE_SIPS,
      T::TOTH,
      T::SIPS_2,
      T::EXTENDED_SIPS_2,
      T::TOTH_2,
  };
  return kinds;
}
bool IsTemperatureDependent(IsothermType t) {
  using T = IsothermType;
  switch (t) {
    case T::HENRY:
    case T::SINGLE_SITE_LANGMUIR_P:
    case T::DUAL_SITE_LANGMUIR_P:
    case T::SINGLE_SITE_LANGMUIR_C:
    case T::DUAL_SITE_LANGMUIR_C:
    case T::EXTENDED_SINGLE_SITE_LANGMUIR_P:
    case T::EXTENDED_DUAL_SITE_LANGMUIR_P:
    case T::SIPS_2:
    case T::EXTENDED_SIPS_2:
    case T::TOTH_2:
      return true;
    default:
      return false;
  }
}
bool PointUsable(const Point& p) {
  return std::isfinite(p.P) && std::isfinite(p.q) && std::isfinite(p.T) && std::isfinite(p.y) &&
         p.P >= 0.0 && p.T > 0.0 && p.y > 0.0;
}
DataSummary Summarise(const std::vector<Point>& pts) {
  const double R = constants::kIdealGasConstant;
  DataSummary s;
  double t_sum = 0.0;
  std::vector<double> temps;
  bool any_pos = false, any_slope = false, any_c = false;
  for (const auto& p : pts) {
    if (!PointUsable(p)) continue;
    s.n++;
    t_sum += p.T;
    const double pp = p.P * p.y;
    const double c = pp / (R * p.T);
    s.q_max = std::max(s.q_max, p.q);
    s.p_max = std::max(s.p_max, pp);
    s.c_max = std::max(s.c_max, c);
    if (pp > 0.0) {
      s.p_min_pos = any_pos ? std::min(s.p_min_pos, pp) : pp;
      any_pos = true;
      s.c_min_pos = any_c ? std::min(s.c_min_pos, c) : c;
      any_c = true;
      if (p.q > 0.0) {
        const double sl = p.q / pp;
        s.slope_min = any_slope ? std::min(s.slope_min, sl) : sl;
        s.slope_max = any_slope ? std::max(s.slope_max, sl) : sl;
        any_slope = true;
      }
    }
    bool seen = false;
    for (double t : temps) if (std::abs(t - p.T) < 0.5) { seen = true; break; }
    if (!seen) temps.push_back(p.T);
  }
  s.n_temperatures = static_cast<int>(temps.size());
  s.t_mean = s.n > 0 ? t_sum / s.n : 0.0;
  return s;
}
double ReferenceTemperature(const DataSummary& s) {
  return std::round(s.t_mean * 100.0) / 100.0;
}
bool DefaultSpecs(IsothermType type, const std::vector<Point>& pts,
                  std::vector<ParamSpec>& out, double& t_ref, std::string& err) {
  using T = IsothermType;
  out.clear();
  const DataSummary s = Summarise(pts);
  if (s.n == 0) { err = "no usable points (need finite p >= 0, q, T > 0)"; return false; }
  if (!(s.q_max > 0.0)) { err = "all loadings are <= 0; nothing to fit"; return false; }
  if (!(s.p_min_pos > 0.0)) { err = "need at least one point with pressure > 0"; return false; }
  t_ref = ReferenceTemperature(s);
  const bool multi_t = s.n_temperatures >= 2;
  const char* single_t_note =
      "single temperature: not identifiable from these data, held fixed (enter a known value if you have one)";

  auto geo = [](double a, double b) { return std::sqrt(a * b); };
  auto capacity = [&](const char* name, double lo_frac) {
    ParamSpec p;
    p.name = name;
    p.scale = Scale::Log;
    p.lower = s.q_max * lo_frac;
    p.upper = s.q_max * 1e3;
    p.value = s.q_max * 1.2;
    out.push_back(p);
  };
  // Affinity [1/Pa] (or [m3/mol] for the concentration forms).
  auto affinity = [&](const char* name, const char* energy, bool conc) {
    ParamSpec p;
    p.name = name;
    p.scale = Scale::Log;
    const double xmax = conc ? s.c_max : s.p_max;
    const double xmin = conc ? s.c_min_pos : s.p_min_pos;
    p.lower = 1e-3 / xmax;
    p.upper = 1e3 / xmin;
    p.value = 1.0 / geo(xmin, xmax);
    if (energy) {
      p.pair_energy = energy;
      p.note = "bounds apply to the value at T_ref, " + std::string(name) + "*exp(" + energy + "/(R*T_ref))";
    }
    out.push_back(p);
  };
  auto energy = [&](const char* name) {
    ParamSpec p;
    p.name = name;
    p.scale = Scale::Linear;
    p.lower = 0.0;
    p.upper = 1.0e5;
    if (multi_t) {
      p.value = 2.0e4;
    } else {
      p.value = 0.0;
      p.fixed = true;
      p.note = single_t_note;
    }
    out.push_back(p);
  };
  auto exponent = [&](const char* name, double lo, double hi, double init) {
    ParamSpec p;
    p.name = name;
    p.scale = Scale::Linear;
    p.lower = lo;
    p.upper = hi;
    p.value = init;
    out.push_back(p);
  };
  auto tdep = [&](const char* name, double lo, double hi) {
    ParamSpec p;
    p.name = name;
    p.scale = Scale::Linear;
    p.lower = lo;
    p.upper = hi;
    p.value = 0.0;
    if (!multi_t) { p.fixed = true; p.note = single_t_note; }
    out.push_back(p);
  };
  auto reference_t = [&](const char* name) {
    ParamSpec p;
    p.name = name;
    p.scale = Scale::Linear;
    p.lower = p.upper = p.value = t_ref;
    p.fixed = true;
    p.note = "reference temperature (mean data temperature)";
    out.push_back(p);
  };

  switch (type) {
    case T::HENRY: {
      ParamSpec k;
      k.name = "k";
      k.scale = Scale::Log;
      k.lower = s.slope_min * 1e-3;
      k.upper = s.slope_max * 1e3;
      k.value = geo(s.slope_min, s.slope_max);
      k.pair_energy = "Q";
      k.note = "bounds apply to the value at T_ref, k*exp(Q/(R*T_ref))";
      out.push_back(k);
      energy("Q");
      break;
    }
    case T::SINGLE_SITE_LANGMUIR_P:
    case T::EXTENDED_SINGLE_SITE_LANGMUIR_P:
      capacity("m1", 0.5);
      affinity("b1", "q1", false);
      energy("q1");
      break;
    case T::SINGLE_SITE_LANGMUIR_C:
      capacity("m1", 0.5);
      affinity("b1", "q1", true);
      energy("q1");
      break;
    case T::DUAL_SITE_LANGMUIR_P:
    case T::EXTENDED_DUAL_SITE_LANGMUIR_P:
      capacity("m1", 1e-3);
      affinity("b1", "q1", false);
      energy("q1");
      capacity("m2", 1e-3);
      affinity("b2", "q2", false);
      energy("q2");
      break;
    case T::DUAL_SITE_LANGMUIR_C:
      capacity("m1", 1e-3);
      affinity("b1", "q1", true);
      energy("q1");
      capacity("m2", 1e-3);
      affinity("b2", "q2", true);
      energy("q2");
      break;
    case T::FREUNDLICH: {
      // q = (a*P)^(1/v)  =>  a = q^v / P; bound it over the exponent range.
      double lo = std::numeric_limits<double>::infinity(), hi = 0.0;
      for (const auto& p : pts) {
        if (!PointUsable(p) || p.P * p.y <= 0.0 || p.q <= 0.0) continue;
        for (double v : {0.2, 1.0, 5.0}) {
          const double a = std::pow(p.q, v) / (p.P * p.y);
          lo = std::min(lo, a);
          hi = std::max(hi, a);
        }
      }
      ParamSpec a;
      a.name = "a";
      a.scale = Scale::Log;
      a.lower = lo * 1e-3;
      a.upper = hi * 1e3;
      a.value = geo(lo, hi);
      out.push_back(a);
      exponent("v", 0.2, 5.0, 1.0);
      break;
    }
    case T::SIPS:
      capacity("q_sat", 0.5);
      affinity("b", nullptr, false);
      exponent("v", 0.2, 5.0, 1.0);
      break;
    case T::DUAL_SITE_SIPS:
      capacity("q_sat1", 1e-3);
      affinity("b1", nullptr, false);
      exponent("v1", 0.2, 5.0, 1.0);
      capacity("q_sat2", 1e-3);
      affinity("b2", nullptr, false);
      exponent("v2", 0.2, 5.0, 1.0);
      break;
    case T::TOTH:
      capacity("q_sat", 0.5);
      affinity("b", nullptr, false);
      exponent("v", 0.05, 2.0, 1.0);
      break;
    case T::SIPS_2:
    case T::EXTENDED_SIPS_2:
      capacity("A1", 0.5);
      tdep("A2", -10.0, 10.0);
      affinity("B1", "B2", false);
      energy("B2");
      exponent("C1", 0.2, 5.0, 1.0);
      tdep("C2", -10.0, 10.0);
      reference_t("T_ref");
      break;
    case T::TOTH_2:
      capacity("q_sat", 0.5);
      tdep("x", -20.0, 20.0);
      reference_t("T_0");
      affinity("b0", nullptr, false);  // b0 is already the affinity at T_0
      out.back().note = "affinity at T_0";
      energy("dH0");
      exponent("w0", 0.05, 2.0, 1.0);
      tdep("alpha", -5.0, 5.0);
      break;
    default:
      err = std::string("no default fit set-up for isotherm kind '") + IsothermKindString(type) +
            "'; supply every parameter with bounds";
      return false;
  }
  // Initial guesses inside their bounds.
  for (auto& p : out) {
    if (p.upper < p.lower) std::swap(p.lower, p.upper);
    p.value = std::clamp(p.value, p.lower, p.upper);
  }
  return true;
}

// ---------------------------------------------------------------------------
// Fit
// ---------------------------------------------------------------------------

namespace detail {

// Maps the free-parameter search vector z to model parameter values.
class Mapper {
 public:
  Mapper(Model& m, const std::vector<ParamSpec>& specs, double t_ref)
      : m_(m), specs_(specs), t_ref_(t_ref) {
    for (size_t i = 0; i < specs_.size(); ++i)
      if (!specs_[i].fixed) free_.push_back(i);
  }

  const std::vector<size_t>& free() const { return free_; }

  // Search-space coordinate of a natural-unit value for spec i.
  double ToZ(size_t i, double natural_eff) const {
    return specs_[i].scale == Scale::Log ? std::log10(natural_eff) : natural_eff;
  }
  double FromZ(size_t i, double z) const {
    return specs_[i].scale == Scale::Log ? std::pow(10.0, z) : z;
  }

  double LowerZ(size_t i) const { return ToZ(i, specs_[i].lower); }
  double UpperZ(size_t i) const { return ToZ(i, specs_[i].upper); }

  // Effective (T_ref) value of a paired parameter from its natural value.
  double ToEffective(size_t i, double natural) const {
    const auto& s = specs_[i];
    if (s.pair_energy.empty()) return natural;
    return natural * std::exp(EnergyValueSpec(s.pair_energy) / (R() * t_ref_));
  }

  // Apply a search vector: fixed values, then free non-paired, then paired.
  void Apply(const dlib::matrix<double, 0, 1>& z) {
    for (size_t i = 0; i < specs_.size(); ++i)
      if (specs_[i].fixed) m_.Set(specs_[i].name, specs_[i].value);
    for (size_t k = 0; k < free_.size(); ++k) {
      const size_t i = free_[k];
      if (specs_[i].pair_energy.empty()) m_.Set(specs_[i].name, FromZ(i, z(static_cast<long>(k))));
    }
    for (size_t k = 0; k < free_.size(); ++k) {
      const size_t i = free_[k];
      if (specs_[i].pair_energy.empty()) continue;
      const double eff = FromZ(i, z(static_cast<long>(k)));
      const double e = m_.Get(specs_[i].pair_energy);
      m_.Set(specs_[i].name, eff * std::exp(-e / (R() * t_ref_)));
    }
  }

  dlib::matrix<double, 0, 1> InitialZ() const {
    dlib::matrix<double, 0, 1> z(static_cast<long>(free_.size()));
    for (size_t k = 0; k < free_.size(); ++k) {
      const size_t i = free_[k];
      const double eff = ToEffective(i, specs_[i].value);
      z(static_cast<long>(k)) = std::clamp(ToZ(i, eff), LowerZ(i), UpperZ(i));
    }
    return z;
  }

  static double R() { return constants::kIdealGasConstant; }

 private:
  // The initial/held value of an energy parameter (for the initial guess).
  double EnergyValueSpec(const std::string& name) const {
    for (const auto& s : specs_) if (s.name == name) return s.value;
    return 0.0;
  }

  Model& m_;
  const std::vector<ParamSpec>& specs_;
  double t_ref_;
  std::vector<size_t> free_;
};

}  // namespace detail
FitResult Fit(IsothermType type, const std::vector<Point>& points_in,
              std::vector<ParamSpec> specs, double t_ref,
              const FitOptions& opt) {
  FitResult res;
  Model model(type);
  if (!model.valid()) { res.error = "cannot build an isotherm of this kind"; return res; }

  std::vector<Point> pts;
  for (const auto& p : points_in) if (PointUsable(p)) pts.push_back(p);
  const int n = static_cast<int>(pts.size());
  res.n = n;
  if (n == 0) { res.error = "no usable points (need finite p >= 0, q, T > 0)"; return res; }
  if (!(t_ref > 0.0)) t_ref = ReferenceTemperature(Summarise(pts));
  res.t_ref = t_ref;

  // Every model parameter must be specified exactly once; a parameter the
  // caller left out is held at the engine default.
  const auto names = model.Names();
  for (const auto& s : specs) {
    if (!model.Field(s.name)) { res.error = "unknown parameter '" + s.name + "' for this isotherm kind"; return res; }
  }
  for (const auto& nm : names) {
    bool have = false;
    for (const auto& s : specs) if (s.name == nm) { have = true; break; }
    if (!have) {
      ParamSpec p;
      p.name = nm;
      p.value = p.lower = p.upper = model.Get(nm);
      p.fixed = true;
      specs.push_back(p);
    }
  }
  for (auto& s : specs) {
    if (!std::isfinite(s.value)) { res.error = "parameter '" + s.name + "' has a non-finite value"; return res; }
    if (s.fixed) continue;
    if (!std::isfinite(s.lower) || !std::isfinite(s.upper)) { res.error = "parameter '" + s.name + "' needs finite bounds"; return res; }
    if (s.upper < s.lower) std::swap(s.lower, s.upper);
    if (s.scale == Scale::Log && !(s.lower > 0.0)) { res.error = "parameter '" + s.name + "' is searched in log space and needs bounds > 0"; return res; }
    if (!(s.upper > s.lower)) { s.fixed = true; s.value = s.lower; continue; }  // degenerate box => fixed
    if (!s.pair_energy.empty()) {
      bool ok = false;
      for (const auto& e : specs) if (e.name == s.pair_energy) { ok = true; break; }
      if (!ok) s.pair_energy.clear();
    }
  }

  detail::Mapper map(model, specs, t_ref);
  const auto& free = map.free();
  const long p = static_cast<long>(free.size());
  if (n < p) {
    res.error = "need at least " + std::to_string(p) + " points to fit " + std::to_string(p) + " free parameters";
    return res;
  }

  double q_mean = 0.0;
  for (const auto& pt : pts) q_mean += pt.q;
  q_mean /= n;
  double sst = 0.0;
  for (const auto& pt : pts) sst += (pt.q - q_mean) * (pt.q - q_mean);
  // Finite stand-in for a non-finite model (keeps dlib's model well-conditioned).
  const double bad = 1e6 * (sst + 1.0) + 1e6;

  // The search runs in the unit box u in [0,1]^p, z = lo + u*(hi - lo), so
  // every free parameter has the same scale for both stages. Every evaluation
  // updates the incumbent, so a stage that stops early never loses progress.
  dlib::matrix<double, 0, 1> lo(p), hi(p);
  for (long k = 0; k < p; ++k) {
    lo(k) = map.LowerZ(free[static_cast<size_t>(k)]);
    hi(k) = map.UpperZ(free[static_cast<size_t>(k)]);
  }
  auto to_z = [&](const Eigen::VectorXd& u) {
    dlib::matrix<double, 0, 1> z(p);
    for (long k = 0; k < p; ++k) z(k) = lo(k) + std::clamp(u(k), 0.0, 1.0) * (hi(k) - lo(k));
    return z;
  };

  int evals = 0;
  dlib::matrix<double, 0, 1> best_z = map.InitialZ();
  double best = std::numeric_limits<double>::infinity();
  std::vector<std::pair<double, Eigen::VectorXd>> history;  // (sse, u) of every evaluation
  Eigen::VectorXd r_buf(n);
  // Residuals q_model - q at u (into r); returns SSE (finite stand-in if not).
  auto residuals = [&](const Eigen::VectorXd& u, Eigen::VectorXd& r) -> double {
    ++evals;
    const auto z = to_z(u);
    map.Apply(z);
    double sse = 0.0;
    for (int i = 0; i < n; ++i) {
      r(i) = model.Eval(pts[static_cast<size_t>(i)].T, pts[static_cast<size_t>(i)].P, pts[static_cast<size_t>(i)].y) -
             pts[static_cast<size_t>(i)].q;
      sse += r(i) * r(i);
    }
    sse = std::isfinite(sse) ? std::min(sse, bad) : bad;
    if (sse < best) { best = sse; best_z = z; }
    return sse;
  };
  auto sse_u = [&](const Eigen::VectorXd& u) {
    const double v = residuals(u, r_buf);
    history.emplace_back(v, u);
    return v;
  };

  if (p > 0) {
    Eigen::VectorXd u0(p);
    const auto z0 = map.InitialZ();
    for (long k = 0; k < p; ++k) u0(k) = hi(k) > lo(k) ? (z0(k) - lo(k)) / (hi(k) - lo(k)) : 0.0;
    sse_u(u0);

    // 1. Global: dlib MaxLIPO + trust region over the unit box.
    const int n_global = opt.max_global_evals > 0 ? opt.max_global_evals : static_cast<int>(100 + 80 * p);
    dlib::matrix<double, 0, 1> ulo(p), uhi(p);
    for (long k = 0; k < p; ++k) { ulo(k) = 0.0; uhi(k) = 1.0; }
    try {
      (void)dlib::find_min_global(
          [&](const dlib::matrix<double, 0, 1>& x) {
            Eigen::VectorXd u(p);
            for (long k = 0; k < p; ++k) u(k) = x(k);
            return sse_u(u);
          },
          ulo, uhi, dlib::max_function_calls(static_cast<size_t>(n_global)));
    } catch (const std::exception& e) {
      res.warning += std::string("global search stopped early: ") + e.what() + ". ";
    }

    // 2. Local: bounded Levenberg-Marquardt (projected onto the box, forward-
    //    difference Jacobian) from the best few distinct points of the global
    //    stage. LM converges quadratically on the near-zero-residual problems
    //    typical of isotherm data, where a derivative-free polish stalls.
    std::sort(history.begin(), history.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    std::vector<Eigen::VectorXd> starts;
    for (const auto& [v, u] : history) {
      (void)v;
      bool distinct = true;
      for (const auto& s : starts) if ((s - u).cwiseAbs().maxCoeff() < 0.02) { distinct = false; break; }
      if (distinct) starts.push_back(u);
      if (starts.size() >= 6) break;
    }
    const int lm_budget = std::max(200, opt.max_local_evals);
    Eigen::VectorXd r(n), r_new(n), r_h(n);
    for (const auto& start : starts) {
      Eigen::VectorXd u = start;
      double f = residuals(u, r);
      double lambda = 1e-3;
      const int eval_cap = evals + lm_budget;
      for (int it = 0; it < 200 && evals < eval_cap && f > 0.0; ++it) {
        Eigen::MatrixXd J(n, p);
        for (long k = 0; k < p; ++k) {
          const double h = 1e-7;
          Eigen::VectorXd uh = u;
          const double step = (u(k) + h <= 1.0) ? h : -h;
          uh(k) += step;
          residuals(uh, r_h);
          J.col(k) = (r_h - r) / step;
        }
        const Eigen::MatrixXd A = J.transpose() * J;
        const Eigen::VectorXd g = J.transpose() * r;
        const double dfloor = std::max(1e-30, 1e-12 * A.diagonal().maxCoeff());
        bool improved = false;
        double f_new = f;
        Eigen::VectorXd u_new = u;
        while (lambda < 1e12) {
          Eigen::MatrixXd M = A;
          for (long k = 0; k < p; ++k) M(k, k) += lambda * std::max(A(k, k), dfloor);
          const Eigen::VectorXd delta = M.ldlt().solve(-g);
          if (!delta.allFinite()) { lambda *= 10.0; continue; }
          u_new = (u + delta).cwiseMax(0.0).cwiseMin(1.0);
          f_new = residuals(u_new, r_new);
          if (f_new < f) { improved = true; break; }
          lambda *= 4.0;
        }
        if (!improved) break;
        const double gain = f - f_new;
        const double move = (u_new - u).cwiseAbs().maxCoeff();
        u = u_new;
        r = r_new;
        f = f_new;
        lambda = std::max(lambda / 3.0, 1e-12);
        if (gain <= 1e-14 * f || move < 1e-13) break;
      }
    }
  }

  // Final state + statistics (best_z is the incumbent from every stage).
  const dlib::matrix<double, 0, 1> z_final = best_z;
  map.Apply(z_final);
  res.evaluations = evals;
  res.q_fit.resize(pts.size());
  double sse = 0.0;
  for (size_t i = 0; i < pts.size(); ++i) {
    res.q_fit[i] = model.Eval(pts[i].T, pts[i].P, pts[i].y);
    const double r = res.q_fit[i] - pts[i].q;
    sse += r * r;
  }
  if (!std::isfinite(sse)) { res.error = "the model is not finite at the best point found"; return res; }
  res.sse = sse;
  res.rmse = std::sqrt(sse / n);
  res.r2 = sst > 0.0 ? 1.0 - sse / sst : (sse == 0.0 ? 1.0 : 0.0);
  res.dof = n - static_cast<int>(p);

  // Linearised standard errors in natural units.
  std::vector<double> se(static_cast<size_t>(p), std::numeric_limits<double>::quiet_NaN());
  if (p > 0 && res.dof > 0) {
    std::vector<double> theta(static_cast<size_t>(p));
    for (long k = 0; k < p; ++k) theta[static_cast<size_t>(k)] = model.Get(specs[free[static_cast<size_t>(k)]].name);
    Eigen::MatrixXd J(n, p);
    bool finite = true;
    for (long k = 0; k < p; ++k) {
      const std::string& nm = specs[free[static_cast<size_t>(k)]].name;
      const double t0 = theta[static_cast<size_t>(k)];
      const double h = std::max(std::abs(t0) * 1e-6, 1e-12);
      model.Set(nm, t0 + h);
      std::vector<double> up(pts.size());
      for (size_t i = 0; i < pts.size(); ++i) up[i] = model.Eval(pts[i].T, pts[i].P, pts[i].y);
      model.Set(nm, t0 - h);
      for (size_t i = 0; i < pts.size(); ++i) {
        const double d = (up[i] - model.Eval(pts[i].T, pts[i].P, pts[i].y)) / (2.0 * h);
        J(static_cast<long>(i), k) = d;
        if (!std::isfinite(d)) finite = false;
      }
      model.Set(nm, t0);
    }
    if (finite) {
      Eigen::VectorXd scale(p);
      for (long k = 0; k < p; ++k) {
        const double nrm = J.col(k).norm();
        scale(k) = nrm > 0.0 ? 1.0 / nrm : 0.0;
      }
      bool any_zero = false;
      for (long k = 0; k < p; ++k) if (scale(k) == 0.0) any_zero = true;
      if (!any_zero) {
        const Eigen::MatrixXd Js = J * scale.asDiagonal();
        const Eigen::MatrixXd A = Js.transpose() * Js;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
        const double emax = es.eigenvalues().maxCoeff();
        const double emin = es.eigenvalues().minCoeff();
        if (es.info() == Eigen::Success && emin > emax * 1e-13 && emin > 0.0) {
          const Eigen::MatrixXd Ainv =
              es.eigenvectors() * es.eigenvalues().cwiseInverse().asDiagonal() * es.eigenvectors().transpose();
          const double s2 = sse / res.dof;
          for (long k = 0; k < p; ++k) se[static_cast<size_t>(k)] = std::sqrt(std::max(0.0, s2 * Ainv(k, k))) * scale(k);
        } else {
          res.warning += "Parameters are strongly correlated (J^T J near-singular); standard errors are not reported. ";
        }
      }
    }
  }

  // Bound hits.
  for (long k = 0; k < p; ++k) {
    const size_t i = free[static_cast<size_t>(k)];
    const double z = z_final(k), l = map.LowerZ(i), u = map.UpperZ(i);
    const double tol = 1e-6 * std::max(1.0, u - l);
    if (z - l < tol) res.warning += specs[i].name + " is at its lower bound; consider widening it. ";
    else if (u - z < tol) res.warning += specs[i].name + " is at its upper bound; consider widening it. ";
  }

  for (const auto& nm : names) {
    ParamResult pr;
    pr.name = nm;
    pr.value = model.Get(nm);
    for (size_t i = 0; i < specs.size(); ++i) {
      if (specs[i].name != nm) continue;
      pr.fixed = specs[i].fixed;
      for (long k = 0; k < p; ++k) if (free[static_cast<size_t>(k)] == i) pr.std_error = se[static_cast<size_t>(k)];
    }
    res.params.push_back(pr);
  }
  res.ok = true;
  return res;
}
bool Evaluate(IsothermType type, const std::vector<std::pair<std::string, double>>& values,
              const std::vector<Point>& pts, std::vector<double>& q, std::string& err) {
  Model m(type);
  if (!m.valid()) { err = "cannot build an isotherm of this kind"; return false; }
  for (const auto& [nm, v] : values) {
    if (!m.Set(nm, v)) { err = "unknown parameter '" + nm + "' for this isotherm kind"; return false; }
  }
  q.resize(pts.size());
  for (size_t i = 0; i < pts.size(); ++i) q[i] = m.Eval(pts[i].T, pts[i].P, pts[i].y);
  return true;
}

}  // namespace isofit
}  // namespace adsim
