/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Cereal includes
#include <cereal/access.hpp>

#include "adsim/math/expression_param.h"


namespace adsim {

  class Polynomial5 {
  public:
    ExpressionParam a{0.0};
    ExpressionParam b{0.0};
    ExpressionParam c{0.0};
    ExpressionParam d{0.0};
    ExpressionParam e{0.0};
    ExpressionParam f{0.0};

    // y = a + bx + cx^2 + dx^3 + ex^4 + fx^5
    double Evaluate(double x) const {
      return (((((f.Get()*x + e.Get())*x + d.Get())*x + c.Get())*x + b.Get())*x + a.Get());
    }

    void EvaluateExpressions(const EvalContext& ctx) {
      a.Evaluate(ctx); b.Evaluate(ctx); c.Evaluate(ctx);
      d.Evaluate(ctx); e.Evaluate(ctx); f.Evaluate(ctx);
    }

  private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::make_nvp("a", a),
				cereal::make_nvp("b", b),
				cereal::make_nvp("c", c),
				cereal::make_nvp("d", d),
				cereal::make_nvp("e", e),
				cereal::make_nvp("f", f)
			);
		}
  };

  static double IntegrateVariableInTime(const std::vector<double>& y, const std::vector<double>& t) {
    const size_t n = std::min(y.size(), t.size());
    if (n < 2) return 0.0;

    // Keep only finite samples.
    std::vector<std::pair<double, double>> pts;
    pts.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      if (std::isfinite(y[i]) && std::isfinite(t[i])) {
        pts.emplace_back(t[i], y[i]); // (time, value)
      }
    }
    if (pts.size() < 2) return 0.0;

    // Sort by time and coalesce duplicate timestamps (average their values).
    std::sort(pts.begin(), pts.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    std::vector<std::pair<double, double>> uniq;
    uniq.reserve(pts.size());
    for (const auto& p : pts) {
      if (uniq.empty() || p.first != uniq.back().first) {
        uniq.push_back(p);
      } else {
        uniq.back().second = 0.5 * (uniq.back().second + p.second);
      }
    }
    if (uniq.size() < 2) return 0.0;

    // Trapezoidal rule on strictly increasing times; ignore non-positive dt.
    double area = 0.0;
    for (size_t i = 1; i < uniq.size(); ++i) {
      const double dt = uniq[i].first - uniq[i - 1].first;
      if (dt <= 0.0) continue;
      area += 0.5 * (uniq[i - 1].second + uniq[i].second) * dt;
    }
    return area; // Units: (units of y) * (units of t)
  }

} // end adsim namespace
