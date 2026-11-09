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

	class Valve {
	public:
		double GetMolarFlowRate(double upstream_pressure, double downstream_pressure,
			double molar_mass, double upstream_temperature, double heat_capacity) const;

		// Pressure drop that a check valve passes, smoothed over a band of
		// kCheckValveBand pascals: 0 for dp <= 0 (reverse flow is blocked),
		// exactly dp for dp >= kCheckValveBand, and the C1 cubic
		// 2 dp^2 / band - dp^3 / band^2 in between (monotonic, below dp).
		// A hard max(dp, 0) would put a kink in the boundary residual where
		// the flow decays through zero, and IDA's Newton iteration then
		// fails to converge there; the cubic keeps the residual
		// differentiable with a continuous slope, and passes the full flow
		// away from the closure.
		static double CheckValvePressureDrop(double dp);
		static constexpr double kCheckValveBand = 10.0;	// [Pa]

		ExpressionParam cv_{100};	// [USGPM]
		ExpressionParam xt_{0.53};	// [-]

		void EvaluateExpressions(const EvalContext& ctx) {
			cv_.Evaluate(ctx); xt_.Evaluate(ctx);
		}

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::make_nvp("CV", cv_),
				cereal::make_nvp("Xt", xt_)
			);
		}
	};

} // end adsim namespace

// Class versions live next to the inline serialize templates: every
// translation unit that saves or loads a deck must see the same numbers.
CEREAL_CLASS_VERSION(adsim::Valve, 1);
