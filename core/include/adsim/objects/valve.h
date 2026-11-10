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

		// Opening fraction at `time` [s] from the start of the step: it
		// rises from opening_start_ to 1 over opening_ramp_time_ along the
		// C1 cubic s + (1 - s)(3x^2 - 2x^3), x = time / ramp, with a zero
		// slope at both ends of the ramp (a kink in the valve flow would
		// stall IDA's Newton iteration, as for the check valve). A ramp
		// time <= 0 or a start >= 1 gives 1: the valve is fully open for the
		// whole step. The flow of the valve law is cv_ times this fraction.
		double OpeningFraction(double time) const;
		static double OpeningFraction(double time, double start, double ramp_time);

		ExpressionParam cv_{100};	// [USGPM]
		ExpressionParam xt_{0.53};	// [-]
		ExpressionParam opening_start_{1.0};		// [-] opening at the start of the step, in (0, 1]
		ExpressionParam opening_ramp_time_{0.0};	// [s] time to open fully; 0: open at once

		void EvaluateExpressions(const EvalContext& ctx) {
			cv_.Evaluate(ctx); xt_.Evaluate(ctx);
			opening_start_.Evaluate(ctx); opening_ramp_time_.Evaluate(ctx);
		}

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::make_nvp("CV", cv_),
				cereal::make_nvp("Xt", xt_)
			);
			// Version 2: the opening ramp. Older decks keep the defaults (fully open).
			if (version >= 2) {
				ar(
					cereal::make_nvp("Opening_Start", opening_start_),
					cereal::make_nvp("Opening_Ramp_Time", opening_ramp_time_)
				);
			}
		}
	};

} // end adsim namespace

// Class versions live next to the inline serialize templates: every
// translation unit that saves or loads a deck must see the same numbers.
CEREAL_CLASS_VERSION(adsim::Valve, 2);
