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
