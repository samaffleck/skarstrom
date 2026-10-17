/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <memory>

// Cereal includes
#include <cereal/types/memory.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>

#include "adsim/math/expression_param.h"


namespace adsim {

  class Equations;

	static const char* WallTemperatureBoundaryConditionList[] = {
		"Zero flux",
		"Fixed"
	};

	enum class WallTemperatureBoundaryConditionType {
		ZERO_FLUX,
		FIXED
	};

	class IWallTemperatureBoundaryCondition {
	public:
		virtual ~IWallTemperatureBoundaryCondition() = default;
		virtual void Update(Equations& state) = 0;
		virtual WallTemperatureBoundaryConditionType GetModelType() = 0;
		virtual void EvaluateExpressions(const EvalContext& ctx) {}

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {}
	};

  class WallTemperatureBoundaryConditionFixed : public IWallTemperatureBoundaryCondition {
	public:
    void Update(Equations& state) override;
		WallTemperatureBoundaryConditionType GetModelType() override { return WallTemperatureBoundaryConditionType::FIXED; }
		void EvaluateExpressions(const EvalContext& ctx) override { temperature_.Evaluate(ctx); }

		ExpressionParam temperature_{293.0};

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IWallTemperatureBoundaryCondition>(this),
				cereal::make_nvp("Fixed_Wall_Temperature", temperature_)
			);
		}
	};

	class WallTemperatureBoundaryConditionZeroFlux : public IWallTemperatureBoundaryCondition {
	public:
		void Update(Equations& state) override;
		WallTemperatureBoundaryConditionType GetModelType() override { return WallTemperatureBoundaryConditionType::ZERO_FLUX; }

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(cereal::virtual_base_class<IWallTemperatureBoundaryCondition>(this));
		}
	};

	// Factory Functions
	inline std::unique_ptr<IWallTemperatureBoundaryCondition> CreateWallTemperatureBoundaryCondition(WallTemperatureBoundaryConditionType type) {
		switch (type) {
		case WallTemperatureBoundaryConditionType::FIXED: return std::make_unique<WallTemperatureBoundaryConditionFixed>();
		case WallTemperatureBoundaryConditionType::ZERO_FLUX: return std::make_unique<WallTemperatureBoundaryConditionZeroFlux>();
		default: return nullptr;
		}
	}

} // end adsim namespace
