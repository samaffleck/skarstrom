/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <memory>

// adsim includes
#include "adsim/boundary_conditions/boundary_value.h"


namespace adsim {

  class Equations;

	static const char* BoundaryLocationList[] = {
		"Top",		// 0
		"Bottom"	// 1
	};

	enum class BoundaryLocation {
		TOP = 0,
		BOTTOM = 1
	};

	static const char* FlowTypeList[] = {
		"Velocity",				// 0
		"Volumetric Flow Rate",	// 1
		"Molar Flow Rate",		// 2
		"Mass Flow Rate"		// 3
	};

	enum class FlowType {
		VELOCITY = 0,
		VOLUME = 1,
		MOLAR = 2,
		MASS = 3
	};

	class IFlowBoundaryValue {
	public:
		virtual ~IFlowBoundaryValue() = default;
		virtual FlowType GetType() = 0;
		virtual double GetVelocity(const Equations& state, double time, BoundaryLocation loc) = 0;
		virtual double GetVolumeFlow(const Equations& state, double time, BoundaryLocation loc) = 0;
		virtual double GetMolarFlow(const Equations& state, double time, BoundaryLocation loc) = 0;
		virtual void Initialise(double step_time) {}
		virtual void EvaluateExpressions(const EvalContext& ctx) {}

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar) {}
	};

	class VelocityBoundaryValue : public IFlowBoundaryValue {
	public:
		VelocityBoundaryValue() { profile_ = std::make_unique<ConstantBoundaryValue>(0.1); }
		FlowType GetType() override { return FlowType::VELOCITY; }
		double GetVelocity(const Equations& state, double time, BoundaryLocation loc) override;
		double GetVolumeFlow(const Equations& state, double time, BoundaryLocation loc) override;
		double GetMolarFlow(const Equations& state, double time, BoundaryLocation loc) override;
		void Initialise(double step_time) override {}
		void EvaluateExpressions(const EvalContext& ctx) override { if (profile_) profile_->EvaluateExpressions(ctx); }

		std::unique_ptr<IBoundaryValue> profile_;

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(cereal::virtual_base_class<IFlowBoundaryValue>(this));
			if (version >= 2) {
				ar(cereal::make_nvp("Profile", profile_));
			} else {
				double velocity = 0.1;
				ar(cereal::make_nvp("Velocity", velocity));
				profile_ = std::make_unique<ConstantBoundaryValue>(velocity);
			}
		}
	};

	class VolumetricFlowBoundaryValue : public IFlowBoundaryValue {
	public:
		VolumetricFlowBoundaryValue() { profile_ = std::make_unique<ConstantBoundaryValue>(1e-5); }
		FlowType GetType() override { return FlowType::VOLUME; }
    double GetVelocity(const Equations& state, double time, BoundaryLocation loc) override;
		double GetVolumeFlow(const Equations& state, double time, BoundaryLocation loc) override;
		double GetMolarFlow(const Equations& state, double time, BoundaryLocation loc) override;
		void Initialise(double step_time) override {}
		void EvaluateExpressions(const EvalContext& ctx) override { if (profile_) profile_->EvaluateExpressions(ctx); }

		std::unique_ptr<IBoundaryValue> profile_;

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(cereal::virtual_base_class<IFlowBoundaryValue>(this));
			if (version >= 2) {
				ar(cereal::make_nvp("Profile", profile_));
			} else {
				double volume_flow = 1e-5;
				ar(cereal::make_nvp("Volume_Flow", volume_flow));
				profile_ = std::make_unique<ConstantBoundaryValue>(volume_flow);
			}
		}
	};

	class MolarFlowBoundaryValue : public IFlowBoundaryValue {
	public:
		MolarFlowBoundaryValue() { profile_ = std::make_unique<ConstantBoundaryValue>(1.0); }
		FlowType GetType() override { return FlowType::MOLAR; }
    double GetVelocity(const Equations& state, double time, BoundaryLocation loc) override;
		double GetVolumeFlow(const Equations& state, double time, BoundaryLocation loc) override;
		double GetMolarFlow(const Equations& state, double time, BoundaryLocation loc) override;
		void Initialise(double step_time) override {}
		void EvaluateExpressions(const EvalContext& ctx) override { if (profile_) profile_->EvaluateExpressions(ctx); }

		std::unique_ptr<IBoundaryValue> profile_;

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(cereal::virtual_base_class<IFlowBoundaryValue>(this));
			if (version >= 2) {
				ar(cereal::make_nvp("Profile", profile_));
			} else {
				double molar_flow = 1.0;
				ar(cereal::make_nvp("Molar_Flow", molar_flow));
				profile_ = std::make_unique<ConstantBoundaryValue>(molar_flow);
			}
		}
	};

	class MassFlowBoundaryValue : public IFlowBoundaryValue {
	public:
		MassFlowBoundaryValue() { profile_ = std::make_unique<ConstantBoundaryValue>(1e-3); }
		FlowType GetType() override { return FlowType::MASS; }
    double GetVelocity(const Equations& state, double time, BoundaryLocation loc) override;
		double GetVolumeFlow(const Equations& state, double time, BoundaryLocation loc) override;
		double GetMolarFlow(const Equations& state, double time, BoundaryLocation loc) override;
		void Initialise(double step_time) override {}
		void EvaluateExpressions(const EvalContext& ctx) override { if (profile_) profile_->EvaluateExpressions(ctx); }

		std::unique_ptr<IBoundaryValue> profile_;

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(cereal::virtual_base_class<IFlowBoundaryValue>(this));
			if (version >= 2) {
				ar(cereal::make_nvp("Profile", profile_));
			} else {
				double mass_flow = 1e-3;
				ar(cereal::make_nvp("Mass_Flow", mass_flow));
				profile_ = std::make_unique<ConstantBoundaryValue>(mass_flow);
			}
		}
	};

	/* Factory function */
	inline std::unique_ptr<IFlowBoundaryValue> CreateFlowBoundaryValue(FlowType type) {
		switch (type) {
		case FlowType::VELOCITY: return std::make_unique<VelocityBoundaryValue>();
		case FlowType::VOLUME:   return std::make_unique<VolumetricFlowBoundaryValue>();
		case FlowType::MOLAR:    return std::make_unique<MolarFlowBoundaryValue>();
		case FlowType::MASS:     return std::make_unique<MassFlowBoundaryValue>();
		default: return nullptr;
		}
	}

} // end adsim namespace
