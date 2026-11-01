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
		void EvaluateExpressions(const EvalContext& ctx) { if (profile_) profile_->EvaluateExpressions(ctx); }

		// Time profile of the flow, in the unit of the concrete kind (velocity,
		// volume, molar or mass flow). Constant unless the user sets a ramp.
		std::unique_ptr<IBoundaryValue> profile_;

	protected:
		IFlowBoundaryValue() = default;
		explicit IFlowBoundaryValue(double default_value)
			: profile_(std::make_unique<ConstantBoundaryValue>(default_value)) {}

		// Archive layout shared by every kind: the profile from version 2 on; a
		// version-1 archive holds a plain number under the kind's legacy name,
		// which loads as a constant profile. `self` is the concrete kind, so the
		// base-class relation cereal records is the concrete one.
		template <class Archive, class Derived>
		static void SerializeProfile(Archive& ar, Derived* self, std::uint32_t const version, const char* legacy_name) {
			ar(cereal::virtual_base_class<IFlowBoundaryValue>(self));
			if (version >= 2) {
				ar(cereal::make_nvp("Profile", self->profile_));
			} else {
				double value = 0.0;
				ar(cereal::make_nvp(legacy_name, value));
				self->profile_ = std::make_unique<ConstantBoundaryValue>(value);
			}
		}

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar) {}
	};

	class VelocityBoundaryValue : public IFlowBoundaryValue {
	public:
		VelocityBoundaryValue() : IFlowBoundaryValue(0.1) {}
		FlowType GetType() override { return FlowType::VELOCITY; }
		double GetVelocity(const Equations& state, double time, BoundaryLocation loc) override;
		double GetVolumeFlow(const Equations& state, double time, BoundaryLocation loc) override;
		double GetMolarFlow(const Equations& state, double time, BoundaryLocation loc) override;

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) { SerializeProfile(ar, this, version, "Velocity"); }
	};

	class VolumetricFlowBoundaryValue : public IFlowBoundaryValue {
	public:
		VolumetricFlowBoundaryValue() : IFlowBoundaryValue(1e-5) {}
		FlowType GetType() override { return FlowType::VOLUME; }
		double GetVelocity(const Equations& state, double time, BoundaryLocation loc) override;
		double GetVolumeFlow(const Equations& state, double time, BoundaryLocation loc) override;
		double GetMolarFlow(const Equations& state, double time, BoundaryLocation loc) override;

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) { SerializeProfile(ar, this, version, "Volume_Flow"); }
	};

	class MolarFlowBoundaryValue : public IFlowBoundaryValue {
	public:
		MolarFlowBoundaryValue() : IFlowBoundaryValue(1.0) {}
		FlowType GetType() override { return FlowType::MOLAR; }
		double GetVelocity(const Equations& state, double time, BoundaryLocation loc) override;
		double GetVolumeFlow(const Equations& state, double time, BoundaryLocation loc) override;
		double GetMolarFlow(const Equations& state, double time, BoundaryLocation loc) override;

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) { SerializeProfile(ar, this, version, "Molar_Flow"); }
	};

	class MassFlowBoundaryValue : public IFlowBoundaryValue {
	public:
		MassFlowBoundaryValue() : IFlowBoundaryValue(1e-3) {}
		FlowType GetType() override { return FlowType::MASS; }
		double GetVelocity(const Equations& state, double time, BoundaryLocation loc) override;
		double GetVolumeFlow(const Equations& state, double time, BoundaryLocation loc) override;
		double GetMolarFlow(const Equations& state, double time, BoundaryLocation loc) override;

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) { SerializeProfile(ar, this, version, "Mass_Flow"); }
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
