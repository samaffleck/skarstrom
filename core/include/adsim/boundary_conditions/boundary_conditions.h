/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>

// Eigen
#include <Eigen/Dense>

// Cereal includes
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>

// adsim includes
#include "adsim/boundary_conditions/boundary_value.h"
#include "adsim/boundary_conditions/flow_boundary_value.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/connection.h"
#include "adsim/objects/compressor.h"


namespace adsim {
	// Forward declerations
	struct OutletConditions;
	class State;
	class Fluid;
	class Reactor;

	static const char* BoundaryTypeList[] = {
		"Pressure",	// 0 
		"Flow",		// 1
    "Closed" // 2
	};

	enum class BoundaryType {
		PRESSURE = 0,
		FLOW = 1,
    CLOSED = 2
	};

	class IBoundaryCondition {
	public:
		virtual ~IBoundaryCondition() = default;
		virtual void Initialise(Equations& state, double step_time) = 0;
		// Non-const: inlet checks also reconcile the inlet mole-fraction map with the
		// fluid's components (see ReconcileMoleFractionMap).
		virtual bool PerformChecks(std::string& error_msg, const Fluid& fluid) = 0;
		virtual void Apply(Equations& state, double time) = 0;
    virtual void Log(const Equations& state, double time) = 0;
		virtual void ClearLog() = 0;
		virtual BoundaryType GetType() = 0;
		virtual void EvaluateExpressions(const EvalContext& ctx) {}

		BoundaryLocation location_ = BoundaryLocation::BOTTOM;
    BoundaryDirection direction_ = BoundaryDirection::OUTLET; 
	
	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::make_nvp("Boundary_Condition_Location", location_),
				cereal::make_nvp("Boundary_Condition_Direction", direction_)
			);
		}
	};

	class PressureBoundaryCondition : public IBoundaryCondition {
	public:
		PressureBoundaryCondition();

		void Initialise(Equations& state, double step_time) override;
		void Apply(Equations& state, double time) override;
    bool PerformChecks(std::string& error_msg, const Fluid& fluid) override;
		void Log(const Equations& state, double time) override;
		void ClearLog() override;
		BoundaryType GetType() override { return BoundaryType::PRESSURE; }

		// Stream-side inlet temperature / mole fraction (buffer tank, connected
		// stream or user value; see SelectInletTemperature). Shared by Apply and
		// Log so the residuals and the logged stream always use the same source.
		// Isothermal runs bypass the selection (every stream-side gas is at the
		// isothermal temperature); otherwise the user value follows the optional
		// temperature profile (see IBoundaryStream::AssumedInletTemperature).
		double AssumedInletTemperature(double time) const {
			if (isothermal_temperature_) return *isothermal_temperature_;
			return EvaluateInletTemperatureProfile(temperature_, temperature_profile_.get(), time);
		}
		double InletTemperature(double time) const {
			if (isothermal_temperature_) return *isothermal_temperature_;
			return SelectInletTemperature(buffer_tank_, connection_, AssumedInletTemperature(time), time);
		}
		double InletMoleFraction(double time, const std::string& component) const {
			return SelectInletMoleFraction(buffer_tank_, connection_, mole_fraction_, component, time);
		}

		void EvaluateExpressions(const EvalContext& ctx) override {
			temperature_.Evaluate(ctx);
			if (temperature_profile_) temperature_profile_->EvaluateExpressions(ctx);
			if (pressure_) pressure_->EvaluateExpressions(ctx);
			EvaluateExpressionMap(mole_fraction_, ctx);
			if (pump_) pump_->EvaluateExpressions(ctx);
			buffer_tank_.EvaluateExpressions(ctx);
		}

		ExpressionParam temperature_{293.0};
		// Stable endpoint identity (see IBoundaryStream::endpoint_id_).
		size_t endpoint_id_ = 0;
		// Optional inlet-temperature ramp; nullptr (default / pre-v3) = constant
		// temperature_. Same semantics as IBoundaryStream::temperature_profile_.
		std::unique_ptr<IBoundaryValue> temperature_profile_{};
		// Isothermal-run override (transient, never serialised; set in Initialise).
		std::optional<double> isothermal_temperature_{};
		ExpressionMap mole_fraction_{};
		std::unique_ptr<IBoundaryValue> pressure_;
    StreamData stream_data_{};
    Connection connection_{};
    std::unique_ptr<IPump> pump_{};
    CycleEndpointSummary cycle_summary_{};
    BufferTankSource buffer_tank_{};
    double step_time_cached_ = 0.0;

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IBoundaryCondition>(this),
				cereal::make_nvp("Inlet_Temperature", temperature_)
			);
			// version < 5: the mole fractions were plain numbers.
			SerializeExpressionMap(ar, "Inlet_Mole_Fraction", mole_fraction_, version < 5);
			ar(
				cereal::make_nvp("Boundary_Pressure", pressure_),
        cereal::make_nvp("Boundary_Condition_Connection", connection_),
        cereal::make_nvp("Pump", pump_)
			);
			if (version >= 2) {
				ar(
					cereal::make_nvp("Buffer_Tank_Source", buffer_tank_),
					cereal::make_nvp("Cycle_Endpoint_Summary", cycle_summary_)
				);
			}
			if (version >= 3) {
				ar(cereal::make_nvp("Inlet_Temperature_Profile", temperature_profile_));
			}
			// version < 3: temperature_profile_ stays nullptr (constant temperature_).
			if (version >= 4) {
				ar(cereal::make_nvp("Endpoint_ID", endpoint_id_));
			}
			// version < 4: endpoint_id_ stays 0 and Cycle::EnsureIds assigns it on load.
		}
	};

	class FlowBoundaryCondition : public IBoundaryCondition {
	public:
		void Initialise(Equations& state, double step_time) override;
		void Apply(Equations& state, double time) override;
    bool PerformChecks(std::string& error_msg, const Fluid& fluid) override;
		void Log(const Equations& state, double time) override;
		void ClearLog() override;
		BoundaryType GetType() override { return BoundaryType::FLOW; }
		void EvaluateExpressions(const EvalContext& ctx) override {
			for (auto& s : streams_) if (s) s->EvaluateExpressions(ctx);
		}

		void AddStream(const std::string& stream_name);
		void RemoveStream(const std::string& stream_name);

		std::vector<std::unique_ptr<IBoundaryStream>> streams_{};

	private:
		double GetTotalMolarFlow(const Equations& state, double time) const;
		double GetTotalComponentMolarFlow(const Equations& state, const std::string& component, double time) const;
		double GetTotalEnthalpyFlow(const Equations& state, double time) const;

		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IBoundaryCondition>(this),
				cereal::make_nvp("Boundary_Stream_List", streams_)
			);
		}
	};

	class ClosedBoundaryCondition : public IBoundaryCondition {
	public:
    ClosedBoundaryCondition() { direction_ = BoundaryDirection::OUTLET; }  
		void Initialise(Equations& state, double step_time) override { /* empty */};
		void Apply(Equations& state, double time) override;
    bool PerformChecks(std::string& error_msg, const Fluid& fluid) override;
		void Log(const Equations& state, double time) override;
		void ClearLog() override;
		BoundaryType GetType() override { return BoundaryType::CLOSED; }
		
	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IBoundaryCondition>(this)
			);
		}
	};

	// Factory Functions
	inline std::unique_ptr<IBoundaryCondition> CreateBoundaryCondition(BoundaryType type) {
		switch (type) {
		case BoundaryType::PRESSURE: return std::make_unique<PressureBoundaryCondition>();
		case BoundaryType::FLOW: return std::make_unique<FlowBoundaryCondition>();
		case BoundaryType::CLOSED: return std::make_unique<ClosedBoundaryCondition>();
		default: return nullptr;
		}
	}

} // end adsim namespace
