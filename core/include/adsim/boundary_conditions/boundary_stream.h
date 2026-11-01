/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <algorithm>
#include <string>
#include <unordered_map>
#include <memory>
#include <vector>

// adsim includes
#include "adsim/boundary_conditions/flow_boundary_value.h"
#include "adsim/boundary_conditions/boundary_value.h"
#include "adsim/boundary_conditions/outlet_conditions.h"
#include "adsim/boundary_conditions/connection.h"
#include "adsim/boundary_conditions/cycle_endpoint_summary.h"
#include "adsim/boundary_conditions/buffer_tank_source.h"
#include "adsim/math/expression_map.h"
#include "adsim/objects/fluid.h"
#include "adsim/objects/valve.h"
#include "adsim/objects/compressor.h"


namespace adsim {

  class Equations;

	static const char* BoundaryStreamTypeList[] = {
		"None",	// 0
		"Valve",	// 1
		"Pump"	// 2
	};

  // Component Type
	enum class BoundaryStreamType {
		NONE = 0,
		VALVE = 1,
		PUMP = 2
	};

	static const char* BoundaryDirectionList[] = {
		"Inlet",	// 0 
		"Outlet"	// 1
	};

	enum class BoundaryDirection {
		INLET = 0,
		OUTLET = 1
	};

	// Make an inlet mole-fraction map keyed exactly by the fluid's components:
	// every component gets an entry (0.0 when missing, e.g. a component added to
	// the fluid after the map was last edited) and entries for components no
	// longer in the fluid are erased. Run by the boundary checks so the residual
	// never meets a missing key.
	inline void ReconcileMoleFractionMap(ExpressionMap& mole_fraction,
	                                     const std::vector<std::string>& components) {
		for (auto it = mole_fraction.begin(); it != mole_fraction.end();) {
			if (std::find(components.begin(), components.end(), it->first) == components.end())
				it = mole_fraction.erase(it);
			else
				++it;
		}
		for (const auto& comp : components)
			mole_fraction.try_emplace(comp, 0.0);
	}

	// Pump type a boundary needs for its flow direction: an INLET pushes gas
	// into the bed (compressor), an OUTLET draws gas out of it (vacuum pump).
	inline PumpType PumpTypeForDirection(BoundaryDirection direction) {
		return direction == BoundaryDirection::INLET ? PumpType::COMPRESSOR : PumpType::VACUUM_PUMP;
	}

	// Replaces `pump` with a default pump of the type required by `direction`,
	// but only when it is null or of the wrong type. A pump that already matches
	// the direction is left untouched so its user settings (use_pump_,
	// efficiency model, polytropic index, energy model) survive. Returns true if
	// the pump was replaced.
	inline bool EnsurePumpMatchesDirection(std::unique_ptr<IPump>& pump, BoundaryDirection direction) {
		const PumpType required = PumpTypeForDirection(direction);
		if (pump && pump->GetType() == required) return false;
		pump = CreatePump(required);
		return true;
	}

	class StreamData {
	public:
		void Clear();
		void Log(
      double molar_flow, 
      double enthalpy_flow, 
      double temperature_boundary,
      double pressure_boundary,
      double power,
      double time,
			std::unordered_map<std::string, double> component_molar_flow,
			std::unordered_map<std::string, double> mole_fraction_boundary
    );
    double GetMolarFlow(double time) const;
		double GetEnthalpyFlow(double time) const;
		double GetTemperature(double time) const;
		double GetPressure(double time) const;
		double GetComponentFlow(double time, const std::string& component) const;
		double GetMoleFraction(double time, const std::string& component) const;
    double GetAveragePurity(const std::string& component) const;
		
		std::vector<double> molar_flow_{};
		std::vector<double> enthalpy_flow_{};
		std::vector<double> temperature_boundary_{};
		std::vector<double> pressure_boundary_{};
		std::vector<double> electrical_power_{};
		std::vector<double> time_series_{};
		std::unordered_map<std::string, std::vector<double>> component_molar_flow_{};
		std::unordered_map<std::string, std::vector<double>> mole_fraction_boundary_{};
	};

	class IBoundaryStream {
	public:
		virtual ~IBoundaryStream() = default;
		virtual BoundaryStreamType GetType() = 0;
		virtual void Initialise(Equations& state, double step_time, BoundaryLocation location, BoundaryDirection dir) = 0;
		virtual bool PerformSpecificChecks(std::string& error_msg, const Fluid& fluid) = 0;
    virtual double GetPower(double outlet_pressure, double molar_flow, double time) = 0;

    virtual double GetMolarFlow(const Equations& state,double time, BoundaryLocation loc) = 0;
    // Molar flow of one component and the enthalpy flow (molar flow times the
    // gas temperature) through this stream, with the sign of GetMolarFlow. An
    // inlet carries the stream-side composition and the temperature from
    // BedInletTemperature; an outlet carries the bed-end gas (the average of
    // the ghost and first interior cell). The same for every stream kind.
    double GetComponentMolarFlow(const Equations& state, double time, BoundaryLocation loc, const std::string& component);
    double GetEnthalpyFlow(const Equations& state, double time, BoundaryLocation loc);
		virtual void EvaluateExpressions(const EvalContext& ctx) {
			temperature_.Evaluate(ctx);
			if (temperature_profile_) temperature_profile_->EvaluateExpressions(ctx);
			EvaluateExpressionMap(mole_fraction, ctx);
			buffer_tank_.EvaluateExpressions(ctx);
		}

    void Log(const Equations& state, double time, BoundaryLocation location);
    bool PerformChecks(std::string& error_msg, const Fluid& fluid);
		void SetDirection(BoundaryDirection dir) { direction_ = dir; }

    // The user-specified inlet temperature at `time` (seconds since step
    // start): temperature_ when the profile is constant, otherwise the ramp
    // (see EvaluateInletTemperatureProfile). In an isothermal run it is the
    // isothermal temperature instead. This is the first-cycle / unconnected
    // value; InletTemperature applies the buffer-tank / connection precedence.
    double AssumedInletTemperature(double time) const {
      if (isothermal_temperature_) return *isothermal_temperature_;
      return EvaluateInletTemperatureProfile(temperature_, temperature_profile_.get(), time);
    }

    // Stream-side inlet temperature / mole fraction (buffer tank, connected
    // stream or user value; see SelectInletTemperature). Used by every stream
    // kind's flow functions and by Log so they always agree. Isothermal runs
    // bypass the selection: every stream-side gas is at the isothermal
    // temperature.
    double InletTemperature(double time) const {
      if (isothermal_temperature_) return *isothermal_temperature_;
      return SelectInletTemperature(buffer_tank_, connection_, AssumedInletTemperature(time), time);
    }
    double InletMoleFraction(double time, const std::string& component) const {
      return SelectInletMoleFraction(buffer_tank_, connection_, mole_fraction, component, time);
    }

		StreamData stream_data_{};
		ExpressionMap mole_fraction{}; // Inlet mole fraction (also acts as first-cycle assumption when buffer_tank_ is active)
		std::string name_{};
		// Stable identity of this endpoint within the cycle (unique across every
		// stream and pressure boundary of every step; 0 = not yet assigned, see
		// Cycle::EnsureIds). Survives renames; KPI endpoint references, buffer-tank
		// sources and parameter paths ("endpoint[<id>]") refer to it.
		size_t endpoint_id_ = 0;
		ExpressionParam temperature_{293.0}; // Inlet temperature (also first-cycle assumption when buffer_tank_ is active)
    // Optional inlet-temperature time profile (linear / exponential ramp over
    // the step, reusing the boundary-value profile classes). nullptr — the
    // default, and what every pre-v3 archive loads as — means constant
    // temperature_. A ConstantBoundaryValue is treated the same way.
    std::unique_ptr<IBoundaryValue> temperature_profile_{};
    // Isothermal-run override (transient, never serialised): set from the
    // solver state by the parent boundary condition's Initialise.
    std::optional<double> isothermal_temperature_{};
    Connection connection_;
    // Aggregates of this endpoint's most-recently-completed run, used by any
    // consumer that has buffer_tank_.summary_ pointing here.
    CycleEndpointSummary cycle_summary_{};
    // When active, the inlet pulls its temperature and composition from
    // another step's cycle_summary_ rather than from temperature_ and
    // mole_fraction, for every stream kind. A direct-flow BoundaryStream also
    // takes its molar flow from the summary (instead of boundary_flow_); a
    // valve or pump keeps its pressure/power-driven flow, evaluated at the
    // tank's temperature. Mutually exclusive with connection_.
    BufferTankSource buffer_tank_{};

    // Last step duration seen by Initialise. Buffer-tank inlets use this to
    // turn `total_moles` from the source into a constant molar flow rate.
    double step_time_cached_ = 0.0;

  protected:
    // Temperature of the gas an INLET delivers to the bed at `time`: the
    // stream-side InletTemperature, raised by compression heating in the
    // pump-driven stream kinds (they override this).
    virtual double BedInletTemperature(const Equations& state, double time, BoundaryLocation loc) {
      (void)state; (void)loc;
      return InletTemperature(time);
    }

    // Set during initialisation to be equal to parent boundary condition direction
		BoundaryDirection direction_ = BoundaryDirection::OUTLET;

  private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::make_nvp("Stream_Name", name_),
				cereal::make_nvp("Boundary_Temperature", temperature_)
			);
			// version < 5: the mole fractions were plain numbers.
			SerializeExpressionMap(ar, "Boundary_Mole_Fraction", mole_fraction, version < 5);
			ar(cereal::make_nvp("Boundary_Connection", connection_));

			if (version >= 2) {
				ar(
					cereal::make_nvp("Buffer_Tank_Source", buffer_tank_),
					cereal::make_nvp("Cycle_Endpoint_Summary", cycle_summary_)
				);
			}
			if (version >= 3) {
				ar(cereal::make_nvp("Temperature_Profile", temperature_profile_));
			}
			// version < 3: temperature_profile_ stays nullptr (constant temperature_).
			if (version >= 4) {
				ar(cereal::make_nvp("Endpoint_ID", endpoint_id_));
			}
			// version < 4: endpoint_id_ stays 0 and Cycle::EnsureIds assigns it on load.
		}
	};

  // ------------------------------------------------------
  // Boundary condition with a valve.  
  // ------------------------------------------------------
	class ValveBoundaryStream : public IBoundaryStream {
	public:
		ValveBoundaryStream();

		BoundaryStreamType GetType() override { return BoundaryStreamType::VALVE; }
		void Initialise(Equations& state, double step_time, BoundaryLocation location, BoundaryDirection dir) override;
		bool PerformSpecificChecks(std::string& error_msg, const Fluid& fluid) override;
    double GetPower(double outlet_pressure, double molar_flow, double time) override { return 0.0; } // No compressor

    double GetMolarFlow(const Equations& state,double time, BoundaryLocation loc) override;
		void EvaluateExpressions(const EvalContext& ctx) override {
			IBoundaryStream::EvaluateExpressions(ctx);
			if (boundary_value_) boundary_value_->EvaluateExpressions(ctx);
			valve_.EvaluateExpressions(ctx);
		}

    std::unique_ptr<IBoundaryValue> boundary_value_; // Pressure
		Valve valve_{};

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IBoundaryStream>(this),
				cereal::make_nvp("Valve", valve_),
				cereal::make_nvp("Boundary_Value", boundary_value_)
			);
		}
	};

	class BoundaryStream : public IBoundaryStream {
	public:
		BoundaryStream();
		BoundaryStreamType GetType() override { return BoundaryStreamType::NONE; }
		void Initialise(Equations& state, double step_time, BoundaryLocation location, BoundaryDirection dir) override {
			direction_ = dir;
			step_time_cached_ = step_time;
			if (boundary_flow_) boundary_flow_->Initialise(step_time);
		}
		bool PerformSpecificChecks(std::string& error_msg, const Fluid& fluid) override;
    double GetPower(double outlet_pressure, double molar_flow, double time) override;

    double GetMolarFlow(const Equations& state,double time, BoundaryLocation loc) override;
		void EvaluateExpressions(const EvalContext& ctx) override {
			IBoundaryStream::EvaluateExpressions(ctx);
			if (boundary_flow_) boundary_flow_->EvaluateExpressions(ctx);
			if (pump_) pump_->EvaluateExpressions(ctx);
		}

		std::unique_ptr<IFlowBoundaryValue> boundary_flow_;
    std::unique_ptr<IPump> pump_{};

	protected:
		double BedInletTemperature(const Equations& state, double time, BoundaryLocation loc) override;

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IBoundaryStream>(this),
				cereal::make_nvp("Boundary_Flow", boundary_flow_),
				cereal::make_nvp("Pump", pump_)
			);
		}
	};

	class PumpBoundaryStream : public IBoundaryStream {
	public:
		PumpBoundaryStream();
		BoundaryStreamType GetType() override { return BoundaryStreamType::PUMP; }
		void Initialise(Equations& state, double step_time, BoundaryLocation location, BoundaryDirection dir) override { direction_ = dir; step_time_cached_ = step_time; }
		bool PerformSpecificChecks(std::string& error_msg, const Fluid& fluid) override;
    double GetPower(double outlet_pressure, double molar_flow, double time) override;

    double GetMolarFlow(const Equations& state,double time, BoundaryLocation loc) override;
		void EvaluateExpressions(const EvalContext& ctx) override {
			IBoundaryStream::EvaluateExpressions(ctx);
			if (power_) power_->EvaluateExpressions(ctx);
			if (pump_) pump_->EvaluateExpressions(ctx);
		}

    std::unique_ptr<IPumpPowerModel> power_{};
    std::unique_ptr<IPump> pump_{};

	protected:
		double BedInletTemperature(const Equations& state, double time, BoundaryLocation loc) override;

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IBoundaryStream>(this),
				cereal::make_nvp("Power_Model", power_),
				cereal::make_nvp("Pump", pump_)
			);
		}
	};

	/* Factory function */
	inline std::unique_ptr<IBoundaryStream> CreateStream(BoundaryStreamType type) {
		switch (type) {
      case BoundaryStreamType::NONE: return std::make_unique<BoundaryStream>();
      case BoundaryStreamType::VALVE: return std::make_unique<ValveBoundaryStream>();
      case BoundaryStreamType::PUMP: return std::make_unique<PumpBoundaryStream>();
      default: return nullptr;
		}
	}

} // end adsim namespace
