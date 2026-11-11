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

	// What one endpoint (a flow stream or a pressure boundary) carries at one
	// instant, evaluated from the bed state: the values the stream log records
	// (StreamData::Log) and the KPI quadratures integrate. molar_flow and the
	// component flows are signed along the bed axis (positive upwards);
	// temperature is the stream-side temperature, power the pump or
	// compressor power [W] (0 without one).
	struct EndpointSample {
		double molar_flow = 0.0;         // [mol/s]
		double enthalpy_flow = 0.0;      // [mol K/s] molar flow times temperature
		double temperature = 0.0;        // [K]
		double pressure = 0.0;           // [Pa] bed-end face pressure
		double power = 0.0;              // [W]
		double end_cell_pressure = 0.0;  // [Pa]
		std::unordered_map<std::string, double> component_molar_flow{};  // [mol/s]
		std::unordered_map<std::string, double> mole_fraction{};         // [-]
	};

	// Mixing at an outlet flow boundary with several streams
	// (FlowBoundaryCondition::Mixing). When some streams leave the bed and
	// others flow back in through the same end (a vent that back-feeds a
	// light-product draw), the gas the leaving streams carry is the mix of
	// the bed's net outflow (bed-end gas) and the inflowing streams' own gas,
	// so every component is conserved at the node. Without it, every stream
	// carried the bed-end gas and the back-fed gas became a copy of it.
	struct NodeMixing {
		bool active = false;  // false: every stream keeps the bed-end gas
		// Composition of the gas the leaving streams carry: on the face basis
		// of the component flows (BedEndFaceMoleFraction) and on the bed-end
		// basis of the logged composition (ghost and end cell average).
		std::unordered_map<std::string, double> y_face{};
		std::unordered_map<std::string, double> y_logged{};
	};

	// Time integrals of one endpoint over one step, integrated by IDAS as
	// quadrature states alongside the solution (SolverData's KPI quadratures)
	// instead of the trapezoid over the logged samples. Flows use the KPI
	// convention: "in" is gas entering the bed through the endpoint, "out" gas
	// leaving it, both >= 0 (the net flow into the bed is in - out). `valid`
	// is false when the quadratures were not integrated for the step (turned
	// off, or the step failed).
	struct EndpointIntegrals {
		bool valid = false;
		double moles_in = 0.0;                    // [mol]
		double moles_out = 0.0;                   // [mol]
		std::unordered_map<std::string, double> component_in{};   // [mol]
		std::unordered_map<std::string, double> component_out{};  // [mol]
		double energy = 0.0;                      // [J] integral of the pump / compressor power
		double abs_axial_moles = 0.0;             // [mol] integral of |molar flow|
		double temperature_weighted_moles = 0.0;  // [K mol] integral of T |molar flow|

		// Net moles entering the bed (negative when gas left).
		double NetMolesIn() const { return moles_in - moles_out; }
		double NetComponentMolesIn(const std::string& component) const {
			const auto in = component_in.find(component);
			const auto out = component_out.find(component);
			return (in == component_in.end() ? 0.0 : in->second) - (out == component_out.end() ? 0.0 : out->second);
		}
	};

	// The same summary from the endpoint's KPI quadratures (StreamData::
	// integrals_, integrated by IDAS alongside the solution) instead of the
	// trapezoid over its logged samples; `location` is the bed end the endpoint
	// sits on. Not valid when the integrals are not. Buffer-tank consumers still
	// read ComputeSummary; this is its quadrature counterpart for comparison.
	CycleEndpointSummary ComputeQuadratureSummary(const EndpointIntegrals& integrals,
	                                              BoundaryLocation location, int cycle_index);

	class StreamData {
	public:
		void Clear();
		// Appends one sample (time from the step start).
		void Log(const EndpointSample& sample, double time);
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
    // Molar flow out of the bed the stream was logged on: positive when gas
    // left that bed, negative when it flowed back in (e.g. a valve outlet
    // whose bed fell below the downstream pressure). GetMolarFlow is signed
    // along the bed axis instead (positive upwards), so this flips it for a
    // stream at the bottom.
    double GetOutflow(double time) const;
		double GetEnthalpyFlow(double time) const;
		double GetTemperature(double time) const;
		double GetPressure(double time) const;
		// Pressure of the bed's end cell (the last packed cell at the stream's
		// bed end) at `time`; what a connected valve reads with
		// IDASettings::valve_bed_pressure_ = EndCell. Falls back to GetPressure
		// when the end cell was not logged.
		double GetEndCellPressure(double time) const;
		void LogEndCellPressure(double pressure) { end_cell_pressure_.emplace_back(pressure); }
		double GetComponentFlow(double time, const std::string& component) const;
		double GetMoleFraction(double time, const std::string& component) const;
    double GetAveragePurity(const std::string& component) const;
		
		std::vector<double> molar_flow_{};   // signed along the bed axis: positive upwards
		std::vector<double> enthalpy_flow_{};
		std::vector<double> temperature_boundary_{};
		std::vector<double> pressure_boundary_{};
		std::vector<double> end_cell_pressure_{};  // one entry per logged time (see LogEndCellPressure)
		std::vector<double> electrical_power_{};
		std::vector<double> time_series_{};
		std::unordered_map<std::string, std::vector<double>> component_molar_flow_{};
		std::unordered_map<std::string, std::vector<double>> mole_fraction_boundary_{};
		// The bed end the stream was logged at (set by IBoundaryStream::Log and PressureBoundaryCondition::Log).
		BoundaryLocation location_ = BoundaryLocation::TOP;
		// The step's quadrature integrals of this endpoint, set when the step
		// ends (not serialised; cleared with the log).
		EndpointIntegrals integrals_{};
		// True for the log of a valve stream (set by its constructor, not
		// serialised, kept by Clear). A connected valve sets its flow from its
		// own bed's pressure and its partner's logged pressure, of either sign
		// (Valve::GetMolarFlowRate). A set-flow stream connected to one must
		// replay that flow as it is, both ways: a one-sided closure on the
		// replay lets the valve draw back gas that the set-flow side never
		// gives, and the connection makes gas. The replay has no closure now
		// (BoundaryStream::GetMolarFlow), so nothing reads the flag; it marks
		// the case for any closure added there later.
		bool from_valve_ = false;
	};

	// Mole fraction of `component` on the bed-end face at `loc`: the face
	// value the species balance convects through that face (the upwind or
	// high-order reconstruction, Equations::mole_fraction_face), so an
	// outlet's component flow n * y counts the moles the balance removes.
	// Falls back to the bed-end gas (ghost and end cell average) when the
	// face values are not set up.
	double BedEndFaceMoleFraction(const Equations& state, BoundaryLocation loc, const std::string& component);

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
    // BedInletTemperature. At an outlet the component flow carries the gas
    // of the bed-end face as the species balance convects it
    // (BedEndFaceMoleFraction), so the KPIs count the moles the balance
    // removes; the enthalpy flow carries the bed-end gas temperature (the
    // average of the ghost and first interior cell). The same for every
    // stream kind.
    double GetComponentMolarFlow(const Equations& state, double time, BoundaryLocation loc, const std::string& component);
    double GetEnthalpyFlow(const Equations& state, double time, BoundaryLocation loc);
		virtual void EvaluateExpressions(const EvalContext& ctx) {
			temperature_.Evaluate(ctx);
			if (temperature_profile_) temperature_profile_->EvaluateExpressions(ctx);
			EvaluateExpressionMap(mole_fraction, ctx);
			buffer_tank_.EvaluateExpressions(ctx);
		}

    void Log(const Equations& state, double time, BoundaryLocation location);
    // The values Log records at `time` (seconds since the step start), without
    // recording them.
    // `mix`: the outlet node's mixing (FlowBoundaryCondition::Mixing); null
    // or inactive keeps the bed-end gas for every stream.
    void Log(const Equations& state, double time, BoundaryLocation location, const NodeMixing* mix);
    EndpointSample Sample(const Equations& state, double time, BoundaryLocation location,
                          const NodeMixing* mix = nullptr);
    // The gas this stream brings in when it flows back into the bed through
    // an outlet: its own composition (the connected partner's, a buffer
    // tank's or its mole fractions; see InletMoleFraction), normalised. False
    // when that composition is all zero ("no backflow composition"): the
    // back-fed gas is then the bed-end gas, as without mixing.
    bool BackflowComposition(double time, const std::vector<std::string>& components,
                             std::unordered_map<std::string, double>& y) const;
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
		// Check valve: when true, an unconnected valve passes gas only in its
		// own direction (out of the bed through an outlet, into the bed
		// through an inlet). The closure is smoothed over
		// Valve::kCheckValveBand; while it is shut the stream carries no flow,
		// so the end acts as a closed end for it. Ignored on a connected
		// stream: flow both ways through a connection (pressure equalisation,
		// light-product pressurisation) is expected. Off by default, which
		// keeps the plain valve law (flow of either sign).
		bool check_valve_ = false;

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IBoundaryStream>(this),
				cereal::make_nvp("Valve", valve_),
				cereal::make_nvp("Boundary_Value", boundary_value_)
			);
			// v2: the check valve (off for older decks).
			if (version >= 2) ar(cereal::make_nvp("Check_Valve", check_valve_));
			else check_valve_ = false;
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
			ultimate_pressure_.Evaluate(ctx);
		}

    std::unique_ptr<IPumpPowerModel> power_{};
    std::unique_ptr<IPump> pump_{};
    // Ultimate (blank-off) pressure of a vacuum pump [Pa]: the suction
    // pressure at which the pump stops delivering gas. 0 turns it off (the
    // flow follows the power alone). Used on an OUTLET (vacuum pump) only;
    // the flow is scaled by UltimatePressureFactor.
    ExpressionParam ultimate_pressure_{0.0};

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
			// Version 2 added the ultimate pressure; older decks have none (off).
			if (version >= 2) ar(cereal::make_nvp("Ultimate_Pressure", ultimate_pressure_));
		}
	};

	// Factor on a vacuum pump's flow for its ultimate pressure p_ult at the
	// suction pressure P1: the rotary-pump speed law u = 1 - p_ult / P1
	// (S_eff = S (1 - p_ult / P1)), with the corner at u = 0 rounded so the
	// factor is C1 (IDA's Newton iteration needs a smooth residual):
	//   f = 0                               for u <= 0 (P1 <= p_ult)
	//   f = 2 u^2 / d - u^3 / d^2           for 0 < u < d
	//   f = u                               for u >= d
	// with d = kUltimatePressureBlend. f is non-decreasing in P1, so the bed
	// pressure approaches p_ult from above and is held there. 1 when p_ult <= 0
	// (no ultimate pressure).
	inline constexpr double kUltimatePressureBlend = 0.02;
	inline double UltimatePressureFactor(double P1, double p_ult) {
		if (!(p_ult > 0.0)) return 1.0;
		if (!(P1 > p_ult)) return 0.0;
		const double u = 1.0 - p_ult / P1;
		const double d = kUltimatePressureBlend;
		if (u >= d) return u;
		return u * u * (2.0 / d - u / (d * d));
	}

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

// Class versions live next to the inline serialize templates: every
// translation unit that saves or loads a deck must see the same numbers.
CEREAL_CLASS_VERSION(adsim::IBoundaryStream, 5);
CEREAL_CLASS_VERSION(adsim::BoundaryStream, 1);
CEREAL_CLASS_VERSION(adsim::ValveBoundaryStream, 2);
CEREAL_CLASS_VERSION(adsim::PumpBoundaryStream, 2);
