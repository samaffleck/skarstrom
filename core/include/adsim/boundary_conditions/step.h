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
#include <memory>
#include <map>

// Cereal includes
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>

// adsim includes
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/outlet_conditions.h"
#include "adsim/objects/reactor_equations.h"
#include "adsim/math/expression_param.h"

namespace adsim {

	// A heating or cooling medium on the wall's outer face during one step: the
	// fluid in a jacket around the column, as in temperature swing and direct
	// air capture plants. While the step runs, the wall exchanges heat with the
	// medium instead of the surroundings:
	//
	//   q_o = h_m (T_w - T_m(t))      [W/m2 of outer wall surface]
	//
	// so h_m replaces each layer's external heat transfer coefficient and T_m
	// replaces Bed::ambient_temperature_. T_m is set, not solved (no energy
	// balance on the medium): a constant or a ramp from the step start, with the
	// same profile classes and semantics as an inlet temperature (see
	// EvaluateInletTemperatureProfile). A disabled medium (the default) leaves
	// the step on the deck's ambient exchange. One medium covers the whole bed.
	//
	// The same medium can also run through an internal heat exchanger (tubes
	// or plates inside the packing) that heats or cools the bed directly:
	//
	//   q_x = a_x h_x (T_m(t) - T_bed)  [W/m3 of bed]
	//
	// with a_x the exchanger's heat-transfer area per bed volume and h_x its
	// coefficient. It acts in every bed cell (not in void regions), alongside
	// the jacket; a_x = 0 (the default) turns it off, and h_m = 0 leaves an
	// exchanger without a jacket. Under LTNE the exchanger's heat is shared
	// between the gas (eb) and the solid (1 - eb), as the wall's is.
	//
	// The class has its own cereal version, so options for the medium extend
	// it without changing Step.
	class HeatingMedium {
	public:
		HeatingMedium() = default;
		HeatingMedium(const HeatingMedium& other) { *this = other; }
		HeatingMedium& operator=(const HeatingMedium& other) {
			if (this == &other) return *this;
			enabled_ = other.enabled_;
			temperature_ = other.temperature_;
			heat_transfer_coefficient_ = other.heat_transfer_coefficient_;
			exchanger_area_per_volume_ = other.exchanger_area_per_volume_;
			exchanger_heat_transfer_coefficient_ = other.exchanger_heat_transfer_coefficient_;
			temperature_profile_ = CloneProfile(other.temperature_profile_.get());
			return *this;
		}
		HeatingMedium(HeatingMedium&&) = default;
		HeatingMedium& operator=(HeatingMedium&&) = default;

		bool enabled_ = false;
		ExpressionParam temperature_{293.0};              // [K] T_m when the profile is constant
		// Optional ramp of T_m over the step; nullptr or a ConstantBoundaryValue
		// means the constant temperature_.
		std::unique_ptr<IBoundaryValue> temperature_profile_{};
		ExpressionParam heat_transfer_coefficient_{100.0}; // [W/(m2 K)] h_m, outer wall surface
		// Internal heat exchanger (off while a_x = 0).
		ExpressionParam exchanger_area_per_volume_{0.0};             // [m2/m3] a_x, per bed volume
		ExpressionParam exchanger_heat_transfer_coefficient_{100.0}; // [W/(m2 K)] h_x

		// a_x h_x [W/(m3 K)]: the exchanger's volumetric conductance (0 when off).
		double ExchangerConductance() const {
			return exchanger_area_per_volume_.Get() * exchanger_heat_transfer_coefficient_.Get();
		}

		// T_m [K] at `time` seconds after the step start.
		double Temperature(double time) const {
			return EvaluateInletTemperatureProfile(temperature_, temperature_profile_.get(), time);
		}

		// An enabled medium needs a positive temperature (at the profile's start
		// and end), finite non-negative h_m, a_x and h_x, and a valid
		// exponential factor.
		// A disabled one is not checked. Appends to error_message on failure.
		bool PerformChecks(const std::string& step_name, std::string& error_message) const;

		void EvaluateExpressions(const EvalContext& ctx) {
			temperature_.Evaluate(ctx);
			heat_transfer_coefficient_.Evaluate(ctx);
			exchanger_area_per_volume_.Evaluate(ctx);
			exchanger_heat_transfer_coefficient_.Evaluate(ctx);
			if (temperature_profile_) temperature_profile_->EvaluateExpressions(ctx);
		}

	private:
		static std::unique_ptr<IBoundaryValue> CloneProfile(const IBoundaryValue* v) {
			if (!v) return nullptr;
			if (auto* c = dynamic_cast<const ConstantBoundaryValue*>(v)) return std::make_unique<ConstantBoundaryValue>(*c);
			if (auto* l = dynamic_cast<const LinearBoundaryValue*>(v)) return std::make_unique<LinearBoundaryValue>(*l);
			if (auto* e = dynamic_cast<const ExponentialBoundaryValue*>(v)) return std::make_unique<ExponentialBoundaryValue>(*e);
			if (auto* g = dynamic_cast<const LogLinearBoundaryValue*>(v)) return std::make_unique<LogLinearBoundaryValue>(*g);
			if (auto* d = dynamic_cast<const DoubleExponentialBoundaryValue*>(v)) return std::make_unique<DoubleExponentialBoundaryValue>(*d);
			return nullptr;
		}

		friend class cereal::access;
		template<class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::make_nvp("Enabled", enabled_),
				cereal::make_nvp("Temperature", temperature_),
				cereal::make_nvp("Temperature_Profile", temperature_profile_),
				cereal::make_nvp("Heat_Transfer_Coefficient", heat_transfer_coefficient_)
			);
			if (version >= 1) {
				ar(
					cereal::make_nvp("Exchanger_Area_Per_Volume", exchanger_area_per_volume_),
					cereal::make_nvp("Exchanger_Heat_Transfer_Coefficient", exchanger_heat_transfer_coefficient_)
				);
			}
			// version 0: no internal exchanger (a_x = 0).
		}
	};

	class Step {
	public:
		Step();

		void Initialise(Equations& state);
    void LogOutletConditions(Equations& state, double time);
    void ClearOutletConditions();

    // Walks both boundary conditions, computes a CycleEndpointSummary from each
    // endpoint's stream_data_ (the just-finished step's logged time series)
    // and stores it on the endpoint's cycle_summary_. Called at the end of
    // every IntegrateStep so any later step in the same or subsequent cycle
    // can read the most-recent run's aggregates via BufferTankSource.
    void UpdateCycleSummary(int cycle_index);

    // Invalidates every endpoint's cycle_summary_ (is_valid = false). Called
    // by both solvers at the start of a run so a summary left over from a
    // previous run (or loaded from a saved deck) is never consumed: until the
    // source step completes in THIS run, buffer-tank consumers use their
    // first-cycle assumption values.
    void ClearCycleSummary();

    // Returns the summary attached to the named endpoint, or nullptr if no
    // such endpoint exists on this step. For FlowBoundaryCondition streams,
    // the lookup key is the stream's name_; for PressureBoundaryCondition
    // it is kPressureBoundaryEndpointName.
    CycleEndpointSummary* FindSummaryForEndpoint(const std::string& endpoint_name);

    // Names of every endpoint exposed by this step's two boundary conditions.
    // Used by the GUI dropdown and by BufferTankSource resolution.
    std::vector<std::string> GetEndpointNames() const;

    // One endpoint exposed by a step: its name plus the bed end it sits on.
    struct Endpoint {
      std::string name;
      BoundaryLocation location = BoundaryLocation::BOTTOM;
      size_t id = 0;  // the stream's / pressure boundary's endpoint_id_
    };

    // Every endpoint (flow streams + pressure boundary) with its location.
    // Flow BC -> one Endpoint per stream (name = stream name); pressure BC ->
    // a single Endpoint named PressureEndpointName(step_name_, location) (the
    // BoundaryLog key); closed BC -> none.
    // Used by the KPI config to resolve an endpoint's bed end.
    std::vector<Endpoint> GetEndpoints() const;

    // KPI / BoundaryLog name of a step's pressure-boundary endpoint:
    // "<step> Pressure Boundary <Top|Bottom>". It embeds the step name, so a
    // step rename must rewrite KPI refs (KPIConfig::OnStepRenamed).
    static std::string PressureEndpointName(const std::string& step_name, BoundaryLocation location) {
      return step_name + " Pressure Boundary " + (location == BoundaryLocation::TOP ? "Top" : "Bottom");
    }

		std::string	step_name_ = "Step";
		// Stable identity within the cycle (0 = not yet assigned; see
		// Cycle::EnsureIds). The name is a display label: references to a step
		// (KPI endpoints, buffer-tank sources, parameter paths "step[<id>]")
		// use this ID, so renaming or reordering steps does not break them.
		size_t id_ = 0;
		std::shared_ptr<IBoundaryCondition> top_boundary_condition_ = nullptr;
		std::shared_ptr<IBoundaryCondition> bottom_boundary_condition_ = nullptr;
		ExpressionParam step_time_{100.0};      // [s] Absolute step duration (also the resolved value in FRACTION mode)
		ExpressionParam step_fraction_{0.0};    // [-] Step duration as fraction of total cycle time (0-1, used in FRACTION mode)
		std::map<size_t, EnergySource> energy_sources_; // LayerID -> EnergySource
		HeatingMedium heating_medium_{};  // jacket on the wall's outer face (off by default)

		// The step's medium when it is enabled, else nullptr (ambient exchange).
		const HeatingMedium* ActiveHeatingMedium() const {
			return heating_medium_.enabled_ ? &heating_medium_ : nullptr;
		}

		void EvaluateExpressions(const EvalContext& ctx) {
			step_time_.Evaluate(ctx);
			step_fraction_.Evaluate(ctx);
			for (auto& [id, src] : energy_sources_) src.EvaluateExpressions(ctx);
			heating_medium_.EvaluateExpressions(ctx);
		}

	private:
		friend class cereal::access;
		template<class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::make_nvp("Step_Time", step_time_),
				cereal::make_nvp("Top_Boundary_Condition", top_boundary_condition_),
				cereal::make_nvp("Bottom_Boundary_Condition", bottom_boundary_condition_),
				cereal::make_nvp("Step_Name", step_name_)
			);
			if (version >= 2) {
				ar(cereal::make_nvp("Energy_Sources", energy_sources_));
			}
			if (version >= 3) {
				ar(cereal::make_nvp("Step_Fraction", step_fraction_));
			}
			if (version >= 5) {
				ar(cereal::make_nvp("ID", id_));
			}
			// version < 5: id_ stays 0 and Cycle::EnsureIds assigns it on load.
			if (version >= 6) {
				ar(cereal::make_nvp("Heating_Medium", heating_medium_));
			}
			// version < 6: no medium (heating_medium_ disabled), as before.
		}
	};

} // end adsim namespace

// Class versions live next to the inline serialize templates: every
// translation unit that saves or loads a deck must see the same numbers.
CEREAL_CLASS_VERSION(adsim::HeatingMedium, 1);
CEREAL_CLASS_VERSION(adsim::Step, 6);
