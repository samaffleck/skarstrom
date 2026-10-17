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

		void EvaluateExpressions(const EvalContext& ctx) {
			step_time_.Evaluate(ctx);
			step_fraction_.Evaluate(ctx);
			for (auto& [id, src] : energy_sources_) src.EvaluateExpressions(ctx);
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
		}
	};

} // end adsim namespace
