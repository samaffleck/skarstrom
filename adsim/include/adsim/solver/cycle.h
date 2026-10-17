/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <optional>
#include <vector>

// Cereal includes
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>

// adsim includes
#include "adsim/boundary_conditions/step.h"
#include "adsim/boundary_conditions/connection.h"
#include "adsim/objects/fluid.h"
#include "adsim/math/expression_param.h"


namespace adsim {

  typedef std::vector<std::pair<ConnectionEndpoint, ConnectionEndpoint>> EndpointPairList;

  // Free functions
  EndpointPairList GetConnectedEndpoints(std::vector<Step>& step_list);

  // Disconnects `conn` together with its partner: every connected endpoint in
  // `step_list` that shares conn's pair id. A Connection holds a raw pointer
  // into its partner's StreamData, so this must run before an endpoint is
  // destroyed. No-op when `conn` is not connected.
  void DisconnectWithPartner(std::vector<Step>& step_list, Connection& conn);

  // Disconnects every endpoint carried by `bc` (its pressure endpoint, or each
  // of its flow streams) together with its partner. Null-safe.
  void DisconnectBoundaryCondition(std::vector<Step>& step_list, IBoundaryCondition* bc);

  enum class StepTimeMode {
    ABSOLUTE,  // Each step stores its duration in seconds
    FRACTION   // Each step stores a fraction (0-1) of total_cycle_time_
  };

  static const char* StepTimeModeList[] = {
    "Absolute",
    "Fraction"
  };

	class Cycle {
	public:
		void AddStep(const std::string& step_name);
		void AddStep(const std::string& step_name, size_t position);
		// Disconnects the step's endpoints from their partners, then erases it.
		void RemoveStep(const std::string& step_name);

		// Safe replace/remove of connectable endpoints: each first disconnects the
		// outgoing endpoint(s) from their partners so no surviving Connection is
		// left CONNECTED with a pointer into destroyed StreamData.
		void ReplaceBoundaryCondition(std::shared_ptr<IBoundaryCondition>& slot,
		                              std::shared_ptr<IBoundaryCondition> replacement);
		void ReplaceStream(std::unique_ptr<IBoundaryStream>& slot,
		                   std::unique_ptr<IBoundaryStream> replacement);
		// Erases bc.streams_[index]; no-op when out of range.
		void RemoveStream(FlowBoundaryCondition& bc, size_t index);

		// Deep-copies the step at `index` and inserts the copy immediately after it.
		// The duplicate gets a unique name ("<name> (copy)", "<name> (copy 2)", ...)
		// and all of its connections are cleared so it starts fully detached.
		// Returns the index of the newly inserted step, or static_cast<size_t>(-1)
		// if `index` is out of range.
		size_t DuplicateStep(size_t index);
		bool PerformChecks(const Fluid& fluid, std::string& error_message);
		double GetTotalSimulationTime() const;
		double GetCycleTime() const;
		Step& GetStep(const std::string& step_name);
		Step* FindStep(const std::string& step_name);
		void ResolveStepTimes();

		// Walks every consumer with an active BufferTankSource and binds its
		// summary_ pointer to the source step's matching endpoint summary.
		// Logs (via error_message) any unresolved references but does not
		// throw — PerformChecks is the gate that hard-fails on missing sources.
		void ResolveBufferTankSources();

		// A BufferTankSource names its source by (step name, endpoint name), so
		// renaming either must rewrite every source that names it, on any step.
		// Callers should use UserData::RenameStep / RenameStream, which decide
		// whether the references follow the renamed step/stream.
		void RenameBufferTankSourceStep(const std::string& old_name, const std::string& new_name);
		void RenameBufferTankSourceEndpoint(const std::string& step_name,
		                                    const std::string& old_endpoint,
		                                    const std::string& new_endpoint);

		// ---- Identity ---------------------------------------------------------
		// Steps and endpoints (flow streams and pressure boundaries) carry
		// stable IDs (Step::id_, endpoint_id_) so references survive renames
		// and reordering. EnsureIds gives every object whose ID is 0 (a deck
		// saved before IDs existed, or an object created in code) the next
		// free ID, and re-numbers a repeated ID (a copied step) on its later
		// occurrence. Assignment walks the steps in order, each step's bottom
		// boundary before its top, streams in order, so an old deck always
		// gets the same IDs: steps 1..N, endpoints 1..M. Called on load and by
		// every structural edit here; idempotent.
		void EnsureIds();

		// One endpoint of the cycle, located by ID.
		struct EndpointHandle {
			size_t id = 0;
			size_t step_index = 0;
			Step* step = nullptr;
			IBoundaryCondition* boundary = nullptr;
			IBoundaryStream* stream = nullptr;              // flow stream, or
			PressureBoundaryCondition* pressure = nullptr;  // pressure boundary
			BoundaryLocation location = BoundaryLocation::BOTTOM;
			// KPI / BoundaryLog name: the stream name, or
			// Step::PressureEndpointName(step, location) for a pressure boundary.
			std::string log_name;
			// Buffer-tank / cycle-summary name: the stream name, or
			// kPressureBoundaryEndpointName.
			std::string summary_name;
		};
		// Every endpoint, in ID-assignment order.
		std::vector<EndpointHandle> Endpoints();
		std::optional<EndpointHandle> FindEndpoint(size_t endpoint_id);
		// Resolves every active buffer tank's source: by ID when set (and
		// refreshes the source names from it), otherwise by name (and stores the
		// IDs). A source that resolves neither way is left for PerformChecks to
		// report.
		void ResolveBufferTankReferences();
		Step* FindStepById(size_t step_id);
		// Index of the step with this ID, or -1.
		int StepIndexById(size_t step_id) const;

		int number_of_cycles_ = 1;
		StepTimeMode step_time_mode_ = StepTimeMode::ABSOLUTE;
		ExpressionParam total_cycle_time_{100.0}; // [s] Used in FRACTION mode
		std::vector<Step> step_list_;

		void EvaluateExpressions(const EvalContext& ctx);

	private:
		friend class cereal::access;
		template<class Archive>
		void save(Archive& ar, const std::uint32_t version) const {
      ar(
        cereal::make_nvp("Number_Of_Cycles", number_of_cycles_),
        cereal::make_nvp("Step_List", step_list_)
      );
      if (version >= 2) {
        ar(
          cereal::make_nvp("Step_Time_Mode", static_cast<int>(step_time_mode_)),
          cereal::make_nvp("Total_Cycle_Time", total_cycle_time_)
        );
      }
    }

    template<class Archive>
    void load(Archive& ar, const std::uint32_t version) {
      ar(
        cereal::make_nvp("Number_Of_Cycles", number_of_cycles_),
        cereal::make_nvp("Step_List", step_list_)
      );
      if (version >= 2) {
        int mode = 0;
        ar(
          cereal::make_nvp("Step_Time_Mode", mode),
          cereal::make_nvp("Total_Cycle_Time", total_cycle_time_)
        );
        step_time_mode_ = static_cast<StepTimeMode>(mode);
      }

      // After all Steps/Connections are deserialized, rebuild the raw-pointer links.
      EnsureIds();
      RebuildConnectionsFromIDs();
      ResolveBufferTankSources();
    }

    void RebuildConnectionsFromIDs() {
      // We have to get all endpoints BEFORE calling DisconnectKeepID!
      auto connected_endpoints = GetConnectedEndpoints(step_list_);

      for (Step& step : step_list_) {
        auto clear_bc = [&](IBoundaryCondition* bc) {
          if (auto* flow = dynamic_cast<FlowBoundaryCondition*>(bc)) {
            for (auto& stream : flow->streams_) {
              stream->connection_.DisconnectKeepID();        // sets NOT_CONNECTED, keeps pair_id_
            }
          } else if (auto* p = dynamic_cast<PressureBoundaryCondition*>(bc)) {
            p->connection_.DisconnectKeepID();               // sets NOT_CONNECTED, keeps pair_id_
          }
        };
        clear_bc(step.bottom_boundary_condition_.get());
        clear_bc(step.top_boundary_condition_.get());
      }

      for (auto& [e1, e2] : connected_endpoints) {
        MakeMutualConnection(e1, e2);
      }
    }

	};

  inline std::vector<ConnectionEndpoint> GetAvailableEndpoints(std::vector<Step>& step_list) {
    std::vector<ConnectionEndpoint> endpoint_list{};

    auto append_endpoint_list = [&endpoint_list](Step& step, IBoundaryCondition* bc) {
      if (auto* flow = dynamic_cast<FlowBoundaryCondition*>(bc)) {
        for (auto& stream : flow->streams_) {
          if (!stream->connection_.CanConnect()) continue;
          ConnectionEndpoint new_endpoint;
          new_endpoint.conn = &stream->connection_;
          new_endpoint.step = &step; 
          new_endpoint.stream = stream.get();
          endpoint_list.emplace_back(new_endpoint);
        }
      } else if (auto* pressure = dynamic_cast<PressureBoundaryCondition*>(bc)) {
        if (pressure->connection_.CanConnect()) {
          ConnectionEndpoint new_endpoint;
          new_endpoint.conn = &pressure->connection_;
          new_endpoint.step = &step;
          new_endpoint.p_boundary = pressure;
          endpoint_list.emplace_back(new_endpoint);
        }
      }
    };

    for (Step& step : step_list) {
      append_endpoint_list(step, step.bottom_boundary_condition_.get());
      append_endpoint_list(step, step.top_boundary_condition_.get());
    }

    return endpoint_list;
  } 

  inline std::vector<IBoundaryStream*> GetAllBoundaryStreams(const std::vector<Step>& step_list) {
    std::vector<IBoundaryStream*> streams;

    auto check_boundary = [&](IBoundaryCondition* bc) {
      if (auto* fbc = dynamic_cast<FlowBoundaryCondition*>(bc)) {
        for (auto& stream : fbc->streams_) {
          streams.push_back(stream.get());
        }
      }
    };

    for (const auto& step : step_list) {
      check_boundary(step.bottom_boundary_condition_.get());
      check_boundary(step.top_boundary_condition_.get());
    }

    return streams;
  }

  struct BoundaryData {
    const StreamData* stream_data = nullptr;
    std::string boundary_name{};
  };

  inline std::vector<BoundaryData> GetAllBoundaryStreams(const Step& step) {
    std::vector<BoundaryData> streams;

    auto check_boundary = [&](const IBoundaryCondition* bc) {
      if (!bc) return;

      if (auto fbc = dynamic_cast<const FlowBoundaryCondition*>(bc)) {
        for (const auto& stream : fbc->streams_) {
          BoundaryData data;
          data.stream_data = &stream->stream_data_;
          data.boundary_name = stream->name_;
          streams.emplace_back(data);
        }
      } else if (auto pbc = dynamic_cast<const PressureBoundaryCondition*>(bc)) {
        BoundaryData data;
        data.stream_data = &pbc->stream_data_;
        std::string location = BoundaryLocationList[(int)bc->location_];
        data.boundary_name = step.step_name_ + " Pressure Boundary " + location; 
        streams.emplace_back(data);
      }
    };

    check_boundary(step.bottom_boundary_condition_.get());
    check_boundary(step.top_boundary_condition_.get());
    return streams;
  }

  inline void SaveCycleToJSON(const std::string& file_name, const Cycle& cycle) {
    std::ofstream file(file_name);
    if (file.is_open()) {
      cereal::JSONOutputArchive archive(file);
      archive(cycle);
    }
  }

  // Deserialises into a temporary and only replaces `cycle` on success, so a
  // malformed file throws and leaves `cycle` unchanged (cereal loads in place).
  // Moving is safe for the connection / buffer-tank raw pointers: they target
  // boundary conditions held by shared_ptr, and moving step_list_ steals the
  // Step buffer without relocating Steps.
  inline void LoadCycleFromJSON(const std::string& file_name, Cycle& cycle) {
    std::ifstream ifs(file_name);
    if (!ifs.is_open()) {
      throw std::runtime_error("LoadCycleFromJSON: could not open file " + file_name);
    }
    Cycle loaded;
    {
      cereal::JSONInputArchive archive(ifs);
      archive(loaded);
    }
    cycle = std::move(loaded);
  }

} // end adsim namespace
