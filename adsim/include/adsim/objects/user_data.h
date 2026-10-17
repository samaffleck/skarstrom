/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// AdSim includes
#include "adsim/objects/reactor.h"
#include "adsim/objects/domain.h"
#include "adsim/equations/equation.h"
#include "adsim/solver/ida_settings.h"
#include "adsim/solver/cycle.h"
#include "adsim/solver/experimental_data.h"
#include "adsim/solver/kpi_config.h"
#include "adsim/solver/cost_estimator.h"
#include "adsim/solver/solver_profile.h"
#include "adsim/math/user_variables.h"

// Cereal includes
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>

// STL includes
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>


namespace adsim {

  class SimulationInfo {
  public:
    std::string title{};
    std::string description{};
    std::string author{};
    size_t sim_version = 1;

  private: 
    friend class ::cereal::access;
    template<class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::make_nvp("Title", title),
        cereal::make_nvp("Description", description),
        cereal::make_nvp("Author", author),
        cereal::make_nvp("Version", sim_version)
      );
    }
  };

  // Stores all of our simulation data.
  // This is the top level data structure which contains everything to be serialised
  class UserData {
  public:
    bool PerformChecks(std::string& err_msg);
    void Initialise();

    // Single entry-point for refreshing every cached ExpressionParam value.
    // Resolves variable_registry dependencies (fixed-point iteration), then
    // cascades through the BC tree and other containers. Returns true on
    // success; on failure (cycle in deps, parse error) writes a message to
    // err_msg. Idempotent.
    bool EvaluateAllExpressions(std::string& err_msg);
    // Convenience overload (errors silently swallowed; for callers that don't care).
    void EvaluateAllExpressions() { std::string ignored; EvaluateAllExpressions(ignored); }

    // Every parameter (and user variable) whose CACHED value is inf/NaN, one
    // human-readable line each, e.g. `Bed geometry: "1/0" = inf`. Read-only:
    // nothing is re-evaluated, so call EvaluateAllExpressions first for fresh
    // values. Expressions that use a run-time state variable (T, P, c, y, v)
    // are skipped: they are evaluated per cell during the run, and their
    // placeholder value outside it is meaningless. PerformChecks fails when
    // this is non-empty.
    std::vector<std::string> FindNonFiniteExpressions();

    // Structural cycle edits. Each performs the Cycle operation AND remaps the
    // step-index-keyed KPI endpoint refs so they keep pointing at the same
    // logical endpoint. GUI and WASM edits must go through these rather than
    // editing cycle.step_list_ directly.
    //
    // Duplicates step `index` (see Cycle::DuplicateStep); the copy is
    // inserted at index + 1 and carries no KPI roles. Returns the new index,
    // or static_cast<size_t>(-1) if `index` is out of range.
    size_t DuplicateStep(size_t index);
    // Removes step `index` (disconnecting its partners) and drops its KPI refs.
    // No-op when out of range.
    void RemoveStep(size_t index);
    // Exchanges steps `a` and `b`; KPI refs follow their step. No-op when
    // either index is out of range.
    void SwapSteps(size_t a, size_t b);
    // Renames step `index`, rewriting KPI refs to its pressure-boundary
    // endpoint ("<step> Pressure Boundary Top|Bottom") and every buffer-tank
    // source (on any step) that names it as its source step. No-op when out
    // of range.
    void RenameStep(size_t index, const std::string& new_name);
    // Exchanges step `index`'s top and bottom boundary conditions (the "Swap
    // Boundary Conditions" button), re-stamps their location_, and remaps the
    // step's KPI refs to the other bed end (pressure names Top<->Bottom,
    // every ref's location flipped). No-op when out of range.
    void SwapStepBoundaries(size_t index);
    // Renames stream `stream_index` of the flow boundary on `side` of step
    // `step_index`, rewriting the KPI refs and buffer-tank sources (on any
    // step) that name it. No-op when the address does not name a stream.
    void RenameStream(size_t step_index, BoundaryLocation side, size_t stream_index,
                      const std::string& new_name);

    // Removes adsorbent `id` from the library AND erases it from every layer
    // composition (and its cost entry), so no reference to it dangles. A layer
    // left with an empty composition then fails Bed::PerformChecks. GUI and
    // WASM edits must remove adsorbents through this rather than calling
    // adsorbent_directory.Remove directly. No-op for an unknown id.
    void RemoveAdsorbent(AdsorbentID id);

    AdsorbentDirectory    adsorbent_directory;
    Bed                   bed;
    Fluid                 fluid;
    Equations             state;
    InitialConditions     ic;
    Cycle                 cycle;
    IDASettings           settings;
    ExperimentalData      exp_data_;
    KPIConfig             kpi_config;
    CostEstimator         cost_estimator;
    SimulationInfo        info;
    UserVariableRegistry  variable_registry;
    
    // Time offset for boundary conditions (time relative to step start)
    double step_start_time_ = 0.0;

    // Incremented on every load/deserialisation to invalidate cached references
    size_t load_generation_ = 0;

    // Transient wall-clock profile of the current solve (not serialised).
    SolverProfile profile_;

    // Jacobian builder of the live integrator, if a custom one is installed
    // (set by SolverData::Initialise; not serialised).
    BandedJacobian* jacobian_ = nullptr;

    // Message of the first exception caught at the IDA residual boundary (see
    // UpdateResiduals) since the last TakeResidualError(); empty if none. Not
    // serialised. The solver appends it to its failure log.
    std::string residual_error_;
    std::string TakeResidualError() { return std::exchange(residual_error_, std::string{}); }

  private: 
    friend class ::cereal::access;
    template<class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::make_nvp("Adsorbent_Directory", adsorbent_directory),
        cereal::make_nvp("Bed", bed),
        cereal::make_nvp("Fluid", fluid),
        cereal::make_nvp("Initial_Conditions", ic),
        cereal::make_nvp("Cycle", cycle),
        cereal::make_nvp("Solver_Settings", settings),
        cereal::make_nvp("Experimental_Data", exp_data_),
        cereal::make_nvp("System_Information", info)
      );
      if (version >= 1) {
        ar(cereal::make_nvp("KPI_Config", kpi_config));
      }
      if (version >= 2) {
        ar(cereal::make_nvp("Variable_Registry", variable_registry));
      }
      if (version >= 3) {
        ar(cereal::make_nvp("Cost_Estimator", cost_estimator));
      }
    }
  };

  // IDA residual callback. It is called from SUNDIALS C code, which an exception
  // must never unwind through (that can terminate the process or leave IDA's
  // internal state corrupt), so every exception is caught here, its message is
  // stored in UserData::residual_error_ for the solver log, and IDA gets an
  // unrecoverable error code (-1): it then returns IDA_RES_FAIL /
  // IDA_FIRST_RES_FAIL and the run fails cleanly with that message.
  inline int UpdateResiduals(sunrealtype t, N_Vector y, N_Vector yp, N_Vector res, void* user_data) {
    UserData* udata = static_cast<UserData*>(user_data);
    auto record = [udata](const char* what) {
      if (udata->residual_error_.empty()) udata->residual_error_ = what;
    };
    try {
      sunrealtype* y_data = N_VGetArrayPointer(y);
      sunrealtype* yp_data = N_VGetArrayPointer(yp);
      sunrealtype* res_data = N_VGetArrayPointer(res);
      auto& state = udata->state;

      ProfileTimer timer(udata->profile_.enabled, udata->profile_.residual_ns);
      ++udata->profile_.residual_calls;

      // Update order matters!
      state.UnpackBanded(y_data, yp_data);

      state.UpdateIndependentEquations(udata);
      state.UpdateDifferentialEquations(udata);

      double t_relative = t - udata->step_start_time_;
      state.UpdateBoundaryConditions(udata, t_relative);

      state.PackResidualsBanded(res_data);
    } catch (const std::exception& e) {
      record(e.what());
      return -1;
    } catch (...) {
      record("unknown exception");
      return -1;
    }

    return 0;
  }

	inline void SaveToJson(UserData& data, const std::string& file_name) {
		data.info.sim_version++;
    
    std::ofstream file(file_name);
    if (file.is_open()) {
      cereal::JSONOutputArchive archive(file);
      archive(data);
    }
	}

  // Shared tail of every UserData load. `deserialise(UserData&)` reads a deck
  // into a FRESH UserData; only once it (and the post-load fix-ups) succeed is
  // the result move-assigned into `data`. cereal deserialises in place, so
  // loading straight into `data` would leave it half-overwritten when a
  // malformed or incompatible deck throws part-way through (G-04). This gives
  // the strong exception guarantee: on throw, `data` is untouched.
  //
  // Move-assigning is safe for the self-referencing raw pointers: connection
  // `StreamData*` / BufferTankSource `summary_` point into boundary conditions
  // held by shared_ptr (heap, not moved), and std::vector move-assignment
  // steals the Step buffer without relocating Steps. Solver-side caches in
  // `state` are rebuilt by Equations::Initialise before every run.
  template<class Deserialise>
  inline void LoadUserDataTransactional(UserData& data, Deserialise&& deserialise) {
    UserData loaded;
    deserialise(loaded);
    loaded.kpi_config.FinaliseLegacyLoad(loaded.fluid.components_);
    // Give KPI endpoint refs saved before IDs existed their IDs (Cycle's load
    // already assigned step / endpoint IDs), and refresh cached names.
    loaded.kpi_config.ResolveIds(loaded.cycle.step_list_);
    // Keep the generation counter monotonic across loads so GUI caches keyed
    // on it are always invalidated (a fresh UserData starts at 0).
    loaded.load_generation_ = data.load_generation_ + 1;
    data = std::move(loaded);
  }

	// Throws on failure and leaves `data` unchanged (see LoadUserDataTransactional).
	inline void LoadFromJson(UserData& data, const std::string& file_name) {
		std::ifstream ifs(file_name);
		if (!ifs.is_open()) {
			throw std::runtime_error("LoadFromJson: could not open file " + file_name);
		}
		LoadUserDataTransactional(data, [&](UserData& loaded) {
			cereal::JSONInputArchive archive(ifs);
			archive(loaded);
		});
	}

  inline void SaveToString(UserData& data, std::string& out_string) {
    data.info.sim_version++;

    std::ostringstream oss;
    {
      cereal::JSONOutputArchive archive(oss);
      // using an NVP ensures the root is an object with a stable key
      archive(cereal::make_nvp("simulation", data));
    }
    out_string = oss.str();
  }

  // Throws on failure and leaves `data` unchanged (see LoadUserDataTransactional).
  inline void LoadFromString(UserData& data, const std::string& in_string) {
    if (in_string.empty()) {
        throw std::invalid_argument("LoadFromString: input JSON string is empty");
    }
    LoadUserDataTransactional(data, [&](UserData& loaded) {
      std::istringstream iss(in_string);
      cereal::JSONInputArchive archive(iss);
      archive(cereal::make_nvp("simulation", loaded));
    });
  }

} // end adsim namespace

CEREAL_CLASS_VERSION(adsim::UserData, 3);
