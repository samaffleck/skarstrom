/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

#include <string>
#include <unordered_map>

#include <cereal/types/string.hpp>
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>

#include "adsim/boundary_conditions/cycle_endpoint_summary.h"
#include "adsim/math/expression_map.h"


namespace adsim {

  class Connection;

  // Models a buffer tank between two steps: a consumer inlet pulls a
  // user-specified percentage of the source endpoint's most-recently-completed
  // total moles, distributed uniformly over the consumer step's duration, and
  // adopts the source's molar-flow-weighted average temperature and composition.
  //
  // Until the source step has completed at least one run (`summary_->is_valid`
  // is false), the consumer falls back to its user-specified "first-cycle
  // assumption" values (`temperature_`, `mole_fraction`, `boundary_flow_`).
  class BufferTankSource {
  public:
    void Activate(const std::string& source_step_name, const std::string& source_endpoint_name);
    void Deactivate();
    bool IsActive() const noexcept { return active_; }
    bool IsReady()  const noexcept { return active_ && summary_ != nullptr && summary_->is_valid; }
    bool PerformChecks(std::string& error_msg) const;

    // Consumer-side accessors. All assume IsReady(); callers should branch on
    // IsReady() and fall back to direct values when it returns false.
    double GetMolarFlow(double consumer_step_time) const;          // [mol/s]
    double GetTemperature() const;                                 // [K]
    double GetMoleFraction(const std::string& component) const;    // [-]

    // Resolved at solver start by Cycle::ResolveBufferTankSources from
    // (source_step_name_, source_endpoint_name_); not serialised.
    CycleEndpointSummary* summary_ = nullptr;

    ExpressionParam percentage_{100.0}; // 0..100, validated in PerformChecks
    // The source endpoint. The IDs (Step::id_, endpoint_id_) are the
    // reference; the names are display labels that
    // Cycle::ResolveBufferTankReferences refreshes from the IDs, so renaming
    // the source step or stream does not break the tank. A source set by name
    // only (IDs 0: an old deck, Activate, the web document) is resolved by
    // name once and gets its IDs then.
    std::string source_step_name_{};
    std::string source_endpoint_name_{};
    size_t source_step_id_ = 0;
    size_t source_endpoint_id_ = 0;

    void EvaluateExpressions(const EvalContext& ctx) { percentage_.Evaluate(ctx); }

  private:
    bool active_ = false;

    friend class cereal::access;

    template <class Archive>
    void save(Archive& ar, std::uint32_t const /*version*/) const {
      ar(
        cereal::make_nvp("Active", active_),
        cereal::make_nvp("Percentage", percentage_),
        cereal::make_nvp("Source_Step_Name", source_step_name_),
        cereal::make_nvp("Source_Endpoint_Name", source_endpoint_name_),
        cereal::make_nvp("Source_Step_ID", source_step_id_),
        cereal::make_nvp("Source_Endpoint_ID", source_endpoint_id_)
      );
    }

    template <class Archive>
    void load(Archive& ar, std::uint32_t const version) {
      LoadFields(ar, version);
      if (version >= 3) {
        ar(
          cereal::make_nvp("Source_Step_ID", source_step_id_),
          cereal::make_nvp("Source_Endpoint_ID", source_endpoint_id_)
        );
      }
      // version < 3: IDs stay 0 and are resolved from the names on load.
    }

    template <class Archive>
    void LoadFields(Archive& ar, std::uint32_t const version) {
      if (version < 2) {
        // v1: Percentage was a plain double.
        double pct = 100.0;
        ar(
          cereal::make_nvp("Active", active_),
          cereal::make_nvp("Percentage", pct),
          cereal::make_nvp("Source_Step_Name", source_step_name_),
          cereal::make_nvp("Source_Endpoint_Name", source_endpoint_name_)
        );
        percentage_ = pct; // ExpressionParam::operator=(double) re-stringifies
      } else {
        ar(
          cereal::make_nvp("Active", active_),
          cereal::make_nvp("Percentage", percentage_),
          cereal::make_nvp("Source_Step_Name", source_step_name_),
          cereal::make_nvp("Source_Endpoint_Name", source_endpoint_name_)
        );
      }
    }
  };

  // Stream-side inlet temperature / composition of an inlet endpoint (every
  // IBoundaryStream kind and PressureBoundaryCondition), shared by the
  // residuals and Log() so the two can never pick different sources:
  //   1. an active buffer tank: its source summary once ready, otherwise the
  //      user-specified first-cycle assumption;
  //   2. an interacting-step connection with logged data: the connected
  //      stream's value at `time`;
  //   3. otherwise the user-specified value (also the first-cycle assumption
  //      of a connection that has no data yet).
  // A component missing from `assumed_mole_fraction` reads as 0.
  double SelectInletTemperature(const BufferTankSource& buffer_tank,
                                const Connection& connection,
                                double assumed_temperature,
                                double time);
  double SelectInletMoleFraction(const BufferTankSource& buffer_tank,
                                 const Connection& connection,
                                 const ExpressionMap& assumed_mole_fraction,
                                 const std::string& component,
                                 double time);

  // Canonical endpoint name used when the source is a PressureBoundaryCondition.
  // FlowBoundaryCondition streams use their stream->name_ as-is.
  inline constexpr const char* kPressureBoundaryEndpointName = "Pressure Boundary";

} // end adsim namespace

CEREAL_CLASS_VERSION(adsim::BufferTankSource, 3);
