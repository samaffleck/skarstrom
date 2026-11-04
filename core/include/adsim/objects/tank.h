/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Dynamic tank: a well-mixed, constant-volume gas vessel whose species moles
// and temperature are unknowns of the bed's DAE in every step (see
// docs/tanks.md). Bed endpoints name a tank through their TankPort; the tank's
// own external streams (make-up, vent valve, back-pressure regulator,
// pressure hold) connect it to the outside.
//
// Balances (n_i the species moles, N = sum n_i, P = N R T / V):
//   dn_i/dt = sum over ports and external streams of the molar flow of i
//             into the tank,
//   d/dt( sum_i n_i c_v,i T ) = sum_p sum_i F_p,i c_p,i T_p + Q,
// with T_p the temperature of the gas that crosses p (the tank's own for gas
// leaving it) and Q the heat from the ambient or the wall node.

#include <memory>
#include <string>
#include <vector>

#include <cereal/access.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/math/expression_map.h"
#include "adsim/math/expression_param.h"
#include "adsim/objects/tank_state.h"

namespace adsim {

  class Fluid;

  enum class TankThermalMode {
    ISOTHERMAL = 0,  // gas held at initial_temperature_ (or the run's isothermal temperature)
    ADIABATIC = 1,   // no heat exchange
    HEAT_LOSS = 2,   // Q = UA (T_amb - T)
    WALL = 3         // wall node: m_w c_w dT_w/dt = h_in A (T - T_w) - h_out A (T_w - T_amb)
  };

  static const char* TankThermalModeList[] = { "isothermal", "adiabatic", "heat_loss", "wall" };

  enum class TankStreamKind {
    MAKE_UP = 0,    // inflow at a set molar flow or through a linear valve from a header
    VALVE = 1,      // linear valve to a downstream pressure, optional check valve
    REGULATOR = 2,  // back-pressure regulator, opening ramps from 0 at P_set to 1 at P_full
    HOLD = 3        // overflow at a held pressure: a regulator with P_full = P_set (1 + band)
  };

  static const char* TankStreamKindList[] = { "make_up", "valve", "regulator", "hold" };

  enum class TankMakeUpMode {
    SET_FLOW = 0,  // F = molar_flow_
    VALVE = 1      // F = Cv max(P_header - P, 0) / 1e5 (no backflow into the header)
  };

  namespace detail {
    // Smoothed max(x, 0): (x + w)^2 / (4 w) on |x| < w, exactly max(x, 0)
    // outside; continuously differentiable and never negative.
    double SmoothPositive(double x, double width);
  }

  // Pressure width [Pa] over which a tank stream's max(P_up - P_down, 0)
  // (check valve, regulator, make-up valve) is smoothed.
  inline constexpr double kTankValveSmoothingPressure = 1.0;
  // Width, as a fraction of the opening band, over which the regulator's
  // clamp((P - P_set) / (P_full - P_set), 0, 1) is smoothed.
  inline constexpr double kTankRegulatorSmoothing = 0.01;

  // An external stream of a tank. Every stream logs a StreamData per step
  // (molar flow positive OUT of the tank, see StreamData::GetOutflow) and has
  // an endpoint ID in the same family as the bed endpoints (registry paths
  // "endpoint[<id>]").
  class ITankStream {
  public:
    virtual ~ITankStream() = default;
    virtual TankStreamKind GetKind() const = 0;
    // Molar flow into the tank [mol/s] (negative: out of it) at tank state s.
    virtual double FlowIn(const TankLiveState& s) const = 0;
    virtual void EvaluateExpressions(const EvalContext& ctx) {
      temperature_.Evaluate(ctx);
      EvaluateExpressionMap(mole_fraction_, ctx);
    }
    // Non-const: reconciles the composition map with the fluid's components.
    virtual bool PerformChecks(std::string& error_msg, const Fluid& fluid);
    // True when gas can flow into the tank through this stream, so its own
    // composition and temperature must be valid.
    virtual bool CanFlowIn() const { return false; }

    std::string name_{"Stream"};
    // Stable identity, shared with the bed endpoints' ID space (0 = not yet
    // assigned; see UserData::SyncTankEndpointIds).
    size_t endpoint_id_ = 0;
    // Gas that enters the tank through this stream: the make-up gas, or the
    // gas beyond a valve without a check valve when it flows back.
    ExpressionMap mole_fraction_{};
    ExpressionParam temperature_{293.15};
    // Transient log of the current step (time since the step start); never
    // serialised.
    StreamData stream_data_{};

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/) {
      ar(
        cereal::make_nvp("Name", name_),
        cereal::make_nvp("Endpoint_ID", endpoint_id_),
        cereal::make_nvp("Mole_Fraction", mole_fraction_),
        cereal::make_nvp("Temperature", temperature_)
      );
    }
  };

  class TankMakeUp : public ITankStream {
  public:
    TankStreamKind GetKind() const override { return TankStreamKind::MAKE_UP; }
    double FlowIn(const TankLiveState& s) const override;
    bool CanFlowIn() const override { return true; }
    bool PerformChecks(std::string& error_msg, const Fluid& fluid) override;
    void EvaluateExpressions(const EvalContext& ctx) override {
      ITankStream::EvaluateExpressions(ctx);
      molar_flow_.Evaluate(ctx);
      cv_.Evaluate(ctx);
      header_pressure_.Evaluate(ctx);
    }

    TankMakeUpMode mode_ = TankMakeUpMode::SET_FLOW;
    ExpressionParam molar_flow_{0.0};          // [mol/s] (SET_FLOW)
    ExpressionParam cv_{1.0};                  // [mol/(s bar)] (VALVE)
    ExpressionParam header_pressure_{101325.0};  // [Pa] (VALVE)

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/) {
      ar(
        cereal::virtual_base_class<ITankStream>(this),
        cereal::make_nvp("Mode", mode_),
        cereal::make_nvp("Molar_Flow", molar_flow_),
        cereal::make_nvp("Cv", cv_),
        cereal::make_nvp("Header_Pressure", header_pressure_)
      );
    }
  };

  class TankValve : public ITankStream {
  public:
    TankStreamKind GetKind() const override { return TankStreamKind::VALVE; }
    double FlowIn(const TankLiveState& s) const override;
    bool CanFlowIn() const override { return !check_valve_; }
    bool PerformChecks(std::string& error_msg, const Fluid& fluid) override;
    void EvaluateExpressions(const EvalContext& ctx) override {
      ITankStream::EvaluateExpressions(ctx);
      cv_.Evaluate(ctx);
      downstream_pressure_.Evaluate(ctx);
    }

    ExpressionParam cv_{1.0};                      // [mol/(s bar)]
    ExpressionParam downstream_pressure_{101325.0};  // [Pa]
    bool check_valve_ = true;                      // true: no flow back into the tank

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/) {
      ar(
        cereal::virtual_base_class<ITankStream>(this),
        cereal::make_nvp("Cv", cv_),
        cereal::make_nvp("Downstream_Pressure", downstream_pressure_),
        cereal::make_nvp("Check_Valve", check_valve_)
      );
    }
  };

  // Back-pressure regulator (toPSAil's pressure-driven product mode):
  //   F_out = Cv clamp((P - P_set) / (P_full - P_set), 0, 1) max(P - P_down, 0) / 1e5,
  // with the clamp and the max smoothed (kTankRegulatorSmoothing,
  // kTankValveSmoothingPressure).
  class TankRegulator : public ITankStream {
  public:
    TankStreamKind GetKind() const override { return TankStreamKind::REGULATOR; }
    double FlowIn(const TankLiveState& s) const override;
    bool PerformChecks(std::string& error_msg, const Fluid& fluid) override;
    void EvaluateExpressions(const EvalContext& ctx) override {
      ITankStream::EvaluateExpressions(ctx);
      cv_.Evaluate(ctx);
      set_pressure_.Evaluate(ctx);
      full_open_pressure_.Evaluate(ctx);
      downstream_pressure_.Evaluate(ctx);
    }
    // The regulator law (shared with TankHold).
    static double Outflow(double P, double cv, double P_set, double P_full, double P_down);

    ExpressionParam cv_{1.0};                        // [mol/(s bar)] at full opening
    ExpressionParam set_pressure_{101325.0};         // [Pa] opening starts
    ExpressionParam full_open_pressure_{111325.0};   // [Pa] fully open
    ExpressionParam downstream_pressure_{101325.0};  // [Pa]

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/) {
      ar(
        cereal::virtual_base_class<ITankStream>(this),
        cereal::make_nvp("Cv", cv_),
        cereal::make_nvp("Set_Pressure", set_pressure_),
        cereal::make_nvp("Full_Open_Pressure", full_open_pressure_),
        cereal::make_nvp("Downstream_Pressure", downstream_pressure_)
      );
    }
  };

  // Overflow at a held pressure (toPSAil's flow-driven product mode): a
  // regulator with P_full = P_set (1 + band_). While the net inflow is
  // positive and below Cv (P_full - P_down) / 1e5 the tank pressure stays
  // within band_ P_set above the set point and the draw-off equals the net
  // inflow; while it is negative the tank pressure falls (nothing flows in).
  class TankHold : public ITankStream {
  public:
    TankStreamKind GetKind() const override { return TankStreamKind::HOLD; }
    double FlowIn(const TankLiveState& s) const override;
    bool PerformChecks(std::string& error_msg, const Fluid& fluid) override;
    void EvaluateExpressions(const EvalContext& ctx) override {
      ITankStream::EvaluateExpressions(ctx);
      set_pressure_.Evaluate(ctx);
      band_.Evaluate(ctx);
      cv_.Evaluate(ctx);
      downstream_pressure_.Evaluate(ctx);
    }

    static constexpr double kDefaultBand = 1e-4;
    ExpressionParam set_pressure_{101325.0};         // [Pa]
    ExpressionParam band_{kDefaultBand};             // [-] (P_full - P_set) / P_set
    ExpressionParam cv_{10.0};                       // [mol/(s bar)] at full opening
    ExpressionParam downstream_pressure_{101325.0};  // [Pa]

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/) {
      ar(
        cereal::virtual_base_class<ITankStream>(this),
        cereal::make_nvp("Set_Pressure", set_pressure_),
        cereal::make_nvp("Band", band_),
        cereal::make_nvp("Cv", cv_),
        cereal::make_nvp("Downstream_Pressure", downstream_pressure_)
      );
    }
  };

  std::unique_ptr<ITankStream> CreateTankStream(TankStreamKind kind);

  class Tank {
  public:
    Tank() = default;
    explicit Tank(std::string name) : name_(std::move(name)) {}

    // Checks the tank against the fluid; reconciles the composition maps with
    // the fluid's components.
    bool PerformChecks(std::string& error_msg, const Fluid& fluid);
    void EvaluateExpressions(const EvalContext& ctx);
    // Initial moles of each component, sum_i = P0 V / (R T0).
    double InitialTotalMoles() const;
    ITankStream* FindStream(size_t endpoint_id);

    std::string name_{"Tank"};
    size_t id_ = 0;  // stable ID (registry paths "tank[<id>]"); set by TankDirectory::Add
    ExpressionParam volume_{0.01};                   // [m3]
    ExpressionParam initial_pressure_{101325.0};     // [Pa]
    ExpressionParam initial_temperature_{293.15};    // [K]
    ExpressionMap initial_mole_fraction_{};          // [-]
    // Number of identical beds sharing the tank, staggered by T_c / N. Only
    // N = 1 (no replay of other beds) is implemented; PerformChecks rejects more.
    int beds_ = 1;
    TankThermalMode thermal_ = TankThermalMode::ISOTHERMAL;
    ExpressionParam ambient_temperature_{293.15};  // [K] (HEAT_LOSS, WALL)
    ExpressionParam ua_{0.0};                      // [W/K] (HEAT_LOSS)
    ExpressionParam wall_area_{1.0};               // [m2] (WALL)
    ExpressionParam h_in_{10.0};                   // [W/(m2 K)] gas to wall (WALL)
    ExpressionParam h_out_{10.0};                  // [W/(m2 K)] wall to ambient (WALL)
    ExpressionParam wall_heat_capacity_{1000.0};   // [J/K] m_w c_w (WALL)
    ExpressionParam initial_wall_temperature_{293.15};  // [K] (WALL)
    std::vector<std::unique_ptr<ITankStream>> streams_{};

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/) {
      ar(
        cereal::make_nvp("Name", name_),
        cereal::make_nvp("ID", id_),
        cereal::make_nvp("Volume", volume_),
        cereal::make_nvp("Initial_Pressure", initial_pressure_),
        cereal::make_nvp("Initial_Temperature", initial_temperature_),
        cereal::make_nvp("Initial_Mole_Fraction", initial_mole_fraction_),
        cereal::make_nvp("Beds", beds_),
        cereal::make_nvp("Thermal_Mode", thermal_),
        cereal::make_nvp("Ambient_Temperature", ambient_temperature_),
        cereal::make_nvp("UA", ua_),
        cereal::make_nvp("Wall_Area", wall_area_),
        cereal::make_nvp("H_In", h_in_),
        cereal::make_nvp("H_Out", h_out_),
        cereal::make_nvp("Wall_Heat_Capacity", wall_heat_capacity_),
        cereal::make_nvp("Initial_Wall_Temperature", initial_wall_temperature_),
        cereal::make_nvp("Streams", streams_)
      );
    }
  };

  // The deck's tanks, in creation order, keyed by a stable ID.
  class TankDirectory {
  public:
    // Adds a tank and returns its ID (never reused within the directory).
    size_t Add(const std::string& name = "Tank");
    void Remove(size_t id);
    Tank* Find(size_t id);
    const Tank* Find(size_t id) const;
    Tank* FindByName(const std::string& name);
    bool Empty() const { return tanks_.empty(); }
    size_t Size() const { return tanks_.size(); }

    std::vector<Tank> tanks_{};

  private:
    size_t next_id_ = 0;

    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/) {
      ar(
        cereal::make_nvp("Tanks", tanks_),
        cereal::make_nvp("Next_ID", next_id_)
      );
    }
  };

} // end adsim namespace

CEREAL_CLASS_VERSION(adsim::ITankStream, 1);
CEREAL_CLASS_VERSION(adsim::TankMakeUp, 1);
CEREAL_CLASS_VERSION(adsim::TankValve, 1);
CEREAL_CLASS_VERSION(adsim::TankRegulator, 1);
CEREAL_CLASS_VERSION(adsim::TankHold, 1);
CEREAL_CLASS_VERSION(adsim::Tank, 1);
CEREAL_CLASS_VERSION(adsim::TankDirectory, 1);
