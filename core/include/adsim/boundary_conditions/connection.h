/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

#include "adsim/boundary_conditions/flow_boundary_value.h"

// STL includes
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <atomic>

// Cereal includes
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>

namespace adsim {
  
  class Step;
  class StreamData;
  class IBoundaryStream;
  class PressureBoundaryCondition;

  static const char* ValveConnectionList[] = {
		"Not Connected",	// 0
		"Connected"		// 1 
	};

	enum class ValveConnection {
		NOT_CONNECTED = 0,
		CONNECTED = 1
	};

  class Connection {
  public:
    void MakeConnection(size_t id, 
      IBoundaryStream& connected_stream);
    void MakeConnection(size_t id,
      PressureBoundaryCondition& connected_pressure_boundary);
    void Disconnect();
    void DisconnectKeepID();
    bool PerformChecks(std::string& error_msg) const;
    bool IsConnected() const noexcept { return type_ == ValveConnection::CONNECTED; }
    bool CanConnect() const noexcept { return !IsConnected(); }
    const size_t GetID() const noexcept { return pair_id_; } 

    StreamData* stream_ = nullptr;
    
  private:
    ValveConnection type_ = ValveConnection::NOT_CONNECTED;
    size_t pair_id_ = 0; 

    friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::make_nvp("Connection_Type", type_),
				cereal::make_nvp("Pair_ID", pair_id_)
			);
		}
  };

  struct ConnectionEndpoint {
    Connection* conn = nullptr;
    Step* step = nullptr;
    IBoundaryStream* stream = nullptr;
    PressureBoundaryCondition* p_boundary = nullptr;
  };

  // Gas passes through a connection only from the higher to the lower
  // pressure.
  //
  // A connected set-flow stream (BoundaryStream, no device) replays the
  // flow its partner logged. This weight, of the pressure difference
  // dp = upstream - downstream [Pa] at the same instant, scales the
  // replayed flow: 0 for dp <= 0, 1 for dp >= kConnectionPressureBand, and
  // the C1 cubic 3x^2 - 2x^3 (x = dp / band) between. A hard switch would
  // leave a kink in the residual that stalls Newton as the flow decays
  // through zero.
  //
  // A connected pressure boundary caps its imposed pressure one band inside
  // the partner's logged pressure (ConnectionCappedPressure), so the
  // partner's weight is 1 once the replay has settled and the moles through
  // the connection balance.
  //
  // The band is wide on purpose (1 kPa, about 1 % of an atmosphere). Both
  // closures read the partner's previous run, so each cycle answers the
  // other's last one. A narrow band makes the closure on/off against a
  // pressure that moves by more than the band from one cycle to the next,
  // and the pair hunts instead of settling (measured on Perez2019_D1 with
  // its vents' check valves on: no cyclic steady state in 250 cycles and
  // the light-product moles out of balance by 44 % with 10 Pa; balanced to
  // 1 % with 1 kPa).
  inline constexpr double kConnectionPressureBand = 1000.0;  // [Pa]
  double ConnectionFlowWeight(double dp);

  // The pressure a connected pressure boundary imposes: `target` (its own
  // profile) kept on the right side of `partner` (the partner's logged
  // pressure at the same instant). An INLET takes gas from the partner, so
  // it cannot rise above it (`inlet` true: the result is at most
  // `partner`); an OUTLET gives gas to it, so it cannot fall below it. The
  // target is returned exactly while it is at least kConnectionPressureBand
  // inside that bound, the bound at or beyond it, and a C1 cubic blend
  // between (never past the bound).
  double ConnectionCappedPressure(double target, double partner, bool inlet);

  size_t GenerateID();
  ConnectionEndpoint MakeEndpoint(Step& step, IBoundaryStream& stream);
  ConnectionEndpoint MakeEndpoint(Step& step, PressureBoundaryCondition& pressure);
  void MutualDisconnect(ConnectionEndpoint& a, ConnectionEndpoint& b);
  bool MakeMutualConnection(ConnectionEndpoint& a, ConnectionEndpoint& b);

} // end adsim namespace

// Here rather than in the .cpp, so that every archive that serialises this
// class sees the version (a translation unit without it would use 0).
CEREAL_CLASS_VERSION(adsim::Connection, 1);
