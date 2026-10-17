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

  size_t GenerateID();
  ConnectionEndpoint MakeEndpoint(Step& step, IBoundaryStream& stream);
  ConnectionEndpoint MakeEndpoint(Step& step, PressureBoundaryCondition& pressure);
  void MutualDisconnect(ConnectionEndpoint& a, ConnectionEndpoint& b);
  bool MakeMutualConnection(ConnectionEndpoint& a, ConnectionEndpoint& b);

} // end adsim namespace
