/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/boundary_conditions/connection.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/boundary_conditions/boundary_conditions.h"


namespace adsim {

  void Connection::MakeConnection(
    size_t id, 
    IBoundaryStream& connected_stream  
  ) {
    type_ = ValveConnection::CONNECTED;
    stream_ = &connected_stream.stream_data_;
    pair_id_ = id;
  }

  void Connection::MakeConnection(
    size_t id, 
    PressureBoundaryCondition& connected_pressure_boundary  
  ) {
    type_ = ValveConnection::CONNECTED;
    stream_ = &connected_pressure_boundary.stream_data_;
    pair_id_ = id;
  }

  // Keeps id (only used during serialisation)
  void Connection::DisconnectKeepID() {
    type_ = ValveConnection::NOT_CONNECTED;
    stream_ = nullptr;
  }

  // Resets id
  void Connection::Disconnect() {
    type_ = ValveConnection::NOT_CONNECTED;
    stream_ = nullptr;
    pair_id_ = 0;
  }

  bool Connection::PerformChecks(std::string& error_msg) const {
    if (type_ == ValveConnection::CONNECTED) {
      if (stream_ == nullptr) {
        error_msg += "No connected boundary selected.\n";
        return false;
      }
    }

    return true;
  }

  // Functions
  size_t GenerateID() {
    static std::atomic<size_t> counter{1};
    return counter.fetch_add(1, std::memory_order_relaxed);
  }
  
  ConnectionEndpoint MakeEndpoint(Step& step, IBoundaryStream& stream) {
    ConnectionEndpoint endpoint;
    endpoint.conn = &stream.connection_;
    endpoint.step = &step;
    endpoint.stream = &stream;
    return endpoint;
  }

  ConnectionEndpoint MakeEndpoint(Step& step, PressureBoundaryCondition& pressure) {
    ConnectionEndpoint endpoint;
    endpoint.conn = &pressure.connection_;
    endpoint.step = &step;
    endpoint.p_boundary = &pressure;
    return endpoint;
  }

  void MutualDisconnect(ConnectionEndpoint& a, ConnectionEndpoint& b) {
    if (a.conn) {
      a.conn->Disconnect();
    }
    if (b.conn) {
      b.conn->Disconnect();
    }
  }

  bool MakeMutualConnection(ConnectionEndpoint& a, ConnectionEndpoint& b) {
    assert(a.conn);
    assert(b.conn);
    if (!a.conn->CanConnect() || !b.conn->CanConnect()) // Can't connect to an endpoint that already has a connection 
      return false; 
    
    size_t id = GenerateID();

    if (a.stream && b.stream) {
      a.conn->MakeConnection(id, *b.stream);
      b.conn->MakeConnection(id, *a.stream);
      return true;
    }
    if (a.stream && b.p_boundary) {
      a.conn->MakeConnection(id, *b.p_boundary);
      b.conn->MakeConnection(id, *a.stream);
      return true;
    }
    if (a.p_boundary && b.stream) {
      a.conn->MakeConnection(id, *b.stream);
      b.conn->MakeConnection(id, *a.p_boundary);
      return true;
    }
    return false; // pressure-pressure not supported
  }

} // end adsim namespace


CEREAL_CLASS_VERSION(adsim::Connection, 1);
