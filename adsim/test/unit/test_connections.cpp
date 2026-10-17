// STL includes
#include <iostream>
#include <cmath>
#include <memory>
#include <set>

// ads includes
#include "adsim/objects/fluid.h"
#include "adsim/boundary_conditions/connection.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/equations/equation.h"
#include "adsim/objects/compressor.h"
#include "adsim/solver/cycle.h"

// Gtest includes
#include <gtest/gtest.h>

// CoolProp includes
#include "CoolProp.h"


namespace adsim {

  // Helper function to create a test step with flow boundary
  std::unique_ptr<Step> CreateTestStep(const std::string& name) {
    auto step = std::make_unique<Step>();
    step->step_name_ = name;
    
    // Create flow boundary with valve stream
    auto flowBC = std::make_unique<FlowBoundaryCondition>();
    flowBC->AddStream("TestStream");
    step->top_boundary_condition_ = std::move(flowBC);
    
    // Create pressure boundary
    auto pressureBC = std::make_unique<PressureBoundaryCondition>();
    step->bottom_boundary_condition_ = std::move(pressureBC);
    
    return step;
  }

  // Helper function to get boundary stream from step
  IBoundaryStream* GetBoundaryStream(Step& step) {
    auto* flowBC = dynamic_cast<FlowBoundaryCondition*>(step.top_boundary_condition_.get());
    if (flowBC && !flowBC->streams_.empty()) {
      return flowBC->streams_[0].get();
    }
    return nullptr;
  }

  // Helper function to get pressure boundary from step
  PressureBoundaryCondition* GetPressureBoundary(Step& step) {
    return dynamic_cast<PressureBoundaryCondition*>(step.bottom_boundary_condition_.get());
  }

  TEST(ConnectionTests, BasicConnectionCreation) {
    auto step1 = CreateTestStep("Step1");
    auto step2 = CreateTestStep("Step2");
    
    auto* stream1 = GetBoundaryStream(*step1);
    auto* stream2 = GetBoundaryStream(*step2);
    
    ASSERT_NE(stream1, nullptr);
    ASSERT_NE(stream2, nullptr);
    
    // Initially not connected
    EXPECT_FALSE(stream1->connection_.IsConnected());
    EXPECT_FALSE(stream2->connection_.IsConnected());
    EXPECT_TRUE(stream1->connection_.CanConnect());
    EXPECT_TRUE(stream2->connection_.CanConnect());
    
    // Create endpoints and connect
    auto endpoint1 = MakeEndpoint(*step1, *stream1);
    auto endpoint2 = MakeEndpoint(*step2, *stream2);
    
    bool success = MakeMutualConnection(endpoint1, endpoint2);
    EXPECT_TRUE(success);
    
    // Both should now be connected
    EXPECT_TRUE(stream1->connection_.IsConnected());
    EXPECT_TRUE(stream2->connection_.IsConnected());
    EXPECT_FALSE(stream1->connection_.CanConnect());
    EXPECT_FALSE(stream2->connection_.CanConnect());
    
    // Both should have the same pair ID
    EXPECT_EQ(stream1->connection_.GetID(), stream2->connection_.GetID());
    EXPECT_NE(stream1->connection_.GetID(), 0);
  }

  TEST(ConnectionTests, StreamToPressureConnection) {
    auto step1 = CreateTestStep("Step1");
    auto step2 = CreateTestStep("Step2");
    
    auto* stream1 = GetBoundaryStream(*step1);
    auto* pressure2 = GetPressureBoundary(*step2);
    
    ASSERT_NE(stream1, nullptr);
    ASSERT_NE(pressure2, nullptr);
    
    // Create endpoints and connect
    auto endpoint1 = MakeEndpoint(*step1, *stream1);
    auto endpoint2 = MakeEndpoint(*step2, *pressure2);
    
    bool success = MakeMutualConnection(endpoint1, endpoint2);
    EXPECT_TRUE(success);
    
    // Both should be connected with same pair ID
    EXPECT_TRUE(stream1->connection_.IsConnected());
    EXPECT_TRUE(pressure2->connection_.IsConnected());
    EXPECT_EQ(stream1->connection_.GetID(), pressure2->connection_.GetID());
  }

  TEST(ConnectionTests, PressureToPressureConnectionFails) {
    auto step1 = CreateTestStep("Step1");
    auto step2 = CreateTestStep("Step2");
    
    auto* pressure1 = GetPressureBoundary(*step1);
    auto* pressure2 = GetPressureBoundary(*step2);
    
    ASSERT_NE(pressure1, nullptr);
    ASSERT_NE(pressure2, nullptr);
    
    // Create endpoints and try to connect
    auto endpoint1 = MakeEndpoint(*step1, *pressure1);
    auto endpoint2 = MakeEndpoint(*step2, *pressure2);
    
    bool success = MakeMutualConnection(endpoint1, endpoint2);
    EXPECT_FALSE(success); // Should fail
    
    // Neither should be connected
    EXPECT_FALSE(pressure1->connection_.IsConnected());
    EXPECT_FALSE(pressure2->connection_.IsConnected());
  }

  TEST(ConnectionTests, AlreadyConnectedFails) {
    auto step1 = CreateTestStep("Step1");
    auto step2 = CreateTestStep("Step2");
    auto step3 = CreateTestStep("Step3");
    
    auto* stream1 = GetBoundaryStream(*step1);
    auto* stream2 = GetBoundaryStream(*step2);
    auto* stream3 = GetBoundaryStream(*step3);
    
    // Connect step1 to step2
    auto endpoint1 = MakeEndpoint(*step1, *stream1);
    auto endpoint2 = MakeEndpoint(*step2, *stream2);
    bool success = MakeMutualConnection(endpoint1, endpoint2);
    EXPECT_TRUE(success);
    
    // Try to connect step1 to step3 (should fail)
    auto endpoint3 = MakeEndpoint(*step3, *stream3);
    bool success2 = MakeMutualConnection(endpoint1, endpoint3);
    EXPECT_FALSE(success2);
    
    // step1 should still be connected to step2, not step3
    EXPECT_TRUE(stream1->connection_.IsConnected());
    EXPECT_TRUE(stream2->connection_.IsConnected());
    EXPECT_FALSE(stream3->connection_.IsConnected());
    EXPECT_EQ(stream1->connection_.GetID(), stream2->connection_.GetID());
  }

  TEST(ConnectionTests, MutualDisconnect) {
    auto step1 = CreateTestStep("Step1");
    auto step2 = CreateTestStep("Step2");
    
    auto* stream1 = GetBoundaryStream(*step1);
    auto* stream2 = GetBoundaryStream(*step2);
    
    // Connect them
    auto endpoint1 = MakeEndpoint(*step1, *stream1);
    auto endpoint2 = MakeEndpoint(*step2, *stream2);
    bool success = MakeMutualConnection(endpoint1, endpoint2);
    EXPECT_TRUE(success);
    
    // Verify connected
    EXPECT_TRUE(stream1->connection_.IsConnected());
    EXPECT_TRUE(stream2->connection_.IsConnected());
    
    // Disconnect
    MutualDisconnect(endpoint1, endpoint2);
    
    // Both should be disconnected
    EXPECT_FALSE(stream1->connection_.IsConnected());
    EXPECT_FALSE(stream2->connection_.IsConnected());
    EXPECT_TRUE(stream1->connection_.CanConnect());
    EXPECT_TRUE(stream2->connection_.CanConnect());
    EXPECT_EQ(stream1->connection_.GetID(), 0);
    EXPECT_EQ(stream2->connection_.GetID(), 0);
  }

  TEST(ConnectionTests, UniquePairIDs) {
    auto step1 = CreateTestStep("Step1");
    auto step2 = CreateTestStep("Step2");
    auto step3 = CreateTestStep("Step3");
    auto step4 = CreateTestStep("Step4");
    
    auto* stream1 = GetBoundaryStream(*step1);
    auto* stream2 = GetBoundaryStream(*step2);
    auto* stream3 = GetBoundaryStream(*step3);
    auto* stream4 = GetBoundaryStream(*step4);
    
    // Create two separate connections
    auto endpoint1 = MakeEndpoint(*step1, *stream1);
    auto endpoint2 = MakeEndpoint(*step2, *stream2);
    auto endpoint3 = MakeEndpoint(*step3, *stream3);
    auto endpoint4 = MakeEndpoint(*step4, *stream4);
    
    bool success1 = MakeMutualConnection(endpoint1, endpoint2);
    bool success2 = MakeMutualConnection(endpoint3, endpoint4);
    
    EXPECT_TRUE(success1);
    EXPECT_TRUE(success2);
    
    // Each pair should have different IDs
    EXPECT_NE(stream1->connection_.GetID(), stream3->connection_.GetID());
    EXPECT_NE(stream2->connection_.GetID(), stream4->connection_.GetID());
    
    // Within each pair, IDs should be the same
    EXPECT_EQ(stream1->connection_.GetID(), stream2->connection_.GetID());
    EXPECT_EQ(stream3->connection_.GetID(), stream4->connection_.GetID());
  }

  TEST(ConnectionTests, StressTestMultipleConnections) {
    const int num_steps = 10;
    std::vector<std::unique_ptr<Step>> steps;
    
    // Create many steps
    for (int i = 0; i < num_steps; ++i) {
      steps.push_back(CreateTestStep("Step" + std::to_string(i)));
    }
    
    // Connect every other step (0-1, 2-3, 4-5, etc.)
    for (int i = 0; i < num_steps - 1; i += 2) {
      auto* stream1 = GetBoundaryStream(*steps[i]);
      auto* stream2 = GetBoundaryStream(*steps[i + 1]);
      
      auto endpoint1 = MakeEndpoint(*steps[i], *stream1);
      auto endpoint2 = MakeEndpoint(*steps[i + 1], *stream2);
      
      bool success = MakeMutualConnection(endpoint1, endpoint2);
      EXPECT_TRUE(success);
      
      // Verify connection
      EXPECT_TRUE(stream1->connection_.IsConnected());
      EXPECT_TRUE(stream2->connection_.IsConnected());
      EXPECT_EQ(stream1->connection_.GetID(), stream2->connection_.GetID());
    }
    
    // Verify all connections have unique IDs
    std::set<size_t> used_ids;
    for (const auto& step : steps) {
      auto* stream = GetBoundaryStream(*step);
      if (stream->connection_.IsConnected()) {
        size_t id = stream->connection_.GetID();
        EXPECT_NE(id, 0);
        used_ids.insert(id);
      }
    }
    
    // Should have num_steps/2 unique connection IDs
    EXPECT_EQ(used_ids.size(), num_steps / 2);
  }

  TEST(ConnectionTests, DisconnectAndReconnect) {
    auto step1 = CreateTestStep("Step1");
    auto step2 = CreateTestStep("Step2");
    auto step3 = CreateTestStep("Step3");
    
    auto* stream1 = GetBoundaryStream(*step1);
    auto* stream2 = GetBoundaryStream(*step2);
    auto* stream3 = GetBoundaryStream(*step3);
    
    // Connect step1 to step2
    auto endpoint1 = MakeEndpoint(*step1, *stream1);
    auto endpoint2 = MakeEndpoint(*step2, *stream2);
    bool success1 = MakeMutualConnection(endpoint1, endpoint2);
    EXPECT_TRUE(success1);
    
    size_t original_id = stream1->connection_.GetID();
    
    // Disconnect
    MutualDisconnect(endpoint1, endpoint2);
    EXPECT_FALSE(stream1->connection_.IsConnected());
    EXPECT_FALSE(stream2->connection_.IsConnected());
    
    // Now connect step1 to step3
    auto endpoint3 = MakeEndpoint(*step3, *stream3);
    bool success2 = MakeMutualConnection(endpoint1, endpoint3);
    EXPECT_TRUE(success2);
    
    // step1 should be connected to step3, not step2
    EXPECT_TRUE(stream1->connection_.IsConnected());
    EXPECT_TRUE(stream3->connection_.IsConnected());
    EXPECT_FALSE(stream2->connection_.IsConnected());
    EXPECT_EQ(stream1->connection_.GetID(), stream3->connection_.GetID());
    EXPECT_NE(stream1->connection_.GetID(), original_id); // Should have new ID
  }

  TEST(ConnectionTests, CycleIntegration) {
    Cycle cycle;
    cycle.AddStep("Step1");
    cycle.AddStep("Step2");
    
    auto& step1 = cycle.GetStep("Step1");
    auto& step2 = cycle.GetStep("Step2");
    
    // Set up boundary conditions
    step1.top_boundary_condition_ = std::make_unique<FlowBoundaryCondition>();
    step1.bottom_boundary_condition_ = std::make_unique<PressureBoundaryCondition>();
    step2.top_boundary_condition_ = std::make_unique<FlowBoundaryCondition>();
    step2.bottom_boundary_condition_ = std::make_unique<PressureBoundaryCondition>();
    
    // Add streams
    auto* flow1 = dynamic_cast<FlowBoundaryCondition*>(step1.top_boundary_condition_.get());
    auto* flow2 = dynamic_cast<FlowBoundaryCondition*>(step2.top_boundary_condition_.get());
    flow1->AddStream("Stream1");
    flow2->AddStream("Stream2");
    
    // Get all endpoints
    auto endpoints = GetAvailableEndpoints(cycle.step_list_);
    EXPECT_EQ(endpoints.size(), 4); // 2 steps × 2 boundaries each
    
    // Find stream endpoints
    ConnectionEndpoint* stream1_endpoint = nullptr;
    ConnectionEndpoint* stream2_endpoint = nullptr;
    
    for (auto& ep : endpoints) {
      if (ep.stream && ep.step->step_name_ == "Step1") {
        stream1_endpoint = &ep;
      } else if (ep.stream && ep.step->step_name_ == "Step2") {
        stream2_endpoint = &ep;
      }
    }
    
    ASSERT_NE(stream1_endpoint, nullptr);
    ASSERT_NE(stream2_endpoint, nullptr);
    
    // Connect them
    bool success = MakeMutualConnection(*stream1_endpoint, *stream2_endpoint);
    EXPECT_TRUE(success);
    
    // Verify connection
    EXPECT_TRUE(stream1_endpoint->conn->IsConnected());
    EXPECT_TRUE(stream2_endpoint->conn->IsConnected());
    EXPECT_EQ(stream1_endpoint->conn->GetID(), stream2_endpoint->conn->GetID());
  }

  // --------------------------------------------------------------------------
  // E-06: connected pump streams must not rewrite the user's suction/discharge
  // expression, and a connected vacuum-pump outlet must use the bed-side
  // temperature (as the unconnected branch does).
  // --------------------------------------------------------------------------
  namespace {

    // Bed state with uniform pressure / fluid temperature over a few cells.
    Equations MakeUniformBedState(double pressure, double temperature, size_t n_cells = 4) {
      Equations state;
      state.pressure.var.assign(n_cells, pressure);
      state.fluid_temperature.var.assign(n_cells, temperature);
      return state;
    }

    // Gives `source` logged data at a constant pressure / temperature.
    void LogConstantStream(IBoundaryStream& source, double pressure, double temperature) {
      source.stream_data_.Clear();
      for (double t : {0.0, 5.0, 10.0}) {
        source.stream_data_.Log(1.0, temperature, temperature, pressure, 0.0, t, {}, {});
      }
    }

    void SetConstantPower(PumpBoundaryStream& stream, double power) {
      auto* model = dynamic_cast<ConstantPumpPowerModel*>(stream.power_.get());
      ASSERT_NE(model, nullptr);
      model->power_ = power;
    }

  } // namespace

  TEST(ConnectedPumpStreamTest, CompressorKeepsSuctionPressureExpression) {
    PumpBoundaryStream stream;
    stream.SetDirection(BoundaryDirection::INLET);
    stream.temperature_ = 300.0;
    SetConstantPower(stream, 2000.0);
    auto* comp = dynamic_cast<Compressor*>(stream.pump_.get());
    ASSERT_NE(comp, nullptr);
    comp->suction_pressure_.SetExpr("0.5 * 2e5");
    stream.EvaluateExpressions(EvalContext{});
    ASSERT_DOUBLE_EQ(comp->suction_pressure_.Get(), 1e5);

    const double P_connected = 1.3e5;
    PumpBoundaryStream source;
    LogConstantStream(source, P_connected, 310.0);
    stream.connection_.MakeConnection(1, source);

    const double P_bed = 3e5;
    Equations state = MakeUniformBedState(P_bed, 320.0);
    const double n = stream.GetMolarFlow(state, 2.0, BoundaryLocation::BOTTOM);
    stream.GetEnthalpyFlow(state, 2.0, BoundaryLocation::BOTTOM);

    // The user's expression and its value are untouched.
    EXPECT_EQ(comp->suction_pressure_.GetExpr(), "0.5 * 2e5");
    EXPECT_DOUBLE_EQ(comp->suction_pressure_.Get(), 1e5);

    // The flow uses the connected stream's pressure as the suction pressure.
    const double expected = comp->CalculateMolarFlow(P_bed, 2000.0, 310.0, P_connected);
    EXPECT_GT(expected, 0.0);
    EXPECT_NEAR(n, expected, 1e-12 * expected);
    EXPECT_GT(std::abs(n - comp->CalculateMolarFlow(P_bed, 2000.0, 310.0)), 1e-6 * expected);
  }

  TEST(ConnectedPumpStreamTest, CompressorInletUsesConnectedSuctionTemperature) {
    // The connected compressor's suction temperature is the connected stream's
    // temperature, for both the molar flow and the enthalpy flow (consistency).
    PumpBoundaryStream stream;
    stream.SetDirection(BoundaryDirection::INLET);
    stream.temperature_ = 293.0; // First-cycle assumption only: must NOT be used
    SetConstantPower(stream, 2000.0);
    auto* comp = dynamic_cast<Compressor*>(stream.pump_.get());
    ASSERT_NE(comp, nullptr);
    comp->suction_pressure_ = 1e5;

    const double P_connected = 1.2e5;
    const double T_connected = 360.0;
    PumpBoundaryStream source;
    LogConstantStream(source, P_connected, T_connected);
    stream.connection_.MakeConnection(1, source);

    const double P_bed = 3e5;
    Equations state = MakeUniformBedState(P_bed, 320.0);

    // The suction volumetric flow Q1 is fixed by power and pressures, and
    // n = Q1 * P1 / (R * T1), so n scales as 1/T1: the flow at T_connected is
    // the flow at the (unused) inlet param scaled by 293 / T_connected.
    const double n_at_param = comp->CalculateMolarFlow(P_bed, 2000.0, 293.0, P_connected);
    const double expected = n_at_param * 293.0 / T_connected;
    ASSERT_GT(expected, 0.0);
    EXPECT_NEAR(expected, comp->CalculateMolarFlow(P_bed, 2000.0, T_connected, P_connected), 1e-12 * expected);

    const double n_bottom = stream.GetMolarFlow(state, 2.0, BoundaryLocation::BOTTOM);
    const double n_top = stream.GetMolarFlow(state, 2.0, BoundaryLocation::TOP);
    EXPECT_NEAR(n_bottom, expected, 1e-10 * expected);
    EXPECT_NEAR(n_top, -expected, 1e-10 * expected);

    // Isothermal (default) compressor: enthalpy flow = n * T_connected, i.e. the
    // flow and the enthalpy are evaluated at the same suction temperature.
    const double H = stream.GetEnthalpyFlow(state, 2.0, BoundaryLocation::BOTTOM);
    EXPECT_NEAR(H, expected * T_connected, 1e-10 * expected * T_connected);
  }

  TEST(ConnectedPumpStreamTest, CompressorOutletTemperatureUsesConnectedSuctionPressure) {
    PumpBoundaryStream stream;
    stream.SetDirection(BoundaryDirection::INLET);
    SetConstantPower(stream, 2000.0);
    auto* comp = dynamic_cast<Compressor*>(stream.pump_.get());
    ASSERT_NE(comp, nullptr);
    comp->energy_eq_ = CreateCompressorEnergyModel(CompressorEnergyType::ADIABATIC);
    comp->suction_pressure_ = 1e5;

    const double P_connected = 1.5e5;
    const double T_connected = 310.0;
    PumpBoundaryStream source;
    LogConstantStream(source, P_connected, T_connected);
    stream.connection_.MakeConnection(1, source);

    const double P_bed = 3e5;
    Equations state = MakeUniformBedState(P_bed, 320.0);
    const double n = stream.GetMolarFlow(state, 2.0, BoundaryLocation::BOTTOM);
    const double H = stream.GetEnthalpyFlow(state, 2.0, BoundaryLocation::BOTTOM);
    ASSERT_GT(n, 0.0);

    // Adiabatic compression from the connected stream's (P, T) to the bed pressure.
    const double k = comp->poly_index_.Get();
    const double T_out = T_connected * std::pow(P_bed / P_connected, (k - 1.0) / k);
    EXPECT_NEAR(H / n, T_out, 1e-9 * T_out);
    EXPECT_DOUBLE_EQ(comp->suction_pressure_.Get(), 1e5);
  }

  TEST(ConnectedPumpStreamTest, VacuumPumpKeepsDischargeExpressionAndUsesBedTemperature) {
    const double P_discharge = 1e5;
    const double P_bed = 2e4;
    const double T_bed = 350.0;
    const double power = 1500.0;

    auto make_vacuum_stream = [&]() {
      auto s = std::make_unique<PumpBoundaryStream>();
      s->SetDirection(BoundaryDirection::OUTLET);
      s->pump_ = CreatePump(PumpType::VACUUM_PUMP);
      s->pump_->use_pump_ = true;
      s->temperature_ = 293.0; // Inlet param: must NOT be used for the outlet flow
      SetConstantPower(*s, power);
      return s;
    };

    // Unconnected reference with the same discharge pressure.
    auto reference = make_vacuum_stream();
    dynamic_cast<VacuumPump*>(reference->pump_.get())->discharge_pressure_ = P_discharge;

    // Connected stream whose user discharge expression differs from the connected pressure.
    auto stream = make_vacuum_stream();
    auto* vp = dynamic_cast<VacuumPump*>(stream->pump_.get());
    ASSERT_NE(vp, nullptr);
    vp->discharge_pressure_.SetExpr("1.2e5");
    stream->EvaluateExpressions(EvalContext{});

    PumpBoundaryStream source;
    LogConstantStream(source, P_discharge, 300.0);
    stream->connection_.MakeConnection(1, source);

    Equations state = MakeUniformBedState(P_bed, T_bed);
    // Hand calculation at the bed temperature and the connected discharge pressure.
    const double expected = vp->CalculateMolarFlow(P_bed, power, T_bed, P_discharge);
    ASSERT_GT(expected, 0.0);

    for (auto loc : {BoundaryLocation::TOP, BoundaryLocation::BOTTOM}) {
      const double n_conn = stream->GetMolarFlow(state, 2.0, loc);
      const double n_ref = reference->GetMolarFlow(state, 2.0, loc);
      EXPECT_NEAR(n_conn, n_ref, 1e-12 * std::abs(n_ref));
      EXPECT_NEAR(std::abs(n_conn), expected, 1e-12 * expected);
    }

    EXPECT_EQ(vp->discharge_pressure_.GetExpr(), "1.2e5");
    EXPECT_DOUBLE_EQ(vp->discharge_pressure_.Get(), 1.2e5);
  }

  // --- G-01: replacing/removing an endpoint must not orphan its partner ---

  namespace {
    // Default steps (top: pressure OUTLET, bottom: flow INLET "Stream 1") with
    // mole fractions set so Cycle::PerformChecks passes on its own.
    Cycle MakeTwoStepCycle() {
      Cycle cycle;
      cycle.AddStep("A");
      cycle.AddStep("B");
      for (auto& step : cycle.step_list_) {
        auto* flow = dynamic_cast<FlowBoundaryCondition*>(step.bottom_boundary_condition_.get());
        flow->streams_[0]->mole_fraction["N2"] = 1.0;
      }
      return cycle;
    }

    FlowBoundaryCondition& BottomFlow(Step& step) {
      return *dynamic_cast<FlowBoundaryCondition*>(step.bottom_boundary_condition_.get());
    }

    PressureBoundaryCondition& TopPressure(Step& step) {
      return *dynamic_cast<PressureBoundaryCondition*>(step.top_boundary_condition_.get());
    }

    // Connects A's bottom stream to B's top pressure boundary.
    bool ConnectAStreamToBPressure(Cycle& cycle) {
      auto& a = cycle.GetStep("A");
      auto& b = cycle.GetStep("B");
      auto ea = MakeEndpoint(a, *BottomFlow(a).streams_[0]);
      auto eb = MakeEndpoint(b, TopPressure(b));
      return MakeMutualConnection(ea, eb);
    }

    Fluid MakeN2Fluid() {
      Fluid fluid;
      fluid.components_ = {"N2"};
      return fluid;
    }
  } // anonymous namespace

  TEST(ConnectionTests, ConnectedCyclePassesChecks) {
    Cycle cycle = MakeTwoStepCycle();
    ASSERT_TRUE(ConnectAStreamToBPressure(cycle));
    std::string err;
    EXPECT_TRUE(cycle.PerformChecks(MakeN2Fluid(), err)) << err;
    EXPECT_TRUE(BottomFlow(cycle.GetStep("A")).streams_[0]->connection_.IsConnected());
  }

  TEST(ConnectionTests, RemoveStepDisconnectsPartner) {
    Cycle cycle = MakeTwoStepCycle();
    ASSERT_TRUE(ConnectAStreamToBPressure(cycle));

    cycle.RemoveStep("B");

    auto& conn = BottomFlow(cycle.GetStep("A")).streams_[0]->connection_;
    EXPECT_FALSE(conn.IsConnected());
    EXPECT_EQ(conn.stream_, nullptr);
    EXPECT_EQ(conn.GetID(), 0u);
  }

  TEST(ConnectionTests, ReplaceBoundaryConditionDisconnectsPartner) {
    Cycle cycle = MakeTwoStepCycle();
    ASSERT_TRUE(ConnectAStreamToBPressure(cycle));

    auto& b = cycle.GetStep("B");
    cycle.ReplaceBoundaryCondition(b.top_boundary_condition_,
                                   CreateBoundaryCondition(BoundaryType::CLOSED));
    b.top_boundary_condition_->location_ = BoundaryLocation::TOP;

    auto& conn = BottomFlow(cycle.GetStep("A")).streams_[0]->connection_;
    EXPECT_FALSE(conn.IsConnected());
    EXPECT_EQ(conn.stream_, nullptr);

    std::string err;
    EXPECT_TRUE(cycle.PerformChecks(MakeN2Fluid(), err)) << err;
  }

  TEST(ConnectionTests, ReplaceFlowBoundaryDisconnectsEveryStreamPartner) {
    Cycle cycle = MakeTwoStepCycle();
    cycle.AddStep("C");
    auto& a = cycle.GetStep("A");
    auto& b = cycle.GetStep("B");
    auto& c = cycle.GetStep("C");
    BottomFlow(a).AddStream("Stream 2");

    auto ea1 = MakeEndpoint(a, *BottomFlow(a).streams_[0]);
    auto ea2 = MakeEndpoint(a, *BottomFlow(a).streams_[1]);
    auto eb = MakeEndpoint(b, TopPressure(b));
    auto ec = MakeEndpoint(c, TopPressure(c));
    ASSERT_TRUE(MakeMutualConnection(ea1, eb));
    ASSERT_TRUE(MakeMutualConnection(ea2, ec));

    cycle.ReplaceBoundaryCondition(a.bottom_boundary_condition_,
                                   CreateBoundaryCondition(BoundaryType::PRESSURE));

    EXPECT_FALSE(TopPressure(b).connection_.IsConnected());
    EXPECT_FALSE(TopPressure(c).connection_.IsConnected());
    EXPECT_EQ(TopPressure(b).connection_.stream_, nullptr);
    EXPECT_EQ(TopPressure(c).connection_.stream_, nullptr);
  }

  TEST(ConnectionTests, ReplaceStreamDisconnectsPartnerOnly) {
    Cycle cycle = MakeTwoStepCycle();
    cycle.AddStep("C");
    auto& a = cycle.GetStep("A");
    auto& b = cycle.GetStep("B");
    auto& c = cycle.GetStep("C");

    // A.stream <-> B.stream (the pair being broken) and C.stream <-> A.pressure
    // (an unrelated pair that must survive).
    auto ea = MakeEndpoint(a, *BottomFlow(a).streams_[0]);
    auto eb = MakeEndpoint(b, *BottomFlow(b).streams_[0]);
    auto ec = MakeEndpoint(c, *BottomFlow(c).streams_[0]);
    auto eap = MakeEndpoint(a, TopPressure(a));
    ASSERT_TRUE(MakeMutualConnection(ea, eb));
    ASSERT_TRUE(MakeMutualConnection(ec, eap));

    auto& slot = BottomFlow(b).streams_[0];
    cycle.ReplaceStream(slot, CreateStream(BoundaryStreamType::VALVE));

    EXPECT_EQ(slot->GetType(), BoundaryStreamType::VALVE);
    EXPECT_FALSE(slot->connection_.IsConnected());
    EXPECT_FALSE(BottomFlow(a).streams_[0]->connection_.IsConnected());
    EXPECT_EQ(BottomFlow(a).streams_[0]->connection_.stream_, nullptr);

    EXPECT_TRUE(BottomFlow(c).streams_[0]->connection_.IsConnected());
    EXPECT_TRUE(TopPressure(a).connection_.IsConnected());
    EXPECT_EQ(BottomFlow(c).streams_[0]->connection_.stream_, &TopPressure(a).stream_data_);
  }

  TEST(ConnectionTests, RemoveStreamByIndexDisconnectsPartner) {
    Cycle cycle = MakeTwoStepCycle();
    ASSERT_TRUE(ConnectAStreamToBPressure(cycle));
    auto& flow = BottomFlow(cycle.GetStep("A"));
    flow.AddStream("Stream 2");

    cycle.RemoveStream(flow, 0);

    ASSERT_EQ(flow.streams_.size(), 1u);
    EXPECT_EQ(flow.streams_[0]->name_, "Stream 2");
    auto& conn = TopPressure(cycle.GetStep("B")).connection_;
    EXPECT_FALSE(conn.IsConnected());
    EXPECT_EQ(conn.stream_, nullptr);

    cycle.RemoveStream(flow, 5); // out of range: no-op
    EXPECT_EQ(flow.streams_.size(), 1u);
  }

  TEST(ConnectionTests, PerformChecksRejectsOrphanedConnection) {
    Cycle cycle = MakeTwoStepCycle();
    ASSERT_TRUE(ConnectAStreamToBPressure(cycle));
    auto& b = cycle.GetStep("B");

    // Replace B's pressure boundary WITHOUT disconnecting (the pre-fix GUI
    // path). Keep the old object alive so the test never touches freed memory.
    auto old_pressure = b.top_boundary_condition_;
    auto fresh = std::make_shared<PressureBoundaryCondition>();
    fresh->location_ = BoundaryLocation::TOP;
    fresh->direction_ = BoundaryDirection::OUTLET;
    b.top_boundary_condition_ = fresh;

    std::string err;
    EXPECT_FALSE(cycle.PerformChecks(MakeN2Fluid(), err));
    EXPECT_NE(err.find("no valid partner"), std::string::npos) << err;
    EXPECT_NE(err.find("step 'A' endpoint 'Stream 1'"), std::string::npos) << err;
  }

  TEST(ConnectionTests, PerformChecksRejectsPartnerWithStalePointer) {
    Cycle cycle = MakeTwoStepCycle();
    ASSERT_TRUE(ConnectAStreamToBPressure(cycle));
    auto& b = cycle.GetStep("B");

    // A copy of the pressure boundary carries the same pair id, so the pair is
    // still found by id, but A's stream_ points into the old object's data.
    auto old_pressure = std::static_pointer_cast<PressureBoundaryCondition>(b.top_boundary_condition_);
    auto copy = std::make_shared<PressureBoundaryCondition>();
    copy->location_ = BoundaryLocation::TOP;
    copy->direction_ = BoundaryDirection::OUTLET;
    copy->connection_ = old_pressure->connection_;
    b.top_boundary_condition_ = copy;

    std::string err;
    EXPECT_FALSE(cycle.PerformChecks(MakeN2Fluid(), err));
    EXPECT_NE(err.find("no valid partner"), std::string::npos) << err;
  }

} // end adsim namespace
