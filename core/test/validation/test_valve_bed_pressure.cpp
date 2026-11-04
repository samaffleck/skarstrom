/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Where the linear valve law reads the bed: the bed-end face (default) or the
// end cell (IDASettings::valve_bed_pressure_ = EndCell, as toPSAil, which
// attaches its valves to the last well-mixed cell of the column).
//
// Steady flow through an inert, isothermal bed. A set molar flow F enters at
// the bottom and leaves at the top through a linear valve,
//
//     F = Cv (P_bed - P_down) / 1e5,
//
// so at steady state the bed-side pressure the valve reads is exactly
//
//     P_bed = P_down + 1e5 F / Cv,
//
// whatever the bed's flow resistance. With the face setting P_bed is the
// bed-end face pressure (the logged stream pressure); with the end-cell
// setting it is the pressure of the last cell. The two differ by the Darcy
// drop over half a cell, which the test makes large (fine particles, coarse
// grid) so that reading the wrong one fails by orders of magnitude more than
// the tolerance.
//
// Connected valves (two steps of the single-bed model, as in pressure
// equalisation): with the end-cell setting each side reads the other's end
// cell, so at the cyclic steady state the donor's outflow obeys
//
//     F = Cv (P_end,donor - P_end,receiver) / 1e5
//
// with both end-cell pressures from the equation log, while the face
// pressures are a full half-cell drop apart on each side.
//
// Neglected terms: none in the steady balance. The through-flow bed starts at
// the steady pressure; the connected steps start from each other's end state
// and relax with a time constant of about 4 s (measured), so after an 80 s
// step the residual is exp(-20) of the 1 bar initial offset. The replayed
// partner pressure converges geometrically from cycle to cycle (factor ~0.33
// measured); after 18 cycles it moves by less than 0.01 Pa per cycle.

#include "validation_helpers.h"

namespace adsim {

  using namespace adsim::validation;

  namespace {

    constexpr double kT = 298.0;
    constexpr double kEps = 0.4;
    constexpr double kL = 1.0, kD = 0.1;
    constexpr double kF = 0.05;        // [mol/s]
    constexpr double kCv = 0.05;       // [mol/s per bar]: 1 bar across the valve
    constexpr double kPDown = 1.0e5;   // [Pa]
    constexpr int kCells = 10;

    void BuildInertBed(UserData& data, ValveBedPressure where) {
      data.bed.wall_.thickness_ = 0.005;
      data.fluid.AddComponent("Nitrogen");
      data.fluid.AddComponent("Oxygen");
      auto ads = data.adsorbent_directory.Add("Packing");
      auto lid = data.bed.layers_.Add("Packing");
      auto& layer = data.bed.layers_.Get(lid);
      layer.momentum_model_ = MomentumModel::KarmenCozeny;
      layer.domain_.length = kL;
      layer.domain_.n_cells = kCells;
      layer.eb_ = kEps;
      layer.particle_diameter_ = 2.0e-4;  // fine: a half-cell drop of ~2 kPa
      layer.composition_[ads] = 1.0;
      std::get<CylinderBed>(data.bed.shape_).diameter = kD;
      auto& adsorbent = data.adsorbent_directory.Get(ads);
      adsorbent.ep_ = 0.0;
      adsorbent.solid_density_ = 1000.0;
      adsorbent.heat_capacity_ = 900.0;
      data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
      data.ic.pressure = kPDown + 1e5 * kF / kCv;
      data.ic.temperature = kT;
      data.ic.wall_temperature = kT;
      data.ic.velocity = 0.0;
      data.ic.mole_fraction["Nitrogen"] = 0.79;
      data.ic.mole_fraction["Oxygen"] = 0.21;
      data.settings.reltol = 1e-9;
      data.settings.max_step = 1.0;
      data.settings.valve_bed_pressure_ = where;
    }

    template <class Map>
    void SetAir(Map& y) {
      y["Nitrogen"] = 0.79;
      y["Oxygen"] = 0.21;
    }

    std::unique_ptr<ConstantBoundaryValue> Constant(double v) {
      auto bv = std::make_unique<ConstantBoundaryValue>();
      bv->constant_value_ = v;
      return bv;
    }

    // Bottom inlet with a set molar flow of air.
    void FeedBottom(Step& step, double F) {
      auto stream = std::make_unique<BoundaryStream>();
      auto flow = std::make_unique<MolarFlowBoundaryValue>();
      static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = F;
      stream->boundary_flow_ = std::move(flow);
      stream->temperature_ = kT;
      SetAir(stream->mole_fraction);
      auto inlet = std::make_unique<FlowBoundaryCondition>();
      inlet->direction_ = BoundaryDirection::INLET;
      inlet->location_ = BoundaryLocation::BOTTOM;
      inlet->streams_.push_back(std::move(stream));
      step.bottom_boundary_condition_ = std::move(inlet);
    }

    // Top valve; returns the stream so the test can read its log.
    ValveBoundaryStream* TopValve(Step& step, BoundaryDirection dir, double P_other) {
      auto valve = std::make_unique<ValveBoundaryStream>();
      valve->boundary_value_ = Constant(P_other);
      valve->valve_.cv_ = kCv;
      valve->temperature_ = kT;
      SetAir(valve->mole_fraction);
      auto* raw = valve.get();
      auto bc = std::make_unique<FlowBoundaryCondition>();
      bc->direction_ = dir;
      bc->location_ = BoundaryLocation::TOP;
      bc->streams_.push_back(std::move(valve));
      step.top_boundary_condition_ = std::move(bc);
      return raw;
    }

    // Pressure of the top end cell (the last interior cell of the log) at the
    // end of a step.
    double TopEndCellPressure(const EquationLog& log, int cycle, int step) {
      const auto P = log.GetVariable(cycle, step, VariableType::PRESSURE);
      const auto& last = P.values.back();
      return last[last.size() - 2];
    }

    struct ThroughFlow { double P_face, P_end, F_out; };

    ThroughFlow RunThroughFlow(ValveBedPressure where) {
      UserData data;
      BuildInertBed(data, where);
      data.cycle.AddStep("Through");
      Step& step = data.cycle.GetStep("Through");
      step.step_time_ = 40.0;
      FeedBottom(step, kF);
      auto* valve = TopValve(step, BoundaryDirection::OUTLET, kPDown);
      data.cycle.number_of_cycles_ = 1;

      CoupledSolver solver{data};
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
      const auto& sd = valve->stream_data_;
      return {sd.pressure_boundary_.back(), TopEndCellPressure(solver.GetLog(), 0, 0),
              sd.GetOutflow(sd.time_series_.back())};
    }

    struct ConnectedFlow {
      double F_donor, F_receiver;
      double P_end_donor, P_end_receiver;
      double P_face_donor, P_face_receiver;
    };

    ConnectedFlow RunConnected(ValveBedPressure where) {
      UserData data;
      BuildInertBed(data, where);
      const double t_step = 80.0;
      data.cycle.AddStep("Donor");
      Step& donor = data.cycle.GetStep("Donor");
      donor.step_time_ = t_step;
      FeedBottom(donor, kF);
      auto* out = TopValve(donor, BoundaryDirection::OUTLET, kPDown);  // first-cycle assumption

      data.cycle.AddStep("Receiver");
      Step& receiver = data.cycle.GetStep("Receiver");
      receiver.step_time_ = t_step;
      auto* in = TopValve(receiver, BoundaryDirection::INLET, kPDown + 1e5 * kF / kCv);
      {
        auto bc = std::make_unique<PressureBoundaryCondition>();
        bc->direction_ = BoundaryDirection::OUTLET;
        bc->location_ = BoundaryLocation::BOTTOM;
        bc->pressure_ = Constant(kPDown);
        SetAir(bc->mole_fraction_);
        receiver.bottom_boundary_condition_ = std::move(bc);
      }
      {
        auto a = MakeEndpoint(data.cycle.GetStep("Donor"), *out);
        auto b = MakeEndpoint(data.cycle.GetStep("Receiver"), *in);
        EXPECT_TRUE(MakeMutualConnection(a, b));
      }
      data.cycle.number_of_cycles_ = 18;

      CoupledSolver solver{data};
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
      const auto& log = solver.GetLog();
      const int last = *log.sys_data.cycle_index_set_.rbegin();
      // GUARD: the replayed coupling has converged (the previous cycle ended
      // the donor step at the same pressure).
      EXPECT_NEAR(TopEndCellPressure(log, last - 1, 0), TopEndCellPressure(log, last, 0), 0.01)
          << "Not at the cyclic steady state";
      const auto& so = out->stream_data_;
      const auto& si = in->stream_data_;
      return {so.GetOutflow(so.time_series_.back()), -si.GetOutflow(si.time_series_.back()),
              TopEndCellPressure(log, last, 0), TopEndCellPressure(log, last, 1),
              so.pressure_boundary_.back(), si.pressure_boundary_.back()};
    }

  } // namespace

  TEST(ValveBedPressure, SteadyThroughFlowReadsTheChosenPressure) {
    const double P_valve = kPDown + 1e5 * kF / kCv;  // the bed side of the valve law
    const auto face = RunThroughFlow(ValveBedPressure::Face);
    const auto cell = RunThroughFlow(ValveBedPressure::EndCell);
    std::cout << "[ValveBedPressure] face run: P_face " << face.P_face << " P_end " << face.P_end
              << " F " << face.F_out << "; end-cell run: P_face " << cell.P_face << " P_end "
              << cell.P_end << " F " << cell.F_out << " (valve side " << P_valve << " Pa)\n";

    // Steady state: what enters leaves.
    EXPECT_NEAR(face.F_out, kF, 1e-6 * kF);
    EXPECT_NEAR(cell.F_out, kF, 1e-6 * kF);
    // GUARD: the end cell sits a large half-cell drop above the face, so the
    // assertions below cannot pass with the wrong pressure.
    ASSERT_GT(face.P_end - face.P_face, 500.0);
    ASSERT_GT(cell.P_end - cell.P_face, 500.0);
    // The valve reads the face (default) or the end cell (option).
    EXPECT_NEAR(face.P_face, P_valve, 1.0);
    EXPECT_NEAR(cell.P_end, P_valve, 1.0);
  }

  TEST(ValveBedPressure, ConnectedValvesReadBothEndCells) {
    const double dP_valve = 1e5 * kF / kCv;
    const auto face = RunConnected(ValveBedPressure::Face);
    const auto cell = RunConnected(ValveBedPressure::EndCell);
    for (const auto* r : {&face, &cell}) {
      std::cout << "[ValveBedPressure] connected (" << (r == &face ? "face" : "end cell")
                << "): F donor " << r->F_donor << " receiver " << r->F_receiver
                << "; end cells " << r->P_end_donor << " / " << r->P_end_receiver
                << "; faces " << r->P_face_donor << " / " << r->P_face_receiver << "\n";
      // The connection conserves moles and both steps are steady.
      EXPECT_NEAR(r->F_donor, kF, 1e-6 * kF);
      EXPECT_NEAR(r->F_receiver, kF, 1e-6 * kF);
    }
    // GUARD: the end-cell and face differences across the valve are far apart
    // (a half-cell drop on each side).
    ASSERT_GT((cell.P_end_donor - cell.P_end_receiver) - (cell.P_face_donor - cell.P_face_receiver), 1000.0);
    // Face setting: the faces are 1e5 F / Cv apart; end-cell setting: the end cells.
    EXPECT_NEAR(face.P_face_donor - face.P_face_receiver, dP_valve, 1.0);
    EXPECT_NEAR(cell.P_end_donor - cell.P_end_receiver, dP_valve, 1.0);
  }

} // end adsim namespace
