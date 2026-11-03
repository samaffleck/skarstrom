// test_bed_series_solver.cpp
//
// Simulation::Series and Simulation::Profile on real runs:
//   * Top / Bottom are the physical boundary values: the pressure, temperature
//     and mole fraction of the bed-end endpoint's log (StreamSeries) at inlets
//     and outlets, for flow and pressure boundaries, and the face value at a
//     closed end; never the ghost cell.
//   * Profiles are at the cell centres (Positions) without the ghost cells,
//     or at the faces (FacePositions) for face variables.
//   * The time base is the logged time from the start of the cycle.

#include <cmath>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "../validation/breakthrough_case.h"
#include "adsim/simulation.h"
#include "adsim/solver/bed_series.h"

namespace adsim {

namespace {

using validation::BreakthroughCase;

const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";

constexpr double kFeedTemperature = 313.15;
constexpr double kBedTemperature = 298.0;

// The plain breakthrough case (validation/breakthrough_case.h), adiabatic and
// adsorbing, with the feed 15 K warmer than the bed so the inlet ghost cell
// (2 T_feed - T_1) is far from the feed temperature.
Simulation Breakthrough(bool feed_from_top) {
  BreakthroughCase c;
  c.n_cells = 30;
  c.run_time = 150.0;
  c.max_step = 1.0;
  c.reltol = 1e-6;
  c.adsorbing = true;
  c.y_feed = 0.2;
  c.feed_temperature = kFeedTemperature;
  c.bed_temperature = kBedTemperature;
  c.feed_from_top = feed_from_top;
  UserData data;
  validation::Build(data, c);
  Simulation sim;
  sim.SetData(std::move(data));
  return sim;
}

// The log name of the one endpoint at this end of the step, if it has one.
std::optional<std::string> EndpointAt(const Simulation& sim, int step, BoundaryLocation loc) {
  std::optional<std::string> name;
  int count = 0;
  for (const auto& s : sim.Streams())
    if (s.step == step && s.location == loc) {
      name = s.name;
      ++count;
    }
  return count == 1 ? name : std::nullopt;
}

// The raw log rows of one variable over one step of the last cycle.
VariableLogEntry Raw(const Simulation& sim, int step, const BedVariable& var) {
  return GetBedVariable(sim.Log(), sim.LastCycle(), step, var);
}

// IBoundaryCondition::GetType is non-const.
IBoundaryCondition& BC(Simulation& sim, int step, BoundaryLocation loc) {
  Step& s = sim.Data().cycle.step_list_.at(static_cast<size_t>(step));
  return loc == BoundaryLocation::TOP ? *s.top_boundary_condition_ : *s.bottom_boundary_condition_;
}

SeriesLocation Where(BoundaryLocation loc) {
  return loc == BoundaryLocation::TOP ? SeriesLocation::Top : SeriesLocation::Bottom;
}

// Checks, step by step over the last cycle, that the bed-end series of y, T
// and P equal the endpoint logs where the end has one endpoint, and the face
// values at a closed end. Counts what it saw.
struct Coverage {
  int flow_inlets = 0, flow_outlets = 0, pressure_inlets = 0, pressure_outlets = 0, closed = 0;
};

Coverage CheckBedEndsAgainstLogs(Simulation& sim, const std::string& component) {
  Coverage seen;
  const int steps = static_cast<int>(sim.Data().cycle.step_list_.size());
  for (BoundaryLocation loc : {BoundaryLocation::BOTTOM, BoundaryLocation::TOP}) {
    const TimeSeries P = sim.Series("pressure", "", Where(loc));
    const TimeSeries T = sim.Series("temperature", "", Where(loc));
    const TimeSeries y = sim.Series("mole_fraction", component, Where(loc));
    size_t offset = 0;
    for (int s = 0; s < steps; ++s) {
      const VariableLogEntry rawP = Raw(sim, s, BedVariable::Bulk(VariableType::PRESSURE));
      const VariableLogEntry rawT = Raw(sim, s, BedVariable::Bulk(VariableType::TEMPERATURE));
      const VariableLogEntry rawY = Raw(sim, s, BedVariable::Component(ComponentVariableType::MOLE_FRACTION, component));
      const size_t n = rawP.times.size();
      IBoundaryCondition& bc = BC(sim, s, loc);
      const std::string where = "step " + std::to_string(s) + (loc == BoundaryLocation::TOP ? " top" : " bottom");
      if (bc.GetType() == BoundaryType::CLOSED) {
        ++seen.closed;
        for (size_t k = 0; k < n; ++k) {
          const auto& p = rawP.values[k];
          const size_t g = loc == BoundaryLocation::TOP ? p.size() - 1 : 0;
          const size_t e = loc == BoundaryLocation::TOP ? p.size() - 2 : 1;
          EXPECT_DOUBLE_EQ(P.values[offset + k], 0.5 * (p[g] + p[e])) << where;
          EXPECT_DOUBLE_EQ(T.values[offset + k], 0.5 * (rawT.values[k][g] + rawT.values[k][e])) << where;
          EXPECT_DOUBLE_EQ(y.values[offset + k], rawY.values[k][g]) << where;
        }
      } else {
        if (const auto name = EndpointAt(sim, s, loc)) {
          const TimeSeries sp = sim.StreamSeries(s, *name, "pressure");
          const TimeSeries st = sim.StreamSeries(s, *name, "temperature");
          const TimeSeries sy = sim.StreamSeries(s, *name, "mole_fraction", component);
          EXPECT_EQ(sp.unit, "Pa");
          EXPECT_EQ(sp.time.size(), n) << where;
          for (size_t k = 0; k < n && k < sp.values.size(); ++k) {
            EXPECT_DOUBLE_EQ(P.values[offset + k], sp.values[k]) << where << ", row " << k;
            EXPECT_DOUBLE_EQ(T.values[offset + k], st.values[k]) << where << ", row " << k;
            EXPECT_DOUBLE_EQ(y.values[offset + k], sy.values[k]) << where << ", row " << k;
          }
          const bool inlet = bc.direction_ == BoundaryDirection::INLET;
          if (bc.GetType() == BoundaryType::PRESSURE) ++(inlet ? seen.pressure_inlets : seen.pressure_outlets);
          else ++(inlet ? seen.flow_inlets : seen.flow_outlets);
        }
      }
      offset += n;
    }
    EXPECT_EQ(offset, P.values.size());
  }
  return seen;
}

}  // namespace

TEST(BedSeriesSolver, BreakthroughBedEndsAreTheBoundaryValues) {
  for (bool from_top : {false, true}) {
    Simulation sim = Breakthrough(from_top);
    const RunResult r = sim.Run();
    ASSERT_TRUE(r.Ok()) << r.message;
    const SeriesLocation in = from_top ? SeriesLocation::Top : SeriesLocation::Bottom;
    const SeriesLocation out = from_top ? SeriesLocation::Bottom : SeriesLocation::Top;

    // Inlet: the feed temperature and composition, not the ghost cell.
    const TimeSeries T_in = sim.Series("temperature", "", in);
    const TimeSeries y_in = sim.Series("mole_fraction", "Oxygen", in);
    for (double v : T_in.values) EXPECT_NEAR(v, kFeedTemperature, 1e-9);
    for (double v : y_in.values) EXPECT_NEAR(v, 0.2, 1e-12);
    // The inlet ghost really is elsewhere (2 T_feed - T_1, with T_1 still
    // near the bed temperature early on): the old reading.
    const VariableLogEntry rawT = Raw(sim, 0, BedVariable::Bulk(VariableType::TEMPERATURE));
    const double ghost = from_top ? rawT.values[1].back() : rawT.values[1].front();
    EXPECT_GT(std::abs(ghost - kFeedTemperature), 5.0);

    // Outlet: the set pressure of the pressure boundary.
    const TimeSeries P_out = sim.Series("pressure", "", out);
    for (double v : P_out.values) EXPECT_NEAR(v, validation::kCasePressure, 1e-5 * validation::kCasePressure);

    const Coverage seen = CheckBedEndsAgainstLogs(sim, "Oxygen");
    EXPECT_EQ(seen.flow_inlets, 1);
    EXPECT_EQ(seen.pressure_outlets, 1);
  }
}

TEST(BedSeriesSolver, CycleBedEndsMatchTheEndpointLogs) {
  Simulation sim = Simulation::FromFile((kExamples / "Perez2019_A.json").string());
  sim.SetParameter("cycle.number_of_cycles", 2);
  const RunResult r = sim.Run();
  ASSERT_TRUE(r.Ok()) << r.message;
  const Coverage seen = CheckBedEndsAgainstLogs(sim, "CarbonDioxide");
  EXPECT_GT(seen.flow_inlets + seen.pressure_inlets, 0);
  EXPECT_GT(seen.flow_outlets + seen.pressure_outlets, 0);
  EXPECT_GT(seen.pressure_inlets + seen.pressure_outlets, 0);
  EXPECT_GT(seen.closed, 0);
}

TEST(BedSeriesSolver, ProfilesAreAtCellCentresOrFaces) {
  Simulation sim = Breakthrough(false);
  ASSERT_TRUE(sim.Run().Ok());
  const std::vector<double> faces = sim.FacePositions();
  const std::vector<double> centres = sim.Positions();
  ASSERT_EQ(faces.size(), 31u);    // 30 cells: faces 0..L
  ASSERT_EQ(centres.size(), 30u);  // one centre per bed cell, no ghosts
  EXPECT_DOUBLE_EQ(faces.front(), 0.0);
  EXPECT_NEAR(faces.back(), 1.0, 1e-12);
  for (size_t i = 0; i < centres.size(); ++i) EXPECT_NEAR(centres[i], (i + 0.5) / 30.0, 1e-12);

  const AxialProfile T = sim.Profile("temperature");
  EXPECT_FALSE(T.faces);
  EXPECT_EQ(T.z, centres);
  const VariableLogEntry rawT = Raw(sim, 0, BedVariable::Bulk(VariableType::TEMPERATURE));
  ASSERT_EQ(T.rows, rawT.values.size());
  ASSERT_EQ(T.cols, rawT.values[0].size() - 2);
  for (size_t t = 0; t < T.rows; ++t)
    for (size_t i = 0; i < T.cols; ++i) EXPECT_EQ(T.values[t * T.cols + i], rawT.values[t][i + 1]);

  const AxialProfile u = sim.Profile("velocity");
  EXPECT_TRUE(u.faces);
  EXPECT_EQ(u.z, faces);
  const VariableLogEntry rawU = Raw(sim, 0, BedVariable::Bulk(VariableType::VELOCITY));
  for (size_t t = 0; t < u.rows; ++t)
    for (size_t i = 0; i < u.cols; ++i) EXPECT_EQ(u.values[t * u.cols + i], rawU.values[t][i]);

  // A z query at a cell centre is that cell; at a face, that face.
  const TimeSeries T5 = sim.Series("temperature", "", SeriesLocation::Position, centres[5]);
  const TimeSeries u0 = sim.Series("velocity", "", SeriesLocation::Bottom);
  const TimeSeries uL = sim.Series("velocity", "", SeriesLocation::Top);
  for (size_t t = 0; t < rawT.values.size(); ++t) {
    EXPECT_NEAR(T5.values[t], rawT.values[t][6], 1e-9 * rawT.values[t][6]);
    EXPECT_EQ(u0.values[t], rawU.values[t].front());
    EXPECT_EQ(uL.values[t], rawU.values[t].back());
  }
}

TEST(BedSeriesSolver, TimeBaseIsUnchanged) {
  Simulation sim = Simulation::FromFile((kExamples / "Perez2019_A.json").string());
  sim.SetParameter("cycle.number_of_cycles", 2);
  ASSERT_TRUE(sim.Run().Ok());
  for (int cycle : {0, 1}) {
    // The logged times of the cycle, from its first logged time.
    std::vector<double> expected;
    const auto& cycle_log = sim.Log().GetSimulationLog().at(cycle);
    const double start = cycle_log.begin()->second.step_data.times.front();
    for (const auto& [step, log] : cycle_log)
      for (double t : log.step_data.times) expected.push_back(t - start);
    for (SeriesLocation where : {SeriesLocation::Top, SeriesLocation::Bottom, SeriesLocation::Position}) {
      EXPECT_EQ(sim.Series("pressure", "", where, 0.3, cycle).time, expected);
      EXPECT_EQ(sim.Series("mole_fraction", "CarbonDioxide", where, 0.3, cycle).time, expected);
      EXPECT_EQ(sim.Series("velocity", "", where, 0.3, cycle).time, expected);
    }
  }
}

}  // namespace adsim
