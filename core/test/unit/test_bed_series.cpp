// test_bed_series.cpp
//
// Bed variables read back from the equation log at the right positions
// (adsim/solver/bed_series.h). A hand-built log on a non-uniform mesh holds a
// linear cell profile with garbage in the ghost cells and a linear face
// profile: every query must return the line at the queried position, never a
// ghost value, and the bed ends must be the face values.

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/equations/equation.h"
#include "adsim/solver/bed_series.h"
#include "adsim/solver/parameterisation.h"

namespace adsim {

namespace {

constexpr double kGhost = 1.0e6;  // never a bed value

// Ghost, three cells of 0.1 m, two of 0.25 m, ghost: N = 7 cells, faces at
// 0, 0.1, 0.2, 0.3, 0.55, 0.8 (positions[i] is the face between cells i and
// i+1), cell centres 0.05, 0.15, 0.25, 0.425, 0.675.
DomainVectors Mesh() {
  DomainVectors d;
  d.cell_width = {0.1, 0.1, 0.1, 0.1, 0.25, 0.25, 0.25};
  const size_t N = d.cell_width.size();
  d.positions.assign(N - 1, 0.0);
  for (size_t i = 1; i + 1 < N; ++i) d.positions[i] = d.positions[i - 1] + d.cell_width[i];
  d.bed_cell_indices = {1, 2, 3, 4, 5};
  return d;
}

double CellLine(double z, double t) { return 2.0 + 3.0 * z + t; }
double FaceLine(double z, double t) { return -1.0 + 5.0 * z - 0.5 * t; }
double FractionLine(double z) { return 0.1 + 0.5 * z; }

// The time is absolute solver time: the cycle starts at t = 100 s.
EquationLog LinearLog() {
  const DomainVectors mesh = Mesh();
  const std::vector<double> centres = CellCentrePositions(mesh);

  Equations state;
  state.domain_vectors = mesh;
  state.mole_fraction.emplace("A", Equation("Mole Fraction", CellType::CellCentre, EquationType::Differential, {}));

  EquationLog log;
  log.sys_data.step_list_ = {"Feed", "Blowdown"};
  const std::vector<std::pair<int, std::vector<double>>> steps = {{0, {100.0, 101.0, 102.5}}, {1, {102.5, 104.0}}};
  for (const auto& [step, times] : steps) {
    for (double t : times) {
      std::vector<double> cells(mesh.cell_width.size(), kGhost);
      std::vector<double> fractions(mesh.cell_width.size(), kGhost);
      for (size_t i = 0; i < centres.size(); ++i) {
        cells[i + 1] = CellLine(centres[i], t);
        fractions[i + 1] = FractionLine(centres[i]);
      }
      // Ghosts as the boundary conditions set them: the temperature mirror
      // 2 T_face - T_1 at the bottom, a copy of the end cell at the top.
      cells.front() = 2.0 * CellLine(0.0, t) - cells[1];
      cells.back() = cells[cells.size() - 2];
      fractions.front() = FractionLine(0.0);  // a mole-fraction ghost holds the face value
      fractions.back() = fractions[fractions.size() - 2];
      std::vector<double> faces(mesh.positions.size());
      for (size_t i = 0; i < faces.size(); ++i) faces[i] = FaceLine(mesh.positions[i], t);
      state.fluid_temperature.var = cells;
      state.velocity.var = faces;
      state.mole_fraction.at("A").var = fractions;
      log.LogState(state, t, step, 0);
    }
  }
  return log;
}

std::vector<double> AllTimes() { return {100.0, 101.0, 102.5, 102.5, 104.0}; }

}  // namespace

TEST(BedSeries, CellCentresAndFaceRows) {
  const DomainVectors mesh = Mesh();
  const std::vector<double> centres = CellCentrePositions(mesh);
  const std::vector<double> expected = {0.05, 0.15, 0.25, 0.425, 0.675};
  ASSERT_EQ(centres.size(), expected.size());
  for (size_t i = 0; i < centres.size(); ++i) EXPECT_NEAR(centres[i], expected[i], 1e-12);
  EXPECT_TRUE(IsFaceRow(mesh, 6));
  EXPECT_FALSE(IsFaceRow(mesh, 7));
}

TEST(BedSeries, LinearCellProfileReadsBackAtItsPositions) {
  const EquationLog log = LinearLog();
  const std::vector<double> t_abs = AllTimes();
  // Interior points (between two cell centres), the cell centres themselves,
  // and points in the half cells at the two ends.
  for (double z : {0.0, 0.02, 0.05, 0.12, 0.3, 0.425, 0.5, 0.675, 0.75, 0.8}) {
    std::vector<double> t, v;
    ASSERT_TRUE(ExtractBulkTimeSeries(log, VariableType::TEMPERATURE, ParamSpatialLoc::Position, z, t, v));
    ASSERT_EQ(v.size(), t_abs.size());
    for (size_t k = 0; k < v.size(); ++k) {
      if (z <= 0.675) {
        EXPECT_NEAR(v[k], CellLine(z, t_abs[k]), 1e-9) << "z = " << z << ", row " << k;
      } else {
        // Above the last centre the top face value (the end cell, zero
        // gradient) holds: the line is flat there.
        EXPECT_NEAR(v[k], CellLine(0.675, t_abs[k]), 1e-9) << "z = " << z << ", row " << k;
      }
    }
  }
}

TEST(BedSeries, BedEndsAreFaceValuesNotGhosts) {
  const EquationLog log = LinearLog();
  const std::vector<double> t_abs = AllTimes();
  std::vector<double> t, bottom, top;
  ASSERT_TRUE(ExtractBulkTimeSeries(log, VariableType::TEMPERATURE, ParamSpatialLoc::Bottom, 0.0, t, bottom));
  ASSERT_TRUE(ExtractBulkTimeSeries(log, VariableType::TEMPERATURE, ParamSpatialLoc::Top, 0.0, t, top));
  for (size_t k = 0; k < t_abs.size(); ++k) {
    EXPECT_NEAR(bottom[k], CellLine(0.0, t_abs[k]), 1e-9);  // mean of the mirror ghost and cell 1
    EXPECT_NEAR(top[k], CellLine(0.675, t_abs[k]), 1e-9);   // zero-gradient outlet: the end cell
  }

  std::vector<double> y_bottom, y_top, y_mid;
  ASSERT_TRUE(ExtractComponentTimeSeries(log, ComponentVariableType::MOLE_FRACTION, "A", ParamSpatialLoc::Bottom,
                                         0.0, t, y_bottom));
  ASSERT_TRUE(ExtractComponentTimeSeries(log, ComponentVariableType::MOLE_FRACTION, "A", ParamSpatialLoc::Top, 0.0,
                                         t, y_top));
  ASSERT_TRUE(ExtractComponentTimeSeries(log, ComponentVariableType::MOLE_FRACTION, "A",
                                         ParamSpatialLoc::Position, 0.02, t, y_mid));
  for (size_t k = 0; k < t_abs.size(); ++k) {
    EXPECT_NEAR(y_bottom[k], FractionLine(0.0), 1e-12);
    EXPECT_NEAR(y_top[k], FractionLine(0.675), 1e-12);
    EXPECT_NEAR(y_mid[k], FractionLine(0.02), 1e-12);
  }
}

TEST(BedSeries, LinearFaceProfileReadsBackAtTheFaces) {
  const EquationLog log = LinearLog();
  const std::vector<double> t_abs = AllTimes();
  for (double z : {0.0, 0.1, 0.17, 0.3, 0.4, 0.8}) {
    std::vector<double> t, v;
    ASSERT_TRUE(ExtractBulkTimeSeries(log, VariableType::VELOCITY, ParamSpatialLoc::Position, z, t, v));
    for (size_t k = 0; k < v.size(); ++k) EXPECT_NEAR(v[k], FaceLine(z, t_abs[k]), 1e-9) << "z = " << z;
  }
  std::vector<double> t, bottom, top;
  ASSERT_TRUE(ExtractBulkTimeSeries(log, VariableType::VELOCITY, ParamSpatialLoc::Bottom, 0.0, t, bottom));
  ASSERT_TRUE(ExtractBulkTimeSeries(log, VariableType::VELOCITY, ParamSpatialLoc::Top, 0.0, t, top));
  for (size_t k = 0; k < t_abs.size(); ++k) {
    EXPECT_NEAR(bottom[k], FaceLine(0.0, t_abs[k]), 1e-12);
    EXPECT_NEAR(top[k], FaceLine(0.8, t_abs[k]), 1e-12);
  }
}

TEST(BedSeries, TimeBaseIsTheLoggedTimeFromTheCycleStart) {
  const EquationLog log = LinearLog();
  for (ParamSpatialLoc where : {ParamSpatialLoc::Bottom, ParamSpatialLoc::Top, ParamSpatialLoc::Position}) {
    std::vector<double> t, v;
    ASSERT_TRUE(ExtractBulkTimeSeries(log, VariableType::TEMPERATURE, where, 0.3, t, v));
    EXPECT_EQ(t, (std::vector<double>{0.0, 1.0, 2.5, 2.5, 4.0}));
    ASSERT_TRUE(ExtractBulkTimeSeries(log, VariableType::VELOCITY, where, 0.3, t, v));
    EXPECT_EQ(t, (std::vector<double>{0.0, 1.0, 2.5, 2.5, 4.0}));
  }
}

}  // namespace adsim
