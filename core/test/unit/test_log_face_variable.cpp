// test_log_face_variable.cpp
//
// Cell-centre values interpolated to the cell faces by the equation log
// (EquationLog::GetFaceVariable / GetComponentFaceVariable). A hand-built log
// on a non-uniform mesh (a void-like thin cell, a fine layer and a coarse
// layer) holds a linear cell profile: every interior face must read back the
// line exactly, which needs each side weighted by the OTHER cell's width.
// The two boundary faces must be the boundary conditions' face value, the
// plain mean of the ghost and the first/last interior cell, whatever width
// the ghost carries (it copies the first/last bed cell, not its neighbour).

#include <vector>

#include <gtest/gtest.h>

#include "adsim/equations/equation.h"
#include "adsim/solver/log.h"

namespace adsim {

namespace {

// Ghost (width of the first bed cell), a thin 0.02 m void-like cell, two
// cells of 0.1 m, two of 0.25 m, ghost: N = 7 cells. Faces at 0, 0.02, 0.12,
// 0.22, 0.47, 0.72 (positions[i] is the face between cells i and i+1).
DomainVectors Mesh() {
  DomainVectors d;
  d.cell_width = {0.1, 0.02, 0.1, 0.1, 0.25, 0.25, 0.25};
  const size_t N = d.cell_width.size();
  d.positions.assign(N - 1, 0.0);
  for (size_t i = 1; i + 1 < N; ++i) d.positions[i] = d.positions[i - 1] + d.cell_width[i];
  return d;
}

// Interior cell centres, from the face positions.
std::vector<double> Centres(const DomainVectors& d) {
  std::vector<double> z(d.cell_width.size(), 0.0);
  for (size_t i = 1; i + 1 < z.size(); ++i) z[i] = d.positions[i - 1] + 0.5 * d.cell_width[i];
  return z;
}

double Line(double z, double t) { return 2.0e5 + 3.0e4 * z + 10.0 * t; }
double FractionLine(double z) { return 0.1 + 0.5 * z; }

constexpr double kGhostBottom = 7.0;  // arbitrary ghost values: the BC sets them
constexpr double kGhostTop = 11.0;

std::vector<double> Row(const DomainVectors& d, double t, bool fraction) {
  const std::vector<double> z = Centres(d);
  std::vector<double> row(z.size());
  for (size_t i = 1; i + 1 < z.size(); ++i) row[i] = fraction ? FractionLine(z[i]) : Line(z[i], t);
  row.front() = kGhostBottom;
  row.back() = kGhostTop;
  return row;
}

const std::vector<double> kTimes = {0.0, 1.5, 4.0};

EquationLog LinearLog() {
  const DomainVectors mesh = Mesh();
  Equations state;
  state.domain_vectors = mesh;
  state.mole_fraction.emplace("A", Equation("Mole Fraction", CellType::CellCentre, EquationType::Differential, {}));

  EquationLog log;
  log.sys_data.step_list_ = {"Feed"};
  for (double t : kTimes) {
    state.pressure.var = Row(mesh, t, false);
    state.velocity.var.assign(mesh.positions.size(), 0.0);
    state.mole_fraction.at("A").var = Row(mesh, t, true);
    log.LogState(state, t, 0, 0);
  }
  return log;
}

}  // namespace

TEST(EquationLogFaces, LinearProfileIsExactAtInteriorFacesOfANonUniformMesh) {
  const DomainVectors mesh = Mesh();
  const EquationLog log = LinearLog();
  const VariableLogEntry P = log.GetFaceVariable(0, 0, VariableType::PRESSURE);
  ASSERT_EQ(P.values.size(), kTimes.size());
  for (size_t k = 0; k < kTimes.size(); ++k) {
    const std::vector<double>& row = P.values[k];
    ASSERT_EQ(row.size(), mesh.positions.size());
    // Interior faces 1..N-3: between two interior cells, including the
    // thin/fine and fine/coarse faces where the widths differ.
    for (size_t f = 1; f + 1 < row.size(); ++f) {
      EXPECT_NEAR(row[f], Line(mesh.positions[f], kTimes[k]), 1e-6) << "face " << f << ", row " << k;
    }
  }
}

TEST(EquationLogFaces, BoundaryFacesAreTheMeanOfGhostAndEndCell) {
  const DomainVectors mesh = Mesh();
  const EquationLog log = LinearLog();
  const VariableLogEntry P = log.GetFaceVariable(0, 0, VariableType::PRESSURE);
  for (size_t k = 0; k < kTimes.size(); ++k) {
    const std::vector<double> cells = Row(mesh, kTimes[k], false);
    const std::vector<double>& row = P.values[k];
    EXPECT_NEAR(row.front(), 0.5 * (cells[0] + cells[1]), 1e-6);
    EXPECT_NEAR(row.back(), 0.5 * (cells[cells.size() - 2] + cells.back()), 1e-6);
  }
}

TEST(EquationLogFaces, ComponentFaceVariableUsesTheSameWeights) {
  const DomainVectors mesh = Mesh();
  const EquationLog log = LinearLog();
  const VariableLogEntry y = log.GetComponentFaceVariable(0, 0, ComponentVariableType::MOLE_FRACTION, "A");
  ASSERT_EQ(y.values.size(), kTimes.size());
  for (const std::vector<double>& row : y.values) {
    ASSERT_EQ(row.size(), mesh.positions.size());
    for (size_t f = 1; f + 1 < row.size(); ++f) {
      EXPECT_NEAR(row[f], FractionLine(mesh.positions[f]), 1e-12) << "face " << f;
    }
  }
}

TEST(EquationLogFaces, FaceCentredVariablesAreReturnedAsLogged) {
  const DomainVectors mesh = Mesh();
  const EquationLog log = LinearLog();
  const VariableLogEntry u = log.GetFaceVariable(0, 0, VariableType::VELOCITY);
  ASSERT_EQ(u.values.size(), kTimes.size());
  for (const std::vector<double>& row : u.values) EXPECT_EQ(row.size(), mesh.positions.size());
}

}  // namespace adsim
