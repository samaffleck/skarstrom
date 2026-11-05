// test_fine_grid_breakthrough.cpp
//
// Regression test for a breakthrough deck that stalled on a fine grid: the
// example deck Ward2022_Q50ccm (pure CO2 at 50 ccm into a helium-filled bed of
// activated carbon, Ergun momentum, pressure outlet at 1 bar) ran on 10 to
// 160 cells but failed on 320 cells with every advection scheme ("too much
// work": 50 000 steps, 24 000 Newton convergence failures, 72 000 Jacobians)
// as the CO2 front left the column (t ~ 2735 s), the step size collapsing to
// 1e-6 s with the largest error in the outlet-face velocity.
//
// The cause was the difference-quotient Jacobian, not the model. IDA's DQ
// increment for a column is floored at reltol |y| + abstol, which for the
// pressure is 0.1 Pa at 1 bar. The pressure outlet imposes the outlet-face
// velocity as the explicit Ergun root u = ErgunVelocity(2 (P_end - P_bc), A, B),
// which is linear in the end-cell pressure only on the scale A^2 / (4 B), and
// that scale shrinks with the cell width (A, B ~ dz): 0.12 Pa on the deck's
// 60 cells, 0.02 Pa on 320. A 0.1 Pa secant across that bend is half the true
// slope on 320 cells, the modified Newton iteration overshoots the end-cell
// pressure about twofold, and it diverges unless the step is tiny. The
// pressure columns now use sqrt(eps) |P| (1.5e-3 Pa at 1 bar; see Increment in
// solver/jacobian.cpp).
//
// The test runs the deck as shipped on 80, 160 and 320 cells and checks:
//   * every run completes, with Newton convergence failures on fewer than 5 %
//     of the steps (before the fix: 48 % on 320 cells, and the run failed);
//   * the 50 % breakthrough time converges at first order (upwind): the
//     320 - 160 difference is at most 0.75 of the 160 - 80 difference (ideal
//     first order: 0.5), and below 0.5 % of the front time;
//   * van Leer completes on 320 cells as well and agrees with upwind's front
//     time there to 0.5 %.
// Guard: on 320 cells the outlet face's linearity scale A^2 / (4 B), from the
// Ergun coefficients with the logged gas state, is less than a quarter of the
// old increment reltol P, so the run is in the regime that stalled.
// Runtime: about 8 s on a Release build.

#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/simulation.h"

namespace adsim {

namespace {

const std::filesystem::path kDeck =
    std::filesystem::current_path() / ".." / "examples" / "decks" / "Ward2022_Q50ccm.json";

struct FineGridRun {
  RunResult result;
  double t50 = std::nan("");   // [s] outlet y_CO2 = 0.5, by linear interpolation
  double linearity_scale = 0;  // [Pa] A^2 / (4 B) of the outlet face
  double old_increment = 0;    // [Pa] reltol P_out + abstol, IDA's DQ floor for P
};

double CrossingTime(const TimeSeries& s, double level) {
  for (size_t k = 1; k < s.values.size(); ++k) {
    const double a = s.values[k - 1], b = s.values[k];
    if (a < level && b >= level)
      return s.time[k - 1] + (level - a) / (b - a) * (s.time[k] - s.time[k - 1]);
  }
  return std::nan("");
}

FineGridRun RunWard(int n_cells, AdvectionScheme scheme) {
  Simulation sim = Simulation::FromFile(kDeck.string());
  sim.SetParameter("layer[1].cells", n_cells);
  sim.Data().settings.advection_scheme_ = scheme;

  FineGridRun r;
  r.result = sim.Run();
  if (!r.result.Ok()) return r;
  r.t50 = CrossingTime(sim.Series("mole_fraction", "CarbonDioxide", SeriesLocation::Top), 0.5);

  // Ergun coefficients of the outlet face (the last bed cell's half-cells in
  // series with the ghost's, which takes the interior gas state): viscous
  // A = mu dz / K with K = dp^2 eb^3 / (150 (1 - eb)^2), inertial
  // B = beta rho dz with beta = 1.75 (1 - eb) / (eb^3 dp).
  const UserData& data = sim.Data();
  const auto& layer = data.bed.layers_.layers_.front();
  const double eb = layer.eb_.Get();
  const double dp = layer.particle_diameter_.Get();
  const double dz = layer.domain_.length.Get() / n_cells;
  const auto& state = data.state;
  const size_t end = state.pressure.var.size() - 2;  // last bed cell
  const double mu = state.viscosity.var[end];
  const double rho = state.fluid_density.var[end];
  const double K = dp * dp * eb * eb * eb / (150.0 * (1.0 - eb) * (1.0 - eb));
  const double beta = 1.75 * (1.0 - eb) / (eb * eb * eb * dp);
  const double A = mu * dz / K;
  const double B = beta * rho * dz;
  r.linearity_scale = A * A / (4.0 * B);
  r.old_increment = data.settings.reltol.Get() * state.pressure.var[end] + state.pressure.abstol;
  return r;
}

void ExpectCompletedCleanly(const FineGridRun& r, const std::string& label) {
  ASSERT_TRUE(r.result.Ok()) << label << ": " << r.result.message;
  const auto& s = r.result.stats;
  EXPECT_LT(s.num_conv_fails, 0.05 * s.num_steps)
      << label << ": " << s.num_conv_fails << " convergence failures in " << s.num_steps << " steps";
  ASSERT_TRUE(std::isfinite(r.t50)) << label << ": the CO2 front did not reach the outlet";
  std::cout << "[FineGridBreakthrough] " << label << ": t50 = " << r.t50 << " s, steps "
            << s.num_steps << ", conv fails " << s.num_conv_fails << ", Jacobians "
            << s.num_jac_evals << ", wall " << r.result.wall_time_s << " s" << std::endl;
}

}  // namespace

TEST(FineGridBreakthrough, Ward2022Q50ccmConvergesTo320Cells) {
  const FineGridRun r80 = RunWard(80, AdvectionScheme::Upwind);
  const FineGridRun r160 = RunWard(160, AdvectionScheme::Upwind);
  const FineGridRun r320 = RunWard(320, AdvectionScheme::Upwind);
  ExpectCompletedCleanly(r80, "upwind, 80 cells");
  ExpectCompletedCleanly(r160, "upwind, 160 cells");
  ExpectCompletedCleanly(r320, "upwind, 320 cells");

  // Guard: the fine grid is in the regime that stalled.
  EXPECT_LT(r320.linearity_scale, 0.25 * r320.old_increment)
      << "outlet face linear over " << r320.linearity_scale << " Pa, old increment "
      << r320.old_increment << " Pa";
  std::cout << "[FineGridBreakthrough] 320 cells: outlet face linear over A^2/(4B) = "
            << r320.linearity_scale << " Pa, IDA's pressure increment " << r320.old_increment << " Pa"
            << std::endl;

  const double d_coarse = std::abs(r160.t50 - r80.t50);
  const double d_fine = std::abs(r320.t50 - r160.t50);
  EXPECT_LE(d_fine, 0.75 * d_coarse)
      << "t50 on 80/160/320 cells: " << r80.t50 << ", " << r160.t50 << ", " << r320.t50 << " s";
  EXPECT_LT(d_fine, 0.005 * r320.t50)
      << "t50 on 160/320 cells: " << r160.t50 << ", " << r320.t50 << " s";

  const FineGridRun vl320 = RunWard(320, AdvectionScheme::VanLeer);
  ExpectCompletedCleanly(vl320, "van Leer, 320 cells");
  EXPECT_LT(std::abs(vl320.t50 - r320.t50), 0.005 * r320.t50)
      << "t50 on 320 cells: upwind " << r320.t50 << " s, van Leer " << vl320.t50 << " s";
}

}  // namespace adsim
