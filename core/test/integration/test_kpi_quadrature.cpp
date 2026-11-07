// test_kpi_quadrature.cpp
//
// KPI quadratures (adsim/solver/kpi_quadrature.h): the endpoint flows, the
// pump energy and the face fluxes integrated by IDAS as quadrature states
// beside the solution. With quadrature error control off (the default) they
// must leave the run untouched; the quadrature KPIs must agree with the
// trapezoid KPIs; the step mole balance (inventory difference against the
// integrated flows) must close; and the buffer-tank summary rebuilt from the
// quadratures must match the one from the logged samples.

#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/simulation.h"
#include "adsim/solver/kpis.h"
#include "adsim/solver/kpi_quadrature.h"

namespace adsim {

namespace {

const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";

// Perez 2019 A (four-step CO2 VSA with a vacuum pump, CO2 tracked) on 20
// cells for three cycles: a fast-tier run (< 1 s) that exercises feed,
// blowdown and evacuation endpoints, one-open-end steps and pump energy.
Simulation ShortPerez(int cycles = 3) {
  Simulation sim = Simulation::FromFile((kExamples / "Perez2019_A.json").string());
  sim.SetParameter("cycle.number_of_cycles", cycles);
  sim.SetParameter("css.method", 0);
  for (auto& layer : sim.Data().bed.layers_.layers_) layer.domain_.n_cells = 20;
  return sim;
}

void ExpectRelNear(double got, double want, double rel, const char* what) {
  ASSERT_TRUE(std::isfinite(got)) << what;
  ASSERT_TRUE(std::isfinite(want)) << what;
  EXPECT_LE(std::abs(got - want), rel * std::abs(want)) << what << ": " << got << " vs " << want;
}

}  // namespace

// Error control off: the quadratures do not reach IDA's step control, so the
// solution, the logs and the reported KPIs are bit-identical to a run
// without them.
TEST(KpiQuadrature, LeavesTheRunUnchanged) {
  Simulation with = ShortPerez();
  Simulation without = ShortPerez();
  without.Data().settings.kpi_quadratures_ = false;
  const RunResult a = with.Run();
  const RunResult b = without.Run();
  ASSERT_TRUE(a.Ok()) << a.message;
  ASSERT_TRUE(b.Ok()) << b.message;
  EXPECT_EQ(a.stats.num_steps, b.stats.num_steps);
  EXPECT_EQ(a.stats.num_res_evals, b.stats.num_res_evals);
  EXPECT_EQ(a.stats.num_nonlin_iters, b.stats.num_nonlin_iters);
  EXPECT_EQ(a.stats.num_err_test_fails, b.stats.num_err_test_fails);

  const auto ka = with.KPIs();
  const auto kb = without.KPIs();
  ASSERT_EQ(ka.size(), kb.size());
  for (size_t i = 0; i < ka.size(); ++i) {
    EXPECT_EQ(ka[i].purity, kb[i].purity);
    EXPECT_EQ(ka[i].recovery, kb[i].recovery);
    EXPECT_EQ(ka[i].productivity_mass, kb[i].productivity_mass);
    EXPECT_EQ(ka[i].productivity_volume, kb[i].productivity_volume);
    EXPECT_EQ(ka[i].energy_consumption, kb[i].energy_consumption);
    EXPECT_EQ(ka[i].production_rate, kb[i].production_rate);
    // Without the quadratures the quadrature KPIs are undefined.
    EXPECT_TRUE(std::isnan(kb[i].quadrature_purity));
    EXPECT_TRUE(std::isnan(kb[i].quadrature_mass_balance_error));
  }
  const TimeSeries sa = with.Series("mole_fraction", "CarbonDioxide", SeriesLocation::Top);
  const TimeSeries sb = without.Series("mole_fraction", "CarbonDioxide", SeriesLocation::Top);
  EXPECT_EQ(sa.time, sb.time);
  EXPECT_EQ(sa.values, sb.values);
}

// The quadrature KPIs integrate the same integrands as the trapezoid over the
// logged samples (200+ per step), so the two agree to the trapezoid's own
// error; both close the cycle mole balance.
TEST(KpiQuadrature, AgreesWithTheTrapezoidKpis) {
  Simulation sim = ShortPerez();
  const RunResult r = sim.Run();
  ASSERT_TRUE(r.Ok()) << r.message;
  const auto kpis = sim.KPIs();
  ASSERT_EQ(kpis.size(), 3u);
  for (const auto& k : kpis) {
    SCOPED_TRACE("cycle " + std::to_string(k.cycle));
    ASSERT_TRUE(k.valid);
    ExpectRelNear(k.quadrature_purity, k.purity, 1e-4, "purity");
    ExpectRelNear(k.quadrature_recovery, k.recovery, 1e-4, "recovery");
    ExpectRelNear(k.quadrature_productivity_mass, k.productivity_mass, 1e-4, "productivity");
    ExpectRelNear(k.quadrature_productivity_volume, k.productivity_volume, 1e-4, "productivity (volume)");
    ExpectRelNear(k.quadrature_production_rate, k.production_rate, 1e-4, "production rate");
    ExpectRelNear(k.quadrature_energy_consumption, k.energy_consumption, 1e-3, "energy");
    ASSERT_TRUE(std::isfinite(k.mass_balance_error));
    ASSERT_TRUE(std::isfinite(k.quadrature_mass_balance_error));
    EXPECT_LT(std::abs(k.quadrature_mass_balance_error), 1e-4);
    EXPECT_LT(std::abs(k.mass_balance_error), 1e-3);
  }
}

// Per-step accumulation form of the mole balance: what entered through the
// endpoints (and through the faces) equals the change of the column's gas and
// adsorbed inventory. Perez 2019 A has one-open-end steps (blowdown and
// evacuation), where this checks the open end's outflow on its own.
TEST(KpiQuadrature, StepInventoryBalanceCloses) {
  Simulation sim = ShortPerez(1);
  const RunResult r = sim.Run();
  ASSERT_TRUE(r.Ok()) << r.message;
  const BoundaryLog& b_log = sim.Boundaries();
  const int n_steps = static_cast<int>(sim.Data().cycle.step_list_.size());
  int one_open_end = 0;
  for (int s = 0; s < n_steps; ++s) {
    SCOPED_TRACE("step " + std::to_string(s));
    const StepBalance* b = b_log.GetStepBalance(0, s);
    ASSERT_NE(b, nullptr);
    ASSERT_TRUE(b->valid);
    if (b->open_ends == 1) ++one_open_end;
    const double scale = std::max(std::abs(b->TotalInventoryChange()), b->TotalThroughput());
    ASSERT_GT(scale, 0.0);
    EXPECT_LT(std::abs(b->EndpointDiscrepancy()) / scale, 1e-5);
    EXPECT_LT(std::abs(b->FaceDiscrepancy()) / scale, 1e-5);
    EXPECT_LT(std::abs(b->RelativeEndpointDiscrepancy()), 1e-5);
  }
  EXPECT_GE(one_open_end, 1);
}

// The buffer-tank summary (total moles, flow-weighted temperature and
// composition) rebuilt from an endpoint's quadratures matches the one from
// its logged samples.
TEST(KpiQuadrature, BufferSummaryFromTheQuadratures) {
  Simulation sim = ShortPerez(1);
  const RunResult r = sim.Run();
  ASSERT_TRUE(r.Ok()) << r.message;
  const BoundaryLog& b_log = sim.Boundaries();
  int checked = 0;
  for (int s = 0; s < static_cast<int>(b_log.sys_data.step_list_.size()); ++s) {
    for (const auto& name : b_log.stream_names.at(b_log.sys_data.step_list_[s])) {
      const StreamData* d = b_log.GetStreamData(0, s, name);
      ASSERT_NE(d, nullptr);
      ASSERT_TRUE(d->integrals_.valid);
      const CycleEndpointSummary trap = ComputeSummary(*d, 0);
      const CycleEndpointSummary quad = ComputeQuadratureSummary(d->integrals_, d->location_, 0);
      if (!trap.is_valid || std::abs(trap.total_moles) < 1e-6) continue;
      SCOPED_TRACE(name);
      ASSERT_TRUE(quad.is_valid);
      ExpectRelNear(quad.total_moles, trap.total_moles, 1e-4, "total moles");
      EXPECT_NEAR(quad.avg_temperature, trap.avg_temperature, 1e-2);
      for (const auto& [comp, y] : trap.avg_mole_fraction) EXPECT_NEAR(quad.avg_mole_fraction.at(comp), y, 1e-4);
      ++checked;
    }
  }
  EXPECT_GE(checked, 2);
}

// The chunked solver (the GUI path) records the same quadratures.
TEST(KpiQuadrature, ChunkedRunRecordsTheSameQuadratures) {
  Simulation blocking = ShortPerez(2);
  ASSERT_TRUE(blocking.Run().Ok());
  Simulation chunked = ShortPerez(2);
  ASSERT_TRUE(chunked.Start());
  RunStatus st = RunStatus::Running;
  for (int i = 0; i < 100000 && st == RunStatus::Running; ++i) st = chunked.Advance(20);
  ASSERT_EQ(st, RunStatus::Completed) << chunked.LastResult().message;
  const auto a = blocking.KPIs();
  const auto b = chunked.KPIs();
  ASSERT_EQ(a.size(), b.size());
  for (size_t i = 0; i < a.size(); ++i) {
    EXPECT_NEAR(b[i].quadrature_purity, a[i].quadrature_purity, 1e-12);
    EXPECT_NEAR(b[i].quadrature_recovery, a[i].quadrature_recovery, 1e-12);
    EXPECT_NEAR(b[i].quadrature_mass_balance_error, a[i].quadrature_mass_balance_error, 1e-12);
  }
}

// With quadrature error control on, the quadratures join IDA's error test;
// the run still completes and the KPIs move by no more than the tolerance.
TEST(KpiQuadrature, ErrorControlOn) {
  Simulation off = ShortPerez(1);
  Simulation on = ShortPerez(1);
  on.Data().settings.quadrature_error_control_ = true;
  const RunResult a = off.Run();
  const RunResult b = on.Run();
  ASSERT_TRUE(a.Ok()) << a.message;
  ASSERT_TRUE(b.Ok()) << b.message;
  EXPECT_GE(b.stats.num_steps, a.stats.num_steps);
  const auto ka = off.KPIs();
  const auto kb = on.KPIs();
  ASSERT_EQ(ka.size(), kb.size());
  for (size_t i = 0; i < ka.size(); ++i) {
    ExpectRelNear(kb[i].quadrature_purity, ka[i].quadrature_purity, 1e-4, "purity");
    ExpectRelNear(kb[i].quadrature_recovery, ka[i].quadrature_recovery, 1e-4, "recovery");
  }
}

}  // namespace adsim
