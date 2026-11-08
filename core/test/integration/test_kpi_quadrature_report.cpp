// test_kpi_quadrature_report.cpp
//
// Opt-in measurement of the KPI quadratures on the example decks (which the
// paper suites write): the reported (quadrature) KPIs against the trapezoid
// KPIs at the last cycle,
// the trapezoid's own error (its change when the log interval is halved), the
// cycle mole-balance error with each, whether the run with the quadratures is
// bit-identical to one without, the cost of the quadrature states, and what
// quadrature error control changes. Skipped unless
//
//   ADSIM_KPI_QUADRATURE_REPORT=<out.csv>   [ADSIM_KPI_QUADRATURE_DECKS=A.json,B.json]
//
// is set; each deck runs four times as shipped (to its CSS or cycle count).

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/simulation.h"
#include "adsim/solver/kpis.h"

namespace adsim {

namespace {

const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";

struct DeckRun {
  RunResult result;
  std::vector<KPIRecord> kpis;  // last cycle
  double mb_trap = std::numeric_limits<double>::quiet_NaN();  // max |cycle mole-balance error| over components
  double mb_quad = std::numeric_limits<double>::quiet_NaN();
  double step_disc = std::numeric_limits<double>::quiet_NaN();  // max |relative step discrepancy|, one open end
  double mb_face = std::numeric_limits<double>::quiet_NaN();  // cycle balance from the face-flux quadratures
  // max |cycle balance| over components from the face-flux quadratures
  double mb_face_comp = std::numeric_limits<double>::quiet_NaN();
};

DeckRun RunOne(const std::string& deck, bool quadratures, int log_points, bool error_control) {
  Simulation sim = Simulation::FromFile((kExamples / deck).string());
  sim.Data().settings.kpi_quadratures_ = quadratures;
  sim.Data().settings.quadrature_error_control_ = error_control;
  sim.Data().settings.max_log_points_per_step_ = log_points;
  DeckRun out;
  out.result = sim.Run();
  if (!sim.HasResults()) return out;
  const int last = sim.LastCycle();
  out.kpis = sim.KPIs(last);
  const auto& eq = sim.Log();
  const auto& b = sim.Boundaries();
  auto worse = [](double acc, double v) { return std::isnan(acc) ? std::abs(v) : std::max(acc, std::abs(v)); };
  for (const auto& c : sim.Data().fluid.components_) {
    const double t = CycleMassBalanceError(eq, b, last, c, KPIIntegration::TRAPEZOID);
    const double q = CycleMassBalanceError(eq, b, last, c, KPIIntegration::QUADRATURE);
    if (std::isfinite(t)) out.mb_trap = worse(out.mb_trap, t);
    if (std::isfinite(q)) out.mb_quad = worse(out.mb_quad, q);
  }
  const int n_steps = static_cast<int>(b.sys_data.step_list_.size());
  for (int s = 0; s < n_steps; ++s) {
    const StepBalance* sb = b.GetStepBalance(last, s);
    if (sb && sb->valid && sb->open_ends == 1) out.step_disc = worse(out.step_disc, sb->RelativeEndpointDiscrepancy());
  }
  // Total-mole cycle balance from the face fluxes the species balance uses.
  const StepBalance* first = b.GetStepBalance(last, 0);
  const StepBalance* end = b.GetStepBalance(last, n_steps - 1);
  if (first && end && first->valid && end->valid) {
    double net = 0.0, entered = 0.0;
    bool ok = true;
    for (int s = 0; s < n_steps; ++s) {
      const StepBalance* sb = b.GetStepBalance(last, s);
      if (!sb || !sb->valid) { ok = false; break; }
      for (size_t c = 0; c < sb->components.size(); ++c) {
        net += sb->FaceNetIn(c);
        entered += std::max(sb->face_in_bottom[c], 0.0) + std::max(sb->face_in_top[c], 0.0);
      }
    }
    double acc = 0.0;
    for (size_t c = 0; c < first->components.size(); ++c) acc += end->inventory_end[c] - first->inventory_start[c];
    if (ok && entered > 0.0) out.mb_face = (net - acc) / entered;
    // Per component, as CycleMassBalanceError but with the face fluxes.
    double inventory = 0.0;
    for (double v : first->inventory_start) inventory += std::abs(v);
    for (size_t c = 0; ok && c < first->components.size(); ++c) {
      double net_c = 0.0, in_c = 0.0;
      for (int s = 0; s < n_steps; ++s) {
        const StepBalance* sb = b.GetStepBalance(last, s);
        net_c += sb->FaceNetIn(c);
        in_c += std::max(sb->face_in_bottom[c], 0.0) + std::max(sb->face_in_top[c], 0.0);
      }
      if (!(in_c > 1e-12 * inventory)) continue;
      const double gen = [&] {
        double g = 0.0;
        for (int s = 0; s < n_steps; ++s) g += eq.ReactionGeneration(last, s, first->components[c]);
        return g;
      }();
      out.mb_face_comp =
          worse(out.mb_face_comp, (net_c + gen - (end->inventory_end[c] - first->inventory_start[c])) / in_c);
    }
  }
  return out;
}

std::string Num(double v) {
  std::ostringstream ss;
  ss << std::setprecision(10) << v;
  return ss.str();
}

double Rel(double a, double b) {
  if (!std::isfinite(a) || !std::isfinite(b)) return std::numeric_limits<double>::quiet_NaN();
  const double s = std::max(std::abs(a), std::abs(b));
  return s > 0.0 ? std::abs(a - b) / s : 0.0;
}

}  // namespace

TEST(KpiQuadratureReport, ExampleDecks) {
  const char* out_path = std::getenv("ADSIM_KPI_QUADRATURE_REPORT");
  if (!out_path) GTEST_SKIP() << "set ADSIM_KPI_QUADRATURE_REPORT=<out.csv> to run";
  std::vector<std::string> decks;
  if (const char* list = std::getenv("ADSIM_KPI_QUADRATURE_DECKS")) {
    std::stringstream ss(list);
    for (std::string d; std::getline(ss, d, ',');) if (!d.empty()) decks.push_back(d);
  } else {
    for (const auto& e : std::filesystem::directory_iterator(kExamples))
      if (e.path().extension() == ".json") decks.push_back(e.path().filename().string());
    std::sort(decks.begin(), decks.end());
  }

  std::ofstream csv(out_path);
  csv << "deck,status,cycles,css_cycle,component,kpi,trapezoid,quadrature,trapezoid_400,rel_quad_vs_trap,"
         "rel_trap_halving,mb_trap,mb_quad,mb_face_total,mb_face_component,step_disc_one_open_end,identical,steps_on,steps_off,res_on,res_off,"
         "wall_on,wall_off,steps_errcon,res_errcon,quad_errcon,rel_errcon\n";
  const char* names[] = {"purity", "recovery", "productivity_mass", "productivity_volume", "energy_consumption",
                         "production_rate"};
  for (const auto& deck : decks) {
    SCOPED_TRACE(deck);
    const DeckRun on = RunOne(deck, true, 200, false);
    const DeckRun off = RunOne(deck, false, 200, false);
    const DeckRun fine = RunOne(deck, true, 400, false);
    const DeckRun ec = RunOne(deck, true, 200, true);

    bool identical = on.result.stats.num_steps == off.result.stats.num_steps
                     && on.result.stats.num_res_evals == off.result.stats.num_res_evals
                     && on.kpis.size() == off.kpis.size();
    for (size_t i = 0; identical && i < on.kpis.size(); ++i) {
      const auto& a = on.kpis[i];
      const auto& b = off.kpis[i];
      auto same = [](double x, double y) { return (std::isnan(x) && std::isnan(y)) || x == y; };
      identical = same(a.trapezoid_purity, b.trapezoid_purity) && same(a.trapezoid_recovery, b.trapezoid_recovery)
                  && same(a.trapezoid_productivity_mass, b.trapezoid_productivity_mass)
                  && same(a.trapezoid_energy_consumption, b.trapezoid_energy_consumption)
                  && same(a.trapezoid_production_rate, b.trapezoid_production_rate);
    }
    EXPECT_TRUE(identical);

    auto row_prefix = [&](const std::string& comp, const std::string& kpi) {
      std::ostringstream ss;
      ss << deck << "," << RunStatusName(on.result.status) << "," << on.result.cycles_run << ","
         << on.result.css_cycle << "," << comp << "," << kpi << ",";
      return ss.str();
    };
    auto row_suffix = [&]() {
      std::ostringstream ss;
      ss << Num(on.mb_trap) << "," << Num(on.mb_quad) << "," << Num(on.mb_face) << "," << Num(on.mb_face_comp) << "," << Num(on.step_disc) << "," << (identical ? 1 : 0) << ","
         << on.result.stats.num_steps << "," << off.result.stats.num_steps << "," << on.result.stats.num_res_evals
         << "," << off.result.stats.num_res_evals << "," << Num(on.result.wall_time_s) << ","
         << Num(off.result.wall_time_s) << "," << ec.result.stats.num_steps << "," << ec.result.stats.num_res_evals;
      return ss.str();
    };
    if (on.kpis.empty()) {
      csv << row_prefix("", "") << ",,,,," << row_suffix() << ",,\n";
    }
    for (size_t i = 0; i < on.kpis.size(); ++i) {
      const KPIRecord& k = on.kpis[i];
      const KPIRecord* kf = i < fine.kpis.size() ? &fine.kpis[i] : nullptr;
      const KPIRecord* ke = i < ec.kpis.size() ? &ec.kpis[i] : nullptr;
      const double trap[] = {k.trapezoid_purity, k.trapezoid_recovery, k.trapezoid_productivity_mass,
                             k.trapezoid_productivity_volume, k.trapezoid_energy_consumption,
                             k.trapezoid_production_rate};
      const double quad[] = {k.purity, k.recovery, k.productivity_mass, k.productivity_volume, k.energy_consumption,
                             k.production_rate};
      const double nan = std::numeric_limits<double>::quiet_NaN();
      const double fine_v[] = {kf ? kf->trapezoid_purity : nan, kf ? kf->trapezoid_recovery : nan,
                               kf ? kf->trapezoid_productivity_mass : nan, kf ? kf->trapezoid_productivity_volume : nan,
                               kf ? kf->trapezoid_energy_consumption : nan, kf ? kf->trapezoid_production_rate : nan};
      const double ec_v[] = {ke ? ke->purity : nan, ke ? ke->recovery : nan, ke ? ke->productivity_mass : nan,
                             ke ? ke->productivity_volume : nan, ke ? ke->energy_consumption : nan,
                             ke ? ke->production_rate : nan};
      for (int j = 0; j < 6; ++j) {
        csv << row_prefix(k.component, names[j]) << Num(trap[j]) << "," << Num(quad[j]) << "," << Num(fine_v[j]) << ","
            << Num(Rel(quad[j], trap[j])) << "," << Num(Rel(fine_v[j], trap[j])) << "," << row_suffix() << ","
            << Num(ec_v[j]) << "," << Num(Rel(ec_v[j], quad[j])) << "\n";
      }
    }
    csv.flush();
  }
}

}  // namespace adsim
