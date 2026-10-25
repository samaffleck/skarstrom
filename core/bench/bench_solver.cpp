// Solver benchmark runner.
//
// Loads one or more saved decks (JSON) through the Simulation facade
// (adsim/simulation.h), solves each serially with profiling
// enabled, and prints wall time, IDA statistics and a wall-time breakdown so
// solver optimisations can be measured and compared run-to-run.
//
// Usage:
//   bench_solver [--repeat N] [--save-dir DIR] [--jacobian ida|band|coloured]
//                [--pattern FILE] deck.json [deck.json ...]
//
// --save-dir writes each deck's CSV logs (first repeat only) to DIR/<deck>/ so
// results can be diffed byte-for-byte against another build.
// --jacobian selects how the Newton Jacobian is built (default: solver default).
// --pattern appends the observed Jacobian coupling of each deck to FILE as
// "deck row_var col_var cell_offset count" lines.
// --unknowns full|ts|c|both selects the reduced-unknown switches (IDASettings
// alias_lte_solid_temperature_ / derived_concentration_; default: solver default).
// --band probe|full sizes the banded matrix from a start-of-run coupling probe
// (IDASettings band_from_probe_) or from the 2 neq each side default.
// --max-step / --nl-coef override the deck's IDA max step [s] and nonlinear
// convergence coefficient (for solver-settings experiments).
// --log low|medium|high prints the solver's log messages to stderr.
// --trace T prints, after every accepted step at t >= T, the step size, error
// norm and the boundary velocities/pressures (for diagnosing stalls).

#include "adsim/simulation.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

namespace {

  struct RunResult {
    int retval = 0;
    double wall_s = 0.0;
    double save_s = 0.0;
    adsim::IntegrationStats stats;
    adsim::SolverProfile profile;
    int colours = 0;
    long audits = 0;
    long failure_audits = 0;
    long misses = 0;
    long mu = 0;
    long ml = 0;
    long out_of_band = 0;
    bool widened = false;
  };

  struct Options {
    int repeat = 1;
    std::string save_dir;
    std::optional<adsim::JacobianMode> jacobian;
    std::string pattern_file;
    std::optional<double> max_step;
    std::optional<double> nl_coef;
    std::optional<adsim::LogLevel> log_level;
    std::optional<double> trace_from;
    std::optional<std::string> unknowns;  // full | ts | c | both
    std::optional<bool> band_probe;
  };

  RunResult RunDeck(const std::string& path, const std::string& save_dir, const Options& opt) {
    RunResult r;
    adsim::Simulation sim = adsim::Simulation::FromFile(path);
    adsim::UserData& data = sim.Data();
    data.profile_.enabled = true;
    if (opt.jacobian) data.settings.jacobian_mode_ = *opt.jacobian;
    data.settings.record_jacobian_pattern_ = !opt.pattern_file.empty();
    if (opt.max_step) data.settings.max_step = *opt.max_step;
    if (opt.nl_coef) data.settings.nonlin_conv_coef = *opt.nl_coef;
    if (opt.band_probe) data.settings.band_from_probe_ = *opt.band_probe;
    if (opt.unknowns) {
      data.settings.alias_lte_solid_temperature_ = (*opt.unknowns == "ts" || *opt.unknowns == "both");
      data.settings.derived_concentration_ = (*opt.unknowns == "c" || *opt.unknowns == "both");
    }

    adsim::RunOptions run;
    if (opt.log_level) {
      run.log_level = *opt.log_level;
      run.on_message = [](adsim::LogLevel, const std::string& msg) {
        std::fprintf(stderr, "%s\n", msg.c_str());
      };
    }
    if (opt.trace_from) {
      const double from = *opt.trace_from;
      run.on_step = [&data, &sim, from](double t, double dt, double err, int cycle, int step) {
        if (t < from) return;
        const auto& u = data.state.velocity.var;
        const auto& P = data.state.pressure.var;
        const auto& T = data.state.fluid_temperature.var;
        const size_t n = P.size();
        std::fprintf(stderr, "TRACE c%d s%d t=%.9f h=%.3e err=%.3f [%s] u0=%.3e u1=%.3e uN=%.3e P0=%.3f P1=%.3f PN=%.3f T1=%.3f\n",
                     cycle + 1, step, t, dt, err, sim.BlockingSolver()->DescribeLargestLocalError().c_str(),
                     u[0], u[1], u[u.size() - 1], P[0], P[1], P[n - 1], T[1]);
      };
    }
    const adsim::RunResult result = sim.Run(run);
    adsim::CoupledSolver& solver = *sim.BlockingSolver();
    r.retval = result.return_code;
    r.wall_s = result.wall_time_s;
    r.stats = result.stats;
    r.profile = data.profile_;
    r.colours = solver.GetJacobian().NumColours();
    r.audits = solver.GetJacobian().NumAudits();
    r.failure_audits = solver.GetJacobian().NumFailureAudits();
    r.misses = solver.GetJacobian().NumMisses();
    r.mu = solver.GetJacobian().UpperBandwidth();
    r.ml = solver.GetJacobian().LowerBandwidth();
    r.out_of_band = solver.GetJacobian().NumOutOfBand();
    r.widened = solver.GetSolverData().BandWidened();
    if (data.settings.jacobian_mode_ == adsim::JacobianMode::ColouredCheck) {
      std::printf("  check: %ld mismatched entries, max |diff|/colmax = %.3e\n",
                  solver.GetJacobian().NumCheckMismatches(), solver.GetJacobian().MaxCheckRelativeError());
      for (const auto& line : solver.GetJacobian().CheckReport()) std::printf("    %s\n", line.c_str());
    }

    if (!opt.pattern_file.empty()) {
      std::ofstream out(opt.pattern_file, std::ios::app);
      const std::string name = std::filesystem::path(path).stem().string();
      for (const auto& [key, count] : solver.GetJacobian().PatternCounts()) {
        out << name << ' ' << std::get<0>(key) << ' ' << std::get<1>(key) << ' '
            << std::get<2>(key) << ' ' << count << '\n';
      }
    }

    if (!save_dir.empty()) {
      std::filesystem::create_directories(save_dir);
      auto s0 = std::chrono::steady_clock::now();
      sim.SaveLogs(save_dir);
      r.save_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - s0).count();
    }
    return r;
  }

  double Pct(std::int64_t part_ns, double wall_s) {
    return wall_s > 0.0 ? 100.0 * (static_cast<double>(part_ns) * 1e-9) / wall_s : 0.0;
  }

  void PrintRun(const std::string& name, int rep, const RunResult& r) {
    const auto& s = r.stats;
    const auto& p = r.profile;
    const long total_res = s.num_res_evals + s.num_lin_res_evals;
    const double jac_share = total_res > 0 ? 100.0 * s.num_lin_res_evals / total_res : 0.0;
    const double us_per_res = p.residual_calls > 0 ? 1e-3 * p.residual_ns / p.residual_calls : 0.0;
    const std::int64_t ida_internal_ns = p.one_step_ns - p.residual_ns;
    const std::int64_t accounted = p.one_step_ns + p.log_ns + p.css_ns;
    const std::int64_t other_ns = static_cast<std::int64_t>(r.wall_s * 1e9) - accounted;

    std::printf("%-36s rep%d rc=%d wall=%8.3fs | steps=%7ld res=%8ld jacRes=%8ld (%4.1f%%) jacs=%6ld "
                "setups=%6ld nni=%7ld etf=%5ld ncf=%5ld | %6.2fus/res | "
                "resid=%4.1f%% idaInt=%4.1f%% log=%4.1f%% css=%4.1f%% other=%4.1f%% | "
                "atMaxStep=%ld/%ld bndETF=%ld/%ld cycSteps | colours=%d audits=%ld (after ncf %ld) misses=%ld band=%ld/%ld oob=%ld widened=%d | save=%.3fs\n",
                name.c_str(), rep, r.retval, r.wall_s,
                s.num_steps, s.num_res_evals, s.num_lin_res_evals, jac_share, s.num_jac_evals,
                s.num_lin_setups, s.num_nonlin_iters, s.num_err_test_fails, s.num_conv_fails,
                us_per_res,
                Pct(p.residual_ns, r.wall_s), Pct(ida_internal_ns, r.wall_s), Pct(p.log_ns, r.wall_s),
                Pct(p.css_ns, r.wall_s), Pct(other_ns, r.wall_s),
                p.steps_at_max_step, s.num_steps, p.boundary_err_test_fails, p.cycle_steps,
                r.colours, r.audits, r.failure_audits, r.misses, r.mu, r.ml, r.out_of_band, r.widened ? 1 : 0, r.save_s);
    std::fflush(stdout);
  }

} // namespace

int main(int argc, char** argv) {
  Options opt;
  std::vector<std::string> decks;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--repeat" && i + 1 < argc) {
      opt.repeat = std::max(1, std::atoi(argv[++i]));
    } else if (arg == "--save-dir" && i + 1 < argc) {
      opt.save_dir = argv[++i];
    } else if (arg == "--jacobian" && i + 1 < argc) {
      std::string mode = argv[++i];
      if (mode == "ida") opt.jacobian = adsim::JacobianMode::IdaInternal;
      else if (mode == "band") opt.jacobian = adsim::JacobianMode::BandDQ;
      else if (mode == "coloured") opt.jacobian = adsim::JacobianMode::Coloured;
      else if (mode == "check") opt.jacobian = adsim::JacobianMode::ColouredCheck;
      else { std::fprintf(stderr, "unknown --jacobian mode: %s\n", mode.c_str()); return 2; }
    } else if (arg == "--max-step" && i + 1 < argc) {
      opt.max_step = std::atof(argv[++i]);
    } else if (arg == "--log" && i + 1 < argc) {
      std::string level = argv[++i];
      opt.log_level = level == "high" ? adsim::LogLevel::HIGH
                    : level == "medium" ? adsim::LogLevel::MEDIUM : adsim::LogLevel::LOW;
    } else if (arg == "--trace" && i + 1 < argc) {
      opt.trace_from = std::atof(argv[++i]);
    } else if (arg == "--nl-coef" && i + 1 < argc) {
      opt.nl_coef = std::atof(argv[++i]);
    } else if (arg == "--pattern" && i + 1 < argc) {
      opt.pattern_file = argv[++i];
    } else if (arg == "--band" && i + 1 < argc) {
      std::string mode = argv[++i];
      if (mode == "probe") opt.band_probe = true;
      else if (mode == "full") opt.band_probe = false;
      else { std::fprintf(stderr, "unknown --band mode: %s\n", mode.c_str()); return 2; }
    } else if (arg == "--unknowns" && i + 1 < argc) {
      opt.unknowns = argv[++i];
      if (*opt.unknowns != "full" && *opt.unknowns != "ts" && *opt.unknowns != "c" && *opt.unknowns != "both") {
        std::fprintf(stderr, "unknown --unknowns mode: %s\n", opt.unknowns->c_str()); return 2;
      }
    } else {
      decks.push_back(arg);
    }
  }
  const int repeat = opt.repeat;
  const std::string& save_dir = opt.save_dir;

  if (decks.empty()) {
    std::fprintf(stderr, "usage: %s [--repeat N] [--save-dir DIR] deck.json [deck.json ...]\n", argv[0]);
    return 2;
  }

  int failures = 0;
  for (const auto& deck : decks) {
    const std::string name = std::filesystem::path(deck).stem().string();
    std::vector<double> walls;
    for (int rep = 0; rep < repeat; ++rep) {
      const std::string dir =
          (rep == 0 && !save_dir.empty()) ? (std::filesystem::path(save_dir) / name).string() : "";
      RunResult r = RunDeck(deck, dir, opt);
      PrintRun(name, rep, r);
      walls.push_back(r.wall_s);
      if (r.retval < 0) ++failures;
    }
    std::sort(walls.begin(), walls.end());
    std::printf("%-36s MEDIAN wall=%.3fs over %d run(s)\n\n", name.c_str(), walls[walls.size() / 2], repeat);
  }
  return failures == 0 ? 0 : 1;
}
