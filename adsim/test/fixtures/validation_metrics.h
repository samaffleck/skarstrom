#pragma once

// ============================================================================
//  Quantitative paper-validation helpers.
//
//  Every paper test in adsim/test/papers/ compares the simulation with numbers
//  published in (or digitised from) the paper, using an explicit error metric,
//  and records each comparison as a row of a machine-readable summary:
//
//      paper, case, quantity, simulated, reference, error, metric, tolerance,
//      pass, units, source, note
//
//  * Reference values live in adsim/test/papers/<paper>/reference_values.csv
//    (one row per published number, with the table / figure / page it was
//    transcribed from). Digitised curves live in adsim/data/<paper>/.
//  * Each test case writes ONE file, <dir>/<paper>__<case>.json, so cases
//    running in parallel never race on a shared file. <dir> is
//    $ADSIM_VALIDATION_DIR (CMake sets it to <build>/validation for ctest) or,
//    if unset, ./Results/validation.
//  * scripts/docs/gen_validation_page.py merges those files into
//    adsim/test/papers/validation_summary.json and generates
//    docs/validation.md.
//
//  Tolerances are set in the tests from the MEASURED error (rounded up
//  modestly) with a comment giving the measured value — they are regression
//  bounds that document how well the model reproduces the paper, not
//  acceptance criteria chosen in advance.
// ============================================================================

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/solver/coupled_solver.h"
#include "adsim/solver/log.h"

namespace adsim::validation {

  // ---------------------------------------------------------------------------
  // Small utilities
  // ---------------------------------------------------------------------------

  struct Series {
    std::vector<double> t;
    std::vector<double> y;
    bool empty() const { return t.empty(); }
  };

  // Sort by time and drop duplicate / non-increasing time stamps (keeping the
  // last sample at a repeated time). Digitised curves and step-boundary log
  // samples can both contain repeated abscissae.
  inline Series MakeStrictlyIncreasing(Series s) {
    std::vector<size_t> idx(s.t.size());
    for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
    std::stable_sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return s.t[a] < s.t[b]; });
    Series out;
    for (size_t k : idx) {
      if (!std::isfinite(s.t[k]) || !std::isfinite(s.y[k])) continue;
      if (!out.t.empty() && !(s.t[k] > out.t.back())) {
        out.y.back() = s.y[k];
        continue;
      }
      out.t.push_back(s.t[k]);
      out.y.push_back(s.y[k]);
    }
    return out;
  }

  inline double Interp(const Series& s, double x) {
    if (s.t.empty()) return std::numeric_limits<double>::quiet_NaN();
    if (x <= s.t.front()) return s.y.front();
    if (x >= s.t.back()) return s.y.back();
    auto it = std::upper_bound(s.t.begin(), s.t.end(), x);
    size_t i1 = static_cast<size_t>(std::distance(s.t.begin(), it));
    size_t i0 = i1 - 1;
    double w = (x - s.t[i0]) / (s.t[i1] - s.t[i0]);
    return s.y[i0] + w * (s.y[i1] - s.y[i0]);
  }

  // Root-mean-square difference between a simulated and an experimental curve,
  // evaluated at the experimental abscissae that fall inside the simulated time
  // window (simulation linearly interpolated). Returns NaN if there is no overlap.
  inline double CurveRMSE(const Series& sim_in, const Series& exp_in, size_t* n_points = nullptr) {
    const Series sim = MakeStrictlyIncreasing(sim_in);
    const Series exp = MakeStrictlyIncreasing(exp_in);
    if (sim.t.size() < 2 || exp.t.empty()) return std::numeric_limits<double>::quiet_NaN();
    double sum = 0.0;
    size_t n = 0;
    for (size_t i = 0; i < exp.t.size(); ++i) {
      if (exp.t[i] < sim.t.front() || exp.t[i] > sim.t.back()) continue;
      double e = Interp(sim, exp.t[i]) - exp.y[i];
      sum += e * e;
      ++n;
    }
    if (n_points) *n_points = n;
    if (n == 0) return std::numeric_limits<double>::quiet_NaN();
    return std::sqrt(sum / static_cast<double>(n));
  }

  // First time at which the curve rises to `threshold` (linear interpolation
  // between samples). NaN if it never does.
  inline double CrossingTime(const Series& s_in, double threshold) {
    const Series s = MakeStrictlyIncreasing(s_in);
    for (size_t i = 0; i < s.t.size(); ++i) {
      if (s.y[i] >= threshold) {
        if (i == 0) return s.t[0];
        double w = (threshold - s.y[i - 1]) / (s.y[i] - s.y[i - 1]);
        return s.t[i - 1] + w * (s.t[i] - s.t[i - 1]);
      }
    }
    return std::numeric_limits<double>::quiet_NaN();
  }

  inline double MaxValue(const Series& s) {
    double m = -std::numeric_limits<double>::infinity();
    for (double v : s.y) m = std::max(m, v);
    return m;
  }

  // Walk up from the CWD (ctest runs in adsim/) to find a repo-relative path.
  inline std::filesystem::path ResolveRepoPath(const std::string& relative) {
    namespace fs = std::filesystem;
    for (fs::path dir = fs::current_path();; dir = dir.parent_path()) {
      if (fs::exists(dir / relative)) return dir / relative;
      if (fs::exists(dir / "adsim" / relative)) return dir / "adsim" / relative;
      if (dir == dir.parent_path()) break;
    }
    return relative;
  }

  // Minimal RFC-4180-ish CSV line splitter (double-quoted fields may contain
  // commas; "" inside quotes is a literal quote).
  inline std::vector<std::string> SplitCsvLine(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool quoted = false;
    for (size_t i = 0; i < line.size(); ++i) {
      char c = line[i];
      if (quoted) {
        if (c == '"') {
          if (i + 1 < line.size() && line[i + 1] == '"') { cur += '"'; ++i; }
          else quoted = false;
        } else {
          cur += c;
        }
      } else if (c == '"') {
        quoted = true;
      } else if (c == ',') {
        out.push_back(cur);
        cur.clear();
      } else if (c != '\r') {
        cur += c;
      }
    }
    out.push_back(cur);
    return out;
  }

  inline std::string StripBom(std::string s) {
    if (s.size() >= 3 && static_cast<unsigned char>(s[0]) == 0xEF &&
        static_cast<unsigned char>(s[1]) == 0xBB && static_cast<unsigned char>(s[2]) == 0xBF)
      s.erase(0, 3);
    return s;
  }

  // Load a two-column digitised CSV (header row + "x,y" rows), e.g.
  // adsim/data/streb_2021/A3/A3_mole_fraction_CO2.csv. `relative` is resolved
  // like the experimental-data loader (relative to adsim/).
  inline Series LoadSeriesCsv(const std::string& relative) {
    Series s;
    std::ifstream f(ResolveRepoPath(relative));
    if (!f) {
      ADD_FAILURE() << "Could not open digitised data file " << relative;
      return s;
    }
    std::string line;
    bool header = true;
    while (std::getline(f, line)) {
      line = StripBom(line);
      if (header) { header = false; continue; }
      if (line.empty()) continue;
      auto cols = SplitCsvLine(line);
      if (cols.size() < 2) continue;
      try {
        s.t.push_back(std::stod(cols[0]));
        s.y.push_back(std::stod(cols[1]));
      } catch (...) {
        // skip non-numeric rows
      }
    }
    return MakeStrictlyIncreasing(s);
  }

  // ---------------------------------------------------------------------------
  // Published reference values (adsim/test/papers/<paper>/reference_values.csv)
  //   columns: case,quantity,value,unit,source,note
  // ---------------------------------------------------------------------------

  struct RefValue {
    double value = std::numeric_limits<double>::quiet_NaN();
    std::string unit;
    std::string source;
    std::string note;
  };

  inline RefValue Reference(const std::string& paper, const std::string& case_name,
                            const std::string& quantity) {
    const std::string rel = "test/papers/" + paper + "/reference_values.csv";
    std::ifstream f(ResolveRepoPath(rel));
    if (!f) {
      ADD_FAILURE() << "Missing reference table " << rel;
      return {};
    }
    std::string line;
    bool header = true;
    while (std::getline(f, line)) {
      line = StripBom(line);
      if (header) { header = false; continue; }
      if (line.empty() || line[0] == '#') continue;
      auto c = SplitCsvLine(line);
      if (c.size() < 5) continue;
      if (c[0] == case_name && c[1] == quantity) {
        RefValue r;
        r.value = std::stod(c[2]);
        r.unit = c[3];
        r.source = c[4];
        if (c.size() > 5) r.note = c[5];
        return r;
      }
    }
    ADD_FAILURE() << "No reference value '" << quantity << "' for case '" << case_name
                  << "' in " << rel;
    return {};
  }

  // ---------------------------------------------------------------------------
  // Simulation series extraction
  // ---------------------------------------------------------------------------

  // Outlet (or inlet) face mole fraction of `component` during one step.
  inline Series FaceMoleFraction(const CoupledSolver& solver, const std::string& component,
                                 int cycle_index, int step_index, bool top = true) {
    const auto e = solver.GetLog().GetComponentFaceVariable(
        cycle_index, step_index, ComponentVariableType::MOLE_FRACTION, component);
    Series s;
    for (size_t i = 0; i < e.times.size() && i < e.values.size(); ++i) {
      if (e.values[i].empty()) continue;
      s.t.push_back(e.times[i]);
      s.y.push_back(top ? e.values[i].back() : e.values[i].front());
    }
    return MakeStrictlyIncreasing(s);
  }

  // A cell variable (e.g. gas temperature) sampled at axial position z [m] of a
  // single-layer bed of length L with n_cells cells, linearly interpolated
  // between cell centres. Handles logs that include one ghost cell at each end.
  inline Series CellVariableAt(const CoupledSolver& solver, VariableType type, double z,
                               double L, int n_cells, int cycle_index = 0, int step_index = 0) {
    const auto e = solver.GetLog().GetVariable(cycle_index, step_index, type);
    Series s;
    if (e.values.empty()) return s;
    const int n_total = static_cast<int>(e.values.front().size());
    const int offset = (n_total == n_cells + 2) ? 1 : 0;
    const double dz = L / n_cells;
    double xi = z / dz - 0.5;  // fractional cell-centre index
    xi = std::clamp(xi, 0.0, static_cast<double>(n_cells - 1));
    const int i0 = static_cast<int>(std::floor(xi));
    const int i1 = std::min(i0 + 1, n_cells - 1);
    const double w = xi - i0;
    for (size_t k = 0; k < e.times.size() && k < e.values.size(); ++k) {
      const auto& v = e.values[k];
      if (static_cast<int>(v.size()) < offset + n_cells) continue;
      s.t.push_back(e.times[k]);
      s.y.push_back((1.0 - w) * v[offset + i0] + w * v[offset + i1]);
    }
    return MakeStrictlyIncreasing(s);
  }

  // ---------------------------------------------------------------------------
  // Per-case report
  // ---------------------------------------------------------------------------

  struct Metric {
    std::string quantity;
    double simulated = 0.0;
    double reference = 0.0;
    double error = 0.0;
    std::string metric;     // e.g. "abs_error", "rel_error", "rmse", "nrmse"
    double tolerance = std::numeric_limits<double>::quiet_NaN();  // NaN => informational only
    bool pass = true;
    std::string units;      // units of simulated/reference
    std::string error_units;
    std::string source;     // table / figure / page of the reference
    std::string note;
  };

  inline std::string JsonEscape(const std::string& s) {
    std::string o;
    for (char c : s) {
      switch (c) {
        case '"': o += "\\\""; break;
        case '\\': o += "\\\\"; break;
        case '\n': o += "\\n"; break;
        case '\t': o += "\\t"; break;
        case '\r': break;
        default:
          if (static_cast<unsigned char>(c) < 0x20) o += ' ';
          else o += c;
      }
    }
    return o;
  }

  inline std::string JsonNumber(double v) {
    if (!std::isfinite(v)) return "null";
    std::ostringstream ss;
    ss.precision(10);
    ss << v;
    return ss.str();
  }

  inline std::filesystem::path OutputDirectory() {
    if (const char* env = std::getenv("ADSIM_VALIDATION_DIR"); env && *env)
      return std::filesystem::path(env);
    return std::filesystem::current_path() / "Results" / "validation";
  }

  class CaseReport {
  public:
    CaseReport(std::string paper, std::string case_name)
      : paper_(std::move(paper)), case_(std::move(case_name)) {}

    ~CaseReport() { Write(); }

    CaseReport(const CaseReport&) = delete;
    CaseReport& operator=(const CaseReport&) = delete;

    const std::string& paper() const { return paper_; }
    const std::string& case_name() const { return case_; }

    // |sim - ref| <= tol (same units as the quantity, e.g. percentage points).
    void CheckAbs(const std::string& quantity, double sim, const RefValue& ref, double tol,
                  const std::string& note = "") {
      Metric m = Base(quantity, sim, ref, note);
      m.metric = "abs_error";
      m.error = std::abs(sim - ref.value);
      m.error_units = ref.unit;
      Finish(m, tol);
    }

    // |sim - ref| / |ref| <= tol (fraction, reported in %).
    void CheckRel(const std::string& quantity, double sim, const RefValue& ref, double tol,
                  const std::string& note = "") {
      Metric m = Base(quantity, sim, ref, note);
      m.metric = "rel_error";
      m.error = std::abs(sim - ref.value) / std::abs(ref.value);
      m.error_units = "-";
      Finish(m, tol);
    }

    // Curve RMSE (absolute, in the curve's units) against a digitised series.
    // `scale` > 0 turns it into a normalised RMSE (RMSE / scale), e.g. the feed
    // mole fraction of that component.
    double CheckCurve(const std::string& quantity, const Series& sim, const Series& exp,
                      double tol, const std::string& units, const std::string& source,
                      double scale = 0.0, const std::string& note = "") {
      size_t n = 0;
      const double rmse = CurveRMSE(sim, exp, &n);
      Metric m;
      m.quantity = quantity;
      m.units = units;
      m.source = source;
      m.note = note.empty() ? ("RMSE over " + std::to_string(n) + " digitised points")
                            : note + "; RMSE over " + std::to_string(n) + " digitised points";
      if (scale > 0.0) {
        m.metric = "nrmse";
        m.error = rmse / scale;
        m.simulated = rmse;
        m.reference = scale;
        m.error_units = "-";
        m.note += "; normalised by " + JsonNumber(scale) + " " + units;
      } else {
        m.metric = "rmse";
        m.error = rmse;
        m.simulated = rmse;
        m.reference = 0.0;
        m.error_units = units;
      }
      EXPECT_GT(n, 0u) << paper_ << "/" << case_ << ": no overlap between simulation and data for " << quantity;
      Finish(m, tol);
      return m.error;
    }

    // Relative error of the time at which the curve first reaches `level`
    // (e.g. half the feed mole fraction): |t_sim - t_exp| / t_exp, with t_exp
    // interpolated from the digitised points. A front-position metric that is
    // less dominated by the steepness of the front than the curve RMSE.
    double CheckCrossingTime(const std::string& quantity, const Series& sim, const Series& exp,
                             double level, double tol, const std::string& source) {
      const double ts = CrossingTime(sim, level);
      const double te = CrossingTime(exp, level);
      Metric m;
      m.quantity = quantity;
      m.simulated = ts;
      m.reference = te;
      m.units = "s";
      m.metric = "rel_error";
      m.error = std::abs(ts - te) / te;
      m.error_units = "-";
      m.source = source;
      m.note = "first time the curve reaches " + JsonNumber(level) +
               " (reference time interpolated from the digitised points)";
      Finish(m, tol);
      return m.error;
    }

    // Assert on an error computed by the caller (e.g. an RMSE pooled over
    // several thermocouples). `simulated`/`reference` are stored as given.
    void CheckError(const std::string& quantity, double error, const std::string& metric, double tol,
                    const std::string& error_units, const std::string& source, const std::string& note,
                    double simulated = std::numeric_limits<double>::quiet_NaN(),
                    double reference = std::numeric_limits<double>::quiet_NaN()) {
      Metric m;
      m.quantity = quantity;
      m.error = error;
      m.metric = metric;
      m.units = error_units;
      m.error_units = error_units;
      m.source = source;
      m.note = note;
      m.simulated = std::isfinite(simulated) ? simulated : error;
      m.reference = std::isfinite(reference) ? reference : 0.0;
      Finish(m, tol);
    }

    // Record a comparison without asserting on it (e.g. a reference whose
    // definition could not be confirmed from the paper, or a quantity the deck
    // cannot compute). Still shown in the generated validation page.
    void Info(const std::string& quantity, double sim, const RefValue& ref, const std::string& metric,
              const std::string& note) {
      Metric m = Base(quantity, sim, ref, note);
      m.metric = metric;
      if (metric == "rel_error") {
        m.error = std::abs(sim - ref.value) / std::abs(ref.value);
        m.error_units = "-";
      } else {
        m.error = std::abs(sim - ref.value);
        m.error_units = ref.unit;
      }
      m.tolerance = std::numeric_limits<double>::quiet_NaN();
      m.pass = true;
      Log(m, false);
      metrics_.push_back(m);
    }

    void Write() {
      if (written_ || metrics_.empty()) return;
      written_ = true;
      namespace fs = std::filesystem;
      std::error_code ec;
      const fs::path dir = OutputDirectory();
      fs::create_directories(dir, ec);
      std::string safe_case = case_;
      for (char& c : safe_case)
        if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-')) c = '_';
      std::ofstream f(dir / (paper_ + "__" + safe_case + ".json"));
      if (!f) return;
      f << "{\n  \"paper\": \"" << JsonEscape(paper_) << "\",\n  \"case\": \"" << JsonEscape(case_)
        << "\",\n  \"metrics\": [\n";
      for (size_t i = 0; i < metrics_.size(); ++i) {
        const auto& m = metrics_[i];
        f << "    {\"quantity\": \"" << JsonEscape(m.quantity) << "\", \"simulated\": " << JsonNumber(m.simulated)
          << ", \"reference\": " << JsonNumber(m.reference) << ", \"error\": " << JsonNumber(m.error)
          << ", \"metric\": \"" << JsonEscape(m.metric) << "\", \"tolerance\": " << JsonNumber(m.tolerance)
          << ", \"pass\": " << (m.pass ? "true" : "false") << ", \"units\": \"" << JsonEscape(m.units)
          << "\", \"error_units\": \"" << JsonEscape(m.error_units) << "\", \"source\": \""
          << JsonEscape(m.source) << "\", \"note\": \"" << JsonEscape(m.note) << "\"}"
          << (i + 1 < metrics_.size() ? ",\n" : "\n");
      }
      f << "  ]\n}\n";
    }

  private:
    Metric Base(const std::string& quantity, double sim, const RefValue& ref, const std::string& note) const {
      Metric m;
      m.quantity = quantity;
      m.simulated = sim;
      m.reference = ref.value;
      m.units = ref.unit;
      m.source = ref.source;
      m.note = note.empty() ? ref.note : (ref.note.empty() ? note : ref.note + "; " + note);
      return m;
    }

    void Finish(Metric& m, double tol) {
      m.tolerance = tol;
      m.pass = std::isfinite(m.error) && m.error <= tol;
      Log(m, true);
      metrics_.push_back(m);
      EXPECT_TRUE(m.pass) << "[validation] " << paper_ << " / " << case_ << " / " << m.quantity << ": "
                          << m.metric << " = " << m.error << " exceeds tolerance " << tol
                          << " (simulated " << m.simulated << ", reference " << m.reference << ")";
    }

    void Log(const Metric& m, bool asserted) const {
      std::cout << "[validation] " << paper_ << " | " << case_ << " | " << m.quantity << " | sim="
                << m.simulated << " ref=" << m.reference << " | " << m.metric << "=" << m.error;
      if (asserted) std::cout << " (tol " << m.tolerance << (m.pass ? ", pass" : ", FAIL") << ")";
      else std::cout << " (info)";
      std::cout << std::endl;
    }

    std::string paper_;
    std::string case_;
    std::vector<Metric> metrics_;
    bool written_ = false;
  };

} // namespace adsim::validation
