// Python bindings for the Skarstrom engine (module skarstrom._core).
//
// A thin layer over adsim::Simulation (adsim/simulation.h) and the isotherm
// fitter (adsim/solver/isotherm_fit.h). Series and profiles come
// back as NumPy arrays. The pure-Python package (skarstrom/__init__.py) wraps
// these classes with a friendlier surface; nothing here is part of the public
// API except through it.

#include <algorithm>
#include <atomic>
#include <map>
#include <chrono>
#include <cmath>
#include <type_traits>
#include <optional>
#include <string>
#include <thread>
#include <vector>

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "adsim/simulation.h"
#include "adsim/version.h"

#include "adsim/equations/isotherm_kinds.h"
#include "adsim/equations/iast.h"
#include "adsim/math/unit_parser.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/deck_edit.h"
#include "adsim/objects/deck_view.h"
#include "adsim/objects/manifests.h"
#include "adsim/objects/model_kinds.h"
#include "adsim/solver/calibration.h"
#include "adsim/solver/cyclic_steady_state.h"
#include "adsim/solver/isotherm_fit.h"
#include "adsim/solver/optimisation.h"
#include "adsim/solver/parameterisation.h"
#include "adsim/solver/sensitivity.h"
#include "adsim/solver/state_cache.h"

namespace py = pybind11;

namespace {

py::array_t<double> ToArray(const std::vector<double>& v) {
  py::array_t<double> a(static_cast<py::ssize_t>(v.size()));
  std::copy(v.begin(), v.end(), a.mutable_data());
  return a;
}

py::array_t<double> ToArray2D(const std::vector<double>& v, size_t rows, size_t cols) {
  py::array_t<double> a({static_cast<py::ssize_t>(rows), static_cast<py::ssize_t>(cols)});
  std::copy(v.begin(), v.end(), a.mutable_data());
  return a;
}

adsim::LogLevel ParseLogLevel(const std::string& s) {
  if (s == "low") return adsim::LogLevel::LOW;
  if (s == "medium") return adsim::LogLevel::MEDIUM;
  if (s == "high") return adsim::LogLevel::HIGH;
  throw py::value_error("log_level must be 'low', 'medium' or 'high', got '" + s + "'");
}

const char* LogLevelName(adsim::LogLevel l) {
  switch (l) {
    case adsim::LogLevel::LOW: return "low";
    case adsim::LogLevel::MEDIUM: return "medium";
    case adsim::LogLevel::HIGH: return "high";
  }
  return "low";
}

adsim::SeriesLocation ParseLocation(const std::string& s) {
  if (s == "top") return adsim::SeriesLocation::Top;
  if (s == "bottom") return adsim::SeriesLocation::Bottom;
  if (s == "position") return adsim::SeriesLocation::Position;
  throw py::value_error("location must be 'top', 'bottom' or 'position', got '" + s + "'");
}

adsim::IsothermType ParseIsothermKind(const std::string& kind) {
  adsim::IsothermType t{};
  if (!adsim::KindFromString(kind, t))
    throw py::value_error("unknown isotherm kind '" + kind + "' (see skarstrom.isotherm_kinds())");
  return t;
}

std::vector<adsim::isofit::Point> MakePoints(const std::vector<double>& P, const std::vector<double>& q,
                                                  const std::vector<double>& T, const std::vector<double>& y,
                                                  const std::vector<double>& q_w = {}) {
  const size_t n = P.size();
  if (!q.empty() && q.size() != n) throw py::value_error("pressure and loading must have the same length");
  if (T.size() != n && T.size() != 1) throw py::value_error("temperature must be a scalar or match pressure");
  if (!y.empty() && y.size() != n && y.size() != 1)
    throw py::value_error("mole_fraction must be a scalar or match pressure");
  if (!q_w.empty() && q_w.size() != n && q_w.size() != 1)
    throw py::value_error("water_loading must be a scalar or match pressure");
  std::vector<adsim::isofit::Point> pts(n);
  for (size_t i = 0; i < n; ++i) {
    pts[i].P = P[i];
    pts[i].q = q.empty() ? 0.0 : q[i];
    pts[i].T = T.size() == 1 ? T[0] : T[i];
    pts[i].y = y.empty() ? 1.0 : (y.size() == 1 ? y[0] : y[i]);
    pts[i].q_w = q_w.empty() ? 0.0 : (q_w.size() == 1 ? q_w[0] : q_w[i]);
  }
  return pts;
}

py::dict SpecToDict(const adsim::isofit::ParamSpec& s) {
  py::dict d;
  d["name"] = s.name;
  d["value"] = s.value;
  d["lower"] = s.lower;
  d["upper"] = s.upper;
  d["fixed"] = s.fixed;
  d["log_scale"] = s.scale == adsim::isofit::Scale::Log;
  d["paired_energy"] = s.pair_energy;
  d["paired_energy_sign"] = s.pair_sign;
  d["note"] = s.note;
  return d;
}

adsim::isofit::ParamSpec SpecFromDict(const py::dict& d) {
  adsim::isofit::ParamSpec s;
  s.name = d["name"].cast<std::string>();
  s.value = d["value"].cast<double>();
  s.lower = d["lower"].cast<double>();
  s.upper = d["upper"].cast<double>();
  s.fixed = d.contains("fixed") ? d["fixed"].cast<bool>() : false;
  s.scale = (d.contains("log_scale") && d["log_scale"].cast<bool>()) ? adsim::isofit::Scale::Log
                                                                    : adsim::isofit::Scale::Linear;
  if (d.contains("paired_energy")) s.pair_energy = d["paired_energy"].cast<std::string>();
  if (d.contains("paired_energy_sign")) s.pair_sign = d["paired_energy_sign"].cast<double>();
  return s;
}


// ---- sensitivity (adsim/solver/sensitivity.h) ------------------------------------
// Inputs, plans and policies cross as plain dicts; skarstrom/sensitivity.py
// builds the public API on top.

namespace sa = adsim::sensitivity;

double DictNumber(const py::dict& d, const char* key) {
  if (!d.contains(key) || d[key].is_none()) return sa::kNaN;
  return d[key].cast<double>();
}

sa::Input InputFromDict(const py::dict& d) {
  sa::Input in;
  in.path = d["path"].cast<std::string>();
  in.distribution = sa::DistributionFromName(d.contains("distribution") ? d["distribution"].cast<std::string>()
                                                                        : std::string("uniform"));
  in.lower = DictNumber(d, "lower");
  in.upper = DictNumber(d, "upper");
  in.mode = DictNumber(d, "mode");
  in.mean = DictNumber(d, "mean");
  in.sd = DictNumber(d, "sd");
  return in;
}

std::vector<sa::Input> InputsFromList(const py::list& l) {
  std::vector<sa::Input> out;
  for (auto item : l) out.push_back(InputFromDict(item.cast<py::dict>()));
  return out;
}

std::vector<sa::Output> OutputsFromList(const std::vector<std::string>& names) {
  std::vector<sa::Output> out;
  for (const auto& n : names) out.push_back(sa::ParseOutput(n));
  return out;
}

sa::RunPolicy PolicyFromDict(const py::dict& d) {
  sa::RunPolicy p;
  if (d.contains("require_css")) p.require_css = d["require_css"].cast<bool>();
  if (d.contains("cycles") && !d["cycles"].is_none()) p.cycles = d["cycles"].cast<int>();
  if (d.contains("settings")) {
    for (auto item : d["settings"].cast<py::dict>()) {
      p.settings.emplace_back(item.first.cast<std::string>(), py::str(item.second).cast<std::string>());
    }
  }
  return p;
}

py::dict PlanToDict(const sa::Plan& plan) {
  py::dict d;
  d["method"] = sa::MethodName(plan.method);
  d["seed"] = plan.seed;
  d["trajectories"] = plan.trajectories;
  d["levels"] = plan.levels;
  d["batch"] = plan.batch;
  d["changed"] = plan.changed;
  d["u"] = plan.u;
  d["values"] = plan.values;
  return d;
}

py::object NumberOrNone(double v) { return std::isnan(v) ? py::object(py::none()) : py::object(py::float_(v)); }

py::dict InputToDict(const sa::Input& in) {
  py::dict d;
  d["path"] = in.path;
  d["distribution"] = sa::DistributionName(in.distribution);
  d["lower"] = NumberOrNone(in.lower);
  d["upper"] = NumberOrNone(in.upper);
  d["mode"] = NumberOrNone(in.mode);
  d["mean"] = NumberOrNone(in.mean);
  d["sd"] = NumberOrNone(in.sd);
  return d;
}

// The study saved with a deck (UserData::sensitivity_study) as a dict:
// {inputs, outputs, policy: {require_css, cycles (None: the deck's), settings},
//  seed, trajectories, levels, threshold, runs, batch, max_runs, rel_tol, bootstrap}.
py::dict StudyToDict(const sa::Study& st) {
  py::dict d;
  py::list inputs;
  for (const auto& in : st.inputs) inputs.append(InputToDict(in));
  d["inputs"] = inputs;
  d["outputs"] = st.outputs;
  py::dict policy, settings;
  policy["require_css"] = st.policy.require_css;
  policy["cycles"] = st.policy.cycles > 0 ? py::object(py::int_(st.policy.cycles)) : py::object(py::none());
  for (const auto& [path, text] : st.policy.settings) settings[py::str(path)] = text;
  policy["settings"] = settings;
  d["policy"] = policy;
  d["seed"] = st.seed;
  d["trajectories"] = st.trajectories;
  d["levels"] = st.levels;
  d["threshold"] = st.threshold;
  d["runs"] = st.runs;
  d["batch"] = st.batch;
  d["max_runs"] = st.max_runs;
  d["rel_tol"] = st.rel_tol;
  d["bootstrap"] = st.bootstrap;
  return d;
}

// The inverse of StudyToDict; a missing key keeps the default.
sa::Study StudyFromDict(const py::dict& d) {
  sa::Study st;
  if (d.contains("inputs")) st.inputs = InputsFromList(d["inputs"].cast<py::list>());
  if (d.contains("outputs") && !d["outputs"].is_none()) st.outputs = d["outputs"].cast<std::vector<std::string>>();
  if (d.contains("policy") && !d["policy"].is_none()) st.policy = PolicyFromDict(d["policy"].cast<py::dict>());
  auto get = [&d](const char* key, auto& out) {
    if (d.contains(key) && !d[key].is_none()) out = d[key].cast<std::decay_t<decltype(out)>>();
  };
  get("seed", st.seed);
  get("trajectories", st.trajectories);
  get("levels", st.levels);
  get("threshold", st.threshold);
  get("runs", st.runs);
  get("batch", st.batch);
  get("max_runs", st.max_runs);
  get("rel_tol", st.rel_tol);
  get("bootstrap", st.bootstrap);
  return st;
}

py::dict IntervalToDict(const sa::Interval& i) {
  py::dict d;
  d["low"] = i.low;
  d["high"] = i.high;
  return d;
}

// ---- optimisation (adsim/solver/optimisation.h) ------------------------------

adsim::KPIType KPITypeFromName(const std::string& name) {
  if (name == "purity") return adsim::KPIType::PURITY;
  if (name == "recovery") return adsim::KPIType::RECOVERY;
  if (name == "productivity_mass") return adsim::KPIType::PRODUCTIVITY_MASS;
  if (name == "productivity_volume") return adsim::KPIType::PRODUCTIVITY_VOLUME;
  if (name == "energy_consumption" || name == "electrical_energy") return adsim::KPIType::ENERGY_CONSUMPTION;
  if (name == "production_rate") return adsim::KPIType::PRODUCTION_RATE;
  if (name == "thermal_energy") return adsim::KPIType::THERMAL_ENERGY;
  if (name == "cooling_duty") return adsim::KPIType::COOLING_DUTY;
  if (name == "end_pressure_bottom") return adsim::KPIType::END_PRESSURE_BOTTOM;
  if (name == "end_pressure_top") return adsim::KPIType::END_PRESSURE_TOP;
  if (name == "end_pressure_min") return adsim::KPIType::END_PRESSURE_MIN;
  if (name == "end_pressure_mean") return adsim::KPIType::END_PRESSURE_MEAN;
  if (name == "end_pressure_max") return adsim::KPIType::END_PRESSURE_MAX;
  throw std::invalid_argument("unknown KPI '" + name +
                              "' (purity, recovery, productivity_mass, productivity_volume, energy_consumption, "
                              "production_rate, thermal_energy, cooling_duty; for a step[...]: end_pressure_bottom, "
                              "end_pressure_top, end_pressure_min, end_pressure_mean, end_pressure_max)");
}

const char* KPITypeKey(adsim::KPIType t) {
  switch (t) {
    case adsim::KPIType::PURITY: return "purity";
    case adsim::KPIType::RECOVERY: return "recovery";
    case adsim::KPIType::PRODUCTIVITY_MASS: return "productivity_mass";
    case adsim::KPIType::PRODUCTIVITY_VOLUME: return "productivity_volume";
    case adsim::KPIType::ENERGY_CONSUMPTION: return "energy_consumption";
    case adsim::KPIType::PRODUCTION_RATE: return "production_rate";
    case adsim::KPIType::THERMAL_ENERGY: return "thermal_energy";
    case adsim::KPIType::COOLING_DUTY: return "cooling_duty";
    case adsim::KPIType::LCOC: return "lcoc";
    case adsim::KPIType::END_PRESSURE_BOTTOM: return "end_pressure_bottom";
    case adsim::KPIType::END_PRESSURE_TOP: return "end_pressure_top";
    case adsim::KPIType::END_PRESSURE_MIN: return "end_pressure_min";
    case adsim::KPIType::END_PRESSURE_MEAN: return "end_pressure_mean";
    case adsim::KPIType::END_PRESSURE_MAX: return "end_pressure_max";
  }
  return "unknown";
}

py::dict CachedRunInfoToDict(const adsim::CachedRunInfo& i) {
  py::dict d;
  d["warm"] = i.warm;
  d["distance"] = i.distance;
  d["start_cycles"] = i.start_cycles;
  d["skipped"] = i.skipped;
  d["retried_cold"] = i.retried_cold;
  d["cached"] = i.cached;
  d["note"] = i.note;
  return d;
}

py::dict EvaluationToDict(const adsim::OptimizationEvaluation& e) {
  py::dict d;
  d["index"] = e.index;
  d["x"] = e.x;
  d["cost"] = e.cost;
  py::dict kpis;
  for (const auto& [component, values] : e.kpis) {
    py::dict v;
    for (const auto& [type, value] : values) v[KPITypeKey(type)] = value;
    kpis[py::str(component)] = v;
  }
  d["kpis"] = kpis;
  d["status"] = adsim::RunStatusName(e.status);
  d["message"] = e.message;
  d["accepted"] = e.accepted;
  d["warm"] = e.warm;
  d["distance"] = e.distance;
  d["note"] = e.note;
  d["css_reached"] = e.css_reached;
  d["cycles_run"] = e.cycles_run;
  d["wall_time_s"] = e.wall_time_s;
  return d;
}

}  // namespace

// A blocking run from Python: the GIL is released while the solver works and
// callbacks re-acquire it. Ctrl-C (or any pending Python signal) cancels the
// run and is re-raised afterwards. `go` runs the simulation with the options.
adsim::RunResult RunFromPython(std::optional<int> cycles, py::object progress, const std::string& log_level,
                               py::object on_message,
                               const std::function<adsim::RunResult(const adsim::RunOptions&)>& go) {
  adsim::RunOptions opt;
  opt.cycles = cycles;
  opt.log_level = ParseLogLevel(log_level);
  bool python_error = false;
  std::optional<py::error_already_set> pending;
  if (!progress.is_none()) {
    opt.progress = [&](const adsim::RunProgress& p) {
      py::gil_scoped_acquire gil;
      if (pending) return false;
      try {
        py::dict d;
        d["cycle"] = p.cycle;
        d["step"] = p.step;
        d["total_cycles"] = p.total_cycles;
        d["total_steps"] = p.total_steps;
        d["fraction"] = p.fraction;
        py::object r = progress(d);
        return r.is_none() || r.cast<bool>();
      } catch (py::error_already_set& e) {
        pending.emplace(std::move(e));
        python_error = true;
        return false;
      }
    };
  }
  if (!on_message.is_none()) {
    opt.on_message = [&](adsim::LogLevel level, const std::string& msg) {
      py::gil_scoped_acquire gil;
      if (pending) return;
      try {
        on_message(LogLevelName(level), msg);
      } catch (py::error_already_set& e) {
        pending.emplace(std::move(e));
        python_error = true;
      }
    };
  }
  auto last_check = std::chrono::steady_clock::now();
  opt.cancel = [&]() {
    if (python_error) return true;
    const auto now = std::chrono::steady_clock::now();
    if (now - last_check < std::chrono::milliseconds(100)) return false;
    last_check = now;
    py::gil_scoped_acquire gil;
    if (PyErr_CheckSignals() != 0) {
      pending.emplace();  // fetches the KeyboardInterrupt (or other) error
      python_error = true;
      return true;
    }
    return false;
  };
  adsim::RunResult r;
  {
    py::gil_scoped_release release;
    r = go(opt);
  }
  if (pending) throw std::move(*pending);
  return r;
}

namespace {

// ---- experimental datasets (adsim/solver/experimental_data.h) ----------------

size_t TagIndex(const std::vector<std::string>& tags, const std::string& tag, const char* what) {
  for (size_t i = 0; i < tags.size(); ++i)
    if (tags[i] == tag) return i;
  std::string valid;
  for (const auto& t : tags) valid += (valid.empty() ? "" : ", ") + t;
  throw py::value_error(std::string("unknown ") + what + " '" + tag + "'; one of: " + valid);
}

py::dict BindingToDict(const adsim::ObservableBinding& b) {
  py::dict d;
  d["kind"] = adsim::ObservableKindTags()[static_cast<size_t>(b.kind)];
  d["variable"] = b.variable;
  d["component"] = b.component;
  d["location"] = adsim::ObservableLocationTags()[static_cast<size_t>(b.location)];
  d["z"] = b.z;
  d["cycle"] = b.cycle;
  d["cycle_offset"] = b.cycle_offset;
  d["x_axis"] = b.XAxis() == adsim::ObservableXAxis::Cycle ? "cycle" : "time";
  d["x_unit"] = b.x_unit;
  d["x_scale"] = b.x_scale;
  d["y_unit"] = b.y_unit;
  d["si_unit"] = adsim::ObservableSIUnit(b);
  d["uncertainty"] = b.uncertainty.empty() ? py::object(py::none()) : py::object(ToArray(b.uncertainty));
  d["weight"] = b.weight;
  return d;
}

// Missing (or None) keys keep the defaults; unknown keys are refused.
adsim::ObservableBinding BindingFromDict(const py::dict& d) {
  static const char* known[] = {"kind", "variable", "component", "location", "z", "cycle", "cycle_offset",
                                "x_axis", "x_unit", "x_scale", "y_unit", "si_unit", "uncertainty", "weight"};
  for (const auto& item : d) {
    const std::string key = py::str(item.first);
    if (std::find(std::begin(known), std::end(known), key) == std::end(known))
      throw py::value_error("unknown observable field '" + key + "'");
  }
  adsim::ObservableBinding b;
  auto get = [&](const char* key) -> py::object {
    return d.contains(key) && !d[key].is_none() ? py::object(d[key]) : py::object(py::none());
  };
  if (auto v = get("kind"); !v.is_none())
    b.kind = static_cast<adsim::ObservableKind>(
        TagIndex(adsim::ObservableKindTags(), v.cast<std::string>(), "observable kind"));
  if (auto v = get("variable"); !v.is_none()) b.variable = v.cast<std::string>();
  if (auto v = get("component"); !v.is_none()) b.component = v.cast<std::string>();
  if (auto v = get("location"); !v.is_none())
    b.location = static_cast<adsim::ObservableLocation>(
        TagIndex(adsim::ObservableLocationTags(), v.cast<std::string>(), "observable location"));
  if (auto v = get("z"); !v.is_none()) b.z = v.cast<double>();
  if (auto v = get("cycle"); !v.is_none()) b.cycle = v.cast<int>();
  if (auto v = get("cycle_offset"); !v.is_none()) b.cycle_offset = v.cast<int>();
  if (auto v = get("x_unit"); !v.is_none()) b.x_unit = v.cast<std::string>();
  if (auto v = get("x_scale"); !v.is_none()) b.x_scale = v.cast<double>();
  if (auto v = get("y_unit"); !v.is_none()) b.y_unit = v.cast<std::string>();
  if (auto v = get("uncertainty"); !v.is_none()) b.uncertainty = v.cast<std::vector<double>>();
  if (auto v = get("weight"); !v.is_none()) b.weight = v.cast<double>();
  return b;
}

py::dict DatasetToDict(const adsim::XYData& xy) {
  py::dict d;
  d["name"] = xy.name;
  d["x"] = ToArray(xy.x);
  d["y"] = ToArray(xy.y);
  d["x_label"] = xy.x_label;
  d["y_label"] = xy.y_label;
  d["observable"] = BindingToDict(xy.observable);
  return d;
}

// A dataset from its points and an observable dict, checked against the deck's
// components: the points need not be stored in the deck.
adsim::XYData DatasetFromArgs(const adsim::Simulation& s, const std::string& name, std::vector<double> x,
                              std::vector<double> y, const py::dict& observable) {
  adsim::XYData xy;
  xy.name = name;
  xy.x = std::move(x);
  xy.y = std::move(y);
  if (xy.x.size() != xy.y.size()) throw py::value_error("x and y differ in length");
  xy.observable = BindingFromDict(observable);
  std::string err;
  if (!adsim::NormaliseObservableBinding(xy.observable, xy, s.Data().fluid.components_, err))
    throw py::value_error(err);
  return xy;
}

// The simulated counterpart of a bound dataset from the last run.
adsim::SimulatedCounterpart Simulated(const adsim::Simulation& s, const adsim::XYData& xy) {
  if (!xy.observable.Bound()) throw py::value_error("dataset '" + xy.name + "' is not bound to an observable");
  std::vector<adsim::KPIRecord> kpis;
  if (xy.observable.kind == adsim::ObservableKind::Kpi) kpis = s.KPIs();
  const adsim::BedBoundarySources sources{&s.Boundaries(), &s.Data().cycle.step_list_};
  return adsim::SimulatedSeries(xy, s.Log(), kpis, sources);
}

}  // namespace

// ---- calibration (adsim/solver/calibration.h) --------------------------------

namespace cal = adsim::calibration;

// A flag Python can raise from any thread (a progress callback, another
// thread, a GUI) to stop a calibration: the fit returns its best point.
struct CalibrationCancel {
  std::atomic<bool> flag{false};
  void Cancel() { flag = true; }
  bool Cancelled() const { return flag.load(); }
  void Reset() { flag = false; }
};

namespace {

double OptionalDouble(const py::dict& d, const char* key) {
  if (!d.contains(key) || d[key].is_none()) return cal::kNaN;
  return d[key].cast<double>();
}

cal::Scale ScaleFromName(const std::string& s) {
  if (s == "linear") return cal::Scale::Linear;
  if (s == "log") return cal::Scale::Log;
  throw std::invalid_argument("a parameter's scale is 'linear' or 'log'; got '" + s + "'");
}

cal::Case CaseFromDict(const py::dict& d) {
  cal::Case c;
  if (d.contains("name") && !d["name"].is_none()) c.name = d["name"].cast<std::string>();
  c.deck_json = d["deck_json"].cast<std::string>();
  if (d.contains("datasets")) {
    for (auto item : d["datasets"].cast<py::list>()) {
      auto t = item.cast<py::tuple>();
      const double w = t[1].is_none() ? cal::kNaN : t[1].cast<double>();
      if (py::isinstance<py::int_>(t[0])) {
        c.datasets.push_back(cal::Dataset(t[0].cast<int>(), w));
      } else {
        c.datasets.push_back(cal::Dataset(t[0].cast<std::string>(), w));
      }
    }
  }
  c.policy.require_css = d.contains("require_css") && d["require_css"].cast<bool>();
  if (d.contains("cycles") && !d["cycles"].is_none()) c.policy.cycles = d["cycles"].cast<int>();
  if (d.contains("settings")) {
    for (auto item : d["settings"].cast<py::dict>())
      c.policy.settings.emplace_back(item.first.cast<std::string>(), py::str(item.second).cast<std::string>());
  }
  return c;
}

cal::Parameter CalParameterFromDict(const py::dict& d) {
  cal::Parameter p;
  p.path = d["path"].cast<std::string>();
  p.lower = OptionalDouble(d, "lower");
  p.upper = OptionalDouble(d, "upper");
  if (d.contains("scale")) p.scale = ScaleFromName(d["scale"].cast<std::string>());
  if (d.contains("shared")) p.sharing = d["shared"].cast<bool>() ? cal::Sharing::Shared : cal::Sharing::PerCase;
  p.initial = OptionalDouble(d, "initial");
  if (d.contains("case_bounds") && !d["case_bounds"].is_none()) {
    for (auto item : d["case_bounds"].cast<py::list>()) {
      cal::Bounds b;
      if (!item.is_none()) {
        auto t = item.cast<py::tuple>();
        b.lower = t[0].is_none() ? cal::kNaN : t[0].cast<double>();
        b.upper = t[1].is_none() ? cal::kNaN : t[1].cast<double>();
      }
      p.case_bounds.push_back(b);
    }
  }
  if (d.contains("case_initial") && !d["case_initial"].is_none()) {
    for (auto item : d["case_initial"].cast<py::list>())
      p.case_initial.push_back(item.is_none() ? cal::kNaN : item.cast<double>());
  }
  return p;
}

cal::Problem ProblemFromLists(const py::list& cases, const py::list& parameters) {
  std::vector<cal::Case> cs;
  for (auto item : cases) cs.push_back(CaseFromDict(item.cast<py::dict>()));
  std::vector<cal::Parameter> ps;
  for (auto item : parameters) ps.push_back(CalParameterFromDict(item.cast<py::dict>()));
  return cal::Problem(std::move(cs), std::move(ps));
}

cal::Options CalOptionsFromDict(const py::dict& d) {
  cal::Options o;
  if (d.contains("method")) {
    const std::string m = d["method"].cast<std::string>();
    if (m == "least_squares") o.solver = cal::Solver::LeastSquares;
    else if (m == "global") o.solver = cal::Solver::Global;
    else throw std::invalid_argument("method is 'least_squares' or 'global'; got '" + m + "'");
  }
  if (d.contains("refine")) o.refine = d["refine"].cast<bool>();
  if (d.contains("max_simulations")) o.max_simulations = d["max_simulations"].cast<int>();
  if (d.contains("max_iterations")) o.max_iterations = d["max_iterations"].cast<int>();
  if (d.contains("ftol")) o.ftol = d["ftol"].cast<double>();
  if (d.contains("xtol")) o.xtol = d["xtol"].cast<double>();
  if (d.contains("differences")) {
    const std::string s = d["differences"].cast<std::string>();
    if (s == "forward") o.difference = cal::Difference::Forward;
    else if (s == "central") o.difference = cal::Difference::Central;
    else throw std::invalid_argument("differences is 'forward' or 'central'; got '" + s + "'");
  }
  if (d.contains("fd_step")) o.fd_step = d["fd_step"].cast<double>();
  if (d.contains("global_evaluations")) o.global_evaluations = d["global_evaluations"].cast<int>();
  if (d.contains("threads") && !d["threads"].is_none()) o.threads = d["threads"].cast<int>();
  if (d.contains("warm_start")) o.warm_start = d["warm_start"].cast<bool>();
  if (d.contains("cache_capacity")) o.cache_capacity = d["cache_capacity"].cast<size_t>();
  if (d.contains("standard_errors")) o.standard_errors = d["standard_errors"].cast<bool>();
  return o;
}

py::array_t<double> NumpyMatrix(const std::vector<std::vector<double>>& m) {
  const py::ssize_t n = static_cast<py::ssize_t>(m.size());
  const py::ssize_t k = n > 0 ? static_cast<py::ssize_t>(m[0].size()) : 0;
  py::array_t<double> a({n, k});
  auto r = a.mutable_unchecked<2>();
  for (py::ssize_t i = 0; i < n; ++i)
    for (py::ssize_t j = 0; j < k; ++j) r(i, j) = m[i][j];
  return a;
}

py::dict CalResultToDict(const cal::Result& r) {
  py::list params;
  for (const auto& p : r.parameters) {
    py::dict d;
    d["path"] = p.path;
    d["label"] = p.label;
    d["case_index"] = p.case_index;
    d["estimate"] = p.estimate;
    d["standard_error"] = p.standard_error;
    d["initial"] = p.initial;
    d["lower"] = p.lower;
    d["upper"] = p.upper;
    d["at_lower"] = p.at_lower;
    d["at_upper"] = p.at_upper;
    d["scale"] = cal::ScaleName(p.scale);
    params.append(d);
  }
  py::list datasets;
  for (const auto& f : r.datasets) {
    py::dict d;
    d["case_index"] = f.case_index;
    d["case_name"] = f.case_name;
    d["name"] = f.name;
    d["weight"] = f.weight;
    d["points"] = f.points;
    d["skipped"] = f.skipped;
    d["rmse"] = f.rmse;
    d["nrmse"] = f.nrmse;
    d["max_abs_error"] = f.max_abs_error;
    d["objective"] = f.objective;
    d["x"] = ToArray(f.x);
    d["measured"] = ToArray(f.measured);
    d["simulated"] = ToArray(f.simulated);
    d["residuals"] = ToArray(f.residuals);
    d["weighted_residuals"] = ToArray(f.weighted_residuals);
    d["series_x"] = ToArray(f.series.x);
    d["series_y"] = ToArray(f.series.y);
    d["series_si_unit"] = f.series.si_unit;
    datasets.append(d);
  }
  py::list trace;
  for (const auto& t : r.trace) {
    py::dict d;
    d["solver"] = cal::SolverName(t.solver);
    d["iteration"] = t.iteration;
    d["objective"] = t.objective;
    d["values"] = ToArray(t.values);
    d["simulations"] = t.simulations;
    trace.append(d);
  }
  py::dict d;
  d["parameters"] = params;
  d["correlation"] = NumpyMatrix(r.correlation);
  d["datasets"] = datasets;
  d["trace"] = trace;
  d["objective"] = r.objective;
  d["sigma2"] = r.sigma2;
  d["residuals"] = r.residuals;
  d["simulations"] = r.simulations;
  d["failed_simulations"] = r.failed_simulations;
  d["iterations"] = r.iterations;
  d["wall_time_s"] = r.wall_time_s;
  d["solver"] = cal::SolverName(r.solver);
  d["refined"] = r.refined;
  d["stop_reason"] = cal::StopReasonName(r.stop);
  d["message"] = r.message;
  d["notes"] = r.notes;
  d["cancelled"] = r.cancelled;
  d["ok"] = r.Ok();
  d["calibrated_decks"] = r.calibrated_decks;
  return d;
}

py::list UnknownsToList(const cal::Problem& p) {
  py::list out;
  for (const auto& u : p.Unknowns()) {
    py::dict d;
    d["parameter"] = u.parameter;
    d["case_index"] = u.case_index;
    d["path"] = u.path;
    d["label"] = u.label;
    d["lower"] = u.lower;
    d["upper"] = u.upper;
    d["initial"] = u.initial;
    d["scale"] = cal::ScaleName(u.scale);
    out.append(d);
  }
  return out;
}

}  // namespace

PYBIND11_MODULE(_core, m) {
  m.doc() = "Skarstrom PSA/TSA simulation engine (compiled core). Use the 'skarstrom' package.";
  m.attr("__version__") = ADSIM_VERSION_STRING;
  // CSSRecord.distance_ratio above which a CSS stop is reported as far from the CSS.
  m.attr("CSS_DISTANCE_WARNING_RATIO") = adsim::CyclicSteadyState::kDistanceWarningRatio;

  // ---- run status / records ---------------------------------------------
  py::class_<adsim::RunHealthWarning>(m, "RunHealthWarning")
      .def_readonly("code", &adsim::RunHealthWarning::code)
      .def_readonly("value", &adsim::RunHealthWarning::value)
      .def_readonly("threshold", &adsim::RunHealthWarning::threshold)
      .def_readonly("cycle", &adsim::RunHealthWarning::cycle)
      .def_readonly("step", &adsim::RunHealthWarning::step)
      .def_readonly("step_name", &adsim::RunHealthWarning::step_name)
      .def_readonly("component", &adsim::RunHealthWarning::component)
      .def_readonly("suspects", &adsim::RunHealthWarning::suspects)
      .def_readonly("hint", &adsim::RunHealthWarning::hint)
      .def_readonly("message", &adsim::RunHealthWarning::message)
      .def("as_dict", [](const adsim::RunHealthWarning& w) {
        py::dict d;
        d["code"] = w.code;
        d["value"] = w.value;
        d["threshold"] = w.threshold;
        d["cycle"] = w.cycle;
        d["step"] = w.step;
        d["step_name"] = w.step_name;
        d["component"] = w.component;
        d["suspects"] = w.suspects;
        d["hint"] = w.hint;
        d["message"] = w.message;
        return d;
      })
      .def("__repr__", [](const adsim::RunHealthWarning& w) { return w.message; });

  py::class_<adsim::RunHealthCycle>(m, "RunHealthCycle")
      .def_readonly("cycle", &adsim::RunHealthCycle::cycle)
      .def_readonly("mole_fraction_sum_error", &adsim::RunHealthCycle::mole_fraction_sum_error)
      .def_readonly("mole_fraction_sum_step", &adsim::RunHealthCycle::mole_fraction_sum_step)
      .def_readonly("state_change", &adsim::RunHealthCycle::state_change)
      .def_readonly("mass_balance_error", &adsim::RunHealthCycle::mass_balance_error)
      .def_readonly("steps", &adsim::RunHealthCycle::steps)
      .def_readonly("rejected", &adsim::RunHealthCycle::rejected);

  py::class_<adsim::RunHealthReport>(m, "RunHealthReport")
      .def_readonly("components", &adsim::RunHealthReport::components)
      .def_readonly("cycles", &adsim::RunHealthReport::cycles)
      .def_readonly("warnings", &adsim::RunHealthReport::warnings)
      .def_property_readonly("healthy", &adsim::RunHealthReport::Healthy);

  m.attr("RUN_HEALTH_CODES") = [] {
    std::vector<std::string> codes;
    for (const char* c : adsim::kRunHealthCodeList) codes.emplace_back(c);
    return codes;
  }();

  py::class_<adsim::RunResult>(m, "RunResult")
      .def_property_readonly("status", [](const adsim::RunResult& r) { return adsim::RunStatusName(r.status); })
      .def_readonly("message", &adsim::RunResult::message)
      .def_readonly("return_code", &adsim::RunResult::return_code)
      .def_readonly("css_reached", &adsim::RunResult::css_reached)
      .def_readonly("css_cycle", &adsim::RunResult::css_cycle)
      .def_readonly("cycles_run", &adsim::RunResult::cycles_run)
      .def_readonly("wall_time_s", &adsim::RunResult::wall_time_s)
      .def_property_readonly("ok", &adsim::RunResult::Ok)
      .def_readonly("healthy", &adsim::RunResult::healthy,
                    "False when a run-health check fired (the status and results are unchanged)")
      .def_readonly("warnings", &adsim::RunResult::health_warnings,
                    "The run-health warnings (RunHealthWarning: code, value, threshold, cycle, step, suspects, hint)")
      .def_property_readonly("stats", [](const adsim::RunResult& r) {
        py::dict d;
        d["steps"] = r.stats.num_steps;
        d["residual_evaluations"] = r.stats.num_res_evals;
        d["nonlinear_iterations"] = r.stats.num_nonlin_iters;
        d["error_test_failures"] = r.stats.num_err_test_fails;
        d["convergence_failures"] = r.stats.num_conv_fails;
        d["jacobian_residual_evaluations"] = r.stats.num_lin_res_evals;
        d["jacobian_evaluations"] = r.stats.num_jac_evals;
        d["linear_setups"] = r.stats.num_lin_setups;
        d["recovered_steps"] = r.stats.recovered_steps;
        return d;
      })
      .def("__repr__", [](const adsim::RunResult& r) {
        return "<RunResult " + std::string(adsim::RunStatusName(r.status)) + ": " + r.message + ">";
      });

  py::class_<adsim::StepEndPressure>(m, "StepEndPressure")
      .def_readonly("cycle", &adsim::StepEndPressure::cycle)
      .def_readonly("step", &adsim::StepEndPressure::step)
      .def_readonly("step_id", &adsim::StepEndPressure::step_id)
      .def_readonly("step_name", &adsim::StepEndPressure::step_name)
      .def_readonly("time", &adsim::StepEndPressure::time)
      .def_readonly("bottom", &adsim::StepEndPressure::bottom)
      .def_readonly("top", &adsim::StepEndPressure::top)
      .def_readonly("bed_min", &adsim::StepEndPressure::bed_min)
      .def_readonly("bed_mean", &adsim::StepEndPressure::bed_mean)
      .def_readonly("bed_max", &adsim::StepEndPressure::bed_max)
      .def("as_dict", [](const adsim::StepEndPressure& r) {
        py::dict d;
        d["cycle"] = r.cycle;
        d["step"] = r.step;
        d["step_id"] = r.step_id;
        d["step_name"] = r.step_name;
        d["time"] = r.time;
        d["bottom"] = r.bottom;
        d["top"] = r.top;
        d["bed_min"] = r.bed_min;
        d["bed_mean"] = r.bed_mean;
        d["bed_max"] = r.bed_max;
        return d;
      })
      .def("__repr__", [](const adsim::StepEndPressure& r) {
        return "StepEndPressure(cycle=" + std::to_string(r.cycle) + ", step='" + r.step_name +
               "', bottom=" + std::to_string(r.bottom) + ", top=" + std::to_string(r.top) + ")";
      });

  py::class_<adsim::KPIRecord>(m, "KPIRecord")
      .def_readonly("cycle", &adsim::KPIRecord::cycle)
      .def_readonly("component", &adsim::KPIRecord::component)
      .def_readonly("valid", &adsim::KPIRecord::valid)
      .def_readonly("purity", &adsim::KPIRecord::purity)
      .def_readonly("recovery", &adsim::KPIRecord::recovery)
      .def_readonly("productivity_mass", &adsim::KPIRecord::productivity_mass)
      .def_readonly("productivity_volume", &adsim::KPIRecord::productivity_volume)
      .def_readonly("energy_consumption", &adsim::KPIRecord::energy_consumption)
      .def_readonly("production_rate", &adsim::KPIRecord::production_rate)
      .def_readonly("electrical_energy", &adsim::KPIRecord::electrical_energy)
      .def_readonly("thermal_energy", &adsim::KPIRecord::thermal_energy)
      .def_readonly("cooling_duty", &adsim::KPIRecord::cooling_duty)
      .def_readonly("step_thermal_energy", &adsim::KPIRecord::step_thermal_energy)
      .def_readonly("step_cooling_duty", &adsim::KPIRecord::step_cooling_duty)
      .def_readonly("trapezoid_purity", &adsim::KPIRecord::trapezoid_purity)
      .def_readonly("trapezoid_recovery", &adsim::KPIRecord::trapezoid_recovery)
      .def_readonly("trapezoid_productivity_mass", &adsim::KPIRecord::trapezoid_productivity_mass)
      .def_readonly("trapezoid_productivity_volume", &adsim::KPIRecord::trapezoid_productivity_volume)
      .def_readonly("trapezoid_energy_consumption", &adsim::KPIRecord::trapezoid_energy_consumption)
      .def_readonly("trapezoid_production_rate", &adsim::KPIRecord::trapezoid_production_rate)
      .def_readonly("mass_balance_error", &adsim::KPIRecord::mass_balance_error)
      .def_readonly("trapezoid_mass_balance_error", &adsim::KPIRecord::trapezoid_mass_balance_error)
      .def_readonly("uncounted_inflow", &adsim::KPIRecord::uncounted_inflow)
      .def_readonly("uncounted_inflow_fraction", &adsim::KPIRecord::uncounted_inflow_fraction)
      .def_readonly("reaction_kpis", &adsim::KPIRecord::reaction_kpis)
      .def_readonly("reaction_generation", &adsim::KPIRecord::reaction_generation)
      .def_readonly("conversion", &adsim::KPIRecord::conversion)
      .def_readonly("selectivity", &adsim::KPIRecord::selectivity)
      .def_readonly("reaction_yield", &adsim::KPIRecord::reaction_yield)
      .def_readonly("space_time_yield", &adsim::KPIRecord::space_time_yield)
      .def_readonly("space_time_yield_unit", &adsim::KPIRecord::space_time_yield_unit)
      .def_readonly("peak_temperature", &adsim::KPIRecord::peak_temperature)
      .def_readonly("peak_temperature_position", &adsim::KPIRecord::peak_temperature_position)
      .def_readonly("peak_temperature_time", &adsim::KPIRecord::peak_temperature_time)
      .def("as_dict", [](const adsim::KPIRecord& k) {
        py::dict d;
        d["cycle"] = k.cycle;
        d["component"] = k.component;
        d["valid"] = k.valid;
        d["purity"] = k.purity;
        d["recovery"] = k.recovery;
        d["productivity_mass"] = k.productivity_mass;
        d["productivity_volume"] = k.productivity_volume;
        d["energy_consumption"] = k.energy_consumption;
        d["production_rate"] = k.production_rate;
        d["electrical_energy"] = k.electrical_energy;
        d["thermal_energy"] = k.thermal_energy;
        d["cooling_duty"] = k.cooling_duty;
        d["step_thermal_energy"] = k.step_thermal_energy;
        d["step_cooling_duty"] = k.step_cooling_duty;
        d["trapezoid_purity"] = k.trapezoid_purity;
        d["trapezoid_recovery"] = k.trapezoid_recovery;
        d["trapezoid_productivity_mass"] = k.trapezoid_productivity_mass;
        d["trapezoid_productivity_volume"] = k.trapezoid_productivity_volume;
        d["trapezoid_energy_consumption"] = k.trapezoid_energy_consumption;
        d["trapezoid_production_rate"] = k.trapezoid_production_rate;
        d["mass_balance_error"] = k.mass_balance_error;
        d["trapezoid_mass_balance_error"] = k.trapezoid_mass_balance_error;
        d["uncounted_inflow"] = k.uncounted_inflow;
        d["uncounted_inflow_fraction"] = k.uncounted_inflow_fraction;
        // The reaction KPIs only on the records that carry them, so a deck
        // without reactions keeps its keys.
        if (k.reaction_kpis) {
          d["reaction_kpis"] = k.reaction_kpis;
          d["reaction_generation"] = k.reaction_generation;
          d["conversion"] = k.conversion;
          d["selectivity"] = k.selectivity;
          d["reaction_yield"] = k.reaction_yield;
          d["space_time_yield"] = k.space_time_yield;
          d["space_time_yield_unit"] = k.space_time_yield_unit;
          d["peak_temperature"] = k.peak_temperature;
          d["peak_temperature_position"] = k.peak_temperature_position;
          d["peak_temperature_time"] = k.peak_temperature_time;
        }
        return d;
      })
      .def("__repr__", [](const adsim::KPIRecord& k) {
        return "<KPIRecord cycle=" + std::to_string(k.cycle) + " " + k.component +
               " purity=" + std::to_string(k.purity) + " recovery=" + std::to_string(k.recovery) + ">";
      });

  py::class_<adsim::CSSRecord>(m, "CSSRecord")
      .def_readonly("cycle", &adsim::CSSRecord::cycle)
      .def_readonly("converged", &adsim::CSSRecord::converged)
      .def_readonly("error", &adsim::CSSRecord::error)
      .def_readonly("tolerance", &adsim::CSSRecord::tolerance)
      .def_readonly("method", &adsim::CSSRecord::method)
      .def_readonly("state_change", &adsim::CSSRecord::state_change)
      .def_readonly("state_tolerance", &adsim::CSSRecord::state_tolerance)
      .def_readonly("passing_cycles", &adsim::CSSRecord::passing_cycles)
      .def_readonly("required_cycles", &adsim::CSSRecord::required_cycles)
      .def_readonly("distance_estimated", &adsim::CSSRecord::distance_estimated,
                    "True when the remaining distance to the cyclic steady state was estimated")
      .def_readonly("convergence_rate", &adsim::CSSRecord::convergence_rate,
                    "Estimated cycle-to-cycle convergence rate rho in [0, 0.999] (nan if not estimated)")
      .def_readonly("distance", &adsim::CSSRecord::distance,
                    "Estimated remaining distance to the CSS, last state change * rho / (1 - rho), "
                    "in the norm of state_change (nan if not estimated)")
      .def_readonly("distance_ratio", &adsim::CSSRecord::distance_ratio,
                    "distance / state_tolerance (nan if not estimated or the state test is off)")
      .def_readonly("distance_note", &adsim::CSSRecord::distance_note,
                    "Why the distance was not estimated (empty when it was)")
      .def_readonly("newton_iterations", &adsim::CSSRecord::newton_iterations,
                    "Newton-Krylov acceleration: Newton steps of the solve made after this cycle "
                    "(-1 when there was none)")
      .def_readonly("newton_cycles", &adsim::CSSRecord::newton_cycles,
                    "Newton-Krylov acceleration: extra cycles (cycle equivalents) that solve ran")
      .def_readonly("newton_residuals", &adsim::CSSRecord::newton_residuals,
                    "Newton-Krylov acceleration: the residual |G(x) - x| at each Newton iterate, the start "
                    "first, in the norm of state_change")
      .def_readonly("newton_krylov_vectors", &adsim::CSSRecord::newton_krylov_vectors,
                    "Newton-Krylov acceleration: Krylov vectors (one cycle each) per Newton step")
      .def_readonly("floquet_multiplier", &adsim::CSSRecord::floquet_multiplier,
                    "Newton-Krylov acceleration: the dominant Floquet multiplier of the cycle map, estimated "
                    "from the last Krylov basis (nan if not estimated)")
      .def_readonly("newton_outcome", &adsim::CSSRecord::newton_outcome,
                    "Newton-Krylov acceleration: why the solve ended")
      .def("__repr__", [](const adsim::CSSRecord& c) {
        return "<CSSRecord cycle=" + std::to_string(c.cycle) + (c.converged ? " converged" : "") +
               " error=" + std::to_string(c.error) + ">";
      });

  py::class_<adsim::StepOperatingSummary>(m, "StepOperatingSummary")
      .def_readonly("step", &adsim::StepOperatingSummary::step)
      .def_readonly("step_name", &adsim::StepOperatingSummary::step_name)
      .def_readonly("max_pressure", &adsim::StepOperatingSummary::max_pressure)
      .def_readonly("min_pressure", &adsim::StepOperatingSummary::min_pressure)
      .def_readonly("max_wall_temperature", &adsim::StepOperatingSummary::max_wall_temperature)
      .def_readonly("max_approach_to_fluidisation", &adsim::StepOperatingSummary::max_approach_to_fluidisation)
      .def_readonly("fluidisation_position", &adsim::StepOperatingSummary::fluidisation_position)
      .def_readonly("fluidisation_time", &adsim::StepOperatingSummary::fluidisation_time)
      .def("__repr__", [](const adsim::StepOperatingSummary& s) {
        return "<StepOperatingSummary step=" + std::to_string(s.step) + " '" + s.step_name + "'>";
      });

  py::class_<adsim::OperatingWarning>(m, "OperatingWarning")
      .def_property_readonly("limit", [](const adsim::OperatingWarning& w) {
        return std::string(adsim::kOperatingLimitList[static_cast<int>(w.limit)]);
      })
      .def_readonly("step", &adsim::OperatingWarning::step)
      .def_readonly("value", &adsim::OperatingWarning::value)
      .def_readonly("limit_value", &adsim::OperatingWarning::limit_value)
      .def_readonly("message", &adsim::OperatingWarning::message)
      .def("__repr__", [](const adsim::OperatingWarning& w) { return w.message; });

  py::class_<adsim::OperatingLimitsReport>(m, "OperatingLimitsReport")
      .def_readonly("cycle", &adsim::OperatingLimitsReport::cycle)
      .def_readonly("design_checked", &adsim::OperatingLimitsReport::design_checked)
      .def_readonly("steps", &adsim::OperatingLimitsReport::steps)
      .def_readonly("warnings", &adsim::OperatingLimitsReport::warnings);

  py::class_<adsim::WallDesignResult>(m, "WallDesignResult")
      .def_readonly("valid", &adsim::WallDesignResult::valid)
      .def_readonly("message", &adsim::WallDesignResult::message)
      .def_readonly("notes", &adsim::WallDesignResult::notes)
      .def_readonly("allowable_stress", &adsim::WallDesignResult::allowable_stress)
      .def_readonly("internal_thickness", &adsim::WallDesignResult::internal_thickness)
      .def_readonly("external_thickness", &adsim::WallDesignResult::external_thickness)
      .def_readonly("required_thickness", &adsim::WallDesignResult::required_thickness)
      .def_readonly("external_governs", &adsim::WallDesignResult::external_governs)
      .def_readonly("thick_shell", &adsim::WallDesignResult::thick_shell);

  py::class_<adsim::TimeSeries>(m, "TimeSeries")
      .def_readonly("variable", &adsim::TimeSeries::variable)
      .def_readonly("component", &adsim::TimeSeries::component)
      .def_readonly("adsorbent", &adsim::TimeSeries::adsorbent)
      .def_readonly("reaction", &adsim::TimeSeries::reaction)
      .def_readonly("unit", &adsim::TimeSeries::unit)
      .def_property_readonly("time", [](const adsim::TimeSeries& t) { return ToArray(t.time); })
      .def_property_readonly("values", [](const adsim::TimeSeries& t) { return ToArray(t.values); })
      .def("__len__", [](const adsim::TimeSeries& t) { return t.time.size(); })
      .def("__repr__", [](const adsim::TimeSeries& t) {
        return "<TimeSeries " + t.variable + (t.component.empty() ? "" : "[" + t.component + "]") +
               (t.reaction.empty() ? "" : " " + t.reaction) +
               (t.adsorbent.empty() ? "" : " (" + t.adsorbent + ")") + " [" +
               t.unit + "], " + std::to_string(t.time.size()) + " points>";
      });

  py::class_<adsim::AxialProfile>(m, "AxialProfile")
      .def_readonly("variable", &adsim::AxialProfile::variable)
      .def_readonly("component", &adsim::AxialProfile::component)
      .def_readonly("adsorbent", &adsim::AxialProfile::adsorbent)
      .def_readonly("reaction", &adsim::AxialProfile::reaction)
      .def_readonly("unit", &adsim::AxialProfile::unit)
      .def_readonly("cycle", &adsim::AxialProfile::cycle)
      .def_readonly("step", &adsim::AxialProfile::step)
      .def_readonly("faces", &adsim::AxialProfile::faces)
      .def_property_readonly("z", [](const adsim::AxialProfile& p) { return ToArray(p.z); })
      .def_property_readonly("time", [](const adsim::AxialProfile& p) { return ToArray(p.time); })
      .def_property_readonly("values",
                             [](const adsim::AxialProfile& p) { return ToArray2D(p.values, p.rows, p.cols); })
      .def("__repr__", [](const adsim::AxialProfile& p) {
        return "<AxialProfile " + p.variable + (p.component.empty() ? "" : "[" + p.component + "]") +
               (p.reaction.empty() ? "" : " " + p.reaction) +
               (p.adsorbent.empty() ? "" : " (" + p.adsorbent + ")") + " [" +
               p.unit + "], cycle " + std::to_string(p.cycle) + " step " + std::to_string(p.step) + ", " +
               std::to_string(p.rows) + " x " + std::to_string(p.cols) + ">";
      });

  // ---- CycleState (adsim/solver/cycle_state.h) -----------------------------
  py::class_<adsim::CycleState>(m, "CycleState")
      .def_readonly("components", &adsim::CycleState::components)
      .def_readonly("cycles", &adsim::CycleState::cycles)
      .def_property_readonly("layers", [](const adsim::CycleState& c) {
        py::list out;
        for (const auto& l : c.layers) {
          py::dict d;
          d["id"] = l.id;
          d["name"] = l.name;
          d["cells"] = l.cells;
          out.append(d);
        }
        return out;
      })
      .def_property_readonly("steps", [](const adsim::CycleState& c) {
        py::list out;
        for (const auto& st : c.steps) {
          py::dict d;
          d["id"] = st.id;
          d["name"] = st.name;
          d["time"] = st.time;
          d["endpoint_ids"] = st.endpoint_ids;
          out.append(d);
        }
        return out;
      })
      .def_property_readonly("fields", [](const adsim::CycleState& c) {
        py::list out;
        for (const auto& f : c.fields) {
          py::dict d;
          d["name"] = f.name;
          d["component"] = f.component;
          d["adsorbent"] = f.adsorbent;
          d["differential"] = f.differential;
          d["index"] = f.index;
          d["cell"] = f.cell;
          out.append(d);
        }
        return out;
      })
      .def_property_readonly("endpoints", [](const adsim::CycleState& c) {
        py::list out;
        for (const auto& e : c.endpoints) {
          py::dict d;
          d["endpoint_id"] = e.endpoint_id;
          d["step_id"] = e.step_id;
          d["name"] = e.name;
          d["history_points"] = e.history.time_series_.size();
          d["summary_valid"] = e.summary.is_valid;
          d["summary_total_moles"] = e.summary.total_moles;
          out.append(d);
        }
        return out;
      })
      .def_property("y", [](const adsim::CycleState& c) { return ToArray(c.y); },
                    [](adsim::CycleState& c, const std::vector<double>& v) {
                      if (v.size() != c.y.size())
                        throw std::invalid_argument("y: expected " + std::to_string(c.y.size()) + " values, got "
                                                    + std::to_string(v.size()));
                      c.y = v;
                    })
      .def_property_readonly("yp", [](const adsim::CycleState& c) { return ToArray(c.yp); })
      .def_property_readonly("differential_index", [](const adsim::CycleState& c) { return c.differential_index; })
      .def_property("differential", [](const adsim::CycleState& c) { return ToArray(c.Differential()); },
                    &adsim::CycleState::SetDifferential)
      .def_property_readonly("has_history", &adsim::CycleState::HasHistory)
      .def_property_readonly("size", &adsim::CycleState::Size)
      .def("field", [](const adsim::CycleState& c, const std::string& name, const std::string& component,
                       size_t adsorbent) { return ToArray(c.FieldValues(name, component, adsorbent)); },
           py::arg("name"), py::arg("component") = "", py::arg("adsorbent") = 0)
      .def("copy", [](const adsim::CycleState& c) { return c; })
      .def("__copy__", [](const adsim::CycleState& c) { return c; })
      .def("__deepcopy__", [](const adsim::CycleState& c, py::dict) { return c; }, py::arg("memo"))
      .def("__len__", &adsim::CycleState::Size)
      .def("__repr__", [](const adsim::CycleState& c) {
        return "<CycleState " + std::to_string(c.Size()) + " unknowns, " + std::to_string(c.endpoints.size())
               + " endpoints, after " + std::to_string(c.cycles) + " cycles>";
      });

  // ---- Simulation ----------------------------------------------------------
  py::class_<adsim::Simulation>(m, "Simulation")
      .def(py::init<>())
      .def_static("from_file", &adsim::Simulation::FromFile, py::arg("path"))
      .def_static("from_json", &adsim::Simulation::FromJson, py::arg("text"))
      .def("load", &adsim::Simulation::Load, py::arg("path"))
      .def("load_json", &adsim::Simulation::LoadJson, py::arg("text"))
      .def("save", &adsim::Simulation::Save, py::arg("path"))
      .def("to_json", &adsim::Simulation::ToJson)
      .def("sensitivity_study", [](const adsim::Simulation& s) { return StudyToDict(s.Data().sensitivity_study); })
      .def("set_sensitivity_study", [](adsim::Simulation& s, const py::dict& study) {
        s.Data().sensitivity_study = StudyFromDict(study);
      }, py::arg("study"))
      .def("check", [](adsim::Simulation& s) {
        std::string err;
        const bool ok = s.Check(err);
        return py::make_tuple(ok, err);
      })

      // Parameters: a list of dicts (the registry's ParameterInfo + value).
      .def("parameters", [](adsim::Simulation& s, bool tunable_only) {
        py::list out;
        for (auto& p : s.Parameters().All()) {
          const auto& i = p.Info();
          if (tunable_only && !i.tunable) continue;
          py::dict d;
          d["path"] = i.path;
          d["alias"] = i.alias;
          d["category"] = i.category;
          d["subcategory"] = i.subcategory;
          d["name"] = i.name;
          d["unit"] = i.unit;
          d["unit_checked"] = i.unit_checked;
          d["choices"] = i.choices;
          d["kind"] = adsim::ParameterKindName(i.kind);
          d["tunable"] = i.tunable;
          d["lower"] = i.lower;
          d["upper"] = i.upper;
          if (i.kind == adsim::ParameterKind::String) d["value"] = p.GetText();
          else if (i.kind == adsim::ParameterKind::Boolean) d["value"] = p.GetValue() != 0.0;
          else d["value"] = p.GetValue();
          d["text"] = p.GetText();
          out.append(d);
        }
        return out;
      }, py::arg("tunable_only") = false)
      .def("get", [](adsim::Simulation& s, const std::string& path) -> py::object {
        auto& p = s.Parameters().Get(path);
        if (p.Info().kind == adsim::ParameterKind::String) return py::str(p.GetText());
        if (p.Info().kind == adsim::ParameterKind::Integer) return py::int_(static_cast<long long>(p.GetValue()));
        if (p.Info().kind == adsim::ParameterKind::Boolean) return py::bool_(p.GetValue() != 0.0);
        return py::float_(p.GetValue());
      }, py::arg("path"))
      .def("get_text", &adsim::Simulation::GetParameterText, py::arg("path"))
      .def("set_value", &adsim::Simulation::SetParameter, py::arg("path"), py::arg("value"))
      .def("set_text", &adsim::Simulation::SetParameterText, py::arg("path"), py::arg("text"))

      // Blocking run. The GIL is released while the solver works; callbacks
      // re-acquire it. Ctrl-C (or any pending Python signal) cancels the run
      // and is re-raised afterwards.
      .def("run", [](adsim::Simulation& s, std::optional<int> cycles, py::object progress,
                     const std::string& log_level, py::object on_message, const adsim::CycleState* initial_state) {
        // Copied: the run replaces the final state, which `initial_state` may be.
        std::optional<adsim::CycleState> start;
        if (initial_state) start = *initial_state;
        return RunFromPython(cycles, progress, log_level, on_message, [&s, &start](const adsim::RunOptions& o) {
          adsim::RunOptions opt = o;
          if (start) opt.initial_state = &*start;
          return s.Run(opt);
        });
      }, py::arg("cycles") = py::none(), py::arg("progress") = py::none(), py::arg("log_level") = "low",
         py::arg("on_message") = py::none(), py::arg("initial_state") = py::none())
      .def("run_cycle", [](adsim::Simulation& s, const adsim::CycleState* state, py::object progress,
                           const std::string& log_level, py::object on_message) {
        return RunFromPython(std::nullopt, progress, log_level, on_message, [&s, state](const adsim::RunOptions& o) {
          return state ? s.RunCycle(*state, o) : s.RunCycle(o);
        });
      }, py::arg("state") = py::none(), py::arg("progress") = py::none(), py::arg("log_level") = "low",
         py::arg("on_message") = py::none())
      .def_property_readonly("has_final_state", &adsim::Simulation::HasFinalState)
      .def("final_state", [](const adsim::Simulation& s) { return s.FinalState(); })
      .def("set_initial_state", [](adsim::Simulation& s, const adsim::CycleState& state) { s.SetInitialState(state); },
           py::arg("state"))
      .def("clear_initial_state", &adsim::Simulation::ClearInitialState)
      .def_property_readonly("initial_state", [](const adsim::Simulation& s) -> py::object {
        if (!s.HasInitialState()) return py::none();
        return py::cast(s.InitialState());
      })
      .def("check_state", &adsim::Simulation::CheckState, py::arg("state"))
      .def("prepare_state", [](adsim::Simulation& s, const adsim::CycleState& state) {
        std::string note;
        adsim::CycleState out = s.PrepareState(state, &note);
        return py::make_tuple(std::move(out), note);
      }, py::arg("state"))
      .def("request_stop", &adsim::Simulation::RequestStop)
      .def_property_readonly("last_result", &adsim::Simulation::LastResult)
      .def_property_readonly("messages", [](const adsim::Simulation& s) {
        py::list out;
        for (const auto& m : s.Messages()) out.append(py::make_tuple(LogLevelName(m.level), m.text));
        return out;
      })

      // Results.
      .def_property_readonly("has_results", &adsim::Simulation::HasResults)
      .def_property_readonly("cycles", &adsim::Simulation::Cycles)
      .def_property_readonly("last_cycle", &adsim::Simulation::LastCycle)
      .def_property_readonly("components", &adsim::Simulation::Components)
      .def_property_readonly("steps", &adsim::Simulation::StepNames)
      .def_property_readonly("positions", [](const adsim::Simulation& s) { return ToArray(s.Positions()); })
      .def_property_readonly("face_positions", [](const adsim::Simulation& s) { return ToArray(s.FacePositions()); })
      .def("streams", [](const adsim::Simulation& s) {
        py::list out;
        for (const auto& st : s.Streams()) {
          py::dict d;
          d["step"] = st.step;
          d["step_name"] = st.step_name;
          d["name"] = st.name;
          d["endpoint_id"] = st.endpoint_id;
          d["location"] = st.location == adsim::BoundaryLocation::TOP ? "top" : "bottom";
          out.append(d);
        }
        return out;
      })
      .def("kpis", [](const adsim::Simulation& s, std::optional<int> cycle) {
        return cycle ? s.KPIs(*cycle) : s.KPIs();
      }, py::arg("cycle") = py::none())
      .def("step_end_pressures", &adsim::Simulation::StepEndPressures, py::arg("cycle") = -1)
      .def("css_history", &adsim::Simulation::CSSHistory)
      .def("operating_limits", &adsim::Simulation::OperatingLimits, py::arg("cycle") = -1)
      .def("run_health", &adsim::Simulation::RunHealth, py::return_value_policy::copy)
      .def("wall_design", &adsim::Simulation::WallDesign)
      .def("series", [](const adsim::Simulation& s, const std::string& variable, const std::string& component,
                        const std::string& location, double z, int cycle, const std::string& adsorbent,
                        const std::string& reaction) {
        return s.Series(variable, component, ParseLocation(location), z, cycle, adsorbent, reaction);
      }, py::arg("variable"), py::arg("component") = "", py::arg("location") = "top", py::arg("z") = 0.0,
         py::arg("cycle") = -1, py::arg("adsorbent") = "", py::arg("reaction") = "")
      .def("profile", &adsim::Simulation::Profile, py::arg("variable"), py::arg("component") = "",
           py::arg("cycle") = -1, py::arg("step") = 0, py::arg("adsorbent") = "", py::arg("reaction") = "")
      .def("series_reactions", [](const adsim::Simulation& s) {
        py::list out;
        for (const auto& r : s.SeriesReactions()) {
          py::dict d;
          d["id"] = r.id;
          d["name"] = r.name;
          d["adsorbent_id"] = r.adsorbent_id;
          d["adsorbent"] = r.adsorbent;
          d["variable_name"] = r.variable_name;
          out.append(d);
        }
        return out;
      })
      .def("series_adsorbents", [](const adsim::Simulation& s) {
        py::list out;
        for (const auto& a : s.SeriesAdsorbents()) {
          py::dict d;
          d["id"] = a.id;
          d["label"] = a.label;
          d["logged"] = a.logged;
          out.append(d);
        }
        return out;
      })
      .def("stream_series", &adsim::Simulation::StreamSeries, py::arg("step"), py::arg("stream"),
           py::arg("variable"), py::arg("component") = "", py::arg("cycle") = -1)
      .def("endpoint_series", &adsim::Simulation::EndpointSeries, py::arg("endpoint_id"), py::arg("variable"),
           py::arg("component") = "", py::arg("cycle") = -1)
      .def("save_logs", &adsim::Simulation::SaveLogs, py::arg("directory"))
      // ---- experimental data stored in the deck (adsim/solver/experimental_data.h) ----
      .def("experimental_data", [](adsim::Simulation& s) {
        py::list out;
        for (const auto& xy : s.Data().exp_data_.GetAllData()) out.append(DatasetToDict(xy));
        return out;
      })
      .def("add_experimental_data", [](adsim::Simulation& s, const std::string& name, std::vector<double> x,
                                       std::vector<double> y, const std::string& x_label, const std::string& y_label,
                                       const py::dict& observable) {
        adsim::XYData xy;
        xy.name = name;
        xy.x = std::move(x);
        xy.y = std::move(y);
        xy.x_label = x_label;
        xy.y_label = y_label;
        xy.observable = BindingFromDict(observable);
        size_t index = 0;
        std::string err;
        if (!adsim::AddExperimentalData(s.Data(), std::move(xy), index, err)) throw py::value_error(err);
        return index;
      }, py::arg("name"), py::arg("x"), py::arg("y"), py::arg("x_label") = "", py::arg("y_label") = "",
         py::arg("observable") = py::dict())
      .def("set_experimental_data_binding", [](adsim::Simulation& s, size_t index, const py::dict& observable) {
        std::string err;
        if (!adsim::SetExperimentalDataBinding(s.Data(), index, BindingFromDict(observable), err))
          throw py::value_error(err);
      }, py::arg("index"), py::arg("observable"))
      .def("remove_experimental_data", [](adsim::Simulation& s, size_t index) {
        std::string err;
        if (!adsim::RemoveExperimentalData(s.Data(), index, err)) throw py::index_error(err);
      }, py::arg("index"))
      // The simulated counterpart of a dataset (stored: by index; or given:
      // name, x, y, observable) and its comparison (parameterisation.h).
      .def("simulated_experimental_data", [](const adsim::Simulation& s, size_t index) {
        const auto& all = s.Data().exp_data_.GetAllData();
        if (index >= all.size()) throw py::index_error("no experimental dataset " + std::to_string(index));
        const auto sim = Simulated(s, all[index]);
        return py::make_tuple(ToArray(sim.x), ToArray(sim.y));
      }, py::arg("index"))
      .def("simulated_dataset", [](const adsim::Simulation& s, const std::string& name, std::vector<double> x,
                                   std::vector<double> y, const py::dict& observable) {
        const auto xy = DatasetFromArgs(s, name, std::move(x), std::move(y), observable);
        const auto sim = Simulated(s, xy);
        return py::make_tuple(ToArray(sim.x), ToArray(sim.y));
      }, py::arg("name"), py::arg("x"), py::arg("y"), py::arg("observable"))
      .def("compare_dataset", [](const adsim::Simulation& s, const std::string& name, std::vector<double> x,
                                 std::vector<double> y, const py::dict& observable) {
        const auto xy = DatasetFromArgs(s, name, std::move(x), std::move(y), observable);
        const auto c = adsim::CompareWithSimulation(xy, Simulated(s, xy));
        py::dict d;
        d["rmse"] = c.rmse;
        d["nrmse"] = c.nrmse;
        d["max_abs_error"] = c.max_abs_error;
        d["points"] = c.points;
        d["weighted_sse"] = c.weighted_sse;
        return d;
      }, py::arg("name"), py::arg("x"), py::arg("y"), py::arg("observable"))
      // ---- structural deck edits (adsim/objects/deck_edit.h) ----------------
      .def("add_component", [](adsim::Simulation& s, const std::string& name) {
        std::string err;
        if (!adsim::AddComponentToDeck(s.Data(), name, err)) throw py::value_error(err);
        s.RefreshParameters();
      }, py::arg("name"))
      .def("remove_component", [](adsim::Simulation& s, const std::string& name) {
        adsim::RemoveComponentFromDeck(s.Data(), name);
        s.RefreshParameters();
      }, py::arg("name"))
      .def("add_library_adsorbent", [](adsim::Simulation& s, const std::string& library_dir, const std::string& key) {
        adsim::AdsorbentID id{};
        std::string err;
        if (!adsim::AddLibraryAdsorbent(s.Data(), library_dir, key, id, err)) throw py::value_error(err);
        s.RefreshParameters();
        return static_cast<long long>(id);
      }, py::arg("library_dir"), py::arg("key"))
      // ---- building a deck from scratch (adsim/objects/deck_builder.h) ----------
      // Each call changes the structure, so the registry is refreshed after it.
      .def("build", [](adsim::Simulation& s, const std::string& op, const py::list& args) -> py::object {
        namespace b = adsim::build;
        auto& d = s.Data();
        auto str = [&](size_t i) { return args[i].cast<std::string>(); };
        auto id = [&](size_t i) { return args[i].cast<size_t>(); };
        py::object out = py::none();
        if (op == "add_component") b::AddComponent(d, str(0));
        else if (op == "add_adsorbent") out = py::int_(b::AddAdsorbent(d, str(0)));
        else if (op == "set_isotherm") b::SetIsotherm(d, id(0), str(1), args[2].cast<std::vector<std::string>>());
        else if (op == "set_mixing_rule") b::SetMixingRule(d, id(0), str(1));
        else if (op == "set_rate_model") b::SetRateModel(d, id(0), str(1), str(2), str(3));
        else if (op == "set_bed_shape") b::SetBedShape(d, str(0));
        else if (op == "set_wall_material") b::SetWallMaterial(d, str(0));
        else if (op == "set_wall_boundary") b::SetWallBoundary(d, str(0));
        else if (op == "set_thermal_equilibrium") b::SetThermalEquilibrium(d, str(0));
        else if (op == "set_particle_heat_transfer") b::SetParticleHeatTransfer(d, str(0));
        else if (op == "add_layer") out = py::int_(b::AddLayer(d, str(0)));
        else if (op == "set_layer_composition")
          b::SetLayerComposition(d, id(0), args[1].cast<std::vector<std::pair<adsim::AdsorbentID, double>>>());
        else if (op == "set_momentum") b::SetMomentum(d, id(0), str(1));
        else if (op == "set_film_correlation") b::SetFilmCorrelation(d, id(0), str(1));
        else if (op == "set_dispersion") b::SetDispersion(d, id(0), str(1));
        else if (op == "set_heat_dispersion") b::SetHeatDispersion(d, id(0), str(1));
        else if (op == "set_internal_heat_transfer") b::SetInternalHeatTransfer(d, id(0), str(1));
        else if (op == "set_external_heat_transfer") b::SetExternalHeatTransfer(d, id(0), str(1));
        else if (op == "add_step") out = py::int_(b::AddStep(d, str(0)));
        else if (op == "set_closed") b::SetClosed(d, id(0), str(1));
        else if (op == "set_pressure_boundary") out = py::int_(b::SetPressureBoundary(d, id(0), str(1), str(2), str(3)));
        else if (op == "set_flow_boundary") b::SetFlowBoundary(d, id(0), str(1), str(2));
        else if (op == "add_stream") out = py::int_(b::AddStream(d, id(0), str(1), str(2), str(3), str(4)));
        else if (op == "add_valve_stream") out = py::int_(b::AddValveStream(d, id(0), str(1), str(2), str(3)));
        else if (op == "add_pump_stream")
          out = py::int_(b::AddPumpStream(d, id(0), str(1), str(2), str(3), str(4), str(5)));
        else if (op == "set_buffer_tank") b::SetBufferTank(d, id(0), id(1));
        else if (op == "set_temperature_profile") b::SetTemperatureProfile(d, id(0), str(1));
        else if (op == "set_pump") b::SetPump(d, id(0), str(1), str(2), str(3));
        else if (op == "connect") b::Connect(d, id(0), id(1));
        else if (op == "set_energy_source") b::SetEnergySource(d, id(0), id(1), args[2].cast<bool>());
        else if (op == "set_kpi")
          b::SetKPI(d, str(0), args[1].cast<std::vector<size_t>>(), args[2].cast<std::vector<size_t>>(),
                    args.size() > 3 ? args[3].cast<std::vector<size_t>>() : std::vector<size_t>{});
        else if (op == "add_reaction")  // structure only; values through the registry
          out = py::int_(b::AddReaction(d, id(0), str(1), args[2].cast<std::map<std::string, double>>(),
                                        args.size() > 3 ? str(3) : "power_law"));
        else if (op == "set_reaction_kind") b::SetReactionKind(d, id(0), id(1), str(2));
        else if (op == "set_reaction_stoichiometry")
          b::SetReactionStoichiometry(d, id(0), id(1), args[2].cast<std::map<std::string, double>>());
        else if (op == "set_reaction_basis") b::SetReactionBasis(d, id(0), id(1), str(2));
        else if (op == "set_effectiveness")
          b::SetEffectiveness(d, id(0), id(1), str(2), args.size() > 3 ? str(3) : std::string{});
        else if (op == "add_adsorption_term") out = py::int_(b::AddAdsorptionTerm(d, id(0), id(1)));
        else if (op == "remove_adsorption_term") b::RemoveAdsorptionTerm(d, id(0), id(1), id(2));
        else if (op == "remove_reaction") b::RemoveReaction(d, id(0), id(1));
        else if (op == "set_reaction_kpi")
          b::SetReactionKPI(d, str(0), args[1].cast<std::vector<std::string>>(), args.size() > 2 ? str(2)
                                                                                                 : "catalyst_mass");
        else if (op == "add_variable") b::AddVariable(d, str(0), str(1), str(2), str(3));
        else throw py::value_error("unknown build operation '" + op + "'");
        s.RefreshParameters();
        return out;
      }, py::arg("op"), py::arg("args"))
      // ---- structural cascades (adsim/objects/deck_edit.h) -------------------
      // Kinds are named by their tags (adsim/objects/model_kinds.h), sides
      // "bottom" / "top". Refused edits raise ValueError with the engine's
      // message and leave the deck unchanged.
      .def("edit", [](adsim::Simulation& s, const std::string& op, const py::list& args) -> py::object {
        using namespace adsim;
        auto& d = s.Data();
        auto str = [&](size_t i) { return args[i].cast<std::string>(); };
        auto id = [&](size_t i) { return args[i].cast<size_t>(); };
        auto side = [&](size_t i) {
          const std::string v = str(i);
          if (v == "bottom") return BoundaryLocation::BOTTOM;
          if (v == "top") return BoundaryLocation::TOP;
          throw py::value_error("unknown side '" + v + "'; one of: bottom, top");
        };
        auto kind = [&](ModelFamily family, size_t i) {
          const std::string tag = str(i);
          const auto& tags = KindTags(family);
          for (size_t k = 0; k < tags.size(); ++k)
            if (tags[k] == tag) return k;
          std::string valid;
          for (const auto& t : tags) valid += (valid.empty() ? "" : ", ") + t;
          throw py::value_error(std::string("unknown ") + ModelFamilyName(family) + " kind '" + tag + "'; one of: " + valid);
        };
        std::string err;
        bool ok = true;
        py::object out = py::none();
        size_t new_id = 0;  // set by the call, read after it (not as a second argument)
        auto ret_id = [&](bool good) { ok = good; if (good) out = py::int_(new_id); };
        if (op == "add_layer") out = py::int_(AddLayer(d, str(0)));
        else if (op == "remove_layer") ok = RemoveLayer(d, id(0), err);
        else if (op == "move_layer") ok = MoveLayer(d, id(0), args[1].cast<int>(), err);
        else if (op == "add_step") out = py::int_(AddStep(d, str(0)));
        else if (op == "remove_step") ok = RemoveStep(d, id(0), err);
        else if (op == "duplicate_step") ret_id(DuplicateStep(d, id(0), new_id, err));
        else if (op == "move_step") ok = MoveStep(d, id(0), args[1].cast<int>(), err);
        else if (op == "add_stream")
          ret_id(AddStream(d, id(0), side(1), static_cast<BoundaryStreamType>(kind(ModelFamily::Stream, 2)), new_id, err));
        else if (op == "remove_stream") ok = RemoveStream(d, id(0), side(1), id(2), err);
        else if (op == "connect") ok = ConnectEndpoints(d, id(0), id(1), err);
        else if (op == "ensure_ids") EnsureIds(d);
        else if (op == "repair_connections") out = py::int_(RepairConnections(d));
        else if (op == "prune_buffer_tanks") out = py::int_(PruneBufferTanks(d));
        else if (op == "sync_kpi_config") { SyncKpiConfigToComponents(d); PruneKpiConfig(d); }
        else if (op == "fill_gaps") FillDeckGaps(d);
        else if (op == "add_experimental_data") {  // name, x, y, x_label, y_label, observable dict
          XYData xy;
          xy.name = str(0);
          xy.x = args[1].cast<std::vector<double>>();
          xy.y = args[2].cast<std::vector<double>>();
          if (args.size() > 3) xy.x_label = str(3);
          if (args.size() > 4) xy.y_label = str(4);
          if (args.size() > 5) xy.observable = BindingFromDict(args[5].cast<py::dict>());
          ret_id(AddExperimentalData(d, std::move(xy), new_id, err));
        }
        else if (op == "set_experimental_data_binding")
          ok = SetExperimentalDataBinding(d, id(0), BindingFromDict(args[1].cast<py::dict>()), err);
        else if (op == "remove_experimental_data") ok = RemoveExperimentalData(d, id(0), err);
        else if (op == "set_isotherm_kind") {
          IsothermType t{};
          if (!KindFromString(str(2), t)) throw py::value_error("unknown isotherm kind '" + str(2) + "'");
          ok = SetIsothermKind(d, id(0), str(1), t, err);
        }
        else if (op == "set_mass_transfer_kind")
          ok = SetMassTransferKind(d, id(0), str(1), static_cast<MassTransferModelType>(kind(ModelFamily::RateModel, 2)),
                                   static_cast<MassTransferCoefficientType>(kind(ModelFamily::RateCoefficient, 3)), err);
        else if (op == "set_mass_transfer_resistance" || op == "set_capacity_factor_mode") {
          const std::string where = str(2);
          if (where != "fast" && where != "slow") throw py::value_error("unknown site '" + where + "'; one of: fast, slow");
          const MassTransferSite site = where == "slow" ? MassTransferSite::Slow : MassTransferSite::Fast;
          if (op == "set_mass_transfer_resistance") {
            ok = SetMassTransferResistance(d, id(0), str(1), site, str(3), args[4].cast<bool>(), err);
          } else {
            const std::string m = str(3);
            int mode = -1;
            for (int i = 0; i < 4; ++i)
              if (m == CapacityFactorModeList[i]) mode = i;
            if (mode < 0) throw py::value_error("unknown capacity-factor mode '" + m + "'; one of: none, slope, secant, fixed");
            ok = SetCapacityFactorMode(d, id(0), str(1), site, static_cast<CapacityFactorMode>(mode), err);
          }
        }
        else if (op == "set_dispersion_kind")
          ok = SetDispersionKind(d, id(0), static_cast<AxialDispersionType>(kind(ModelFamily::Dispersion, 1)), err);
        else if (op == "set_particle_shape") ok = SetParticleShape(d, id(0), str(1), err);
        else if (op == "set_heat_dispersion_kind")
          ok = SetHeatDispersionKind(d, id(0), static_cast<AxialHeatDispersionType>(kind(ModelFamily::HeatDispersion, 1)), err);
        else if (op == "set_internal_heat_transfer_kind")
          ok = SetInternalHeatTransferKind(d, id(0),
                                           static_cast<InternalHeatTransferType>(kind(ModelFamily::InternalHeatTransfer, 1)), err);
        else if (op == "set_external_heat_transfer_kind")
          ok = SetExternalHeatTransferKind(d, id(0),
                                           static_cast<ExternalHeatTransferType>(kind(ModelFamily::ExternalHeatTransfer, 1)), err);
        else if (op == "set_particle_heat_transfer_kind")
          ok = SetParticleHeatTransferKind(d, static_cast<ParticleHeatTransferType>(kind(ModelFamily::ParticleHeatTransfer, 0)), err);
        else if (op == "set_boundary_kind")
          ok = SetBoundaryConditionKind(d, id(0), side(1), static_cast<BoundaryType>(kind(ModelFamily::Boundary, 2)), err);
        else if (op == "set_stream_kind")
          ok = SetStreamKind(d, id(0), side(1), id(2), static_cast<BoundaryStreamType>(kind(ModelFamily::Stream, 3)), err);
        else if (op == "set_flow_kind")
          ok = SetFlowKind(d, id(0), static_cast<FlowType>(kind(ModelFamily::Flow, 1)), err);
        else if (op == "set_profile_kind") {
          const std::string t = str(1);
          ProfileTarget target;
          if (t == "pressure") target = ProfileTarget::Pressure;
          else if (t == "flow") target = ProfileTarget::Flow;
          else if (t == "temperature") target = ProfileTarget::Temperature;
          else throw py::value_error("unknown profile target '" + t + "'; one of: pressure, flow, temperature");
          ok = SetProfileKind(d, id(0), target, static_cast<BoundaryValueType>(kind(ModelFamily::Profile, 2)), err);
        }
        else if (op == "set_pump_kind") ok = SetPumpKind(d, id(0), kind(ModelFamily::Pump, 1), err);
        else if (op == "set_pump_efficiency_kind")
          ok = SetPumpEfficiencyKind(d, id(0), static_cast<PumpEfficiencyType>(kind(ModelFamily::PumpEfficiency, 1)), err);
        else if (op == "set_pump_power_kind")
          ok = SetPumpPowerKind(d, id(0), static_cast<PumpPowerType>(kind(ModelFamily::PumpPower, 1)), err);
        else throw py::value_error("unknown edit operation '" + op + "'");
        s.RefreshParameters();
        if (!ok) throw py::value_error(err);
        return out;
      }, py::arg("op"), py::arg("args"))
      // The deck's structure in builder vocabulary: which model kind each
      // object uses, how the steps' ends are set up, the connections and the
      // KPI configuration. Values are read from the parameter registry.
      // The deck view (adsim/objects/deck_view.h) as JSON text.
      .def("deck_view", [](adsim::Simulation& s, bool pretty) {
        adsim::DeckViewOptions o;
        o.pretty = pretty;
        return adsim::DeckView(s.Data(), o);
      }, py::arg("pretty") = false)
      .def("describe", [](adsim::Simulation& s) {
        namespace b = adsim::build;
        using namespace adsim;
        auto& d = s.Data();
        d.cycle.EnsureIds();
        d.cycle.ResolveBufferTankReferences();
        auto kind = [](const std::vector<std::string>& kinds, int i) -> std::string {
          return i >= 0 && i < static_cast<int>(kinds.size()) ? kinds[static_cast<size_t>(i)] : std::string("?");
        };
        py::dict out;
        out["components"] = d.fluid.components_;

        py::list adsorbents;
        std::vector<AdsorbentID> ids;
        for (const auto& [id, ads] : d.adsorbent_directory.GetAll()) ids.push_back(id);
        std::sort(ids.begin(), ids.end());
        for (AdsorbentID id : ids) {
          auto& ads = d.adsorbent_directory.Get(id);
          auto& src = ads.adsorption_source_;
          py::dict a;
          a["id"] = id;
          a["name"] = ads.name_;
          std::map<std::string, std::vector<std::string>> by_kind;
          for (const auto& c : d.fluid.components_) {
            auto it = src.IsothermMap().find(c);
            if (it != src.IsothermMap().end()) by_kind[IsothermKindString(it->second)].push_back(c);
          }
          py::list isotherms;
          for (auto& [k, comps] : by_kind) {
            py::dict iso;
            iso["kind"] = k;
            iso["components"] = comps;
            isotherms.append(iso);
          }
          a["isotherms"] = isotherms;
          a["mixing_rule"] = kind(b::MixingRuleKinds(), static_cast<int>(src.mixing_rule_));
          py::dict rates;
          for (const auto& c : d.fluid.components_) {
            auto it = src.mass_transfer_models_.find(c);
            if (it == src.mass_transfer_models_.end() || !it->second || !it->second->mtc_model_) continue;
            py::dict r;
            r["model"] = kind(b::RateModelKinds(), static_cast<int>(it->second->GetModelType()));
            // The preset the switches match (a former kind's tag), or
            // series_resistance with the switches and mode (custom).
            auto site = [](adsim::IMassTransferCoefficientModel* m) {
              py::dict out;
              auto* g = dynamic_cast<adsim::MassTransferCoefficientSeriesResistance*>(m);
              if (!g) return out;
              const std::string preset = g->PresetTag();
              out["coefficient"] = preset == "custom" ? std::string("series_resistance") : preset;
              py::list on;
              for (const auto& t : adsim::ResistanceTags())
                if (g->Enabled(t)) on.append(t);
              out["resistances"] = on;
              out["capacity_factor_mode"] = std::string(adsim::CapacityFactorModeList[static_cast<int>(g->CapacityMode())]);
              return out;
            };
            py::dict fast = site(it->second->mtc_model_.get());
            std::string coefficient = "series_resistance";
            if (fast.contains("coefficient")) coefficient = fast["coefficient"].cast<std::string>();
            else if (const auto* k = adsim::KindAt(adsim::ModelFamily::RateCoefficient,
                                                   static_cast<size_t>(it->second->mtc_model_->GetModelType())))
              coefficient = k->tag;  // an older class built by hand
            r["coefficient"] = coefficient;
            // The switches and mode only where the coefficient does not imply
            // them: a custom model, or a mode the preset does not set; the
            // slow site only where it differs from the fast one.
            if (fast.contains("coefficient")) {
              const std::string mode = fast["capacity_factor_mode"].cast<std::string>();
              const std::string preset_mode = coefficient == "macropore_capacity" ? "secant" : "none";
              if (coefficient == "series_resistance") r["resistances"] = fast["resistances"];
              if (coefficient == "series_resistance" || mode != preset_mode) r["capacity_factor_mode"] = mode;
            }
            py::dict slow = site(it->second->SlowCoefficientModel());
            if (slow.contains("coefficient") && fast.contains("coefficient") &&
                (!slow["resistances"].equal(fast["resistances"]) ||
                 !slow["capacity_factor_mode"].equal(fast["capacity_factor_mode"])))
              r["slow"] = slow;
            rates[py::str(c)] = r;
          }
          a["rate_models"] = rates;
          py::list reactions;
          for (const auto& r : ads.reaction_source_.reactions_) {
            py::dict rx;
            rx["id"] = r.id_;
            rx["name"] = r.name_;
            py::dict nu;
            for (const auto& [c, v] : r.nu_) nu[py::str(c)] = v;
            rx["nu"] = nu;
            rx["kind"] = r.rate_ ? std::string(ReactionRateKindString(r.rate_->GetType())) : std::string("power_law");
            rx["basis"] = kind(b::ReactionBasisKinds(), static_cast<int>(r.basis_));
            rx["effectiveness"] = std::string(EffectivenessModeString(r.effectiveness_.mode_));
            rx["key_component"] = r.effectiveness_.key_component_;
            const auto* lhhw = r.rate_ ? dynamic_cast<const LhhwRate*>(r.rate_.get()) : nullptr;
            rx["adsorption_terms"] = lhhw ? lhhw->adsorption_.size() : size_t{0};
            reactions.append(rx);
          }
          a["reactions"] = reactions;
          adsorbents.append(a);
        }
        out["adsorbents"] = adsorbents;

        py::dict bed;
        bed["shape"] = kind(b::BedShapeKinds(), static_cast<int>(d.bed.shape_.index()));
        bed["wall_material"] = kind(b::WallMaterialKinds(), static_cast<int>(d.bed.wall_.material));
        auto& eqs = d.bed.equations_;
        bed["wall_boundary"] = eqs.wall_temperature_boundary_condition_
            ? kind(b::WallBoundaryKinds(), static_cast<int>(eqs.wall_temperature_boundary_condition_->GetModelType()))
            : std::string("zero_flux");
        bed["thermal_equilibrium"] = kind(b::ThermalEquilibriumKinds(), static_cast<int>(eqs.thermal_equilibrium_mode_));
        bed["particle_heat_transfer"] = eqs.particle_heat_transfer_model_
            ? kind(b::ParticleHeatTransferKinds(), static_cast<int>(eqs.particle_heat_transfer_model_->GetModelType()))
            : std::string("constant");
        out["bed"] = bed;

        py::list layers;
        for (auto& layer : d.bed.layers_.layers_) {
          py::dict l;
          l["id"] = layer.id_;
          l["name"] = layer.name_;
          l["momentum"] = kind(b::MomentumKinds(), static_cast<int>(layer.momentum_model_));
          l["film_correlation"] = kind(b::FilmCorrelationKinds(), static_cast<int>(layer.film_correlation_));
          auto& le = layer.equations_;
          l["dispersion"] = kind(b::DispersionKinds(), static_cast<int>(le.dispersion_model_->GetModelType()));
          l["heat_dispersion"] = kind(b::HeatDispersionKinds(), static_cast<int>(le.heat_dispersion_model_->GetModelType()));
          l["internal_heat_transfer"] =
              kind(b::InternalHeatTransferKinds(), static_cast<int>(le.internal_heat_transfer_model_->GetModelType()));
          l["external_heat_transfer"] =
              kind(b::ExternalHeatTransferKinds(), static_cast<int>(le.external_heat_transfer_model_->GetModelType()));
          py::list comp;
          // The evaluated fractions; a fraction that is a formula also has its
          // text in composition_expressions (by adsorbent id).
          py::dict comp_expr;
          for (const auto& [id, f] : layer.composition_) {
            comp.append(py::make_tuple(id, f.Get()));
            if (f.HasExpression()) comp_expr[py::int_(id)] = f.GetExpr();
          }
          l["composition"] = comp;
          l["composition_expressions"] = comp_expr;
          layers.append(l);
        }
        out["layers"] = layers;

        auto pump_of = [&](IPump* p) -> py::object {
          if (!p || !p->use_pump_) return py::none();
          py::dict q;
          q["kind"] = kind(b::PumpKinds(), static_cast<int>(p->GetType()) + 1);
          q["energy"] = p->energy_eq_ ? kind(b::PumpEnergyKinds(), static_cast<int>(p->energy_eq_->GetType()))
                                      : std::string("isothermal");
          q["efficiency"] = p->efficiency_model_
              ? kind(b::PumpEfficiencyKinds(), static_cast<int>(p->efficiency_model_->GetType()))
              : std::string("constant");
          return q;
        };
        auto tank_of = [](const BufferTankSource& t) -> py::object {
          if (!t.IsActive() || t.source_endpoint_id_ == 0) return py::none();
          return py::int_(t.source_endpoint_id_);
        };
        auto tprofile = [&](IBoundaryValue* v) -> std::string {
          return v ? kind(b::ProfileKinds(), static_cast<int>(v->GetType())) : std::string("constant");
        };
        auto end_of = [&](IBoundaryCondition* bc) -> py::dict {
          py::dict e;
          const char* dir = bc && bc->direction_ == BoundaryDirection::INLET ? "inlet" : "outlet";
          if (auto* p = dynamic_cast<PressureBoundaryCondition*>(bc)) {
            e["kind"] = "pressure";
            e["direction"] = dir;
            e["endpoint"] = p->endpoint_id_;
            e["profile"] = tprofile(p->pressure_.get());
            e["temperature_profile"] = tprofile(p->temperature_profile_.get());
            e["pump"] = pump_of(p->pump_.get());
            e["buffer_tank"] = tank_of(p->buffer_tank_);
          } else if (auto* f = dynamic_cast<FlowBoundaryCondition*>(bc)) {
            e["kind"] = "flow";
            e["direction"] = dir;
            py::list streams;
            for (auto& st : f->streams_) {
              if (!st) continue;
              py::dict sd;
              sd["name"] = st->name_;
              sd["endpoint"] = st->endpoint_id_;
              sd["temperature_profile"] = tprofile(st->temperature_profile_.get());
              sd["buffer_tank"] = tank_of(st->buffer_tank_);
              if (auto* bs = dynamic_cast<BoundaryStream*>(st.get())) {
                sd["type"] = "stream";
                sd["flow"] = kind(b::FlowKinds(), static_cast<int>(bs->boundary_flow_->GetType()));
                sd["profile"] = tprofile(bs->boundary_flow_->profile_.get());
                sd["pump"] = pump_of(bs->pump_.get());
              } else if (auto* vs = dynamic_cast<ValveBoundaryStream*>(st.get())) {
                sd["type"] = "valve";
                sd["profile"] = tprofile(vs->boundary_value_.get());
                sd["pump"] = py::none();
              } else if (auto* ps = dynamic_cast<PumpBoundaryStream*>(st.get())) {
                sd["type"] = "pump";
                sd["power"] = ps->power_ ? kind(b::PumpPowerKinds(), static_cast<int>(ps->power_->GetType()))
                                         : std::string("constant");
                sd["pump"] = pump_of(ps->pump_.get());
              } else {
                sd["type"] = "unsupported";
              }
              streams.append(sd);
            }
            e["streams"] = streams;
          } else {
            e["kind"] = "closed";
          }
          return e;
        };
        py::list steps;
        for (auto& st : d.cycle.step_list_) {
          py::dict sd;
          sd["id"] = st.id_;
          sd["name"] = st.step_name_;
          sd["bottom"] = end_of(st.bottom_boundary_condition_.get());
          sd["top"] = end_of(st.top_boundary_condition_.get());
          py::list energy;
          for (const auto& [layer_id, src] : st.energy_sources_) energy.append(layer_id);
          sd["energy_sources"] = energy;
          steps.append(sd);
        }
        out["steps"] = steps;

        // Connected pairs, by endpoint ID (each pair once).
        std::map<size_t, std::vector<size_t>> by_pair;
        for (auto& h : d.cycle.Endpoints()) {
          const Connection* c = h.stream ? &h.stream->connection_ : (h.pressure ? &h.pressure->connection_ : nullptr);
          if (c && c->IsConnected()) by_pair[c->GetID()].push_back(h.id);
        }
        py::list connections;
        for (auto& [pair, eps] : by_pair)
          if (eps.size() == 2) connections.append(py::make_tuple(eps[0], eps[1]));
        out["connections"] = connections;

        py::list kpis;
        for (auto& cfg : d.kpi_config.per_component) {
          if (!cfg.tracked) continue;
          auto ids_of = [&](std::vector<KPIEndpointRef> refs) {
            std::vector<size_t> v;
            for (auto& r : refs)
              if (r.Resolve(d.cycle.step_list_)) v.push_back(r.endpoint_id);
            return v;
          };
          py::dict k;
          k["component"] = cfg.product_component;
          k["feed"] = ids_of(cfg.feed_endpoints);
          k["product"] = ids_of(cfg.product_endpoints);
          k["recycle"] = ids_of(cfg.recycle_endpoints);
          kpis.append(k);
        }
        out["kpis"] = kpis;
        {
          const auto& rx = d.kpi_config.reaction;
          py::dict r;
          r["key_reactant"] = rx.key_reactant;
          r["products"] = rx.products;
          r["basis"] = rx.basis == SpaceTimeYieldBasis::BED_VOLUME ? "bed_volume" : "catalyst_mass";
          out["reaction_kpi"] = r;
        }

        py::list variables;
        for (const auto& v : d.variable_registry.GetAll()) {
          py::dict vd;
          vd["name"] = v.name;
          vd["category"] = b::VariableCategoryKind(v.unit_category);
          vd["formula"] = v.value_expr;
          vd["description"] = v.description;
          variables.append(vd);
        }
        out["variables"] = variables;
        return out;
      });

  // Every model kind's descriptor (adsim/objects/model_kinds.h): family ->
  // [{tag, label, parameters: [{name, unit, default, lower, upper, tunable,
  // label}]}], isotherms under "isotherm" with their extra flags.
  m.def("model_kinds", []() {
    auto params = [](const adsim::KindDescriptor& k) {
      py::list ps;
      for (const auto& p : k.parameters) {
        py::dict d;
        d["name"] = p.name;
        d["unit"] = p.unit;
        d["default"] = p.default_value;
        d["lower"] = p.lower;
        d["upper"] = p.upper;
        d["tunable"] = p.tunable;
        d["label"] = p.Label();
        ps.append(d);
      }
      return ps;
    };
    // The families from the engine's manifest (one definition: legacy flags,
    // the generic coefficient's resistances, capacity modes and presets).
    py::dict out = py::module_::import("json").attr("loads")(adsim::ModelKindManifest());
    py::list isotherms;
    for (const auto& k : adsim::IsothermKindDescriptors()) {
      py::dict d;
      d["tag"] = k.tag;
      d["label"] = k.label;
      d["parameters"] = params(k);
      d["competitive"] = k.competitive;
      d["water_coupled"] = k.water_coupled;
      d["fittable"] = k.fittable;
      d["iast"] = k.iast != adsim::IastFamily::None;
      d["paired"] = k.has_pair ? py::object(py::str(adsim::IsothermKindString(k.pair))) : py::object(py::none());
      isotherms.append(d);
    }
    out["isotherm"] = isotherms;
    return out;
  });

  m.def("build_kinds", []() {
    namespace b = adsim::build;
    py::dict k;
    k["isotherm"] = b::IsothermKinds();
    k["variable_category"] = b::VariableCategoryKinds();
    k["rate_model"] = b::RateModelKinds();
    k["rate_coefficient"] = b::RateCoefficientKinds();
    k["rate_coefficient_preset"] = b::RateCoefficientPresets();
    k["mixing_rule"] = b::MixingRuleKinds();
    k["reaction_rate"] = b::ReactionRateKinds();
    k["reaction_basis"] = b::ReactionBasisKinds();
    k["effectiveness"] = b::EffectivenessKinds();
    k["momentum"] = b::MomentumKinds();
    k["film_correlation"] = b::FilmCorrelationKinds();
    k["particle_shape"] = adsim::ParticleShapeKinds();
    k["dispersion"] = b::DispersionKinds();
    k["heat_dispersion"] = b::HeatDispersionKinds();
    k["internal_heat_transfer"] = b::InternalHeatTransferKinds();
    k["external_heat_transfer"] = b::ExternalHeatTransferKinds();
    k["particle_heat_transfer"] = b::ParticleHeatTransferKinds();
    k["thermal_equilibrium"] = b::ThermalEquilibriumKinds();
    k["bed_shape"] = b::BedShapeKinds();
    k["wall_material"] = b::WallMaterialKinds();
    k["wall_boundary"] = b::WallBoundaryKinds();
    k["profile"] = b::ProfileKinds();
    k["flow"] = b::FlowKinds();
    k["pump"] = b::PumpKinds();
    k["pump_energy"] = b::PumpEnergyKinds();
    k["pump_efficiency"] = b::PumpEfficiencyKinds();
    k["pump_power"] = b::PumpPowerKinds();
    return k;
  });

  // Units: the engine's unit parser (adsim/math/unit_parser.h), so Python
  // accepts exactly the units decks and the CLI do.
  m.def("parse_unit", [](const std::string& text) {
    adsim::UnitInfo u;
    std::string why;
    if (!adsim::ParseUnitExpression(text, u, &why)) throw py::value_error("unknown unit '" + text + "': " + why);
    py::dict d;
    d["factor"] = u.factor;
    d["offset"] = u.offset;
    d["dimension"] = py::make_tuple(u.dimension.length, u.dimension.mass, u.dimension.time,
                                    u.dimension.temperature, u.dimension.amount);
    d["dimension_name"] = adsim::DimensionName(u.dimension);
    return d;
  }, py::arg("unit"));
  m.def("convert_unit", &adsim::ConvertUnit, py::arg("value"), py::arg("from_unit"), py::arg("to_unit"));
  m.def("split_quantity", [](const std::string& text) -> py::object {
    double number = 0.0;
    std::string unit;
    if (!adsim::SplitQuantity(text, number, unit)) return py::none();
    return py::make_tuple(number, unit);
  }, py::arg("text"));

  // Deck view schema and display manifests (adsim/objects/manifests.h), as JSON text.
  m.def("deck_view_schema", []() { return adsim::DeckViewSchema(); });
  m.def("kpi_manifest", []() { return adsim::KpiManifest(); });
  m.def("unit_manifest", []() { return adsim::UnitManifest(); });
  m.def("warning_manifest", []() { return adsim::WarningManifest(); });
  m.def("component_manifest", []() { return adsim::ComponentManifest(); });
  m.def("series_manifest", []() { return adsim::SeriesManifest(); });
  m.def("example_deck_manifest", [](const std::string& manifest, const std::string& decks) {
    return adsim::ExampleDeckManifest(manifest, decks);
  }, py::arg("manifest_path"), py::arg("decks_dir"));
  m.def("check_unit", [](const std::string& text, const std::string& category) {
    const adsim::UnitCheck c = adsim::CheckUnit(text, category);
    py::dict d;
    d["ok"] = c.ok;
    d["category"] = c.category;
    d["dimension"] = c.dimension;
    d["factor"] = c.factor;
    d["offset"] = c.offset;
    d["message"] = c.message;
    return d;
  }, py::arg("unit"), py::arg("category") = "");

  m.def("bulk_variables", &adsim::BulkVariableNames);
  m.def("component_variables", &adsim::ComponentVariableNames);
  m.def("stream_variables", &adsim::StreamVariableNames);
  m.def("stream_component_variables", &adsim::StreamComponentVariableNames);
  m.def("per_adsorbent_variables", &adsim::PerAdsorbentVariableNames);
  m.def("reaction_variables", &adsim::ReactionVariableNames);

  // ---- gas properties -------------------------------------------------------
  // Ideal-gas heat capacity C_p0(T) [J/(mol K)] of the NASA polynomials
  // (ideal_gas_heat_capacity.h), as the default gas heat capacity mode uses it.
  m.def("ideal_gas_heat_capacity", [](const std::string& component, double temperature) {
    const adsim::NasaPolynomial* p = adsim::Fluid::IdealGasCoeffs(component);
    if (!p) throw py::value_error("no ideal-gas heat capacity polynomial for '" + component + "'");
    return adsim::Fluid::EvalIdealGasHeatCapacity(*p, temperature);
  }, py::arg("component"), py::arg("temperature"));
  // Provenance of the polynomials: the table's reference and, per component, the
  // TM-4513 species label and data-source code (or Air's mixing note).
  m.def("ideal_gas_heat_capacity_reference", []() { return std::string(adsim::kIdealGasHeatCapacityReference); });
  m.def("ideal_gas_heat_capacity_source", [](const std::string& component) {
    const adsim::NasaPolynomial* p = adsim::Fluid::IdealGasCoeffs(component);
    if (!p) throw py::value_error("no ideal-gas heat capacity polynomial for '" + component + "'");
    py::dict d;
    d["species"] = std::string(p->species);
    d["source"] = std::string(p->source);
    return d;
  }, py::arg("component"));
  m.def("components_with_held_heat_capacity", []() {
    std::vector<std::string> out;
    for (const auto& c : adsim::Components::componentNames)
      if (!adsim::Fluid::IdealGasCoeffs(c)) out.push_back(c);
    return out;
  });

  // ---- isotherms -------------------------------------------------------------
  m.def("isotherm_kinds", []() {
    std::vector<std::string> out;
    for (auto t : adsim::isofit::FittableKinds()) out.emplace_back(adsim::IsothermKindString(t));
    return out;
  });
  m.def("isotherm_parameters", [](const std::string& kind) {
    adsim::isofit::Model model(ParseIsothermKind(kind));
    if (!model.valid()) throw py::value_error("cannot build an isotherm of kind '" + kind + "'");
    return model.Names();
  }, py::arg("kind"));
  // `component`: the fluid (CoolProp name) whose saturation pressure the GAB
  // kinds use; empty means water. Other kinds ignore it. `water_loading`
  // [mol/kg]: the adsorbed water loading the water-coupled kinds read
  // (empty: 0, dry).
  m.def("evaluate_isotherm", [](const std::string& kind, const std::map<std::string, double>& params,
                                const std::vector<double>& pressure, const std::vector<double>& temperature,
                                const std::vector<double>& mole_fraction, const std::string& component,
                                const std::vector<double>& water_loading) {
    const auto pts = MakePoints(pressure, {}, temperature, mole_fraction, water_loading);
    std::vector<std::pair<std::string, double>> values(params.begin(), params.end());
    std::vector<double> q;
    std::string err;
    if (!adsim::isofit::Evaluate(ParseIsothermKind(kind), values, pts, q, err, component)) throw py::value_error(err);
    return ToArray(q);
  }, py::arg("kind"), py::arg("params"), py::arg("pressure"), py::arg("temperature"),
     py::arg("mole_fraction") = std::vector<double>{}, py::arg("component") = std::string(),
     py::arg("water_loading") = std::vector<double>{});
  // IAST mixture loadings from several pure (kind, params) isotherms at one
  // (T, P, y): the engine's adsim/equations/iast.h.
  m.def("evaluate_iast", [](const std::vector<std::string>& kinds,
                            const std::vector<std::map<std::string, double>>& params,
                            double pressure, double temperature, const std::vector<double>& mole_fraction,
                            const std::vector<std::string>& components) {
    if (kinds.size() != params.size() || kinds.size() != mole_fraction.size())
      throw py::value_error("kinds, params and mole_fraction must have one entry per component");
    if (!components.empty() && components.size() != kinds.size())
      throw py::value_error("components must be empty or have one entry per component");
    std::vector<std::unique_ptr<adsim::IIsotherm>> owned;
    std::vector<adsim::IIsotherm*> isos;
    std::vector<std::string> names;
    for (size_t i = 0; i < kinds.size(); ++i) {
      // The isotherm's component name: the GAB kinds look up its saturation
      // pressure (water unless `components` names the fluid).
      const auto type = ParseIsothermKind(kinds[i]);
      const bool gab = type == adsim::IsothermType::GAB || type == adsim::IsothermType::GAB_YOUNG;
      const std::string name = !components.empty() && !components[i].empty() ? components[i]
                               : gab ? std::string("Water") : "component_" + std::to_string(i);
      auto iso = adsim::CreateIsotherm(type, name);
      if (!iso) throw py::value_error("cannot build an isotherm of kind '" + kinds[i] + "'");
      auto fields = adsim::IsothermFields(iso.get(), name);
      for (const auto& [nm, v] : params[i]) {
        bool found = false;
        for (auto& f : fields) if (f.first == nm) { f.second->Get() = v; found = true; }
        if (!found) throw py::value_error("unknown parameter '" + nm + "' for isotherm kind '" + kinds[i] + "'");
      }
      isos.push_back(iso.get());
      names.push_back(name);
      owned.push_back(std::move(iso));
    }
    std::vector<double> q;
    const auto r = adsim::EvalIastMixture(isos, names, temperature, pressure, mole_fraction, q);
    if (!r.converged) throw py::value_error("the IAST solve did not converge at this state");
    return ToArray(q);
  }, py::arg("kinds"), py::arg("params"), py::arg("pressure"), py::arg("temperature"), py::arg("mole_fraction"),
     py::arg("components") = std::vector<std::string>{});
  m.def("default_fit_specs", [](const std::string& kind, const std::vector<double>& pressure,
                                const std::vector<double>& loading, const std::vector<double>& temperature,
                                const std::vector<double>& water_loading) {
    const auto pts = MakePoints(pressure, loading, temperature, {}, water_loading);
    std::vector<adsim::isofit::ParamSpec> specs;
    double t_ref = 0.0;
    std::string err;
    if (!adsim::isofit::DefaultSpecs(ParseIsothermKind(kind), pts, specs, t_ref, err))
      throw py::value_error(err);
    py::list out;
    for (const auto& s : specs) out.append(SpecToDict(s));
    return py::make_tuple(out, t_ref);
  }, py::arg("kind"), py::arg("pressure"), py::arg("loading"), py::arg("temperature"),
     py::arg("water_loading") = std::vector<double>{});
  m.def("fit_isotherm", [](const std::string& kind, const std::vector<double>& pressure,
                           const std::vector<double>& loading, const std::vector<double>& temperature,
                           const std::vector<double>& mole_fraction, py::object specs_in, double t_ref,
                           int max_global_evals, const std::string& component,
                           const std::vector<double>& water_loading) {
    const auto type = ParseIsothermKind(kind);
    const auto pts = MakePoints(pressure, loading, temperature, mole_fraction, water_loading);
    std::vector<adsim::isofit::ParamSpec> specs;
    if (specs_in.is_none()) {
      std::string err;
      double tr = 0.0;
      if (!adsim::isofit::DefaultSpecs(type, pts, specs, tr, err)) throw py::value_error(err);
      if (!(t_ref > 0.0)) t_ref = tr;
    } else {
      for (auto item : specs_in.cast<py::list>()) specs.push_back(SpecFromDict(item.cast<py::dict>()));
    }
    adsim::isofit::FitOptions opt;
    opt.max_global_evals = max_global_evals;
    opt.component = component;
    adsim::isofit::FitResult r;
    {
      py::gil_scoped_release release;
      r = adsim::isofit::Fit(type, pts, specs, t_ref, opt);
    }
    if (!r.ok) throw py::value_error("isotherm fit failed: " + r.error);
    py::dict d;
    py::dict params;
    py::dict errors;
    py::list fixed;
    for (const auto& p : r.params) {
      params[py::str(p.name)] = p.value;
      errors[py::str(p.name)] = p.std_error;
      if (p.fixed) fixed.append(p.name);
    }
    d["kind"] = kind;
    d["params"] = params;
    d["std_errors"] = errors;
    d["fixed"] = fixed;
    d["sse"] = r.sse;
    d["rmse"] = r.rmse;
    d["r2"] = r.r2;
    d["n"] = r.n;
    d["dof"] = r.dof;
    d["evaluations"] = r.evaluations;
    d["t_ref"] = r.t_ref;
    d["warning"] = r.warning;
    d["q_fit"] = ToArray(r.q_fit);
    return d;
  }, py::arg("kind"), py::arg("pressure"), py::arg("loading"), py::arg("temperature"),
     py::arg("mole_fraction") = std::vector<double>{}, py::arg("specs") = py::none(), py::arg("t_ref") = 0.0,
     py::arg("max_global_evals") = 0, py::arg("component") = std::string(),
     py::arg("water_loading") = std::vector<double>{});

  m.def("sensitivity_morris_plan", [](const py::list& inputs, int trajectories, int levels, std::uint64_t seed) {
    return PlanToDict(sa::MorrisPlan(InputsFromList(inputs), trajectories, levels, seed));
  }, py::arg("inputs"), py::arg("trajectories"), py::arg("levels"), py::arg("seed"));
  m.def("sensitivity_monte_carlo_plan", [](const py::list& inputs, int runs, std::uint64_t seed, int batch) {
    return PlanToDict(sa::MonteCarloPlan(InputsFromList(inputs), runs, seed, batch));
  }, py::arg("inputs"), py::arg("runs"), py::arg("seed"), py::arg("batch"));
  m.def("sensitivity_quantile", [](const py::dict& input, const std::vector<double>& u) {
    const sa::Input in = InputFromDict(input);
    in.Validate();
    std::vector<double> out;
    for (double x : u) out.push_back(in.Quantile(x));
    return ToArray(out);
  }, py::arg("input"), py::arg("u"));
  m.def("sensitivity_kpi_names", &sa::KPINames);
  m.def("step_output_names", &adsim::StepOutputNames);
  // ---- warm-start cache and optimisation (state_cache.h, optimisation.h) ----
  py::class_<adsim::CycleStateCache>(m, "StateCache")
      .def(py::init<size_t>(), py::arg("capacity") = adsim::CycleStateCache::kDefaultCapacity)
      .def(py::init<std::vector<double>, std::vector<double>, size_t>(), py::arg("lower"), py::arg("upper"),
           py::arg("capacity") = adsim::CycleStateCache::kDefaultCapacity)
      .def("normalise", &adsim::CycleStateCache::Normalise, py::arg("x"))
      .def("add", [](adsim::CycleStateCache& c, const std::vector<double>& x, const adsim::CycleState& state) {
        c.Add(x, state);
      }, py::arg("x"), py::arg("state"))
      .def("nearest", [](const adsim::CycleStateCache& c, const std::vector<double>& x) -> py::object {
        double d = 0.0;
        const adsim::CycleState* s = c.Nearest(x, &d);
        if (!s) return py::none();
        return py::make_tuple(*s, d);
      }, py::arg("x"))
      .def("clear", &adsim::CycleStateCache::Clear)
      .def("__len__", &adsim::CycleStateCache::Size)
      .def_property_readonly("size", &adsim::CycleStateCache::Size)
      .def_property_readonly("capacity", &adsim::CycleStateCache::Capacity)
      .def_property_readonly("dimension", &adsim::CycleStateCache::Dimension)
      .def_property_readonly("lower", &adsim::CycleStateCache::Lower)
      .def_property_readonly("upper", &adsim::CycleStateCache::Upper);

  m.def("run_with_state_cache", [](adsim::Simulation& s, adsim::CycleStateCache* cache, const std::vector<double>& x,
                                   std::optional<int> cycles, py::object progress, const std::string& log_level,
                                   py::object on_message) {
    adsim::CachedRunInfo info;
    adsim::RunResult r = RunFromPython(cycles, progress, log_level, on_message,
                                       [&](const adsim::RunOptions& o) {
                                         return adsim::RunWithStateCache(s, cache, x, o, &info);
                                       });
    return py::make_tuple(r, CachedRunInfoToDict(info));
  }, py::arg("sim"), py::arg("cache"), py::arg("x"), py::arg("cycles") = py::none(),
     py::arg("progress") = py::none(), py::arg("log_level") = "low", py::arg("on_message") = py::none());

  m.def("optimise", [](const std::string& deck_json, const py::list& variables, const py::list& objectives,
                       const py::list& constraints, const py::dict& options, py::object on_evaluation) {
    std::vector<adsim::OptimizationVariable> vars;
    for (auto item : variables) {
      auto t = item.cast<py::tuple>();
      vars.push_back({t[0].cast<std::string>(), t[1].cast<double>(), t[2].cast<double>()});
    }
    adsim::OptimizationCostFunction cost;
    for (auto item : objectives) {
      auto d = item.cast<py::dict>();
      const std::string dir = d["direction"].cast<std::string>();
      if (dir != "min" && dir != "max") throw std::invalid_argument("an objective's direction is 'min' or 'max'");
      cost.AddObjective(d["component"].cast<std::string>(), KPITypeFromName(d["kpi"].cast<std::string>()),
                        d["weight"].cast<double>(),
                        dir == "min" ? adsim::ObjectiveDirection::MINIMIZE : adsim::ObjectiveDirection::MAXIMIZE);
    }
    for (auto item : constraints) {
      auto d = item.cast<py::dict>();
      const std::string sense = d["sense"].cast<std::string>();
      if (sense != ">=" && sense != "<=") throw std::invalid_argument("a constraint's sense is '>=' or '<='");
      cost.AddConstraint(d["component"].cast<std::string>(), KPITypeFromName(d["kpi"].cast<std::string>()),
                         d["target"].cast<double>(), d["penalty"].cast<double>(),
                         sense == ">=" ? adsim::ConstraintSense::GREATER_OR_EQUAL
                                       : adsim::ConstraintSense::LESS_OR_EQUAL);
    }
    adsim::SimulationOptimizerSettings st;
    if (options.contains("evaluations")) {
      const int total = options["evaluations"].cast<int>();
      st.n_init_samples = options.contains("init_samples") ? options["init_samples"].cast<int>() : 5;
      st.n_init_samples = std::max(1, std::min(st.n_init_samples, total));
      st.n_iterations = std::max(0, total - st.n_init_samples);
    }
    if (options.contains("warm_start")) st.warm_start = options["warm_start"].cast<bool>();
    if (options.contains("rerun_cold")) st.rerun_cold = options["rerun_cold"].cast<bool>();
    if (options.contains("require_css")) st.require_css = options["require_css"].cast<bool>();
    if (options.contains("cycles") && !options["cycles"].is_none()) st.cycles = options["cycles"].cast<int>();
    if (options.contains("cache_capacity")) st.cache_capacity = options["cache_capacity"].cast<size_t>();
    if (options.contains("settings")) {
      for (auto item : options["settings"].cast<py::dict>())
        st.settings.emplace_back(item.first.cast<std::string>(), py::str(item.second).cast<std::string>());
    }
    std::optional<py::error_already_set> pending;
    auto last_check = std::chrono::steady_clock::now();
    st.cancel = [&]() {
      if (pending) return true;
      const auto now = std::chrono::steady_clock::now();
      if (now - last_check < std::chrono::milliseconds(100)) return false;
      last_check = now;
      py::gil_scoped_acquire gil;
      if (PyErr_CheckSignals() != 0) {
        pending.emplace();
        return true;
      }
      return false;
    };
    if (!on_evaluation.is_none()) {
      st.on_evaluation = [&](const adsim::OptimizationEvaluation& e) {
        py::gil_scoped_acquire gil;
        if (pending) return;
        try {
          on_evaluation(EvaluationToDict(e));
        } catch (py::error_already_set& err) {
          pending.emplace(std::move(err));
        }
      };
    }
    adsim::CheckOptimization(deck_json, vars, cost, st);
    adsim::SimulationOptimizationResult r;
    {
      py::gil_scoped_release release;
      r = adsim::OptimizeSimulation(deck_json, vars, cost, st);
    }
    if (pending) throw std::move(*pending);
    py::dict d;
    d["x"] = r.x;
    d["cost"] = r.cost;
    d["search_cost"] = r.search_cost;
    d["best_index"] = r.best_index;
    d["rerun_cold"] = r.rerun_cold;
    d["cold"] = r.rerun_cold ? py::object(EvaluationToDict(r.cold)) : py::object(py::none());
    py::list evs;
    for (const auto& e : r.evaluations) evs.append(EvaluationToDict(e));
    d["evaluations"] = evs;
    d["cancelled"] = r.cancelled;
    return d;
  }, py::arg("deck_json"), py::arg("variables"), py::arg("objectives"), py::arg("constraints"),
     py::arg("options"), py::arg("on_evaluation") = py::none());

  m.def("sensitivity_default_outputs", [](const std::string& deck_json) {
    adsim::Simulation sim = adsim::Simulation::FromJson(deck_json);
    std::vector<std::string> out;
    for (const auto& o : sa::DefaultOutputs(sim.Data())) out.push_back(o.Name());
    return out;
  }, py::arg("deck_json"));
  m.def("sensitivity_check", [](const std::string& deck_json, const py::list& inputs,
                                const std::vector<std::string>& outputs, const py::dict& policy) {
    sa::CheckStudy(deck_json, InputsFromList(inputs), OutputsFromList(outputs), PolicyFromDict(policy));
  }, py::arg("deck_json"), py::arg("inputs"), py::arg("outputs"), py::arg("policy"));
  // One sample, GIL released; Ctrl-C cancels the run and is re-raised.
  m.def("sensitivity_evaluate", [](const std::string& deck_json, const py::list& inputs,
                                   const std::vector<double>& values, const std::vector<std::string>& outputs,
                                   const py::dict& policy, adsim::CycleStateCache* cache,
                                   const std::vector<double>& key) {
    const auto in = InputsFromList(inputs);
    const auto out = OutputsFromList(outputs);
    const auto pol = PolicyFromDict(policy);
    std::optional<py::error_already_set> pending;
    auto last_check = std::chrono::steady_clock::now();
    adsim::CancelPoll cancel = [&]() {
      if (pending) return true;
      const auto now = std::chrono::steady_clock::now();
      if (now - last_check < std::chrono::milliseconds(100)) return false;
      last_check = now;
      py::gil_scoped_acquire gil;
      if (PyErr_CheckSignals() != 0) {
        pending.emplace();
        return true;
      }
      return false;
    };
    sa::SampleOutcome r;
    {
      py::gil_scoped_release release;
      r = sa::EvaluateSample(deck_json, in, values, out, pol, cancel, cache, key);
    }
    if (pending) throw std::move(*pending);
    py::dict d;
    d["status"] = adsim::RunStatusName(r.status);
    d["message"] = r.message;
    d["accepted"] = r.accepted;
    d["outputs"] = r.outputs;
    d["wall_time_s"] = r.wall_time_s;
    d["css_reached"] = r.css_reached;
    d["cycles_run"] = r.cycles_run;
    d["warm"] = r.warm;
    d["distance"] = r.distance;
    d["note"] = r.note;
    return d;
  }, py::arg("deck_json"), py::arg("inputs"), py::arg("values"), py::arg("outputs"), py::arg("policy"),
     py::arg("cache") = nullptr, py::arg("key") = std::vector<double>{});
  m.def("sensitivity_analyse_morris", [](const py::list& inputs, const py::dict& plan_in,
                                         const std::vector<std::string>& names,
                                         const std::vector<std::vector<double>>& y, int bootstrap) {
    sa::Plan plan;
    plan.method = sa::Method::Morris;
    plan.inputs = InputsFromList(inputs);
    plan.seed = plan_in["seed"].cast<std::uint64_t>();
    plan.trajectories = plan_in["trajectories"].cast<int>();
    plan.levels = plan_in["levels"].cast<int>();
    plan.changed = plan_in["changed"].cast<std::vector<int>>();
    plan.u = plan_in["u"].cast<std::vector<std::vector<double>>>();
    plan.values = plan_in["values"].cast<std::vector<std::vector<double>>>();
    const sa::MorrisResult r = sa::AnalyseMorris(plan, names, y, bootstrap);
    py::list outputs;
    for (const auto& o : r.outputs) {
      py::list effects;
      for (const auto& e : o.effects) {
        py::dict d;
        d["input"] = e.input;
        d["mu"] = e.mu;
        d["mu_star"] = e.mu_star;
        d["sigma"] = e.sigma;
        d["mu_star_interval"] = IntervalToDict(e.mu_star_interval);
        d["effects"] = e.effects;
        effects.append(d);
      }
      py::dict d;
      d["output"] = o.output;
      d["missing"] = o.missing;
      d["effects"] = effects;
      outputs.append(d);
    }
    py::dict d;
    d["runs"] = r.runs;
    d["outputs"] = outputs;
    return d;
  }, py::arg("inputs"), py::arg("plan"), py::arg("names"), py::arg("y"), py::arg("bootstrap") = 1000);
  m.def("sensitivity_analyse_monte_carlo", [](const py::list& inputs, const std::vector<std::vector<double>>& values,
                                              const std::vector<std::string>& names,
                                              const std::vector<std::vector<double>>& y, std::uint64_t seed,
                                              int bootstrap, double rel_tol) {
    const sa::MonteCarloResult r =
        sa::AnalyseMonteCarlo(InputsFromList(inputs), values, names, y, seed, bootstrap, rel_tol);
    py::list outputs;
    for (const auto& s : r.outputs) {
      py::list shares;
      for (const auto& x : s.shares) {
        py::dict d;
        d["input"] = x.input;
        d["spearman"] = x.spearman;
        d["srrc"] = x.srrc;
        d["share"] = x.share;
        shares.append(d);
      }
      py::dict d;
      d["output"] = s.output;
      d["valid"] = s.valid;
      d["missing"] = s.missing;
      d["mean"] = s.mean;
      d["sd"] = s.sd;
      d["min"] = s.min;
      d["max"] = s.max;
      d["p05"] = s.p05;
      d["p50"] = s.p50;
      d["p95"] = s.p95;
      d["p05_interval"] = IntervalToDict(s.p05_interval);
      d["p50_interval"] = IntervalToDict(s.p50_interval);
      d["p95_interval"] = IntervalToDict(s.p95_interval);
      d["converged"] = s.converged;
      d["r2"] = s.r2;
      d["shares"] = shares;
      outputs.append(d);
    }
    py::dict d;
    d["runs"] = r.runs;
    d["converged"] = r.converged;
    d["outputs"] = outputs;
    return d;
  }, py::arg("inputs"), py::arg("values"), py::arg("names"), py::arg("y"), py::arg("seed"),
     py::arg("bootstrap") = 1000, py::arg("rel_tol") = 0.1);

  // ---- calibration ---------------------------------------------------------------
  py::class_<CalibrationCancel>(m, "CalibrationCancel")
      .def(py::init<>())
      .def("cancel", &CalibrationCancel::Cancel)
      .def("reset", &CalibrationCancel::Reset)
      .def_property_readonly("cancelled", &CalibrationCancel::Cancelled)
      .def("__repr__", [](const CalibrationCancel& c) {
        return std::string("<CalibrationCancel ") + (c.Cancelled() ? "cancelled>" : "armed>");
      });
  // Validates the problem (no runs) and returns its unknowns, in order.
  m.def("calibration_check", [](const py::list& cases, const py::list& parameters) {
    return UnknownsToList(ProblemFromLists(cases, parameters));
  }, py::arg("cases"), py::arg("parameters"));
  // A blocking fit. The GIL is released while the engine's threads run;
  // `progress(dict)` re-acquires it (returning False cancels). Ctrl-C, or
  // `cancel.cancel()` from any thread, stops the fit at its best point so far;
  // a KeyboardInterrupt (or an error raised by `progress`) is re-raised after.
  m.def("calibrate", [](const py::list& cases, const py::list& parameters, const py::dict& options,
                        py::object progress, CalibrationCancel* token) {
    const cal::Problem problem = ProblemFromLists(cases, parameters);
    cal::Options opt = CalOptionsFromDict(options);
    std::atomic<bool> stop{false};
    std::optional<py::error_already_set> pending;  // touched with the GIL held
    const auto caller = std::this_thread::get_id();
    auto last_check = std::chrono::steady_clock::now();  // the calling thread only
    opt.cancel = [&]() {
      if (stop.load() || (token != nullptr && token->Cancelled())) return true;
      // Python handles signals only on its main thread: poll from the thread
      // that called (it takes a share of every batch's runs).
      if (std::this_thread::get_id() != caller) return false;
      const auto now = std::chrono::steady_clock::now();
      if (now - last_check < std::chrono::milliseconds(100)) return false;
      last_check = now;
      py::gil_scoped_acquire gil;
      if (PyErr_CheckSignals() != 0) {
        if (!pending) pending.emplace();
        else PyErr_Clear();
        stop = true;
        return true;
      }
      return false;
    };
    if (!progress.is_none()) {
      opt.progress = [&](const cal::Progress& p) {
        py::gil_scoped_acquire gil;
        if (pending) return;
        try {
          py::dict d;
          d["solver"] = cal::SolverName(p.solver);
          d["iteration"] = p.iteration;
          d["objective"] = p.objective;
          d["values"] = ToArray(p.values);
          d["simulations"] = p.simulations;
          py::object r = progress(d);
          if (!r.is_none() && !r.cast<bool>()) stop = true;
        } catch (py::error_already_set& e) {
          pending.emplace(std::move(e));
          stop = true;
        }
      };
    }
    cal::Result r;
    {
      py::gil_scoped_release release;
      r = cal::Calibrate(problem, opt);
    }
    if (pending) throw std::move(*pending);
    return CalResultToDict(r);
  }, py::arg("cases"), py::arg("parameters"), py::arg("options"), py::arg("progress") = py::none(),
     py::arg("cancel") = nullptr);
  // The objective and residuals at a point (per unknown, SI units): one run per case.
  m.def("calibration_evaluate", [](const py::list& cases, const py::list& parameters,
                                   const std::vector<double>& values, const py::dict& options) {
    const cal::Problem problem = ProblemFromLists(cases, parameters);
    const cal::Options opt = CalOptionsFromDict(options);
    cal::Evaluation e;
    {
      py::gil_scoped_release release;
      e = cal::Evaluate(problem, values, opt);
    }
    py::dict d;
    d["ok"] = e.ok;
    d["objective"] = e.objective;
    d["residuals"] = ToArray(e.residuals);
    d["message"] = e.message;
    return d;
  }, py::arg("cases"), py::arg("parameters"), py::arg("values"), py::arg("options"));
}
