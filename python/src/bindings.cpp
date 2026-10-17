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
#include <optional>
#include <string>
#include <vector>

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "adsim/simulation.h"
#include "adsim/version.h"

#include "adsim/equations/isotherm_kinds.h"
#include "adsim/math/unit_parser.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/deck_edit.h"
#include "adsim/solver/isotherm_fit.h"

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
                                                  const std::vector<double>& T, const std::vector<double>& y) {
  const size_t n = P.size();
  if (!q.empty() && q.size() != n) throw py::value_error("pressure and loading must have the same length");
  if (T.size() != n && T.size() != 1) throw py::value_error("temperature must be a scalar or match pressure");
  if (!y.empty() && y.size() != n && y.size() != 1)
    throw py::value_error("mole_fraction must be a scalar or match pressure");
  std::vector<adsim::isofit::Point> pts(n);
  for (size_t i = 0; i < n; ++i) {
    pts[i].P = P[i];
    pts[i].q = q.empty() ? 0.0 : q[i];
    pts[i].T = T.size() == 1 ? T[0] : T[i];
    pts[i].y = y.empty() ? 1.0 : (y.size() == 1 ? y[0] : y[i]);
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
  return s;
}

}  // namespace

PYBIND11_MODULE(_core, m) {
  m.doc() = "Skarstrom PSA/TSA simulation engine (compiled core). Use the 'skarstrom' package.";
  m.attr("__version__") = ADSIM_VERSION_STRING;

  // ---- run status / records ---------------------------------------------
  py::class_<adsim::RunResult>(m, "RunResult")
      .def_property_readonly("status", [](const adsim::RunResult& r) { return adsim::RunStatusName(r.status); })
      .def_readonly("message", &adsim::RunResult::message)
      .def_readonly("return_code", &adsim::RunResult::return_code)
      .def_readonly("css_reached", &adsim::RunResult::css_reached)
      .def_readonly("css_cycle", &adsim::RunResult::css_cycle)
      .def_readonly("cycles_run", &adsim::RunResult::cycles_run)
      .def_readonly("wall_time_s", &adsim::RunResult::wall_time_s)
      .def_property_readonly("ok", &adsim::RunResult::Ok)
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
        return d;
      })
      .def("__repr__", [](const adsim::RunResult& r) {
        return "<RunResult " + std::string(adsim::RunStatusName(r.status)) + ": " + r.message + ">";
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
      .def("__repr__", [](const adsim::CSSRecord& c) {
        return "<CSSRecord cycle=" + std::to_string(c.cycle) + (c.converged ? " converged" : "") +
               " error=" + std::to_string(c.error) + ">";
      });

  py::class_<adsim::TimeSeries>(m, "TimeSeries")
      .def_readonly("variable", &adsim::TimeSeries::variable)
      .def_readonly("component", &adsim::TimeSeries::component)
      .def_readonly("unit", &adsim::TimeSeries::unit)
      .def_property_readonly("time", [](const adsim::TimeSeries& t) { return ToArray(t.time); })
      .def_property_readonly("values", [](const adsim::TimeSeries& t) { return ToArray(t.values); })
      .def("__len__", [](const adsim::TimeSeries& t) { return t.time.size(); })
      .def("__repr__", [](const adsim::TimeSeries& t) {
        return "<TimeSeries " + t.variable + (t.component.empty() ? "" : "[" + t.component + "]") + " [" +
               t.unit + "], " + std::to_string(t.time.size()) + " points>";
      });

  py::class_<adsim::AxialProfile>(m, "AxialProfile")
      .def_readonly("variable", &adsim::AxialProfile::variable)
      .def_readonly("component", &adsim::AxialProfile::component)
      .def_readonly("unit", &adsim::AxialProfile::unit)
      .def_readonly("cycle", &adsim::AxialProfile::cycle)
      .def_readonly("step", &adsim::AxialProfile::step)
      .def_property_readonly("z", [](const adsim::AxialProfile& p) { return ToArray(p.z); })
      .def_property_readonly("time", [](const adsim::AxialProfile& p) { return ToArray(p.time); })
      .def_property_readonly("values",
                             [](const adsim::AxialProfile& p) { return ToArray2D(p.values, p.rows, p.cols); })
      .def("__repr__", [](const adsim::AxialProfile& p) {
        return "<AxialProfile " + p.variable + (p.component.empty() ? "" : "[" + p.component + "]") + " [" +
               p.unit + "], cycle " + std::to_string(p.cycle) + " step " + std::to_string(p.step) + ", " +
               std::to_string(p.rows) + " x " + std::to_string(p.cols) + ">";
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
                     const std::string& log_level, py::object on_message) {
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
          r = s.Run(opt);
        }
        if (pending) throw std::move(*pending);
        return r;
      }, py::arg("cycles") = py::none(), py::arg("progress") = py::none(), py::arg("log_level") = "low",
         py::arg("on_message") = py::none())
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
      .def("css_history", &adsim::Simulation::CSSHistory)
      .def("series", [](const adsim::Simulation& s, const std::string& variable, const std::string& component,
                        const std::string& location, double z, int cycle) {
        return s.Series(variable, component, ParseLocation(location), z, cycle);
      }, py::arg("variable"), py::arg("component") = "", py::arg("location") = "top", py::arg("z") = 0.0,
         py::arg("cycle") = -1)
      .def("profile", &adsim::Simulation::Profile, py::arg("variable"), py::arg("component") = "",
           py::arg("cycle") = -1, py::arg("step") = 0)
      .def("stream_series", &adsim::Simulation::StreamSeries, py::arg("step"), py::arg("stream"),
           py::arg("variable"), py::arg("component") = "", py::arg("cycle") = -1)
      .def("endpoint_series", &adsim::Simulation::EndpointSeries, py::arg("endpoint_id"), py::arg("variable"),
           py::arg("component") = "", py::arg("cycle") = -1)
      .def("save_logs", &adsim::Simulation::SaveLogs, py::arg("directory"))
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
        else if (op == "set_temperature_profile") b::SetTemperatureProfile(d, id(0), str(1));
        else if (op == "set_pump") b::SetPump(d, id(0), str(1), str(2), str(3));
        else if (op == "connect") b::Connect(d, id(0), id(1));
        else if (op == "set_energy_source") b::SetEnergySource(d, id(0), id(1), args[2].cast<bool>());
        else if (op == "set_kpi")
          b::SetKPI(d, str(0), args[1].cast<std::vector<size_t>>(), args[2].cast<std::vector<size_t>>());
        else if (op == "add_variable") b::AddVariable(d, str(0), str(1), str(2), str(3));
        else throw py::value_error("unknown build operation '" + op + "'");
        s.RefreshParameters();
        return out;
      }, py::arg("op"), py::arg("args"))
      // The deck's structure in builder vocabulary: which model kind each
      // object uses, how the steps' ends are set up, the connections and the
      // KPI configuration. Values are read from the parameter registry.
      .def("describe", [](adsim::Simulation& s) {
        namespace b = adsim::build;
        using namespace adsim;
        auto& d = s.Data();
        d.cycle.EnsureIds();
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
          py::dict rates;
          for (const auto& c : d.fluid.components_) {
            auto it = src.mass_transfer_models_.find(c);
            if (it == src.mass_transfer_models_.end() || !it->second || !it->second->mtc_model_) continue;
            py::dict r;
            r["model"] = kind(b::RateModelKinds(), static_cast<int>(it->second->GetModelType()));
            r["coefficient"] = kind(b::RateCoefficientKinds(), static_cast<int>(it->second->mtc_model_->GetModelType()));
            rates[py::str(c)] = r;
          }
          a["rate_models"] = rates;
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
          auto& le = layer.equations_;
          l["dispersion"] = kind(b::DispersionKinds(), static_cast<int>(le.dispersion_model_->GetModelType()));
          l["heat_dispersion"] = kind(b::HeatDispersionKinds(), static_cast<int>(le.heat_dispersion_model_->GetModelType()));
          l["internal_heat_transfer"] =
              kind(b::InternalHeatTransferKinds(), static_cast<int>(le.internal_heat_transfer_model_->GetModelType()));
          l["external_heat_transfer"] =
              kind(b::ExternalHeatTransferKinds(), static_cast<int>(le.external_heat_transfer_model_->GetModelType()));
          py::list comp;
          for (const auto& [id, f] : layer.composition_) comp.append(py::make_tuple(id, f));
          l["composition"] = comp;
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
              if (auto* bs = dynamic_cast<BoundaryStream*>(st.get())) {
                sd["type"] = "stream";
                sd["flow"] = kind(b::FlowKinds(), static_cast<int>(bs->boundary_flow_->GetType()));
                IBoundaryValue* prof = nullptr;
                if (auto* fv = dynamic_cast<VelocityBoundaryValue*>(bs->boundary_flow_.get())) prof = fv->profile_.get();
                else if (auto* fv = dynamic_cast<VolumetricFlowBoundaryValue*>(bs->boundary_flow_.get())) prof = fv->profile_.get();
                else if (auto* fv = dynamic_cast<MolarFlowBoundaryValue*>(bs->boundary_flow_.get())) prof = fv->profile_.get();
                else if (auto* fv = dynamic_cast<MassFlowBoundaryValue*>(bs->boundary_flow_.get())) prof = fv->profile_.get();
                sd["profile"] = tprofile(prof);
                sd["pump"] = pump_of(bs->pump_.get());
              } else if (auto* vs = dynamic_cast<ValveBoundaryStream*>(st.get())) {
                sd["type"] = "valve";
                sd["profile"] = tprofile(vs->boundary_value_.get());
                sd["pump"] = py::none();
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
          kpis.append(k);
        }
        out["kpis"] = kpis;

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

  m.def("build_kinds", []() {
    namespace b = adsim::build;
    py::dict k;
    k["isotherm"] = b::IsothermKinds();
    k["variable_category"] = b::VariableCategoryKinds();
    k["rate_model"] = b::RateModelKinds();
    k["rate_coefficient"] = b::RateCoefficientKinds();
    k["momentum"] = b::MomentumKinds();
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

  m.def("bulk_variables", &adsim::BulkVariableNames);
  m.def("component_variables", &adsim::ComponentVariableNames);
  m.def("stream_variables", &adsim::StreamVariableNames);
  m.def("stream_component_variables", &adsim::StreamComponentVariableNames);

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
  m.def("evaluate_isotherm", [](const std::string& kind, const std::map<std::string, double>& params,
                                const std::vector<double>& pressure, const std::vector<double>& temperature,
                                const std::vector<double>& mole_fraction) {
    const auto pts = MakePoints(pressure, {}, temperature, mole_fraction);
    std::vector<std::pair<std::string, double>> values(params.begin(), params.end());
    std::vector<double> q;
    std::string err;
    if (!adsim::isofit::Evaluate(ParseIsothermKind(kind), values, pts, q, err)) throw py::value_error(err);
    return ToArray(q);
  }, py::arg("kind"), py::arg("params"), py::arg("pressure"), py::arg("temperature"),
     py::arg("mole_fraction") = std::vector<double>{});
  m.def("default_fit_specs", [](const std::string& kind, const std::vector<double>& pressure,
                                const std::vector<double>& loading, const std::vector<double>& temperature) {
    const auto pts = MakePoints(pressure, loading, temperature, {});
    std::vector<adsim::isofit::ParamSpec> specs;
    double t_ref = 0.0;
    std::string err;
    if (!adsim::isofit::DefaultSpecs(ParseIsothermKind(kind), pts, specs, t_ref, err))
      throw py::value_error(err);
    py::list out;
    for (const auto& s : specs) out.append(SpecToDict(s));
    return py::make_tuple(out, t_ref);
  }, py::arg("kind"), py::arg("pressure"), py::arg("loading"), py::arg("temperature"));
  m.def("fit_isotherm", [](const std::string& kind, const std::vector<double>& pressure,
                           const std::vector<double>& loading, const std::vector<double>& temperature,
                           const std::vector<double>& mole_fraction, py::object specs_in, double t_ref,
                           int max_global_evals) {
    const auto type = ParseIsothermKind(kind);
    const auto pts = MakePoints(pressure, loading, temperature, mole_fraction);
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
     py::arg("max_global_evals") = 0);
}
