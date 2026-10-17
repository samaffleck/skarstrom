/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "adsim/simulation.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <cereal/archives/json.hpp>

#include "adsim/solver/kpis.h"
#include "adsim/solver/parameterisation.h"

namespace adsim {

  const char* RunStatusName(RunStatus status) {
    switch (status) {
      case RunStatus::NotRun:     return "not_run";
      case RunStatus::Running:    return "running";
      case RunStatus::Completed:  return "completed";
      case RunStatus::ReachedCSS: return "reached_css";
      case RunStatus::Failed:     return "failed";
      case RunStatus::Cancelled:  return "cancelled";
    }
    return "not_run";
  }

  // -------------------------------------------------------------------------
  // Variable names
  // -------------------------------------------------------------------------

  namespace {

    struct BulkName { const char* name; VariableType type; };
    struct ComponentName { const char* name; ComponentVariableType type; };
    struct StreamName { const char* name; BoundaryVariableType type; };
    struct StreamComponentName { const char* name; BoundaryComponentVariableType type; };

    // snake_case of VariableTypeList, in enum order.
    constexpr BulkName kBulk[] = {
      {"pressure", VariableType::PRESSURE},
      {"temperature", VariableType::TEMPERATURE},
      {"wall_temperature", VariableType::WALL_TEMPERATURE},
      {"velocity", VariableType::VELOCITY},
      {"volumetric_flow_rate", VariableType::VOLUMETRIC_FLOW_RATE},
      {"molar_flow_rate", VariableType::MOLAR_FLOW_RATE},
      {"total_mole_fraction", VariableType::TOTAL_MOLE_FRACTION},
      {"approach_to_fluidisation", VariableType::APPROACH_TO_FLUIDISATION},
      {"viscosity", VariableType::VISCOSITY},
      {"fluid_density", VariableType::FLUID_DENSITY},
      {"mass_source", VariableType::MASS_SOURCE},
      {"energy_source", VariableType::ENERGY_SOURCE},
      {"molar_concentration", VariableType::MOLAR_CONCENTRATION},
      {"heat_capacity", VariableType::HEAT_CAPACITY},
      {"molar_mass", VariableType::MOLAR_MASS},
      {"thermal_conductivity", VariableType::THERMAL_CONDUCTIVITY},
      {"heat_dispersion", VariableType::HEAT_DISPERSION},
      {"schmidt_number", VariableType::SCHMIDT_NUMBER},
      {"reynolds_number", VariableType::REYNOLDS_NUMBER},
      {"sherwood_number", VariableType::SHERWOOD_NUMBER},
      {"pressure_gradient", VariableType::PRESSURE_GRADIENT},
      {"internal_heat_transfer_coefficient", VariableType::INTERNAL_HEAT_TRANSFER_COEFFICIENT},
      {"external_heat_transfer_coefficient", VariableType::EXTERNAL_HEAT_TRANSFER_COEFFICIENT},
      {"axial_dispersion_coefficient", VariableType::AXIAL_DISPERSION_COEFFICIENT},
      {"solid_temperature", VariableType::SOLID_TEMPERATURE},
      {"particle_heat_transfer_coefficient", VariableType::PARTICLE_HEAT_TRANSFER_COEFFICIENT},
    };
    static_assert(std::size(kBulk) == std::size(VariableTypeList),
                  "kBulk must name every VariableType (keep it in VariableTypeList order)");

    constexpr ComponentName kComponent[] = {
      {"mole_fraction", ComponentVariableType::MOLE_FRACTION},
      {"adsorbed_concentration", ComponentVariableType::ADSORBED_CONCENTRATION},
      {"equilibrium_adsorbed_concentration", ComponentVariableType::EQUILIBRIUM_ADSORBED_CONCENTRATION},
      {"component_mass_source", ComponentVariableType::COMPONENT_MASS_SOURCE},
      {"component_energy_source", ComponentVariableType::COMPONENT_ENERGY_SOURCE},
      {"component_molar_concentration", ComponentVariableType::COMPONENT_MOLAR_CONCENTRATION},
      {"mass_transfer_coefficient", ComponentVariableType::MASS_TRANSFER_COEFFICIENT},
      {"molecular_diffusion", ComponentVariableType::MOLECULAR_DIFFUSION},
    };
    static_assert(std::size(kComponent) == std::size(ComponentVariableTypeList),
                  "kComponent must name every ComponentVariableType");

    constexpr StreamName kStream[] = {
      {"temperature", BoundaryVariableType::TEMPERATURE},
      {"molar_flow_rate", BoundaryVariableType::MOLAR_FLOW_RATE},
      {"compressor_power", BoundaryVariableType::COMPRESSOR_POWER},
    };
    constexpr StreamComponentName kStreamComponent[] = {
      {"mole_fraction", BoundaryComponentVariableType::MOLE_FRACTION},
      {"molar_flow", BoundaryComponentVariableType::MOLAR_FLOW},
    };

    template <class Table>
    std::vector<std::string> Names(const Table& table) {
      std::vector<std::string> out;
      for (const auto& e : table) out.emplace_back(e.name);
      return out;
    }

    template <class Table>
    auto Lookup(const Table& table, const std::string& name) -> std::optional<decltype(table[0].type)> {
      for (const auto& e : table)
        if (name == e.name) return e.type;
      return std::nullopt;
    }

    std::string Join(const std::vector<std::string>& v) {
      std::string s;
      for (size_t i = 0; i < v.size(); ++i) {
        if (i) s += ", ";
        s += v[i];
      }
      return s;
    }

    // "[Pa]" -> "Pa"
    std::string StripBrackets(const std::string& unit) {
      if (unit.size() >= 2 && unit.front() == '[' && unit.back() == ']') return unit.substr(1, unit.size() - 2);
      return unit;
    }

  }  // namespace

  std::vector<std::string> BulkVariableNames() { return Names(kBulk); }
  std::vector<std::string> ComponentVariableNames() { return Names(kComponent); }
  std::vector<std::string> StreamVariableNames() { return Names(kStream); }
  std::vector<std::string> StreamComponentVariableNames() { return Names(kStreamComponent); }

  // -------------------------------------------------------------------------
  // Loading
  // -------------------------------------------------------------------------

  void LoadFromJsonText(UserData& data, const std::string& json) {
    if (json.empty()) throw std::invalid_argument("deck JSON is empty");
    LoadUserDataTransactional(data, [&](UserData& loaded) {
      try {
        // String format root ("simulation"). The NVP lookup fails at the first
        // read, before `loaded` is touched, when the root is the file format.
        std::istringstream iss(json);
        cereal::JSONInputArchive ar(iss);
        ar(cereal::make_nvp("simulation", loaded));
      } catch (const cereal::Exception&) {
        loaded = UserData{};
        std::istringstream iss(json);
        cereal::JSONInputArchive ar(iss);
        ar(loaded);
      }
    });
  }

  // -------------------------------------------------------------------------
  // Simulation: deck
  // -------------------------------------------------------------------------

  Simulation::Simulation() : data_(std::make_unique<UserData>()) {}
  Simulation::~Simulation() = default;
  Simulation::Simulation(Simulation&&) noexcept = default;
  Simulation& Simulation::operator=(Simulation&&) noexcept = default;

  Simulation Simulation::FromFile(const std::string& path) {
    Simulation s;
    s.Load(path);
    return s;
  }

  Simulation Simulation::FromJson(const std::string& json) {
    Simulation s;
    s.LoadJson(json);
    return s;
  }

  void Simulation::ResetResults() {
    if (chunked_ && chunked_->IsActive()) chunked_->Stop();
    solver_.reset();
    chunked_.reset();
    registry_.reset();
    result_ = RunResult{};
    messages_.clear();
  }

  void Simulation::Load(const std::string& path) {
    auto fresh = std::make_unique<UserData>();
    LoadFromJson(*fresh, path);
    fresh->EvaluateAllExpressions();
    ResetResults();
    fresh->load_generation_ = data_->load_generation_ + 1;
    data_ = std::move(fresh);
  }

  void Simulation::LoadJson(const std::string& json) {
    auto fresh = std::make_unique<UserData>();
    LoadFromJsonText(*fresh, json);
    fresh->EvaluateAllExpressions();
    ResetResults();
    fresh->load_generation_ = data_->load_generation_ + 1;
    data_ = std::move(fresh);
  }

  void Simulation::SetData(UserData&& data) {
    if (IsRunning()) throw std::logic_error("cannot replace the deck while a run is in progress");
    // Move-assigned in place: the deck keeps its address, so a finished
    // solver (which holds a reference to it) and its results stay valid.
    *data_ = std::move(data);
    registry_.reset();
  }

  void Simulation::Save(const std::string& path) {
    std::ofstream file(path);
    if (!file.is_open()) throw std::runtime_error("cannot write " + path);
    data_->info.sim_version++;
    cereal::JSONOutputArchive archive(file);
    archive(*data_);
  }

  std::string Simulation::ToJson() {
    std::string out;
    SaveToString(*data_, out);
    return out;
  }

  bool Simulation::Check(std::string& error) {
    error.clear();
    if (!data_->EvaluateAllExpressions(error)) return false;
    return data_->PerformChecks(error);
  }

  // -------------------------------------------------------------------------
  // Simulation: parameters
  // -------------------------------------------------------------------------

  ParameterRegistry& Simulation::Parameters() {
    if (!registry_) registry_ = std::make_unique<ParameterRegistry>(*data_);
    return *registry_;
  }

  void Simulation::RefreshParameters() {
    if (registry_) registry_->Refresh();
    else registry_ = std::make_unique<ParameterRegistry>(*data_);
  }

  double Simulation::GetParameter(std::string_view path) {
    return Parameters().Get(path).GetValue();
  }

  std::string Simulation::GetParameterText(std::string_view path) {
    return Parameters().Get(path).GetText();
  }

  void Simulation::SetParameter(std::string_view path, double value) {
    if (IsRunning()) throw std::logic_error("cannot change parameters while a run is in progress");
    Parameter& p = Parameters().Get(path);
    p.SetValue(value);
    data_->EvaluateAllExpressions();
    // Enumeration depends on some values (e.g. a component's mole fraction
    // never does, but the step-time mode or a pump switch could); refresh so
    // later lookups see the current structure.
    registry_->Refresh();
  }

  void Simulation::SetParameterText(std::string_view path, const std::string& text) {
    if (IsRunning()) throw std::logic_error("cannot change parameters while a run is in progress");
    Parameter& p = Parameters().Get(path);
    const ParameterKind kind = p.Info().kind;
    const std::string resolved = p.Path();
    const std::string previous = p.GetText();
    p.SetText(text);
    std::string err;
    if (!data_->EvaluateAllExpressions(err)) {
      // Restore the old text so the deck stays usable.
      p.SetText(previous);
      data_->EvaluateAllExpressions();
      throw std::invalid_argument("cannot set " + resolved + " = '" + text + "': " + err);
    }
    if (kind == ParameterKind::Expression) {
      const ExpressionParam* e = p.ExpressionPointer();
      if (e && !e->LastError().empty()) {
        const std::string why = e->LastError();
        p.SetText(previous);
        data_->EvaluateAllExpressions();
        throw std::invalid_argument("cannot set " + resolved + " = '" + text + "': " + why);
      }
    }
    registry_->Refresh();
  }

  // -------------------------------------------------------------------------
  // Simulation: run
  // -------------------------------------------------------------------------

  RunResult Simulation::Run(const RunOptions& options) {
    if (IsRunning()) throw std::logic_error("a chunked run is already in progress");
    ResetResults();
    UserData& data = *data_;
    if (options.cycles) data.cycle.number_of_cycles_ = std::max(1, *options.cycles);
    {
      std::string err;
      if (!data.EvaluateAllExpressions(err)) {
        result_.status = RunStatus::Failed;
        result_.return_code = -1;
        result_.message = "Expression evaluation failed: " + err;
        return result_;
      }
    }

    solver_ = std::make_unique<CoupledSolver>(data);
    CoupledSolver& solver = *solver_;
    solver.SetLogLevel(options.log_level);
    const LogLevel keep = options.log_level;
    auto forward = options.on_message;
    solver.SetMessageCallback([this, keep, forward](LogLevel level, const std::string& msg) {
      if (static_cast<int>(level) > static_cast<int>(keep)) return;
      messages_.push_back({level, msg});
      if (forward) forward(level, msg);
    });

    bool cancelled = false;
    if (options.progress) {
      auto progress = options.progress;
      solver.SetProgressCallback([&solver, &cancelled, progress](int c, int s, int tc, int ts) {
        RunProgress p{c, s, tc, ts, 0.0};
        const double total = static_cast<double>(tc) * static_cast<double>(ts);
        p.fraction = total > 0.0 ? (static_cast<double>(c) * ts + s + 1) / total : 0.0;
        if (!progress(p)) {
          cancelled = true;
          solver.RequestStop();
        }
      });
    }
    if (options.cancel || options.on_step) {
      auto cancel = options.cancel;
      auto on_step = options.on_step;
      solver.SetMetricsCallback([&solver, &cancelled, cancel, on_step](double t, double dt, double err, int c, int st) {
        if (on_step) on_step(t, dt, err, c, st);
        if (cancel && !cancelled && cancel()) {
          cancelled = true;
          solver.RequestStop();
        }
      });
    }

    const auto t0 = std::chrono::steady_clock::now();
    int rc = -1;
    std::string exception_text;
    try {
      rc = solver.IntegrateCycle(data.cycle, false);
    } catch (const std::exception& e) {
      exception_text = e.what();
    }
    result_.wall_time_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    // The callbacks capture locals of this call; the solver outlives it.
    solver.SetProgressCallback({});
    solver.SetMetricsCallback({});
    result_.return_code = rc;
    result_.stats = solver.GetStatistics();
    result_.css_cycle = solver.GetCSSCycleIndex();
    result_.css_reached = solver.ReachedCSS();
    const auto& cycles = solver.GetLog().sys_data.cycle_index_set_;
    result_.cycles_run = cycles.empty() ? 0 : *cycles.rbegin() + 1;

    if (!exception_text.empty()) {
      result_.status = RunStatus::Failed;
      result_.message = "solver threw: " + exception_text;
    } else if (rc == 99 || cancelled || solver.IsStopRequested()) {
      result_.status = RunStatus::Cancelled;
      result_.message = "run cancelled";
    } else if (rc < 0) {
      result_.status = RunStatus::Failed;
      // The last LOW message holds the reason (checks, initialisation or the
      // failing step).
      result_.message = "run failed";
      for (auto it = messages_.rbegin(); it != messages_.rend(); ++it) {
        if (it->text.rfind("Solver failed after", 0) == 0) continue;
        result_.message = it->text;
        break;
      }
    } else if (result_.css_reached) {
      result_.status = RunStatus::ReachedCSS;
      result_.message = "cyclic steady state reached at cycle " + std::to_string(result_.css_cycle + 1);
    } else {
      result_.status = RunStatus::Completed;
      result_.message = "completed " + std::to_string(result_.cycles_run) + " cycle(s)";
    }
    return result_;
  }

  void Simulation::RequestStop() {
    if (solver_) solver_->RequestStop();
    if (chunked_) chunked_->Stop();
  }

  bool Simulation::Start(LogLevel log_level, SolverMessageCallback on_message) {
    if (IsRunning()) throw std::logic_error("a chunked run is already in progress");
    ResetResults();
    chunked_log_level_ = log_level;
    chunked_on_message_ = std::move(on_message);
    chunked_ = std::make_unique<CoupledSolverChuncked>(*data_);
    chunked_->SetLogLevel(log_level);
    chunked_->SetMessageCallback([this](LogLevel level, const std::string& msg) {
      messages_.push_back({level, msg});
      if (chunked_on_message_) chunked_on_message_(level, msg);
    });
    chunked_->Start();
    const auto st = chunked_->GetStatus();
    if (st.GetStatus() == SolverStatus::FAILED) {
      result_.status = RunStatus::Failed;
      result_.message = st.GetMessage();
      return false;
    }
    result_.status = RunStatus::Running;
    return true;
  }

  RunStatus Simulation::Advance(int max_time_ms) {
    if (!chunked_) return result_.status;
    if (chunked_->IsActive()) {
      const auto t0 = std::chrono::steady_clock::now();
      chunked_->Integrate(max_time_ms);
      result_.wall_time_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }
    const auto st = chunked_->GetStatus();
    const auto& cycles = chunked_->GetLog().sys_data.cycle_index_set_;
    result_.cycles_run = cycles.empty() ? 0 : *cycles.rbegin() + 1;
    switch (st.GetStatus()) {
      case SolverStatus::IN_PROGRESS:
        result_.status = RunStatus::Running;
        break;
      case SolverStatus::COMPLETED:
        result_.status = RunStatus::Completed;
        result_.message = st.GetMessage();
        break;
      case SolverStatus::AT_CSS:
        result_.status = RunStatus::ReachedCSS;
        result_.message = st.GetMessage();
        result_.css_reached = true;
        result_.css_cycle = cycles.empty() ? -1 : *cycles.rbegin();
        break;
      case SolverStatus::FAILED:
        result_.status = RunStatus::Failed;
        result_.message = st.GetMessage();
        break;
      case SolverStatus::NOT_STARTED:
        result_.status = RunStatus::Cancelled;
        result_.message = st.GetMessage();
        break;
    }
    return result_.status;
  }

  void Simulation::Stop() {
    if (chunked_ && chunked_->IsActive()) {
      chunked_->Stop();
      result_.status = RunStatus::Cancelled;
      result_.message = chunked_->GetStatus().GetMessage();
    }
    if (solver_) solver_->RequestStop();
  }

  bool Simulation::IsRunning() const { return chunked_ && chunked_->IsActive(); }

  double Simulation::ProgressPercent() const {
    if (chunked_) return chunked_->GetProgressPercent();
    return result_.Ok() ? 100.0 : 0.0;
  }

  // -------------------------------------------------------------------------
  // Simulation: results
  // -------------------------------------------------------------------------

  namespace {
    const EquationLog& EmptyEquationLog() { static const EquationLog log{}; return log; }
    const BoundaryLog& EmptyBoundaryLog() { static const BoundaryLog log{}; return log; }
  }

  bool Simulation::HasResults() const { return !Log().GetSimulationLog().empty(); }

  const EquationLog& Simulation::Log() const {
    if (solver_) return solver_->GetLog();
    if (chunked_) return chunked_->GetLog();
    return EmptyEquationLog();
  }

  const BoundaryLog& Simulation::Boundaries() const {
    if (solver_) return solver_->GetBoundaryLog();
    if (chunked_) return chunked_->GetBoundaryLog();
    return EmptyBoundaryLog();
  }

  std::vector<int> Simulation::Cycles() const {
    std::vector<int> out;
    for (const auto& [c, steps] : Log().GetSimulationLog()) {
      (void)steps;
      out.push_back(c);
    }
    return out;
  }

  int Simulation::LastCycle() const {
    const auto& log = Log().GetSimulationLog();
    return log.empty() ? -1 : log.rbegin()->first;
  }

  int Simulation::ResolveCycle(int cycle) const {
    const int c = cycle < 0 ? LastCycle() : cycle;
    if (c < 0) throw std::logic_error("no results: run the simulation first");
    if (Log().GetSimulationLog().find(c) == Log().GetSimulationLog().end())
      throw std::out_of_range("cycle " + std::to_string(c) + " is not in the log (logged cycles are the last " +
                              std::to_string(data_->settings.max_cycles_to_keep_) + ")");
    return c;
  }

  std::vector<std::string> Simulation::Components() const { return data_->fluid.components_; }

  std::vector<std::string> Simulation::StepNames() const {
    std::vector<std::string> out;
    for (const auto& s : data_->cycle.step_list_) out.push_back(s.step_name_);
    return out;
  }

  std::vector<StreamInfo> Simulation::Streams() const {
    std::vector<StreamInfo> out;
    // Endpoints() only reads; it is non-const because its handles are mutable.
    for (const auto& h : const_cast<Cycle&>(data_->cycle).Endpoints()) {
      StreamInfo s;
      s.step = static_cast<int>(h.step_index);
      s.step_name = h.step->step_name_;
      s.name = h.log_name;
      s.endpoint_id = h.id;
      s.location = h.location;
      out.push_back(std::move(s));
    }
    return out;
  }

  std::vector<double> Simulation::Positions() const { return Log().GetDomain().positions; }

  std::vector<KPIRecord> Simulation::KPIs(int cycle) const {
    std::vector<KPIRecord> out;
    if (!HasResults()) return out;
    const int c = ResolveCycle(cycle);
    const auto& eq_log = Log();
    const auto& b_log = Boundaries();
    // Re-bind the endpoint refs by ID first, so a step or stream renamed
    // since they were last resolved still counts (on a copy: this is const).
    KPIConfig config = data_->kpi_config;
    config.ResolveIds(data_->cycle.step_list_);
    for (const auto& cfg : config.per_component) {
      if (!cfg.tracked) continue;
      KPIParams p(eq_log, b_log);
      p.cycle_index_ = c;
      p.product_comp_ = cfg.product_component;
      p.feed_endpoints = cfg.feed_endpoints;
      p.product_endpoints = cfg.product_endpoints;
      KPIRecord r;
      r.cycle = c;
      r.component = cfg.product_component;
      r.valid = p.Check();
      if (r.valid) {
        constexpr double nan = std::numeric_limits<double>::quiet_NaN();
        bool ok = true;
        auto value = [&](double v) { return (ok && std::isfinite(v)) ? v : nan; };
        double v = CalculatePurity(p, &ok);                     r.purity = value(v);
        v = CalculateRecovery(p, &ok);                          r.recovery = value(v);
        v = CalculateProductivityMass(p, &ok);                  r.productivity_mass = value(v);
        v = CalculateProductivityVolume(p, &ok);                r.productivity_volume = value(v);
        v = CalculateSpecificEnergyConsumption(p, &ok);         r.energy_consumption = value(v);
        v = CalculateProductionRate(p, &ok);                    r.production_rate = value(v);
      }
      out.push_back(std::move(r));
    }
    return out;
  }

  std::vector<KPIRecord> Simulation::KPIs() const {
    std::vector<KPIRecord> out;
    for (int c : Cycles()) {
      auto rows = KPIs(c);
      out.insert(out.end(), rows.begin(), rows.end());
    }
    return out;
  }

  std::vector<CSSRecord> Simulation::CSSHistory() {
    std::vector<CSSRecord> out;
    for (int c : Cycles()) {
      CSSResult r = data_->settings.css.CheckCSS(Log(), c);
      out.push_back({c, r.converged, r.error, r.tolerance, r.method_name});
    }
    return out;
  }

  TimeSeries Simulation::Series(const std::string& variable, const std::string& component,
                                SeriesLocation where, double z, int cycle) const {
    const int c = ResolveCycle(cycle);
    const ParamSpatialLoc loc = where == SeriesLocation::Top      ? ParamSpatialLoc::Top
                              : where == SeriesLocation::Bottom   ? ParamSpatialLoc::Bottom
                                                                  : ParamSpatialLoc::Position;
    TimeSeries ts;
    ts.variable = variable;
    ts.component = component;
    if (auto bulk = Lookup(kBulk, variable)) {
      if (!component.empty())
        throw std::invalid_argument("'" + variable + "' is a bulk variable; do not pass a component");
      ExtractBulkTimeSeries(Log(), *bulk, loc, z, ts.time, ts.values, c);
      ts.unit = StripBrackets(GetVariableUnits(*bulk));
    } else if (auto comp = Lookup(kComponent, variable)) {
      if (component.empty())
        throw std::invalid_argument("'" + variable + "' is a component variable; pass a component");
      ExtractComponentTimeSeries(Log(), *comp, component, loc, z, ts.time, ts.values, c);
      ts.unit = StripBrackets(GetComponentVariableUnits(*comp));
    } else {
      throw std::invalid_argument("unknown variable '" + variable + "'. Bulk: " + Join(BulkVariableNames()) +
                                  ". Component: " + Join(ComponentVariableNames()));
    }
    return ts;
  }

  AxialProfile Simulation::Profile(const std::string& variable, const std::string& component,
                                   int cycle, int step) const {
    const int c = ResolveCycle(cycle);
    VariableLogEntry e;
    AxialProfile p;
    if (auto bulk = Lookup(kBulk, variable)) {
      if (!component.empty())
        throw std::invalid_argument("'" + variable + "' is a bulk variable; do not pass a component");
      e = Log().GetVariable(c, step, *bulk);
      p.unit = StripBrackets(GetVariableUnits(*bulk));
    } else if (auto comp = Lookup(kComponent, variable)) {
      if (component.empty())
        throw std::invalid_argument("'" + variable + "' is a component variable; pass a component");
      e = Log().GetComponentVariable(c, step, *comp, component);
      p.unit = StripBrackets(GetComponentVariableUnits(*comp));
    } else {
      throw std::invalid_argument("unknown variable '" + variable + "'. Bulk: " + Join(BulkVariableNames()) +
                                  ". Component: " + Join(ComponentVariableNames()));
    }
    if (e.times.empty())
      throw std::out_of_range("no log data for '" + variable + "' at cycle " + std::to_string(c) + ", step " +
                              std::to_string(step));
    p.variable = variable;
    p.component = component;
    p.cycle = c;
    p.step = step;
    p.z = Log().GetDomain().positions;
    p.time = e.times;
    p.rows = e.values.size();
    p.cols = p.z.size();
    p.values.assign(p.rows * p.cols, std::numeric_limits<double>::quiet_NaN());
    for (size_t t = 0; t < p.rows; ++t) {
      const auto& row = e.values[t];
      for (size_t i = 0; i < p.cols && i < row.size(); ++i) p.values[t * p.cols + i] = row[i];
    }
    return p;
  }

  TimeSeries Simulation::StreamSeries(int step, const std::string& stream, const std::string& variable,
                                      const std::string& component, int cycle) const {
    const int c = ResolveCycle(cycle);
    BoundaryVariableData d;
    if (component.empty()) {
      auto v = Lookup(kStream, variable);
      if (!v)
        throw std::invalid_argument("unknown stream variable '" + variable + "'. Bulk: " +
                                    Join(StreamVariableNames()) + ". Per component (pass a component): " +
                                    Join(StreamComponentVariableNames()));
      d = Boundaries().GetVariableData(c, step, stream, *v);
    } else {
      auto v = Lookup(kStreamComponent, variable);
      if (!v)
        throw std::invalid_argument("unknown per-component stream variable '" + variable + "'. Expected: " +
                                    Join(StreamComponentVariableNames()));
      d = Boundaries().GetComponentVariableData(c, step, stream, *v, component);
    }
    if (d.time.empty()) {
      std::string known;
      for (const auto& s : Streams())
        if (s.step == step) known += (known.empty() ? "" : ", ") + ("'" + s.name + "'");
      throw std::out_of_range("no log data for stream '" + stream + "' on step " + std::to_string(step) +
                              " at cycle " + std::to_string(c) + " (endpoints on this step: " + known + ")");
    }
    TimeSeries ts;
    ts.variable = variable;
    ts.component = component;
    ts.unit = StripBrackets(d.units);
    ts.time = std::move(d.time);
    ts.values = std::move(d.variable_timeseries);
    return ts;
  }

  TimeSeries Simulation::EndpointSeries(size_t endpoint_id, const std::string& variable,
                                        const std::string& component, int cycle) const {
    for (const auto& s : Streams())
      if (s.endpoint_id == endpoint_id) return StreamSeries(s.step, s.name, variable, component, cycle);
    throw std::out_of_range("no endpoint with ID " + std::to_string(endpoint_id) + " in this deck");
  }

  void Simulation::SaveLogs(const std::string& directory) const { Log().SaveToFile(directory); }

}  // namespace adsim
