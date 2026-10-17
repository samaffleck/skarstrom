/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/solver/log.h"
#include "adsim/solver/cycle.h"
#include "adsim/objects/reactor.h"
#include "adsim/objects/user_data.h"

#include <filesystem>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <charconv>

namespace adsim {

  void EquationLog::Initialise(UserData& user_data) {
    auto& settings = user_data.settings;
    auto& fluid = user_data.fluid;
    auto& bed = user_data.bed;
    auto& cycle = user_data.cycle;
    auto& state = user_data.state;

    Clear();

    domain = state.domain_vectors;
    max_cycles_to_keep_ = std::max(1, settings.max_cycles_to_keep_);
    max_log_points_per_step_ = settings.max_log_points_per_step_;

    sys_data.components_ = fluid.components_;
    sys_data.adsorbent_mass = GetAdsorbentMass(bed, user_data.adsorbent_directory);
    sys_data.molar_mass = fluid.molar_mass_;
    sys_data.step_list_.clear();
    for (const auto& step : cycle.step_list_) {
      sys_data.step_list_.emplace_back(step.step_name_);
    }
    sys_data.total_cycle_time = cycle.GetCycleTime();
  }

  void EquationLog::LogState(const Equations& state, double t, int step_index, int cycle_index) {
    // Check if we're starting a new cycle
    bool is_new_cycle = (current_logged_cycle_ != cycle_index);
    bool is_new_step = is_new_cycle || (current_logged_step_ != step_index);
    if (is_new_cycle) {
      current_logged_cycle_ = cycle_index;
      // Evict old cycles if we're over the limit (sliding window)
      // max_cycles_to_keep_ is always >= 1, so always check
      EvictOldCyclesIfNeeded();
    }
    current_logged_step_ = step_index;

    // The domain snapshot only changes between steps (PopulateEnergySources at
    // step initialisation), so refresh it once per (cycle, step) instead of
    // deep-copying every solver step.
    if (is_new_step) {
      domain = state.domain_vectors;
      sys_data.cycle_index_set_.insert(cycle_index);
      sys_data.step_index_set_.insert(step_index);
    }

    // Log regular variables (coupled equations)
    // Use [] operator to create entries if they don't exist
    auto& sim_log = var_log_[cycle_index][step_index];

    // OPTIMIZATION: Store times once per step, not per variable
    sim_log.step_data.times.push_back(t);
    sim_log.step_data.cycle_index = cycle_index;
    sim_log.step_data.step_index = step_index;

    // Entry metadata (name, units, step name) is constant for the lifetime of
    // an entry — assign it only on first use rather than re-building the
    // strings on every solver step.
    const std::string& step_name = sys_data.step_list_[step_index];

    auto log_var = [&](VariableType type, const std::vector<double>& values) {
      auto& entry = sim_log.variables[type];
      if (entry.values.empty()) {
        entry.variable_name = VariableTypeList[static_cast<size_t>(type)];
        entry.units = GetVariableUnits(type);
        entry.step_name = step_name;
      }
      entry.values.push_back(CopyRow(values));
    };

    auto log_comp_var = [&](ComponentVariableType type, const auto& eq_map) {
      for (const auto& [component, eq] : eq_map) {
        auto& entry = sim_log.component_variables[component][type];
        if (entry.values.empty()) {
          entry.variable_name = component + " " + ComponentVariableTypeList[static_cast<size_t>(type)];
          entry.units = GetComponentVariableUnits(type);
          entry.step_name = step_name;
        }
        entry.values.push_back(CopyRow(eq.var));
      }
    };

    log_var(VariableType::PRESSURE, state.pressure.var);
    log_var(VariableType::TEMPERATURE, state.fluid_temperature.var);
    log_var(VariableType::WALL_TEMPERATURE, state.wall_temperature.var);
    log_var(VariableType::VELOCITY, state.velocity.var);
    log_var(VariableType::MOLAR_CONCENTRATION, state.concentration.var);
    log_var(VariableType::TOTAL_MOLE_FRACTION, state.total_mole_fraction.var);
    log_var(VariableType::FLUID_DENSITY, state.fluid_density.var);
    log_var(VariableType::VISCOSITY, state.viscosity.var);
    log_var(VariableType::HEAT_CAPACITY, state.heat_capacity.var);
    log_var(VariableType::REYNOLDS_NUMBER, state.reynolds_number.var);
    log_var(VariableType::SCHMIDT_NUMBER, state.schmidt_number.var);
    log_var(VariableType::SHERWOOD_NUMBER, state.sherwood_number.var);
    log_var(VariableType::THERMAL_CONDUCTIVITY, state.fluid_thermal_conductivity.var);
    log_var(VariableType::HEAT_DISPERSION, state.axial_thermal_conductivity.var);
    log_var(VariableType::INTERNAL_HEAT_TRANSFER_COEFFICIENT, state.internal_heat_transfer_coefficient.var);
    log_var(VariableType::EXTERNAL_HEAT_TRANSFER_COEFFICIENT, state.external_heat_transfer_coefficient.var);
    log_var(VariableType::MOLAR_MASS, state.molar_mass.var);
    log_var(VariableType::VOLUMETRIC_FLOW_RATE, state.volumetric_flow_rate.var);
    log_var(VariableType::MOLAR_FLOW_RATE, state.molar_flow_rate.var);
    log_var(VariableType::MASS_SOURCE, state.mass_source.var);
    log_var(VariableType::ENERGY_SOURCE, state.energy_source.var);
    log_var(VariableType::APPROACH_TO_FLUIDISATION, state.approach_to_fluidisation.var);
    log_var(VariableType::PRESSURE_GRADIENT, state.pressure_gradient.var);
    log_var(VariableType::AXIAL_DISPERSION_COEFFICIENT, state.axial_dispersion_coefficient.var);
    log_var(VariableType::SOLID_TEMPERATURE, state.solid_temperature.var);

    // Particle HTC only present in LTNE mode
    if (!state.particle_heat_transfer_coefficient.var.empty()) {
      log_var(VariableType::PARTICLE_HEAT_TRANSFER_COEFFICIENT, state.particle_heat_transfer_coefficient.var);
    }

    log_comp_var(ComponentVariableType::COMPONENT_MOLAR_CONCENTRATION, state.molar_concentration);
    log_comp_var(ComponentVariableType::MOLE_FRACTION, state.mole_fraction);
    log_comp_var(ComponentVariableType::MOLECULAR_DIFFUSION, state.molecular_diffusion);
    log_comp_var(ComponentVariableType::COMPONENT_MASS_SOURCE, state.component_mass_source);
    log_comp_var(ComponentVariableType::COMPONENT_ENERGY_SOURCE, state.component_energy_source);

    // Adsorbed-phase quantities are tracked per adsorbent internally. Log the
    // mass-fraction-weighted aggregate per component (sum_k w_k * x_{k,comp}),
    // which reduces to the single adsorbent's value for a single-adsorbent bed.
    // TODO(per-adsorbent curves): also emit one curve per (component, adsorbent);
    // requires the AdsorbentDirectory for names and matching GUI/CSV plumbing.
    const auto& cell_ads = state.domain_vectors.cell_adsorbents;
    const auto& bed_cells = state.domain_vectors.bed_cell_indices;
    auto aggregate_per_component = [&](const auto& nested, const std::string& comp) {
      size_t M = nested.empty() ? 0 : nested.begin()->second.at(comp).var.size();
      // Resolve each adsorbent's series for `comp` once, not per cell.
      static thread_local std::vector<std::pair<AdsorbentID, const std::vector<double>*>> series;
      series.clear();
      for (const auto& [ads_id, comp_map] : nested) series.emplace_back(ads_id, &comp_map.at(comp).var);
      Var1D out = TakeRow();
      out.assign(M, 0.0);
      for (size_t i = 0; i < M; ++i) {
        double s = 0.0;
        for (const auto& ca : cell_ads[i]) {
          const std::vector<double>* x = nullptr;
          for (const auto& [ads_id, ptr] : series)
            if (ads_id == ca.id) { x = ptr; break; }
          if (!x) x = &nested.at(ca.id).at(comp).var;  // throws, as before, for an unknown id
          s += ca.mass_fraction * (*x)[i];
        }
        out[i] = s;
      }
      if (M > 0 && !bed_cells.empty()) {  // mirror ghost cells onto the first/last bed cell
        out[0] = out[bed_cells.front()];
        out[M - 1] = out[bed_cells.back()];
      }
      return out;
    };
    auto log_aggregate = [&](ComponentVariableType type, const auto& nested) {
      for (const auto& [component, eq] : state.mole_fraction) {
        auto& entry = sim_log.component_variables[component][type];
        if (entry.values.empty()) {
          entry.variable_name = component + " " + ComponentVariableTypeList[static_cast<size_t>(type)];
          entry.units = GetComponentVariableUnits(type);
          entry.step_name = step_name;
        }
        entry.values.push_back(aggregate_per_component(nested, component));
      }
    };
    log_aggregate(ComponentVariableType::ADSORBED_CONCENTRATION, state.ads_concentration);
    log_aggregate(ComponentVariableType::EQUILIBRIUM_ADSORBED_CONCENTRATION, state.eq_ads_concentration);
    log_aggregate(ComponentVariableType::MASS_TRANSFER_COEFFICIENT, state.mass_transfer_coefficient);
  }

  // Format doubles with std::to_chars (printf "%.6g" semantics — identical
  // text to the default ostream operator<<) instead of going through the
  // locale-aware num_put machinery, which dominates the cost of writing
  // large result sets.
  static inline void AppendDouble(std::string& out, double v) {
    char buf[40];
    auto res = std::to_chars(buf, buf + sizeof(buf), v, std::chars_format::general, 6);
    out.append(buf, res.ptr);
  }

  static void WriteCsvRows(std::ofstream& file, int cycle_index, int step_index,
                           const std::vector<double>& times,
                           const std::vector<std::vector<double>>& values,
                           int max_points) {
    const size_t n = std::min(times.size(), values.size());
    if (n == 0) return;

    // Decimate the written rows to at most ~max_points per (cycle, step): the
    // solver logs every internal step in memory (needed at full resolution for
    // KPIs, mass-balance closure and CSS), but the CSV is only for plotting, so
    // a uniform subsample keeps the files small and the write fast. The first
    // and last row of the step are always emitted. <= 0 writes every row.
    size_t stride = 1;
    if (max_points > 0 && n > static_cast<size_t>(max_points)) {
      stride = (n + static_cast<size_t>(max_points) - 1) / static_cast<size_t>(max_points);
    }

    const std::string prefix = std::to_string(cycle_index) + "," + std::to_string(step_index) + ",";
    std::string line;
    auto write_row = [&](size_t i) {
      line.clear();
      line += prefix;
      AppendDouble(line, times[i]);
      for (double v : values[i]) {
        line += ',';
        AppendDouble(line, v);
      }
      line += '\n';
      file.write(line.data(), static_cast<std::streamsize>(line.size()));
    };

    size_t last_written = n;  // sentinel
    for (size_t i = 0; i < n; i += stride) {
      write_row(i);
      last_written = i;
    }
    if (last_written != n - 1) write_row(n - 1);  // always keep the final point
  }

  static inline std::string ToFilenameSafe(std::string name) {
    std::replace(name.begin(), name.end(), ' ', '_');
    for (char& c : name) {
      if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-' || c == '.')) c = '_';
    }
    return name;
  }

  std::string EquationLog::SanitizeFilename(const std::string& name) const {
    return ToFilenameSafe(name);
  }

  void EquationLog::SaveToFile(const std::string& directory) const {
    namespace fs = std::filesystem;
    fs::path out_dir(directory);
    if (!out_dir.is_absolute()) out_dir = fs::current_path() / out_dir;
    fs::create_directories(out_dir);

    // First pass: Open all files and write headers
    std::unordered_map<std::string, std::ofstream> files;
    
    for (const auto& [cycle_index, cycle_data] : var_log_) {
      for (const auto& [step_index, step_data] : cycle_data) {
        // Regular variables
        for (const auto& [var_type, entry] : step_data.variables) {
          if (step_data.step_data.times.empty() || entry.values.empty()) continue;
          
          const std::string filename = SanitizeFilename(VariableTypeList[static_cast<size_t>(var_type)]) + ".csv";
          if (files.find(filename) == files.end()) {
            fs::path filepath = out_dir / filename;
            files[filename].open(filepath);
            if (files[filename].is_open()) {
              files[filename] << "Cycle Index,Step Index,time";
              const size_t num_cells = entry.values[0].size();
              for (size_t i = 0; i < num_cells; ++i) files[filename] << ",cell_" << i;
              files[filename] << "\n";
            }
          }
        }

        // Component variables
        for (const auto& [component, comp_vars] : step_data.component_variables) {
          for (const auto& [comp_var_type, entry] : comp_vars) {
            if (step_data.step_data.times.empty() || entry.values.empty()) continue;
            
            const std::string base_name = ComponentVariableTypeList[static_cast<size_t>(comp_var_type)];
            const std::string full_name = component + " " + base_name;
            const std::string filename = SanitizeFilename(full_name) + ".csv";
            if (files.find(filename) == files.end()) {
              fs::path filepath = out_dir / filename;
              files[filename].open(filepath);
              if (files[filename].is_open()) {
                files[filename] << "Cycle Index,Step Index,time";
                const size_t num_cells = entry.values[0].size();
                for (size_t i = 0; i < num_cells; ++i) files[filename] << ",cell_" << i;
                files[filename] << "\n";
              }
            }
          }
        }
      }
    }

    // Second pass: Write all data
    for (const auto& [cycle_index, cycle_data] : var_log_) {
      for (const auto& [step_index, step_data] : cycle_data) {
        const auto& shared_times = step_data.step_data.times;  // Shared times vector
        
        // Regular variables
        for (const auto& [var_type, entry] : step_data.variables) {
          if (shared_times.empty() || entry.values.empty()) continue;
          
          const std::string filename = SanitizeFilename(VariableTypeList[static_cast<size_t>(var_type)]) + ".csv";
          auto& file = files[filename];
          if (!file.is_open()) continue;

          WriteCsvRows(file, cycle_index, step_index, shared_times, entry.values, max_log_points_per_step_);
        }

        // Component variables
        for (const auto& [component, comp_vars] : step_data.component_variables) {
          for (const auto& [comp_var_type, entry] : comp_vars) {
            if (shared_times.empty() || entry.values.empty()) continue;
            
            const std::string base_name = ComponentVariableTypeList[static_cast<size_t>(comp_var_type)];
            const std::string full_name = component + " " + base_name;
            const std::string filename = SanitizeFilename(full_name) + ".csv";
            auto& file = files[filename];
            if (!file.is_open()) continue;

            WriteCsvRows(file, cycle_index, step_index, shared_times, entry.values, max_log_points_per_step_);
          }
        }
      }
    }
  }

  Var1D EquationLog::TakeRow() {
    if (spare_rows_.empty()) return {};
    Var1D row = std::move(spare_rows_.back());
    spare_rows_.pop_back();
    return row;
  }

  Var1D EquationLog::CopyRow(const std::vector<double>& values) {
    if (spare_rows_.empty()) return values;
    Var1D row = TakeRow();
    row.assign(values.begin(), values.end());
    return row;
  }

  void EquationLog::RecycleRows(Log& log) {
    auto recycle = [&](VariableData& data) {
      for (auto& row : data.values) spare_rows_.push_back(std::move(row));
    };
    for (auto& [type, data] : log.variables) recycle(data);
    for (auto& [component, vars] : log.component_variables)
      for (auto& [type, data] : vars) recycle(data);
  }

  void EquationLog::Clear() {
    var_log_.clear();
    spare_rows_.clear();
    spare_rows_.shrink_to_fit();
    current_logged_cycle_ = -1;
    current_logged_step_ = -1;
    sys_data.cycle_index_set_.clear();
    sys_data.step_index_set_.clear();
  }

  void EquationLog::EvictOldCyclesIfNeeded() {
    // Ensure max_cycles_to_keep_ is at least 1
    if (max_cycles_to_keep_ < 1) {
      max_cycles_to_keep_ = 1;
    }
    
    // Count how many cycles we have
    int num_cycles = static_cast<int>(var_log_.size());
    
    // If we're over the limit, remove oldest cycles
    while (num_cycles > max_cycles_to_keep_) {
      // Get the oldest cycle (smallest cycle_index)
      auto oldest_cycle_it = var_log_.begin();
      int oldest_cycle = oldest_cycle_it->first;
      
      // Don't evict the current cycle being logged
      if (oldest_cycle == current_logged_cycle_) {
        break; // Can't evict current cycle, stop
      }
      
      // Remove the cycle, keeping its row buffers for the cycles still to come
      for (auto& [step_index, log] : oldest_cycle_it->second) RecycleRows(log);
      var_log_.erase(oldest_cycle_it);
      sys_data.cycle_index_set_.erase(oldest_cycle);
      num_cycles--;
    }
  }

  const Log& EquationLog::GetLog(int cycle_index, int step_index) const {
    static const Log kEmptyLog;
    auto cycle_it = var_log_.find(cycle_index);
    if (cycle_it == var_log_.end()) {
      return kEmptyLog;
    }
    auto step_it = cycle_it->second.find(step_index);
    if (step_it == cycle_it->second.end()) {
      return kEmptyLog;
    }
    return step_it->second;
  }

  VariableLogEntry EquationLog::GetVariable(int cycle_index, int step_index, VariableType type) const {
    // Check if cycle_index and step_index exist
    auto cycle_it = var_log_.find(cycle_index);
    if (cycle_it == var_log_.end()) {
      return VariableLogEntry{}; // Return empty entry if cycle doesn't exist
    }
    
    auto step_it = cycle_it->second.find(step_index);
    if (step_it == cycle_it->second.end()) {
      return VariableLogEntry{}; // Return empty entry if step doesn't exist
    }
    
    const auto& sim_log = step_it->second;
    auto var_it = sim_log.variables.find(type);
    if (var_it == sim_log.variables.end()) {
      return VariableLogEntry{}; // Return empty entry if variable type doesn't exist
    }
    
    // Reconstruct VariableLogEntry from optimized structure
    VariableLogEntry entry;
    entry.cycle_index = sim_log.step_data.cycle_index;
    entry.step_index = sim_log.step_data.step_index;
    entry.times = sim_log.step_data.times;  // Already double, no conversion needed
    entry.values = var_it->second.values;
    entry.variable_name = var_it->second.variable_name;
    entry.units = var_it->second.units;
    entry.step_name = var_it->second.step_name;
    entry.object_name = var_it->second.object_name;
    
    return entry;
  }

  VariableLogEntry EquationLog::GetComponentVariable(int cycle_index, int step_index,
      ComponentVariableType type, const std::string& component) const 
  {
    // Check if cycle_index and step_index exist
    auto cycle_it = var_log_.find(cycle_index);
    if (cycle_it == var_log_.end()) {
      return VariableLogEntry{}; // Return empty entry if cycle doesn't exist
    }
    
    auto step_it = cycle_it->second.find(step_index);
    if (step_it == cycle_it->second.end()) {
      return VariableLogEntry{}; // Return empty entry if step doesn't exist
    }
    
    const auto& sim_log = step_it->second;
    auto comp_it = sim_log.component_variables.find(component);
    if (comp_it == sim_log.component_variables.end()) {
      return VariableLogEntry{}; // Return empty entry if component doesn't exist
    }
    
    auto type_it = comp_it->second.find(type);
    if (type_it == comp_it->second.end()) {
      return VariableLogEntry{}; // Return empty entry if component variable type doesn't exist
    }
    
    // Reconstruct VariableLogEntry from optimized structure
    VariableLogEntry entry;
    entry.cycle_index = sim_log.step_data.cycle_index;
    entry.step_index = sim_log.step_data.step_index;
    entry.times = sim_log.step_data.times;  // Already double, no conversion needed
    entry.values = type_it->second.values;
    entry.variable_name = type_it->second.variable_name;
    entry.units = type_it->second.units;
    entry.step_name = type_it->second.step_name;
    entry.object_name = type_it->second.object_name;
    
    return entry;
  }

  static std::vector<double> LinearInterpolateCellCenter(const std::vector<double>& data,
    const std::vector<double>& dx)
  {
    assert(data.size() == dx.size());

    std::vector<double> out(data.size() - 1);
    for (int i = 0; i < out.size(); ++i) {
      double rdx1 = dx[i] / (dx[i] + dx[i + 1]);
      double rdx2 = dx[i + 1] / (dx[i] + dx[i + 1]);
      out[i] = rdx1 * data[i] + rdx2 * data[i + 1];
    }
    return out;
  }

  VariableLogEntry EquationLog::GetFaceVariable(int cycle_index, int step_index, VariableType type) const {
    // Get the base variable entry first
    VariableLogEntry cell_centre_entry = GetVariable(cycle_index, step_index, type);
    if (cell_centre_entry.values.empty()) {
      return VariableLogEntry{}; // Return empty if not found
    }
    
    // Check if variable is already at cell faces (e.g., velocity)
    if (cell_centre_entry.values[0].size() == domain.cell_width.size() - 1) {
      // Already face-centered, return as-is
      return cell_centre_entry;
    }
    
    // Interpolate cell-centered values to face values
    const auto& dx = domain.cell_width;
    VariableLogEntry cell_face_entry;
    // Copy all metadata from the original entry
    cell_face_entry.cycle_index = cell_centre_entry.cycle_index;
    cell_face_entry.step_index = cell_centre_entry.step_index;
    cell_face_entry.variable_name = cell_centre_entry.variable_name;
    cell_face_entry.units = cell_centre_entry.units;
    cell_face_entry.step_name = cell_centre_entry.step_name;
    cell_face_entry.object_name = cell_centre_entry.object_name;
    // Copy times and interpolate values
    cell_face_entry.times = cell_centre_entry.times;
    cell_face_entry.values.reserve(cell_centre_entry.values.size());
    
    for (size_t i = 0; i < cell_centre_entry.values.size(); ++i) {
      cell_face_entry.values.push_back(
        LinearInterpolateCellCenter(cell_centre_entry.values[i], dx)
      );
    }

    return cell_face_entry;
  }

  VariableLogEntry EquationLog::GetComponentFaceVariable(int cycle_index, int step_index,
      ComponentVariableType type, const std::string& component) const {
    
    // Get the base component variable entry first
    VariableLogEntry cell_centre_entry = GetComponentVariable(cycle_index, step_index, type, component);
    if (cell_centre_entry.values.empty()) {
      return VariableLogEntry{}; // Return empty if not found
    }
    
    // Check if variable is already at cell faces
    if (cell_centre_entry.values[0].size() == domain.cell_width.size() - 1) {
      return cell_centre_entry;
    }
    
    // Interpolate cell-centered values to face values
    const auto& dx = domain.cell_width;
    VariableLogEntry cell_face_entry;
    // Copy all metadata from the original entry
    cell_face_entry.cycle_index = cell_centre_entry.cycle_index;
    cell_face_entry.step_index = cell_centre_entry.step_index;
    cell_face_entry.variable_name = cell_centre_entry.variable_name;
    cell_face_entry.units = cell_centre_entry.units;
    cell_face_entry.step_name = cell_centre_entry.step_name;
    cell_face_entry.object_name = cell_centre_entry.object_name;
    // Copy times and interpolate values
    cell_face_entry.times = cell_centre_entry.times;
    cell_face_entry.values.reserve(cell_centre_entry.values.size());
    
    for (size_t i = 0; i < cell_centre_entry.values.size(); ++i) {
      cell_face_entry.values.push_back(
        LinearInterpolateCellCenter(cell_centre_entry.values[i], dx)
      );
    }

    return cell_face_entry;
  }

    void BoundaryLog::LogStep(const StreamData& stream, const std::string& stream_name, int cycle_idx, int step_idx) {
    BoundaryLogEntry entry;
    entry.cycle_index = cycle_idx;
    entry.step_index = step_idx;
    entry.stream_data = stream;
    log_[cycle_idx][step_idx][stream_name] = std::move(entry);
  }

  void BoundaryLog::Clear() { 
    log_.clear();
    sys_data.cycle_index_set_.clear();
    sys_data.step_index_set_.clear();
    sys_data.step_list_.clear();
    stream_names.clear();
    current_logged_cycle_ = -1;
  }

  void BoundaryLog::LogEndOfStep(const Step& step, int cycle_idx, int step_idx) {
    // Check if we're starting a new cycle
    bool is_new_cycle = (current_logged_cycle_ != cycle_idx);
    if (is_new_cycle) {
      current_logged_cycle_ = cycle_idx;
      EvictOldCyclesIfNeeded();
    }
    
    auto streams = GetAllBoundaryStreams(step);
    log_[cycle_idx][step_idx];
    sys_data.cycle_index_set_.insert(cycle_idx);
    sys_data.step_index_set_.insert(step_idx);
    for (const auto& stream : streams) {
      LogStep(*stream.stream_data, stream.boundary_name, cycle_idx, step_idx);
    }
  }

  void BoundaryLog::Initialise(UserData& user_data) {
    auto& settings = user_data.settings;
    auto& fluid = user_data.fluid;
    auto& bed = user_data.bed;
    auto& cycle = user_data.cycle;

    max_cycles_to_keep_ = std::max(1, settings.max_cycles_to_keep_);
    
    Clear();
    for (const auto& step : cycle.step_list_) {
      sys_data.step_list_.emplace_back(step.step_name_);
      auto streams = GetAllBoundaryStreams(step); 
      auto& names = stream_names[step.step_name_];
      for (const auto& stream : streams) {
        names.push_back(stream.boundary_name);
      }
    }

    sys_data.components_ = fluid.components_;
    sys_data.adsorbent_mass = GetAdsorbentMass(bed, user_data.adsorbent_directory);
    sys_data.total_cycle_time = cycle.GetTotalSimulationTime();
    sys_data.molar_mass = fluid.molar_mass_;
  }

  void BoundaryLog::EvictOldCyclesIfNeeded() {
    if (max_cycles_to_keep_ < 1) {
      max_cycles_to_keep_ = 1;
    }
    
    // Count how many cycles we have
    int num_cycles = static_cast<int>(log_.size());
    
    // If we're over the limit, remove oldest cycles
    while (num_cycles > max_cycles_to_keep_) {
      // Get the oldest cycle (smallest cycle_index)
      auto oldest_cycle_it = log_.begin();
      int oldest_cycle = oldest_cycle_it->first;
      
      // Don't evict the current cycle being logged
      if (oldest_cycle == current_logged_cycle_) {
        break; // Can't evict current cycle, stop
      }
      
      // Remove the cycle
      log_.erase(oldest_cycle_it);
      sys_data.cycle_index_set_.erase(oldest_cycle);
      num_cycles--;
    }
  }

  BoundaryVariableData BoundaryLog::GetVariableData(
    int cycle_idx, 
    int step_idx, 
    const std::string& stream_name, 
    BoundaryVariableType variable
  ) const {
    BoundaryVariableData out;

    // Check if cycle_index and step_index exist
    auto cycle_it = log_.find(cycle_idx);
    if (cycle_it == log_.end()) {
      return BoundaryVariableData{}; // Return empty entry if cycle doesn't exist
    }
    
    auto step_it = cycle_it->second.find(step_idx);
    if (step_it == cycle_it->second.end()) {
      return BoundaryVariableData{}; // Return empty entry if step doesn't exist
    }
    
    const auto& sim_log = step_it->second;
    auto var_it = sim_log.find(stream_name);
    if (var_it == sim_log.end()) {
      return BoundaryVariableData{}; // Return empty entry if variable type doesn't exist
    }

    const auto& entry = log_.at(cycle_idx).at(step_idx).at(stream_name); 
    out.variable_timeseries = GetBoundaryVariableVector(variable, entry.stream_data);
    out.time = entry.stream_data.time_series_;
    out.units = GetVariableUnits(variable);
    out.variable_name = kBoundaryVariableTypeList[static_cast<int>(variable)];
    out.cycle_idx = cycle_idx;
    // out.step_name = sys_data.step_list_.at(step_idx);
    out.step_name = sys_data.step_list_[step_idx];
    out.object_name = stream_name;
    return out;
  }

  BoundaryVariableData BoundaryLog::GetComponentVariableData(
    int cycle_idx, 
    int step_idx, 
    const std::string& stream_name, 
    BoundaryComponentVariableType variable, 
    const std::string& component
  ) const {
    BoundaryVariableData out;
    
    // Check if cycle_index and step_index exist
    auto cycle_it = log_.find(cycle_idx);
    if (cycle_it == log_.end()) {
      return BoundaryVariableData{}; // Return empty entry if cycle doesn't exist
    }
    
    auto step_it = cycle_it->second.find(step_idx);
    if (step_it == cycle_it->second.end()) {
      return BoundaryVariableData{}; // Return empty entry if step doesn't exist
    }
    
    const auto& sim_log = step_it->second;
    auto var_it = sim_log.find(stream_name);
    if (var_it == sim_log.end()) {
      return BoundaryVariableData{}; // Return empty entry if variable type doesn't exist
    }

    const auto& entry = log_.at(cycle_idx).at(step_idx).at(stream_name); 

    out.variable_timeseries = GetBoundaryComponentVariableVector(variable, entry.stream_data, component);
    out.time = entry.stream_data.time_series_;
    out.units = GetComponentVariableUnits(variable);
    out.variable_name = kBoundaryComponentVariableTypeList[static_cast<int>(variable)];
    out.cycle_idx = cycle_idx;
    // out.step_name = sys_data.step_list_.at(step_idx);
    out.step_name = sys_data.step_list_[step_idx];
    out.object_name = stream_name;
    return out;
  }

  std::vector<double> GetBoundaryVariableVector(BoundaryVariableType variable, const StreamData& stream_data) {
    switch (variable) {
    case BoundaryVariableType::TEMPERATURE: return stream_data.temperature_boundary_;
    case BoundaryVariableType::MOLAR_FLOW_RATE: return stream_data.molar_flow_;
    case BoundaryVariableType::COMPRESSOR_POWER: return stream_data.electrical_power_;
    default: return std::vector<double>();
    }
  }

  [[nodiscard]]
  const std::vector<double>& GetBoundaryComponentVariableVector(
      BoundaryComponentVariableType variable,
      const StreamData& stream_data,
      const std::string& component) noexcept
  {
    static const std::vector<double> kEmpty;

    switch (variable) {
    case BoundaryComponentVariableType::MOLE_FRACTION: {
      const auto it = stream_data.mole_fraction_boundary_.find(component);
      return (it != stream_data.mole_fraction_boundary_.end()) ? it->second : kEmpty;
    }
    case BoundaryComponentVariableType::MOLAR_FLOW: {
      const auto it = stream_data.component_molar_flow_.find(component);
      return (it != stream_data.component_molar_flow_.end()) ? it->second : kEmpty;
    }
    }

    return kEmpty;
  }

  std::string GetVariableUnits(BoundaryVariableType variable) {
    switch (variable) {
      case BoundaryVariableType::TEMPERATURE: return "[K]";
      case BoundaryVariableType::MOLAR_FLOW_RATE: return "[mol/s]";
      case BoundaryVariableType::COMPRESSOR_POWER: return "[W]";
      default: return "[-]";
    }
  }
  
  std::string GetComponentVariableUnits(BoundaryComponentVariableType variable) {
    switch (variable) {
      case BoundaryComponentVariableType::MOLE_FRACTION: return "[mol/mol]";
      case BoundaryComponentVariableType::MOLAR_FLOW: return "[mol/s]";
      default: return "[-]";
    }
  }

  std::string GetVariableUnits(VariableType variable) {
    switch (variable) {
      case VariableType::PRESSURE: return "[Pa]";
      case VariableType::TEMPERATURE: return "[K]";
      case VariableType::WALL_TEMPERATURE: return "[K]";
      case VariableType::VELOCITY: return "[m/s]";
      case VariableType::VOLUMETRIC_FLOW_RATE: return "[m3/s]";
      case VariableType::MOLAR_FLOW_RATE: return "[mol/s]";
      case VariableType::TOTAL_MOLE_FRACTION: return "[mol/mol]";
      case VariableType::APPROACH_TO_FLUIDISATION: return "[-]";
      case VariableType::VISCOSITY: return "[Pa.s]";
      case VariableType::FLUID_DENSITY: return "[kg/m^3]";
      case VariableType::MASS_SOURCE: return "[mol/(m^3.s)]";
      case VariableType::ENERGY_SOURCE: return "[J/(m^3.s)]";
      case VariableType::MOLAR_CONCENTRATION: return "[mol/m^3]";
      case VariableType::HEAT_CAPACITY: return "[J/(kg.K)]";
      case VariableType::MOLAR_MASS: return "[kg/mol]";
      case VariableType::THERMAL_CONDUCTIVITY: return "[W/(m.K)]";
      case VariableType::HEAT_DISPERSION: return "[W/(m.K)]";
      case VariableType::SCHMIDT_NUMBER: return "[-]";
      case VariableType::REYNOLDS_NUMBER: return "[-]";
      case VariableType::SHERWOOD_NUMBER: return "[-]";
      case VariableType::PRESSURE_GRADIENT: return "[Pa/m]";
      case VariableType::INTERNAL_HEAT_TRANSFER_COEFFICIENT: return "[W/(m^2.K)]";
      case VariableType::EXTERNAL_HEAT_TRANSFER_COEFFICIENT: return "[W/(m^2.K)]";
      case VariableType::AXIAL_DISPERSION_COEFFICIENT: return "[m^2/s]";
      case VariableType::SOLID_TEMPERATURE: return "[K]";
      case VariableType::PARTICLE_HEAT_TRANSFER_COEFFICIENT: return "[W/(m^3.K)]";
      default: return "[-]";
    }
  }
  
  std::string GetComponentVariableUnits(ComponentVariableType variable) {
    switch (variable) {
      case ComponentVariableType::MOLE_FRACTION: return "[mol/mol]";
      case ComponentVariableType::ADSORBED_CONCENTRATION: return "[mol/kg]";
      case ComponentVariableType::EQUILIBRIUM_ADSORBED_CONCENTRATION: return "[mol/kg]";
      case ComponentVariableType::COMPONENT_MASS_SOURCE: return "[kg/(m^3.s)]";
      case ComponentVariableType::COMPONENT_ENERGY_SOURCE: return "[J/(m^3.s)]";
      case ComponentVariableType::COMPONENT_MOLAR_CONCENTRATION: return "[mol/m^3]";
      case ComponentVariableType::MASS_TRANSFER_COEFFICIENT: return "[1/s]";
      case ComponentVariableType::MOLECULAR_DIFFUSION: return "[m^2/s]";
      default: return "[-]";
    }
  }

} // end adsim namespace
