// log.h
#pragma once

// Adsim includes
#include "adsim/equations/equation.h"
#include "adsim/boundary_conditions/boundary_stream.h"

// STL includes
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <set>
#include <iterator>

// SUNDIALS includes
#include <ida/ida.h>


#if !defined(NDEBUG) || defined(_DEBUG)
  #include <iostream>
  #define LOG(msg) \
    do { std::cerr << msg << '\n'; } while (0)
#else
  #define LOG(msg) ((void)0)
#endif

namespace adsim {

  class Bed;
  class Cycle;
  class Fluid;

  enum class SystemObject {
    Bed,
    FlowBoundary,
    KPI,
    VoidSpace
  };

  inline constexpr std::string_view kSystemObjectList[] = {
    "Bed",
    "Flow Boundary",
    "KPI",
    "Void Space"
  };

  enum class BoundaryVariableType {
    TEMPERATURE,
    MOLAR_FLOW_RATE,
    COMPRESSOR_POWER
  };
  
  inline constexpr std::string_view kBoundaryVariableTypeList[] = {
    "Temperature", 
    "Molar Flow Rate", 
    "Compressor Power"
  };
  
  enum class BoundaryComponentVariableType {
    MOLE_FRACTION,
    MOLAR_FLOW
  };

  inline constexpr std::string_view kBoundaryComponentVariableTypeList[] = {
    "Mole Fraction",
    "Molar Flow"
  };


  enum class VariableType {
    PRESSURE,
    TEMPERATURE,
    WALL_TEMPERATURE,
    VELOCITY,
    VOLUMETRIC_FLOW_RATE,
    MOLAR_FLOW_RATE,
    TOTAL_MOLE_FRACTION,
    APPROACH_TO_FLUIDISATION,
    VISCOSITY,
    FLUID_DENSITY,
    MASS_SOURCE,
    ENERGY_SOURCE,
    MOLAR_CONCENTRATION,
    HEAT_CAPACITY,
    MOLAR_MASS,
    THERMAL_CONDUCTIVITY,
    HEAT_DISPERSION,
    SCHMIDT_NUMBER,
    REYNOLDS_NUMBER,
    SHERWOOD_NUMBER,
    PRESSURE_GRADIENT,
    INTERNAL_HEAT_TRANSFER_COEFFICIENT,
    EXTERNAL_HEAT_TRANSFER_COEFFICIENT,
    AXIAL_DISPERSION_COEFFICIENT,
    SOLID_TEMPERATURE,
    PARTICLE_HEAT_TRANSFER_COEFFICIENT
  };

  enum class ComponentVariableType {
    MOLE_FRACTION,
    ADSORBED_CONCENTRATION,
    EQUILIBRIUM_ADSORBED_CONCENTRATION,
    COMPONENT_MASS_SOURCE,
    COMPONENT_ENERGY_SOURCE,
    COMPONENT_MOLAR_CONCENTRATION,
    MASS_TRANSFER_COEFFICIENT,
    MOLECULAR_DIFFUSION
  };

  // Variables that only have physical meaning in packed-bed cells; the plotter
  // hides these when the user is viewing a void cell.
  inline constexpr bool IsPorousOnlyVariable(VariableType t) {
    return t == VariableType::SOLID_TEMPERATURE
        || t == VariableType::PARTICLE_HEAT_TRANSFER_COEFFICIENT;
  }

  inline constexpr bool IsPorousOnlyComponentVariable(ComponentVariableType t) {
    return t == ComponentVariableType::ADSORBED_CONCENTRATION
        || t == ComponentVariableType::EQUILIBRIUM_ADSORBED_CONCENTRATION
        || t == ComponentVariableType::MASS_TRANSFER_COEFFICIENT;
  }

  // Order must match the above
  static const char* VariableTypeList[] = {
    "Pressure",
    "Fluid Temperature",
    "Wall Temperature",
    "Velocity",
    "Volumetric Flow Rate",
    "Molar Flow Rate",
    "Total Mole Fraction",
    "Approach To Fluidisation",
    "Viscosity",
    "Fluid Density",
    "Mass Source",
    "Energy Source",
    "Molar Concentration",
    "Heat Capacity",
    "Molar Mass",
    "Thermal Conductivity",
    "Heat Dispersion",
    "Schmidt Number",
    "Reynolds Number",
    "Sherwood Number",
    "Pressure Gradient",
    "Internal Heat Transfer Coefficient",
    "External Heat Transfer Coefficient",
    "Axial Dispersion Coefficient",
    "Solid Temperature",
    "Particle Heat Transfer Coefficient"
  };

  static const char* ComponentVariableTypeList[] = {
    "Mole Fraction",
    "Adsorbed Concentration",
    "Equilibrium Adsorbed Concentration",
    "Component Mass Source",
    "Component Energy Source",
    "Component Molar Concentration",
    "Mass Transfer Coefficient",
    "Molecular Diffusion"
  };

  // Indexed by ComponentVariableType: a spare name would show up in the UI as
  // an option that maps to no enum value.
  static_assert(std::size(ComponentVariableTypeList) ==
                    static_cast<size_t>(ComponentVariableType::MOLECULAR_DIFFUSION) + 1,
                "ComponentVariableTypeList must have one entry per ComponentVariableType");

  std::string GetVariableUnits(VariableType variable);
  std::string GetComponentVariableUnits(ComponentVariableType variable);

  // BoundaryLogEntry stores the data for an entire step
  struct BoundaryLogEntry {
    StreamData stream_data{};
    int cycle_index = 0;
    int step_index = 0;
  };

  struct BoundaryVariableData {
    // [time]
    std::vector<double> variable_timeseries;
    std::vector<double> time;
    std::string variable_name;
    std::string units;
    int cycle_idx;
    std::string step_name;
    std::string object_name;
  };

  struct SystemData {
    std::set<int> cycle_index_set_;
    std::set<int> step_index_set_;
    std::vector<std::string> step_list_;
    std::vector<std::string> components_;
    std::unordered_map<std::string, double> molar_mass;
    double adsorbent_mass;
    double total_cycle_time;
  };

  struct BoundaryLog {
  public:
    void LogEndOfStep(const Step& step, int cycle_idx, int step_idx);
    void Initialise(UserData& user_data);
    BoundaryVariableData GetVariableData(int cycle_idx, int step_idx, const std::string& stream_name, BoundaryVariableType variable) const;
    BoundaryVariableData GetComponentVariableData(int cycle_idx, int step_idx, const std::string& stream_name, BoundaryComponentVariableType variable, const std::string& component) const;
    
    int max_cycles_to_keep_ = 10;
    SystemData sys_data;
    std::unordered_map<std::string, std::vector<std::string>> stream_names{};
    
  private:
    // log_[cycle_idx][step_idx]["Flow_Boundary_Name"]
    std::map<int, std::map<int, std::map<std::string, BoundaryLogEntry>>> log_;
    
    int current_logged_cycle_ = -1;     // Track which cycle we're currently logging
    void LogStep(const StreamData& state, const std::string& stream_name, int cycle_idx, int step_idx);
    void Clear();
    void EvictOldCyclesIfNeeded();
    // The logged entry of one stream, or nullptr when the cycle, step or
    // stream was not logged.
    const BoundaryLogEntry* FindEntry(int cycle_idx, int step_idx, const std::string& stream_name) const;
    // One variable's time series of `entry` with its description filled in.
    BoundaryVariableData MakeVariableData(const BoundaryLogEntry& entry, const std::string& stream_name,
        std::vector<double> values, std::string variable_name, std::string units) const;

  };

  std::vector<double> GetBoundaryVariableVector(
    BoundaryVariableType variable, 
    const StreamData& stream_data);
  
  [[nodiscard]]
  const std::vector<double>& GetBoundaryComponentVariableVector(
      BoundaryComponentVariableType variable,
      const StreamData& stream_data,
      const std::string& component) noexcept;

  std::string GetVariableUnits(BoundaryVariableType variable);
  std::string GetComponentVariableUnits(BoundaryComponentVariableType variable);

  // Memory optimization: Use float instead of double (50% reduction)
  using Var0D = double;
  using Var1D = std::vector<double>;
  // using Var0D = float;
  // using Var1D = std::vector<float>;

  // Step-level shared data (times stored once per step, not per variable)
  struct StepData {
    std::vector<double> times{};  // Shared times vector - stored once per step (as double)
    int cycle_index = 0;
    int step_index = 0;
  };

  // Internal storage structure (optimized)
  struct VariableData {
    std::vector<Var1D> values{};  // Store values as double; times stored in StepData
    std::string variable_name;
    std::string units;
    std::string step_name;
    std::string object_name;
  };

  // Public API structure (reconstructed on access, maintains backward compatibility)
  struct VariableLogEntry {
    int cycle_index = 0;
    int step_index = 0;
    std::vector<double> times{};  // Shared from StepData (as double)
    std::vector<std::vector<double>> values{};

    std::string variable_name;
    std::string units;
    std::string step_name;
    std::string object_name;
  };

  using VariableLog = std::unordered_map<VariableType, VariableData>;
  using ComponentVariableLog = std::unordered_map<std::string, std::unordered_map<ComponentVariableType, VariableData>>;

  struct Log {
    StepData step_data;  // Shared metadata including times
    VariableLog variables;
    ComponentVariableLog component_variables;
  };

  // [cycle_index][step_index]
  using SimulationLog = std::map<int, std::map<int, Log>>;

  class EquationLog {
  public:
    void Initialise(UserData& user_data);   
    void LogState(const Equations& state, double t, int step_index, int cycle_index);
    void SaveToFile(const std::string& directory = "Results") const;
    const SimulationLog& GetSimulationLog() const { return var_log_; }
    const Log& GetLog(int cycle_index, int step_index) const;
    const DomainVectors& GetDomain() const { return domain; }

    VariableLogEntry GetVariable(int cycle_index, int step_index, VariableType type) const;
    VariableLogEntry GetComponentVariable(int cycle_index, int step_index,
      ComponentVariableType type, const std::string& component) const;
    VariableLogEntry GetFaceVariable(int cycle_index, int step_index, VariableType type) const;
    VariableLogEntry GetComponentFaceVariable(int cycle_index, int step_index,
      ComponentVariableType type, const std::string& component) const;

    int max_cycles_to_keep_ = 10;
    // Max rows written per (cycle, step) in SaveToFile (CSV plotting output
    // only; the in-memory log stays full-resolution). <= 0 writes every row.
    int max_log_points_per_step_ = 0;
    SystemData sys_data;

  private:
    void Clear();
    void EvictOldCyclesIfNeeded();
    std::string SanitizeFilename(const std::string& name) const;

    // Row buffers of evicted cycles, reused by LogState so steady-state logging
    // neither allocates nor frees a vector per variable per solver step.
    Var1D TakeRow();
    Var1D CopyRow(const std::vector<double>& values);
    void RecycleRows(Log& log);
    std::vector<Var1D> spare_rows_;
  
    SimulationLog var_log_;
    DomainVectors domain;
    int current_logged_cycle_ = -1;     // Track which cycle we're currently logging
    int current_logged_step_ = -1;      // Track which step we're currently logging

  };

  class Logger {
  public:
    EquationLog eq_log_;
    BoundaryLog b_log_;
    bool logging_enabled_ = false;
    int max_cycles_to_keep_ = 10;

  };

  inline void PrintData(sunrealtype* y_data, size_t N) {
    for (int i = 0; i < N; ++i) {
      std::cout << y_data[i] << " ";
    }
    std::cout << std::endl;
  }

  inline void LogTime(double t) {
    std::cout << "t = " << t << "\n";
  }

} // end adsim namespace
