// log.h
#pragma once

// Adsim includes
#include "adsim/equations/equation.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/solver/kpi_quadrature.h"

// STL includes
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <set>
#include <iterator>
#include <cstddef>
#include <functional>

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
    COMPRESSOR_POWER,
    PRESSURE  // the bed-end face pressure (a pressure boundary's set value)
  };
  
  inline constexpr std::string_view kBoundaryVariableTypeList[] = {
    "Temperature", 
    "Molar Flow Rate", 
    "Compressor Power",
    "Pressure"
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
    PARTICLE_HEAT_TRANSFER_COEFFICIENT,
    // Heat of reaction released per bed volume, Q_rxn (logged only when the
    // deck has reactions; Equations::reaction_heat).
    REACTION_HEAT
  };

  enum class ComponentVariableType {
    MOLE_FRACTION,
    ADSORBED_CONCENTRATION,
    EQUILIBRIUM_ADSORBED_CONCENTRATION,
    COMPONENT_MASS_SOURCE,
    COMPONENT_ENERGY_SOURCE,
    COMPONENT_MOLAR_CONCENTRATION,
    MASS_TRANSFER_COEFFICIENT,
    MOLECULAR_DIFFUSION,
    // Slow-site loading q2 of a dual-rate LDF component (logged only when the
    // deck has one; 0 for the components without a slow site).
    SLOW_ADSORBED_CONCENTRATION,
    // Moles of the component made by the reactions per bed volume and time,
    // S_i (logged only when the deck has reactions). Not part of
    // COMPONENT_MASS_SOURCE, which stays the adsorption source.
    REACTION_SOURCE
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
        || t == ComponentVariableType::SLOW_ADSORBED_CONCENTRATION
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
    "Particle Heat Transfer Coefficient",
    "Reaction Heat"
  };

  static_assert(std::size(VariableTypeList) == static_cast<size_t>(VariableType::REACTION_HEAT) + 1,
                "VariableTypeList must have one entry per VariableType");

  static const char* ComponentVariableTypeList[] = {
    "Mole Fraction",
    "Adsorbed Concentration",
    "Equilibrium Adsorbed Concentration",
    "Component Mass Source",
    "Component Energy Source",
    "Component Molar Concentration",
    "Mass Transfer Coefficient",
    "Molecular Diffusion",
    "Slow Adsorbed Concentration",
    "Reaction Source"
  };

  // Indexed by ComponentVariableType: a spare name would show up in the UI as
  // an option that maps to no enum value.
  static_assert(std::size(ComponentVariableTypeList) ==
                    static_cast<size_t>(ComponentVariableType::REACTION_SOURCE) + 1,
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
    // The logged data of one endpoint over one step (time from the step
    // start), or nullptr when the cycle, step or endpoint was not logged.
    const StreamData* GetStreamData(int cycle_idx, int step_idx, const std::string& stream_name) const;
    // The mole balance of one step (KPI quadratures and column inventory),
    // kept for the same cycles as the stream logs; nullptr when not logged.
    void LogStepBalance(int cycle_idx, int step_idx, StepBalance balance);
    const StepBalance* GetStepBalance(int cycle_idx, int step_idx) const;
    
    int max_cycles_to_keep_ = 1;
    SystemData sys_data;
    std::unordered_map<std::string, std::vector<std::string>> stream_names{};
    
  private:
    // log_[cycle_idx][step_idx]["Flow_Boundary_Name"]
    std::map<int, std::map<int, std::map<std::string, BoundaryLogEntry>>> log_;
    // step_balances_[cycle_idx][step_idx]
    std::map<int, std::map<int, StepBalance>> step_balances_;
    
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

  // Memory budget of the equation log per (cycle, step) [bytes]. LogState
  // records a full snapshot at every accepted solver step; a step that takes
  // very many steps (a stalled or very stiff run) would otherwise grow without
  // bound (about 100 kB per solver step at 320 cells). When a step's rows reach
  // the budget, its interior rows are thinned (see EquationLog::LogState). Set
  // above the largest bundled deck (Streb2021_A1: 16,672 rows x ~65 kB, about
  // 1.07 GB in one step), so no bundled deck is thinned.
  inline constexpr std::size_t kMaxEquationLogBytesPerStep = std::size_t(2) << 30;  // 2 GiB

  // Step-level shared data (times stored once per step, not per variable)
  struct StepData {
    std::vector<double> times{};  // Shared times vector - stored once per step (as double)
    int cycle_index = 0;
    int step_index = 0;
    std::size_t row_bytes = 0;    // bytes of one logged row (every variable + its time)
    int thinnings = 0;            // times the step's rows were halved to stay in budget
    // Budget bookkeeping (EquationLog::ThinIfOverBudget): LogState calls in
    // this step, the stride of the kept rows in calls (1 until the first
    // thinning), and whether the last row is off that stride (kept only
    // because it is the newest, replaced by the next call).
    std::size_t calls = 0;
    std::size_t stride = 1;
    bool last_provisional = false;
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
  // Per-adsorbent series, [adsorbent id][component][type]: the loading of each
  // solid of a mixed layer on its own (see IsPerAdsorbentVariable).
  using AdsorbentVariableLog = std::map<AdsorbentID, ComponentVariableLog>;

  // The component variables that are also logged per adsorbent in a bed with a
  // mixed layer: the adsorbed and equilibrium loading, the mass-transfer
  // coefficient and the slow-site loading [per kg of that solid].
  inline constexpr bool IsPerAdsorbentVariable(ComponentVariableType t) {
    return t == ComponentVariableType::ADSORBED_CONCENTRATION
        || t == ComponentVariableType::EQUILIBRIUM_ADSORBED_CONCENTRATION
        || t == ComponentVariableType::MASS_TRANSFER_COEFFICIENT
        || t == ComponentVariableType::SLOW_ADSORBED_CONCENTRATION;
  }

  // Display name of a per-adsorbent series, e.g. "CarbonDioxide Adsorbed
  // Concentration (13X)": the component series' name with the adsorbent's
  // label in brackets.
  std::string PerAdsorbentVariableName(const std::string& component, ComponentVariableType type,
                                       const std::string& adsorbent_label);

  // Heat put into the column during one step, as the solver applied it, in two
  // separate parts. All of it is 0 in an isothermal run, where neither the
  // energy sources nor the wall balance act.
  //  * Energy sources (heaters and coolers in the bed): constant powers over
  //    the step, set when the step starts.
  //  * The step's heating medium (jacket) on the wall's outer face: the heat
  //    it passes into and out of the wall, integrated over the step. Per cell
  //    of outer surface the flow into the wall is pi d_o dx h_m (T_m(t) - T_w);
  //    cells where it is positive add to medium_heat_in, the others to
  //    medium_heat_out. Integrated with the trapezoidal rule over every
  //    logged point (each accepted solver step and any fill-in point), as the
  //    log is written, so thinning the stored rows does not change it. Both are
  //    0 for a step without a medium.
  //    The medium's internal heat exchanger, when it has one, adds its heat
  //    to the same two totals: per bed cell a_x h_x dV (T_m(t) - T_bed), with
  //    T_bed the gas temperature (LTE) or the eb / (1 - eb) mix of gas and
  //    solid (LTNE), signed and integrated the same way. Its share is also
  //    kept on its own (exchanger_heat_in / exchanger_heat_out).
  struct StepHeatDuty {
    double heating_power = 0.0;    // [W] sum of the positive energy-source powers
    double cooling_power = 0.0;    // [W] magnitude of the negative ones (>= 0)
    double medium_heat_in = 0.0;   // [J] heat the medium gave the column (jacket + exchanger) over the step
    double medium_heat_out = 0.0;  // [J] heat the column gave the medium over the step (>= 0)
    double exchanger_heat_in = 0.0;   // [J] the internal exchanger's share of medium_heat_in
    double exchanger_heat_out = 0.0;  // [J] the internal exchanger's share of medium_heat_out
  };

  // Instantaneous heat flow [W] from the running step's heating medium into the
  // column, split into its positive (in) and negative (out, as a magnitude)
  // parts: over the cells of outer wall surface (jacket) and over the bed
  // cells (internal exchanger). in / out are the totals, exchanger_in /
  // exchanger_out the exchanger's share. All 0 with no medium or in an
  // isothermal run. `time_in_step` is seconds since the step started;
  // `non_local` is true under LTNE (the exchanger's heat is split between gas
  // and solid).
  struct MediumHeatFlow {
    double in = 0.0;
    double out = 0.0;
    double exchanger_in = 0.0;
    double exchanger_out = 0.0;
  };
  MediumHeatFlow HeatingMediumHeatFlow(const Equations& state, double time_in_step, bool non_local);

  // Rate of one reaction on one solid, keyed (adsorbent id, reaction id).
  using ReactionVariableKey = std::pair<AdsorbentID, std::size_t>;
  using ReactionVariableLog = std::map<ReactionVariableKey, VariableData>;

  // Display name of a reaction's rate series, "SMR Rate (Ni catalyst)": the
  // reaction's name, then the label of the solid that carries it.
  std::string ReactionRateVariableName(const std::string& reaction, const std::string& adsorbent_label);
  inline constexpr const char* kReactionRateUnits = "[mol/(kg.s)]";

  struct Log {
    StepData step_data;  // Shared metadata including times
    VariableLog variables;
    ComponentVariableLog component_variables;
    // Logged only for a deck with reactions: eta r of each reaction on each
    // solid in every cell [mol/(kg of that solid s)], 0 where it is absent.
    ReactionVariableLog reaction_variables;
    // Moles of each component made by the reactions in the bed over the step,
    // integral of sum_n dV_n S_i,n dt [mol] (trapezoidal rule over every
    // logged point, as the log is written; empty without reactions). The
    // mass-balance CSS tests count it as generation.
    std::map<std::string, double> reaction_generation;
    // Logged only when some cell holds more than one solid (a mixed layer):
    // otherwise each cell's aggregate in component_variables already is the
    // loading of its one solid.
    AdsorbentVariableLog adsorbent_variables;
    StepHeatDuty heat_duty;  // recorded with the step's first log; evicted with the cycle
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
    // The per-adsorbent series of `type` (IsPerAdsorbentVariable) for one
    // adsorbent: the loading of that solid [per kg of it] in each cell, 0 in
    // the cells without it, ghost cells mirroring the first and last bed cell
    // (as the aggregate). When the run logged no per-adsorbent series (no
    // mixed layer) it is read off the aggregate, which in every cell holding
    // that adsorbent is its loading. Empty when not logged.
    VariableLogEntry GetAdsorbentComponentVariable(int cycle_index, int step_index, ComponentVariableType type,
      const std::string& component, AdsorbentID adsorbent) const;
    VariableLogEntry GetFaceVariable(int cycle_index, int step_index, VariableType type) const;
    // eta r of reaction `reaction` on adsorbent `adsorbent` (Log::
    // reaction_variables); empty when not logged.
    VariableLogEntry GetReactionRate(int cycle_index, int step_index, AdsorbentID adsorbent,
                                     std::size_t reaction) const;
    // Moles of `component` made by the reactions over one step [mol]
    // (Log::reaction_generation); 0 when none or not logged.
    double ReactionGeneration(int cycle_index, int step_index, const std::string& component) const;

    // The adsorbents of the bed (id -> label, the adsorbent's name, made
    // unique with " #<id>" when two share one), set by Initialise.
    const std::map<AdsorbentID, std::string>& AdsorbentLabels() const { return adsorbent_labels_; }
    // True when this run logs per-adsorbent series (some cell holds more than
    // one solid).
    bool LogsPerAdsorbent() const { return per_adsorbent_; }
    // True for an LTNE run (separate solid temperature).
    bool NonLocal() const { return non_local_; }

    // One reaction on one solid whose rate the run logs (GetReactionRate),
    // in the solver's order; set by Initialise, empty without reactions.
    struct LoggedReaction {
      AdsorbentID adsorbent = 0;
      std::size_t reaction = 0;    // ReactionID
      std::string reaction_name;
      std::string adsorbent_label;  // AdsorbentLabels()
    };
    const std::vector<LoggedReaction>& Reactions() const { return reactions_; }
    VariableLogEntry GetComponentFaceVariable(int cycle_index, int step_index,
      ComponentVariableType type, const std::string& component) const;

    int max_cycles_to_keep_ = 1;
    // Max rows written per (cycle, step) in SaveToFile (CSV plotting output
    // only; the in-memory log keeps every solver step, up to the budget
    // below). <= 0 writes every row.
    int max_log_points_per_step_ = 0;
    // Memory budget per (cycle, step) [bytes]; kMaxEquationLogBytesPerStep
    // (runtime only, not a deck field; tests lower it). 0 disables thinning.
    std::size_t max_bytes_per_step_ = kMaxEquationLogBytesPerStep;
    // Receives one message per (cycle, step) the first time that step reaches
    // the budget (the solvers forward it to their LOW-level log). May be unset.
    std::function<void(const std::string&)> message_sink_;
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
    // Keep a step within max_bytes_per_step_: see the definition in log.cpp.
    void ThinIfOverBudget(Log& log);
    std::vector<Var1D> spare_rows_;
  
    SimulationLog var_log_;
    DomainVectors domain;
    int current_logged_cycle_ = -1;     // Track which cycle we're currently logging
    int current_logged_step_ = -1;      // Track which step we're currently logging
    // Previous logged point of the current step, for integrating the medium duty.
    double medium_prev_t_ = 0.0;
    MediumHeatFlow medium_prev_flow_{};
    // sum_n dV_n S_i,n at the previous logged point of the current step
    // (reaction_generation's trapezoid), per component.
    std::map<std::string, double> reaction_prev_rate_;
    bool non_local_ = false;  // LTNE run (set in Initialise)
    bool per_adsorbent_ = false;  // some cell holds more than one solid (set in Initialise / per step)
    std::map<AdsorbentID, std::string> adsorbent_labels_;
    std::vector<LoggedReaction> reactions_;

  };

  class Logger {
  public:
    EquationLog eq_log_;
    BoundaryLog b_log_;
    bool logging_enabled_ = false;
    int max_cycles_to_keep_ = 1;

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
