/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// KPIs of every cycle of a run.
//
// The equation and boundary logs keep only the last settings.max_cycles_to_keep
// cycles (their spatial profiles are large), so KPIs computed from them on
// demand cover only those cycles. The solvers therefore evaluate the KPIs of
// every tracked component when a cycle completes and keep the few numbers in a
// KPIHistory, which is never evicted: a run's KPI history spans all its cycles.
// A KPIRecord is a few dozen doubles per tracked component, so a 2000-cycle run
// holds well under 1 MB. The CSS check of every cycle is kept the same way
// (CSSHistory, cyclic_steady_state.h).

#include <limits>
#include <map>
#include <string>
#include <vector>

#include "adsim/solver/kpi_config.h"

namespace adsim {

  class EquationLog;
  struct BoundaryLog;

  // KPIs of one tracked component at one cycle. A KPI that is not defined for
  // the run (no product drawn, zero denominator) is NaN, and so is every KPI
  // when the component's KPI configuration is incomplete (`valid` false: no
  // feed or no product endpoint).
  struct KPIRecord {
    int cycle = 0;
    std::string component;
    bool valid = false;
    double purity = std::numeric_limits<double>::quiet_NaN();               // [mol/mol]
    double recovery = std::numeric_limits<double>::quiet_NaN();             // [mol/mol]
    double productivity_mass = std::numeric_limits<double>::quiet_NaN();    // [kg/(kg s)]
    double productivity_volume = std::numeric_limits<double>::quiet_NaN();  // [kg/(m3 s)]
    // Specific electrical energy (compressor and vacuum-pump work per kg of
    // net product). `energy_consumption` is the original name, kept for
    // compatibility; `electrical_energy` holds the same value.
    double energy_consumption = std::numeric_limits<double>::quiet_NaN();   // [J/kg]
    double production_rate = std::numeric_limits<double>::quiet_NaN();      // [mol/s]
    double electrical_energy = std::numeric_limits<double>::quiet_NaN();    // [J/kg], = energy_consumption
    // Specific thermal energy: heat supplied to the bed over the cycle (step
    // energy sources; a heating medium when there is one) per kg of net
    // product. Defined and undefined as the electrical energy.
    double thermal_energy = std::numeric_limits<double>::quiet_NaN();       // [J/kg]
    // Specific cooling duty: heat removed over the cycle per kg of net product,
    // as a positive magnitude. Reported for information, not part of the
    // thermal energy.
    double cooling_duty = std::numeric_limits<double>::quiet_NaN();         // [J/kg]
    // Each step's share of thermal_energy and cooling_duty, in the deck's step
    // order (they sum to the cycle values); empty when `valid` is false.
    std::vector<double> step_thermal_energy;                                // [J/kg]
    std::vector<double> step_cooling_duty;                                  // [J/kg]

    // The same KPIs from the KPI quadratures (endpoint flows and pump power
    // integrated by IDAS alongside the solution, kpi_quadrature.h) instead of
    // the trapezoid over the logged samples. Reported beside the KPIs above,
    // which stay the reported values; NaN when the quadratures were off.
    double quadrature_purity = std::numeric_limits<double>::quiet_NaN();
    double quadrature_recovery = std::numeric_limits<double>::quiet_NaN();
    double quadrature_productivity_mass = std::numeric_limits<double>::quiet_NaN();
    double quadrature_productivity_volume = std::numeric_limits<double>::quiet_NaN();
    double quadrature_energy_consumption = std::numeric_limits<double>::quiet_NaN();
    double quadrature_production_rate = std::numeric_limits<double>::quiet_NaN();

    // Cycle mole balance of the component [-]: (moles that entered through
    // every endpoint of every step, net, minus the change of the column's gas
    // and adsorbed inventory over the cycle) over the moles that entered,
    // with the endpoint flows from the trapezoid (mass_balance_error) or the
    // quadratures. Zero for exact conservation; NaN without the quadratures
    // (the inventory is recorded with them).
    double mass_balance_error = std::numeric_limits<double>::quiet_NaN();
    double quadrature_mass_balance_error = std::numeric_limits<double>::quiet_NaN();

    // ---- Reaction KPIs (KPIConfig::reaction) ------------------------------
    // Set on the records of the key reactant A and of each product P when
    // the deck has reactions and a key reactant is configured; otherwise
    // `reaction_kpis` is false and every field below is NaN (empty unit).
    // The key reactant and the products get a record even when their
    // separation KPIs are not tracked (`valid` then says whether those are
    // configured). From the moles G each component's reactions made in the
    // bed over the cycle (EquationLog::ReactionGeneration) and the moles of A
    // fed, F_A,in (A's KPI feed endpoints, or, when it has none, every
    // endpoint's inflow of A):
    //   conversion        X_A = -G_A / F_A,in           (A's record)
    //   selectivity       S_P = G_P / (-G_A)            (P's record)
    //   reaction_yield    Y_P = G_P / F_A,in            (P's record)
    //   space_time_yield  G_P M_P / (t_cycle m_cat), or / (t_cycle V_bed) with
    //                     the bed-volume basis           (P's record)
    // Plain molar ratios, with no stoichiometric normalisation.
    bool reaction_kpis = false;
    double reaction_generation = std::numeric_limits<double>::quiet_NaN();  // G [mol per cycle]
    double conversion = std::numeric_limits<double>::quiet_NaN();           // [mol/mol]
    double selectivity = std::numeric_limits<double>::quiet_NaN();          // [mol/mol]
    double reaction_yield = std::numeric_limits<double>::quiet_NaN();       // [mol/mol]
    double space_time_yield = std::numeric_limits<double>::quiet_NaN();     // [kg/(kg s)] or [kg/(m3 s)]
    std::string space_time_yield_unit;  // "kg/(kg s)" or "kg/(m3 s)" (the basis); empty without
    // Highest bed temperature over the cycle (the solid temperature in an
    // LTNE run), with the centre of its cell [m from the bottom] and the time
    // from the start of the cycle [s]. The same on every reaction record of
    // the cycle.
    double peak_temperature = std::numeric_limits<double>::quiet_NaN();           // [K]
    double peak_temperature_position = std::numeric_limits<double>::quiet_NaN();  // [m]
    double peak_temperature_time = std::numeric_limits<double>::quiet_NaN();      // [s]
  };

  // Highest bed temperature of one logged cycle: the fluid temperature (the
  // solid temperature in an LTNE run) over the packed cells and every logged
  // point of every step. NaN fields when the cycle is not logged.
  struct PeakTemperature {
    double temperature = std::numeric_limits<double>::quiet_NaN();  // [K]
    double position = std::numeric_limits<double>::quiet_NaN();     // [m] cell centre, from the bottom
    double time = std::numeric_limits<double>::quiet_NaN();         // [s] from the start of the cycle
  };
  PeakTemperature CyclePeakTemperature(const EquationLog& eq_log, int cycle);

  enum class KPIIntegration;
  // Cycle mole balance of one component at one logged cycle [-]
  // (KPIRecord::mass_balance_error): the net moles that entered through every
  // endpoint of every step, plus the moles the reactions made in the bed
  // (EquationLog::ReactionGeneration; 0 without reactions), minus the change
  // of the column inventory over the cycle, over the moles that entered. NaN when the cycle's step balances
  // (recorded with the KPI quadratures) or a stream log are missing, or when
  // the component's inflow is round-off (below 1e-12 of the column inventory).
  double CycleMassBalanceError(const EquationLog& eq_log, const BoundaryLog& b_log, int cycle,
                               const std::string& component, KPIIntegration integration);

  // KPIs of every tracked component of `config` at one logged cycle, computed
  // from the logs, plus the records of the reaction KPIs' key reactant and
  // products (KPIRecord::reaction_kpis). The endpoint references are re-bound by ID against `steps`
  // first (on a copy), so a step or stream renamed since they were resolved
  // still counts. Empty when no component is tracked.
  std::vector<KPIRecord> ComputeCycleKPIs(const KPIConfig& config, const std::vector<Step>& steps,
                                          const EquationLog& eq_log, const BoundaryLog& b_log,
                                          int cycle);

  // The KPIs of each completed cycle, recorded when the cycle completed.
  class KPIHistory {
  public:
    void Clear() { records_.clear(); }
    // Compute and keep the KPIs of `cycle` (replaces an earlier record of it).
    void Record(const KPIConfig& config, const std::vector<Step>& steps,
                const EquationLog& eq_log, const BoundaryLog& b_log, int cycle);

    bool Has(int cycle) const { return records_.find(cycle) != records_.end(); }
    // The records of one cycle (empty when it was not recorded).
    const std::vector<KPIRecord>& Cycle(int cycle) const;
    // Recorded cycle indices, ascending.
    std::vector<int> Cycles() const;
    const std::map<int, std::vector<KPIRecord>>& Records() const { return records_; }

  private:
    std::map<int, std::vector<KPIRecord>> records_;
  };

}  // namespace adsim
