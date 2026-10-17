/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Shared helpers for the analytical validation suite.
//
// Everything in here is written from first principles (textbook energy /
// species balances, ideal-gas relations, cylinder geometry) rather than from
// the solver's own discretisation. That is deliberate: these helpers form the
// independent reference the solver is measured against, so they must NOT be
// refactored to share code with core/src.

#include "../fixtures/user_data_test.h"

#include <cmath>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

namespace adsim::validation {

  inline constexpr double kR = constants::kIdealGasConstant;

  // ---------------------------------------------------------------------------
  // Gas property helpers
  // ---------------------------------------------------------------------------

  // Molar isobaric heat capacity [J/(mol-K)] of a mixture, using the same cached
  // mass-basis Cp values the solver caches from the initial conditions.
  // For a single-component gas (the case used by the energy tests) this is
  // exactly Cp_mass * M.
  inline double MixtureMolarCp(const Fluid& fluid,
                               const std::unordered_map<std::string, double>& y) {
    double cp_mass = 0.0;  // [J/(kg-K)]
    double mm = 0.0;       // [kg/mol]
    for (const auto& [component, yi] : y) {
      cp_mass += yi * fluid.heat_capacity_.at(component);
      mm += yi * fluid.molar_mass_.at(component);
    }
    return cp_mass * mm;
  }

  // Molar isochoric heat capacity [J/(mol-K)]. The solver's equation of state is
  // the ideal gas law (C = P/RT), so Cv = Cp - R is the thermodynamically
  // consistent partner of whatever Cp it caches.
  inline double MixtureMolarCv(const Fluid& fluid,
                               const std::unordered_map<std::string, double>& y) {
    return MixtureMolarCp(fluid, y) - kR;
  }

  // ---------------------------------------------------------------------------
  // Log access
  // ---------------------------------------------------------------------------

  // Number of interior (non-ghost) cells in a logged snapshot.
  inline size_t InteriorCellCount(const std::vector<double>& snapshot) {
    return snapshot.size() >= 2 ? snapshot.size() - 2 : 0;
  }

  // Volume-weighted mean of a cell-centred variable over the interior cells.
  inline double InteriorMean(const std::vector<double>& snapshot,
                             const std::vector<double>& volume) {
    double num = 0.0, den = 0.0;
    for (size_t i = 1; i + 1 < snapshot.size(); ++i) {
      num += snapshot[i] * volume[i];
      den += volume[i];
    }
    return den > 0.0 ? num / den : 0.0;
  }

  // Linear interpolation of a per-snapshot scalar series at a requested time.
  inline double InterpolateSeries(const std::vector<double>& times,
                                  const std::vector<double>& values,
                                  double t) {
    return interp_linear(times, values, t);
  }

  // Trapezoidal time integral of a per-snapshot scalar series.
  inline double TimeIntegral(const std::vector<double>& times,
                             const std::vector<double>& values) {
    double integral = 0.0;
    for (size_t i = 1; i < times.size() && i < values.size(); ++i) {
      integral += 0.5 * (values[i] + values[i - 1]) * (times[i] - times[i - 1]);
    }
    return integral;
  }

  // ---------------------------------------------------------------------------
  // Cylinder geometry (independent of DomainVectors' own bookkeeping)
  // ---------------------------------------------------------------------------

  // Wall (annulus) volume of cell i [m3]: pi/4 * (do^2 - di^2) * dx.
  inline double WallCellVolume(const DomainVectors& dom, size_t i) {
    const double di = dom.di[i];
    const double de = dom.de[i];
    return 0.25 * M_PI * (de * de - di * di) * dom.cell_width[i];
  }

  // Internal (bed-side) wetted area of cell i [m2]: pi * di * dx.
  inline double InternalWallArea(const DomainVectors& dom, size_t i) {
    return M_PI * dom.di[i] * dom.cell_width[i];
  }

  // External (ambient-side) area of cell i [m2]: pi * do * dx.
  inline double ExternalWallArea(const DomainVectors& dom, size_t i) {
    return M_PI * dom.de[i] * dom.cell_width[i];
  }

  // ---------------------------------------------------------------------------
  // Energy accounting
  // ---------------------------------------------------------------------------

  // Internal-energy change of the packed bed (gas + solid) between the first and
  // last logged snapshot of a step, for a CLOSED bed [J].
  //
  //   dU = sum_i dV_i * [ et_i * C_i * Cv + rho_b_i * Cps_i ] * (T_end - T_0)
  //
  // The gas term uses Cv (not Cp) because the bed is a rigid, constant-volume
  // control volume: the P-V work term in the solver's temperature equation is
  // exactly what converts Cp into Cv when C is fixed. C_i is evaluated from the
  // initial state and is a constant of the motion for a closed, non-adsorbing
  // bed.
  inline double ClosedBedInternalEnergyChange(const EquationLog& log,
                                              int cycle_index,
                                              int step_index,
                                              double Cv_molar) {
    const auto T = log.GetVariable(cycle_index, step_index, VariableType::TEMPERATURE);
    const auto P = log.GetVariable(cycle_index, step_index, VariableType::PRESSURE);
    if (T.values.empty() || P.values.empty()) return 0.0;

    const auto& dom = log.GetDomain();
    const auto& T0 = T.values.front();
    const auto& Tf = T.values.back();
    const auto& P0 = P.values.front();

    double dU = 0.0;
    for (size_t i = 1; i + 1 < T0.size(); ++i) {
      const double C0 = P0[i] / (kR * T0[i]);
      const double heat_capacity_per_volume =
        dom.et[i] * C0 * Cv_molar + dom.bulk_density[i] * dom.Cps[i];
      dU += dom.volume[i] * heat_capacity_per_volume * (Tf[i] - T0[i]);
    }
    return dU;
  }

  // Internal-energy change of the column wall between the first and last logged
  // snapshot [J]: sum_i V_wall_i * rho_w * Cp_w * (Tw_end - Tw_0).
  inline double WallInternalEnergyChange(const EquationLog& log,
                                         int cycle_index,
                                         int step_index) {
    const auto Tw = log.GetVariable(cycle_index, step_index, VariableType::WALL_TEMPERATURE);
    if (Tw.values.empty()) return 0.0;

    const auto& dom = log.GetDomain();
    const auto& Tw0 = Tw.values.front();
    const auto& Twf = Tw.values.back();

    double dU = 0.0;
    for (size_t i = 1; i + 1 < Tw0.size(); ++i) {
      dU += WallCellVolume(dom, i) * dom.rho_w[i] * dom.Cpw[i] * (Twf[i] - Tw0[i]);
    }
    return dU;
  }

  // Heat lost from the wall to the ambient over a step [J], integrating the
  // textbook Newton-cooling flux over the true outer surface:
  //
  //   Q = int sum_i ( pi * do * dx_i ) * ho_i * (Tw_i - T_amb) dt
  inline double AmbientHeatLoss(const EquationLog& log,
                                int cycle_index,
                                int step_index,
                                double T_ambient) {
    const auto Tw = log.GetVariable(cycle_index, step_index, VariableType::WALL_TEMPERATURE);
    const auto ho = log.GetVariable(cycle_index, step_index,
                                    VariableType::EXTERNAL_HEAT_TRANSFER_COEFFICIENT);
    if (Tw.values.empty() || ho.values.empty()) return 0.0;

    const auto& dom = log.GetDomain();
    std::vector<double> power(Tw.values.size(), 0.0);
    for (size_t n = 0; n < Tw.values.size() && n < ho.values.size(); ++n) {
      double q = 0.0;
      for (size_t i = 1; i + 1 < Tw.values[n].size(); ++i) {
        q += ExternalWallArea(dom, i) * ho.values[n][i] * (Tw.values[n][i] - T_ambient);
      }
      power[n] = q;
    }
    return TimeIntegral(Tw.times, power);
  }

  // Heat transferred from the bed to the wall over a step [J], integrating the
  // textbook flux over the true inner surface:
  //
  //   Q = int sum_i ( pi * di * dx_i ) * hi_i * (T_i - Tw_i) dt
  inline double BedToWallHeatTransfer(const EquationLog& log,
                                      int cycle_index,
                                      int step_index) {
    const auto T = log.GetVariable(cycle_index, step_index, VariableType::TEMPERATURE);
    const auto Tw = log.GetVariable(cycle_index, step_index, VariableType::WALL_TEMPERATURE);
    const auto hi = log.GetVariable(cycle_index, step_index,
                                    VariableType::INTERNAL_HEAT_TRANSFER_COEFFICIENT);
    if (T.values.empty() || Tw.values.empty() || hi.values.empty()) return 0.0;

    const auto& dom = log.GetDomain();
    std::vector<double> power(T.values.size(), 0.0);
    for (size_t n = 0; n < T.values.size() && n < Tw.values.size() && n < hi.values.size(); ++n) {
      double q = 0.0;
      for (size_t i = 1; i + 1 < T.values[n].size(); ++i) {
        q += InternalWallArea(dom, i) * hi.values[n][i] * (T.values[n][i] - Tw.values[n][i]);
      }
      power[n] = q;
    }
    return TimeIntegral(T.times, power);
  }

  // Heat added to the bed by the volumetric energy source over a step [J],
  // integrating the logged ENERGY_SOURCE field (W/m3) over the interior volume
  // and over time. Reading it back from the log rather than from the input makes
  // the energy-balance check fully self-contained.
  inline double EnergySourceHeat(const EquationLog& log,
                                 int cycle_index,
                                 int step_index) {
    const auto Qe = log.GetVariable(cycle_index, step_index, VariableType::ENERGY_SOURCE);
    if (Qe.values.empty()) return 0.0;

    const auto& dom = log.GetDomain();
    std::vector<double> power(Qe.values.size(), 0.0);
    for (size_t n = 0; n < Qe.values.size(); ++n) {
      double p = 0.0;
      for (size_t i = 1; i + 1 < Qe.values[n].size(); ++i) {
        p += Qe.values[n][i] * dom.volume[i];  // [W/m3] * [m3] = [W]
      }
      power[n] = p;
    }
    return TimeIntegral(Qe.times, power);
  }

  // Live per-step energy-closure check for a SEALED bed (both ends closed), the
  // energy analogue of the mole balance in user_data_test.h. With no flow the
  // first law over bed + wall reads, with no approximation,
  //
  //     dU_bed + dU_wall + Q_ambient - Q_source = 0
  //
  // over every logged step: whatever the heater adds, minus what leaves to
  // ambient, equals the internal-energy rise of the bed and wall. Unlike the mass
  // balance this is exact (not discretisation-limited) precisely because the bed
  // is sealed and carries no boundary enthalpy. Returns the largest relative
  // imbalance over all steps; the caller asserts it against a small tolerance.
  //
  // `Cv_molar` is the gas isochoric molar heat capacity (Cv, not Cp, because the
  // rigid sealed vessel converts the Cp in the temperature equation to Cv), and
  // `T_ambient` the ambient temperature the external-flux reference is built from.
  inline double SealedBedEnergyImbalance(const EquationLog& log,
                                         double Cv_molar,
                                         double T_ambient) {
    double worst = 0.0;
    for (int cycle_index : log.sys_data.cycle_index_set_) {
      for (int step_index : log.sys_data.step_index_set_) {
        const double dU_bed = ClosedBedInternalEnergyChange(log, cycle_index, step_index, Cv_molar);
        const double dU_wall = WallInternalEnergyChange(log, cycle_index, step_index);
        const double Q_ambient = AmbientHeatLoss(log, cycle_index, step_index, T_ambient);
        const double Q_source = EnergySourceHeat(log, cycle_index, step_index);

        const double imbalance = dU_bed + dU_wall + Q_ambient - Q_source;
        const double scale = std::max({std::abs(dU_bed), std::abs(dU_wall),
                                       std::abs(Q_ambient), std::abs(Q_source), 1.0});
        worst = std::max(worst, std::abs(imbalance) / scale);
      }
    }
    return worst;
  }

  // Total moles of `component` held in the adsorbed phase at a snapshot [mol].
  inline double AdsorbedMoles(const EquationLog& log,
                              int cycle_index,
                              int step_index,
                              const std::string& component,
                              bool at_end) {
    const auto q = log.GetComponentVariable(cycle_index, step_index,
                                            ComponentVariableType::ADSORBED_CONCENTRATION, component);
    if (q.values.empty()) return 0.0;

    const auto& dom = log.GetDomain();
    const auto& snapshot = at_end ? q.values.back() : q.values.front();

    double moles = 0.0;
    for (size_t i = 1; i + 1 < snapshot.size(); ++i) {
      moles += snapshot[i] * dom.volume[i] * dom.bulk_density[i];
    }
    return moles;
  }

  // Total moles of `component` held in the gas phase at a snapshot [mol].
  inline double GasMoles(const EquationLog& log,
                         int cycle_index,
                         int step_index,
                         const std::string& component,
                         bool at_end) {
    const auto c = log.GetComponentVariable(cycle_index, step_index,
                                            ComponentVariableType::COMPONENT_MOLAR_CONCENTRATION,
                                            component);
    if (c.values.empty()) return 0.0;

    const auto& dom = log.GetDomain();
    const auto& snapshot = at_end ? c.values.back() : c.values.front();

    double moles = 0.0;
    for (size_t i = 1; i + 1 < snapshot.size(); ++i) {
      moles += snapshot[i] * dom.volume[i] * dom.et[i];
    }
    return moles;
  }

  // ---------------------------------------------------------------------------
  // Outlet time series
  // ---------------------------------------------------------------------------

  // Mole fraction of `component` at the top cell face, per logged snapshot.
  inline std::pair<std::vector<double>, std::vector<double>> TopFaceMoleFraction(
      const EquationLog& log, int cycle_index, int step_index, const std::string& component) {
    const auto y = log.GetComponentFaceVariable(cycle_index, step_index,
                                                ComponentVariableType::MOLE_FRACTION, component);
    std::vector<double> series;
    series.reserve(y.values.size());
    for (const auto& snapshot : y.values) {
      series.push_back(snapshot.empty() ? 0.0 : snapshot.back());
    }
    return {y.times, series};
  }

} // namespace adsim::validation
