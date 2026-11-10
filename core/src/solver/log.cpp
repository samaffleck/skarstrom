/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/solver/log.h"
#include "adsim/boundary_conditions/step.h"

#include <sstream>
#include "adsim/solver/cycle.h"
#include "adsim/objects/reactor.h"
#include "adsim/objects/user_data.h"

#include <filesystem>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <charconv>
#include <numbers>

namespace adsim {

  namespace {
    // True when some cell holds more than one solid (a mixed layer).
    bool HasMixedCell(const DomainVectors& domain) {
      for (const auto& cell : domain.cell_adsorbents)
        if (cell.size() > 1) return true;
      return false;
    }
  }  // namespace

  std::string PerAdsorbentVariableName(const std::string& component, ComponentVariableType type,
                                       const std::string& adsorbent_label) {
    return component + " " + ComponentVariableTypeList[static_cast<size_t>(type)] + " (" + adsorbent_label + ")";
  }

  std::string ReactionRateVariableName(const std::string& reaction, const std::string& adsorbent_label) {
    return reaction + " Rate (" + adsorbent_label + ")";
  }

  void EquationLog::Initialise(UserData& user_data) {
    auto& settings = user_data.settings;
    auto& fluid = user_data.fluid;
    auto& bed = user_data.bed;
    auto& cycle = user_data.cycle;
    auto& state = user_data.state;

    Clear();

    domain = state.domain_vectors;
    non_local_ = bed.equations_.thermal_equilibrium_mode_ == ThermalEquilibriumMode::NON_LOCAL;
    per_adsorbent_ = HasMixedCell(domain);

    // Adsorbent labels for the per-adsorbent series: the name, with " #<id>"
    // when another adsorbent has the same name.
    adsorbent_labels_.clear();
    {
      std::map<std::string, int> uses;
      for (const auto& [id, ads] : user_data.adsorbent_directory.GetAll()) ++uses[ads.name_];
      for (const auto& [id, ads] : user_data.adsorbent_directory.GetAll())
        adsorbent_labels_[id] = uses[ads.name_] > 1 ? ads.name_ + " #" + std::to_string(id) : ads.name_;
    }
    // The reactions whose rates are logged (one per active reaction per solid).
    reactions_.clear();
    for (const auto& field : state.reaction_fields) {
      LoggedReaction r;
      r.adsorbent = field->adsorbent;
      r.reaction = field->reaction;
      r.reaction_name = field->reaction_name;
      auto label = adsorbent_labels_.find(field->adsorbent);
      r.adsorbent_label = label == adsorbent_labels_.end() ? std::to_string(field->adsorbent) : label->second;
      reactions_.push_back(std::move(r));
    }
    // Keep at least the cycles the CSS check reads (it is a function of the log).
    // When a cycle is checked the log holds it plus max_cycles_to_keep_ earlier
    // cycles (eviction happens when the next cycle starts), hence the - 1.
    max_cycles_to_keep_ = std::max({1, settings.max_cycles_to_keep_, settings.css.CyclesToLog() - 1});
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

  MediumHeatFlow HeatingMediumHeatFlow(const Equations& state, double time_in_step, bool non_local) {
    MediumHeatFlow flow;
    if (state.IsIsothermal()) return flow;
    const Step* step = state.CurrentStep();
    const HeatingMedium* medium = step ? step->ActiveHeatingMedium() : nullptr;
    if (!medium) return flow;

    // The same outer-face term as the wall residual (UpdateWallTemperature),
    // per unit time and over the true outer area pi d_o dx of each cell, with
    // the coefficient the residual uses: h_m in every bed cell, 0 in a void
    // cell (UpdateExternalHeatTransferCoefficient). Taken from the medium, not
    // from the state's coefficient, which a log point taken before the step's
    // first residual evaluation could still hold from the previous step.
    const double T_m = medium->Temperature(time_in_step);
    const double h_m = medium->heat_transfer_coefficient_;
    const auto& Tw = state.wall_temperature.var;
    const auto& d_o = state.domain_vectors.de;
    const auto& dx = state.domain_vectors.cell_width;
    const auto& is_void = state.domain_vectors.is_void_cell;
    const size_t n = std::min({Tw.size(), d_o.size(), dx.size()});
    for (size_t i = 1; i + 1 < n; ++i) {  // interior cells only
      if (i < is_void.size() && is_void[i]) continue;
      const double q = std::numbers::pi * d_o[i] * dx[i] * h_m * (T_m - Tw[i]);
      if (q > 0.0) flow.in += q;
      else flow.out -= q;
    }

    // The internal exchanger: the same term as the bed balances
    // (UpdateTemperature, and under LTNE UpdateFluidTemperature with the eb
    // share and UpdateSolidTemperature with the 1 - eb share), bed cells only.
    const double ua = medium->ExchangerConductance();
    if (ua > 0.0) {
      const auto& T = state.fluid_temperature.var;
      const auto& Ts = state.solid_temperature.var;
      const auto& dV = state.domain_vectors.volume;
      const auto& eb = state.domain_vectors.eb;
      const size_t m = std::min({T.size(), dV.size(), eb.size(), non_local ? Ts.size() : T.size()});
      for (size_t i = 1; i + 1 < m; ++i) {
        if (i < is_void.size() && is_void[i]) continue;
        const double dT = non_local ? eb[i] * (T_m - T[i]) + (1.0 - eb[i]) * (T_m - Ts[i]) : T_m - T[i];
        const double q = ua * dV[i] * dT;
        if (q > 0.0) { flow.in += q; flow.exchanger_in += q; }
        else { flow.out -= q; flow.exchanger_out -= q; }
      }
    }
    return flow;
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
      per_adsorbent_ = HasMixedCell(domain);
      sys_data.cycle_index_set_.insert(cycle_index);
      sys_data.step_index_set_.insert(step_index);
    }

    // Log regular variables (coupled equations)
    // Use [] operator to create entries if they don't exist
    auto& sim_log = var_log_[cycle_index][step_index];

    // The step's heat duty: the energy-source powers applied over the step
    // (the domain snapshot above), stored once per step with its log.
    if (is_new_step) {
      StepHeatDuty duty;
      if (!state.IsIsothermal()) {
        const auto& q = state.domain_vectors.energy_source_power;  // [W/m3]
        const auto& vol = state.domain_vectors.volume;              // [m3]
        const size_t n = std::min(q.size(), vol.size());
        for (size_t i = 0; i < n; ++i) {
          const double p = q[i] * vol[i];
          if (p > 0.0) duty.heating_power += p;
          else duty.cooling_power -= p;
        }
      }
      sim_log.heat_duty = duty;
    }

    // The heating medium's duty: integrate its heat flow into the wall over
    // the step (trapezoidal rule between consecutive logged points). The first
    // logged point of a step is the step start.
    {
      const double step_start = sim_log.step_data.times.empty() ? t : sim_log.step_data.times.front();
      const MediumHeatFlow flow = HeatingMediumHeatFlow(state, t - step_start, non_local_);
      if (!is_new_step) {
        const double dt = t - medium_prev_t_;
        if (dt > 0.0) {
          sim_log.heat_duty.medium_heat_in += 0.5 * (medium_prev_flow_.in + flow.in) * dt;
          sim_log.heat_duty.medium_heat_out += 0.5 * (medium_prev_flow_.out + flow.out) * dt;
          sim_log.heat_duty.exchanger_heat_in += 0.5 * (medium_prev_flow_.exchanger_in + flow.exchanger_in) * dt;
          sim_log.heat_duty.exchanger_heat_out += 0.5 * (medium_prev_flow_.exchanger_out + flow.exchanger_out) * dt;
        }
      }
      // The reactions' generation, integrated the same way (the previous
      // point's time is medium_prev_t_, read before it is updated).
      if (state.HasReactions()) {
        const auto& dV = state.domain_vectors.volume;
        for (const auto& [component, eq] : state.component_reaction_source) {
          double rate = 0.0;  // [mol/s]
          const size_t n = std::min(eq.var.size(), dV.size());
          for (size_t i = 0; i < n; ++i) rate += dV[i] * eq.var[i];
          double& total = sim_log.reaction_generation[component];
          if (!is_new_step) {
            const double dt = t - medium_prev_t_;
            if (dt > 0.0) total += 0.5 * (reaction_prev_rate_[component] + rate) * dt;
          }
          reaction_prev_rate_[component] = rate;
        }
      }
      medium_prev_t_ = t;
      medium_prev_flow_ = flow;
    }

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

    // Reactions: the heat, each component's source and each reaction's rate
    // on each solid (only for a deck with reactions).
    if (state.HasReactions()) {
      log_var(VariableType::REACTION_HEAT, state.reaction_heat.var);
      log_comp_var(ComponentVariableType::REACTION_SOURCE, state.component_reaction_source);
      for (const auto& field : state.reaction_fields) {
        auto& entry = sim_log.reaction_variables[{field->adsorbent, field->reaction}];
        if (entry.values.empty()) {
          auto label = adsorbent_labels_.find(field->adsorbent);
          entry.variable_name = ReactionRateVariableName(
              field->reaction_name, label == adsorbent_labels_.end() ? std::to_string(field->adsorbent) : label->second);
          entry.units = kReactionRateUnits;
          entry.step_name = step_name;
        }
        entry.values.push_back(CopyRow(field->rate.var));
      }
    }

    // Adsorbed-phase quantities are tracked per adsorbent internally. Log the
    // mass-fraction-weighted aggregate per component (sum_k w_k * x_{k,comp}),
    // which reduces to the single adsorbent's value for a single-adsorbent bed.
    // In a bed with a mixed layer each solid's own series are logged too
    // (below), so that what each one holds can be read and plotted.
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
    if (state.HasSlowSites()) {
      // Slow-site loading of the dual-rate LDF components, aggregated like the
      // loading; an (adsorbent, component) without a slow site adds 0.
      for (const auto& [component, eq] : state.mole_fraction) {
        auto& entry = sim_log.component_variables[component][ComponentVariableType::SLOW_ADSORBED_CONCENTRATION];
        if (entry.values.empty()) {
          entry.variable_name =
              component + " " + ComponentVariableTypeList[static_cast<size_t>(ComponentVariableType::SLOW_ADSORBED_CONCENTRATION)];
          entry.units = GetComponentVariableUnits(ComponentVariableType::SLOW_ADSORBED_CONCENTRATION);
          entry.step_name = step_name;
        }
        const size_t M = state.mole_fraction.at(component).var.size();
        Var1D out = TakeRow();
        out.assign(M, 0.0);
        for (size_t i = 0; i < M; ++i) {
          double s = 0.0;
          for (const auto& ca : cell_ads[i]) {
            auto a = state.ads_concentration_slow.find(ca.id);
            if (a == state.ads_concentration_slow.end()) continue;
            auto c = a->second.find(component);
            if (c == a->second.end()) continue;
            s += ca.mass_fraction * c->second.var[i];
          }
          out[i] = s;
        }
        if (M > 0 && !bed_cells.empty()) {
          out[0] = out[bed_cells.front()];
          out[M - 1] = out[bed_cells.back()];
        }
        entry.values.push_back(std::move(out));
      }
    }

    // Per-adsorbent series: in a bed where some cell holds more than one solid,
    // each solid's own loading [per kg of that solid] beside the aggregate, 0
    // in the cells without it and the ghost cells mirroring the first and last
    // bed cell, as the aggregate. Skipped otherwise: there the aggregate of
    // each cell is the loading of its one solid (GetAdsorbentComponentVariable
    // reads it from there), and the log stays the same size.
    if (per_adsorbent_) {
      // `x` (one adsorbent's values, or nullptr for none) in the cells that
      // hold adsorbent `id`.
      auto masked = [&](AdsorbentID id, const std::vector<double>* x, size_t M) {
        Var1D out = TakeRow();
        out.assign(M, 0.0);
        if (x) {
          for (size_t i = 0; i < M && i < cell_ads.size() && i < x->size(); ++i)
            for (const auto& ca : cell_ads[i])
              if (ca.id == id) { out[i] = (*x)[i]; break; }
        }
        if (M > 0 && !bed_cells.empty()) {
          out[0] = out[bed_cells.front()];
          out[M - 1] = out[bed_cells.back()];
        }
        return out;
      };
      auto log_per_adsorbent = [&](ComponentVariableType type, const auto& nested) {
        for (const auto& [component, y] : state.mole_fraction) {
          const size_t M = y.var.size();
          for (const auto& [id, label] : adsorbent_labels_) {
            auto a = state.ads_concentration.find(id);
            if (a == state.ads_concentration.end()) continue;  // not in the bed
            const std::vector<double>* x = nullptr;
            if (auto n = nested.find(id); n != nested.end())
              if (auto c = n->second.find(component); c != n->second.end()) x = &c->second.var;
            auto& entry = sim_log.adsorbent_variables[id][component][type];
            if (entry.values.empty()) {
              entry.variable_name = PerAdsorbentVariableName(component, type, label);
              entry.units = GetComponentVariableUnits(type);
              entry.step_name = step_name;
            }
            entry.values.push_back(masked(id, x, M));
          }
        }
      };
      log_per_adsorbent(ComponentVariableType::ADSORBED_CONCENTRATION, state.ads_concentration);
      log_per_adsorbent(ComponentVariableType::EQUILIBRIUM_ADSORBED_CONCENTRATION, state.eq_ads_concentration);
      log_per_adsorbent(ComponentVariableType::MASS_TRANSFER_COEFFICIENT, state.mass_transfer_coefficient);
      if (state.HasSlowSites())
        log_per_adsorbent(ComponentVariableType::SLOW_ADSORBED_CONCENTRATION, state.ads_concentration_slow);
    }

    ThinIfOverBudget(sim_log);
  }

  // A row is every variable's snapshot plus its time. A step whose rows would
  // exceed max_bytes_per_step_ (a stalled or very stiff run taking very many
  // solver steps) is thinned: its rows are halved, keeping every other one
  // from the first, and from then on only every `stride`-th LogState call
  // (doubled at each halving) is kept for good, so the rows stay evenly spaced
  // in calls. The newest row is always kept: if it is off the stride it is
  // provisional and the next call replaces it, so the end-of-step state (read
  // by the CSS check and its acceleration) stays exact. The step then holds
  // between half the budget and the budget. Profiles, series and the face-flux
  // integrals of the mass-balance check use these coarser points in that step;
  // the per-stream logs that feed the KPIs are not affected.
  void EquationLog::ThinIfOverBudget(Log& log) {
    auto& sd = log.step_data;
    const std::size_t call = sd.calls++;  // index of the row LogState just appended
    if (max_bytes_per_step_ == 0) return;

    auto for_each_series = [&](auto&& f) {
      for (auto& [type, data] : log.variables) f(data.values);
      for (auto& [component, vars] : log.component_variables)
        for (auto& [type, data] : vars) f(data.values);
      for (auto& [id, comps] : log.adsorbent_variables)
        for (auto& [component, vars] : comps)
          for (auto& [type, data] : vars) f(data.values);
      for (auto& [key, data] : log.reaction_variables) f(data.values);
    };
    if (sd.row_bytes == 0) {
      std::size_t doubles = 1;  // the time
      for_each_series([&](auto& values) { if (!values.empty()) doubles += values.back().size(); });
      sd.row_bytes = doubles * sizeof(double);
    }

    // Drop the previous provisional row (second to last): the new row replaces it.
    auto erase_at = [&](std::size_t i) {
      sd.times.erase(sd.times.begin() + static_cast<std::ptrdiff_t>(i));
      for_each_series([&](auto& values) {
        if (i >= values.size()) return;
        spare_rows_.push_back(std::move(values[i]));
        values.erase(values.begin() + static_cast<std::ptrdiff_t>(i));
      });
    };
    if (sd.last_provisional && sd.times.size() >= 2) erase_at(sd.times.size() - 2);
    sd.last_provisional = sd.stride > 1 && call % sd.stride != 0;

    const std::size_t n = sd.times.size();
    if (n < 3 || n * sd.row_bytes <= max_bytes_per_step_) return;

    // Halve: the rows on the stride are at calls 0, s, 2s, ...; keep every
    // other one (calls 0, 2s, 4s, ...) and the newest row.
    const std::size_t m = sd.last_provisional ? n - 1 : n;  // rows on the stride
    auto keep = [m, n](std::size_t i) { return i == n - 1 || (i < m && i % 2 == 0); };
    std::size_t w = 0;
    for (std::size_t i = 0; i < n; ++i)
      if (keep(i)) sd.times[w++] = sd.times[i];
    sd.times.resize(w);
    for_each_series([&](auto& values) {
      if (values.size() != n) return;  // not logged at every row: leave it
      std::size_t k = 0;
      for (std::size_t i = 0; i < n; ++i) {
        if (keep(i)) {
          if (k != i) values[k] = std::move(values[i]);
          ++k;
        } else {
          spare_rows_.push_back(std::move(values[i]));
        }
      }
      values.resize(k);
    });
    sd.stride *= 2;
    sd.last_provisional = call % sd.stride != 0;

    if (sd.thinnings++ == 0 && message_sink_) {
      const std::string& name = (sd.step_index >= 0 && static_cast<std::size_t>(sd.step_index) < sys_data.step_list_.size())
                                    ? sys_data.step_list_[static_cast<std::size_t>(sd.step_index)]
                                    : std::string("?");
      std::ostringstream ss;
      ss << "  Step \"" << name << "\" (cycle " << sd.cycle_index + 1 << "): the result log reached its memory budget ("
         << (max_bytes_per_step_ >> 20) << " MB, " << n << " solver steps logged); keeping every other logged point"
         << " from here, so this step's profiles and mass-balance integrals use coarser time points.";
      message_sink_(ss.str());
    }
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

        // Per-adsorbent variables (mixed layers), named like
        // "CarbonDioxide_Adsorbed_Concentration__13X_.csv"
        for (const auto& [ads_id, comps] : step_data.adsorbent_variables) {
          for (const auto& [component, comp_vars] : comps) {
            for (const auto& [comp_var_type, entry] : comp_vars) {
              if (step_data.step_data.times.empty() || entry.values.empty()) continue;
              const std::string filename = SanitizeFilename(entry.variable_name) + ".csv";
              if (files.find(filename) == files.end()) {
                files[filename].open(out_dir / filename);
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

        // Reaction rates, named like "SMR_Rate__Ni_catalyst_.csv"
        for (const auto& [key, entry] : step_data.reaction_variables) {
          if (step_data.step_data.times.empty() || entry.values.empty()) continue;
          const std::string filename = SanitizeFilename(entry.variable_name) + ".csv";
          if (files.find(filename) == files.end()) {
            files[filename].open(out_dir / filename);
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

        // Per-adsorbent variables
        for (const auto& [ads_id, comps] : step_data.adsorbent_variables) {
          for (const auto& [component, comp_vars] : comps) {
            for (const auto& [comp_var_type, entry] : comp_vars) {
              if (shared_times.empty() || entry.values.empty()) continue;
              auto& file = files[SanitizeFilename(entry.variable_name) + ".csv"];
              if (!file.is_open()) continue;
              WriteCsvRows(file, cycle_index, step_index, shared_times, entry.values, max_log_points_per_step_);
            }
          }
        }

        // Reaction rates
        for (const auto& [key, entry] : step_data.reaction_variables) {
          if (shared_times.empty() || entry.values.empty()) continue;
          auto& file = files[SanitizeFilename(entry.variable_name) + ".csv"];
          if (!file.is_open()) continue;
          WriteCsvRows(file, cycle_index, step_index, shared_times, entry.values, max_log_points_per_step_);
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
    for (auto& [id, comps] : log.adsorbent_variables)
      for (auto& [component, vars] : comps)
        for (auto& [type, data] : vars) recycle(data);
    for (auto& [key, data] : log.reaction_variables) recycle(data);
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

  namespace {
    // Public entry rebuilt from the compact per-step storage (times are kept
    // once per step, the rows per variable).
    VariableLogEntry MakeEntry(const StepData& step, const VariableData& var) {
      VariableLogEntry entry;
      entry.cycle_index = step.cycle_index;
      entry.step_index = step.step_index;
      entry.times = step.times;
      entry.values = var.values;
      entry.variable_name = var.variable_name;
      entry.units = var.units;
      entry.step_name = var.step_name;
      entry.object_name = var.object_name;
      return entry;
    }

    // Face values of a cell-centred row: face i lies between cells i and i+1,
    // dx[i] / 2 from the centre of cell i and dx[i+1] / 2 from that of cell
    // i+1, so the linear interpolant weights each side by the OTHER cell's
    // width, (y_i dx_{i+1} + y_{i+1} dx_i) / (dx_i + dx_{i+1}), and is exact
    // for a linear profile on a non-uniform mesh. The ghost cells (first and
    // last) mirror their neighbour: the boundary conditions set them so that
    // the boundary face value is the plain mean of the ghost and the end cell
    // (e.g. P_face = (P_0 + P_1) / 2), whatever width the ghost carries (it
    // copies the first/last bed cell, which differs from a void end cell), so
    // those two faces take equal weights.
    std::vector<double> LinearInterpolateCellCenter(const std::vector<double>& data,
      const std::vector<double>& dx)
    {
      assert(data.size() == dx.size());
      const size_t n = data.size();
      if (n < 2) return {};

      std::vector<double> out(n - 1);
      for (size_t i = 0; i + 1 < n; ++i) {
        const bool boundary_face = (i == 0) || (i + 2 == n);
        const double sum = dx[i] + dx[i + 1];
        if (boundary_face || !(sum > 0.0)) {
          out[i] = 0.5 * (data[i] + data[i + 1]);
        } else {
          out[i] = (data[i] * dx[i + 1] + data[i + 1] * dx[i]) / sum;
        }
      }
      return out;
    }

    // The entry with every row interpolated from the cell centres to the cell
    // faces. An entry that is already face-centred (e.g. velocity: one value
    // fewer than there are cells) or empty is returned as it is.
    VariableLogEntry InterpolateToFaces(VariableLogEntry entry, const std::vector<double>& dx) {
      if (entry.values.empty()) return VariableLogEntry{};
      if (entry.values[0].size() == dx.size() - 1) return entry;
      for (auto& row : entry.values) row = LinearInterpolateCellCenter(row, dx);
      return entry;
    }
  }

  VariableLogEntry EquationLog::GetVariable(int cycle_index, int step_index, VariableType type) const {
    const Log& log = GetLog(cycle_index, step_index);
    auto var_it = log.variables.find(type);
    if (var_it == log.variables.end()) return VariableLogEntry{};
    return MakeEntry(log.step_data, var_it->second);
  }

  VariableLogEntry EquationLog::GetReactionRate(int cycle_index, int step_index, AdsorbentID adsorbent,
                                              std::size_t reaction) const {
    const Log& log = GetLog(cycle_index, step_index);
    auto it = log.reaction_variables.find({adsorbent, reaction});
    if (it == log.reaction_variables.end()) return VariableLogEntry{};
    return MakeEntry(log.step_data, it->second);
  }

  double EquationLog::ReactionGeneration(int cycle_index, int step_index, const std::string& component) const {
    const Log& log = GetLog(cycle_index, step_index);
    auto it = log.reaction_generation.find(component);
    return it == log.reaction_generation.end() ? 0.0 : it->second;
  }

  VariableLogEntry EquationLog::GetComponentVariable(int cycle_index, int step_index,
      ComponentVariableType type, const std::string& component) const
  {
    const Log& log = GetLog(cycle_index, step_index);
    auto comp_it = log.component_variables.find(component);
    if (comp_it == log.component_variables.end()) return VariableLogEntry{};
    auto type_it = comp_it->second.find(type);
    if (type_it == comp_it->second.end()) return VariableLogEntry{};
    return MakeEntry(log.step_data, type_it->second);
  }

  VariableLogEntry EquationLog::GetAdsorbentComponentVariable(int cycle_index, int step_index,
      ComponentVariableType type, const std::string& component, AdsorbentID adsorbent) const
  {
    if (!IsPerAdsorbentVariable(type)) return VariableLogEntry{};
    const Log& log = GetLog(cycle_index, step_index);
    if (!log.adsorbent_variables.empty()) {
      auto ads_it = log.adsorbent_variables.find(adsorbent);
      if (ads_it == log.adsorbent_variables.end()) return VariableLogEntry{};
      auto comp_it = ads_it->second.find(component);
      if (comp_it == ads_it->second.end()) return VariableLogEntry{};
      auto type_it = comp_it->second.find(type);
      if (type_it == comp_it->second.end()) return VariableLogEntry{};
      return MakeEntry(log.step_data, type_it->second);
    }

    // No mixed layer: every cell holds at most one solid, so in a cell holding
    // `adsorbent` the aggregate w x is its value x (w = 1 unless the layer's
    // one mass fraction was set otherwise).
    auto label = adsorbent_labels_.find(adsorbent);
    if (label == adsorbent_labels_.end()) return VariableLogEntry{};
    VariableLogEntry e = GetComponentVariable(cycle_index, step_index, type, component);
    if (e.times.empty()) return e;
    const auto& cell_ads = domain.cell_adsorbents;
    const auto& bed_cells = domain.bed_cell_indices;
    for (auto& row : e.values) {
      const size_t M = row.size();
      for (size_t i = 0; i < M; ++i) {
        double x = 0.0;
        if (i < cell_ads.size())
          for (const auto& ca : cell_ads[i])
            if (ca.id == adsorbent && ca.mass_fraction > 0.0) { x = row[i] / ca.mass_fraction; break; }
        row[i] = x;
      }
      if (M > 0 && !bed_cells.empty()) {
        row[0] = row[bed_cells.front()];
        row[M - 1] = row[bed_cells.back()];
      }
    }
    e.variable_name = PerAdsorbentVariableName(component, type, label->second);
    return e;
  }

  VariableLogEntry EquationLog::GetFaceVariable(int cycle_index, int step_index, VariableType type) const {
    return InterpolateToFaces(GetVariable(cycle_index, step_index, type), domain.cell_width);
  }

  VariableLogEntry EquationLog::GetComponentFaceVariable(int cycle_index, int step_index,
      ComponentVariableType type, const std::string& component) const {
    return InterpolateToFaces(GetComponentVariable(cycle_index, step_index, type, component), domain.cell_width);
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
    step_balances_.clear();
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

    // Same floor as the equation log, so both logs hold the same cycles.
    max_cycles_to_keep_ = std::max({1, settings.max_cycles_to_keep_, settings.css.CyclesToLog() - 1});
    
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
      step_balances_.erase(oldest_cycle);
      log_.erase(oldest_cycle_it);
      sys_data.cycle_index_set_.erase(oldest_cycle);
      num_cycles--;
    }
  }

  const BoundaryLogEntry* BoundaryLog::FindEntry(int cycle_idx, int step_idx, const std::string& stream_name) const {
    auto cycle_it = log_.find(cycle_idx);
    if (cycle_it == log_.end()) return nullptr;
    auto step_it = cycle_it->second.find(step_idx);
    if (step_it == cycle_it->second.end()) return nullptr;
    auto stream_it = step_it->second.find(stream_name);
    if (stream_it == step_it->second.end()) return nullptr;
    return &stream_it->second;
  }

  void BoundaryLog::LogStepBalance(int cycle_idx, int step_idx, StepBalance balance) {
    step_balances_[cycle_idx][step_idx] = std::move(balance);
  }

  const StepBalance* BoundaryLog::GetStepBalance(int cycle_idx, int step_idx) const {
    const auto c = step_balances_.find(cycle_idx);
    if (c == step_balances_.end()) return nullptr;
    const auto s = c->second.find(step_idx);
    return s == c->second.end() ? nullptr : &s->second;
  }

  const StreamData* BoundaryLog::GetStreamData(int cycle_idx, int step_idx, const std::string& stream_name) const {
    const BoundaryLogEntry* entry = FindEntry(cycle_idx, step_idx, stream_name);
    return entry ? &entry->stream_data : nullptr;
  }

  BoundaryVariableData BoundaryLog::MakeVariableData(const BoundaryLogEntry& entry, const std::string& stream_name,
      std::vector<double> values, std::string variable_name, std::string units) const {
    BoundaryVariableData out;
    out.variable_timeseries = std::move(values);
    out.time = entry.stream_data.time_series_;
    out.units = std::move(units);
    out.variable_name = std::move(variable_name);
    out.cycle_idx = entry.cycle_index;
    out.step_name = sys_data.step_list_[entry.step_index];
    out.object_name = stream_name;
    return out;
  }

  BoundaryVariableData BoundaryLog::GetVariableData(
    int cycle_idx, 
    int step_idx, 
    const std::string& stream_name, 
    BoundaryVariableType variable
  ) const {
    const BoundaryLogEntry* entry = FindEntry(cycle_idx, step_idx, stream_name);
    if (!entry) return BoundaryVariableData{};
    return MakeVariableData(*entry, stream_name,
      GetBoundaryVariableVector(variable, entry->stream_data),
      std::string(kBoundaryVariableTypeList[static_cast<int>(variable)]), GetVariableUnits(variable));
  }

  BoundaryVariableData BoundaryLog::GetComponentVariableData(
    int cycle_idx, 
    int step_idx, 
    const std::string& stream_name, 
    BoundaryComponentVariableType variable, 
    const std::string& component
  ) const {
    const BoundaryLogEntry* entry = FindEntry(cycle_idx, step_idx, stream_name);
    if (!entry) return BoundaryVariableData{};
    return MakeVariableData(*entry, stream_name,
      GetBoundaryComponentVariableVector(variable, entry->stream_data, component),
      std::string(kBoundaryComponentVariableTypeList[static_cast<int>(variable)]), GetComponentVariableUnits(variable));
  }

  std::vector<double> GetBoundaryVariableVector(BoundaryVariableType variable, const StreamData& stream_data) {
    switch (variable) {
    case BoundaryVariableType::TEMPERATURE: return stream_data.temperature_boundary_;
    case BoundaryVariableType::MOLAR_FLOW_RATE: return stream_data.molar_flow_;
    case BoundaryVariableType::COMPRESSOR_POWER: return stream_data.electrical_power_;
    case BoundaryVariableType::PRESSURE: return stream_data.pressure_boundary_;
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
      case BoundaryVariableType::PRESSURE: return "[Pa]";
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
      case VariableType::REACTION_HEAT: return "[J/(m^3.s)]";
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
      case ComponentVariableType::SLOW_ADSORBED_CONCENTRATION: return "[mol/kg]";
      case ComponentVariableType::REACTION_SOURCE: return "[mol/(m^3.s)]";
      default: return "[-]";
    }
  }

} // end adsim namespace
