// kpis.h
#pragma once

// STL includes
#include <cmath>
#include <string>
#include <vector>

// adsim includes
#include "adsim/solver/log.h"
#include "adsim/solver/kpi_config.h"
#include "adsim/math/general_functions.h"


namespace adsim {

  enum class KPIVariableType {
    PURITY,
    RECOVERY,
    PRODUCTIVITY_MASS,
    PRODUCTIVITY_VOLUME,
    ENERGY_CONSUMPTION,
    PRODUCTION_RATE
  };

  inline constexpr std::string_view kKPIVariableTypeList[] = {
    "Purity",
    "Recovery",
    "Productivity (mass)",
    "Productivity (volume)",
    "Energy Consumption",
    "Production Rate"
  };

  static std::string GetKPIUnits(KPIVariableType variable) {
    switch (variable) {
      case KPIVariableType::PURITY: return "[mol/mol]";
      case KPIVariableType::RECOVERY: return "[mol/mol]";
      case KPIVariableType::PRODUCTIVITY_MASS: return "[kg_gas/(kg_solid·s)]";
      case KPIVariableType::PRODUCTIVITY_VOLUME: return "[kg_gas/(m3·s)]";
      case KPIVariableType::ENERGY_CONSUMPTION: return "[J/kg]";
      case KPIVariableType::PRODUCTION_RATE: return "[mol/s]";
      default: return "[-]";
    }
  }

  class KPIParams {
  public:
    KPIParams(const EquationLog& eq_log, const BoundaryLog& b_log) : eq_log_(eq_log), b_log_(b_log) {}

    bool Check() const {
      if (product_comp_.empty()) {
        return false;
      }

      // Need at least one product endpoint and one feed endpoint to attribute
      // throughput. Recycle is no longer a role.
      if (product_endpoints.empty() || feed_endpoints.empty()) {
        return false;
      }

      const auto& mr = eq_log_.sys_data.molar_mass;
      if (mr.find(product_comp_) == mr.end()) {
        return false;
      }

      return true;
    }

    int cycle_index_{};
    std::string product_comp_{};
    // Endpoints (a flow stream or a pressure boundary on a given step) whose
    // throughput counts as feed / product for this component. Each is summed.
    std::vector<KPIEndpointRef> feed_endpoints{};
    std::vector<KPIEndpointRef> product_endpoints{};
    const EquationLog& eq_log_;
    const BoundaryLog& b_log_;
  };

  // Face-reconstruction of total boundary molar flux (positive = entering the
  // bed). Used by cyclic-steady-state detection, the mass-balance plotter
  // panels and the validation tests. NOT used by the endpoint-keyed feed/product
  // KPIs, which sample per-stream logs.
  static double GetTotalMolesIn(
    const adsim::EquationLog& logger,
    int cycle_index,
    int step_index,
    adsim::BoundaryLocation location
  ) {
    const auto& pressure = logger.GetFaceVariable(cycle_index, step_index, adsim::VariableType::PRESSURE);
    const std::vector<adsim::vec>& P_ = pressure.values;

    const auto& temperature = logger.GetFaceVariable(cycle_index, step_index, adsim::VariableType::TEMPERATURE);
    const std::vector<adsim::vec>& T_ = temperature.values;

    const adsim::vec& time = pressure.times;

    const auto& A = logger.GetDomain().area;
    const auto& velocity = logger.GetFaceVariable(cycle_index, step_index, adsim::VariableType::VELOCITY).values;

    adsim::vec ct(P_.size());
    adsim::vec u(P_.size());

    const double R = adsim::constants::kIdealGasConstant;

    double area = A[0];
    const int Na = A.size() - 1;
    if (location == adsim::BoundaryLocation::TOP) area = A[Na];

    for (int i = 0; i < P_.size(); ++i) {
      const auto& Pvec = P_[i];
      const auto& Tvec = T_[i];
      const auto& uvec = velocity[i];

      if (location == adsim::BoundaryLocation::BOTTOM) {
        ct[i] = Pvec[0] / (R * Tvec[0]); // Cant just average C[0] and C[1]
        u[i] = uvec[0];
      }
      else {
        int N = Pvec.size() - 1;
        ct[i] = Pvec[N] / (R * Tvec[N]);
        u[i] = -uvec[uvec.size() - 1];
      }
    }

    // Compute the time integral of u*C using trapezoidal rule
    double integral = 0.0;
    for (int i = 1; i < time.size(); ++i) {
      double dt = time[i] - time[i - 1];
      double flux_i = u[i] * ct[i];
      double flux_prev = u[i - 1] * ct[i - 1];
      integral += 0.5 * (flux_i + flux_prev) * dt;
    }

    return integral * area;
  }

  // Signed moles entering the bed at one endpoint over the step: positive =
  // entering, negative = leaving. The per-stream log records molar flow with a
  // velocity sign (positive = bottom->top), so a TOP endpoint is negated to
  // share the "entering the bed" convention with the rest of the KPI maths.
  static double EndpointMolesIn(
    const adsim::BoundaryLog& b_log,
    int cycle_index,
    const KPIEndpointRef& ep
  ) {
    adsim::BoundaryVariableData d = b_log.GetVariableData(
      cycle_index, ep.step_index, ep.endpoint_name,
      adsim::BoundaryVariableType::MOLAR_FLOW_RATE);
    const double integral = IntegrateVariableInTime(d.variable_timeseries, d.time);
    return (ep.location == adsim::BoundaryLocation::BOTTOM) ? integral : -integral;
  }

  // As EndpointMolesIn, but for a single component (integrates the endpoint's
  // directly-logged per-component molar flow).
  static double EndpointComponentMolesIn(
    const adsim::BoundaryLog& b_log,
    int cycle_index,
    const KPIEndpointRef& ep,
    const std::string& component
  ) {
    adsim::BoundaryVariableData d = b_log.GetComponentVariableData(
      cycle_index, ep.step_index, ep.endpoint_name,
      adsim::BoundaryComponentVariableType::MOLAR_FLOW, component);
    const double integral = IntegrateVariableInTime(d.variable_timeseries, d.time);
    return (ep.location == adsim::BoundaryLocation::BOTTOM) ? integral : -integral;
  }

  // Face-reconstruction of per-component boundary flux (positive = entering the
  // bed), matching the discrete species-equation flux (upwind convection + axial
  // dispersion). Used by cyclic-steady-state detection, the mass-balance plotter
  // panels and the validation tests. NOT used by the endpoint-keyed feed/product
  // KPIs, which sample per-stream logs instead.
  static double GetTotalComponentMolesIn(
    const adsim::EquationLog& logger,
    int cycle_index,
    int step_index,
    adsim::BoundaryLocation location,
    const std::string& component
  ) {
    // Use the same discrete boundary flux form as the species equation:
    // convective upwind flux + axial-dispersion flux at the boundary face.
    const auto& pressure = logger.GetVariable(cycle_index, step_index, adsim::VariableType::PRESSURE);
    const std::vector<adsim::vec>& P_ = pressure.values;

    const auto& temperature = logger.GetVariable(cycle_index, step_index, adsim::VariableType::TEMPERATURE);
    const std::vector<adsim::vec>& T_ = temperature.values;

    const auto& mole_fraction = logger.GetComponentVariable(cycle_index, step_index, adsim::ComponentVariableType::MOLE_FRACTION, component);
    const std::vector<adsim::vec>& yi_ = mole_fraction.values;

    const auto& axial_dispersion = logger.GetVariable(cycle_index, step_index, adsim::VariableType::AXIAL_DISPERSION_COEFFICIENT);
    const std::vector<adsim::vec>& Dax_ = axial_dispersion.values;

    const adsim::vec& time = pressure.times;

    const auto& domain = logger.GetDomain();
    const auto& A = domain.area;
    const auto& dx = domain.cell_width;
    const auto& eb = domain.eb;
    const auto& velocity = logger.GetFaceVariable(cycle_index, step_index, adsim::VariableType::VELOCITY).values;

    if (P_.empty()) return 0.0;

    // Guard against mismatched sizes (can occur if data is still being written).
    const size_t n_snapshots = P_.size();
    if (T_.size() != n_snapshots || yi_.size() != n_snapshots ||
        Dax_.size() != n_snapshots || velocity.size() != n_snapshots) {
      return 0.0;
    }

    adsim::vec flux_in(n_snapshots, 0.0);

    const double R = adsim::constants::kIdealGasConstant;

    double area = A[0];
    const int Na = A.size() - 1;
    if (location == adsim::BoundaryLocation::TOP) area = A[Na];

    for (size_t i = 0; i < n_snapshots; ++i) {
      const auto& Pvec = P_[i];
      const auto& Tvec = T_[i];
      const auto& uvec = velocity[i];
      const auto& yivec = yi_[i];
      const auto& Daxvec = Dax_[i];

      if (location == adsim::BoundaryLocation::BOTTOM) {
        const double u_face = uvec[0];
        const double yi_face = (u_face > 0.0) ? yivec[0] : yivec[1];
        const double Pw = 0.5 * (Pvec[0] + Pvec[1]);
        const double Tw = 0.5 * (Tvec[0] + Tvec[1]);
        const double convective_in = u_face * Pw * yi_face / (R * Tw);

        // D_w term at n=1 (first interior cell), matching UpdateMoleFraction.
        const double fw = dx[1] / (dx[0] + dx[1]);
        double d_w = 0.0;
        if (Daxvec[0] > 0.0 && Daxvec[1] > 0.0) {
          d_w = 1.0 / ((1.0 - fw) / Daxvec[0] + fw / Daxvec[1]);
        }
        const double D_w_ci = (eb[1] * d_w / (R * dx[1])) * (Pw / Tw * (yivec[1] - yivec[0]));

        // Bottom inflow is F_w - D_w.
        flux_in[i] = convective_in - D_w_ci;
      }
      else {
        const int Nf = static_cast<int>(uvec.size()) - 1;
        const int Np = static_cast<int>(Pvec.size()) - 1;      // top ghost cell index
        const int n = Np - 1;                                   // top interior cell index
        const double u_face = uvec[Nf];
        const double yi_face = (u_face > 0.0) ? yivec[n] : yivec[n + 1];
        const double Pe = 0.5 * (Pvec[n] + Pvec[n + 1]);
        const double Te = 0.5 * (Tvec[n] + Tvec[n + 1]);
        const double convective_in = -u_face * Pe * yi_face / (R * Te);

        // D_e term at n=N-2 (last interior cell), matching UpdateMoleFraction.
        const double fe = dx[n + 1] / (dx[n] + dx[n + 1]);
        double d_e = 0.0;
        if (Daxvec[n] > 0.0 && Daxvec[n + 1] > 0.0) {
          d_e = 1.0 / ((1.0 - fe) / Daxvec[n] + fe / Daxvec[n + 1]);
        }
        const double D_e_ci = (eb[n] * d_e / (R * dx[n])) * (Pe / Te * (yivec[n + 1] - yivec[n]));

        // Top inflow is -F_e + D_e.
        flux_in[i] = convective_in + D_e_ci;
      }
    }

    // Time-integrate net boundary flux into the domain.
    double integral = 0.0;
    for (int i = 1; i < time.size(); ++i) {
      double dt = time[i] - time[i - 1];
      double flux_i = flux_in[i];
      double flux_prev = flux_in[i - 1];
      integral += 0.5 * (flux_i + flux_prev) * dt;
    }

    return integral * area;
  }

  static double CalculateAccumulation_Volume(
    const adsim::EquationLog& logger,
    int cycle_index,
    int step_index,
    adsim::VariableType variable_type
  ) {
    const auto& variable = logger.GetVariable(cycle_index, step_index, variable_type);
    if (variable.values.empty()) {
      return 0.0;
    }

    const auto& domain = logger.GetDomain();
    const auto& volume = domain.volume;
    const auto& et = domain.et;

    // Get initial and final values
    const auto& values_t0 = variable.values[0];
    const auto& values_tend = variable.values.back();

    // Integrate over volume: sum(variable * volume * void_fraction) for each cell
    double accumulation = 0.0;
    for (size_t i = 1; i < values_t0.size() - 1; ++i) { // Skip ghost cells
      double value_start = values_t0[i];
      double value_end = values_tend[i];
      accumulation += (value_end - value_start) * volume[i] * et[i];
    }

    return accumulation;
  }

  static double CalculateComponentAccumulation_Volume(
    const adsim::EquationLog& logger,
    int cycle_index,
    int step_index,
    adsim::ComponentVariableType variable_type,
    const std::string& component
  ) {
    const auto& variable = logger.GetComponentVariable(cycle_index, step_index, variable_type, component);
    if (variable.values.empty()) {
      return 0.0;
    }

    const auto& domain = logger.GetDomain();
    const auto& volume = domain.volume;
    const auto& et = domain.et;

    const auto& values_t0 = variable.values[0];
    const auto& values_tend = variable.values.back();

    double accumulation = 0.0;
    for (size_t i = 1; i < values_t0.size() - 1; ++i) { // Skip ghost cells
      double value_start = values_t0[i];
      double value_end = values_tend[i];
      accumulation += (value_end - value_start) * volume[i] * et[i];
    }

    return accumulation;
  }

  static double CalculateComponentAccumulation_Mass(
    const adsim::EquationLog& logger,
    int cycle_index,
    int step_index,
    adsim::ComponentVariableType variable_type,
    const std::string& component
  ) {
    const auto& variable = logger.GetComponentVariable(cycle_index, step_index, variable_type, component);
    if (variable.values.empty()) {
      return 0.0;
    }

    const auto& domain = logger.GetDomain();
    const auto& volume = domain.volume;
    const auto& bulk_density = domain.bulk_density;

    const auto& values_t0 = variable.values[0];
    const auto& values_tend = variable.values.back();

    double accumulation = 0.0;
    for (size_t i = 1; i < values_t0.size() - 1; ++i) { // Skip ghost cells
      double value_start = values_t0[i];
      double value_end = values_tend[i];
      accumulation += (value_end - value_start) * volume[i] * bulk_density[i];
    }

    return accumulation;
  }

  // Net moles leaving the bed across the selected product endpoints. Signed sum:
  // an outlet (net leaving) adds; a product-quality stream re-entering the bed
  // (e.g. light-product pressurisation) subtracts. This is how a "recycle" is
  // accounted — there is no separate recycle role. For a pure outlet the signed
  // sum equals the gross throughput, so ordinary product endpoints are unaffected.
  static double CalculateTotalMolesInProduct(const KPIParams& params) {
    if (!params.Check()) {
      return 0.0;
    }

    double moles_in_product = 0.0;
    for (const auto& ep : params.product_endpoints) {
      double moles_in = EndpointMolesIn(params.b_log_, params.cycle_index_, ep);
      moles_in_product += -moles_in;  // leaving = positive product
    }

    return moles_in_product;
  }

  static double CalculateComponentMolesInProduct(const KPIParams& params) {
    if (!params.Check()) {
      return 0.0;
    }

    double product_moles = 0.0;
    for (const auto& ep : params.product_endpoints) {
      double moles_in = EndpointComponentMolesIn(
        params.b_log_, params.cycle_index_, ep, params.product_comp_);
      product_moles += -moles_in;  // component leaving = positive product
    }

    return product_moles;
  }

  static double CalculateComponentMolesInFeed(const KPIParams& params) {
    if (!params.Check()) {
      return 0.0;
    }

    // Net component moles entering the bed across the selected feed endpoints.
    double component_moles_in_feed = 0.0;
    for (const auto& ep : params.feed_endpoints) {
      double moles_in = EndpointComponentMolesIn(
        params.b_log_, params.cycle_index_, ep, params.product_comp_);
      component_moles_in_feed += moles_in;  // entering = positive feed
    }

    return component_moles_in_feed;
  }

  // A component amount indistinguishable from zero at double precision relative
  // to the total it is compared with (e.g. a trace component whose integrated
  // flow is a few ulp negative) is treated as exactly zero, so a ratio KPI never
  // comes out as a meaningless -1e-22.
  inline double ZeroIfRoundoff(double amount, double reference) {
    return std::abs(amount) <= 1e-12 * std::abs(reference) ? 0.0 : amount;
  }

  // -- Undefined-KPI convention --------------------------------------------
  // Every ratio KPI below (purity, recovery, productivity, production rate,
  // specific energy) is UNDEFINED when its denominator is not a positive,
  // finite number — e.g. no net product drawn in the selected cycle (purity,
  // specific energy), no net feed (recovery), or a zero cycle time / adsorbent
  // mass / bed volume (productivities, production rate) — or when
  // KPIParams::Check() fails. An undefined KPI returns 0.0, never NaN/Inf, and,
  // when the optional `defined` out-parameter is supplied, sets *defined to
  // false (true otherwise).
  //
  // Display / report callers can ignore `defined`. Callers that must not take
  // the placeholder 0 at face value MUST check it: in particular the optimiser
  // cost, where an undefined "Energy Consumption" of 0 would look ideal under
  // "minimise" (OptimizationCostFunction::ComputeCost treats a non-finite KPI
  // value as undefined and penalises it, so feed it NaN when !defined), and the
  // web KPI emitter, which reports an undefined KPI as null.

  // numerator / denominator, or 0 with ok=false when the ratio is undefined (or
  // `ok` is already false, so a chain of ratios propagates "undefined").
  inline double KPIRatio(double numerator, double denominator, bool& ok) {
    if (!ok || !std::isfinite(numerator) || !std::isfinite(denominator) || !(denominator > 0.0)) {
      ok = false;
      return 0.0;
    }
    const double ratio = numerator / denominator;
    if (!std::isfinite(ratio)) {  // overflow for a denormal denominator
      ok = false;
      return 0.0;
    }
    return ratio;
  }

  // Writes the definedness flag (if requested) and returns `value`, or 0 when undefined.
  inline double KPIResult(double value, bool ok, bool* defined) {
    if (defined) *defined = ok;
    return ok ? value : 0.0;
  }

  static double CalculatePurity(const KPIParams& params, bool* defined = nullptr) {

    // -- Validation -- 

    if (!params.Check()) {
      return KPIResult(0.0, false, defined);
    }

    // -- Start calculation --

    bool ok = true;
    double total_moles_in_product = CalculateTotalMolesInProduct(params);
    double component_moles_in_product =
        ZeroIfRoundoff(CalculateComponentMolesInProduct(params), total_moles_in_product);
    double purity = KPIRatio(component_moles_in_product, total_moles_in_product, ok);

    return KPIResult(purity, ok, defined);
  }

  static double CalculateRecovery(const KPIParams& params, bool* defined = nullptr) {
    
    // -- Validation --
    
    if (!params.Check()) {
      return KPIResult(0.0, false, defined);
    }

    // -- Start calculation --

    bool ok = true;
    double moles_in_feed = CalculateComponentMolesInFeed(params);
    double moles_in_product = ZeroIfRoundoff(CalculateComponentMolesInProduct(params), moles_in_feed);
    double recovery = KPIRatio(moles_in_product, moles_in_feed, ok);

    return KPIResult(recovery, ok, defined);
  }

  // Returns [kg product/s]
  static double CalculateProductivity(const KPIParams& params, bool* defined = nullptr) {
    
    // -- Validation --
    
    if (!params.Check()) {
      return KPIResult(0.0, false, defined);
    }
    
    const auto& product_component = params.product_comp_;
    const auto& eq_log = params.eq_log_;
    const auto& mr = eq_log.sys_data.molar_mass;
    
    // -- Start calculation --

    bool ok = true;
    double product_moles = CalculateComponentMolesInProduct(params);    
    double cycle_time = eq_log.sys_data.total_cycle_time; // [s]
    double product_mass = product_moles * mr.at(product_component); // [kg/s]
    double productivity = KPIRatio(product_mass, cycle_time, ok); // [kg product/(s)]

    return KPIResult(productivity, ok, defined);
  }

  static double CalculateProductivityMass(const KPIParams& params, bool* defined = nullptr) {
    const auto& eq_log = params.eq_log_;
    double mass = eq_log.sys_data.adsorbent_mass; // [kg]
    bool ok = true;
    double productivity = CalculateProductivity(params, &ok);
    double mass_productivity = KPIRatio(productivity, mass, ok); // [kg product/(kg-solid.s)]

    return KPIResult(mass_productivity, ok, defined);
  }

  static double CalculateProductivityVolume(const KPIParams& params, bool* defined = nullptr) {
    const auto& eq_log = params.eq_log_;
    auto& cell_volumes = eq_log.GetDomain().volume;

    double volume = 0.0;
    for (double dv : cell_volumes) {
      volume += dv;
    }

    bool ok = true;
    double productivity = CalculateProductivity(params, &ok);
    double volume_productivity = KPIRatio(productivity, volume, ok); // [kg product/(m3.s)]

    return KPIResult(volume_productivity, ok, defined);
  }

  static double CalculateProductionRate(const KPIParams& params, bool* defined = nullptr) {
    // -- Validation --
    if (!params.Check()) {
      return KPIResult(0.0, false, defined);
    }

    const auto& eq_log = params.eq_log_;

    // -- Start calculation --
    bool ok = true;
    const double moles_out = CalculateComponentMolesInProduct(params);
    const double cycle_time = eq_log.sys_data.total_cycle_time; // [s]
    const double production_rate = KPIRatio(moles_out, cycle_time, ok); // [mol/s]

    return KPIResult(production_rate, ok, defined);
  }

  static double CalculateStepEnergy(const KPIParams& params, int step_idx) {
    
    // -- Validation --
    
    if (!params.Check()) {
      return 0.0;
    }
    
    int cycle_index = params.cycle_index_;
    const auto& product_component = params.product_comp_;
    const auto& eq_log = params.eq_log_;
    const auto& b_log = params.b_log_;
    const auto& sys_data = eq_log.sys_data;
    const auto& step_index_set = sys_data.step_index_set_;
    const auto& step_list = sys_data.step_list_;

    const std::string& step_name = step_list[step_idx];

    auto it = b_log.stream_names.find(step_name);
    if (it == b_log.stream_names.end()) {
      return 0.0;
    }

    auto& stream_name_list = b_log.stream_names.at(step_name);

    double step_energy = 0.0;

    for (const auto& stream_name : stream_name_list) {
      adsim::BoundaryVariableData data = b_log.GetVariableData(
        cycle_index, 
        step_idx,
        stream_name, 
        adsim::BoundaryVariableType::COMPRESSOR_POWER
      );

      double stream_energy = IntegrateVariableInTime(data.variable_timeseries, data.time);
      step_energy += stream_energy;
    }

    return step_energy; // [J] 
  }

  static double CalculateSpecificEnergyConsumption(const KPIParams& params, bool* defined = nullptr) {
    
    // -- Validation --
    
    if (!params.Check()) {
      return KPIResult(0.0, false, defined);
    }
    
    int cycle_index = params.cycle_index_;
    const auto& product_component = params.product_comp_;
    const auto& eq_log = params.eq_log_;
    const auto& b_log = params.b_log_;
    const auto& sys_data = eq_log.sys_data;
    const auto& step_index_set = sys_data.step_index_set_;
    const auto& step_list = sys_data.step_list_;
    const auto& mr = sys_data.molar_mass;
    
    double total_energy = 0.0; // [J]
    
    size_t n_steps = step_index_set.size();

    for (size_t s = 0; s < n_steps; ++s) {
      total_energy += CalculateStepEnergy(params, s);
    }

    bool ok = true;
    double product_moles = CalculateComponentMolesInProduct(params);    
    double product_mass = product_moles * mr.at(product_component); // [kg_gas]
    double total_specific_energy = KPIRatio(total_energy, product_mass, ok); // [J/kg_gas]
    
    return KPIResult(total_specific_energy, ok, defined);
  }

}
