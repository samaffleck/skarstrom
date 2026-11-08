// kpis.h
#pragma once

// STL includes
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

// adsim includes
#include "adsim/solver/log.h"
#include "adsim/solver/kpi_config.h"
#include "adsim/math/general_functions.h"


namespace adsim {

  // ENERGY_CONSUMPTION is the specific electrical energy: compressor and
  // vacuum-pump work over the cycle per kg of net product. ELECTRICAL_ENERGY
  // is the same enumerator under its explicit name; the old name is kept so
  // existing code and stored indices keep their meaning. New KPIs (such as a
  // thermal energy) are appended at the end so stored indices do not move.
  enum class KPIVariableType {
    PURITY,
    RECOVERY,
    PRODUCTIVITY_MASS,
    PRODUCTIVITY_VOLUME,
    ENERGY_CONSUMPTION,
    PRODUCTION_RATE,
    THERMAL_ENERGY,
    COOLING_DUTY,
    // Aliases last: an enumerator after them would take the alias's value + 1.
    ELECTRICAL_ENERGY = ENERGY_CONSUMPTION
  };

  inline constexpr std::string_view kKPIVariableTypeList[] = {
    "Purity",
    "Recovery",
    "Productivity (mass)",
    "Productivity (volume)",
    "Electrical Energy",
    "Production Rate",
    "Thermal Energy",
    "Cooling Duty"
  };

  static std::string GetKPIUnits(KPIVariableType variable) {
    switch (variable) {
      case KPIVariableType::PURITY: return "[mol/mol]";
      case KPIVariableType::RECOVERY: return "[mol/mol]";
      case KPIVariableType::PRODUCTIVITY_MASS: return "[kg_gas/(kg_solid·s)]";
      case KPIVariableType::PRODUCTIVITY_VOLUME: return "[kg_gas/(m3·s)]";
      case KPIVariableType::ENERGY_CONSUMPTION: return "[J/kg]";
      case KPIVariableType::PRODUCTION_RATE: return "[mol/s]";
      case KPIVariableType::THERMAL_ENERGY: return "[J/kg]";
      case KPIVariableType::COOLING_DUTY: return "[J/kg]";
      default: return "[-]";
    }
  }

  // Where the KPI integrals come from: the IDAS quadratures integrated
  // alongside the solution (StreamData::integrals_, kpi_quadrature.h; the
  // default and the reported KPIs) or the trapezoid over the logged boundary
  // samples (KPIRecord::trapezoid_*). The KPI formulas are the same; only the
  // time integrals of the endpoint flows and the pump power differ. An
  // endpoint whose quadratures were not integrated (turned off with
  // IDASettings::kpi_quadratures_, or its step failed) falls back to the
  // trapezoid, so QUADRATURE always gives a value where TRAPEZOID does.
  enum class KPIIntegration { TRAPEZOID, QUADRATURE };

  // Volume [m3] of the packed bed: the cells the domain lists as packed
  // (DomainVectors::bed_cell_indices: no ghost cells, no void end cells), or
  // every interior cell when it lists none. What the per-bed-volume KPIs
  // (productivity_volume, the space-time yield on the bed-volume basis)
  // divide by.
  inline double PackedBedVolume(const DomainVectors& domain) {
    double volume = 0.0;
    if (!domain.bed_cell_indices.empty()) {
      for (int i : domain.bed_cell_indices)
        if (i >= 0 && static_cast<size_t>(i) < domain.volume.size()) volume += domain.volume[static_cast<size_t>(i)];
    } else {
      for (size_t i = 1; i + 1 < domain.volume.size(); ++i) volume += domain.volume[i];
    }
    return volume;
  }

  class KPIParams {
  public:
    KPIParams(const EquationLog& eq_log, const BoundaryLog& b_log) : eq_log_(eq_log), b_log_(b_log) {}

    bool Check() const {
      if (product_comp_.empty()) {
        return false;
      }

      // Need at least one product endpoint and one feed endpoint to attribute
      // throughput. Recycle endpoints are optional.
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
    // throughput counts as feed / product / recycle for this component. Each
    // is summed. Recycle (product fed back into the bed) is subtracted from
    // the product for the net KPIs; purity uses the product alone.
    std::vector<KPIEndpointRef> feed_endpoints{};
    std::vector<KPIEndpointRef> product_endpoints{};
    std::vector<KPIEndpointRef> recycle_endpoints{};
    KPIIntegration integration_ = KPIIntegration::QUADRATURE;
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
  // (The default, QUADRATURE, reads the endpoint's quadrature integrals and
  // falls back to the trapezoid when they were not integrated.)
  // The quadrature integrals of one endpoint over its step, or nullptr when
  // the endpoint was not logged or its quadratures were not integrated.
  inline const EndpointIntegrals* EndpointQuadrature(const adsim::BoundaryLog& b_log, int cycle_index,
                                                     int step_index, const std::string& endpoint_name) {
    const StreamData* d = b_log.GetStreamData(cycle_index, step_index, endpoint_name);
    return (d && d->integrals_.valid) ? &d->integrals_ : nullptr;
  }

  static double EndpointMolesIn(
    const adsim::BoundaryLog& b_log,
    int cycle_index,
    const KPIEndpointRef& ep,
    KPIIntegration integration = KPIIntegration::QUADRATURE
  ) {
    if (integration == KPIIntegration::QUADRATURE) {
      if (const EndpointIntegrals* q = EndpointQuadrature(b_log, cycle_index, ep.step_index, ep.endpoint_name))
        return q->NetMolesIn();
    }
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
    const std::string& component,
    KPIIntegration integration = KPIIntegration::QUADRATURE
  ) {
    if (integration == KPIIntegration::QUADRATURE) {
      if (const EndpointIntegrals* q = EndpointQuadrature(b_log, cycle_index, ep.step_index, ep.endpoint_name))
        return q->NetComponentMolesIn(component);
    }
    adsim::BoundaryVariableData d = b_log.GetComponentVariableData(
      cycle_index, ep.step_index, ep.endpoint_name,
      adsim::BoundaryComponentVariableType::MOLAR_FLOW, component);
    const double integral = IntegrateVariableInTime(d.variable_timeseries, d.time);
    return (ep.location == adsim::BoundaryLocation::BOTTOM) ? integral : -integral;
  }

  // Face-reconstruction of per-component boundary flux (positive = entering the
  // bed), matching the discrete species-equation flux at the boundary face:
  // upwind convection only, since no axial dispersion crosses a boundary face
  // of the bed (Danckwerts; see UpdateMoleFraction). Used by cyclic-steady-state detection, the mass-balance plotter
  // panels and the validation tests. NOT used by the endpoint-keyed feed/product
  // KPIs, which sample per-stream logs instead.
  static double GetTotalComponentMolesIn(
    const adsim::EquationLog& logger,
    int cycle_index,
    int step_index,
    adsim::BoundaryLocation location,
    const std::string& component
  ) {
    // Use the same discrete boundary flux form as the species equation: the
    // convective upwind flux through the boundary face. (The residual carries
    // no dispersive flux there: at an inlet the total flux is u C y_in, at an
    // outlet or closed end the ghost mirrors the first cell.)
    const auto& pressure = logger.GetVariable(cycle_index, step_index, adsim::VariableType::PRESSURE);
    const std::vector<adsim::vec>& P_ = pressure.values;

    const auto& temperature = logger.GetVariable(cycle_index, step_index, adsim::VariableType::TEMPERATURE);
    const std::vector<adsim::vec>& T_ = temperature.values;

    const auto& mole_fraction = logger.GetComponentVariable(cycle_index, step_index, adsim::ComponentVariableType::MOLE_FRACTION, component);
    const std::vector<adsim::vec>& yi_ = mole_fraction.values;


    const adsim::vec& time = pressure.times;

    const auto& domain = logger.GetDomain();
    const auto& A = domain.area;
    const auto& velocity = logger.GetFaceVariable(cycle_index, step_index, adsim::VariableType::VELOCITY).values;

    if (P_.empty()) return 0.0;

    // Guard against mismatched sizes (can occur if data is still being written).
    const size_t n_snapshots = P_.size();
    if (T_.size() != n_snapshots || yi_.size() != n_snapshots ||
        velocity.size() != n_snapshots) {
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

      if (location == adsim::BoundaryLocation::BOTTOM) {
        const double u_face = uvec[0];
        const double yi_face = (u_face > 0.0) ? yivec[0] : yivec[1];
        const double Pw = 0.5 * (Pvec[0] + Pvec[1]);
        const double Tw = 0.5 * (Tvec[0] + Tvec[1]);
        flux_in[i] = u_face * Pw * yi_face / (R * Tw);
      }
      else {
        const int Nf = static_cast<int>(uvec.size()) - 1;
        const int Np = static_cast<int>(Pvec.size()) - 1;      // top ghost cell index
        const int n = Np - 1;                                   // top interior cell index
        const double u_face = uvec[Nf];
        const double yi_face = (u_face > 0.0) ? yivec[n] : yivec[n + 1];
        const double Pe = 0.5 * (Pvec[n] + Pvec[n + 1]);
        const double Te = 0.5 * (Tvec[n] + Tvec[n + 1]);
        flux_in[i] = -u_face * Pe * yi_face / (R * Te);
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

  // Amount accumulated in the bed over the logged step: the change of
  // sum_i value[i] * volume[i] * factor[i] over the interior cells between the
  // first and the last logged instant (0 when nothing was logged). The factor
  // turns the per-cell value into an amount: the void fraction for a gas
  // concentration, the bulk density for an adsorbed loading.
  static double AccumulatedInBed(
    const adsim::VariableLogEntry& variable,
    const std::vector<double>& volume,
    const std::vector<double>& factor
  ) {
    if (variable.values.empty()) {
      return 0.0;
    }
    const auto& values_t0 = variable.values[0];
    const auto& values_tend = variable.values.back();

    double accumulation = 0.0;
    for (size_t i = 1; i + 1 < values_t0.size(); ++i) { // Skip ghost cells
      accumulation += (values_tend[i] - values_t0[i]) * volume[i] * factor[i];
    }
    return accumulation;
  }

  // Gas-phase accumulation [mol] of a component: molar concentration times the
  // void volume of each cell.
  static double CalculateComponentAccumulation_Volume(
    const adsim::EquationLog& logger,
    int cycle_index,
    int step_index,
    adsim::ComponentVariableType variable_type,
    const std::string& component
  ) {
    const auto& domain = logger.GetDomain();
    return AccumulatedInBed(logger.GetComponentVariable(cycle_index, step_index, variable_type, component),
                            domain.volume, domain.et);
  }

  // Adsorbed-phase accumulation [mol] of a component: loading [mol/kg] times
  // the adsorbent mass of each cell.
  static double CalculateComponentAccumulation_Mass(
    const adsim::EquationLog& logger,
    int cycle_index,
    int step_index,
    adsim::ComponentVariableType variable_type,
    const std::string& component
  ) {
    const auto& domain = logger.GetDomain();
    return AccumulatedInBed(logger.GetComponentVariable(cycle_index, step_index, variable_type, component),
                            domain.volume, domain.bulk_density);
  }

  // Net moles leaving the bed across the selected product endpoints. Signed sum:
  // an outlet (net leaving) adds; flow back into the bed through a product
  // endpoint subtracts. For a pure outlet the signed sum equals the gross
  // throughput. Product fed back through other endpoints (a purge or a
  // light-product pressurisation) belongs in the recycle role, see
  // CalculateNetComponentMolesInProduct.
  static double CalculateTotalMolesInProduct(const KPIParams& params) {
    if (!params.Check()) {
      return 0.0;
    }

    double moles_in_product = 0.0;
    for (const auto& ep : params.product_endpoints) {
      double moles_in = EndpointMolesIn(params.b_log_, params.cycle_index_, ep, params.integration_);
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
        params.b_log_, params.cycle_index_, ep, params.product_comp_, params.integration_);
      product_moles += -moles_in;  // component leaving = positive product
    }

    return product_moles;
  }

  // Component moles the recycle endpoints carry into the bed (entering =
  // positive).
  static double CalculateComponentMolesRecycled(const KPIParams& params) {
    if (!params.Check()) {
      return 0.0;
    }

    double recycled = 0.0;
    for (const auto& ep : params.recycle_endpoints) {
      recycled += EndpointComponentMolesIn(params.b_log_, params.cycle_index_, ep, params.product_comp_,
                                           params.integration_);
    }

    return recycled;
  }

  // Product the process delivers: product drawn minus product recycled.
  // Recovery, productivity, production rate and specific energy use this.
  static double CalculateNetComponentMolesInProduct(const KPIParams& params) {
    return CalculateComponentMolesInProduct(params) - CalculateComponentMolesRecycled(params);
  }

  static double CalculateComponentMolesInFeed(const KPIParams& params) {
    if (!params.Check()) {
      return 0.0;
    }

    // Net component moles entering the bed across the selected feed endpoints.
    double component_moles_in_feed = 0.0;
    for (const auto& ep : params.feed_endpoints) {
      double moles_in = EndpointComponentMolesIn(
        params.b_log_, params.cycle_index_, ep, params.product_comp_, params.integration_);
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
    double moles_in_product = ZeroIfRoundoff(CalculateNetComponentMolesInProduct(params), moles_in_feed);
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
    double product_moles = CalculateNetComponentMolesInProduct(params);    
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
    // The packed bed's volume: no ghost cells (and no void end cells).
    const double volume = PackedBedVolume(params.eq_log_.GetDomain());

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
    const double moles_out = CalculateNetComponentMolesInProduct(params);
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
      if (params.integration_ == KPIIntegration::QUADRATURE) {
        if (const EndpointIntegrals* q = EndpointQuadrature(b_log, cycle_index, step_idx, stream_name)) {
          step_energy += q->energy;
          continue;
        }
      }
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

  // Specific electrical energy [J/kg]: the time-integrated compressor and
  // vacuum-pump power of every step over the net product mass. Step energy
  // sources (heaters, coolers) are not included. This is the
  // `energy_consumption` KPI, also reported as `electrical_energy`.
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
    double product_moles = CalculateNetComponentMolesInProduct(params);    
    double product_mass = product_moles * mr.at(product_component); // [kg_gas]
    double total_specific_energy = KPIRatio(total_energy, product_mass, ok); // [J/kg_gas]
    
    return KPIResult(total_specific_energy, ok, defined);
  }

  // The specific electrical energy under its explicit name.
  static double CalculateSpecificElectricalEnergy(const KPIParams& params, bool* defined = nullptr) {
    return CalculateSpecificEnergyConsumption(params, defined);
  }

  // ---- Thermal energy and cooling duty -------------------------------------
  //
  // The thermal energy is the heat supplied to the bed over the cycle; the
  // cooling duty is the heat removed, reported separately (as a positive
  // magnitude) and not subtracted from the thermal energy. Both come from the
  // heat duty the solver logged for each step (StepHeatDuty, zero in an
  // isothermal run): the step energy sources' powers times the step's
  // duration, plus the heat the step's heating medium passed into (thermal
  // energy) and out of (cooling duty) the column, through the jacket and any
  // internal heat exchanger, integrated over the step. The
  // electrical energy (compressor and pump work) is separate.

  // Duration [s] of a logged step: from its first logged time (the step start)
  // to the next step's start, or to the cycle's end for the last step. 0 when
  // the step was not logged.
  inline double LoggedStepDuration(const EquationLog& eq_log, int cycle_index, int step_idx) {
    const auto& times = eq_log.GetLog(cycle_index, step_idx).step_data.times;
    if (times.empty()) return 0.0;
    const double start = times.front();
    const auto& next = eq_log.GetLog(cycle_index, step_idx + 1).step_data.times;
    double end = 0.0;
    if (!next.empty()) {
      end = next.front();
    } else {
      const auto& sim_log = eq_log.GetSimulationLog();
      const auto it = sim_log.find(cycle_index);
      if (it == sim_log.end()) return 0.0;
      double cycle_start = start;  // the first logged time of the cycle
      for (const auto& entry : it->second) {
        const auto& t = entry.second.step_data.times;
        if (!t.empty()) { cycle_start = t.front(); break; }
      }
      end = cycle_start + eq_log.sys_data.total_cycle_time;
    }
    return std::max(0.0, end - start);
  }

  struct StepThermalDuty {
    double heating = 0.0;  // [J] heat supplied over the step (energy sources + medium)
    double cooling = 0.0;  // [J] heat removed over the step (>= 0; energy sources + medium)
    double medium_heating = 0.0;  // [J] the heating medium's share of `heating`
    double medium_cooling = 0.0;  // [J] the heating medium's share of `cooling`
  };

  // Heat supplied and removed over one step of the params' cycle [J]. Needs no
  // KPI configuration: it does not depend on the product.
  static StepThermalDuty CalculateStepThermalDuty(const KPIParams& params, int step_idx) {
    const auto& eq_log = params.eq_log_;
    const StepHeatDuty& duty = eq_log.GetLog(params.cycle_index_, step_idx).heat_duty;
    const double dt = LoggedStepDuration(eq_log, params.cycle_index_, step_idx);
    return {duty.heating_power * dt + duty.medium_heat_in, duty.cooling_power * dt + duty.medium_heat_out,
            duty.medium_heat_in, duty.medium_heat_out};
  }

  // Heat supplied and removed over the whole cycle [J].
  static StepThermalDuty CalculateCycleThermalDuty(const KPIParams& params) {
    StepThermalDuty total;
    const int n_steps = static_cast<int>(params.eq_log_.sys_data.step_list_.size());
    for (int s = 0; s < n_steps; ++s) {
      const StepThermalDuty d = CalculateStepThermalDuty(params, s);
      total.heating += d.heating;
      total.cooling += d.cooling;
      total.medium_heating += d.medium_heating;
      total.medium_cooling += d.medium_cooling;
    }
    return total;
  }

  // Net product mass of the params' cycle [kg]; `ok` false when undefined.
  static double NetProductMass(const KPIParams& params, bool& ok) {
    if (!params.Check()) { ok = false; return 0.0; }
    const double moles = CalculateNetComponentMolesInProduct(params);
    const double mass = moles * params.eq_log_.sys_data.molar_mass.at(params.product_comp_);
    if (!std::isfinite(mass) || !(mass > 0.0)) ok = false;
    return mass;
  }

  // Specific thermal energy [J/kg]: heat supplied over the cycle per kg of net
  // product. Defined and undefined exactly as the electrical energy.
  static double CalculateSpecificThermalEnergy(const KPIParams& params, bool* defined = nullptr) {
    if (!params.Check()) return KPIResult(0.0, false, defined);
    bool ok = true;
    const double mass = NetProductMass(params, ok);
    const double value = KPIRatio(CalculateCycleThermalDuty(params).heating, mass, ok);
    return KPIResult(value, ok, defined);
  }

  // Specific cooling duty [J/kg]: heat removed over the cycle per kg of net
  // product, as a positive magnitude. For information; not part of the
  // thermal energy.
  static double CalculateSpecificCoolingDuty(const KPIParams& params, bool* defined = nullptr) {
    if (!params.Check()) return KPIResult(0.0, false, defined);
    bool ok = true;
    const double mass = NetProductMass(params, ok);
    const double value = KPIRatio(CalculateCycleThermalDuty(params).cooling, mass, ok);
    return KPIResult(value, ok, defined);
  }

}
