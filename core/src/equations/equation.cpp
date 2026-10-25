/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/equations/equation.h"
#include "adsim/objects/user_data.h"
#include "adsim/boundary_conditions/step.h"

#include <algorithm>


namespace adsim {

  namespace {
    // Transport properties some selected closure model reads during the
    // residual (see transport_properties.h), closed over their dependencies:
    // Sh reads Sc and Re, Sc reads Dm.
    unsigned TransportPropertiesRead(const Bed& bed, const AdsorbentDirectory& dir) {
      unsigned reads = kNoTransportProperties;
      const auto& particle = bed.equations_.particle_heat_transfer_model_;
      if (bed.equations_.thermal_equilibrium_mode_ == ThermalEquilibriumMode::NON_LOCAL && particle)
        reads |= particle->ReadsTransportProperties();
      for (const auto& layer : bed.layers_.layers_) {
        const auto& eq = layer.equations_;
        if (eq.dispersion_model_) reads |= eq.dispersion_model_->ReadsTransportProperties();
        if (eq.heat_dispersion_model_) reads |= eq.heat_dispersion_model_->ReadsTransportProperties();
        if (eq.internal_heat_transfer_model_) reads |= eq.internal_heat_transfer_model_->ReadsTransportProperties();
        if (eq.external_heat_transfer_model_) reads |= eq.external_heat_transfer_model_->ReadsTransportProperties();
        for (const auto& [ads_id, fraction] : layer.composition_) {
          for (const auto& [component, model] : dir.Get(ads_id).adsorption_source_.mass_transfer_models_)
            if (model && model->mtc_model_) reads |= model->mtc_model_->ReadsTransportProperties();
        }
      }
      if (reads & kSherwoodNumber) reads |= kSchmidtNumber | kReynoldsNumber;
      if (reads & kSchmidtNumber) reads |= kMolecularDiffusion;
      return reads;
    }
  } // namespace

  // IndependentEquation implementations
  void IndependentEquation::Resize(size_t n_cells) {
    if (cell_type == CellType::CellFace) 
      n_cells--;
    var.resize(n_cells);
  }

  void IndependentEquation::Update(UserData* udata) {
    update_func(udata);
  }

  size_t IndependentEquation::Size() const {
    return var.size();
  }

  // Equation implementations
  void Equation::SetInitialConditions(double x) {
    for (size_t i = 0; i < var.size(); i++) {
      var[i] = x;
      ddt_var[i] = 0.0;
      res_var[i] = 0.0;
    }
  }

  void Equation::SetSmoothInitialConditions(double x_inlet, double x_outlet) {
    for (size_t i = 0; i < var.size(); i++) {
      double xi = static_cast<double>(i) / (var.size() - 1);
      var[i] = x_inlet + (x_outlet - x_inlet) * xi;
      ddt_var[i] = 0.0;
      res_var[i] = 0.0;
    }
  }

  void Equation::Resize(size_t n_cells) {
    if (cell_type == CellType::CellFace)
      n_cells--;

    var.resize(n_cells);
    ddt_var.resize(n_cells);
    res_var.resize(n_cells);
  }

  bool Equation::IsActiveAt(size_t cell_idx, const std::vector<bool>* is_void_mask) const {
    if (cell_idx >= var.size()) return false;
    if (domain == EquationDomain::PorousMediaOnly
        && is_void_mask
        && cell_idx < is_void_mask->size()
        && (*is_void_mask)[cell_idx]) {
      return false;
    }
    if (active_mask
        && cell_idx < active_mask->size()
        && !(*active_mask)[cell_idx]) {
      return false;
    }
    return true;
  }

  void Equation::SetID(size_t& start_index, sunrealtype* id) {
    double x = 1.0;  // 1.0 = differential equation
    if (equation_type == EquationType::Algebraic) x = 0.0;  // 0.0 = algebraic constraint

    size_t N = Size();
    size_t end_idx = start_index + N;
    
    id[start_index] = 0.0; // Ghost cell
    for (size_t i = start_index + 1; i < end_idx - 1; i++) {
      id[i] = x;
    }
    id[end_idx - 1] = 0.0; // Ghost cell
    
    start_index += N;
  }

  double Equation::GetIDValue(size_t cell_idx) const {
    // Ghost cells (first and last) always have ID = 0.0
    size_t N = Size();
    if (cell_idx == 0 || cell_idx == N - 1) {
      return 0.0;
    }
    
    // Interior cells: differential = 1.0, algebraic = 0.0
    if (equation_type == EquationType::Algebraic) {
      return 0.0;
    }
    return 1.0;
  }

  void Equation::Pack(size_t& start_index, sunrealtype* y_data, sunrealtype* yp_data) const {
    size_t N = Size();
    for (size_t i = 0; i < N; i++) {
      y_data[start_index + i] = var[i];
      yp_data[start_index + i] = ddt_var[i];
    }
    start_index += N;
  }
  
  void Equation::PackResiduals(size_t& start_index, sunrealtype* res_data) const {
    size_t N = Size();
    for (size_t i = 0; i < N; i++) {
      res_data[start_index + i] = res_var[i];
    }
    start_index += N;
  }

  void Equation::Unpack(size_t& start_index, const sunrealtype* y_data, const sunrealtype* yp_data) {
    size_t N = Size();
    for (size_t i = 0; i < N; i++) {
      var[i] = y_data[start_index + i];
      ddt_var[i] = yp_data[start_index + i];
    }
    start_index += N;
  }
  
  void Equations::Initialise(const Fluid& fluid, const Bed& bed, const AdsorbentDirectory& dir) {
    domain_vectors.Initialise(bed, dir, fluid);
    equations_.clear();
    independent_equations_.clear();
    banded_slots_.clear();

    // Build the set of adsorbents present in the bed and their per-cell presence
    // masks. A per-adsorbent loading equation is active only in the cells of the
    // layers whose composition contains that adsorbent, plus the adjacent ghost
    // cell when that adsorbent occupies the first/last bed cell (matching the
    // ghost-mirroring residual in UpdateAdsorbedConcentration). For a
    // single-adsorbent bed this reproduces the legacy active-cell set exactly,
    // so the packed state vector is byte-for-byte unchanged.
    adsorbent_ids_.clear();
    adsorbent_cell_mask_.clear();
    {
      const int N = bed.GetNumberOfCells();
      for (const auto& layer : bed.layers_.layers_) {
        for (const auto& [ads_id, frac] : layer.composition_) {
          auto& mask = adsorbent_cell_mask_[ads_id];
          if (mask.empty()) mask.assign(N, false);
          for (int i = layer.domain_.start_index; i <= layer.domain_.end_index; ++i) {
            if (i >= 0 && i < N) mask[i] = true;
          }
        }
      }
      const int first_bed = bed.GetFirstBedCellIndex();
      const int last_bed = bed.GetLastBedCellIndex();
      for (auto& [ads_id, mask] : adsorbent_cell_mask_) {
        if (first_bed >= 0 && first_bed < N && mask[first_bed]) mask[0] = true;
        if (last_bed >= 0 && last_bed < N && mask[last_bed]) mask[N - 1] = true;
        adsorbent_ids_.push_back(ads_id);
      }
      std::sort(adsorbent_ids_.begin(), adsorbent_ids_.end());
    }

    const auto thermal_mode = bed.equations_.thermal_equilibrium_mode_;

    // Set solid_temperature equation type based on thermal mode
    if (thermal_mode == ThermalEquilibriumMode::LOCAL) {
      solid_temperature.SetEquationType(EquationType::Algebraic);
    } else {
      solid_temperature.SetEquationType(EquationType::Differential);
    }

    // Solid temperature only has physical meaning inside the packed bed.
    solid_temperature.SetDomain(EquationDomain::PorousMediaOnly);

    // Isothermal runs drop the three temperature unknowns from the DAE: they
    // are not registered, so they are never packed, unpacked, given residuals
    // or Jacobian columns, and HoldIsothermalTemperature keeps their values at
    // T_iso (which is what every residual, isotherm and log reads).
    isothermal_ = bed.equations_.IsIsothermal();
    isothermal_temperature_ = 0.0;

    // Reduced unknowns (IDASettings): in LTE the solid temperature is the
    // algebraic constraint Ts = Tf, so it can be a copy of Tf instead of an
    // unknown; the gas concentration is C = P / (R T) in every cell (its time
    // derivative appears in no residual), so it can be derived from P and T.
    // Unregistered fields keep their storage and are refreshed from the
    // unpacked state by SyncDerivedFields.
    solid_temperature_aliased_ = alias_lte_solid_temperature_ && !isothermal_
                                 && thermal_mode == ThermalEquilibriumMode::LOCAL;

    // Registration order matters!
    if (!derived_concentration_) equations_.emplace_back(&concentration);
    equations_.emplace_back(&pressure);
    if (!isothermal_) {
      equations_.emplace_back(&fluid_temperature);
      if (!solid_temperature_aliased_) equations_.emplace_back(&solid_temperature);
      equations_.emplace_back(&wall_temperature);
    }
    equations_.emplace_back(&velocity);

    mole_fraction.clear();
    for (const auto& component : fluid.components_) {
      auto [it, inserted] = mole_fraction.emplace(
        std::piecewise_construct,
        std::forward_as_tuple(component),
        std::forward_as_tuple(
          component + " Mole Fraction",
          CellType::CellCentre,
          EquationType::Differential, 
          [component](UserData* udata) { UpdateMoleFraction(udata, component); },
          1e-8
        )
      );
      equations_.emplace_back(&it->second);
    }

    ads_concentration.clear();
    for (AdsorbentID ads_id : adsorbent_ids_) {
      const std::vector<bool>* mask = GetAdsorbentMask(ads_id);
      auto& comp_map = ads_concentration[ads_id];
      for (const auto& component : fluid.components_) {
        auto [it, inserted] = comp_map.emplace(
          std::piecewise_construct,
          std::forward_as_tuple(component),
          std::forward_as_tuple(
            component + " Adsorbed Concentration",
            CellType::CellCentre,
            EquationType::Differential,
            [ads_id, component](UserData* udata) { UpdateAdsorbedConcentration(udata, ads_id, component); },
            1e-8
          )
        );
        it->second.SetDomain(EquationDomain::PorousMediaOnly);
        it->second.SetActiveMask(mask);
        equations_.emplace_back(&it->second);
      }
    }

    molar_concentration.clear();
    for (const auto& component : fluid.components_) {
      auto [it, inserted] = molar_concentration.emplace(
        std::piecewise_construct,
        std::forward_as_tuple(component),
        std::forward_as_tuple(
          component + " Molar Concentration", 
          CellType::CellCentre, 
          [component](UserData* udata) { UpdateMolarConcentration(udata, component); }
        )
      );
      independent_equations_.emplace_back(&it->second);
    }

    molecular_diffusion.clear();
    for (const auto& component : fluid.components_) {
      auto [it, inserted] = molecular_diffusion.emplace(
        std::piecewise_construct,
        std::forward_as_tuple(component),
        std::forward_as_tuple(
          component + " Molecular Diffusion Coefficient",
          CellType::CellCentre,
          [component](UserData* udata) { UpdateMolecularDiffusionCoefficient(udata, component); }
        )
      );
      independent_equations_.emplace_back(&it->second);
    }

    independent_equations_.emplace_back(&total_mole_fraction);
    independent_equations_.emplace_back(&fluid_density);
    independent_equations_.emplace_back(&viscosity);
    independent_equations_.emplace_back(&heat_capacity);
    independent_equations_.emplace_back(&reynolds_number);
    independent_equations_.emplace_back(&schmidt_number);
    independent_equations_.emplace_back(&sherwood_number);
    independent_equations_.emplace_back(&fluid_thermal_conductivity);
    independent_equations_.emplace_back(&axial_thermal_conductivity);
    independent_equations_.emplace_back(&internal_heat_transfer_coefficient);
    independent_equations_.emplace_back(&external_heat_transfer_coefficient);
    if (thermal_mode == ThermalEquilibriumMode::NON_LOCAL) {
      independent_equations_.emplace_back(&particle_heat_transfer_coefficient);
    }
    independent_equations_.emplace_back(&molar_mass);
    independent_equations_.emplace_back(&volumetric_flow_rate);
    independent_equations_.emplace_back(&molar_flow_rate);
    independent_equations_.emplace_back(&approach_to_fluidisation);
    independent_equations_.emplace_back(&pressure_gradient);
    independent_equations_.emplace_back(&axial_dispersion_coefficient);

    mass_transfer_coefficient.clear();
    for (AdsorbentID ads_id : adsorbent_ids_) {
      auto& comp_map = mass_transfer_coefficient[ads_id];
      for (const auto& component : fluid.components_) {
        auto [it, inserted] = comp_map.emplace(
          std::piecewise_construct,
          std::forward_as_tuple(component),
          std::forward_as_tuple(
            component + " Mass Transfer Coefficient",
            CellType::CellCentre,
            [ads_id, component](UserData* udata) { UpdateMassTransferCoefficient(udata, ads_id, component); }
          )
        );
        independent_equations_.emplace_back(&it->second);
      }
    }

    eq_ads_concentration.clear();
    for (AdsorbentID ads_id : adsorbent_ids_) {
      auto& comp_map = eq_ads_concentration[ads_id];
      for (const auto& component : fluid.components_) {
        auto [it, inserted] = comp_map.emplace(
          std::piecewise_construct,
          std::forward_as_tuple(component),
          std::forward_as_tuple(
            component + " Equilibrium Adsorbed Concentration",
            CellType::CellCentre,
            [ads_id, component](UserData* udata) { UpdateEquilibriumAdsorbedConcentration(udata, ads_id, component); }
          )
        );
        independent_equations_.emplace_back(&it->second);
      }
    }

    component_mass_source.clear();
    for (const auto& component : fluid.components_) {
      auto [it, inserted] = component_mass_source.emplace(
        std::piecewise_construct,
        std::forward_as_tuple(component),
        std::forward_as_tuple(
          component + " Component Mass Source",
          CellType::CellCentre,
          [component](UserData* udata) { UpdateComponentMassSource(udata, component); }
        )
      );
      independent_equations_.emplace_back(&it->second);
    }

    component_energy_source.clear();
    for (const auto& component : fluid.components_) {
      auto [it, inserted] = component_energy_source.emplace(
        std::piecewise_construct,
        std::forward_as_tuple(component),
        std::forward_as_tuple(
          component + " Component Energy Source",
          CellType::CellCentre,
          [component](UserData* udata) { UpdateComponentEnergySource(udata, component); }
        )
      );
      independent_equations_.emplace_back(&it->second);
    }

    // The totals sum the per-component sources, so register them after those.
    independent_equations_.emplace_back(&mass_source);
    independent_equations_.emplace_back(&energy_source);

    mole_fraction_face.clear();
    for (const auto& component : fluid.components_) {
      auto [it, inserted] = mole_fraction_face.emplace(
        std::piecewise_construct,
        std::forward_as_tuple(component),
        std::forward_as_tuple(
          component + " Mole Fraction Face",
          CellType::CellFace,
          [component](UserData* udata) { UpdateMoleFractionFace(udata, component); }
        )
      );
      independent_equations_.emplace_back(&it->second);
    }

    // Logging-only fields: no residual reads them, so they are skipped on
    // residual evaluations and computed only when the state is logged. Ordered
    // so each reads fields already updated in the same pass (pressure gradient
    // before approach to fluidisation, component sources before their totals).
    // Transport properties no selected model reads are only needed for the
    // log. Dm before Sc before Sh: each reads the one before.
    const unsigned transport_reads = TransportPropertiesRead(bed, dir);
    diagnostic_equations_.clear();
    if (!(transport_reads & kMolecularDiffusion))
      for (const auto& component : fluid.components_)
        diagnostic_equations_.push_back(&molecular_diffusion.at(component));
    if (!(transport_reads & kReynoldsNumber)) diagnostic_equations_.push_back(&reynolds_number);
    if (!(transport_reads & kSchmidtNumber)) diagnostic_equations_.push_back(&schmidt_number);
    if (!(transport_reads & kSherwoodNumber)) diagnostic_equations_.push_back(&sherwood_number);
    if (!(transport_reads & kGasThermalConductivity)) diagnostic_equations_.push_back(&fluid_thermal_conductivity);

    diagnostic_equations_.push_back(&volumetric_flow_rate);
    diagnostic_equations_.push_back(&molar_flow_rate);
    for (const auto& component : fluid.components_)
      diagnostic_equations_.push_back(&component_mass_source.at(component));
    diagnostic_equations_.push_back(&mass_source);
    for (const auto& component : fluid.components_)
      diagnostic_equations_.push_back(&component_energy_source.at(component));
    diagnostic_equations_.push_back(&energy_source);
    diagnostic_equations_.push_back(&pressure_gradient);
    diagnostic_equations_.push_back(&approach_to_fluidisation);
    diagnostic_equations_.push_back(&total_mole_fraction);

    residual_independent_equations_.clear();
    for (auto* eq : independent_equations_) {
      if (std::find(diagnostic_equations_.begin(), diagnostic_equations_.end(), eq) == diagnostic_equations_.end())
        residual_independent_equations_.push_back(eq);
    }

    // Cache the adsorbed-phase terms consumed by pressure and energy residuals.
    // Iterate the maps in the same order as the old per-residual construction:
    // summation order is observable in floating-point results.
    adsorbent_slot_.clear();
    loading_rate_terms_.clear();
    energy_terms_.clear();
    loading_rate_terms_.reserve(ads_concentration.size());
    energy_terms_.reserve(ads_concentration.size());
    for (const auto& [ads_id, comp_map] : ads_concentration) {
      const size_t slot = loading_rate_terms_.size();
      adsorbent_slot_.emplace(ads_id, slot);
      auto& rates = loading_rate_terms_.emplace_back();
      auto& terms = energy_terms_.emplace_back();
      rates.reserve(comp_map.size());
      terms.reserve(comp_map.size());
      const auto& dH = domain_vectors.dH_ads.at(ads_id);
      for (const auto& [component, qi] : comp_map) {
        rates.push_back(&qi.ddt_var);
        const double Cpa = fluid.adsorbed_heat_capacity_mode_ == AdsorbedHeatCapacityMode::GAS_PHASE
          ? fluid.heat_capacity_.at(component) * fluid.molar_mass_.at(component)
          : fluid.adsorbed_heat_capacity_.at(component).Get();
        terms.push_back({ Cpa, &qi.var, &qi.ddt_var, &dH.at(component) });
      }
    }
  }

  void Equation::Update(UserData* udata) { 
    update_func(udata); 
  }

  size_t Equation::Size() const { 
    return var.size(); 
  }

  void Equations::UpdateDifferentialEquations(UserData* udata) {
    for (const auto& eq : equations_) {
      eq->Update(udata);
    }
  }

  void Equations::Resize(size_t n_cells) {
    banded_slots_.clear();
    for (const auto& eq : equations_) {
      eq->Resize(n_cells);
    }
    if (isothermal_) {
      // Not registered (see Initialise) but still read by the residuals,
      // isotherms and logs, so they need their storage.
      fluid_temperature.Resize(n_cells);
      solid_temperature.Resize(n_cells);
      wall_temperature.Resize(n_cells);
    }
    if (solid_temperature_aliased_) solid_temperature.Resize(n_cells);
    if (derived_concentration_) concentration.Resize(n_cells);

    for (const auto& eq : independent_equations_) {
      eq->Resize(n_cells);
    }
  }

  void Equations::HoldIsothermalTemperature(double T) {
    if (!isothermal_) return;
    isothermal_temperature_ = T;
    fluid_temperature.SetInitialConditions(T);
    solid_temperature.SetInitialConditions(T);
    wall_temperature.SetInitialConditions(T);
    // Keep the gas concentration consistent with the equation of state at the
    // isothermal temperature (it may differ from the IC temperature).
    for (size_t i = 0; i < concentration.var.size() && i < pressure.var.size(); ++i)
      concentration.var[i] = pressure.var[i] / (R * T);
  }

  int Equations::GetSystemSize() {
    const auto* void_mask = &domain_vectors.is_void_cell;
    int size = 0;
    for (const auto& eq : equations_) {
      const size_t N = eq->Size();
      for (size_t i = 0; i < N; ++i) {
        if (eq->IsActiveAt(i, void_mask)) size++;
      }
    }
    return size;
  }

  void Equations::SetInitialConditions(InitialConditions& ic, ThermalEquilibriumMode mode) {
    double C_0 = ic.pressure / (R * ic.temperature);

    concentration.SetInitialConditions(C_0);
    pressure.SetInitialConditions(ic.pressure);
    fluid_temperature.SetInitialConditions(ic.temperature);
    solid_temperature.SetInitialConditions(ic.temperature);
    wall_temperature.SetInitialConditions(ic.wall_temperature);
    velocity.SetInitialConditions(ic.velocity);
    
    for (auto& [component, eq] : mole_fraction){
      double yi_0 = 0.0;
      if (ic.mole_fraction.find(component) != ic.mole_fraction.end()) {
        yi_0 = ic.mole_fraction.at(component).Get();
      }
      eq.SetInitialConditions(yi_0);
    }
    
    // Initialize adsorbed concentration from initial conditions, per adsorbent.
    // Note: Will be updated to equilibrium value after isotherm calculation.
    for (auto& [ads_id, comp_map] : ads_concentration) {
      for (auto& [component, eq] : comp_map) {
        double qi_0 = 0.0;
        if (ic.adsorbed_concentration.find(component) != ic.adsorbed_concentration.end()) {
          qi_0 = ic.adsorbed_concentration.at(component);
        }
        eq.SetInitialConditions(qi_0);
      }
    }
  }

  void Equations::SetAdsorbedPhaseInitialConditions() {
    // Initialize each adsorbent's loading to its own equilibrium value (from its
    // own isotherm). Each adsorbent starts in equilibrium with the gas-phase IC.
    for (auto& [ads_id, comp_map] : ads_concentration) {
      for (auto& [component, ads_eq] : comp_map) {
        auto& qi_eq = ads_eq.var;
        const auto& qi_sat_eq = eq_ads_concentration.at(ads_id).at(component).var;
        for (size_t i = 0; i < qi_eq.size(); ++i) {
          qi_eq[i] = qi_sat_eq[i];
          ads_eq.ddt_var[i] = 0.0;  // Initial time derivative is zero
        }
      }
    }
  }
  
  void Equations::SetID(sunrealtype* id) {
    assert(id);
    size_t start_index = 0;
    for (const auto& eq : equations_) {
      eq->SetID(start_index, id);
    }
  }

  void Equations::Pack(sunrealtype* y_data, sunrealtype* yp_data) const {
    assert(y_data);
    assert(yp_data);

    size_t start_index = 0;
    for (const auto& eq : equations_) {
      eq->Pack(start_index, y_data, yp_data);
    }
  }
  
  void Equations::PackResiduals(sunrealtype* res_data) const {
    assert(res_data);
    size_t start_index = 0;
    for (const auto& eq : equations_) {
      eq->PackResiduals(start_index, res_data);
    }
  }
  
  void Equations::Unpack(const sunrealtype* y_data, const sunrealtype* yp_data) {
    assert(y_data);
    assert(yp_data);

    size_t start_index = 0;
    for (const auto& eq : equations_) {
      eq->Unpack(start_index, y_data, yp_data);
    }
    SyncDerivedFields();
  }

  void Equations::SyncDerivedFields() {
    if (solid_temperature_aliased_) {
      solid_temperature.var = fluid_temperature.var;
      solid_temperature.ddt_var = fluid_temperature.ddt_var;
    }
    if (derived_concentration_) {
      // Same expression (and temperature floor) as the former residual
      // C - P / (R T); dC/dt by the chain rule, for anything that logs it.
      const auto& P = pressure.var;
      const auto& dP = pressure.ddt_var;
      const auto& T = fluid_temperature.var;
      const auto& dT = fluid_temperature.ddt_var;
      auto& C = concentration.var;
      auto& dC = concentration.ddt_var;
      const size_t n = std::min(C.size(), std::min(P.size(), T.size()));
      for (size_t i = 0; i < n; ++i) {
        const double Tc = detail::DenomFloor(T[i], constants::kMinDenomTemperature);
        const double inv_RT = 1.0 / (R * Tc);
        C[i] = P[i] * inv_RT;
        dC[i] = (dP[i] - P[i] * dT[i] / Tc) * inv_RT;
      }
    }
  }

  void Equations::SetIDBanded(sunrealtype* id) {
    assert(id);
    const auto* void_mask = &domain_vectors.is_void_cell;

    size_t max_size = 0;
    for (const auto& eq : equations_) {
      max_size = std::max(max_size, eq->Size());
    }

    size_t global_index = 0;
    for (size_t cell_idx = 0; cell_idx < max_size; ++cell_idx) {
      for (const auto& eq : equations_) {
        if (!eq->IsActiveAt(cell_idx, void_mask)) continue;
        id[global_index] = eq->GetIDValue(cell_idx);
        global_index++;
      }
    }
  }

  const std::vector<Equations::BandedSlot>& Equations::BandedSlots() const {
    if (!banded_slots_.empty()) return banded_slots_;

    const auto* void_mask = &domain_vectors.is_void_cell;
    size_t max_size = 0;
    for (const auto& eq : equations_) {
      max_size = std::max(max_size, eq->Size());
    }

    for (size_t cell_idx = 0; cell_idx < max_size; ++cell_idx) {
      for (const auto& eq : equations_) {
        if (!eq->IsActiveAt(cell_idx, void_mask)) continue;
        banded_slots_.push_back({ &eq->var[cell_idx], &eq->ddt_var[cell_idx], &eq->res_var[cell_idx], eq, cell_idx });
      }
    }
    return banded_slots_;
  }

  std::string Equations::DescribeBandedIndex(size_t g) const {
    const auto& slots = BandedSlots();
    if (g >= slots.size()) return "?";
    return slots[g].eq->name + " @ cell " + std::to_string(slots[g].cell);
  }

  void Equations::PackBanded(sunrealtype* y_data, sunrealtype* yp_data) const {
    assert(y_data);
    assert(yp_data);
    const auto& slots = BandedSlots();
    for (size_t g = 0; g < slots.size(); ++g) {
      y_data[g] = *slots[g].var;
      yp_data[g] = *slots[g].ddt_var;
    }
  }

  void Equations::PackResidualsBanded(sunrealtype* res_data) const {
    assert(res_data);
    const auto& slots = BandedSlots();
    for (size_t g = 0; g < slots.size(); ++g) {
      res_data[g] = *slots[g].res_var;
    }
  }

  void Equations::UnpackBanded(const sunrealtype* y_data, const sunrealtype* yp_data) {
    assert(y_data);
    assert(yp_data);
    const auto& slots = BandedSlots();
    for (size_t g = 0; g < slots.size(); ++g) {
      *slots[g].var = y_data[g];
      *slots[g].ddt_var = yp_data[g];
    }
    SyncDerivedFields();
  }
  
  void Equations::UpdateIndependentEquations(UserData* udata) {
    if (!udata) return;

    for (const auto& eq : residual_independent_equations_) {
      eq->Update(udata);
    }
  }

  void Equations::UpdateDiagnosticEquations(UserData* udata) {
    if (!udata) return;

    for (const auto& eq : diagnostic_equations_) {
      eq->Update(udata);
    }
  }

  void Equations::UpdateBoundaryConditions(UserData* udata, double time) {
    assert(current_step_);
    if (!current_step_) return;  // assert is compiled out in Release; don't deref
    auto& state = udata->state;
    current_step_->bottom_boundary_condition_->Apply(state, time);
    current_step_->top_boundary_condition_->Apply(state, time);
  }

  void Equations::SetToleranceVector(sunrealtype* abstol_data) const {
    size_t idx = 0;
    for (const auto& eq : equations_) {
      for (size_t i = 0; i < eq->Size(); i++) {
        abstol_data[idx++] = eq->abstol;
      }
    }
  }

  void Equations::SetToleranceVectorBanded(sunrealtype* abstol_data) const {
    const auto* void_mask = &domain_vectors.is_void_cell;

    size_t max_size = 0;
    for (const auto& eq : equations_) {
      max_size = std::max(max_size, eq->Size());
    }

    size_t global_index = 0;
    for (size_t cell_idx = 0; cell_idx < max_size; ++cell_idx) {
      for (const auto& eq : equations_) {
        if (!eq->IsActiveAt(cell_idx, void_mask)) continue;
        abstol_data[global_index++] = eq->abstol;
      }
    }
  }

} // End adsim namespace
