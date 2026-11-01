/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <vector>
#include <functional>
#include <unordered_map>
#include <string>
#include <optional>

// SUNDIALS includes
#include <ida/ida.h>

// AdSim includes
#include "adsim/solver/initial_conditions.h"
#include "adsim/math/constants.h"
#include "adsim/equations/adsorbent_equations.h"
#include "adsim/objects/domain.h"
#include "adsim/objects/reactor_equations.h"


namespace adsim {

  const double R = adsim::constants::kIdealGasConstant;

  enum class EquationType {
    Differential,
    Algebraic
  };

  enum class CellType {
    CellCentre,
    CellFace
  };

  enum class EquationDomain {
    AllCells,         // packed at every cell including voids (default)
    PorousMediaOnly,  // packed only where adsorbent is present; skipped at void cells
  };

  // Forward decleration
  class UserData;
  class Step;

  using vec = std::vector<double>;
  using update_fun = std::function<void(UserData* udata)>;

  class IndependentEquation {
  public:
    IndependentEquation(
      std::string name_,
      CellType cell_type_,
      update_fun update_func_
    ) : name(name_), cell_type(cell_type_), update_func(update_func_) {}

    std::string name;
    vec var;  // Vector of values
    
    void Resize(size_t n_cells);
    void Update(UserData* udata);
    size_t Size() const;

  private:
    CellType cell_type;
    update_fun update_func;
  };

  class Equation {
  public:
    Equation(
      std::string name_, 
      CellType cell_type_, 
      EquationType equation_type_,
      update_fun update_func_, 
      double abstol_ = 1e-3
    ) : name(name_), cell_type(cell_type_), equation_type(equation_type_), abstol(abstol_), update_func(update_func_) {}

    std::string name;
    double abstol;  // Absolute tolerance for this equation type
    vec var;      // Vector
    vec ddt_var;  // Time derivative of vector
    vec res_var;  // Residual vector

    void SetInitialConditions(double x);
    void SetSmoothInitialConditions(double x_inlet, double x_outlet);
    void Resize(size_t n_cells);
    void SetID(size_t& start_index, sunrealtype* id);
    double GetIDValue(size_t cell_idx) const; // Returns ID value (0.0 or 1.0) for a specific cell
    void Pack(size_t& start_index, sunrealtype* y_data, sunrealtype* yp_data) const;
    void PackResiduals(size_t& start_index, sunrealtype* res_data) const;
    void Unpack(size_t& start_index, const sunrealtype* y_data, const sunrealtype* yp_data);
    void Update(UserData* udata);
    void SetEquationType(EquationType type) { equation_type = type; }
    EquationType GetEquationType() const { return equation_type; }
    void SetDomain(EquationDomain d) { domain = d; }
    EquationDomain GetDomain() const { return domain; }
    // Optional per-equation active-cell mask. When set, this equation only
    // contributes a state-vector entry where (*active_mask)[cell] is true (in
    // addition to the PorousMediaOnly void rule). Used to restrict a
    // per-adsorbent loading equation to the cells of the layers that contain
    // that adsorbent. The pointee must outlive the equation (owned by Equations).
    void SetActiveMask(const std::vector<bool>* m) { active_mask = m; }
    // Returns true if this equation contributes a state-vector entry at cell_idx.
    // PorousMediaOnly equations are inactive at cells flagged by is_void_mask;
    // equations with an active_mask are inactive where the mask is false.
    bool IsActiveAt(size_t cell_idx, const std::vector<bool>* is_void_mask) const;
    size_t Size() const;

  private:
    CellType cell_type;
    EquationType equation_type;
    EquationDomain domain = EquationDomain::AllCells;
    const std::vector<bool>* active_mask = nullptr;
    update_fun update_func;

  };

  using ComponentEquation = std::unordered_map<Component, Equation>;
  using ComponentIndependentEquation = std::unordered_map<Component, IndependentEquation>;
  // Per-adsorbent nesting: AdsorbentID -> (component -> equation). Adsorbed-phase
  // quantities (loading, equilibrium loading, mass-transfer coefficient) are kept
  // separate per adsorbent so each adsorbent retains its own isotherm and kinetics.
  using AdsorbentComponentEquation = std::unordered_map<AdsorbentID, ComponentEquation>;
  using AdsorbentComponentIndependentEquation = std::unordered_map<AdsorbentID, ComponentIndependentEquation>;

  struct AdsorbentEnergyTerm {
    double heat_capacity;
    const vec* loading;
    const vec* loading_rate;
    const vec* heat_of_adsorption;
  };

  class Equations {
  public:    
    // Must be called before Initialise (UserData::Initialise does, from the
    // IDASettings flags of the same name). See IDASettings.
    void SetReducedUnknowns(bool alias_lte_solid_temperature, bool derived_concentration) {
      alias_lte_solid_temperature_ = alias_lte_solid_temperature;
      derived_concentration_ = derived_concentration;
    }
    // Must be called before Initialise. See IDASettings::pressure_differential_.
    void SetPressureDifferential(bool on) { pressure_differential_ = on; }
    // True when Ts is not a DAE unknown but a copy of Tf (LTE, non-isothermal,
    // flag on) / when C is not a DAE unknown but P / (R T) (flag on).
    bool SolidTemperatureAliased() const { return solid_temperature_aliased_; }
    bool ConcentrationDerived() const { return derived_concentration_; }
    void Initialise(const Fluid& fluid, const Bed& bed, const AdsorbentDirectory& dir);
    void UpdateDifferentialEquations(UserData* udata);
    // Fields the residual depends on (called on every residual evaluation).
    void UpdateIndependentEquations(UserData* udata);
    // Logging-only fields (flow rates, source totals, pressure gradient, ...);
    // nothing in the residual reads them. Call after UpdateIndependentEquations.
    void UpdateDiagnosticEquations(UserData* udata);
    void UpdateBoundaryConditions(UserData* udata, double time);
    void SetCurrentStep(Step* step) { current_step_ = step; }
    void Resize(size_t n_cells);
    int GetSystemSize();
    void SetInitialConditions(InitialConditions& ic, ThermalEquilibriumMode mode);
    // Isothermal runs (ReactorEquations::energy_balance_mode_ == ISOTHERMAL):
    // Initialise leaves the fluid / solid / wall temperature equations out of
    // the DAE, and this pins their values to T (time derivative 0) for the
    // whole run. Call after SetInitialConditions. No-op when not isothermal.
    void HoldIsothermalTemperature(double T);
    bool IsIsothermal() const { return isothermal_; }
    // The isothermal temperature, or nullopt for a non-isothermal run.
    std::optional<double> IsothermalTemperature() const {
      return isothermal_ ? std::optional<double>(isothermal_temperature_) : std::nullopt;
    }
    void SetAdsorbedPhaseInitialConditions();
    void SetID(sunrealtype* id);
    void Pack(sunrealtype* y_data, sunrealtype* yp_data) const;
    void PackResiduals(sunrealtype* res_data) const;
    void Unpack(const sunrealtype* y_data, const sunrealtype* yp_data);
    void SetIDBanded(sunrealtype* id);
    void PackBanded(sunrealtype* y_data, sunrealtype* yp_data) const;
    void PackResidualsBanded(sunrealtype* res_data) const;
    void UnpackBanded(const sunrealtype* y_data, const sunrealtype* yp_data);
    void SetToleranceVector(sunrealtype* abstol_data) const;
    void SetToleranceVectorBanded(sunrealtype* abstol_data) const;
    // "<equation name> @ cell <i>" for entry g of the banded state vector.
    std::string DescribeBandedIndex(size_t g) const;
    const std::vector<Equation*>& GetEquations() const { return equations_; }
    const std::vector<IndependentEquation*>& GetIndependentEquations() const { return independent_equations_; }
    // Adsorbents that appear in at least one layer composition, in ascending ID
    // order (deterministic for banded layout / logging). Populated in Initialise().
    const std::vector<AdsorbentID>& GetAdsorbentIDs() const { return adsorbent_ids_; }
    const std::vector<const vec*>& GetLoadingRateTerms(AdsorbentID id) const {
      return loading_rate_terms_.at(adsorbent_slot_.at(id));
    }
    const std::vector<AdsorbentEnergyTerm>& GetEnergyTerms(AdsorbentID id) const {
      return energy_terms_.at(adsorbent_slot_.at(id));
    }
    // Per-adsorbent active-cell mask (true where the adsorbent is present),
    // or nullptr if the adsorbent is unknown. Stable after Initialise().
    const std::vector<bool>* GetAdsorbentMask(AdsorbentID id) const {
      auto it = adsorbent_cell_mask_.find(id);
      return it == adsorbent_cell_mask_.end() ? nullptr : &it->second;
    }
            
    // Coupled Differential Algebraic Equations
    Equation concentration{ "Concentration", CellType::CellCentre, EquationType::Differential, UpdateConcentration, 1e-8 };
    Equation pressure{ "Pressure", CellType::CellCentre, EquationType::Algebraic, UpdatePressure, 1e-6 };
    Equation fluid_temperature{ "Fluid Temperature", CellType::CellCentre, EquationType::Differential, UpdateTemperature, 1e-4 };
    Equation solid_temperature{ "Solid Temperature", CellType::CellCentre, EquationType::Differential, UpdateSolidTemperature, 1e-4 };
    Equation wall_temperature{ "Wall Temperature", CellType::CellCentre, EquationType::Differential, UpdateWallTemperature, 1e-4 };
    Equation velocity{ "Velocity", CellType::CellFace, EquationType::Algebraic, UpdateVelocity, 1e-6 };
    ComponentEquation mole_fraction;
    AdsorbentComponentEquation ads_concentration; // Actual adsorbed conc (rate equation), per adsorbent

    // Independent equations - not used in solver
    IndependentEquation total_mole_fraction{ "Total Mole Fraction", CellType::CellCentre, UpdateTotalMoleFraction };
    IndependentEquation fluid_density{ "Fluid Density", CellType::CellCentre, UpdateFluidDensity };
    IndependentEquation viscosity{ "Viscosity", CellType::CellCentre, UpdateViscosity };
    IndependentEquation heat_capacity{ "Heat Capacity", CellType::CellCentre, UpdateHeatCapacity };
    IndependentEquation sherwood_number{ "Sherwood Number", CellType::CellCentre, UpdateSherwoodNumber };
    IndependentEquation reynolds_number{ "Reynolds Number", CellType::CellCentre, UpdateReynoldsNumber };
    IndependentEquation schmidt_number{ "Schmidt Number", CellType::CellCentre, UpdateSchmidtNumber };
    IndependentEquation fluid_thermal_conductivity{ "Fluid Thermal Conductivity", CellType::CellCentre, UpdateFluidThermalConductivity };
    IndependentEquation axial_thermal_conductivity{ "Axial Thermal Conductivity", CellType::CellCentre, UpdateAxialThermalConductivity };
    IndependentEquation internal_heat_transfer_coefficient{ "Internal Heat Transfer Coefficient", CellType::CellCentre, UpdateInternalHeatTransferCoefficient };
    IndependentEquation external_heat_transfer_coefficient{ "External Heat Transfer Coefficient", CellType::CellCentre, UpdateExternalHeatTransferCoefficient };
    IndependentEquation particle_heat_transfer_coefficient{ "Particle Heat Transfer Coefficient", CellType::CellCentre, UpdateParticleHeatTransferCoefficient };
    IndependentEquation molar_mass{ "Molar Mass", CellType::CellCentre, UpdateMolarMass };
    IndependentEquation volumetric_flow_rate{ "Volumetric Flow Rate", CellType::CellFace, UpdateVolumetricFlowRate };
    IndependentEquation molar_flow_rate{ "Molar Flow Rate", CellType::CellFace, UpdateMolarFlowRate };
    IndependentEquation mass_source{ "Mass Source", CellType::CellCentre, UpdateMassSource };
    IndependentEquation energy_source{ "Energy Source", CellType::CellCentre, UpdateEnergySource };
    IndependentEquation approach_to_fluidisation{ "Approach To Fluidisation", CellType::CellFace, UpdateApproachToFluidisation };
    IndependentEquation pressure_gradient{ "Pressure Gradient", CellType::CellCentre, UpdatePressureGradient };
    IndependentEquation axial_dispersion_coefficient{ "Axial Dispersion Coefficient", CellType::CellCentre, UpdateAxialDispersionCoefficient };
    ComponentIndependentEquation molar_concentration;
    AdsorbentComponentIndependentEquation mass_transfer_coefficient; // per adsorbent
    ComponentIndependentEquation molecular_diffusion;
    AdsorbentComponentIndependentEquation eq_ads_concentration; // Equilibrium adsorbed conc (from isotherm), per adsorbent
    ComponentIndependentEquation component_mass_source;
    ComponentIndependentEquation component_energy_source;
    ComponentIndependentEquation mole_fraction_face;
    
    // Domain vector
    DomainVectors domain_vectors;

    // Component-independent terms of the mole-fraction residual, filled once
    // per residual evaluation by the first component's UpdateMoleFraction and
    // read by every component's cell loop (see there). Array-of-structs so the
    // cell loop carries two base pointers, not one per term.
    struct MoleFractionFaceTerms {   // face f sits between cells f and f+1
      double P;          // 0.5 * (P[f] + P[f+1])
      double T;          // DenomFloor(0.5 * (T[f] + T[f+1]))
      double P_over_T;   // P / T
      double d;          // harmonic-mean axial dispersion across the face
    };
    struct MoleFractionCellTerms {
      double Pc, Tc;              // DenomFloor'ed cell P and T
      double dx_et;               // dx * et
      double D_e_coef, D_w_coef;  // (eb * Tc * d_e) / (Pc * dx), west likewise
      double Tc_over_Pc;          // Tc / Pc
      double R_Tc_dx;             // R * Tc * dx
    };
    std::vector<MoleFractionFaceTerms> mole_fraction_face_terms;
    std::vector<MoleFractionCellTerms> mole_fraction_cell_terms;

  private:
    std::vector<Equation*> equations_;
    // All independent fields in registration order (what gets logged), split
    // into the residual-required ones and the logging-only diagnostics.
    std::vector<IndependentEquation*> independent_equations_;
    std::vector<IndependentEquation*> residual_independent_equations_;
    std::vector<IndependentEquation*> diagnostic_equations_;
    // Adsorbents present in the bed (ascending ID) and their per-cell presence
    // masks. Built once in Initialise(); pointers into adsorbent_cell_mask_ are
    // handed to the per-adsorbent loading equations via SetActiveMask().
    std::vector<AdsorbentID> adsorbent_ids_;
    std::unordered_map<AdsorbentID, std::vector<bool>> adsorbent_cell_mask_;
    // Stable pointer tables for residual loops. They point to vector objects,
    // whose storage can resize without invalidating the table.
    std::unordered_map<AdsorbentID, size_t> adsorbent_slot_;
    std::vector<std::vector<const vec*>> loading_rate_terms_;
    std::vector<std::vector<AdsorbentEnergyTerm>> energy_terms_;
    Step* current_step_ = nullptr;
    // Isothermal run: temperature equations not registered in equations_.
    bool isothermal_ = false;
    double isothermal_temperature_ = 0.0;
    // Reduced unknowns (see SetReducedUnknowns). solid_temperature_aliased_ is
    // the effective value: the flag, and LTE, and not isothermal.
    bool alias_lte_solid_temperature_ = true;
    bool derived_concentration_ = true;
    bool solid_temperature_aliased_ = false;
    bool pressure_differential_ = true;
    // Fills the unregistered fields from the unpacked state: Ts (and dTs/dt)
    // from Tf when aliased, C (and dC/dt) from P and T when derived.
    void SyncDerivedFields();

    // Banded state layout (cell-major, active entries only), resolved once so
    // the per-residual pack/unpack is a straight copy. Built on first use;
    // cleared whenever Initialise/Resize can change the layout.
    struct BandedSlot { double* var; double* ddt_var; const double* res_var; const Equation* eq; size_t cell; };
    mutable std::vector<BandedSlot> banded_slots_;
    const std::vector<BandedSlot>& BandedSlots() const;

  };

} // end adsim namespace
