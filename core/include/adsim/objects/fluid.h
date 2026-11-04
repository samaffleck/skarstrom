/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include <functional>
#include <fstream>
#include <cctype>
#include <iostream>

// Eigen includes
#include <Eigen/Core>

// Cereal includes
#include <cereal/types/vector.hpp>
#include <cereal/types/unordered_map.hpp>

#include "adsim/math/expression_param.h"
#include "adsim/objects/viscosity_coefficients.h"
#include "adsim/objects/thermal_conductivity_coefficients.h"
#include <cereal/types/memory.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>

// CoolProp includes
#include "CoolProp.h"


namespace adsim {

  // Contains component data that can not be found in coolprop
  struct ComponentData  {
    // Air as default params
    double eps = 78.6; // = e/k_b [=] [K]
    double sigma = 3.711; // = [Angstrom]
  };

  namespace Components{
    // NOTE: these are `inline` (C++17) so the whole program shares ONE instance.
    // With `static` each translation unit got its own copy, so eps/sigma edits
    // made through the wasm contract layer never reached the copy read by
    // binary_diffusion_model.cpp.
    inline std::vector<std::string> componentNames = {
      "Air", 
      "Ammonia", 
      "Argon", 
      "Benzene",
      "Butane", 
      "CarbonDioxide", 
      "CarbonMonoxide",
      "Cyclohexane",
      "Ethane", 
      "Ethanol",
      "Ethylene", 
      "Fluorine",
      "Helium", 
      "Hexane", 
      "Hydrogen", 
      "HydrogenSulfide", 
      "Isobutane", 
      "Isopentane", 
      "Krypton",
      "Methane", 
      "Methanol", 
      "Neon", 
      "Nitrogen", 
      "NitrousOxide", 
      "Octane", 
      "Oxygen",
      "Pentane", 
      "Propane", 
      "Propylene", 
      "SulfurDioxide", 
      "SulfurHexafluoride", 
      "Toluene", 
      "Water", 
      "Xenon"
    };

    // 1. Poling BE, Prausnitz JM, O’Connell JP. Properties of Gases and Liquids, Fifth Edition. Fifth edition. New York, N.Y: McGraw-Hill Education; 2020. 1 p. (McGraw-Hill’s AccessEngineering).
    //
    // These Lennard-Jones parameters are fixed physical constants, so the table
    // is `const`: it must never be mutated at runtime. Keeping it immutable is
    // what prevents per-run edits from leaking between decks (the table is
    // process-global and is not part of any deck's serialized state).
    inline const std::unordered_map<std::string, ComponentData> componentData = {

      { "Acetone",            { 560.2, 4.600 } },  // CH3COCH3
      { "Acetylene",          { 231.8, 4.033 } },  // C2H2
      { "Air",                {  78.6, 3.711 } },
      { "Ammonia",            { 558.3, 2.900 } },  // NH3
      { "Argon",              {  93.3, 3.542 } },  // Ar
      { "Arsine",             { 259.8, 4.145 } },  // AsH3

      { "Benzene",            { 412.3, 5.349 } },  // C6H6
      { "BoronChloride",      { 337.7, 5.127 } },  // BCl3
      { "BoronFluoride",      { 186.3, 4.198 } },  // BF3
      { "Bromine",            { 507.9, 4.296 } },  // Br2
      { "Butane",             { 531.4, 4.687 } },  // n-C4H10

      { "CarbonDisulfide",    { 467.0, 4.483 } },  // CS2
      { "CarbonDioxide",      { 195.2, 3.941 } },  // CO2
      { "CarbonMonoxide",     {  91.7, 3.690 } },  // CO
      { "CarbonTetrachloride",{ 322.7, 5.947 } },  // CCl4
      { "CarbonTetrafluoride",{ 134.0, 4.662 } },  // CF4
      { "CarbonylSulfide",    { 336.0, 4.130 } },  // COS
      { "Chloroform",         { 340.2, 5.389 } },  // CHCl3
      { "Chlorine",           { 316.0, 4.217 } },  // Cl2
      { "Cyclohexane",        { 297.1, 6.182 } },  // C6H12
      { "Cyclopropane",       { 248.9, 4.807 } },  // C3H6
      { "Cyanogen",           { 348.6, 4.361 } },  // C2N2

      { "DimethylPropane",    { 193.4, 6.464 } },  // 2,2-dimethylpropane, C(CH3)4

      { "Ethane",             { 215.7, 4.443 } },  // C2H6
      { "Ethanol",            { 362.6, 4.530 } },  // C2H5OH
      { "EthylAcetate",       { 521.3, 5.205 } },  // CH3COOC2H5
      { "EthylChloride",      { 300.0, 4.898 } },  // C2H5Cl
      { "EthylEther",         { 313.8, 5.678 } },  // CH3OCH2CH3
      { "Ethylene",           { 224.7, 4.163 } },  // C2H4

      { "Fluorine",           { 112.6, 3.357 } },  // F2

      { "Helium",             { 10.22, 2.551 } }, // He
      { "Hexane",             { 399.3, 5.949 } },  // n-C6H14
      { "Hydrogen",           {  59.7, 2.827 } },  // H2
      { "HydrogenBromide",    { 449.0, 3.353 } },  // HBr
      { "HydrogenChloride",   { 344.7, 3.339 } },  // HCl
      { "HydrogenCyanide",    { 569.1, 3.630 } },  // HCN
      { "HydrogenFluoride",   { 330.0, 3.148 } },  // HF
      { "HydrogenIodide",     { 288.7, 4.211 } },  // HI
      { "HydrogenPeroxide",   { 289.3, 4.196 } },  // H2O2
      { "HydrogenSulfide",    { 301.1, 3.623 } },  // H2S

      { "Iodine",             { 474.2, 5.160 } },  // I2
      { "Isobutane",          { 330.1, 5.278 } },  // iso-C4H10
      { "Isopentane",         { 341.1, 5.784 } },  // (using n-C5H12 parameters)

      { "Krypton",            { 178.9, 3.655 } },  // Kr

      { "Mercury",            { 750.0, 2.969 } },  // Hg
      { "MercuricBromide",    { 686.2, 5.080 } },  // HgBr2
      { "MercuricChloride",   { 750.0, 4.550 } },  // HgCl2
      { "MercuricIodide",     { 695.6, 5.625 } },  // HgI2
      { "Methane",            { 148.6, 3.758 } },  // CH4
      { "Methanol",           { 481.8, 3.626 } },  // CH3OH
      { "MethylAcetate",      { 469.8, 4.936 } },  // CH3COOCH3
      { "MethylAcetylene",    { 251.8, 4.761 } },  // CH3CCH
      { "MethylBorate",       { 396.7, 5.503 } },  // B(OCH3)3
      { "MethylBromide",      { 449.2, 4.118 } },  // CH3Br
      { "MethylChloride",     { 350.0, 4.182 } },  // CH3Cl
      { "MethylEther",        { 395.0, 4.307 } },  // CH3OCH3

      { "Neon",               {  32.8, 2.820 } },  // Ne
      { "NitricOxide",        { 116.7, 3.492 } },  // NO
      { "Nitrogen",           {  71.4, 3.798 } },  // N2
      { "NitrosylChloride",   { 395.3, 4.112 } },  // NOCl
      { "NitrousOxide",       { 232.4, 3.828 } },  // N2O

      { "Octane",             { 320.0, 7.45  } },  // user value (not in table)
      { "Oxygen",             { 106.7, 3.467 } },  // O2

      { "Pentane",            { 341.1, 5.784 } },  // n-C5H12
      { "Phosphine",          { 251.5, 3.981 } },  // PH3
      { "Propane",            { 237.1, 5.118 } },  // C3H8
      { "PropylAlcohol",      { 576.7, 4.549 } },  // n-C3H7OH
      { "Propylene",          { 298.9, 4.678 } },  // CH3CHCH3

      { "SiliconHydride",     { 207.6, 4.084 } },  // SiH4
      { "SiliconTetrafluoride",{171.9, 4.880 } },  // SiF4
      { "StannicBromide",     { 563.7, 6.388 } },  // SnBr4

      { "SulfurDioxide",      { 335.4, 4.112 } },  // SO2
      { "SulfurHexafluoride", { 222.1, 5.128 } },  // SF6

      { "Toluene",            { 285.7, 5.93  } },  // user value (not in table)

      { "UraniumHexafluoride",{ 236.8, 5.967 } },  // UF6

      { "Water",              { 809.1, 2.641 } },  // H2O

      { "Xenon",              { 231.0, 4.047 } }   // Xe
    };

    // Read-only lookup into the constant table. Unknown components fall back to
    // the ComponentData defaults (air), matching the previous operator[]
    // behaviour without inserting into (mutating) the table.
    inline const ComponentData& GetComponentData(const std::string& name) {
      static const ComponentData fallback{};
      auto it = componentData.find(name);
      return it != componentData.end() ? it->second : fallback;
    }

  }

  using Component = std::string;

  enum class AdsorbedHeatCapacityMode {
    CONSTANT,   // Use constant user-defined value from adsorbed_heat_capacity_ [J/(mol-K)]
    GAS_PHASE   // Equal to the component's gas-phase molar heat capacity: Cpa = Cp,i * M_i [J/(mol-K)]
  };

  static const char* AdsorbedHeatCapacityModeList[] = {
    "Constant",
    "Gas phase"
  };

  // Selects which of the two self-consistent forms of the column energy balance
  // to use. Both keep the adsorbed-phase heat capacity in the accumulation term;
  // they differ by a single adsorption source term, paired with the definition of
  // the supplied heat of adsorption:
  //   ISOSTERIC       - heat of adsorption is the isosteric heat q_st; source is
  //                     rho_b * dH * dq/dt (Casas/Mazzotti, Ruthven). Default.
  //   INTERNAL_ENERGY - heat of adsorption is the internal-energy convention used
  //                     by Haghpanah 2013 / Ward & Pini; source additionally carries
  //                     the adsorbed-phase sensible-enthalpy transfer term, i.e.
  //                     rho_b * (dH - Cpa*T) * dq/dt.
  // The two are related by q_st = dH_internal_energy - Cpa*T, so the convention must
  // match the value stored in the adsorbent heat_of_adsorption_.
  //   ISOSTERIC_RT    - the supplied value is -dU, the internal energy of
  //                     adsorption, and the isosteric heat follows the ideal-gas
  //                     relation q_st = -dU + R*T at the local temperature of the
  //                     adsorbed phase; source is rho_b * (-dU + R*T) * dq/dt
  //                     (PEESEgroup / Yancy-Caballero 2020, whose energy balance
  //                     writes the source as -(dU - R*T)).
  enum class HeatOfAdsorptionConvention {
    ISOSTERIC,
    INTERNAL_ENERGY,
    ISOSTERIC_RT
  };

  static const char* HeatOfAdsorptionConventionList[] = {
    "Isosteric",
    "Internal energy",
    "Isosteric from -dU + RT"
  };

  // Source of the pure-component gas heat capacities C_p,i.
  //   COOLPROP - CoolProp at the initial-condition pressure and temperature
  //              (default), held for the whole run.
  //   CONSTANT - the user's molar value per component, gas_heat_capacity_
  //              [J/(mol-K)] (e.g. the constant 30.7 J/(mol-K) of
  //              Yancy-Caballero 2020 and Haghpanah 2013).
  enum class GasHeatCapacityMode {
    COOLPROP,
    CONSTANT
  };

  // Display names (front ends); the registry names are "coolprop,constant".
  static const char* GasHeatCapacityModeList[] = {
    "CoolProp",
    "User defined"
  };

  // Source of the gas mixture viscosity.
  //   POLYNOMIAL - mole-fraction average of the pure-component mu_i(T)
  //                polynomials (default; viscosity_coefficients.h).
  //   CONSTANT   - the user's mixture value viscosity_ [Pa-s], the same in
  //                every cell at every temperature and composition.
  enum class ViscosityMode {
    POLYNOMIAL,
    CONSTANT
  };

  // Display names (front ends); the registry names are "polynomial,constant".
  // The polynomials are fitted to CoolProp's dilute-gas values (Chapman-Enskog
  // for the few components CoolProp has no viscosity model for).
  static const char* ViscosityModeList[] = {
    "CoolProp",
    "User defined"
  };

  class InitialConditions;

  class Fluid {
  public:
    // Seed of a component's constant gas heat capacity [J/(mol-K)]: 7/2 R, the
    // ideal diatomic gas. Only used once the CONSTANT mode is chosen.
    static constexpr double kDefaultGasHeatCapacity = 29.1;

    void AddComponent(const std::string& component_name);
    void RemoveComponent(const std::string& component_name);
    void Initialise(InitialConditions& ic);
    bool PerformChecks(std::string& error_message) const;

    // Refresh the cached value of every ExpressionParam owned by the fluid
    // (currently the per-component adsorbed heat capacities). Part of the
    // UserData::EvaluateAllExpressions cascade so the expression string stays
    // the source of truth. Defined inline (EvalContext is complete here).
    void EvaluateExpressions(const EvalContext& ctx) {
      for (auto& [c, p] : adsorbed_heat_capacity_) p.Evaluate(ctx);
      for (auto& [c, p] : gas_heat_capacity_) p.Evaluate(ctx);
      viscosity_.Evaluate(ctx);
    }

    // Gas mixture viscosity [Pa-s] of the CONSTANT viscosity mode, or a
    // negative value when the viscosity comes from the polynomials.
    double ConstantViscosity() const {
      return viscosity_mode_ == ViscosityMode::CONSTANT ? viscosity_.Get() : -1.0;
    }

    // Pure-component gas-phase viscosity mu_i(T) [Pa-s] from a hard-coded degree-4
    // polynomial per component (see viscosity_coefficients.h — generated offline
    // from CoolProp / Chapman-Enskog, NOT refit at runtime). Calling CoolProp per
    // cell per residual is ~20 us/call and would add minutes to a solve; a
    // polynomial eval is ~1 ns. Gas viscosity is ~pressure-independent over
    // PSA/TSA swings (<~1% across 0.1-10 bar), so the fit is temperature-only.
    static constexpr double kViscTmin = 150.0; // [K] fit/eval window (no extrapolation)
    static constexpr double kViscTmax = 650.0; // [K]

    // Convenience getter (find + eval). For hot loops, resolve the coefficients
    // once with ViscosityCoeffs() and call EvalViscosity() so the per-cell inner
    // loop does no string hashing / map lookups.
    static double GetComponentViscosity(const std::string& component, double T);
    static const std::array<double, 5>* ViscosityCoeffs(const std::string& component) {
      auto it = kViscosityCoefficients.find(component);
      return it == kViscosityCoefficients.end() ? nullptr : &it->second;
    }
    static double EvalViscosity(const std::array<double, 5>& k, double T) {
      const double Tc = T < kViscTmin ? kViscTmin : (T > kViscTmax ? kViscTmax : T);
      double mu = k[4];                                     // Horner, rising powers
      for (int i = 3; i >= 0; --i) mu = mu * Tc + k[i];
      return mu > 0.0 ? mu : 1.0e-6;                        // guard fit undershoot at edges
    }

    // Pure-component gas-phase thermal conductivity k_i(T) [W/(m.K)], same
    // hard-coded-polynomial scheme as viscosity (see thermal_conductivity_
    // coefficients.h — CoolProp fit, or Eucken relation where CoolProp lacks a
    // conductivity model). Same [kViscTmin, kViscTmax] evaluation window.
    static double GetComponentThermalConductivity(const std::string& component, double T);
    static const std::array<double, 5>* ThermalConductivityCoeffs(const std::string& component) {
      auto it = kThermalConductivityCoefficients.find(component);
      return it == kThermalConductivityCoefficients.end() ? nullptr : &it->second;
    }
    static double EvalThermalConductivity(const std::array<double, 5>& k, double T) {
      const double Tc = T < kViscTmin ? kViscTmin : (T > kViscTmax ? kViscTmax : T);
      double kg = k[4];                                     // Horner, rising powers
      for (int i = 3; i >= 0; --i) kg = kg * Tc + k[i];
      return kg > 0.0 ? kg : 1.0e-4;                        // guard fit undershoot at edges
    }

    // Gas-phase heat capacity of a pure component on a molar basis,
    // C_p,i * M_i [J/(mol.K)], from the cached per-kilogram heat_capacity_.
    // Requires Initialise().
    double MolarHeatCapacity(const std::string& component) const {
      return heat_capacity_.at(component) * molar_mass_.at(component);
    }

    // Molar heat capacity Cpa of `component` in the adsorbed phase [J/(mol.K)]:
    // the gas-phase molar value (GAS_PHASE) or the user's constant (CONSTANT).
    double AdsorbedPhaseHeatCapacity(const std::string& component) const {
      return adsorbed_heat_capacity_mode_ == AdsorbedHeatCapacityMode::GAS_PHASE
        ? MolarHeatCapacity(component)
        : adsorbed_heat_capacity_.at(component).Get();
    }

    // Mass-specific heat capacity of a gas mixture [J/(kg.K)]:
    //   C_p,g = sum_i y_i M_i C_p,i / sum_i y_i M_i,
    // the molar mixture heat capacity divided by the mixture molar mass. Mole
    // fractions mix MOLAR heat capacities; mixing the per-kilogram C_p,i by mole
    // fraction is wrong unless every M_i is equal. Takes the two mole-fraction
    // weighted sums, sum_i y_i C_p,i M_i and sum_i y_i M_i, so the solver's hot
    // loop (UpdateHeatCapacity) and MixtureHeatCapacity below share one rule.
    // Returns 0 when the molar-mass sum is not positive (no gas).
    static double MixtureHeatCapacityFromSums(double sum_y_molar_cp, double sum_y_molar_mass) {
      return sum_y_molar_mass > 0.0 ? sum_y_molar_cp / sum_y_molar_mass : 0.0;
    }

    // C_p,g [J/(kg.K)] of the mixture with mole fractions y (component -> y_i;
    // components absent from y count as 0). Requires Initialise().
    double MixtureHeatCapacity(const std::unordered_map<std::string, double>& y) const {
      double sum_cp = 0.0, sum_m = 0.0;
      for (const auto& [component, yi] : y) {
        sum_cp += yi * MolarHeatCapacity(component);
        sum_m += yi * molar_mass_.at(component);
      }
      return MixtureHeatCapacityFromSums(sum_cp, sum_m);
    }

    std::vector<std::string> components_{};
    std::unordered_map<std::string, double> molar_mass_{};
    std::unordered_map<std::string, ExpressionParam> adsorbed_heat_capacity_{}; // [J/(mol-K)] - used when mode is CONSTANT
    std::unordered_map<std::string, double> heat_capacity_{}; // [J/(kg-K)] - pure-component, cached by Initialise (CoolProp or gas_heat_capacity_); mix with MixtureHeatCapacity
    AdsorbedHeatCapacityMode adsorbed_heat_capacity_mode_ = AdsorbedHeatCapacityMode::CONSTANT;
    HeatOfAdsorptionConvention heat_of_adsorption_convention_ = HeatOfAdsorptionConvention::ISOSTERIC;
    GasHeatCapacityMode gas_heat_capacity_mode_ = GasHeatCapacityMode::COOLPROP;
    std::unordered_map<std::string, ExpressionParam> gas_heat_capacity_{}; // [J/(mol-K)] - used when gas_heat_capacity_mode_ is CONSTANT
    ViscosityMode viscosity_mode_ = ViscosityMode::POLYNOMIAL;
    ExpressionParam viscosity_{1.8e-5}; // [Pa-s] - mixture value, used when viscosity_mode_ is CONSTANT


  private:
    bool IsComponentInList(const std::string& component_name) const;
    bool IsComponentValid(const std::string& component_name) const;

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::make_nvp("Components", components_),
        cereal::make_nvp("Adsorbed_Heat_Capacity", adsorbed_heat_capacity_)
      );
      if (version >= 2) {
        ar(cereal::make_nvp("Adsorbed_Heat_Capacity_Mode", adsorbed_heat_capacity_mode_));
      }
      // Legacy v1 files carry no mode field; the CONSTANT default is retained.
      if (version >= 3) {
        ar(cereal::make_nvp("Heat_Of_Adsorption_Convention", heat_of_adsorption_convention_));
      }
      // Legacy v1/v2 files carry no convention field; the ISOSTERIC default is retained.
      if (version >= 4) {
        ar(cereal::make_nvp("Gas_Heat_Capacity_Mode", gas_heat_capacity_mode_),
           cereal::make_nvp("Gas_Heat_Capacity", gas_heat_capacity_),
           cereal::make_nvp("Viscosity_Mode", viscosity_mode_),
           cereal::make_nvp("Viscosity", viscosity_));
      }
      // Legacy v1-v3 files: CoolProp heat capacities and polynomial viscosities.
      // On load, seed an editable constant heat capacity for every component
      // without one (as AddComponentToDeck does), so the map lists them all.
      if constexpr (Archive::is_loading::value) {
        for (const auto& c : components_) gas_heat_capacity_.try_emplace(c, kDefaultGasHeatCapacity);
      }
    }
  };
  
  inline static double GetHeatCapacity(const std::string& component, double P, double T) {
    double cp = 0.0; // failure sentinel, consistent with the other CoolProp getters below
    try {
      cp = CoolProp::PropsSI("C", "P", P, "T", T, component); // [J/(kg-K)]
    }
    catch (const CoolProp::CoolPropBaseError& e) {
      std::cerr << "CoolProp Error: " << e.what() << std::endl;
    }

    return cp;
  }
  
  inline static double GetViscosity(const std::string& component, double P, double T) {
    double vis = 0.0;
    try {
      vis = CoolProp::PropsSI("V", "P", P, "T", T, component); 
    }
    catch (const CoolProp::CoolPropBaseError& e) {
      std::cerr << "CoolProp Error: " << e.what() << std::endl;
    }

    return vis;
  }

  inline static double GetThermalConductivity(const std::string& component, double P, double T) {
    double k = 0.0;
    try {
      k = CoolProp::PropsSI("CONDUCTIVITY", "P", P, "T", T, component); 
    }
    catch (const CoolProp::CoolPropBaseError& e) {
      std::cerr << "CoolProp Error: " << e.what() << std::endl;
    }

    return k;
  }

  inline static double GetMolarMass(const std::string& component) {
    double mr = 0.0;
    try {
      mr = CoolProp::Props1SI(component, "molar_mass");
    }
    catch (const CoolProp::CoolPropBaseError& e) {
      std::cerr << "CoolProp Error: " << e.what() << std::endl;
    }

    return mr;
  }

} // end adsim namespace
