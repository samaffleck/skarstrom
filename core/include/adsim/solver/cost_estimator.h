#pragma once

// STL includes
#include <cmath>
#include <map>
#include <numbers>
#include <unordered_map>
#include <vector>

// Cereal includes
#include <cereal/types/vector.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/cereal.hpp>

// Adsim includes
#include "adsim/objects/reactor.h"
#include "adsim/objects/wall.h"
#include "adsim/solver/kpis.h"


namespace adsim {

  // ===========================================================================
  // Pipe nominal diameter (DN) — used by both pipework and valve costing
  // Inside diameter and wall thickness: ASME B36.19 Schedule 10S (small bore)
  //                                     ASME B36.10 Schedule 40 (DN40 and above)
  // ===========================================================================

  enum PipeDN {
    DN6,  DN8,  DN10, DN15,  DN20,  DN25,  DN32,
    DN40, DN50, DN65, DN80,  DN100, DN125, DN150,
    DN200, DN250, DN300
  };

  // Inside diameter [mm]
  inline const std::map<PipeDN, double> kInsideDiameter_DN = {
    {DN6,    6.2}, {DN8,    8.8}, {DN10,  12.5}, {DN15,  15.8},
    {DN20,  20.9}, {DN25,  26.6}, {DN32,  34.8}, {DN40,  41.9},
    {DN50,  52.5}, {DN65,  66.9}, {DN80,  82.5}, {DN100, 102.3},
    {DN125, 128.2},{DN150, 154.1},{DN200, 202.7}, {DN250, 254.5},
    {DN300, 304.8}
  };

  // Wall thickness [mm]
  inline const std::map<PipeDN, double> kWallThickness_DN = {
    {DN6,  2.00}, {DN8,  2.35}, {DN10, 2.35}, {DN15, 2.77},
    {DN20, 2.77}, {DN25, 3.05}, {DN32, 3.38}, {DN40, 3.56},
    {DN50, 3.91}, {DN65, 4.19}, {DN80, 5.49}, {DN100, 6.02},
    {DN125, 6.35},{DN150, 6.35},{DN200, 6.35}, {DN250, 6.35},
    {DN300, 7.11}
  };

  // ===========================================================================
  // Valve costs
  // Base costs (SS316L, PN16) indexed by (ValveType, PipeDN).
  // Material and pressure rating adjustments applied via factor tables.
  // Reference: Sinnott & Towler; Matches & Sümer industrial cost data (~2024).
  // ===========================================================================

  enum class ValveType {
    Ball,
    Globe,
    Butterfly,
    CheckSwing,
    CheckWafer,
    NeedleAngle
  };

  inline constexpr const char* kValveTypeList[] = {
    "Ball",
    "Globe",
    "Butterfly",
    "Check (swing)",
    "Check (wafer)",
    "Needle/Angle"
  };

  enum class ValvePressureRating {
    PN6, PN10, PN16, PN25, PN40
  };

  inline constexpr const char* kValvePressureRatingList[] = {
    "PN6", "PN10", "PN16", "PN25", "PN40"
  };

  // Base purchase cost [£] per valve — SS316L, PN16, flanged end
  // Rows: ValveType index; Cols: PipeDN index (DN6..DN300)
  // Note: butterfly valves only practical from DN50 upward
  struct ValveCostKey {
    ValveType type;
    PipeDN    dn;
    bool operator==(const ValveCostKey& o) const {
      return type == o.type && dn == o.dn;
    }
  };

  struct ValveCostKeyHash {
    std::size_t operator()(const ValveCostKey& k) const {
      return std::hash<int>{}(static_cast<int>(k.type)) ^ (std::hash<int>{}(static_cast<int>(k.dn)) << 8);
    }
  };

  // clang-format off
  inline const std::unordered_map<ValveCostKey, double, ValveCostKeyHash> kValveBaseCost = {
  //  {type,                 DN}     cost [£]
    {{ValveType::Ball,       DN6},    120},
    {{ValveType::Ball,       DN8},    140},
    {{ValveType::Ball,       DN10},   160},
    {{ValveType::Ball,       DN15},   180},
    {{ValveType::Ball,       DN20},   220},
    {{ValveType::Ball,       DN25},   280},
    {{ValveType::Ball,       DN32},   360},
    {{ValveType::Ball,       DN40},   450},
    {{ValveType::Ball,       DN50},   600},
    {{ValveType::Ball,       DN65},   900},
    {{ValveType::Ball,       DN80},  1200},
    {{ValveType::Ball,       DN100}, 1800},
    {{ValveType::Ball,       DN125}, 2800},
    {{ValveType::Ball,       DN150}, 3800},
    {{ValveType::Ball,       DN200}, 6500},
    {{ValveType::Ball,       DN250},10000},
    {{ValveType::Ball,       DN300},14000},

    {{ValveType::Globe,      DN6},    220},
    {{ValveType::Globe,      DN8},    260},
    {{ValveType::Globe,      DN10},   300},
    {{ValveType::Globe,      DN15},   360},
    {{ValveType::Globe,      DN20},   460},
    {{ValveType::Globe,      DN25},   600},
    {{ValveType::Globe,      DN32},   800},
    {{ValveType::Globe,      DN40},  1000},
    {{ValveType::Globe,      DN50},  1400},
    {{ValveType::Globe,      DN65},  2000},
    {{ValveType::Globe,      DN80},  2800},
    {{ValveType::Globe,      DN100}, 4200},
    {{ValveType::Globe,      DN125}, 6500},
    {{ValveType::Globe,      DN150}, 9000},
    {{ValveType::Globe,      DN200},15000},
    {{ValveType::Globe,      DN250},24000},
    {{ValveType::Globe,      DN300},34000},

    {{ValveType::Butterfly,  DN50},   300},
    {{ValveType::Butterfly,  DN65},   420},
    {{ValveType::Butterfly,  DN80},   560},
    {{ValveType::Butterfly,  DN100},  800},
    {{ValveType::Butterfly,  DN125}, 1100},
    {{ValveType::Butterfly,  DN150}, 1500},
    {{ValveType::Butterfly,  DN200}, 2500},
    {{ValveType::Butterfly,  DN250}, 4000},
    {{ValveType::Butterfly,  DN300}, 5500},

    {{ValveType::CheckSwing, DN15},   180},
    {{ValveType::CheckSwing, DN20},   220},
    {{ValveType::CheckSwing, DN25},   280},
    {{ValveType::CheckSwing, DN32},   360},
    {{ValveType::CheckSwing, DN40},   460},
    {{ValveType::CheckSwing, DN50},   620},
    {{ValveType::CheckSwing, DN65},   920},
    {{ValveType::CheckSwing, DN80},  1250},
    {{ValveType::CheckSwing, DN100}, 1900},
    {{ValveType::CheckSwing, DN150}, 4000},
    {{ValveType::CheckSwing, DN200}, 7000},

    {{ValveType::CheckWafer, DN50},   200},
    {{ValveType::CheckWafer, DN65},   280},
    {{ValveType::CheckWafer, DN80},   360},
    {{ValveType::CheckWafer, DN100},  500},
    {{ValveType::CheckWafer, DN150},  900},
    {{ValveType::CheckWafer, DN200}, 1500},

    {{ValveType::NeedleAngle, DN6},   180},
    {{ValveType::NeedleAngle, DN8},   210},
    {{ValveType::NeedleAngle, DN10},  250},
    {{ValveType::NeedleAngle, DN15},  320},
    {{ValveType::NeedleAngle, DN20},  420},
    {{ValveType::NeedleAngle, DN25},  560},
  };
  // clang-format on

  // Material cost factor relative to SS316L baseline
  inline const std::unordered_map<WallMaterial, double> kValveMaterialFactor = {
    {WallMaterial::SS316L,            1.00},
    {WallMaterial::SS304,             0.80},
    {WallMaterial::SS304L,            0.85},
    {WallMaterial::CarbonSteel_A106,  0.45},
    {WallMaterial::CarbonSteel_A53,   0.42},
    {WallMaterial::Hastelloy_C276,    4.50},
    {WallMaterial::Inconel_625,       4.00},
    {WallMaterial::Aluminium_6061,    0.70},
    {WallMaterial::Titanium_Gr2,      3.00},
    {WallMaterial::Custom,            1.00}
  };

  // Pressure rating cost factor relative to PN16 baseline
  inline const std::unordered_map<ValvePressureRating, double> kValvePressureFactor = {
    {ValvePressureRating::PN6,  0.80},
    {ValvePressureRating::PN10, 0.90},
    {ValvePressureRating::PN16, 1.00},
    {ValvePressureRating::PN25, 1.25},
    {ValvePressureRating::PN40, 1.60}
  };

  struct ValveInstance {
    ValveType           type     = ValveType::Ball;
    PipeDN              dn       = PipeDN::DN25;
    WallMaterial        material = WallMaterial::SS316L;
    ValvePressureRating rating   = ValvePressureRating::PN16;
    int                 count    = 1;

    double CalculateCost() const {
      auto it = kValveBaseCost.find({type, dn});
      if (it == kValveBaseCost.end()) return 0.0;
      const double base      = it->second;
      const double mat_f     = kValveMaterialFactor.at(material);
      const double press_f   = kValvePressureFactor.at(rating);
      return base * mat_f * press_f * count;
    }

    template<class Archive>
    void serialize(Archive& ar) {
      ar(cereal::make_nvp("type", type),
         cereal::make_nvp("dn", dn),
         cereal::make_nvp("material", material),
         cereal::make_nvp("rating", rating),
         cereal::make_nvp("count", count));
    }
  };

  class ValveCost {
  public:
    std::vector<ValveInstance> valves;

    double CalculateCost() const {
      double cost = 0.0;
      for (const auto& v : valves) cost += v.CalculateCost();
      return cost;
    }

  private:
    friend class cereal::access;
    template<class Archive>
    void serialize(Archive& ar) {
      ar(cereal::make_nvp("valves", valves));
    }
  };

  // ===========================================================================
  // Pipework costs
  // ===========================================================================

  // Fabricated pipe cost [£/kg] — approximate market rates for process piping
  inline const std::unordered_map<WallMaterial, double> kPipeMaterialCost = {
    {WallMaterial::SS316L,            8.0},
    {WallMaterial::SS304,             5.5},
    {WallMaterial::SS304L,            6.0},
    {WallMaterial::CarbonSteel_A106,  1.5},
    {WallMaterial::CarbonSteel_A53,   1.4},
    {WallMaterial::Hastelloy_C276,   65.0},
    {WallMaterial::Inconel_625,      55.0},
    {WallMaterial::Aluminium_6061,    5.0},
    {WallMaterial::Titanium_Gr2,     40.0},
    {WallMaterial::Custom,            8.0}
  };

  class PipeSegmentCost {
  public:
    double       length   = 1.0;             // [m]
    PipeDN       pipe_dn  = PipeDN::DN25;
    WallMaterial material = WallMaterial::SS316L;

    double CalculateCost() const {
      const double d_i = kInsideDiameter_DN.at(pipe_dn) * 1e-3; // [m]
      const double t_w = kWallThickness_DN.at(pipe_dn)  * 1e-3; // [m]
      const double d_o = d_i + 2.0 * t_w;
      const double cross_area = std::numbers::pi / 4.0 * (d_o * d_o - d_i * d_i); // [m²]
      const double volume     = length * cross_area;                                // [m³]
      const double rho        = GetWallMaterialProperties(material).density;
      const double cost_per_kg = kPipeMaterialCost.at(material);
      return volume * rho * cost_per_kg; // [£]
    }

    template<class Archive>
    void serialize(Archive& ar) {
      ar(cereal::make_nvp("length", length),
         cereal::make_nvp("pipe_dn", pipe_dn),
         cereal::make_nvp("material", material));
    }
  };

  class PipeworkCost {
  public:
    double CalculatePipeworkCost() const {
      double cost = 0.0;
      for (const auto& segment : pipe_segments) {
        cost += segment.CalculateCost();
      }
      return cost;
    }

    std::vector<PipeSegmentCost> pipe_segments;

  private:
    friend class cereal::access;
    template<class Archive>
    void serialize(Archive& ar) {
      ar(cereal::make_nvp("pipe_segments", pipe_segments));
    }
  };

  // ===========================================================================
  // Vessel cost
  // Sinnott & Towler "Chemical Engineering Design" — mass-based correlation
  // ===========================================================================

  // Fabricated vessel material cost [£/kg] — includes forming + machining allowance
  inline const std::unordered_map<WallMaterial, double> kVesselMaterialCost = {
    {WallMaterial::SS316L,            15.0},
    {WallMaterial::SS304,             10.0},
    {WallMaterial::SS304L,            11.0},
    {WallMaterial::CarbonSteel_A106,   3.5},
    {WallMaterial::CarbonSteel_A53,    3.2},
    {WallMaterial::Hastelloy_C276,    90.0},
    {WallMaterial::Inconel_625,       80.0},
    {WallMaterial::Aluminium_6061,     8.0},
    {WallMaterial::Titanium_Gr2,      60.0},
    {WallMaterial::Custom,            15.0}
  };

  class VesselCost {
  public:
    double installation_factor = 4.0;

    double CalculateCost(const Bed& bed) const {
      const double D   = bed.GetHydraulicDiameter(); // [m]
      const double t_w = bed.wall_.thickness_;        // [m]
      const WallMaterial mat = bed.wall_.material;
      const double rho = GetWallMaterialProperties(mat).density; // [kg/m³]

      // Total column length = sum of layer lengths
      double L = 0.0;
      for (const auto& layer : bed.layers_.layers_) {
        L += layer.domain_.length;
      }

      // Shell mass: cylinder
      const double mass_shell = std::numbers::pi * D * L * t_w * rho; // [kg]

      // Head mass: 2× standard 2:1 semi-ellipsoidal heads.
      // Equivalent area factor ≈ 1.09 relative to a flat disc (Sinnott & Towler).
      const double mass_heads = 2.0 * 1.09 * (std::numbers::pi / 4.0) * D * D * t_w * rho; // [kg]

      const double vessel_mass = mass_shell + mass_heads;
      const double cost_per_kg = kVesselMaterialCost.at(mat);

      return vessel_mass * cost_per_kg * installation_factor; // [£]
    }

  private:
    friend class cereal::access;
    template<class Archive>
    void serialize(Archive& ar) {
      ar(cereal::make_nvp("installation_factor", installation_factor));
    }
  };

  // ===========================================================================
  // Rotating equipment costs — power-law correlation
  // Peters, Timmerhaus & West "Plant Design and Economics for Chemical Engineers"
  //   C = C_ref × (P / P_ref)^exponent
  //
  // Use CompressorUnit / VacuumPumpUnit for individual machines; collect them
  // in CompressorBank / VacuumPumpBank to model mixed-size configurations.
  // ===========================================================================

  // A rotating machine costed by the power law
  //   cost = reference_cost * (rated_power / reference_power)^cost_exponent * count.
  // CompressorUnit and VacuumPumpUnit differ only in their defaults.
  struct PoweredUnit {
    std::string label;
    double rated_power_kW;        // Design shaft power [kW]
    double reference_cost;        // Reference purchase cost [£]
    double reference_power_kW;
    double cost_exponent = 0.67;
    int    count         = 1;     // Number of identical units

    double CalculateCost() const {
      if (rated_power_kW <= 0.0) return 0.0;
      return reference_cost
           * std::pow(rated_power_kW / reference_power_kW, cost_exponent)
           * count;
    }

    template<class Archive>
    void serialize(Archive& ar) {
      ar(cereal::make_nvp("label", label),
         cereal::make_nvp("rated_power_kW", rated_power_kW),
         cereal::make_nvp("reference_cost", reference_cost),
         cereal::make_nvp("reference_power_kW", reference_power_kW),
         cereal::make_nvp("cost_exponent", cost_exponent),
         cereal::make_nvp("count", count));
    }
  };

  struct CompressorUnit : PoweredUnit {
    CompressorUnit() : PoweredUnit{"Compressor", 100.0, 30000, 100.0} {}
  };

  struct VacuumPumpUnit : PoweredUnit {
    VacuumPumpUnit() : PoweredUnit{"Vacuum pump", 50.0, 20000, 50.0} {}
  };

  // One or more units of the same kind, costed together.
  template <class Unit>
  class EquipmentBank {
  public:
    std::vector<Unit> units = { Unit{} };

    double CalculateCost() const {
      double cost = 0.0;
      for (const auto& u : units) cost += u.CalculateCost();
      return cost;
    }

  private:
    friend class cereal::access;
    template<class Archive>
    void serialize(Archive& ar) {
      ar(cereal::make_nvp("units", units));
    }
  };

  using CompressorBank = EquipmentBank<CompressorUnit>;
  using VacuumPumpBank = EquipmentBank<VacuumPumpUnit>;

  class HeatExchangerCost {
  public:
    bool   include = false;
    double cost    = 0.0; // [£]

    double CalculateCost() const { return include ? cost : 0.0; }

  private:
    friend class cereal::access;
    template<class Archive>
    void serialize(Archive& ar) {
      ar(cereal::make_nvp("include", include),
         cereal::make_nvp("cost", cost));
    }
  };

  // ===========================================================================
  // CostEstimator — Levelized Cost of Capture [£/kg_product]
  //
  // LCOC = (Annualized_CAPEX + Annual_OPEX) / Annual_Production
  //
  // Annualized_CAPEX = Total_CAPEX × CRF
  // CRF = i(1+i)^n / ((1+i)^n - 1)
  // ===========================================================================

  class CostEstimator {
  public:

    // --- Equipment sub-objects ---
    VesselCost        vessel_cost_{};
    CompressorBank    compressors_{};
    VacuumPumpBank    vacuum_pumps_{};
    PipeworkCost      pipework_cost_{};
    ValveCost         valve_cost_{};
    HeatExchangerCost heat_exchanger_cost_{};

    // --- Adsorbent ---
    std::unordered_map<AdsorbentID, double> adsorbent_cost_per_kg{}; // [£/kg]
    double adsorbent_replacement_years = 5.0; // How often adsorbent is fully replaced [years]

    // --- Plant configuration ---
    int    vessel_count            = 1;
    int    train_count             = 1;
    int    operating_days_per_year = 330;

    // --- Financial parameters ---
    double plant_life_years  = 25.0;
    double discount_rate     = 0.08; // [-]

    // --- OPEX rates ---
    double electricity_cost_per_kWh = 0.04;  // [£/kWh]
    double maintenance_fraction     = 0.03;  // fraction of total CAPEX per year
    double overhead_fraction        = 0.02;  // fraction of total CAPEX per year
    double labour_cost_per_year     = 50000; // [£/year]

    // Returns Levelized Cost of Capture [£/kg_product]. Follows the kpis.h
    // undefined-KPI convention: with no net production (or an invalid
    // KPIParams) LCOC is undefined -> returns 0 and sets *defined = false.
    double CalculateCost(const Bed& bed, const AdsorbentDirectory& dir, const KPIParams& params,
                         bool* defined = nullptr) {

      const int n_units = vessel_count * train_count;

      // ---- CAPEX [£] ----
      const double adsorbent_capex   = CalculateAdsorbentFillCost(bed, dir) * n_units;
      const double vessel_capex      = vessel_cost_.CalculateCost(bed) * n_units;
      const double compressor_capex  = compressors_.CalculateCost() * train_count;
      const double vacuum_pump_capex = vacuum_pumps_.CalculateCost() * train_count;
      const double pipework_capex    = pipework_cost_.CalculatePipeworkCost();
      const double valve_capex       = valve_cost_.CalculateCost();
      const double hex_capex         = heat_exchanger_cost_.CalculateCost();

      const double total_capex = adsorbent_capex + vessel_capex + compressor_capex
                               + vacuum_pump_capex + pipework_capex + valve_capex + hex_capex;

      // ---- Annualized CAPEX [£/year] ----
      const double annual_capex = total_capex * CalculateCRF();

      // ---- Annual production [kg/year] ----
      bool ok = true;
      const double productivity_kg_s    = CalculateProductivity(params, &ok); // [kg/s] single vessel
      const double annual_production_kg = productivity_kg_s
                                        * 3600.0 * 24.0 * operating_days_per_year
                                        * n_units;

      if (!ok || !(annual_production_kg > 0.0)) return KPIResult(0.0, false, defined);

      // ---- Annual OPEX [£/year] ----
      const double annual_adsorbent_replacement = CalculateAnnualAdsorbentReplacement(bed, dir) * n_units;

      const double specific_energy_J_per_kg = CalculateSpecificEnergyConsumption(params); // [J/kg]
      const double specific_energy_kWh      = specific_energy_J_per_kg / 3.6e6;           // [kWh/kg]
      const double annual_electricity        = specific_energy_kWh * electricity_cost_per_kWh
                                             * annual_production_kg;                       // [£/year]

      const double annual_maintenance = maintenance_fraction * total_capex;
      const double annual_overhead    = overhead_fraction    * total_capex;

      const double total_annual_opex = annual_adsorbent_replacement
                                     + annual_electricity
                                     + annual_maintenance
                                     + annual_overhead
                                     + labour_cost_per_year;

      // ---- LCOC [£/kg] ----
      const double lcoc = KPIRatio(annual_capex + total_annual_opex, annual_production_kg, ok);
      return KPIResult(lcoc, ok, defined);
    }

    // Cost of a single adsorbent fill across all layers [£]. Prices are keyed
    // by AdsorbentID: each layer's adsorbent mass is split per adsorbent by
    // mass fraction (GetAdsorbentBulkDensities sums to the layer bulk density)
    // and priced at that adsorbent's cost. An adsorbent with no price entry
    // contributes 0.
    double CalculateAdsorbentFillCost(const Bed& bed, const AdsorbentDirectory& dir) const {
      double cost = 0.0;
      const double area = bed.GetArea();
      for (const auto& layer : bed.layers_.layers_) {
        const double volume = layer.domain_.length * area; // [m3 bed]
        for (const auto& [ads_id, bulk_density] : layer.GetAdsorbentBulkDensities(dir)) {
          auto it = adsorbent_cost_per_kg.find(ads_id);
          if (it != adsorbent_cost_per_kg.end()) {
            cost += volume * bulk_density * it->second; // [kg] * [£/kg]
          }
        }
      }
      return cost; // [£]
    }

  private:
    friend class cereal::access;
    template<class Archive>
    void serialize(Archive& ar) {
      ar(cereal::make_nvp("vessel_cost", vessel_cost_),
         cereal::make_nvp("compressors", compressors_),
         cereal::make_nvp("vacuum_pumps", vacuum_pumps_),
         cereal::make_nvp("pipework_cost", pipework_cost_),
         cereal::make_nvp("valve_cost", valve_cost_),
         cereal::make_nvp("heat_exchanger_cost", heat_exchanger_cost_),
         cereal::make_nvp("adsorbent_cost_per_kg", adsorbent_cost_per_kg),
         cereal::make_nvp("adsorbent_replacement_years", adsorbent_replacement_years),
         cereal::make_nvp("vessel_count", vessel_count),
         cereal::make_nvp("train_count", train_count),
         cereal::make_nvp("operating_days_per_year", operating_days_per_year),
         cereal::make_nvp("plant_life_years", plant_life_years),
         cereal::make_nvp("discount_rate", discount_rate),
         cereal::make_nvp("electricity_cost_per_kWh", electricity_cost_per_kWh),
         cereal::make_nvp("maintenance_fraction", maintenance_fraction),
         cereal::make_nvp("overhead_fraction", overhead_fraction),
         cereal::make_nvp("labour_cost_per_year", labour_cost_per_year));
    }

    // Capital Recovery Factor: i(1+i)^n / ((1+i)^n - 1)
    double CalculateCRF() const {
      const double i      = discount_rate;
      const double factor = std::pow(1.0 + i, plant_life_years);
      return (i * factor) / (factor - 1.0);
    }

    // Annual adsorbent replacement cost [£/year]
    double CalculateAnnualAdsorbentReplacement(const Bed& bed, const AdsorbentDirectory& dir) const {
      if (adsorbent_replacement_years <= 0.0) return 0.0;
      return CalculateAdsorbentFillCost(bed, dir) / adsorbent_replacement_years;
    }
  };

} // namespace adsim
