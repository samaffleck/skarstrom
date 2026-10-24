/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <optional>
#include <span>

// Cereal includes
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>

#include "adsim/math/expression_param.h"


namespace adsim {

  enum class WallMaterial {
    Custom,
    SS316L,
    SS304,
    SS304L,
    CarbonSteel_A106,
    CarbonSteel_A53,
    Hastelloy_C276,
    Inconel_625,
    Aluminium_6061,
    Titanium_Gr2
  };

  struct WallMaterialProperties {
    double density;               // [kg/m3]
    double heat_capacity;         // [J/(kg-K)]
    double thermal_conductivity;  // [W/(m-K)]
  };

  // Reference values at ~20°C
  // Sources: ASM Handbook, ASME B31.3, MatWeb
  inline constexpr const char* kWallMaterialList[] = {
    "Custom",
    "SS316L",
    "SS304",
    "SS304L",
    "Carbon Steel (A106)",
    "Carbon Steel (A53)",
    "Hastelloy C-276",
    "Inconel 625",
    "Aluminium 6061",
    "Titanium Gr.2"
  };

  inline WallMaterialProperties GetWallMaterialProperties(WallMaterial m) {
    switch (m) {
      case WallMaterial::SS316L:
        return {7990, 500, 16.3};
      case WallMaterial::SS304:
        return {7900, 500, 16.2};
      case WallMaterial::SS304L:
        return {7900, 500, 16.2};
      case WallMaterial::CarbonSteel_A106:
        return {7850, 490, 51.9};
      case WallMaterial::CarbonSteel_A53:
        return {7850, 490, 51.9};
      case WallMaterial::Hastelloy_C276:
        return {8890, 427, 10.2};
      case WallMaterial::Inconel_625:
        return {8440, 410, 9.8};
      case WallMaterial::Aluminium_6061:
        return {2700, 896, 167.0};
      case WallMaterial::Titanium_Gr2:
        return {4510, 520, 16.4};
      default:
        return {7990, 500, 16.3}; // fallback: SS316L
    }
  }

  // ---------------------------------------------------------------------------
  // Maximum Allowable Stress — ASME Section II Part D
  // ---------------------------------------------------------------------------
  // Temperature-stress data points for piecewise linear interpolation.
  // Source: ASME BPVC Section II Part D, Table 1A (ferrous) and Table 1B
  //         (non-ferrous), ~2023 edition. Values at ~20°C (68°F) through the
  //         maximum listed temperature for each alloy.
  // ---------------------------------------------------------------------------

  struct StressDataPoint {
    double temperature_C;  // [°C]
    double max_stress;     // [Pa]
  };

  // Returns the maximum allowable stress [Pa] for a material at temperature_K
  // using piecewise linear interpolation between ASME tabulated values.
  // Clamps to the lowest tabulated stress if temperature_K exceeds the table range.
  // Returns 0.0 for WallMaterial::Custom (caller must supply the value).
  inline double GetMaxAllowableStress(WallMaterial m, double temperature_K) {
    const double T_C = temperature_K - 273.15;

    auto interpolate = [&](std::span<const StressDataPoint> pts) -> double {
      assert(!pts.empty());
      if (T_C <= pts.front().temperature_C) return pts.front().max_stress;
      if (T_C >= pts.back().temperature_C)  return pts.back().max_stress;
      for (std::size_t i = 1; i < pts.size(); ++i) {
        if (T_C <= pts[i].temperature_C) {
          const double frac = (T_C - pts[i-1].temperature_C)
                            / (pts[i].temperature_C - pts[i-1].temperature_C);
          return pts[i-1].max_stress + frac * (pts[i].max_stress - pts[i-1].max_stress);
        }
      }
      return pts.back().max_stress;
    };

    switch (m) {
      // SS316L — ASME IID Table 1A, UNS S31603
      case WallMaterial::SS316L: {
        static constexpr StressDataPoint pts[] = {
          { -29.0, 115.1e6 }, {  38.0, 115.1e6 }, {  93.0, 109.6e6 },
          { 149.0, 104.8e6 }, { 204.0, 101.4e6 }, { 260.0,  97.2e6 },
          { 316.0,  93.1e6 }, { 371.0,  86.2e6 }, { 427.0,  77.9e6 },
          { 482.0,  62.1e6 }, { 538.0,  48.3e6 }, { 566.0,  41.4e6 }
        };
        return interpolate(pts);
      }
      // SS304 — ASME IID Table 1A, UNS S30400
      case WallMaterial::SS304: {
        static constexpr StressDataPoint pts[] = {
          { -29.0, 129.6e6 }, {  38.0, 129.6e6 }, {  93.0, 119.3e6 },
          { 149.0, 110.3e6 }, { 204.0, 104.1e6 }, { 260.0,  99.3e6 },
          { 316.0,  96.5e6 }, { 371.0,  91.7e6 }, { 427.0,  82.7e6 },
          { 482.0,  71.0e6 }, { 538.0,  57.2e6 }, { 566.0,  48.3e6 }
        };
        return interpolate(pts);
      }
      // SS304L — ASME IID Table 1A, UNS S30403
      case WallMaterial::SS304L: {
        static constexpr StressDataPoint pts[] = {
          { -29.0, 115.1e6 }, {  38.0, 115.1e6 }, {  93.0,  99.3e6 },
          { 149.0,  91.0e6 }, { 204.0,  85.5e6 }, { 260.0,  82.1e6 },
          { 316.0,  78.6e6 }, { 371.0,  75.8e6 }, { 427.0,  71.0e6 },
          { 482.0,  60.0e6 }, { 538.0,  46.9e6 }
        };
        return interpolate(pts);
      }
      // Carbon Steel A106 Grade B — ASME IID Table 1A, SA-106 Gr.B
      case WallMaterial::CarbonSteel_A106: {
        static constexpr StressDataPoint pts[] = {
          { -29.0, 120.7e6 }, { 316.0, 120.7e6 },
          { 371.0, 111.0e6 }, { 427.0,  87.6e6 },
          { 482.0,  54.5e6 }, { 538.0,  27.6e6 }
        };
        return interpolate(pts);
      }
      // Carbon Steel A53 Grade B — same allowable stress as SA-106 Gr.B
      case WallMaterial::CarbonSteel_A53: {
        static constexpr StressDataPoint pts[] = {
          { -29.0, 120.7e6 }, { 316.0, 120.7e6 },
          { 371.0, 111.0e6 }, { 427.0,  87.6e6 },
          { 482.0,  54.5e6 }, { 538.0,  27.6e6 }
        };
        return interpolate(pts);
      }
      // Hastelloy C-276 — ASME IID Table 1B, UNS N10276
      case WallMaterial::Hastelloy_C276: {
        static constexpr StressDataPoint pts[] = {
          { -29.0, 194.4e6 }, {  38.0, 194.4e6 }, {  93.0, 176.5e6 },
          { 149.0, 164.1e6 }, { 204.0, 154.4e6 }, { 260.0, 146.2e6 },
          { 316.0, 139.3e6 }, { 371.0, 133.1e6 }, { 427.0, 128.2e6 },
          { 482.0, 124.1e6 }, { 538.0, 120.7e6 }
        };
        return interpolate(pts);
      }
      // Inconel 625 — ASME IID Table 1B, UNS N06625 (annealed)
      case WallMaterial::Inconel_625: {
        static constexpr StressDataPoint pts[] = {
          { -29.0, 198.6e6 }, {  38.0, 198.6e6 }, {  93.0, 179.3e6 },
          { 149.0, 160.7e6 }, { 204.0, 148.3e6 }, { 260.0, 142.0e6 },
          { 316.0, 137.9e6 }, { 371.0, 133.8e6 }, { 427.0, 129.6e6 },
          { 482.0, 125.5e6 }, { 538.0, 122.7e6 }
        };
        return interpolate(pts);
      }
      // Aluminium 6061-T6 — ASME IID Table 1B, UNS A96061
      // Note: strength drops rapidly above ~150°C
      case WallMaterial::Aluminium_6061: {
        static constexpr StressDataPoint pts[] = {
          { -29.0, 103.4e6 }, {  66.0, 103.4e6 }, {  93.0,  96.5e6 },
          { 121.0,  82.7e6 }, { 149.0,  48.3e6 }, { 177.0,  27.6e6 },
          { 204.0,  17.2e6 }
        };
        return interpolate(pts);
      }
      // Titanium Grade 2 — ASME IID Table 1B, UNS R50400
      case WallMaterial::Titanium_Gr2: {
        static constexpr StressDataPoint pts[] = {
          { -29.0,  98.6e6 }, {  38.0,  98.6e6 }, {  93.0,  83.4e6 },
          { 149.0,  76.5e6 }, { 204.0,  70.3e6 }, { 260.0,  64.8e6 },
          { 316.0,  49.6e6 }, { 371.0,  31.0e6 }
        };
        return interpolate(pts);
      }
      default:
        return 0.0; // Custom: caller must supply
    }
  }

  // ---------------------------------------------------------------------------
  // Elastic modulus — ASME Section II Part D, Tables TM
  // ---------------------------------------------------------------------------
  // Approximate values (converted from the psi tables at 70, 200, 300, ... °F)
  // after ASME BPVC II-D Tables TM-1 (carbon and austenitic steels), TM-2
  // (aluminium alloys), TM-4 (nickel alloys) and TM-5 (titanium). Used only by
  // the external-pressure (buckling) check of the wall design; verify them
  // against the edition in use for a final design.
  // Returns 0.0 for WallMaterial::Custom.
  inline double GetElasticModulus(WallMaterial m, double temperature_K) {
    const double T_C = temperature_K - 273.15;
    auto interpolate = [&](std::span<const StressDataPoint> pts) -> double {
      if (T_C <= pts.front().temperature_C) return pts.front().max_stress;
      if (T_C >= pts.back().temperature_C)  return pts.back().max_stress;
      for (std::size_t i = 1; i < pts.size(); ++i) {
        if (T_C <= pts[i].temperature_C) {
          const double frac = (T_C - pts[i-1].temperature_C)
                            / (pts[i].temperature_C - pts[i-1].temperature_C);
          return pts[i-1].max_stress + frac * (pts[i].max_stress - pts[i-1].max_stress);
        }
      }
      return pts.back().max_stress;
    };
    switch (m) {
      // Austenitic stainless steels (TM-1, 304/316 group)
      case WallMaterial::SS316L:
      case WallMaterial::SS304:
      case WallMaterial::SS304L: {
        static constexpr StressDataPoint pts[] = {
          {  21.0, 195.1e9 }, {  93.0, 190.3e9 }, { 149.0, 186.2e9 }, { 204.0, 182.7e9 },
          { 260.0, 177.9e9 }, { 316.0, 174.4e9 }, { 371.0, 171.0e9 }, { 427.0, 166.2e9 },
          { 482.0, 162.0e9 }, { 538.0, 157.2e9 }
        };
        return interpolate(pts);
      }
      // Carbon steels, C <= 0.30 % (TM-1)
      case WallMaterial::CarbonSteel_A106:
      case WallMaterial::CarbonSteel_A53: {
        static constexpr StressDataPoint pts[] = {
          {  21.0, 202.7e9 }, {  93.0, 198.6e9 }, { 149.0, 195.1e9 }, { 204.0, 192.4e9 },
          { 260.0, 188.2e9 }, { 316.0, 182.7e9 }, { 371.0, 175.8e9 }, { 427.0, 166.9e9 },
          { 482.0, 154.4e9 }, { 538.0, 140.7e9 }
        };
        return interpolate(pts);
      }
      // Hastelloy C-276, UNS N10276 (TM-4)
      case WallMaterial::Hastelloy_C276: {
        static constexpr StressDataPoint pts[] = {
          {  21.0, 205.5e9 }, {  93.0, 201.3e9 }, { 149.0, 198.6e9 }, { 204.0, 195.8e9 },
          { 260.0, 193.1e9 }, { 316.0, 190.3e9 }, { 371.0, 187.5e9 }, { 427.0, 184.8e9 },
          { 482.0, 182.0e9 }, { 538.0, 179.3e9 }
        };
        return interpolate(pts);
      }
      // Inconel 625, UNS N06625 (TM-4)
      case WallMaterial::Inconel_625: {
        static constexpr StressDataPoint pts[] = {
          {  21.0, 207.5e9 }, {  93.0, 203.4e9 }, { 149.0, 200.6e9 }, { 204.0, 197.9e9 },
          { 260.0, 195.1e9 }, { 316.0, 192.4e9 }, { 371.0, 189.6e9 }, { 427.0, 186.8e9 },
          { 482.0, 183.4e9 }, { 538.0, 180.0e9 }
        };
        return interpolate(pts);
      }
      // Aluminium 6061 (TM-2)
      case WallMaterial::Aluminium_6061: {
        static constexpr StressDataPoint pts[] = {
          {  21.0, 68.9e9 }, {  93.0, 66.2e9 }, { 149.0, 63.4e9 }, { 204.0, 60.0e9 }
        };
        return interpolate(pts);
      }
      // Titanium Grade 2, UNS R50400 (TM-5)
      case WallMaterial::Titanium_Gr2: {
        static constexpr StressDataPoint pts[] = {
          {  21.0, 106.9e9 }, {  93.0, 103.4e9 }, { 149.0, 100.7e9 }, { 204.0, 96.5e9 },
          { 260.0,  91.7e9 }, { 316.0,  86.9e9 }, { 371.0,  82.0e9 }
        };
        return interpolate(pts);
      }
      default:
        return 0.0;
    }
  }

  // Austenitic stainless steels and nickel alloys may be given an allowable
  // stress up to 0.9 of the yield strength at temperature (ASME II-D Mandatory
  // Appendix 1); every other material is limited to 2/3 of it.
  inline bool IsAusteniticOrNickelAlloy(WallMaterial m) {
    return m == WallMaterial::SS316L || m == WallMaterial::SS304 || m == WallMaterial::SS304L
        || m == WallMaterial::Hastelloy_C276 || m == WallMaterial::Inconel_625;
  }

  // A lower bound on the yield strength at temperature [Pa], from the
  // allowable stress: II-D sets S <= 2/3 Sy(T) (0.9 Sy(T) for austenitic and
  // nickel alloys), so Sy(T) >= 1.5 S (S / 0.9). The external-pressure check
  // uses it in place of a yield-strength table; a lower yield strength gives a
  // lower allowable external pressure, so the bound is conservative.
  inline double YieldStrengthLowerBound(WallMaterial m, double temperature_K) {
    const double S = GetMaxAllowableStress(m, temperature_K);
    return IsAusteniticOrNickelAlloy(m) ? S / 0.9 : 1.5 * S;
  }

  // The code a wall's minimum thickness is calculated to (see wall_design.h).
  enum class WallDesignCode {
    ASME_VIII_1,  // ASME BPVC Section VIII Division 1, UG-27 (pressure vessel shell)
    ASME_B31_3    // ASME B31.3 304.1.2 (process piping: tube and pipe columns)
  };

  inline constexpr const char* kWallDesignCodeList[] = {
    "ASME BPVC VIII-1 (UG-27)",
    "ASME B31.3 (304.1.2)"
  };


	class BedWall {
	public:
    WallMaterial material = WallMaterial::SS316L;
    ExpressionParam thickness_{0.02}; // [m]

    // Read-only material properties — set via ApplyMaterial()
    ExpressionParam density_{7990};              // [kg/m3]
    ExpressionParam heat_capacity_{500};         // [J/(kg-K)]
    ExpressionParam thermal_conductivity_{16.3}; // [W/(m-K)]

    // Pressure design (see wall_design.h). The design is in use when the
    // design pressure is positive or the wall is in vacuum service; then the
    // solver warns when a run leaves the design envelope.
    WallDesignCode  design_code_ = WallDesignCode::ASME_VIII_1;
    ExpressionParam design_pressure_{0.0};       // [Pa gauge] — maximum internal (operating) pressure
    ExpressionParam design_temperature_{293.0};  // [K]        — for allowable stress lookup
    ExpressionParam weld_joint_efficiency_{1.0}; // E: 1.0 seamless/full RT, 0.85 spot, 0.70 none
    ExpressionParam corrosion_allowance_{0.0};   // [m] — added to the required thickness
    bool            vacuum_service_ = false;     // if true, also design for full vacuum (external pressure)

    void ApplyMaterial(WallMaterial m) {
      material = m;
      if (m != WallMaterial::Custom) {
        const auto props = GetWallMaterialProperties(m);
        density_             = props.density;
        heat_capacity_       = props.heat_capacity;
        thermal_conductivity_ = props.thermal_conductivity;
      }
    }

    bool DesignInUse() const {
      return design_pressure_.Get() > 0.0 || vacuum_service_;
    }

    void EvaluateExpressions(const EvalContext& ctx) {
      thickness_.Evaluate(ctx);
      density_.Evaluate(ctx);
      heat_capacity_.Evaluate(ctx);
      thermal_conductivity_.Evaluate(ctx);
      design_pressure_.Evaluate(ctx);
      design_temperature_.Evaluate(ctx);
      weld_joint_efficiency_.Evaluate(ctx);
      corrosion_allowance_.Evaluate(ctx);
    }

	private:
		friend class cereal::access;
		template<class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
        cereal::make_nvp("Material", material),
				cereal::make_nvp("Density", density_),
				cereal::make_nvp("Specific_Heat_Capacity", heat_capacity_),
				cereal::make_nvp("Thermal_Conductivity", thermal_conductivity_),
				cereal::make_nvp("Wall_Thickness", thickness_),
        cereal::make_nvp("Design_Pressure", design_pressure_),
        cereal::make_nvp("Design_Temperature", design_temperature_),
        cereal::make_nvp("Weld_Joint_Efficiency", weld_joint_efficiency_),
        cereal::make_nvp("Vacuum_Service", vacuum_service_)
			);
      if (version >= 4) {
        ar(
          cereal::make_nvp("Design_Code", design_code_),
          cereal::make_nvp("Corrosion_Allowance", corrosion_allowance_)
        );
      }
		}
	};

} // adsim namespace

// Version 4 adds Design_Code and Corrosion_Allowance. Decks saved before it
// carry version 0 (the version was declared in wall.cpp, where no archive saw
// it) and have every field up to Vacuum_Service.
CEREAL_CLASS_VERSION(adsim::BedWall, 4);
