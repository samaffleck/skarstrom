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
  // ASME BPV VIII-1 Minimum Wall Thickness
  // ---------------------------------------------------------------------------

  struct AsmeWallThicknessResult {
    double min_thickness;    // [m]  — minimum required wall thickness per ASME
    double allowable_stress; // [Pa] — S at design temperature
    bool   is_valid;         // false when P >= S*E (UG-27 formula breaks down)
  };

  // Calculates the minimum required wall thickness for a cylindrical pressure
  // vessel per ASME BPV Section VIII Division 1.
  //
  // Internal pressure (UG-27, circumferential stress):
  //   t = P * R / (S * E - 0.6 * P)
  //
  // Vacuum / external pressure:
  //   A conservative simplification is used: the vessel is designed for full
  //   atmospheric vacuum differential (101 325 Pa) applied as an equivalent
  //   internal pressure using the UG-27 formula. This is conservative relative
  //   to the iterative UG-28 chart procedure and is appropriate for preliminary
  //   sizing. A formal UG-28 analysis is required for final design.
  //
  // Parameters:
  //   inner_diameter        [m]   — vessel inside diameter
  //   design_pressure       [Pa gauge] — maximum allowable working pressure
  //   design_temperature    [K]   — design temperature for S lookup
  //   material              —     wall material (Custom returns is_valid=false)
  //   weld_joint_efficiency —     E: 1.0 seamless/full RT, 0.85 spot, 0.70 none
  //   vacuum_service        —     if true, design for full external vacuum
  inline AsmeWallThicknessResult CalculateAsmeMinWallThickness(
      double inner_diameter,
      double design_pressure,
      double design_temperature,
      WallMaterial material,
      double weld_joint_efficiency,
      bool   vacuum_service)
  {
    const double S = GetMaxAllowableStress(material, design_temperature);
    if (S <= 0.0) {
      return {0.0, S, false};
    }

    const double R = inner_diameter / 2.0;
    const double E = weld_joint_efficiency;

    // For vacuum service use full atmospheric differential as the effective
    // design pressure (conservative simplification — see function comment above).
    const double P = vacuum_service ? 101325.0 : design_pressure;

    const double denominator = S * E - 0.6 * P;
    if (denominator <= 0.0) {
      return {0.0, S, false};
    }

    const double t = P * R / denominator;
    return {t, S, true};
  }


	class BedWall {
	public:
    WallMaterial material = WallMaterial::SS316L;
    ExpressionParam thickness_{0.02}; // [m]

    // Read-only material properties — set via ApplyMaterial()
    ExpressionParam density_{7990};              // [kg/m3]
    ExpressionParam heat_capacity_{500};         // [J/(kg-K)]
    ExpressionParam thermal_conductivity_{16.3}; // [W/(m-K)]

    // Design parameters for ASME minimum wall thickness calculation
    ExpressionParam design_pressure_{0.0};       // [Pa gauge] — user-specified MAWP
    ExpressionParam design_temperature_{293.0};  // [K]        — for allowable stress lookup
    ExpressionParam weld_joint_efficiency_{1.0}; // E: 1.0 seamless/full RT, 0.85 spot, 0.70 none
    bool            vacuum_service_ = false;     // if true, design for external pressure (vacuum)

    void ApplyMaterial(WallMaterial m) {
      material = m;
      if (m != WallMaterial::Custom) {
        const auto props = GetWallMaterialProperties(m);
        density_             = props.density;
        heat_capacity_       = props.heat_capacity;
        thermal_conductivity_ = props.thermal_conductivity;
      }
    }

    // Returns the ASME minimum wall thickness result using the stored design
    // parameters. Caller must ensure the column shape is CylinderBed — the
    // ASME BPV VIII-1 formulae only apply to cylindrical pressure vessels.
    AsmeWallThicknessResult CalculateMinThickness(double inner_diameter) const {
      return CalculateAsmeMinWallThickness(
          inner_diameter,
          design_pressure_.Get(),
          design_temperature_.Get(),
          material,
          weld_joint_efficiency_.Get(),
          vacuum_service_);
    }

    void EvaluateExpressions(const EvalContext& ctx) {
      thickness_.Evaluate(ctx);
      density_.Evaluate(ctx);
      heat_capacity_.Evaluate(ctx);
      thermal_conductivity_.Evaluate(ctx);
      design_pressure_.Evaluate(ctx);
      design_temperature_.Evaluate(ctx);
      weld_joint_efficiency_.Evaluate(ctx);
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
		}
	};

} // adsim namespace
