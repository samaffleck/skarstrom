/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/objects/wall_design.h"
#include "adsim/objects/reactor.h"

#include <cmath>
#include <iomanip>
#include <sstream>

namespace adsim {

  namespace {

    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    // B31.3 Table 304.1.1 coefficient Y below 482 °C; above it Y rises, which
    // thins the wall, so 0.4 is conservative there.
    constexpr double kB31_3_Y = 0.4;
    // B31.3 304.1.2(b) and VIII-1 UG-27 thin-shell limit on P / (S E).
    constexpr double kThinShellLimit = 0.385;

    std::string Celsius(double kelvin) {
      std::ostringstream ss;
      ss << std::fixed << std::setprecision(0) << kelvin - 273.15 << " °C";
      return ss.str();
    }

    // Allowable stress S [Pa] at the design temperature, or 0 (Custom).
    double AllowableStress(const WallDesignInput& in) {
      return GetMaxAllowableStress(in.material, in.design_temperature);
    }

    // Nominal internal-pressure thickness [m] (incl. c), or NaN when the
    // code's formula does not apply; `why` says why.
    double InternalThickness(const WallDesignInput& in, double S, bool& thick, std::string& why) {
      const double P = in.design_pressure;
      const double E = in.joint_efficiency;
      const double c = in.corrosion_allowance;
      const double SE = S * E;
      thick = false;
      if (P <= 0.0) return 0.0;
      if (P >= SE) {
        why = "the design pressure is not below S·E: no wall thickness can hold it";
        return kNaN;
      }
      if (in.code == WallDesignCode::ASME_VIII_1) {
        const double R = 0.5 * in.inner_diameter + c;  // corroded inside radius
        if (P <= kThinShellLimit * SE) return P * R / (SE - 0.6 * P) + c;
        thick = true;  // Appendix 1-2(a)(1)
        const double Z = (SE + P) / (SE - P);
        return R * (std::sqrt(Z) - 1.0) + c;
      }
      // ASME B31.3 304.1.2 eq. (3b), W = 1.
      if (P > kThinShellLimit * SE) {
        why = "B31.3 304.1.2(b): P / (S·E) is above 0.385, so the pipe needs special consideration";
        return kNaN;
      }
      const double t = P * (in.inner_diameter + 2.0 * c) / (2.0 * (SE - P * (1.0 - kB31_3_Y)));
      const double D = in.inner_diameter + 2.0 * (t + c);
      if (t >= D / 6.0) {
        why = "B31.3 304.1.2(b): the pressure design thickness is not below D/6, so the pipe needs special consideration";
        return kNaN;
      }
      return t + c;
    }

  } // namespace

  double AllowableExternalPressure(const WallDesignInput& in, double thickness) {
    const double t = thickness - in.corrosion_allowance;   // corroded thickness
    const double Do = in.inner_diameter + 2.0 * thickness;
    const double L = in.length;
    const double Ey = GetElasticModulus(in.material, in.design_temperature);
    const double Sy = YieldStrengthLowerBound(in.material, in.design_temperature);
    if (!(t > 0.0) || !(Do > 0.0) || !(L > 0.0) || !(Ey > 0.0) || !(Sy > 0.0)) return kNaN;

    // Code Case 2286 / VIII-2 4.4.5.1: elastic hoop buckling stress.
    const double Ro = 0.5 * Do;
    const double Mx = L / std::sqrt(Ro * t);
    const double Dt = Do / t;
    double Ch;
    if (Mx >= 2.0 * std::pow(Dt, 0.94)) Ch = 0.55 / Dt;
    else if (Mx > 13.0)                 Ch = 1.12 * std::pow(Mx, -1.058);
    else if (Mx > 1.5)                  Ch = 0.92 / (Mx - 0.579);
    else                                Ch = 1.0;
    const double Fhe = 1.6 * Ch * Ey / Dt;

    // Inelastic buckling stress.
    const double r = Fhe / Sy;
    double Fic;
    if (r >= 2.439)      Fic = Sy;
    else if (r > 0.552)  Fic = 0.7 * Sy * std::pow(r, 0.4);
    else                 Fic = Fhe;

    // Design factor (4.4.2).
    double FS;
    if (Fic <= 0.55 * Sy) FS = 2.0;
    else if (Fic < Sy)    FS = 2.407 - 0.741 * (Fic / Sy);
    else                  FS = 1.667;

    return 2.0 * (Fic / FS) / Dt;
  }

  double MaxAllowableInternalPressure(const WallDesignInput& in, double thickness) {
    const double S = AllowableStress(in);
    const double SE = S * in.joint_efficiency;
    const double c = in.corrosion_allowance;
    const double t = thickness - c;
    if (!(SE > 0.0) || !(t > 0.0) || !(in.inner_diameter > 0.0)) return kNaN;
    if (in.code == WallDesignCode::ASME_VIII_1) {
      const double R = 0.5 * in.inner_diameter + c;
      if (t <= 0.5 * R) return SE * t / (R + 0.6 * t);            // UG-27(c)(1)
      const double Z = std::pow((R + t) / R, 2);                   // Appendix 1-2(a)(1)
      return SE * (Z - 1.0) / (Z + 1.0);
    }
    // B31.3 eq. (3a) solved for P, with the outside diameter D.
    const double D = in.inner_diameter + 2.0 * thickness;
    if (t >= D / 6.0) return kNaN;
    const double P = 2.0 * t * SE / (D - 2.0 * t * kB31_3_Y);
    return P > kThinShellLimit * SE ? kNaN : P;
  }

  WallDesignResult DesignWall(const WallDesignInput& in) {
    WallDesignResult res;
    const double S = AllowableStress(in);
    res.allowable_stress = S;
    if (!(S > 0.0)) {
      res.message = "the wall material is Custom: choose a listed material, whose ASME allowable stress is known";
      return res;
    }
    if (!(in.inner_diameter > 0.0)) {
      res.message = "the column diameter must be positive";
      return res;
    }
    if (!(in.joint_efficiency > 0.0) || in.joint_efficiency > 1.0) {
      res.message = "the joint efficiency must be in (0, 1]";
      return res;
    }
    if (in.corrosion_allowance < 0.0) {
      res.message = "the corrosion allowance cannot be negative";
      return res;
    }
    if (in.design_pressure <= 0.0 && !in.vacuum_service) {
      res.message = "no design pressure: give a positive design pressure or set vacuum service";
      return res;
    }

    std::string why;
    res.internal_thickness = InternalThickness(in, S, res.thick_shell, why);
    if (std::isnan(res.internal_thickness)) {
      res.message = why;
      return res;
    }
    if (res.thick_shell)
      res.notes.push_back("P is above 0.385 S·E: the thick-shell rule of VIII-1 Appendix 1-2 was used");

    if (in.vacuum_service) {
      if (!(in.length > 0.0)) {
        res.message = "vacuum service needs the column length: add a layer";
        return res;
      }
      // The smallest nominal thickness whose allowable external pressure is
      // at least full vacuum (Pa grows with t): bracket, then bisect.
      const double c = in.corrosion_allowance;
      double lo = c, hi = c + 1e-4;
      const double hi_max = c + in.inner_diameter;
      auto holds = [&](double t) {
        const double pa = AllowableExternalPressure(in, t);
        return std::isfinite(pa) && pa >= kFullVacuumExternalPressure;
      };
      while (!holds(hi) && hi < hi_max) { lo = hi; hi *= 2.0; }
      if (!holds(hi)) {
        res.message = "no wall thinner than the column diameter holds full vacuum";
        return res;
      }
      for (int i = 0; i < 200 && (hi - lo) > 1e-9; ++i) {
        const double mid = 0.5 * (lo + hi);
        (holds(mid) ? hi : lo) = mid;
      }
      res.external_thickness = hi;
      res.notes.push_back("full vacuum by the Code Case 2286 formula method over the packed length; "
                          "heads and stiffening rings are not modelled");
    }

    res.required_thickness = std::max(res.internal_thickness, res.external_thickness);
    res.external_governs = res.external_thickness > res.internal_thickness;
    if (in.code == WallDesignCode::ASME_B31_3) {
      res.notes.push_back("B31.3 with W = 1 and Y = 0.4, and ASME II-D (VIII-1) allowable stresses");
      // W < 1 above 427 °C for carbon steel and 510 °C for other materials.
      const bool carbon = in.material == WallMaterial::CarbonSteel_A106
                       || in.material == WallMaterial::CarbonSteel_A53;
      if (in.design_temperature > 273.15 + (carbon ? 427.0 : 510.0))
        res.notes.push_back("above " + Celsius(273.15 + (carbon ? 427.0 : 510.0))
                            + " B31.3 requires a weld strength reduction factor W < 1; W = 1 was used");
    }
    res.valid = true;
    return res;
  }

  std::optional<WallDesignInput> WallDesignInputFor(const Bed& bed) {
    const auto* cyl = std::get_if<CylinderBed>(&bed.shape_);
    if (!cyl) return std::nullopt;
    const auto& w = bed.wall_;
    WallDesignInput in;
    in.code = w.design_code_;
    in.material = w.material;
    in.inner_diameter = cyl->diameter.Get();
    in.length = 0.0;
    for (const auto& layer : bed.layers_.layers_) in.length += layer.domain_.length.Get();
    in.design_pressure = w.design_pressure_.Get();
    in.design_temperature = w.design_temperature_.Get();
    in.joint_efficiency = w.weld_joint_efficiency_.Get();
    in.corrosion_allowance = w.corrosion_allowance_.Get();
    in.vacuum_service = w.vacuum_service_;
    return in;
  }

  WallDesignResult DesignBedWall(const Bed& bed) {
    const auto in = WallDesignInputFor(bed);
    if (!in) {
      WallDesignResult res;
      res.message = "the wall design formulas are for cylindrical columns";
      return res;
    }
    return DesignWall(*in);
  }

} // namespace adsim
