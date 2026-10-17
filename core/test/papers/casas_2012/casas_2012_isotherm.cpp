// ============================================================================
//  Casas / Schell (2012) — ISOTHERM acceptance test.
//
//  Gate before any process (breakthrough) validation: the ArrheniusSipsIsotherm
//  configured with the Schell 2012 constants must reproduce the paper's
//  equilibrium isotherm. Two checks:
//    1. Implementation: the C++ pure-component loading matches an independent
//       (Python/numpy) reference evaluation of the same closed-form equations.
//    2. Fidelity: the continuous model reproduces the paper's per-temperature
//       Sips isotherm (Table 3) across the full measured range 25–140 °C.
//  See casas_2012.h (SetSchellSipsParams) for the fitted constants and README.md
//  for the derivation.
// ============================================================================

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>

#include "adsim/equations/isotherm.h"
#include "adsim/math/constants.h"

#include <gtest/gtest.h>

#include "casas_2012.h"  // SetSchellSipsParams — single source of truth
#include "../../fixtures/validation_metrics.h"

namespace adsim {

  namespace {
    // Pure-component loading via the paper's per-temperature Sips (Table 3):
    //   q = n∞ · (k·p)^c / (1 + (k·p)^c),  k in [1/Pa], p in [Pa].
    double SipsTable3(double ninf, double k_per_MPa, double c, double p_Pa) {
      const double k = k_per_MPa / 1.0e6;
      const double z = std::pow(k * p_Pa, c);
      return ninf * z / (1.0 + z);
    }

    double LoadingCO2(double T, double p_Pa) {
      ArrheniusSipsIsothermParameters co2, h2;
      SetSchellSipsParams(co2, h2);
      double q = 0.0;
      IsothermFunction::ArrheniusSips(q, T, p_Pa, 1.0, co2);
      return q;
    }
    double LoadingH2(double T, double p_Pa) {
      ArrheniusSipsIsothermParameters co2, h2;
      SetSchellSipsParams(co2, h2);
      double q = 0.0;
      IsothermFunction::ArrheniusSips(q, T, p_Pa, 1.0, h2);
      return q;
    }
  } // namespace

  // 1. Implementation correctness: the C++ isotherm must equal its documented
  //    closed form q = n∞(T)·(k(T)·p)^{c(T)} / (1 + (k(T)·p)^{c(T)}), with
  //    n∞, k, c evaluated from the same stored constants. Computed inline here
  //    to catch any coding error in IsothermFunction::ArrheniusSips.
  TEST(Casas2012Isotherm, MatchesClosedForm) {
    const double R = constants::kIdealGasConstant;
    ArrheniusSipsIsothermParameters co2, h2;
    SetSchellSipsParams(co2, h2);

    auto closed_form = [&](const ArrheniusSipsIsothermParameters& p, double T, double pp) {
      const double n_inf = double(p.n_inf_a) * std::exp(double(p.n_inf_b) / (R * T));
      const double k = double(p.k_A) * std::exp(double(p.k_B) / (R * T));
      const double c = double(p.c_ref) + double(p.c_amp) * (2.0 / M_PI) *
                       std::atan(double(p.c_beta) * (T - double(p.c_Tref)));
      const double z = std::pow(k * pp, c);
      return n_inf * z / (1.0 + z);
    };

    for (double Tc : { 25.0, 45.0, 65.0, 100.0, 140.0 }) {
      const double T = Tc + 273.15;
      for (double p_bar : { 0.5, 1.0, 5.0, 15.0, 35.0 }) {
        const double p = p_bar * 1.0e5;
        double q_co2 = 0.0, q_h2 = 0.0;
        IsothermFunction::ArrheniusSips(q_co2, T, p, 1.0, co2);
        IsothermFunction::ArrheniusSips(q_h2, T, p, 1.0, h2);
        EXPECT_NEAR(q_co2, closed_form(co2, T, p), 1e-9) << "CO2 T=" << Tc << " p=" << p_bar;
        EXPECT_NEAR(q_h2, closed_form(h2, T, p), 1e-9) << "H2 T=" << Tc << " p=" << p_bar;
      }
    }
  }

  // 2. Paper fidelity: the continuous model reproduces the paper's per-temperature
  //    Sips isotherm (Table 3) across 25–140 °C. Tolerances reflect the paper's
  //    own Table-3 → Table-4 smoothing (CO2's k is non-monotonic at 45/65 °C).
  TEST(Casas2012Isotherm, ReproducesTable3AcrossTemperatures) {
    // Table 3 (Schell 2012): n∞ [mol/kg], k [1/MPa], c [-].
    struct Row { double Tc, ninf, k, c; };
    const std::vector<Row> co2 = {
      { 25., 13.27, 0.666, 0.74 }, { 45., 12.11, 0.488, 0.77 },
      { 65.,  9.85, 0.498, 0.88 }, { 100., 8.38, 0.340, 0.93 },
      { 140., 7.22, 0.225, 0.94 },
    };
    const std::vector<Row> h2 = {
      { 25., 6.66, 0.037, 0.96 }, { 45., 6.66, 0.029, 0.96 }, { 65., 6.66, 0.023, 0.96 },
    };

    auto sweep = [](const std::vector<Row>& rows, auto model_fn, double tol_pct) {
      double overall = 0.0;
      for (const auto& row : rows) {
        const double T = row.Tc + 273.15;
        double worst = 0.0;
        for (double p_bar = 0.5; p_bar <= 40.0; p_bar += 0.5) {
          const double p = p_bar * 1.0e5;
          const double q_paper = SipsTable3(row.ninf, row.k, row.c, p);
          const double q_model = model_fn(T, p);
          if (q_paper > 1e-6)
            worst = std::max(worst, std::abs(q_model - q_paper) / q_paper * 100.0);
        }
        EXPECT_LT(worst, tol_pct)
            << "T=" << row.Tc << " C: max deviation " << worst << "% exceeds " << tol_pct << "%";
        overall = std::max(overall, worst);
      }
      return overall;
    };

    const double worst_co2 = sweep(co2, [](double T, double p) { return LoadingCO2(T, p); }, 8.0);
    const double worst_h2  = sweep(h2,  [](double T, double p) { return LoadingH2(T, p); }, 2.0);

    // Record in the validation summary (max relative deviation over 0.5-40 bar
    // and all tabulated temperatures; tolerances are the pre-existing ones).
    validation::CaseReport report("casas_2012", "Isotherm_Table3");
    report.CheckError("CO2 loading vs Table 3 Sips (25-140 C, 0.5-40 bar), max deviation", worst_co2 / 100.0,
                      "max_rel_error", 0.08, "-", "Schell et al. 2012, Table 3",
                      "ArrheniusSipsIsotherm with the casas_2012.h constants vs the paper's per-temperature Sips fits");
    report.CheckError("H2 loading vs Table 3 Sips (25-65 C, 0.5-40 bar), max deviation", worst_h2 / 100.0,
                      "max_rel_error", 0.02, "-", "Schell et al. 2012, Table 3",
                      "ArrheniusSipsIsotherm with the casas_2012.h constants vs the paper's per-temperature Sips fits");
  }

} // namespace adsim
