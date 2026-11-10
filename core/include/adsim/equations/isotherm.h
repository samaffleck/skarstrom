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
#include <map>
#include <unordered_map>
#include <memory>
#include <iterator>

// Eigen includes
#include <Eigen/Dense>

// Cereal includes
#include <cereal/types/vector.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/string.hpp>
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>

// adsim includes
#include "adsim/equations/equation.h"
#include "adsim/math/user_variables.h"
#include "adsim/math/expression_param.h"


namespace adsim {

  // Forward declare
  class IsothermModel;
  class PorousMedia;
  class FluidData;
  class Equations;

  enum class IsothermType {
    SINGLE_SITE_LANGMUIR_P,
    DUAL_SITE_LANGMUIR_P,
    INERT, 
    HENRY,
    SINGLE_SITE_LANGMUIR_C,
    DUAL_SITE_LANGMUIR_C,
    EXTENDED_SINGLE_SITE_LANGMUIR_P,
    EXTENDED_DUAL_SITE_LANGMUIR_P,
    EXTENDED_SINGLE_SITE_LANGMUIR_C,
    EXTENDED_DUAL_SITE_LANGMUIR_C,
    ANTI_LANGMUIR,
    FREUNDLICH,
    SIPS,
    SIPS_2,
    EXTENDED_SIPS_2,
    DUAL_SITE_SIPS,
    REDLICH_PETERSON,
    TOTH,
    TOTH_2,
    ARRHENIUS_SIPS,
    EXTENDED_LANGMUIR_FREUNDLICH,
    GAB,
    GAB_YOUNG,
    TOTH_WADST,
    TOTH_HUMID,
  };

  static const char* isothermModelNames[] = {
    "Single Site Langmuir",      // 0
    "Dual Site Langmuir",        // 1
    "Inert",                                        // 2
    "Henry",                     // 3
    "Single Site Langmuir (Concentration)",   // 4
    "Dual Site Langmuir (Concentration)",      // 5
    "Extended Single Site Langmuir",      // 6
    "Extended Dual Site Langmuir",// 7
    "Extended Single Site Langmuir (Concentration)",      // 8
    "Extended Dual Site Langmuir (Concentration)",      // 9
    "Anti-Langmuir",  
    "Freundlich",
    "Sips",
    "Sips 2",
    "Extended Sips 2",
    "Dual Site Sips",
    "Redlich-Peterson",
    "Toth",
    "Toth 2",
    "Arrhenius Sips (Schell 2012)",
    "Extended Langmuir-Freundlich",
    "GAB",
    "GAB (Young 2021 energies)",
    "Toth WADST (CO2 with water, Young 2021)",
    "Toth humid (CO2 with water, Stampi-Bombelli 2020)"
  };

  // The names are indexed by IsothermType, so inserting a value mid-enum
  // without adding the matching name would silently mislabel every isotherm
  // after it in the GUI and the web contract.
  static_assert(std::size(isothermModelNames) ==
                    static_cast<size_t>(IsothermType::TOTH_HUMID) + 1,
                "isothermModelNames must have one entry per IsothermType");

  class InertIsothermParameters {
  public:
    void EvaluateExpressions(const EvalContext&) {}  // nothing to evaluate

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {}
  };

  class HenryIsothermParameters {
  public:
    ExpressionParam k{1.0e-12};
    ExpressionParam Q{4e4};

    void EvaluateExpressions(const EvalContext& ctx) {
      k.Evaluate(ctx); Q.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("IP_k", k), cereal::make_nvp("IP_Q", Q));
    }
  };

  class AntiLangmuirIsothermParameters {
  public:
    ExpressionParam a{1.0e-6};
    ExpressionParam b{9.0e-6};

    void EvaluateExpressions(const EvalContext& ctx) {
      a.Evaluate(ctx); b.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("IP_a", a), cereal::make_nvp("IP_b", b));
    }
  };

  class FreundlichIsothermParameters {
  public:
    ExpressionParam a{1.0};
    ExpressionParam v{2.0};

    void EvaluateExpressions(const EvalContext& ctx) {
      a.Evaluate(ctx); v.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("IP_a", a), cereal::make_nvp("IP_v", v));
    }
  };

  class SipsIsothermParameters {
  public:
    ExpressionParam q_sat{2.0};
    ExpressionParam b{1.0e-5};
    ExpressionParam v{1.0};

    void EvaluateExpressions(const EvalContext& ctx) {
      q_sat.Evaluate(ctx); b.Evaluate(ctx); v.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("IP_q_sat", q_sat), cereal::make_nvp("IP_b", b),
         cereal::make_nvp("IP_v", v));
    }
  };

  class Sips2IsothermParameters {
  public:
    ExpressionParam A1{4.0};        // [mol/kg]
    ExpressionParam A2{-0.6};       // [-]
    ExpressionParam B1{3.6e-4};     // [1/bar]
    ExpressionParam B2{15.0e3};     // [J/mol]
    ExpressionParam C1{1.0};        // [-]
    ExpressionParam C2{0.0};        // [-]
    ExpressionParam T_ref{298};     // [K]

    void EvaluateExpressions(const EvalContext& ctx) {
      A1.Evaluate(ctx); A2.Evaluate(ctx);
      B1.Evaluate(ctx); B2.Evaluate(ctx);
      C1.Evaluate(ctx); C2.Evaluate(ctx);
      T_ref.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("IP_A1", A1), cereal::make_nvp("IP_A2", A2),
         cereal::make_nvp("IP_B1", B1), cereal::make_nvp("IP_B2", B2),
         cereal::make_nvp("IP_C1", C1), cereal::make_nvp("IP_C2", C2),
         cereal::make_nvp("IP_T_ref", T_ref));
    }
  };

  // Temperature-dependent Sips parameters following Schell/Casas et al. (2012),
  // "Pure and binary adsorption of CO2, H2 and N2 on activated carbon", Adsorption.
  //   n_inf(T) = a · exp(b / (R·T))                      [Arrhenius, Eq. 6]
  //   k(T)     = A · exp(B / (R·T))                      [Arrhenius, Eq. 7]
  //   c(T)     = c_ref + c_amp·(2/π)·arctan(c_beta·(T − c_Tref))   [S-shaped exponent]
  // The competitive multicomponent loading (Ruthven 1984 empirical binary Sips) is
  //   q_i = n_inf_i · (k_i·p_i)^{c_i} / (1 + Σ_j (k_j·p_j)^{c_j}),  p_i = P·y_i.
  // For components with a temperature-independent exponent (H2, N2 in the paper)
  // set c_amp = 0, giving c(T) = c_ref. For a temperature-independent saturation
  // capacity (H2) set b = 0.
  class ArrheniusSipsIsothermParameters {
  public:
    ExpressionParam n_inf_a{1.0};    // a   [mol/kg]
    ExpressionParam n_inf_b{0.0};    // b   [J/mol]
    ExpressionParam k_A{1.0e-8};     // A   [1/Pa]
    ExpressionParam k_B{0.0};        // B   [J/mol]
    ExpressionParam c_ref{1.0};      // c_ref   [-]
    ExpressionParam c_amp{0.0};      // c_amp   [-]
    ExpressionParam c_beta{0.0};     // c_beta  [1/K]
    ExpressionParam c_Tref{298.15};  // c_Tref  [K]

    void EvaluateExpressions(const EvalContext& ctx) {
      n_inf_a.Evaluate(ctx); n_inf_b.Evaluate(ctx);
      k_A.Evaluate(ctx); k_B.Evaluate(ctx);
      c_ref.Evaluate(ctx); c_amp.Evaluate(ctx);
      c_beta.Evaluate(ctx); c_Tref.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("IP_n_inf_a", n_inf_a), cereal::make_nvp("IP_n_inf_b", n_inf_b),
         cereal::make_nvp("IP_k_A", k_A), cereal::make_nvp("IP_k_B", k_B),
         cereal::make_nvp("IP_c_ref", c_ref), cereal::make_nvp("IP_c_amp", c_amp),
         cereal::make_nvp("IP_c_beta", c_beta), cereal::make_nvp("IP_c_Tref", c_Tref));
    }
  };

  // Extended Langmuir-Freundlich with temperature-dependent parameters (the
  // form of Yang, "Gas Separation by Adsorption Processes", 1987, used by
  // toPSAil's calcIsothermExtLangFreu.m):
  //   q_i  = q_s,i z_i / (1 + sum_j z_j),   z_i = b_i (P y_i / P_ref)^{n_i}
  //   q_s,i = k1 + k2 T,   b_i = k3 exp(k4 / T),   1 / n_i = k5 + k6 / T
  // The affinity multiplies the power (b P^n), unlike Sips 2's (b P)^n. P_ref
  // carries the pressure unit k3 was fitted in (1e5 Pa: k3 in bar^-n).
  class ExtendedLangmuirFreundlichIsothermParameters {
  public:
    ExpressionParam k1{4.0};        // [mol/kg]
    ExpressionParam k2{0.0};        // [mol/(kg K)]
    ExpressionParam k3{1.0};        // [-] (per P_ref^n)
    ExpressionParam k4{0.0};        // [K]
    ExpressionParam k5{1.0};        // [-]
    ExpressionParam k6{0.0};        // [K]
    ExpressionParam P_ref{1.0e5};   // [Pa]

    void EvaluateExpressions(const EvalContext& ctx) {
      k1.Evaluate(ctx); k2.Evaluate(ctx); k3.Evaluate(ctx); k4.Evaluate(ctx);
      k5.Evaluate(ctx); k6.Evaluate(ctx); P_ref.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("IP_k1", k1), cereal::make_nvp("IP_k2", k2),
         cereal::make_nvp("IP_k3", k3), cereal::make_nvp("IP_k4", k4),
         cereal::make_nvp("IP_k5", k5), cereal::make_nvp("IP_k6", k6),
         cereal::make_nvp("IP_P_ref", P_ref));
    }
  };

  // GAB (Guggenheim-Anderson-de Boer) isotherm for a condensable vapour, in
  // the relative pressure x = p_i / p_sat,i(T), p_i = P y_i:
  //   q = q_m C K x / ((1 - K x) (1 - K x + C K x))
  // with
  //   q_m(T) = q_m0 exp(dH_m / (R T)),  C(T) = C0 exp(dH_C / (R T)),
  //   K(T)   = K0 exp(dH_K / (R T)).
  // BET is the special case K = 1 (K0 = 1, dH_K = 0). p_sat is the pure
  // component's saturation pressure (SaturationPressure). Near
  // the pole K x -> 1, and for x < 0, the loading is continued as described at
  // GabLoading.
  class GabIsothermParameters {
  public:
    ExpressionParam q_m0{3.0};     // [mol/kg] monolayer capacity (at dH_m = 0: q_m itself)
    ExpressionParam dH_m{0.0};     // [J/mol]
    ExpressionParam C0{1.0};       // [-]
    ExpressionParam dH_C{5000.0};  // [J/mol]
    ExpressionParam K0{0.7};       // [-]
    ExpressionParam dH_K{0.0};     // [J/mol]

    void EvaluateExpressions(const EvalContext& ctx) {
      q_m0.Evaluate(ctx); dH_m.Evaluate(ctx); C0.Evaluate(ctx);
      dH_C.Evaluate(ctx); K0.Evaluate(ctx); dH_K.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("IP_q_m0", q_m0), cereal::make_nvp("IP_dH_m", dH_m),
         cereal::make_nvp("IP_C0", C0), cereal::make_nvp("IP_dH_C", dH_C),
         cereal::make_nvp("IP_K0", K0), cereal::make_nvp("IP_dH_K", dH_K));
    }
  };

  // GAB with the layer energies of Young et al. (2021), "The impact of binary
  // water-CO2 isotherm models on the optimal performance of sorbent-based
  // direct air capture processes", Energy Environ. Sci. 14, 5377 (also used by
  // pyDAC):
  //   C = exp((E1 - E10) / (R T)),   K = exp((E2_9 - E10) / (R T)),
  //   E1 = C_E - exp(D T),  E2_9 = F + G T,  E10 = E10_A + E10_B T,
  // with q_m constant; the loading is GabLoading(q_m, C, K, x) as for GAB. E10
  // is the heat of condensation of the bulk liquid (57220 - 44.38 T J/mol for
  // water), E1 the heat of adsorption of the first layer and E2_9 that of
  // layers 2 to 9. The paper calls C_E "C"; it is renamed so the energy is not
  // confused with the GAB constant C. Low et al. (2023), J. Chem. Eng. Data 68,
  // Eqs. 17-22, use the same form. The defaults are Young's Lewatit VP OC 1065
  // values as reprinted by Low et al. (SI Table S5), the energies in J/mol.
  class GabYoungIsothermParameters {
  public:
    ExpressionParam q_m{3.63};       // [mol/kg]
    ExpressionParam C_E{47110.0};    // [J/mol]
    ExpressionParam D{0.023744};     // [1/K]
    ExpressionParam F{57706.0};      // [J/mol]
    ExpressionParam G{-47.814};      // [J/(mol K)]
    ExpressionParam E10_A{57220.0};  // [J/mol]
    ExpressionParam E10_B{-44.38};   // [J/(mol K)]

    void EvaluateExpressions(const EvalContext& ctx) {
      q_m.Evaluate(ctx); C_E.Evaluate(ctx); D.Evaluate(ctx); F.Evaluate(ctx);
      G.Evaluate(ctx); E10_A.Evaluate(ctx); E10_B.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("IP_q_m", q_m), cereal::make_nvp("IP_C_E", C_E),
         cereal::make_nvp("IP_D", D), cereal::make_nvp("IP_F", F),
         cereal::make_nvp("IP_G", G), cereal::make_nvp("IP_E10_A", E10_A),
         cereal::make_nvp("IP_E10_B", E10_B));
    }
  };

  class DualSiteSipsIsothermParameters {
  public:
    ExpressionParam q_sat1{2.0};
    ExpressionParam b1{1.0e-5};
    ExpressionParam v1{1.0};
    ExpressionParam q_sat2{2.0};
    ExpressionParam b2{1.0e-5};
    ExpressionParam v2{1.0};

    void EvaluateExpressions(const EvalContext& ctx) {
      q_sat1.Evaluate(ctx); b1.Evaluate(ctx); v1.Evaluate(ctx);
      q_sat2.Evaluate(ctx); b2.Evaluate(ctx); v2.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("IP_q_sat1", q_sat1), cereal::make_nvp("IP_b1", b1),
         cereal::make_nvp("IP_v1", v1), cereal::make_nvp("IP_q_sat2", q_sat2),
         cereal::make_nvp("IP_b2", b2), cereal::make_nvp("IP_v2", v2));
    }
  };

  class RedlichPetersonIsothermParameters {
  public:
    ExpressionParam a{2.0};
    ExpressionParam b{1.0e-5};
    ExpressionParam v{1.0};

    void EvaluateExpressions(const EvalContext& ctx) {
      a.Evaluate(ctx); b.Evaluate(ctx); v.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("IP_a", a), cereal::make_nvp("IP_b", b),
         cereal::make_nvp("IP_v", v));
    }
  };

  class TothIsothermParameters {
  public:
    ExpressionParam q_sat{2.0};
    ExpressionParam b{1.0e-5};
    ExpressionParam v{1.0};

    void EvaluateExpressions(const EvalContext& ctx) {
      q_sat.Evaluate(ctx); b.Evaluate(ctx); v.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("IP_q_sat", q_sat), cereal::make_nvp("IP_b", b),
         cereal::make_nvp("IP_v", v));
    }
  };

  class Toth2IsothermParameters {
  public:
    ExpressionParam q_sat{3.715}; // [mol/kg]
    ExpressionParam x{0.0};       // [-]
    ExpressionParam T_0{324.0};   // [K]
    ExpressionParam b0{0.027};    // [1/Pa]
    ExpressionParam dH0{100000};  // [J/mol]
    ExpressionParam w0{0.28};     // [-]
    ExpressionParam alpha{0.446}; // [-]

    void EvaluateExpressions(const EvalContext& ctx) {
      q_sat.Evaluate(ctx); x.Evaluate(ctx); T_0.Evaluate(ctx);
      b0.Evaluate(ctx); dH0.Evaluate(ctx); w0.Evaluate(ctx);
      alpha.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("IP_q_sat", q_sat), cereal::make_nvp("IP_x", x),
         cereal::make_nvp("IP_T_0", T_0), cereal::make_nvp("IP_b0", b0),
         cereal::make_nvp("IP_dH0", dH0), cereal::make_nvp("IP_w0", w0),
         cereal::make_nvp("IP_alpha", alpha));
    }
  };

  // ---- CO2 isotherms that depend on the adsorbed water loading ----------------
  //
  // The CO2 loading of these two kinds depends on the loading q_w [mol/kg] of
  // water adsorbed on the same adsorbent in the same cell: the actual
  // (rate-equation) loading of the component named kCoadsorbedWaterComponent,
  // which has its own isotherm (gab, gab_young, ...) and its own mass-transfer
  // model. Both papers use the adsorbed water loading, not its equilibrium
  // value, so with slow water uptake the CO2 capacity follows the water on the
  // sorbent rather than the gas humidity. Without a water component (or with
  // no water adsorbed) both kinds are their dry Toth isotherm. A negative q_w
  // (a Newton iterate) is taken as 0. They are not IAST-mixable.
  inline constexpr char kCoadsorbedWaterComponent[] = "Water";

  // Weighted-average dual-site Toth (WADST) of Young et al. (2021), "The impact
  // of binary water-CO2 isotherm models on the optimal performance of
  // sorbent-based direct air capture processes", Energy Environ. Sci. 14,
  // 5377, Eq. 14:
  //   q = (1 - exp(-A / q_w)) q_dry(p, T) + exp(-A / q_w) q_wet(p, T),
  // each site a temperature-dependent Toth (Eqs. 1-4),
  //   q_site = q_inf(T) b(T) p / (1 + (b(T) p)^tau(T))^(1 / tau(T)),
  //   q_inf = q_inf0 exp(chi (1 - T / T_0)),  b = b0 exp(dH / (R T)),
  //   tau = tau0 + alpha (1 - T_0 / T),
  // with p = P y the CO2 partial pressure [Pa] and dH = -dH_ads > 0 for
  // exothermic adsorption. The affinity has no reference temperature (the
  // paper's form), so b0 is b at infinite temperature; toth_2 writes the same
  // site as b0' exp(dH / (R T_0) (T_0 / T - 1)) with b0' = b0 exp(dH / (R T_0)).
  // At q_w = 0 the loading is exactly the dry site. The defaults are the
  // paper's Lewatit VP OC 1065 values.
  class TothWadstIsothermParameters {
  public:
    ExpressionParam T_0{298.15};         // [K]
    ExpressionParam q_inf_dry{4.86};     // [mol/kg]
    ExpressionParam chi_dry{0.0};        // [-]
    ExpressionParam b0_dry{2.85e-21};    // [1/Pa]
    ExpressionParam dH_dry{117798.0};    // [J/mol]
    ExpressionParam tau0_dry{0.209};     // [-]
    ExpressionParam alpha_dry{0.523};    // [-]
    ExpressionParam q_inf_wet{9.035};    // [mol/kg]
    ExpressionParam chi_wet{0.0};        // [-]
    ExpressionParam b0_wet{1.230e-18};   // [1/Pa]
    ExpressionParam dH_wet{203687.0};    // [J/mol]
    ExpressionParam tau0_wet{0.053};     // [-]
    ExpressionParam alpha_wet{0.053};    // [-]
    ExpressionParam A{1.532};            // [mol/kg]

    void EvaluateExpressions(const EvalContext& ctx) {
      T_0.Evaluate(ctx);
      q_inf_dry.Evaluate(ctx); chi_dry.Evaluate(ctx); b0_dry.Evaluate(ctx);
      dH_dry.Evaluate(ctx); tau0_dry.Evaluate(ctx); alpha_dry.Evaluate(ctx);
      q_inf_wet.Evaluate(ctx); chi_wet.Evaluate(ctx); b0_wet.Evaluate(ctx);
      dH_wet.Evaluate(ctx); tau0_wet.Evaluate(ctx); alpha_wet.Evaluate(ctx);
      A.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("IP_T_0", T_0),
         cereal::make_nvp("IP_q_inf_dry", q_inf_dry), cereal::make_nvp("IP_chi_dry", chi_dry),
         cereal::make_nvp("IP_b0_dry", b0_dry), cereal::make_nvp("IP_dH_dry", dH_dry),
         cereal::make_nvp("IP_tau0_dry", tau0_dry), cereal::make_nvp("IP_alpha_dry", alpha_dry),
         cereal::make_nvp("IP_q_inf_wet", q_inf_wet), cereal::make_nvp("IP_chi_wet", chi_wet),
         cereal::make_nvp("IP_b0_wet", b0_wet), cereal::make_nvp("IP_dH_wet", dH_wet),
         cereal::make_nvp("IP_tau0_wet", tau0_wet), cereal::make_nvp("IP_alpha_wet", alpha_wet),
         cereal::make_nvp("IP_A", A));
    }
  };

  // Humid Toth of Stampi-Bombelli, Storione and Mazzotti (2020), "On
  // comparing packed beds and monoliths for CO2 capture from air through
  // experiments, theory, and modeling", Adsorption 26, 1183-1197 (corrected
  // 2022), Eqs. 11-14: the temperature-dependent Toth of toth_2 (same names)
  //   q = n_s b p / (1 + (b p)^t)^(1 / t),
  //   n_s = q_sat exp(x (1 - T / T_0)) / (1 - gamma q_w),
  //   b   = b0 exp(dH0 / (R T_0) (T_0 / T - 1)) (1 + beta q_w),
  //   t   = w0 + alpha (1 - T_0 / T),
  // with dH0 > 0 for exothermic adsorption (b falls with T). At q_w = 0 this is
  // toth_2 exactly. The defaults are the paper's corrected APDES-NFC values.
  class TothHumidIsothermParameters {
  public:
    ExpressionParam q_sat{2.38};    // [mol/kg]
    ExpressionParam x{0.0};         // [-]
    ExpressionParam T_0{296.0};     // [K]
    ExpressionParam b0{0.0707};     // [1/Pa]
    ExpressionParam dH0{57047.0};   // [J/mol]
    ExpressionParam w0{0.4148};     // [-]
    ExpressionParam alpha{-1.606};  // [-]
    ExpressionParam gamma{0.0016};  // [kg/mol]
    ExpressionParam beta{59.1};     // [kg/mol]

    void EvaluateExpressions(const EvalContext& ctx) {
      q_sat.Evaluate(ctx); x.Evaluate(ctx); T_0.Evaluate(ctx);
      b0.Evaluate(ctx); dH0.Evaluate(ctx); w0.Evaluate(ctx);
      alpha.Evaluate(ctx); gamma.Evaluate(ctx); beta.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("IP_q_sat", q_sat), cereal::make_nvp("IP_x", x),
         cereal::make_nvp("IP_T_0", T_0), cereal::make_nvp("IP_b0", b0),
         cereal::make_nvp("IP_dH0", dH0), cereal::make_nvp("IP_w0", w0),
         cereal::make_nvp("IP_alpha", alpha), cereal::make_nvp("IP_gamma", gamma),
         cereal::make_nvp("IP_beta", beta));
    }
  };

  class LangmuirIsothermParameters {
  public:
    ExpressionParam m1{3.0};
    ExpressionParam b1{1e-11};
    ExpressionParam q1{4e4};

    void EvaluateExpressions(const EvalContext& ctx) {
      m1.Evaluate(ctx); b1.Evaluate(ctx); q1.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("IP_m1", m1), cereal::make_nvp("IP_b1", b1),
         cereal::make_nvp("IP_q1", q1));
    }
  };

  class LangmuirIsothermParametersConcentration {
  public:
    ExpressionParam m1{3.0};
    ExpressionParam b1{1e-7};
    ExpressionParam q1{4e4};

    void EvaluateExpressions(const EvalContext& ctx) {
      m1.Evaluate(ctx); b1.Evaluate(ctx); q1.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("IP_m1", m1), cereal::make_nvp("IP_b1", b1),
         cereal::make_nvp("IP_q1", q1));
    }
  };

  class DualSiteLangmuirIsothermParameters {
  public:
    ExpressionParam m1{3.0};
    ExpressionParam b1{1e-12};
    ExpressionParam q1{4.7e4};
    ExpressionParam m2{3.0};
    ExpressionParam b2{1e-13};
    ExpressionParam q2{4.2e4};

    void EvaluateExpressions(const EvalContext& ctx) {
      m1.Evaluate(ctx); b1.Evaluate(ctx); q1.Evaluate(ctx);
      m2.Evaluate(ctx); b2.Evaluate(ctx); q2.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("IP_m1", m1), cereal::make_nvp("IP_b1", b1),
         cereal::make_nvp("IP_q1", q1), cereal::make_nvp("IP_m2", m2),
         cereal::make_nvp("IP_b2", b2), cereal::make_nvp("IP_q2", q2));
    }
  };

  class DualSiteLangmuirIsothermParametersConcentration {
  public:
    ExpressionParam m1{3.0};
    ExpressionParam b1{1e-8};
    ExpressionParam q1{4e4};
    ExpressionParam m2{3.0};
    ExpressionParam b2{1e-7};
    ExpressionParam q2{4e4};

    void EvaluateExpressions(const EvalContext& ctx) {
      m1.Evaluate(ctx); b1.Evaluate(ctx); q1.Evaluate(ctx);
      m2.Evaluate(ctx); b2.Evaluate(ctx); q2.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("IP_m1", m1), cereal::make_nvp("IP_b1", b1),
         cereal::make_nvp("IP_q1", q1), cereal::make_nvp("IP_m2", m2),
         cereal::make_nvp("IP_b2", b2), cereal::make_nvp("IP_q2", q2));
    }
  };

  // Namespace for isotherm functions
  namespace IsothermFunction {
    void Henry(double& qi, double T, double P, double yi, const HenryIsothermParameters& ips);
    void Langmuir(double& qi, double T, double P, double yi, const LangmuirIsothermParameters& ips);
    void DualSiteLangmuir(double& qi, double T, double P, double yi, const DualSiteLangmuirIsothermParameters& ips);
    void Langmuir_Concentration(double& qi, double T, double P, double yi, const LangmuirIsothermParametersConcentration& ips);
    void DualSiteLangmuir_Concentration(double& qi, double T, double P, double yi, const DualSiteLangmuirIsothermParametersConcentration& ips);
    void AntiLangmuir(double& qi, double T, double P, double yi, const AntiLangmuirIsothermParameters& ips);
    void Freundlich(double& qi, double T, double P, double yi, const FreundlichIsothermParameters& ips);
    void Sips(double& qi, double T, double P, double yi, const SipsIsothermParameters& ips);
    void Sips2(double& qi, double T, double P, double yi, const Sips2IsothermParameters& ips);
    // Single-component loading of the competitive Sips 2 kind (extended_sips_2):
    // Sips2 with the same continuation below y = 1e-5 as ExtendedSips2Isotherm.
    void ExtendedSips2(double& qi, double T, double P, double yi, const Sips2IsothermParameters& ips);
    void DualSiteSips(double& qi, double T, double P, double yi, const DualSiteSipsIsothermParameters& ips);
    void RedlichPeterson(double& qi, double T, double P, double yi, const RedlichPetersonIsothermParameters& ips);
    void Toth(double& qi, double T, double P, double yi, const TothIsothermParameters& ips);
    void Toth2(double& qi, double T, double P, double yi, const Toth2IsothermParameters& ips);
    // Single-component Arrhenius Sips loading (Schell/Casas 2012). The competitive
    // multicomponent form is handled in ArrheniusSipsIsotherm::Update/UpdateAll;
    // this evaluates the pure-component isotherm q = n_inf·(k·p)^c/(1+(k·p)^c).
    void ArrheniusSips(double& qi, double T, double P, double yi, const ArrheniusSipsIsothermParameters& ips);
    // Single-component extended Langmuir-Freundlich loading q = q_s z / (1 + z).
    void ExtendedLangmuirFreundlich(double& qi, double T, double P, double yi, const ExtendedLangmuirFreundlichIsothermParameters& ips);
    // GAB kinds: p_sat [Pa] is the component's saturation pressure at T
    // (SaturationPressure(component, T)); the relative pressure is P yi / p_sat.
    void Gab(double& qi, double T, double P, double yi, double p_sat, const GabIsothermParameters& ips);
    void GabYoung(double& qi, double T, double P, double yi, double p_sat, const GabYoungIsothermParameters& ips);
    // CO2-with-water kinds: q_w [mol/kg] is the adsorbed water loading
    // (negative values are taken as 0; 0 gives the dry isotherm).
    void TothWadst(double& qi, double T, double P, double yi, double q_w, const TothWadstIsothermParameters& ips);
    void TothHumid(double& qi, double T, double P, double yi, double q_w, const TothHumidIsothermParameters& ips);
  }

  // GAB constants at temperature T: q_m [mol/kg], C [-], K [-].
  struct GabConstants { double q_m, C, K; };
  GabConstants GabConstantsAt(double T, const GabIsothermParameters& ips);
  GabConstants GabConstantsAt(double T, const GabYoungIsothermParameters& ips);

  // The relative pressure above which the GAB loading is continued: the pole
  // of q is at x = 1/K, and the loading is exact for x <= x_c = min(1, kGabMaxKx / K).
  inline constexpr double kGabMaxKx = 0.95;

  // GAB loading [mol/kg] at relative pressure x:
  //   q = q_m C K x / ((1 - K x) (1 - K x + C K x))   for 0 <= x <= x_c,
  // and, so that q stays finite, smooth (C1) and increasing everywhere,
  //   q = q_m C K x                                   for x < 0 (the Henry line),
  //   q = q(x_c) + q'(x_c) (x - x_c)                  for x > x_c,
  // with x_c = min(1, kGabMaxKx / K). With K <= kGabMaxKx (the usual case for
  // water) the GAB is therefore exact up to saturation, x = 1; with K near or
  // above 1 (BET, K = 1) it is exact up to K x = kGabMaxKx. Condensation is not
  // modelled.
  double GabLoading(double q_m, double C, double K, double x);

  // Saturation pressure p_sat(T) [Pa] of a pure component (a CoolProp fluid
  // name), from the vapour-pressure ancillary equation of the CoolProp fluid
  // file, held in the generated table component_properties.h: a closed-form
  // fit of the full Helmholtz-energy saturation curve, about 0.1 us per call
  // (water: within 0.02 % of the IAPWS-95 saturation pressure, 273-647 K). T is
  // clamped to [0.75 T_triple, T_c]; above T_c the result is p_c. NaN for a
  // component the table does not have or without an ancillary. Thread-safe
  // (the table is constant).
  double SaturationPressure(const std::string& component, double T);

  struct IsothermUpdateParams {
    IsothermUpdateParams(
      ComponentEquation& yi_,
      ComponentIndependentEquation& qi_sat_,
      vec& T_,
      vec& P_,
      const std::string& component_,
      int start_index_,
      int end_index_,
      const UserVariableRegistry* registry_ = nullptr,
      const ComponentEquation* q_map_ = nullptr
    ) : yi_map(yi_), qi_sat_map(qi_sat_), T(T_), P(P_), component(component_), start_index(start_index_), end_index(end_index_), registry(registry_), q_map(q_map_) {}
    ~IsothermUpdateParams() = default;

    ComponentEquation& yi_map;
    ComponentIndependentEquation& qi_sat_map;
    vec& T;
    vec& P;
    const std::string& component;
    int start_index;
    int end_index;
    const UserVariableRegistry* registry = nullptr;
    // The adsorbent's actual (rate-equation) loadings, component -> q
    // [mol/kg]. Only the kinds that depend on another component's loading
    // read it (toth_wadst, toth_humid: the water loading); null leaves them
    // at their dry limit.
    const ComponentEquation* q_map = nullptr;
  };

  // Per-component isotherm parameters, keyed by component name.
  //
  // An ordered map, not std::unordered_map: competitive isotherms sum over
  // every component (the 1 + sum_j b_j c_j denominators), so iteration order
  // is summation order. An unordered_map's order depends on the standard
  // library's hash and on insertion history (the order of the entries in the
  // deck file, itself the iteration order of whichever library saved it), so a
  // deck and its re-saved copy could run to different numbers. With FMA
  // contraction (Apple clang on arm64 by default) even a two-component sum
  // depends on the order: fma(x2, c2, x1*c1) != fma(x1, c1, x2*c2). std::map
  // iterates in key order whatever the history; cereal writes both map kinds
  // identically ({key, value} arrays), so decks are unaffected.
  template <class P>
  using ParamMap = std::map<std::string, P>;

  // Keys of a per-component parameter map (used by IIsotherm::GetComponents).
  template <class Map>
  std::vector<std::string> ComponentKeys(const Map& params) {
    std::vector<std::string> keys;
    keys.reserve(params.size());
    for (const auto& [component, p] : params) keys.push_back(component);
    return keys;
  }

  class IIsotherm {
  public:
    virtual ~IIsotherm() = default;
    virtual void Update(IsothermUpdateParams& params) = 0;
    // Update several components owned by this isotherm in one call. The
    // default loops over the per-component Update; competitive (extended)
    // isotherms override this to share the per-cell denominator across
    // components instead of recomputing it for each one.
    virtual void UpdateAll(IsothermUpdateParams& params, const std::vector<std::string>& components) {
      for (const auto& comp : components) {
        IsothermUpdateParams p(params.yi_map, params.qi_sat_map, params.T, params.P,
                               comp, params.start_index, params.end_index, params.registry, params.q_map);
        Update(p);
      }
    }
    virtual void RemoveComponent(const std::string& component) = 0;
    // Components this isotherm holds parameters for (the keys of its params_).
    virtual std::vector<std::string> GetComponents() const = 0;
    virtual IsothermType GetIsothermType() = 0;
    virtual void EvaluateExpressions(const EvalContext& ctx) {}

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar) {}
  };

  // Every concrete isotherm is a per-component map of one parameter class plus
  // the same bookkeeping: construction for one component, component removal,
  // key listing, expression evaluation and serialisation. This base supplies
  // all of it; a concrete isotherm only adds Update (and UpdateAll when it is
  // competitive). The archive layout is unchanged: the parameter map is still
  // written under kParamsNvp directly after the IIsotherm base, and each
  // concrete class keeps its own cereal_class_version.
  inline constexpr char kIsothermParamsNvp[] = "Isotherm_Params";
  inline constexpr char kIsoParamsNvp[] = "ISO_Params";  // the Inert and Henry archives

  template <class Params, IsothermType kType, const char* kParamsNvp = kIsothermParamsNvp>
  class ParamIsotherm : public IIsotherm {
  public:
    ParamIsotherm() = default;
    explicit ParamIsotherm(std::string component) { params_[std::move(component)]; }
    void RemoveComponent(const std::string& component) override { params_.erase(component); }
    std::vector<std::string> GetComponents() const override { return ComponentKeys(params_); }
    IsothermType GetIsothermType() override { return kType; }
    void EvaluateExpressions(const EvalContext& ctx) override {
      for (auto& [c, p] : params_) p.EvaluateExpressions(ctx);
    }

    ParamMap<Params> params_{};

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/) {
      ar(
        cereal::virtual_base_class<IIsotherm>(this),
        cereal::make_nvp(kParamsNvp, params_)
      );
    }
  };

  class InertIsotherm : public ParamIsotherm<InertIsothermParameters, IsothermType::INERT, kIsoParamsNvp> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
  };

  class HenryIsotherm : public ParamIsotherm<HenryIsothermParameters, IsothermType::HENRY, kIsoParamsNvp> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
  };

  class LangmuirIsotherm : public ParamIsotherm<LangmuirIsothermParameters, IsothermType::SINGLE_SITE_LANGMUIR_P> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
  };

  class AntiLangmuirIsotherm : public ParamIsotherm<AntiLangmuirIsothermParameters, IsothermType::ANTI_LANGMUIR> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
  };

  class FreundlichIsotherm : public ParamIsotherm<FreundlichIsothermParameters, IsothermType::FREUNDLICH> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
  };
  
  class SipsIsotherm : public ParamIsotherm<SipsIsothermParameters, IsothermType::SIPS> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
  };

  class Sips2Isotherm : public ParamIsotherm<Sips2IsothermParameters, IsothermType::SIPS_2> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
  };

  // Competitive Sips 2: q_i = q_max,i z_i / (1 + sum_j z_j), z_i = (b_i P y_i)^{s_i}.
  // Below y_i = 1e-5 each z_i is continued by a quadratic in y_i (exact above,
  // zero at y_i = 0, finite slope); see docs/isotherms.md, "Competitive Sips
  // kinds near zero mole fraction". ArrheniusSipsIsotherm does the same.
  class ExtendedSips2Isotherm : public ParamIsotherm<Sips2IsothermParameters, IsothermType::EXTENDED_SIPS_2> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
    void UpdateAll(IsothermUpdateParams& state, const std::vector<std::string>& components) override;
  };
  
  // Competitive Sips with Arrhenius n_inf(T), k(T) and an arctan-shaped exponent
  // c(T), matching Schell/Casas et al. (2012). See ArrheniusSipsIsothermParameters
  // for the exact functional forms. Like ExtendedSips2Isotherm this is inherently
  // competitive: the loading of each component shares a common denominator, and
  // each term is continued below y_i = 1e-5 in the same way.
  class ArrheniusSipsIsotherm : public ParamIsotherm<ArrheniusSipsIsothermParameters, IsothermType::ARRHENIUS_SIPS> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
    void UpdateAll(IsothermUpdateParams& state, const std::vector<std::string>& components) override;
  };

  // Competitive (shared denominator) extended Langmuir-Freundlich; see
  // ExtendedLangmuirFreundlichIsothermParameters.
  class ExtendedLangmuirFreundlichIsotherm : public ParamIsotherm<ExtendedLangmuirFreundlichIsothermParameters, IsothermType::EXTENDED_LANGMUIR_FREUNDLICH> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
    void UpdateAll(IsothermUpdateParams& state, const std::vector<std::string>& components) override;
  };

  // GAB kinds: each component's loading depends on its own relative pressure
  // only (not competitive). See GabIsothermParameters / GabYoungIsothermParameters.
  class GabIsotherm : public ParamIsotherm<GabIsothermParameters, IsothermType::GAB> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
  };

  class GabYoungIsotherm : public ParamIsotherm<GabYoungIsothermParameters, IsothermType::GAB_YOUNG> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
  };

  // CO2 kinds that read the adsorbed water loading (IsothermUpdateParams::
  // q_map, component kCoadsorbedWaterComponent); see TothWadstIsothermParameters
  // and TothHumidIsothermParameters.
  class TothWadstIsotherm : public ParamIsotherm<TothWadstIsothermParameters, IsothermType::TOTH_WADST> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
  };

  class TothHumidIsotherm : public ParamIsotherm<TothHumidIsothermParameters, IsothermType::TOTH_HUMID> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
  };

  class DualSiteSipsIsotherm : public ParamIsotherm<DualSiteSipsIsothermParameters, IsothermType::DUAL_SITE_SIPS> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
  };
  
  class RedlichPetersonIsotherm : public ParamIsotherm<RedlichPetersonIsothermParameters, IsothermType::REDLICH_PETERSON> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
  };
  
  class TothIsotherm : public ParamIsotherm<TothIsothermParameters, IsothermType::TOTH> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
  };

  class Toth2Isotherm : public ParamIsotherm<Toth2IsothermParameters, IsothermType::TOTH_2> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
  };

  class ExtendedLangmuirIsotherm : public ParamIsotherm<LangmuirIsothermParameters, IsothermType::EXTENDED_SINGLE_SITE_LANGMUIR_P> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
    void UpdateAll(IsothermUpdateParams& state, const std::vector<std::string>& components) override;
  };

  class ExtendedLangmuirIsotherm_Concentration : public ParamIsotherm<LangmuirIsothermParametersConcentration, IsothermType::EXTENDED_SINGLE_SITE_LANGMUIR_C> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
    void UpdateAll(IsothermUpdateParams& state, const std::vector<std::string>& components) override;
  };

  class DualSiteLangmuirIsotherm : public ParamIsotherm<DualSiteLangmuirIsothermParameters, IsothermType::DUAL_SITE_LANGMUIR_P> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
  };

  class ExtendedDualSiteLangmuirIsotherm : public ParamIsotherm<DualSiteLangmuirIsothermParameters, IsothermType::EXTENDED_DUAL_SITE_LANGMUIR_P> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
    void UpdateAll(IsothermUpdateParams& state, const std::vector<std::string>& components) override;
  };

  class ExtendedDualSiteLangmuirIsotherm_Concentration : public ParamIsotherm<DualSiteLangmuirIsothermParametersConcentration, IsothermType::EXTENDED_DUAL_SITE_LANGMUIR_C> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
    void UpdateAll(IsothermUpdateParams& state, const std::vector<std::string>& components) override;
  };

  class LangmuirIsotherm_Concentration : public ParamIsotherm<LangmuirIsothermParametersConcentration, IsothermType::SINGLE_SITE_LANGMUIR_C> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
  };

  class DualSiteLangmuirIsotherm_Concentration : public ParamIsotherm<DualSiteLangmuirIsothermParametersConcentration, IsothermType::DUAL_SITE_LANGMUIR_C> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
  };

  // Factory Functions
  inline std::unique_ptr<IIsotherm> CreateIsotherm(IsothermType type, std::string component) {
    switch (type) {
    case IsothermType::INERT: return std::make_unique<InertIsotherm>(component);
    case IsothermType::HENRY: return std::make_unique<HenryIsotherm>(component);
    case IsothermType::SINGLE_SITE_LANGMUIR_P: return std::make_unique<LangmuirIsotherm>(component);
    case IsothermType::SINGLE_SITE_LANGMUIR_C: return std::make_unique<LangmuirIsotherm_Concentration>(component);
    case IsothermType::DUAL_SITE_LANGMUIR_P: return std::make_unique<DualSiteLangmuirIsotherm>(component);
    case IsothermType::DUAL_SITE_LANGMUIR_C: return std::make_unique<DualSiteLangmuirIsotherm_Concentration>(component);
    case IsothermType::EXTENDED_SINGLE_SITE_LANGMUIR_P: return std::make_unique<ExtendedLangmuirIsotherm>(component);
    case IsothermType::EXTENDED_SINGLE_SITE_LANGMUIR_C: return std::make_unique<ExtendedLangmuirIsotherm_Concentration>(component);
    case IsothermType::EXTENDED_DUAL_SITE_LANGMUIR_P: return std::make_unique<ExtendedDualSiteLangmuirIsotherm>(component);
    case IsothermType::EXTENDED_DUAL_SITE_LANGMUIR_C: return std::make_unique<ExtendedDualSiteLangmuirIsotherm_Concentration>(component);
    case IsothermType::ANTI_LANGMUIR: return std::make_unique<AntiLangmuirIsotherm>(component);
    case IsothermType::FREUNDLICH: return std::make_unique<FreundlichIsotherm>(component);
    case IsothermType::SIPS: return std::make_unique<SipsIsotherm>(component);
    case IsothermType::SIPS_2: return std::make_unique<Sips2Isotherm>(component);
    case IsothermType::EXTENDED_SIPS_2: return std::make_unique<ExtendedSips2Isotherm>(component);
    case IsothermType::DUAL_SITE_SIPS: return std::make_unique<DualSiteSipsIsotherm>(component);
    case IsothermType::REDLICH_PETERSON: return std::make_unique<RedlichPetersonIsotherm>(component);
    case IsothermType::TOTH: return std::make_unique<TothIsotherm>(component);
    case IsothermType::TOTH_2: return std::make_unique<Toth2Isotherm>(component);
    case IsothermType::ARRHENIUS_SIPS: return std::make_unique<ArrheniusSipsIsotherm>(component);
    case IsothermType::EXTENDED_LANGMUIR_FREUNDLICH: return std::make_unique<ExtendedLangmuirFreundlichIsotherm>(component);
    case IsothermType::GAB: return std::make_unique<GabIsotherm>(component);
    case IsothermType::GAB_YOUNG: return std::make_unique<GabYoungIsotherm>(component);
    case IsothermType::TOTH_WADST: return std::make_unique<TothWadstIsotherm>(component);
    case IsothermType::TOTH_HUMID: return std::make_unique<TothHumidIsotherm>(component);
    default: return nullptr;
    }
  }

} // end adsim namespace

// Class versions live next to the inline serialize templates: every
// translation unit that saves or loads a deck must see the same numbers.
CEREAL_CLASS_VERSION(adsim::InertIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::HenryIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::LangmuirIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::AntiLangmuirIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::FreundlichIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::SipsIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::Sips2Isotherm, 2);
CEREAL_CLASS_VERSION(adsim::ExtendedSips2Isotherm, 2);
CEREAL_CLASS_VERSION(adsim::DualSiteSipsIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::RedlichPetersonIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::TothIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::Toth2Isotherm, 2);
CEREAL_CLASS_VERSION(adsim::ArrheniusSipsIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::ExtendedLangmuirFreundlichIsotherm, 1);
CEREAL_CLASS_VERSION(adsim::GabIsotherm, 1);
CEREAL_CLASS_VERSION(adsim::GabYoungIsotherm, 1);
CEREAL_CLASS_VERSION(adsim::TothWadstIsotherm, 1);
CEREAL_CLASS_VERSION(adsim::TothHumidIsotherm, 1);
CEREAL_CLASS_VERSION(adsim::ExtendedLangmuirIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::DualSiteLangmuirIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::ExtendedDualSiteLangmuirIsotherm, 2);
CEREAL_CLASS_VERSION(adsim::LangmuirIsotherm_Concentration, 2);
CEREAL_CLASS_VERSION(adsim::ExtendedLangmuirIsotherm_Concentration, 2);
CEREAL_CLASS_VERSION(adsim::DualSiteLangmuirIsotherm_Concentration, 2);
CEREAL_CLASS_VERSION(adsim::ExtendedDualSiteLangmuirIsotherm_Concentration, 2);
