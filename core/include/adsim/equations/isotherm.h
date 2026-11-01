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
    "Arrhenius Sips (Schell 2012)"
  };

  // The names are indexed by IsothermType, so inserting a value mid-enum
  // without adding the matching name would silently mislabel every isotherm
  // after it in the GUI and the web contract.
  static_assert(std::size(isothermModelNames) ==
                    static_cast<size_t>(IsothermType::ARRHENIUS_SIPS) + 1,
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
    void DualSiteSips(double& qi, double T, double P, double yi, const DualSiteSipsIsothermParameters& ips);
    void RedlichPeterson(double& qi, double T, double P, double yi, const RedlichPetersonIsothermParameters& ips);
    void Toth(double& qi, double T, double P, double yi, const TothIsothermParameters& ips);
    void Toth2(double& qi, double T, double P, double yi, const Toth2IsothermParameters& ips);
    // Single-component Arrhenius Sips loading (Schell/Casas 2012). The competitive
    // multicomponent form is handled in ArrheniusSipsIsotherm::Update/UpdateAll;
    // this evaluates the pure-component isotherm q = n_inf·(k·p)^c/(1+(k·p)^c).
    void ArrheniusSips(double& qi, double T, double P, double yi, const ArrheniusSipsIsothermParameters& ips);
  }

  struct IsothermUpdateParams {
    IsothermUpdateParams(
      ComponentEquation& yi_,
      ComponentIndependentEquation& qi_sat_,
      vec& T_,
      vec& P_,
      const std::string& component_,
      int start_index_,
      int end_index_,
      const UserVariableRegistry* registry_ = nullptr
    ) : yi_map(yi_), qi_sat_map(qi_sat_), T(T_), P(P_), component(component_), start_index(start_index_), end_index(end_index_), registry(registry_) {}
    ~IsothermUpdateParams() = default;

    ComponentEquation& yi_map;
    ComponentIndependentEquation& qi_sat_map;
    vec& T;
    vec& P;
    const std::string& component;
    int start_index;
    int end_index;
    const UserVariableRegistry* registry = nullptr;
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
                               comp, params.start_index, params.end_index, params.registry);
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

  class ExtendedSips2Isotherm : public ParamIsotherm<Sips2IsothermParameters, IsothermType::EXTENDED_SIPS_2> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
    void UpdateAll(IsothermUpdateParams& state, const std::vector<std::string>& components) override;
  };
  
  // Competitive Sips with Arrhenius n_inf(T), k(T) and an arctan-shaped exponent
  // c(T), matching Schell/Casas et al. (2012). See ArrheniusSipsIsothermParameters
  // for the exact functional forms. Like ExtendedSips2Isotherm this is inherently
  // competitive: the loading of each component shares a common denominator.
  class ArrheniusSipsIsotherm : public ParamIsotherm<ArrheniusSipsIsothermParameters, IsothermType::ARRHENIUS_SIPS> {
  public:
    using ParamIsotherm::ParamIsotherm;
    void Update(IsothermUpdateParams& state) override;
    void UpdateAll(IsothermUpdateParams& state, const std::vector<std::string>& components) override;
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
    default: return nullptr;
    }
  }

} // end adsim namespace
