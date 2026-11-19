/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <cereal/types/base_class.hpp>
#include <cereal/types/string.hpp>
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>

#include "adsim/equations/transport_properties.h"
#include "adsim/math/expression_param.h"


namespace adsim {

  class Reactor;
  class Wall;
  class BedWall;
  class PorousMedia;
  class Adsorbent;
  class MixedAdsorbentLayer;
  class Equations;

  // Stored values; never renumber.
  enum class InternalHeatTransferType {
    CONSTANT = 0,
    DIXON_1996 = 1,
    LEVA_1947 = 2,
    LI_FINLAYSON_1977 = 3,
    MARTIN_NILLES_1993 = 4
  };

  static const char* InternalHeatTransferModeList[] = {
    "Constant", // 0
    "Dixon (1996)", // 1
    "Leva (1947)", // 2
    "Li and Finlayson (1977)", // 3
    "Martin and Nilles (1993)" // 4
  };

  // ---- Bed-to-wall heat transfer in a 1-D bed model --------------------------
  //
  // A packed-bed wall correlation gives the wall coefficient h_w of a 2-D
  // (radial) pseudo-homogeneous model. A 1-D model has no radial profile, so
  // h_w is lumped in series with the radial conduction resistance of the bed
  // (effective radial conductivity k_e) into one apparent coefficient h_i
  // (Dixon 1996):
  //     1/h_i = 1/h_w + D_c/(6 k_e) (Bi + 3)/(Bi + 4),   Bi = h_w D_c / (2 k_e).
  // The kinds that use it differ in h_w (and k_e). h_w <= 0 gives h_i = 0.
  double LumpedWallHeatTransferCoefficient(double h_w, double k_e, double D_c);

  // Dixon's wall coefficient h_w [W/(m2 K)]: a stagnant term plus the film
  // term (k_g/d_p) 0.0835 Re^0.91. re is the particle Reynolds number (>= 0).
  double DixonWallCoefficient(double k_g, double k_s, double eb, double dp, double D_c, double re);

  // Dixon's effective radial bed conductivity k_e [W/(m K)]: the stagnant
  // term k_g [eb + (1 - eb)/((2/3) k_g/k_s + 0.22 eb^2)] plus the dynamic
  // term |u| rho_g Cp_g d_p / Pe_H, Pe_H = 8.65 (1 + 19.4 (d_p/D_c)^2).
  // u is the superficial velocity [m/s], cp_g the gas heat capacity per unit
  // mass as the energy balance uses it.
  double DixonRadialConductivity(double k_g, double k_s, double eb, double dp, double D_c, double u,
                                 double rho_g, double cp_g);

  // Li and Finlayson (1977) wall Nusselt number h_w d_p / k_g of a packed
  // bed: 0.17 Re^0.79 for spheres, 0.16 Re^0.93 for cylinders (re >= 0).
  double LiFinlaysonWallNusselt(double re, bool cylinders);

  // Martin and Nilles (1993) wall Nusselt number h_w d_p / k_g:
  //     Nu_w = (1.3 + 5 d_p/D_c) k_0/k_g + 0.19 Re^(3/4) Pr^(1/3),
  // with k_0/k_g the stagnant bed conductivity ratio (Zehner-Schlunder).
  double MartinNillesWallNusselt(double re, double pr, double k0_over_kg, double dp_over_D);

  // Effective radial bed conductivity k_r [W/(m K)] of the VDI Heat Atlas
  // model that Martin and Nilles's h_w belongs to: the stagnant term plus the
  // radial dynamic term Pe_0/K_r (Bauer and Schlunder 1978),
  //     k_r = k_g [k_0/k_g + Re Pr / K_r],   K_r = 7 [2 - (1 - 2 d_p/D_c)^2].
  double ZehnerSchlunderRadialConductivity(double k_g, double k0_over_kg, double re, double pr, double dp_over_D);

  // Stored values; never renumber.
  enum class ExternalHeatTransferType {
    CONSTANT = 0,
    HOLMAN = 1,
    CHURCHILL_CHU_VERTICAL = 2,
    CHURCHILL_CHU_HORIZONTAL = 3,
    YANG_VERTICAL_CYLINDER = 4,
    CHURCHILL_BERNSTEIN = 5,
    MIXED_CONVECTION_VERTICAL = 6,
    MIXED_CONVECTION_HORIZONTAL = 7
  };

  static const char* ExternalHeatTransferModeList[] = {
    "Constant",
    "Holman (2008)",
    "Churchill and Chu (1975), vertical column",
    "Churchill and Chu (1975), horizontal column",
    "Yang (1985), vertical cylinder",
    "Churchill and Bernstein (1977), cylinder in crossflow",
    "Mixed convection (Churchill 1977), vertical column",
    "Mixed convection (Churchill 1977), horizontal column"
  };

  class IInternalHeatTransferModel {
  public:
    virtual ~IInternalHeatTransferModel() = default;
    virtual void Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall) = 0;
    virtual InternalHeatTransferType GetModelType() = 0;
    virtual unsigned ReadsTransportProperties() const { return kAllTransportProperties; }
    virtual void EvaluateExpressions(const EvalContext& ctx) {}

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {}
  };

  class IExternalHeatTransferModel {
  public:
    virtual ~IExternalHeatTransferModel() = default;
    // ambient_temperature: the temperature of the surroundings [K].
    virtual void Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall,
                        double ambient_temperature) = 0;
    virtual ExternalHeatTransferType GetModelType() = 0;
    virtual unsigned ReadsTransportProperties() const { return kAllTransportProperties; }
    virtual void EvaluateExpressions(const EvalContext& ctx) {}

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {}
  };

  class ConstantInternalHeatTransferModel : public IInternalHeatTransferModel {
  public:
    void Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall) override;
    InternalHeatTransferType GetModelType() override { return InternalHeatTransferType::CONSTANT; }
    unsigned ReadsTransportProperties() const override { return kNoTransportProperties; }
    void EvaluateExpressions(const EvalContext& ctx) override { heat_transfer_coefficient_.Evaluate(ctx); }

    ExpressionParam heat_transfer_coefficient_{0.0}; // [W/m2]

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::virtual_base_class<IInternalHeatTransferModel>(this),
        cereal::make_nvp("Gas_Wall_Heat_Transfer_Coefficient", heat_transfer_coefficient_)
      );
    }
  };

  class InternalHeatTransferModel_Dixon : public IInternalHeatTransferModel {
  public:
    void Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall) override;
    InternalHeatTransferType GetModelType() override { return InternalHeatTransferType::DIXON_1996; }
    unsigned ReadsTransportProperties() const override { return kReynoldsNumber | kGasThermalConductivity; }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::virtual_base_class<IInternalHeatTransferModel>(this)
      );
    }
  };

  // Li and Finlayson (1977) wall coefficient (LiFinlaysonWallNusselt, by the
  // layer's particle shape), lumped with Dixon's radial conductivity
  // (LumpedWallHeatTransferCoefficient). No parameters.
  class InternalHeatTransferModel_LiFinlayson : public IInternalHeatTransferModel {
  public:
    void Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall) override;
    InternalHeatTransferType GetModelType() override { return InternalHeatTransferType::LI_FINLAYSON_1977; }
    unsigned ReadsTransportProperties() const override { return kReynoldsNumber | kGasThermalConductivity; }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::virtual_base_class<IInternalHeatTransferModel>(this));
    }
  };

  // Martin and Nilles (1993) wall coefficient (MartinNillesWallNusselt, with
  // the Zehner-Schlunder stagnant conductivity), lumped with a radial
  // conductivity (LumpedWallHeatTransferCoefficient). No parameters.
  class InternalHeatTransferModel_MartinNilles : public IInternalHeatTransferModel {
  public:
    void Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall) override;
    InternalHeatTransferType GetModelType() override { return InternalHeatTransferType::MARTIN_NILLES_1993; }
    unsigned ReadsTransportProperties() const override { return kReynoldsNumber | kGasThermalConductivity; }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::virtual_base_class<IInternalHeatTransferModel>(this));
    }
  };

  // Leva (1947) packed-bed wall heat-transfer correlation, as used by Casas
  // et al. (2012) for the internal (gas-to-wall) coefficient h_L:
  //     Nu = h_L·D_t / k_g = eta1 · Re^{eta2} · exp(−exp_coeff · d_p / D_t)
  // where D_t is the column diameter and Re = ρ·u·d_p/μ (particle Reynolds).
  // eta1, eta2 default to Leva's canonical 0.813 and 0.90; Casas fitted them to
  // the reference breakthrough experiment (eta1 = 41.13, eta2 = 0.32).
  class InternalHeatTransferModel_Leva : public IInternalHeatTransferModel {
  public:
    void Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall) override;
    InternalHeatTransferType GetModelType() override { return InternalHeatTransferType::LEVA_1947; }
    unsigned ReadsTransportProperties() const override { return kReynoldsNumber | kGasThermalConductivity; }
    void EvaluateExpressions(const EvalContext& ctx) override {
      eta1.Evaluate(ctx); eta2.Evaluate(ctx); exp_coeff.Evaluate(ctx);
    }

    ExpressionParam eta1{0.813};     // Re coefficient (Leva canonical 0.813)
    ExpressionParam eta2{0.90};      // Re exponent    (Leva canonical 0.90)
    ExpressionParam exp_coeff{6.0};  // geometric decay coefficient exp(−exp_coeff·dp/D_t)

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::virtual_base_class<IInternalHeatTransferModel>(this),
        cereal::make_nvp("Leva_eta1", eta1),
        cereal::make_nvp("Leva_eta2", eta2),
        cereal::make_nvp("Leva_exp_coeff", exp_coeff)
      );
    }
  };

  class ConstantExternalHeatTransferModel : public IExternalHeatTransferModel {
  public:
    void Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall,
                double ambient_temperature) override;
    ExternalHeatTransferType GetModelType() override { return ExternalHeatTransferType::CONSTANT; }
    unsigned ReadsTransportProperties() const override { return kNoTransportProperties; }
    void EvaluateExpressions(const EvalContext& ctx) override { heat_transfer_coefficient_.Evaluate(ctx); }

    ExpressionParam heat_transfer_coefficient_{0.0};

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::virtual_base_class<IExternalHeatTransferModel>(this),
        cereal::make_nvp("Wall_Ambient_Heat_Transfer_Coefficient", heat_transfer_coefficient_)
      );
    }
  };

  // Forced convection, Nu = c Re^m Pr^n on the outer diameter, with the air
  // properties at the film temperature (ambient_air.h).
  class ExternalHeatTransferModel_Holman : public IExternalHeatTransferModel {
  public:
    void Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall,
                double ambient_temperature) override;
    ExternalHeatTransferType GetModelType() override { return ExternalHeatTransferType::HOLMAN; }
    unsigned ReadsTransportProperties() const override { return kNoTransportProperties; }
    void EvaluateExpressions(const EvalContext& ctx) override {
      c.Evaluate(ctx); m.Evaluate(ctx); n.Evaluate(ctx); external_velocity.Evaluate(ctx);
    }

    ExpressionParam c{0.027};             // coefficient
    ExpressionParam m{0.8};               // Power of Re
    ExpressionParam n{0.3333333333};      // Power of Pr
    ExpressionParam external_velocity{0.5}; // External air flow [m/s]

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::virtual_base_class<IExternalHeatTransferModel>(this),
        cereal::make_nvp("Coefficient", c),
        cereal::make_nvp("Power_Coefficient_Re", m),
        cereal::make_nvp("Power_Coefficient_Pr", n),
        cereal::make_nvp("External_Velocity", external_velocity)
      );
    }
  };

  // Natural convection to still air by Churchill and Chu (1975)
  // (ambient_air.h): a vertical column as a vertical plate whose height is the
  // column length, a horizontal column as a horizontal cylinder of the outer
  // wall diameter. Each cell uses its own wall temperature; the air
  // properties are at the film temperature. No parameters.
  class ExternalHeatTransferModel_ChurchillChu : public IExternalHeatTransferModel {
  public:
    explicit ExternalHeatTransferModel_ChurchillChu(bool vertical = true) : vertical_(vertical) {}
    void Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall,
                double ambient_temperature) override;
    ExternalHeatTransferType GetModelType() override {
      return vertical_ ? ExternalHeatTransferType::CHURCHILL_CHU_VERTICAL
                       : ExternalHeatTransferType::CHURCHILL_CHU_HORIZONTAL;
    }
    unsigned ReadsTransportProperties() const override { return kNoTransportProperties; }
    bool Vertical() const { return vertical_; }

  private:
    bool vertical_ = true;

    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      (void)version;
      ar(
        cereal::virtual_base_class<IExternalHeatTransferModel>(this),
        cereal::make_nvp("Vertical", vertical_)
      );
    }
  };

  // Natural convection from a vertical column to still air at any height-to-
  // diameter ratio, by Yang (1985) (ambient_air.h): the column is a vertical
  // cylinder of the bed height and the outer wall diameter. Each cell uses its
  // own wall temperature; the air properties are at the film temperature. No
  // parameters. Unlike churchill_chu_vertical (a plate), it holds for a
  // slender column.
  class ExternalHeatTransferModel_YangCylinder : public IExternalHeatTransferModel {
  public:
    void Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall,
                double ambient_temperature) override;
    ExternalHeatTransferType GetModelType() override { return ExternalHeatTransferType::YANG_VERTICAL_CYLINDER; }
    unsigned ReadsTransportProperties() const override { return kNoTransportProperties; }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      (void)version;
      ar(cereal::virtual_base_class<IExternalHeatTransferModel>(this));
    }
  };

  // Forced convection from a cylinder in crossflow by Churchill and Bernstein
  // (1977) (ambient_air.h): air at external_velocity crosses the column, a
  // cylinder of the outer wall diameter. Each cell uses its own wall
  // temperature; the air properties are at the film temperature. Valid for
  // Re Pr > 0.2.
  class ExternalHeatTransferModel_ChurchillBernstein : public IExternalHeatTransferModel {
  public:
    void Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall,
                double ambient_temperature) override;
    ExternalHeatTransferType GetModelType() override { return ExternalHeatTransferType::CHURCHILL_BERNSTEIN; }
    unsigned ReadsTransportProperties() const override { return kNoTransportProperties; }
    void EvaluateExpressions(const EvalContext& ctx) override { external_velocity.Evaluate(ctx); }

    ExpressionParam external_velocity{0.5}; // External air flow [m/s]

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      (void)version;
      ar(
        cereal::virtual_base_class<IExternalHeatTransferModel>(this),
        cereal::make_nvp("External_Velocity", external_velocity)
      );
    }
  };

  // Mixed convection: air at external_velocity crossing a column that also
  // loses heat by natural convection, h = (h_F^n + h_N^n)^(1/n)
  // (ambient_air.h, MixedConvectionCoefficient). h_N is Yang's vertical
  // cylinder of the bed height (vertical) or Churchill and Chu's horizontal
  // cylinder (horizontal); h_F is Churchill and Bernstein's forced term. Each
  // cell uses its own wall temperature; the air properties are at the film
  // temperature.
  class ExternalHeatTransferModel_MixedConvection : public IExternalHeatTransferModel {
  public:
    explicit ExternalHeatTransferModel_MixedConvection(bool vertical = true) : vertical_(vertical) {}
    void Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall,
                double ambient_temperature) override;
    ExternalHeatTransferType GetModelType() override {
      return vertical_ ? ExternalHeatTransferType::MIXED_CONVECTION_VERTICAL
                       : ExternalHeatTransferType::MIXED_CONVECTION_HORIZONTAL;
    }
    unsigned ReadsTransportProperties() const override { return kNoTransportProperties; }
    void EvaluateExpressions(const EvalContext& ctx) override {
      external_velocity.Evaluate(ctx); exponent.Evaluate(ctx);
    }
    bool Vertical() const { return vertical_; }

    ExpressionParam external_velocity{0.5}; // External air flow [m/s]
    ExpressionParam exponent{3.0};          // Blending exponent n [-]

  private:
    bool vertical_ = true;

    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      (void)version;
      ar(
        cereal::virtual_base_class<IExternalHeatTransferModel>(this),
        cereal::make_nvp("Vertical", vertical_),
        cereal::make_nvp("External_Velocity", external_velocity),
        cereal::make_nvp("Exponent", exponent)
      );
    }
  };

  // ---- Particle Heat Transfer (gas-to-solid) ----

  // Stored values; never renumber.
  enum class ParticleHeatTransferType {
    CONSTANT = 0,
    FILM = 1
  };

  static const char* ParticleHeatTransferModeList[] = {
    "Constant",
    "Film"
  };

  class IParticleHeatTransferModel {
  public:
    virtual ~IParticleHeatTransferModel() = default;
    virtual void Update(Equations& state, const MixedAdsorbentLayer& layer) = 0;
    virtual ParticleHeatTransferType GetModelType() = 0;
    virtual unsigned ReadsTransportProperties() const { return kAllTransportProperties; }
    virtual void EvaluateExpressions(const EvalContext& ctx) {}

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {}
  };

  class ConstantParticleHeatTransferModel : public IParticleHeatTransferModel {
  public:
    void Update(Equations& state, const MixedAdsorbentLayer& layer) override;
    ParticleHeatTransferType GetModelType() override { return ParticleHeatTransferType::CONSTANT; }
    unsigned ReadsTransportProperties() const override { return kNoTransportProperties; }
    void EvaluateExpressions(const EvalContext& ctx) override { heat_transfer_coefficient_.Evaluate(ctx); }

    ExpressionParam heat_transfer_coefficient_{100.0}; // [W/(m2.K)] gas-to-particle

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::virtual_base_class<IParticleHeatTransferModel>(this),
        cereal::make_nvp("Particle_Heat_Transfer_Coefficient", heat_transfer_coefficient_)
      );
    }
  };

  // Gas-solid film coefficient from correlations, h_p = h a_v [W/(m3 K)]:
  //   packed bed:  Nu = h d_p / k_g by the layer's film correlation
  //                (MixedAdsorbentLayer::film_correlation_, film_correlation.h;
  //                default Wakao-Kaguei-Funazkri, 2 + 1.1 Pr^(1/3) Re^0.6;
  //                Re on the superficial velocity and d_p);
  //   structured:  Nu = h D_h / k_g, the channel's fully developed Nu_T
  //                (Shah and London 1978), or Hawthorn's (1974) developing
  //                form when the layer's entrance effects are on
  //                (ChannelFilmNumber; Re on the channel velocity and D_h).
  // a_v is the layer's specific surface (domain vector av).
  class FilmParticleHeatTransferModel : public IParticleHeatTransferModel {
  public:
    void Update(Equations& state, const MixedAdsorbentLayer& layer) override;
    ParticleHeatTransferType GetModelType() override { return ParticleHeatTransferType::FILM; }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      (void)version;
      ar(cereal::virtual_base_class<IParticleHeatTransferModel>(this));
    }
  };

  // Factory Functions
  inline std::unique_ptr<IParticleHeatTransferModel> CreateParticleHeatTransferModel(ParticleHeatTransferType type) {
    switch (type) {
    case ParticleHeatTransferType::CONSTANT: return std::make_unique<ConstantParticleHeatTransferModel>();
    case ParticleHeatTransferType::FILM: return std::make_unique<FilmParticleHeatTransferModel>();
    default: return nullptr;
    }
  }

  inline std::unique_ptr<IInternalHeatTransferModel> CreateInternalHeatTransferModel(InternalHeatTransferType type) {
    switch (type) {
    case InternalHeatTransferType::CONSTANT: return std::make_unique<ConstantInternalHeatTransferModel>();
    case InternalHeatTransferType::DIXON_1996: return std::make_unique<InternalHeatTransferModel_Dixon>();
    case InternalHeatTransferType::LEVA_1947: return std::make_unique<InternalHeatTransferModel_Leva>();
    case InternalHeatTransferType::LI_FINLAYSON_1977: return std::make_unique<InternalHeatTransferModel_LiFinlayson>();
    case InternalHeatTransferType::MARTIN_NILLES_1993: return std::make_unique<InternalHeatTransferModel_MartinNilles>();
    default: return nullptr;
    }
  }
  
  inline std::unique_ptr<IExternalHeatTransferModel> CreateExternalHeatTransferModel(ExternalHeatTransferType type) {
    switch (type) {
    case ExternalHeatTransferType::CONSTANT: return std::make_unique<ConstantExternalHeatTransferModel>();
    case ExternalHeatTransferType::HOLMAN: return std::make_unique<ExternalHeatTransferModel_Holman>();
    case ExternalHeatTransferType::CHURCHILL_CHU_VERTICAL:
      return std::make_unique<ExternalHeatTransferModel_ChurchillChu>(true);
    case ExternalHeatTransferType::CHURCHILL_CHU_HORIZONTAL:
      return std::make_unique<ExternalHeatTransferModel_ChurchillChu>(false);
    case ExternalHeatTransferType::YANG_VERTICAL_CYLINDER:
      return std::make_unique<ExternalHeatTransferModel_YangCylinder>();
    case ExternalHeatTransferType::CHURCHILL_BERNSTEIN:
      return std::make_unique<ExternalHeatTransferModel_ChurchillBernstein>();
    case ExternalHeatTransferType::MIXED_CONVECTION_VERTICAL:
      return std::make_unique<ExternalHeatTransferModel_MixedConvection>(true);
    case ExternalHeatTransferType::MIXED_CONVECTION_HORIZONTAL:
      return std::make_unique<ExternalHeatTransferModel_MixedConvection>(false);
    default: return nullptr;
    }
  }

} // end adsim namespace

// Class versions live next to the inline serialize templates: every
// translation unit that saves or loads a deck must see the same numbers.
CEREAL_CLASS_VERSION(adsim::ConstantParticleHeatTransferModel, 1);
CEREAL_CLASS_VERSION(adsim::ConstantInternalHeatTransferModel, 1);
CEREAL_CLASS_VERSION(adsim::InternalHeatTransferModel_Dixon, 1);
CEREAL_CLASS_VERSION(adsim::InternalHeatTransferModel_Leva, 1);
CEREAL_CLASS_VERSION(adsim::InternalHeatTransferModel_LiFinlayson, 1);
CEREAL_CLASS_VERSION(adsim::InternalHeatTransferModel_MartinNilles, 1);
CEREAL_CLASS_VERSION(adsim::ConstantExternalHeatTransferModel, 1);
CEREAL_CLASS_VERSION(adsim::ExternalHeatTransferModel_Holman, 1);
CEREAL_CLASS_VERSION(adsim::ExternalHeatTransferModel_ChurchillChu, 1);
CEREAL_CLASS_VERSION(adsim::ExternalHeatTransferModel_YangCylinder, 1);
CEREAL_CLASS_VERSION(adsim::ExternalHeatTransferModel_ChurchillBernstein, 1);
CEREAL_CLASS_VERSION(adsim::ExternalHeatTransferModel_MixedConvection, 1);
CEREAL_CLASS_VERSION(adsim::FilmParticleHeatTransferModel, 1);
