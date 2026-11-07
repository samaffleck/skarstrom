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

  enum class InternalHeatTransferType {
    CONSTANT,
    DIXON_1996,
    LEVA_1947
  };

  static const char* InternalHeatTransferModeList[] = {
    "Constant", // 0
    "Dixon (1996)", // 1
    "Leva (1947)" // 2
  };

  enum class ExternalHeatTransferType {
    CONSTANT,
    HOLMAN
  };

  static const char* ExternalHeatTransferModeList[] = {
    "Constant",     
    "Holman (2008)"  
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
    virtual void Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall) = 0;
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
    void Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall) override;
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

  class ExternalHeatTransferModel_Holman : public IExternalHeatTransferModel {
  public:
    void Update(Equations& state, const MixedAdsorbentLayer& layer, const BedWall& wall) override;
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
  //   packed bed:  Nu = h d_p / k_g = 2 + 1.1 Pr^(1/3) Re^0.6 (Wakao and
  //                Funazkri 1978; Re on the superficial velocity and d_p);
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
    default: return nullptr;
    }
  }
  
  inline std::unique_ptr<IExternalHeatTransferModel> CreateExternalHeatTransferModel(ExternalHeatTransferType type) {
    switch (type) {
    case ExternalHeatTransferType::CONSTANT: return std::make_unique<ConstantExternalHeatTransferModel>();
    case ExternalHeatTransferType::HOLMAN: return std::make_unique<ExternalHeatTransferModel_Holman>();
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
CEREAL_CLASS_VERSION(adsim::ConstantExternalHeatTransferModel, 1);
CEREAL_CLASS_VERSION(adsim::ExternalHeatTransferModel_Holman, 1);
CEREAL_CLASS_VERSION(adsim::FilmParticleHeatTransferModel, 1);
