/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <atomic>
#include <limits>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <map>
#include <memory>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <variant>

// Eigen includes
#include <Eigen/Dense>

// Cereal includes
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>

// ads includes
#include "adsim/math/sorted_map_archive.h"
#include "adsim/equations/film_correlation.h"
#include "adsim/objects/fluid.h"
#include "adsim/equations/isotherm.h"
#include "adsim/equations/iast.h"
#include "adsim/equations/mass_transfer_model.h"
#include "adsim/objects/reactor_equations.h"
#include "adsim/objects/reaction.h"
#include "adsim/objects/domain.h"
#include "adsim/equations/equation.h"
#include "adsim/math/expression_param.h"


namespace adsim {

  class MixedAdsorbentLayer;

  // How one adsorbent combines its components' isotherms into mixture loadings.
  //   None: each component's isotherm is evaluated on its own at the partial
  //         pressure (the extended kinds couple components through their own
  //         shared denominators).
  //   IAST: Ideal Adsorbed Solution Theory over the components' pure isotherms
  //         (iast.h); every component keeps its own kind and parameters, and an
  //         extended kind is taken in its pure-component limit.
  // Enumerator values are what the deck stores, so never renumber them.
  enum class MixingRule {
    None = 0,
    IAST = 1
  };

  static const char* MixingRuleList[] = {
    "None",  // 0
    "IAST"   // 1
  };

  class AdsorptionSource {
  public:
    AdsorptionSource() = default;
    AdsorptionSource(AdsorptionSource&&) = default;
    AdsorptionSource& operator=(AdsorptionSource&&) = default;
    AdsorptionSource(const AdsorptionSource&) = delete;
    AdsorptionSource& operator=(const AdsorptionSource&) = delete;

    // Seeds defaults for every fluid component and purges per-component entries
    // (isotherm params, isotherm map, heat of adsorption, rate model) for any
    // component the fluid does not contain.
    void Initialise(const Fluid& fluid);
    // Erases every per-component entry for `component`. Never inserts.
    void RemoveComponent(const std::string& component);
    void SetIsothermType(const std::string& component, IsothermType isotherm_type);
    void AddIsotherm(const std::vector<std::string>& components, std::unique_ptr<IIsotherm> isotherm);
    // Give `components` one isotherm of `type` (created with default
    // parameters if this source has none of that type yet). The components
    // leave their previous isotherm, and an isotherm no component uses any
    // more is removed.
    void AssignIsotherm(const std::vector<std::string>& components, IsothermType type);
    // Which isotherm type each component uses.
    const std::unordered_map<std::string, IsothermType>& IsothermMap() const { return isotherm_map_; }
    void UpdateIsotherm(Equations& state, AdsorbentID ads_id, const std::string& component, int start_index, int end_index, const UserVariableRegistry* registry);
    // Update every component's isotherm in one pass for a single adsorbent.
    // Components sharing a competitive (extended) isotherm reuse its per-cell
    // denominator instead of recomputing it once per component. Under
    // MixingRule::IAST the loadings are the IAST solution over the components'
    // pure isotherms instead (one solve per cell, warm-started from the cell's
    // last solve). Writes into the adsorbent's own equilibrium-loading map
    // (eq_ads_concentration[ads_id]).
    void UpdateIsothermAll(Equations& state, AdsorbentID ads_id, const std::vector<std::string>& components, int start_index, int end_index, const UserVariableRegistry* registry);
    // solid_density [kg/m3] is the owning adsorbent's particle density, which
    // the capacity-factor coefficient needs (0 when unknown).
    void UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const std::string& component, const Fluid& fluid, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry, double solid_density = 0.0);
    // The slow-site coefficient of a dual-rate LDF component, into
    // state.slow_mass_transfer_coefficient (no-op without a slow site).
    void UpdateSlowMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const std::string& component, const Fluid& fluid, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry, double solid_density = 0.0);
    std::unique_ptr<IIsotherm>& GetIsotherm(const std::string& component);

    std::unordered_map<std::string, std::unique_ptr<IMassTransferModel>> mass_transfer_models_{};
    std::unordered_map<std::string, ExpressionParam> heat_of_adsorption_{}; // Heat of adsorption [J/mol]
    MixingRule mixing_rule_ = MixingRule::None;

    void EvaluateExpressions(const EvalContext& ctx) {
      for (auto& [c, mtm] : mass_transfer_models_) if (mtm) mtm->EvaluateExpressions(ctx);
      for (auto& [t, iso] : isotherms_) if (iso) iso->EvaluateExpressions(ctx);
      for (auto& [c, h] : heat_of_adsorption_) h.Evaluate(ctx);
    }
    void CompileForRun(const UserVariableRegistry* registry) {
      for (auto& [c, mtm] : mass_transfer_models_) if (mtm) mtm->CompileForRun(registry);
      // Every run starts the IAST iterations cold, so a rerun of the same deck
      // reproduces the first run bit for bit.
      iast_psi_.clear();
      iast_scratch_ = IastScratch{};
    }

  private:
    std::unique_ptr<IIsotherm>& GetIsotherm(IsothermType type, const std::string& component);
    void UpdateIsothermIast(Equations& state, AdsorbentID ads_id, const std::vector<std::string>& components, int start_index, int end_index);

    std::unordered_map<IsothermType, std::unique_ptr<IIsotherm>> isotherms_{};
    std::unordered_map<std::string, IsothermType> isotherm_map_{};

    // IAST working storage (not serialised): the resolved pure isotherms, the
    // per-cell mixture state and the reduced spreading pressure of each cell's
    // last solve, the warm start of the next one.
    std::vector<PureIsothermRef> iast_refs_{};
    std::vector<PureIsotherm> iast_pure_{};
    std::vector<double> iast_y_{}, iast_q_{}, iast_psi_{};
    std::vector<const std::vector<double>*> iast_y_of_{};
    std::vector<std::vector<double>*> iast_q_of_{};
    IastScratch iast_scratch_{};

    friend class cereal::access;
    template<class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      SerializeMap(ar, "Mass_Transfer_Models", mass_transfer_models_);
      SerializeMap(ar, "Heat_Of_Adsorption", heat_of_adsorption_);
      SerializeMap(ar, "Isotherms", isotherms_);
      SerializeMap(ar, "Isotherm_Map", isotherm_map_);

      // Version 1 adds the mixing rule; older archives mean None.
      if (version >= 1) {
        ar(cereal::make_nvp("Mixing_Rule", mixing_rule_));
      } else {
        mixing_rule_ = MixingRule::None;
      }
    }
  };

  // Momentum balance (pressure drop vs superficial velocity u) of a packed
  // layer. The face velocity solves  -dP/dz = A u + B |u| u  with
  //   KarmenCozeny: A = mu / K,  K = eb^3 dp^2 / (180 (1 - eb)^2),  B = 0
  //                 (Darcy's law with the Kozeny-Carman permeability; the
  //                 legacy model and the one every pre-v3 archive loads as);
  //   Ergun:        A = 150 mu (1 - eb)^2 / (eb^3 dp^2),
  //                 B = 1.75 rho (1 - eb) / (eb^3 dp)
  //                 (viscous Blake-Kozeny + inertial Burke-Plummer terms; the
  //                 PSA-literature standard, default for new layers);
  //   BlakeKozeny:  A = mu / K,  K = eb^3 dp^2 / (150 (1 - eb)^2),  B = 0
  //                 (Darcy's law with the Blake-Kozeny permeability: Ergun's
  //                 viscous term without the inertial one, as used by
  //                 Haghpanah et al. 2013 and Estupinan Perez et al. 2019).
  // Enumerator values are what the deck stores, so never renumber them; new
  // models are appended.
  enum class MomentumModel {
    KarmenCozeny = 0,
    Ergun = 1,
    BlakeKozeny = 2
  };
  // Legacy name: the selector used to choose only the permeability.
  using PermeabilityModel = MomentumModel;

  static const char* MomentumModelList[] = {
    "Darcy (Kozeny-Carman)",  // 0
    "Ergun",                  // 1
    "Darcy (Blake-Kozeny)"    // 2
  };

  class Adsorbent {
  public:
    Adsorbent() = default;
    Adsorbent(const std::string& name);
    ~Adsorbent() = default;
    Adsorbent(Adsorbent&&) = default;
    Adsorbent& operator=(Adsorbent&&) = default;
    Adsorbent(const Adsorbent&) = delete;
    Adsorbent& operator=(const Adsorbent&) = delete;

    std::string name_ = "Default Solid";
    ExpressionParam solid_density_{1200.0};         // Solid density [kg/m3]
    ExpressionParam heat_capacity_{500.0};          // Specific heat capacity [J/(kg-K)]
    ExpressionParam thermal_conductivity_{0.5};     // Thermal conductivity [W/(m-K)]
    ExpressionParam ep_{0.35};                      // Intra-particle void fraction [m3/m3]
    AdsorptionSource adsorption_source_{};
    // Heterogeneous reactions this solid catalyses (objects/reaction.h); empty
    // for a plain sorbent. Rates are per kg of this solid.
    ReactionSource reaction_source_{};

    // Defined out-of-line in porous_media.cpp so it can cascade into
    // adsorption_source_ (isotherms + MTC models) once those are migrated.
    void EvaluateExpressions(const EvalContext& ctx);

  private:
    friend class cereal::access;
    template<class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::make_nvp("Name", name_),
        cereal::make_nvp("Density", solid_density_),
        cereal::make_nvp("Specific_Heat_Capacity", heat_capacity_),
        cereal::make_nvp("Thermal_Conductivity", thermal_conductivity_),
        cereal::make_nvp("Adsorption_Source", adsorption_source_),
        cereal::make_nvp("Particle_Void_Fraction", ep_)
      );
      // Version 5 adds the reactions; older archives have none.
      if (version >= 5) {
        ar(cereal::make_nvp("Reaction_Source", reaction_source_));
      } else if constexpr (Archive::is_loading::value) {
        reaction_source_ = ReactionSource{};
      }
    }
  };

  using AdsorbentID = size_t;

  // Shared material library - stores adsorbent definitions by ID
  class AdsorbentDirectory {
  public:
    AdsorbentID Add(const std::string& name = "New Adsorbent") {
      ++next_id_;
      adsorbents_.emplace(next_id_, Adsorbent(name));
      return next_id_;
    }

    // Insert an adsorbent under a caller-chosen id (used when reconstructing a
    // directory from an external document that owns the ids). Keeps next_id_
    // ahead of every explicit id so subsequent Add() calls never collide.
    AdsorbentID AddWithId(AdsorbentID id, Adsorbent&& adsorbent) {
      adsorbents_.insert_or_assign(id, std::move(adsorbent));
      if (id > next_id_) next_id_ = id;
      return id;
    }

    void Remove(AdsorbentID id) {
      adsorbents_.erase(id);
    }

    Adsorbent& Get(AdsorbentID id) {
      return adsorbents_.at(id);
    }

    const Adsorbent& Get(AdsorbentID id) const {
      return adsorbents_.at(id);
    }

    const std::unordered_map<AdsorbentID, Adsorbent>& GetAll() const {
      return adsorbents_;
    }

    bool Contains(AdsorbentID id) const {
      return adsorbents_.find(id) != adsorbents_.end();
    }

    bool Empty() const {
      return adsorbents_.empty();
    }

    size_t Size() const {
      return adsorbents_.size();
    }

  private:
    std::unordered_map<AdsorbentID, Adsorbent> adsorbents_;
    AdsorbentID next_id_ = 0;

    friend class cereal::access;
    template<class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      SerializeMap(ar, "Adsorbents", adsorbents_);
      ar(cereal::make_nvp("Next_ID", next_id_));
    }
  };

  using LayerID = size_t;

  // What a layer is packed with (a geometry option of MixedAdsorbentLayer):
  //   Packed:   a bed of particles (bed void fraction, particle diameter and
  //             a momentum model).
  //   Monolith: straight parallel channels in a honeycomb (cell density, wall
  //             thickness, channel shape, washcoat).
  //   Laminate: parallel sheets with a gap between them (sheet thickness,
  //             gap, washcoat).
  // Enumerator values are what the deck stores, so never renumber them.
  enum class LayerKind {
    Packed = 0,
    Monolith = 1,
    Laminate = 2
  };

  static const char* LayerKindList[] = {
    "Packed bed",  // 0
    "Monolith",    // 1
    "Laminate"     // 2
  };

  // Cross-section of a monolith channel. Stored; never renumber.
  enum class ChannelShape {
    Square = 0,
    Circular = 1
  };

  static const char* ChannelShapeList[] = {
    "Square",    // 0
    "Circular"   // 1
  };

  // Shape of the particles of a packed layer. The stored size is the particle
  // diameter in both cases (MixedAdsorbentLayer::particle_diameter_):
  //   Sphere:   spheres of diameter d_p (default; archives older than
  //             MixedAdsorbentLayer v7).
  //   Cylinder: extrudates or fibres of diameter d_p = 2 r_f, taken as
  //             infinitely long in the diffusion geometry (the end faces are
  //             not counted).
  // The LDF diffusion geometry (GetDiffusionGeometry) and the specific surface
  // area (ParticleShapeSurfaceFactor) follow the shape; the pressure drop,
  // dispersion and film correlations use d_p as an equivalent-sphere diameter.
  // Enumerator values are what the deck stores, so never renumber them.
  enum class ParticleShape {
    Sphere = 0,
    Cylinder = 1
  };

  static const char* ParticleShapeList[] = {
    "Sphere",    // 0
    "Cylinder"   // 1
  };

  // Factor c of the packed specific surface area a_v = c (1 - eb) / d_p:
  // 6 for spheres, 4 for infinite cylinders of diameter d_p. The solver
  // (domain.cpp) and MixedAdsorbentLayer::GetSpecificSurfaceArea share it.
  inline double ParticleShapeSurfaceFactor(ParticleShape shape) {
    return shape == ParticleShape::Cylinder ? 4.0 : 6.0;
  }

  // Kind tags of ParticleShape, by enumerator value (registry, deck view and
  // deck edits).
  inline const std::vector<std::string>& ParticleShapeKinds() {
    static const std::vector<std::string> k{"sphere", "cylinder"};
    return k;
  }

  // Geometry of a structured layer. Every field is kept whatever the kind,
  // so switching kind loses nothing; only the fields of the layer's kind are
  // used (MixedAdsorbentLayer's geometry getters).
  class StructuredGeometry {
  public:
    ChannelShape channel_shape_ = ChannelShape::Square;   // monolith
    ExpressionParam cell_density_{6.2e5};       // monolith channels per cross-section [1/m2] (400 cpsi)
    ExpressionParam wall_thickness_{1.65e-4};   // monolith wall between channels [m]
    ExpressionParam sheet_thickness_{2e-4};     // laminate sheet [m]
    ExpressionParam gap_{2e-4};                 // laminate gap between sheets [m]
    ExpressionParam washcoat_thickness_{0.0};   // adsorbent layer on the channel walls [m]; 0 = the wall is the adsorbent
    // Film coefficients: false = fully developed laminar Nu and Sh (Shah and
    // London); true = Hawthorn's length-averaged developing-flow correction
    // over the layer length (MixedAdsorbentLayer::GetChannelNusselt).
    bool entrance_effects_ = false;

    void EvaluateExpressions(const EvalContext& ctx) {
      cell_density_.Evaluate(ctx);
      wall_thickness_.Evaluate(ctx);
      sheet_thickness_.Evaluate(ctx);
      gap_.Evaluate(ctx);
      washcoat_thickness_.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template<class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::make_nvp("Channel_Shape", channel_shape_),
        cereal::make_nvp("Cell_Density", cell_density_),
        cereal::make_nvp("Wall_Thickness", wall_thickness_),
        cereal::make_nvp("Sheet_Thickness", sheet_thickness_),
        cereal::make_nvp("Gap", gap_),
        cereal::make_nvp("Washcoat_Thickness", washcoat_thickness_)
      );
      // v2: the developing-flow switch (older archives: fully developed).
      if (version >= 2) ar(cereal::make_nvp("Entrance_Effects", entrance_effects_));
      else entrance_effects_ = false;
    }
  };

  // Read-only geometry of a layer for front ends (LayerGeometrySummary).
  struct LayerGeometry {
    LayerKind kind = LayerKind::Packed;
    double open_fraction = 0.0;           // gas fraction of the cross-section [-]
    double hydraulic_diameter = 0.0;      // [m]
    double specific_surface_area = 0.0;   // gas-solid area per bed volume [1/m]
    double diffusion_length = 0.0;        // [m]: washcoat thickness, else half the wall/sheet (packed: dp / 2)
    double channel_count = std::numeric_limits<double>::quiet_NaN();  // monolith channels across the column; NaN otherwise
    double washcoat_fraction = 0.0;       // washcoat share of the solid cross-section [-]
    double friction_factor_reynolds = 0.0;  // laminar Fanning f Re of the channel; 0 for packed
  };

  class MixedAdsorbentLayer {
  public:
    MixedAdsorbentLayer() = default;
    MixedAdsorbentLayer(const std::string& layer_name) : name_(layer_name) {}
    ~MixedAdsorbentLayer() = default;
    MixedAdsorbentLayer(MixedAdsorbentLayer&&) = default;
    MixedAdsorbentLayer& operator=(MixedAdsorbentLayer&&) = default;

    std::string name_ = "Layer 1";
    Domain domain_{};
    ExpressionParam eb_{0.37};                  // Inter-particle bed voidage [m3/m3]
    ExpressionParam particle_diameter_{1e-3};   // Particle diameter [m] (a cylinder's: 2 r_f)
    // Shape of the packed particles (ParticleShape). Not used by a structured
    // layer; kept for switching back to packed.
    ParticleShape particle_shape_ = ParticleShape::Sphere;
    // Default for a newly constructed layer is Ergun; archives older than
    // MixedAdsorbentLayer v3 load as KarmenCozeny (see serialize).
    MomentumModel momentum_model_ = MomentumModel::Ergun;
    LayerEquations equations_{};
    // Mass fraction w_k of each solid. Each is an expression (a number, or a
    // formula over the deck variables such as "x" and "1 - x" for a two-solid
    // blend), evaluated by EvaluateExpressions; every reader uses the value.
    std::map<AdsorbentID, ExpressionParam> composition_{};
    // Packed (default; archives older than v4) or a structured kind, whose
    // geometry is structured_. For a structured kind momentum_model_ is not
    // used: the channels give laminar Darcy flow (GetPermeability).
    LayerKind kind_ = LayerKind::Packed;
    StructuredGeometry structured_{};
    // Gas-particle film correlation of a packed layer (film_correlation.h):
    // the Nu of the film particle heat-transfer model and the Sh of the
    // Haynes k_f. Not used by a structured layer (its channel Nu = Sh).
    // Archives older than v6 load the default, Wakao-Funazkri.
    FilmCorrelation film_correlation_ = FilmCorrelation::WakaoFunazkri;

    bool IsStructured() const { return kind_ != LayerKind::Packed; }

    void EvaluateExpressions(const EvalContext& ctx) {
      eb_.Evaluate(ctx);
      particle_diameter_.Evaluate(ctx);
      domain_.EvaluateExpressions(ctx);
      equations_.EvaluateExpressions(ctx);
      structured_.EvaluateExpressions(ctx);
      for (auto& [id, w] : composition_) w.Evaluate(ctx);
    }

    // Composition management
    void AddAdsorbent(AdsorbentID id, double mass_fraction = 1.0);
    void RemoveAdsorbent(AdsorbentID id);
    // True when the composition is usable: not empty, no negative fraction,
    // and the fractions sum to 1 within kCompositionSumTolerance (the error
    // rules of CheckLayerComposition that need no adsorbent directory).
    bool ValidateComposition() const;
    // The evaluated mass fractions (w_k by adsorbent id).
    std::map<AdsorbentID, double> GetMassFractions() const;
    // True when a fraction is a formula rather than a plain number (the deck
    // edits that rescale the others refuse such a layer: deck_edit.h).
    bool HasCompositionExpression() const;

    // Computed properties
    // The gas fraction of the bed outside the solid: eb_ for a packed bed,
    // the open frontal area of the channels for a structured layer. Every use
    // of the bed void fraction (gas holdup, bulk density, the dispersion area,
    // the gas/solid split of heat) reads this.
    double GetVoidFraction() const { return IsStructured() ? GetOpenFraction() : eb_.Get(); }

    double GetTotalVoidFraction(const AdsorbentDirectory& dir) const {
      const double e = GetVoidFraction();
      return e + (1 - e) * GetEffectiveParticleVoidFraction(dir);
    }

    // Layer bulk density (1 - eb) * rho_s [kg adsorbent / m3 bed]; equals the
    // sum of GetAdsorbentBulkDensities().
    double GetBulkDensity(const AdsorbentDirectory& dir) const {
      return (1 - GetVoidFraction()) * GetEffectiveSolidDensity(dir);
    }

    // Viscous permeability K [m2] (A = mu / K) and inertial coefficient
    // beta [1/m] (B = beta * rho) of the layer's momentum model. A structured
    // layer: K = OFA D_h^2 / (2 fRe) (laminar channel flow on the superficial
    // velocity) and beta = 0, whatever momentum_model_.
    double GetPermeability() const;
    double GetInertialCoefficient() const;

    // Geometry. Packed: OFA = eb, D_h = (2/3) eb dp / (1 - eb), a_v =
    // 6 (1 - eb) / dp. Structured, with pitch p = 1 / sqrt(cell density) and
    // open width a = p - w - 2 t_wc (monolith):
    //   square:   OFA = a^2 / p^2,          D_h = a,              a_v = 4 a / p^2,  fRe = 14.23
    //   circular: OFA = pi a^2 / (4 p^2),   D_h = a,              a_v = pi a / p^2, fRe = 16
    //   laminate: OFA = (g - 2 t_wc)/(g + t), D_h = 2 (g - 2 t_wc), a_v = 2 / (g + t), fRe = 24
    // Values outside the valid range (a <= 0, ...) are returned as computed;
    // CheckLayerGeometry reports them.
    double GetOpenFraction() const;
    double GetHydraulicDiameter() const;
    double GetSpecificSurfaceArea() const;
    // Laminar Fanning friction factor times Re of a channel (0 for packed).
    double GetFrictionFactorReynolds() const;
    // Diffusion length into the solid [m]: the washcoat thickness, or half
    // the wall (monolith) or sheet (laminate) without one; dp / 2 if packed.
    double GetDiffusionLength() const;
    // Washcoat share of the solid cross-section [-] (0 for packed or no washcoat).
    double GetWashcoatFraction() const;

    // Fully developed laminar Nusselt number of a channel at constant wall
    // temperature, Nu_T (Shah and London 1978): 2.976 square, 3.657 circle,
    // 7.541 parallel plates; also the Sherwood number by the heat-mass
    // analogy. 0 for a packed bed (its correlations depend on Re).
    double GetChannelNusselt() const;
    // Taylor-Aris coefficient kappa of D_ax = D_m + u_ch^2 D_h^2 / (kappa D_m)
    // on the channel's own D_h: 192 for a circular tube (Taylor 1953, Aris
    // 1956), also used for square channels (an approximation); 840 for
    // parallel plates (Aris 1956: 210 on the gap width, D_h = 2 x gap).
    // 0 for a packed bed.
    double GetTaylorArisKappa() const;

    // Diffusion geometry of the linear-driving-force rate coefficients
    // k = c D_eff / l^2 (Glueckauf 1955): a packed bed of spheres c = 60 on
    // l = d_p (15 on the radius); a packed bed of infinite cylinders c = 8 on
    // the radius l = r_f = d_p / 2; a slab (monolith wall, laminate sheet or
    // washcoat) c = 3 on its half-thickness l. volume_to_surface is the
    // solid's V/S for the film resistance: r_p / 3 for a sphere, r_f / 2 for
    // a cylinder, l for a slab. `shape` says which; every shape but the
    // sphere has c on the radius or half-thickness `length`.
    enum class Shape { Sphere, Cylinder, Slab };
    struct DiffusionGeometry {
      double c = 60.0;
      double length = 0.0;             // [m]
      double volume_to_surface = 0.0;  // [m]
      Shape shape = Shape::Sphere;
    };
    DiffusionGeometry GetDiffusionGeometry() const;

    // Effective material properties. composition_ holds mass fractions w_k:
    // rho_s = 1 / sum_k(w_k / rho_k); eps_p and k_s are weighted by the solid
    // volume fractions phi_k = w_k * rho_s / rho_k; Cp_s = sum_k w_k Cp_k.
    double GetEffectiveParticleVoidFraction(const AdsorbentDirectory& dir) const;
    double GetEffectiveSolidDensity(const AdsorbentDirectory& dir) const;
    double GetEffectiveHeatCapacity(const AdsorbentDirectory& dir) const;
    double GetEffectiveThermalConductivity(const AdsorbentDirectory& dir) const;

    // Solid-volume fraction phi_k of each adsorbent (sums to 1).
    std::map<AdsorbentID, double> GetSolidVolumeFractions(const AdsorbentDirectory& dir) const;
    // Per-adsorbent bulk density rho_b,k = (1 - eb) * rho_mix * w_k [kg of k / m3 bed].
    std::map<AdsorbentID, double> GetAdsorbentBulkDensities(const AdsorbentDirectory& dir) const;
    // Mixture particle density 1 / sum_k(w_k / rho_k) with w normalised to sum 1.
    double GetMixtureSolidDensity(const AdsorbentDirectory& dir) const;
    // W = sum_k w_k (1 for a valid composition; not normalised).
    double GetTotalMassFraction() const;

    void Initialise(const AdsorbentDirectory& dir, const Fluid& fluid,
                    int start_index, int end_index);

    LayerID id_{};  // Set by AdsorbentLayers::Add()

  private:
    friend class cereal::access;
    template<class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::make_nvp("Name", name_),
        cereal::make_nvp("Domain", domain_),
        cereal::make_nvp("Bed_Void_Fraction", eb_),
        cereal::make_nvp("Particle_Diameter", particle_diameter_)
      );
      if (version >= 3) {
        ar(cereal::make_nvp("Momentum_Model", momentum_model_));
      } else {
        // Pre-v3 archives could only hold Darcy / Kozeny-Carman (the field had
        // a single enumerator), so they keep that behaviour on load.
        MomentumModel legacy = MomentumModel::KarmenCozeny;
        ar(cereal::make_nvp("Permeability_Model", legacy));
        momentum_model_ = MomentumModel::KarmenCozeny;
      }
      ar(cereal::make_nvp("Layer_Equations", equations_));
      // v5: each mass fraction is an expression. Older archives hold plain
      // numbers, which load as constant expressions with the same value.
      if (version >= 5) {
        ar(cereal::make_nvp("Composition", composition_));
      } else {
        std::map<AdsorbentID, double> numbers;
        ar(cereal::make_nvp("Composition", numbers));
        composition_.clear();
        for (const auto& [id, w] : numbers) composition_[id] = w;
      }
      ar(cereal::make_nvp("ID", id_));
      // v4: the layer kind and its structured geometry (always written).
      // Older archives are packed beds with the default geometry.
      if (version >= 4) {
        ar(
          cereal::make_nvp("Kind", kind_),
          cereal::make_nvp("Structured_Geometry", structured_)
        );
      } else {
        kind_ = LayerKind::Packed;
      }
      // v6: the film correlation (Wakao-Funazkri, the fixed Sh before v6).
      if (version >= 6) ar(cereal::make_nvp("Film_Correlation", film_correlation_));
      else film_correlation_ = FilmCorrelation::WakaoFunazkri;
      // v7: the particle shape (spheres before v7).
      if (version >= 7) ar(cereal::make_nvp("Particle_Shape", particle_shape_));
      else particle_shape_ = ParticleShape::Sphere;
    }
  };

  // Film Nusselt or Sherwood number of a laminar channel. L <= 0: the fully
  // developed value nu_inf (Shah and London 1978, Nu_T). L > 0: Hawthorn's
  // (1974) length-averaged developing-flow correction over a channel of
  // length L, nu_inf (1 + 0.095 Re Pr D_h / L)^0.45, with Pr (heat) or Sc
  // (mass) as `pr` and Re on the channel velocity and D_h.
  inline double ChannelFilmNumber(double nu_inf, double re, double pr, double dh, double length) {
    if (!(length > 0.0)) return nu_inf;
    return nu_inf * std::pow(1.0 + 0.095 * re * pr * dh / length, 0.45);
  }

  // The geometry of a layer for front ends; `cross_section_area` [m2] (the
  // column's, Bed::GetArea) gives the monolith channel count.
  LayerGeometry LayerGeometrySummary(const MixedAdsorbentLayer& layer, double cross_section_area);

  // Errors of a layer's geometry (structured kinds; none for packed): open
  // fraction outside (0, 1), a non-positive hydraulic diameter or diffusion
  // length (the walls and washcoat leave no channel: w + 2 t_wc >= p, or
  // 2 t_wc >= gap), or a non-positive input. Each message starts "Error: ".
  std::vector<std::string> CheckLayerGeometry(const MixedAdsorbentLayer& layer);

  // |sum_k w_k - 1| above which a layer composition is an error.
  inline constexpr double kCompositionSumTolerance = 1e-6;

  // One finding of CheckLayerComposition. Errors stop a run (Bed::PerformChecks
  // reports them); warnings do not (LayerCompositionWarnings, deck_warnings.h).
  struct LayerCompositionIssue {
    enum class Severity { Error, Warning };
    Severity severity = Severity::Error;
    std::string message;  // starts with "Error: " or "Warning: "
  };

  // The rules for a layer's composition (mass fractions w_k), in one place for
  // the deck checks, the run-start warnings and front ends:
  //   Error:   no adsorbent; an id the directory does not hold; a negative
  //            fraction; fractions that do not sum to 1 (beyond
  //            kCompositionSumTolerance). Every effective property scales with
  //            the sum (porous_media.cpp), so a sum of 0.8 would silently mean
  //            80 % of the solid; the fractions are never renormalised.
  //   Warning: a fraction of exactly 0 (that solid takes no part in the run).
  // A solid whose components all use the Inert isotherm is a valid member (an
  // inert diluent or packing). Issues come in the order of the composition map.
  std::vector<LayerCompositionIssue> CheckLayerComposition(const MixedAdsorbentLayer& layer,
                                                           const AdsorbentDirectory& dir);

  // Ordered collection of bed layers (replaces old Adsorbents class)
  class AdsorbentLayers {
  public:
    LayerID Add(const std::string& layer_name = "New Layer");
    void Remove(LayerID id);
    bool Swap(LayerID id1, LayerID id2);

    // Keep next_layer_id_ ahead of an externally-assigned layer id (used when
    // reconstructing layers from a document that owns the ids), so subsequent
    // Add() calls never mint a duplicate.
    void EnsureNextIdAtLeast(LayerID id) {
      if (id > next_layer_id_) next_layer_id_ = id;
    }
    MixedAdsorbentLayer& Get(LayerID id);
    const MixedAdsorbentLayer& Get(LayerID id) const;

    std::vector<MixedAdsorbentLayer> layers_{};  // ordered bottom-to-top

  private:
    LayerID next_layer_id_ = 0;

    friend class cereal::access;
    template<class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::make_nvp("Layers", layers_),
        cereal::make_nvp("Next_Layer_ID", next_layer_id_)
      );
    }
  };

  // Viscous permeability K [m2] of a packed layer: Kozeny-Carman (constant
  // 180) for KarmenCozeny, the Ergun/Blake-Kozeny viscous term (constant 150)
  // for Ergun and BlakeKozeny.
  double GetLayerPermeability(double eb, double dp, MomentumModel model);
  // Inertial (Burke-Plummer) coefficient beta [1/m] such that the inertial
  // pressure gradient is beta * rho * |u| * u: 1.75 (1 - eb) / (eb^3 dp) for
  // Ergun, 0 for KarmenCozeny and BlakeKozeny (which therefore take the
  // linear Darcy face path everywhere).
  double GetLayerInertialCoefficient(double eb, double dp, MomentumModel model);

  // Superficial velocity u solving  g = A u + B |u| u  for a driving pressure
  // drop g (P_upstream - P_downstream, or -dP/dz), with A > 0 the viscous and
  // B >= 0 the inertial coefficient. Closed form
  //   u = sign(g) (-A + sqrt(A^2 + 4 B |g|)) / (2 B),
  // written as 2 g / (A + sqrt(A^2 + 4 B |g|)) so there is no cancellation, it
  // is continuously differentiable through g = 0 (slope 1/A), and it is exactly
  // g / A (Darcy) when B = 0. Returns 0 for a degenerate denominator.
  inline double ErgunVelocity(double g, double A, double B) {
    const double denom = A + std::sqrt(A * A + 4.0 * B * std::abs(g));
    if (!(denom > 0.0)) return 0.0;
    return 2.0 * g / denom;
  }

  // Momentum residual of a face whose velocity u is a DAE unknown, in velocity
  // units:  u - (g - B |u| u) / A,  zero exactly when g = A u + B |u| u, i.e.
  // at u = ErgunVelocity(g, A, B). This implicit form is what the solver uses:
  // it is a C1 polynomial in u with slope (A + 2 B |u|) / A >= 1 and linear in
  // g, so Newton converges robustly, whereas the explicit root has a
  // square-root characteristic (slope 1/A at g = 0 falling to ~1/sqrt(B |g|)
  // within a fraction of a pascal at high Reynolds number) that made Newton
  // fail at step switches of a high-flow PSA cycle. For B = 0 it is exactly
  // the Darcy residual u - g / A.
  // The exception is the boundary face of a pressure or closed boundary
  // condition, which fixes that face's pressure drop: there the condition
  // imposes u - ErgunVelocity(g_bc, A, B) instead, with g_bc taken from the
  // boundary value, because this implicit form's slope in u varies by an
  // order of magnitude when the face velocity jumps at a step switch (see
  // ImposeBoundaryFaceVelocity in boundary_conditions.cpp).
  inline double ErgunMomentumResidual(double u, double g, double A, double B) {
    return u - (g - B * std::abs(u) * u) / A;
  }

  // Writes the file atomically: the JSON goes to a uniquely named temporary
  // file beside it, which then replaces the target in one rename, so a reader
  // never sees a half-written file (the library files are re-saved by tests
  // that run in parallel with tests reading them).
  //
  // Throws std::runtime_error naming the file and the reason when the
  // temporary file cannot be created or written (missing or read-only
  // directory, full disk) or the rename keeps failing; the temporary file is
  // then removed and `file_name` is left as it was. The rename is retried for
  // about 5 s with a back-off from 5 ms to 50 ms: on Windows it fails while
  // another process holds the target open (a test reading the file).
  inline void SaveAdsorbentToJSON(const std::string& file_name, const Adsorbent& adsorbent) {
    static std::atomic<unsigned long> counter{0};
    std::ostringstream tmp_name;
    tmp_name << file_name << ".tmp" << std::hash<std::thread::id>{}(std::this_thread::get_id()) << "_" << counter++;
    const std::string tmp = tmp_name.str();
    auto fail = [&](const std::string& reason) {
      std::error_code ignored;
      std::filesystem::remove(tmp, ignored);
      throw std::runtime_error("SaveAdsorbentToJSON: could not write " + file_name + ": " + reason);
    };
    {
      std::ofstream file(tmp);
      if (!file.is_open()) fail("cannot create the temporary file " + tmp);
      {
        cereal::JSONOutputArchive archive(file);
        archive(adsorbent);
      }  // the archive's destructor completes the JSON before the file closes
      file.close();
      if (file.fail()) fail("writing the temporary file " + tmp + " failed");
    }
    std::error_code ec;
    auto delay = std::chrono::milliseconds(5);
    const auto give_up = std::chrono::steady_clock::now() + std::chrono::seconds(5);
    for (;;) {
      std::filesystem::rename(tmp, file_name, ec);
      if (!ec) return;
      if (std::chrono::steady_clock::now() >= give_up) fail("cannot replace it (" + ec.message() + ")");
      std::this_thread::sleep_for(delay);
      delay = std::min(delay * 2, std::chrono::milliseconds(50));
    }
  }

  inline void LoadAdsorbentFromJSON(const std::string& file_name, Adsorbent& adsorbent) {
    std::ifstream ifs(file_name);
		cereal::JSONInputArchive archive(ifs);
		archive(adsorbent);
  }

} // end adsim namespace

// In the header (not porous_media.cpp) because Adsorbent::serialize branches
// on it: every translation unit that saves or loads a deck must see version 5.
CEREAL_CLASS_VERSION(adsim::Adsorbent, 5);

// Class versions live next to the inline serialize templates: every
// translation unit that saves or loads a deck must see the same numbers.
CEREAL_CLASS_VERSION(adsim::AdsorptionSource, 1);
CEREAL_CLASS_VERSION(adsim::AdsorbentDirectory, 1);
CEREAL_CLASS_VERSION(adsim::StructuredGeometry, 2);
CEREAL_CLASS_VERSION(adsim::MixedAdsorbentLayer, 7);
CEREAL_CLASS_VERSION(adsim::AdsorbentLayers, 1);
