/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <atomic>
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
#include "adsim/objects/fluid.h"
#include "adsim/equations/isotherm.h"
#include "adsim/equations/mass_transfer_model.h"
#include "adsim/objects/reactor_equations.h"
#include "adsim/objects/domain.h"
#include "adsim/equations/equation.h"
#include "adsim/math/expression_param.h"


namespace adsim {

  class MixedAdsorbentLayer;

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
    // denominator instead of recomputing it once per component. Writes into the
    // adsorbent's own equilibrium-loading map (eq_ads_concentration[ads_id]).
    void UpdateIsothermAll(Equations& state, AdsorbentID ads_id, const std::vector<std::string>& components, int start_index, int end_index, const UserVariableRegistry* registry);
    void UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const std::string& component, const Fluid& fluid, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry);
    std::unique_ptr<IIsotherm>& GetIsotherm(const std::string& component);

    std::unordered_map<std::string, std::unique_ptr<IMassTransferModel>> mass_transfer_models_{};
    std::unordered_map<std::string, ExpressionParam> heat_of_adsorption_{}; // Heat of adsorption [J/mol]

    void EvaluateExpressions(const EvalContext& ctx) {
      for (auto& [c, mtm] : mass_transfer_models_) if (mtm) mtm->EvaluateExpressions(ctx);
      for (auto& [t, iso] : isotherms_) if (iso) iso->EvaluateExpressions(ctx);
      for (auto& [c, h] : heat_of_adsorption_) h.Evaluate(ctx);
    }
    void CompileForRun(const UserVariableRegistry* registry) {
      for (auto& [c, mtm] : mass_transfer_models_) if (mtm) mtm->CompileForRun(registry);
    }

  private:
    std::unique_ptr<IIsotherm>& GetIsotherm(IsothermType type, const std::string& component);

    std::unordered_map<IsothermType, std::unique_ptr<IIsotherm>> isotherms_{};
    std::unordered_map<std::string, IsothermType> isotherm_map_{};

    friend class cereal::access;
    template<class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::make_nvp("Mass_Transfer_Models", mass_transfer_models_),
        cereal::make_nvp("Heat_Of_Adsorption", heat_of_adsorption_),
        cereal::make_nvp("Isotherms", isotherms_),
        cereal::make_nvp("Isotherm_Map", isotherm_map_)
      );

      if (version >= 2) {
        // Future expansion
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
  //                 PSA-literature standard, default for new layers).
  // Enumerator values are what the deck stores, so never renumber them.
  enum class MomentumModel {
    KarmenCozeny = 0,
    Ergun = 1
  };
  // Legacy name: the selector used to choose only the permeability.
  using PermeabilityModel = MomentumModel;

  static const char* MomentumModelList[] = {
    "Darcy (Kozeny-Carman)",  // 0
    "Ergun"                   // 1
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

    // Defined out-of-line in porous_media.cpp so it can cascade into
    // adsorption_source_ (isotherms + MTC models) once those are migrated.
    void EvaluateExpressions(const EvalContext& ctx);

  private:
    friend class cereal::access;
    template<class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/) {
      ar(
        cereal::make_nvp("Name", name_),
        cereal::make_nvp("Density", solid_density_),
        cereal::make_nvp("Specific_Heat_Capacity", heat_capacity_),
        cereal::make_nvp("Thermal_Conductivity", thermal_conductivity_),
        cereal::make_nvp("Adsorption_Source", adsorption_source_),
        cereal::make_nvp("Particle_Void_Fraction", ep_)
      );
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
      ar(
        cereal::make_nvp("Adsorbents", adsorbents_),
        cereal::make_nvp("Next_ID", next_id_)
      );
    }
  };

  using LayerID = size_t;

  class StructuredLayer {
  public:
    StructuredLayer() = default;
    ~StructuredLayer() = default;

    std::string name_ = "Layer 1";
    double eb_ = 0.6; // void fraction [m3/m3]: A_void / A total
    double wall_thickness_ = 0.02; // equivelant to particle diameter [m]
    int number_of_channels_ = 100;

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
    ExpressionParam particle_diameter_{1e-3};   // Particle diameter [m]
    // Default for a newly constructed layer is Ergun; archives older than
    // MixedAdsorbentLayer v3 load as KarmenCozeny (see serialize).
    MomentumModel momentum_model_ = MomentumModel::Ergun;
    LayerEquations equations_{};
    std::map<AdsorbentID, double> composition_{};

    void EvaluateExpressions(const EvalContext& ctx) {
      eb_.Evaluate(ctx);
      particle_diameter_.Evaluate(ctx);
      domain_.EvaluateExpressions(ctx);
      equations_.EvaluateExpressions(ctx);
    }

    // Composition management
    void AddAdsorbent(AdsorbentID id, double mass_fraction = 1.0);
    void RemoveAdsorbent(AdsorbentID id);
    bool ValidateComposition() const;

    // Computed properties
    double GetTotalVoidFraction(const AdsorbentDirectory& dir) const {
      const double e = eb_.Get();
      return e + (1 - e) * GetEffectiveParticleVoidFraction(dir);
    }

    // Layer bulk density (1 - eb) * rho_s [kg adsorbent / m3 bed]; equals the
    // sum of GetAdsorbentBulkDensities().
    double GetBulkDensity(const AdsorbentDirectory& dir) const {
      return (1 - eb_.Get()) * GetEffectiveSolidDensity(dir);
    }

    // Viscous permeability K [m2] (A = mu / K) and inertial coefficient
    // beta [1/m] (B = beta * rho) of the layer's momentum model.
    double GetPermeability() const;
    double GetInertialCoefficient() const;

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
      ar(
        cereal::make_nvp("Layer_Equations", equations_),
        cereal::make_nvp("Composition", composition_),
        cereal::make_nvp("ID", id_)
      );
    }
  };

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
  // for Ergun.
  double GetLayerPermeability(double eb, double dp, MomentumModel model);
  // Inertial (Burke-Plummer) coefficient beta [1/m] such that the inertial
  // pressure gradient is beta * rho * |u| * u: 1.75 (1 - eb) / (eb^3 dp) for
  // Ergun, 0 for KarmenCozeny.
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
  inline void SaveAdsorbentToJSON(const std::string& file_name, const Adsorbent& adsorbent) {
    static std::atomic<unsigned long> counter{0};
    std::ostringstream tmp_name;
    tmp_name << file_name << ".tmp" << std::hash<std::thread::id>{}(std::this_thread::get_id()) << "_" << counter++;
    {
      std::ofstream file(tmp_name.str());
      if (!file.is_open()) return;
      cereal::JSONOutputArchive archive(file);
      archive(adsorbent);
    }  // the archive's destructor completes the JSON before the file closes
    std::error_code ec;
    std::filesystem::rename(tmp_name.str(), file_name, ec);
    if (ec) std::filesystem::remove(tmp_name.str(), ec);
  }

  inline void LoadAdsorbentFromJSON(const std::string& file_name, Adsorbent& adsorbent) {
    std::ifstream ifs(file_name);
		cereal::JSONInputArchive archive(ifs);
		archive(adsorbent);
  }

} // end adsim namespace
