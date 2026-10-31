/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/objects/porous_media.h"
#include "adsim/objects/fluid.h"
#include "adsim/objects/reactor.h"
#include "adsim/equations/isotherm.h"
#include "adsim/equations/isotherm_kinds.h"
#include <algorithm>
#include "adsim/math/constants.h"


namespace adsim {

	Adsorbent::Adsorbent(const std::string& name) : name_(name) {}

	void Adsorbent::EvaluateExpressions(const EvalContext& ctx) {
		solid_density_.Evaluate(ctx);
		heat_capacity_.Evaluate(ctx);
		thermal_conductivity_.Evaluate(ctx);
		ep_.Evaluate(ctx);
		adsorption_source_.EvaluateExpressions(ctx);
	}

	void AdsorptionSource::Initialise(const Fluid& fluid) {
		// Purge entries for components the fluid does not contain (e.g. a library
		// adsorbent imported with extra components, or one that missed a
		// Bed::RemoveComponent). Competitive isotherms iterate all of params_ and
		// look each key up in the mole-fraction map, so a stale key would throw.
		const auto& components = fluid.components_;
		auto in_fluid = [&components](const std::string& c) {
			return std::find(components.begin(), components.end(), c) != components.end();
		};
		std::vector<std::string> stale;
		auto collect = [&](const std::string& c) {
			if (!in_fluid(c) && std::find(stale.begin(), stale.end(), c) == stale.end()) {
				stale.push_back(c);
			}
		};
		for (const auto& [c, type] : isotherm_map_) collect(c);
		for (const auto& [c, hoa] : heat_of_adsorption_) collect(c);
		for (const auto& [c, mtm] : mass_transfer_models_) collect(c);
		for (const auto& [type, isotherm] : isotherms_) {
			if (isotherm) for (const auto& c : isotherm->GetComponents()) collect(c);
		}
		for (const auto& c : stale) RemoveComponent(c);

		for (const std::string& component : fluid.components_) {
			auto it_hoa = heat_of_adsorption_.find(component);
			if (it_hoa == heat_of_adsorption_.end()) {
				heat_of_adsorption_[component] = 0.0;
			}

			auto it_mtc = mass_transfer_models_.find(component);
			if (it_mtc == mass_transfer_models_.end()) {
				mass_transfer_models_[component] = std::make_unique<LinearDrivingForce>();
			}

			if (isotherm_map_.find(component) == isotherm_map_.end()) {
				SetIsothermType(component, IsothermType::INERT);
			}

			// Ensure each component only has 1 isotherm (the current one it is pointing to)
			auto isotherm_type = isotherm_map_[component];
			for (auto& [type, isotherm] : isotherms_) {
				if (type != isotherm_type) {
					isotherm->RemoveComponent(component);
				}
			}
		}
	}

	void AdsorptionSource::RemoveComponent(const std::string& component) {
		isotherm_map_.erase(component);
		heat_of_adsorption_.erase(component);
		mass_transfer_models_.erase(component);
		for (auto& [type, isotherm] : isotherms_) {
			if (isotherm) isotherm->RemoveComponent(component);
		}
	}

	void AdsorptionSource::SetIsothermType(const std::string& component, IsothermType isotherm_type) {
		if (isotherms_.find(isotherm_type) == isotherms_.end()) {
			isotherms_[isotherm_type] = CreateIsotherm(isotherm_type, component);
		}

		isotherm_map_[component] = isotherm_type;
	}

	void AdsorptionSource::AssignIsotherm(const std::vector<std::string>& components, IsothermType type) {
		if (components.empty()) return;
		for (const auto& c : components) {
			auto it = isotherm_map_.find(c);
			if (it == isotherm_map_.end() || it->second == type) continue;
			auto old = isotherms_.find(it->second);
			if (old != isotherms_.end() && old->second) old->second->RemoveComponent(c);
		}
		auto& iso = isotherms_[type];
		if (!iso) iso = CreateIsotherm(type, components.front());
		for (const auto& c : components) {
			EnsureIsothermComponent(iso.get(), c);
			isotherm_map_[c] = type;
		}
		for (auto it = isotherms_.begin(); it != isotherms_.end();) {
			const IsothermType t = it->first;
			const bool used = std::any_of(isotherm_map_.begin(), isotherm_map_.end(),
			                              [t](const auto& m) { return m.second == t; });
			it = used ? std::next(it) : isotherms_.erase(it);
		}
	}

	void AdsorptionSource::AddIsotherm(const std::vector<std::string>& components, std::unique_ptr<IIsotherm> isotherm) {
		auto isotherm_type = isotherm->GetIsothermType();
		for (const auto& component : components) {
			isotherm_map_[component] = isotherm_type;
		}
		isotherms_[isotherm_type] = std::move(isotherm);
	}

  void AdsorptionSource::UpdateIsotherm(Equations& state, AdsorbentID ads_id, const std::string& component, int start_index, int end_index, const UserVariableRegistry* registry) {
    auto& isotherm = GetIsotherm(component);
    auto& T_for_isotherm = state.solid_temperature.var;
    // Route the isotherm's qi_sat writes to this adsorbent's own equilibrium-loading
    // map (the inner per-component map), so each adsorbent keeps a separate isotherm.
    IsothermUpdateParams params(state.mole_fraction, state.eq_ads_concentration.at(ads_id), T_for_isotherm, state.pressure.var, component, start_index, end_index, registry);
    isotherm->Update(params);
	}

  void AdsorptionSource::UpdateIsothermAll(Equations& state, AdsorbentID ads_id, const std::vector<std::string>& components, int start_index, int end_index, const UserVariableRegistry* registry) {
    // Group components by isotherm object (preserving component order) so
    // competitive isotherms can share their per-cell denominators.
    std::vector<std::pair<IIsotherm*, std::vector<std::string>>> groups;
    for (const auto& component : components) {
      IIsotherm* isotherm = GetIsotherm(component).get();
      auto it = std::find_if(groups.begin(), groups.end(),
                             [isotherm](const auto& g) { return g.first == isotherm; });
      if (it == groups.end()) {
        groups.emplace_back(isotherm, std::vector<std::string>{ component });
      } else {
        it->second.push_back(component);
      }
    }

    // Route qi_sat writes to this adsorbent's own equilibrium-loading map.
    auto& T_for_isotherm = state.solid_temperature.var;
    for (auto& [isotherm, comps] : groups) {
      IsothermUpdateParams params(state.mole_fraction, state.eq_ads_concentration.at(ads_id), T_for_isotherm, state.pressure.var, comps.front(), start_index, end_index, registry);
      isotherm->UpdateAll(params, comps);
    }
  }

  void AdsorptionSource::UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const std::string& component, const Fluid& fluid, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry, double solid_density) {
    auto& model = mass_transfer_models_.at(component)->mtc_model_;
    model->SetEquilibriumContext(GetIsotherm(component).get(), solid_density);
    model->UpdateMassTransferCoefficient(state, ads_id, fluid, component, layer, registry);
  }

	std::unique_ptr<IIsotherm>& AdsorptionSource::GetIsotherm(const std::string& component) {
		if (isotherm_map_.find(component) == isotherm_map_.end()) {
			SetIsothermType(component, IsothermType::INERT);
		}

		auto isotherm_type = isotherm_map_.at(component);

		return GetIsotherm(isotherm_type, component);
	}

	std::unique_ptr<IIsotherm>& AdsorptionSource::GetIsotherm(IsothermType type, const std::string& component) {
		if (isotherms_.find(type) == isotherms_.end()) {
			isotherms_[type] = CreateIsotherm(type, component);
		}

		return isotherms_.at(type);
	}

	// MixedAdsorbentLayer implementations

	void MixedAdsorbentLayer::AddAdsorbent(AdsorbentID id, double mass_fraction) {
		composition_[id] = mass_fraction;
	}

	void MixedAdsorbentLayer::RemoveAdsorbent(AdsorbentID id) {
		composition_.erase(id);
	}

	bool MixedAdsorbentLayer::ValidateComposition() const {
		if (composition_.empty()) return false;
		double total = 0.0;
		for (const auto& [id, frac] : composition_) {
			if (frac < 0.0 || frac > 1.0) return false;
			total += frac;
		}
		return std::abs(total - 1.0) < 1e-6;
	}

	double MixedAdsorbentLayer::GetPermeability() const {
		return GetLayerPermeability(eb_, particle_diameter_, momentum_model_);
	}

	double MixedAdsorbentLayer::GetInertialCoefficient() const {
		return GetLayerInertialCoefficient(eb_, particle_diameter_, momentum_model_);
	}

	// Mixing rules. composition_ holds MASS fractions w_k, so the extensive
	// (per-volume) properties mix on the solid-volume basis
	//   phi_k = (w_k / rho_k) / sum_j (w_j / rho_j),   sum_k phi_k = 1,
	// giving rho_mix = sum_k phi_k rho_k = 1 / sum_k (w_k / rho_k) when sum_k w_k = 1.
	// Specific (per-mass) properties such as Cp mix on the mass basis directly.
	// Compositions are not normalised: as before, every effective property scales
	// linearly with W = sum_k w_k (W = 1 for a valid composition). For a single
	// adsorbent phi = 1 exactly, so the results are bit-for-bit w * x_k.
	std::map<AdsorbentID, double> MixedAdsorbentLayer::GetSolidVolumeFractions(const AdsorbentDirectory& dir) const {
		std::map<AdsorbentID, double> phi;
		double v_total = 0.0;
		for (const auto& [id, frac] : composition_) {
			const double rho = dir.Get(id).solid_density_;
			const double v = (frac > 0.0 && rho > 0.0) ? frac / rho : 0.0;  // specific solid volume
			phi[id] = v;
			v_total += v;
		}
		for (auto& [id, v] : phi) v = (v_total > 0.0) ? v / v_total : 0.0;
		return phi;
	}

	double MixedAdsorbentLayer::GetTotalMassFraction() const {
		double total = 0.0;
		for (const auto& [id, frac] : composition_) total += frac;
		return total;
	}

	double MixedAdsorbentLayer::GetMixtureSolidDensity(const AdsorbentDirectory& dir) const {
		double density = 0.0;
		for (const auto& [id, phi] : GetSolidVolumeFractions(dir)) {
			density += phi * dir.Get(id).solid_density_;
		}
		return density;
	}

	double MixedAdsorbentLayer::GetEffectiveParticleVoidFraction(const AdsorbentDirectory& dir) const {
		double ep = 0.0;
		for (const auto& [id, phi] : GetSolidVolumeFractions(dir)) {
			ep += phi * dir.Get(id).ep_;
		}
		return GetTotalMassFraction() * ep;
	}

	double MixedAdsorbentLayer::GetEffectiveSolidDensity(const AdsorbentDirectory& dir) const {
		return GetTotalMassFraction() * GetMixtureSolidDensity(dir);
	}

	std::map<AdsorbentID, double> MixedAdsorbentLayer::GetAdsorbentBulkDensities(const AdsorbentDirectory& dir) const {
		// rho_b,k = (1 - eb) * rho_mix * w_k: the mass of adsorbent k per bed volume.
		// These sum to GetBulkDensity() = (1 - eb) * W * rho_mix.
		const double rho_mix = GetMixtureSolidDensity(dir);
		const double solid_fraction = 1.0 - eb_.Get();
		std::map<AdsorbentID, double> rho_b;
		for (const auto& [id, frac] : composition_) {
			rho_b[id] = solid_fraction * (rho_mix * frac);
		}
		return rho_b;
	}

	double MixedAdsorbentLayer::GetEffectiveHeatCapacity(const AdsorbentDirectory& dir) const {
		double cp = 0.0;
		for (const auto& [id, frac] : composition_) {
			cp += frac * dir.Get(id).heat_capacity_;
		}
		return cp;
	}

	double MixedAdsorbentLayer::GetEffectiveThermalConductivity(const AdsorbentDirectory& dir) const {
		if (composition_.empty()) return 0.5;
		double k = 0.0;
		for (const auto& [id, phi] : GetSolidVolumeFractions(dir)) {
			k += phi * dir.Get(id).thermal_conductivity_;
		}
		return GetTotalMassFraction() * k;
	}

	void MixedAdsorbentLayer::Initialise(const AdsorbentDirectory& dir, const Fluid& fluid,
	                                      int start_index, int end_index) {
		domain_.start_index = start_index;
		domain_.end_index = end_index;
		for (auto& [ads_id, frac] : composition_) {
			auto& adsorbent = const_cast<AdsorbentDirectory&>(dir).Get(ads_id);
			adsorbent.adsorption_source_.Initialise(fluid);
		}
	}

	// AdsorbentLayers implementations

	LayerID AdsorbentLayers::Add(const std::string& layer_name) {
		++next_layer_id_;
		layers_.emplace_back(layer_name);
		layers_.back().id_ = next_layer_id_;
		return next_layer_id_;
	}

	void AdsorbentLayers::Remove(LayerID id) {
		auto it = std::remove_if(layers_.begin(), layers_.end(),
			[id](const MixedAdsorbentLayer& layer) {
			return layer.id_ == id;
		});

		if (it != layers_.end()) {
			layers_.erase(it, layers_.end());
		}
	}

	bool AdsorbentLayers::Swap(LayerID id1, LayerID id2) {
		auto it1 = std::find_if(layers_.begin(), layers_.end(),
			[id1](const MixedAdsorbentLayer& layer) {
			return layer.id_ == id1;
		});

		auto it2 = std::find_if(layers_.begin(), layers_.end(),
			[id2](const MixedAdsorbentLayer& layer) {
			return layer.id_ == id2;
		});

		if (it1 != layers_.end() && it2 != layers_.end()) {
			std::swap(*it1, *it2);
			return true;
		}

		return false;
	}

	MixedAdsorbentLayer& AdsorbentLayers::Get(LayerID id) {
		auto it = std::find_if(layers_.begin(), layers_.end(),
			[id](const MixedAdsorbentLayer& layer) {
			return layer.id_ == id;
		});

		if (it != layers_.end()) {
			return *it;
		}

		throw std::runtime_error("Layer with given ID does not exist");
	}

	const MixedAdsorbentLayer& AdsorbentLayers::Get(LayerID id) const {
		auto it = std::find_if(layers_.begin(), layers_.end(),
			[id](const MixedAdsorbentLayer& layer) {
			return layer.id_ == id;
		});

		if (it != layers_.end()) {
			return *it;
		}

		throw std::runtime_error("Layer with given ID does not exist");
	}

  double GetLayerPermeability(double eb, double dp, MomentumModel model) {
    if (model == MomentumModel::KarmenCozeny) {
      return eb * eb * eb * dp * dp / (180 * (1 - eb) * (1 - eb));
    }
    if (model == MomentumModel::Ergun) {
      return eb * eb * eb * dp * dp / (150 * (1 - eb) * (1 - eb));
    }

    return 0.0;
  }

  double GetLayerInertialCoefficient(double eb, double dp, MomentumModel model) {
    if (model == MomentumModel::Ergun) {
      return 1.75 * (1 - eb) / (eb * eb * eb * dp);
    }
    return 0.0;
  }

} // end adsim namespace

CEREAL_CLASS_VERSION(adsim::Adsorbent, 4);
CEREAL_CLASS_VERSION(adsim::AdsorbentDirectory, 1);
CEREAL_CLASS_VERSION(adsim::MixedAdsorbentLayer, 3);
CEREAL_CLASS_VERSION(adsim::AdsorbentLayers, 1);
