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
#include <iomanip>
#include "adsim/math/constants.h"


namespace adsim {

	Adsorbent::Adsorbent(const std::string& name) : name_(name) {}

	void Adsorbent::EvaluateExpressions(const EvalContext& ctx) {
		solid_density_.Evaluate(ctx);
		heat_capacity_.Evaluate(ctx);
		thermal_conductivity_.Evaluate(ctx);
		ep_.Evaluate(ctx);
		adsorption_source_.EvaluateExpressions(ctx);
		reaction_source_.EvaluateExpressions(ctx);
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
			// A dual-rate model's slow site takes the fast site's coefficient kind.
			if (auto* dual = dynamic_cast<DualRateLinearDrivingForce*>(mass_transfer_models_[component].get()))
				dual->SyncSlowSite();

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

  namespace {
    // The adsorbent's actual loadings, for the isotherm kinds that read
    // another component's loading (IsothermUpdateParams::q_map).
    const ComponentEquation* AdsorbedLoadings(const Equations& state, AdsorbentID ads_id) {
      auto it = state.ads_concentration.find(ads_id);
      return it == state.ads_concentration.end() ? nullptr : &it->second;
    }
  }  // namespace

  void AdsorptionSource::UpdateIsotherm(Equations& state, AdsorbentID ads_id, const std::string& component, int start_index, int end_index, const UserVariableRegistry* registry) {
    auto& isotherm = GetIsotherm(component);
    auto& T_for_isotherm = state.solid_temperature.var;
    // Route the isotherm's qi_sat writes to this adsorbent's own equilibrium-loading
    // map (the inner per-component map), so each adsorbent keeps a separate isotherm.
    IsothermUpdateParams params(state.mole_fraction, state.eq_ads_concentration.at(ads_id), T_for_isotherm, state.pressure.var, component, start_index, end_index, registry, AdsorbedLoadings(state, ads_id));
    isotherm->Update(params);
	}

  void AdsorptionSource::UpdateIsothermAll(Equations& state, AdsorbentID ads_id, const std::vector<std::string>& components, int start_index, int end_index, const UserVariableRegistry* registry) {
    if (mixing_rule_ == MixingRule::IAST) {
      UpdateIsothermIast(state, ads_id, components, start_index, end_index);
      return;
    }

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
    const ComponentEquation* q_map = AdsorbedLoadings(state, ads_id);
    for (auto& [isotherm, comps] : groups) {
      IsothermUpdateParams params(state.mole_fraction, state.eq_ads_concentration.at(ads_id), T_for_isotherm, state.pressure.var, comps.front(), start_index, end_index, registry, q_map);
      isotherm->UpdateAll(params, comps);
    }
  }

  void AdsorptionSource::UpdateIsothermIast(Equations& state, AdsorbentID ads_id, const std::vector<std::string>& components, int start_index, int end_index) {
    const int n = static_cast<int>(components.size());
    auto& T = state.solid_temperature.var;
    auto& P = state.pressure.var;
    auto& q_map = state.eq_ads_concentration.at(ads_id);

    // Resolve each component's pure isotherm once per update; the cell loop
    // then only applies the temperature factors.
    iast_refs_.resize(n);
    iast_pure_.resize(n);
    iast_y_.resize(n);
    iast_q_.resize(n);
    iast_y_of_.resize(n);
    iast_q_of_.resize(n);
    for (int i = 0; i < n; ++i) {
      const auto& comp = components[i];
      iast_refs_[i] = ResolvePureIsotherm(GetIsotherm(comp).get(), comp);
      iast_y_of_[i] = &state.mole_fraction.at(comp).var;
      iast_q_of_[i] = &q_map.at(comp).var;
    }
    if (iast_psi_.size() != P.size()) iast_psi_.assign(P.size(), 0.0);

    for (int c = start_index; c <= end_index; ++c) {
      for (int i = 0; i < n; ++i) {
        iast_pure_[i] = PreparePureIsotherm(iast_refs_[i], T[c]);
        iast_y_[i] = (*iast_y_of_[i])[c];
      }
      const IastResult r = SolveIast(iast_pure_.data(), n, P[c], iast_y_.data(), iast_q_.data(),
                                     iast_scratch_, iast_psi_[c]);
      iast_psi_[c] = r.converged ? r.psi : 0.0;
      for (int i = 0; i < n; ++i) (*iast_q_of_[i])[c] = iast_q_[i];
    }
  }

  void AdsorptionSource::UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const std::string& component, const Fluid& fluid, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry, double solid_density) {
    auto& model = mass_transfer_models_.at(component)->mtc_model_;
    model->SetEquilibriumContext(GetIsotherm(component).get(), solid_density, this);
    model->UpdateMassTransferCoefficient(state, ads_id, fluid, component, layer, registry);
  }

  void AdsorptionSource::UpdateSlowMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const std::string& component, const Fluid& fluid, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry, double solid_density) {
    IMassTransferCoefficientModel* model = mass_transfer_models_.at(component)->SlowCoefficientModel();
    auto slow_it = state.slow_mass_transfer_coefficient.find(ads_id);
    if (!model || slow_it == state.slow_mass_transfer_coefficient.end()) return;
    auto& k_slow = slow_it->second.at(component).var;
    // The coefficient models write the state's (fast-site) coefficient field:
    // swap the slow-site field in for the call (O(1), no copy) and back.
    auto& k_fast = state.mass_transfer_coefficient.at(ads_id).at(component).var;
    std::swap(k_fast, k_slow);
    model->SetEquilibriumContext(GetIsotherm(component).get(), solid_density, this);
    model->UpdateMassTransferCoefficient(state, ads_id, fluid, component, layer, registry);
    std::swap(k_fast, k_slow);
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

	std::map<AdsorbentID, double> MixedAdsorbentLayer::GetMassFractions() const {
		std::map<AdsorbentID, double> w;
		for (const auto& [id, frac] : composition_) w[id] = frac.Get();
		return w;
	}

	bool MixedAdsorbentLayer::HasCompositionExpression() const {
		for (const auto& [id, frac] : composition_)
			if (frac.HasExpression()) return true;
		return false;
	}

	bool MixedAdsorbentLayer::ValidateComposition() const {
		if (composition_.empty()) return false;
		double total = 0.0;
		for (const auto& [id, w] : composition_) {
			const double frac = w.Get();
			if (!w.LastError().empty()) return false;
			if (!(frac >= 0.0)) return false;  // also rejects NaN
			total += frac;
		}
		return std::abs(total - 1.0) <= kCompositionSumTolerance;
	}

	std::vector<LayerCompositionIssue> CheckLayerComposition(const MixedAdsorbentLayer& layer,
	                                                         const AdsorbentDirectory& dir) {
		using Severity = LayerCompositionIssue::Severity;
		std::vector<LayerCompositionIssue> issues;
		const std::string ln = "\"" + layer.name_ + "\"";
		auto number = [](double v) {
			std::ostringstream ss;
			ss << std::setprecision(8) << v;
			return ss.str();
		};
		if (layer.composition_.empty()) {
			issues.push_back({Severity::Error, "Error: Layer " + ln +
				" has no adsorbent selected. Each layer must have at least one adsorbent."});
			return issues;
		}
		double total = 0.0;
		bool skip_sum = false;
		for (const auto& [id, w] : layer.composition_) {
			const double frac = w.Get();
			total += frac;
			// A fraction whose expression did not evaluate (an unknown
			// variable, a syntax error) keeps a stale value; report it.
			if (!w.LastError().empty()) {
				issues.push_back({Severity::Error, "Error: Layer " + ln + ": the mass fraction '" + w.GetExpr() +
					"' of adsorbent id " + std::to_string(id) + " cannot be evaluated: " + w.LastError()});
				skip_sum = true;
				continue;
			}
			// A missing id would make every directory lookup during
			// initialisation throw.
			if (!dir.Contains(id)) {
				issues.push_back({Severity::Error, "Error: Layer " + ln + " references adsorbent id " + std::to_string(id) +
					", which does not exist in the adsorbent library. Remove it from the layer composition or select another adsorbent."});
				continue;
			}
			const std::string an = "\"" + dir.Get(id).name_ + "\"";
			// Isotherm kinds that read the adsorbed water loading
			// (IsWaterCoupledIsothermKind): no IAST form, and without a water
			// component they are only their dry isotherm.
			{
				const auto& source = dir.Get(id).adsorption_source_;
				for (const auto& [component, type] : source.IsothermMap()) {
					if (!IsWaterCoupledIsothermKind(type)) continue;
					const std::string kind = std::string("'") + isothermModelNames[static_cast<int>(type)] + "' (" +
						IsothermKindString(type) + ")";
					if (source.mixing_rule_ == MixingRule::IAST)
						issues.push_back({Severity::Error, "Error: adsorbent " + an + " uses the IAST mixing rule, but " +
							component + " has the isotherm " + kind + ", whose loading depends on the adsorbed water loading"
							" and which IAST cannot mix. Set the adsorbent's mixing rule to none."});
					if (source.IsothermMap().find(kCoadsorbedWaterComponent) == source.IsothermMap().end())
						issues.push_back({Severity::Warning, "Warning: adsorbent " + an + " gives " + component +
							" the isotherm " + kind + ", which depends on the loading of the component named " +
							kCoadsorbedWaterComponent + ", but the deck has no such component: the isotherm is"
							" its dry form. Add the component " + kCoadsorbedWaterComponent + " with its own isotherm."});
				}
			}
			// Dual-rate LDF: the fast-site fraction must lie in (0, 1].
			for (const auto& [component, mtm] : dir.Get(id).adsorption_source_.mass_transfer_models_) {
				if (!mtm || mtm->GetModelType() != MassTransferModelType::DUAL_RATE_LDF) continue;
				const double phi = mtm->FastSiteFraction();
				if (!(phi > 0.0 && phi <= 1.0))
					issues.push_back({Severity::Error, "Error: adsorbent " + an + " gives " + component +
						" a dual-rate LDF model with the fast-site fraction (split) " + number(phi) +
						"; it must lie in (0, 1]."});
			}
			if (!(frac >= 0.0)) {
				skip_sum = true;
				issues.push_back({Severity::Error, "Error: Layer " + ln + " gives adsorbent " + an + " the mass fraction " +
					number(frac) + ". Mass fractions must be between 0 and 1."});
			} else if (frac == 0.0) {
				issues.push_back({Severity::Warning, "Warning: adsorbent " + an + " has mass fraction 0 in layer " + ln +
					", so it takes no part in the run. Remove it from the layer's composition."});
			}
		}
		if (!skip_sum && !(std::abs(total - 1.0) <= kCompositionSumTolerance)) {
			// Formulas are shown with their values, so the user sees which
			// variable value broke the sum.
			std::string terms;
			if (layer.HasCompositionExpression()) {
				for (const auto& [id, w] : layer.composition_) {
					terms += terms.empty() ? " (" : ", ";
					terms += w.HasExpression() ? w.GetExpr() + " = " + number(w.Get()) : number(w.Get());
				}
				terms += ")";
			}
			issues.push_back({Severity::Error, "Error: Layer " + ln + ": mass fractions" + terms + " sum to " + number(total) +
				"; they must sum to 1. Change the fractions so they add up to 1 (they are not rescaled automatically)."
				" To vary one solid's share in a study, write the fractions as expressions of a deck variable: for two"
				" solids, set them to x and 1 - x with a deck variable x, and vary x."});
		}
		return issues;
	}

	double MixedAdsorbentLayer::GetPermeability() const {
		if (!IsStructured()) return GetLayerPermeability(eb_, particle_diameter_, momentum_model_);
		// Laminar channel flow: u_channel = D_h^2 / (2 fRe mu) (-dP/dz) with the
		// Fanning fRe, and the superficial velocity is OFA u_channel.
		const double dh = GetHydraulicDiameter();
		return GetOpenFraction() * dh * dh / (2.0 * GetFrictionFactorReynolds());
	}

	double MixedAdsorbentLayer::GetInertialCoefficient() const {
		if (!IsStructured()) return GetLayerInertialCoefficient(eb_, particle_diameter_, momentum_model_);
		return 0.0;  // laminar channels: the linear Darcy path
	}

	namespace {
		constexpr double kPi = 3.14159265358979323846;

		// Monolith pitch p and open channel width a = p - w - 2 t_wc [m].
		void MonolithWidths(const StructuredGeometry& g, double& p, double& a) {
			const double n = g.cell_density_.Get();
			p = n > 0.0 ? 1.0 / std::sqrt(n) : 0.0;
			a = p - g.wall_thickness_.Get() - 2.0 * g.washcoat_thickness_.Get();
		}
	}  // namespace

	double MixedAdsorbentLayer::GetOpenFraction() const {
		const StructuredGeometry& g = structured_;
		switch (kind_) {
			case LayerKind::Packed: return eb_.Get();
			case LayerKind::Monolith: {
				double p, a;
				MonolithWidths(g, p, a);
				if (!(p > 0.0)) return 0.0;
				const double square = a * a / (p * p);
				return g.channel_shape_ == ChannelShape::Circular ? 0.25 * kPi * square : square;
			}
			case LayerKind::Laminate: {
				const double period = g.gap_.Get() + g.sheet_thickness_.Get();
				return period > 0.0 ? (g.gap_.Get() - 2.0 * g.washcoat_thickness_.Get()) / period : 0.0;
			}
		}
		return 0.0;
	}

	double MixedAdsorbentLayer::GetHydraulicDiameter() const {
		const StructuredGeometry& g = structured_;
		switch (kind_) {
			case LayerKind::Packed: {
				const double e = eb_.Get();
				return (2.0 / 3.0) * e * particle_diameter_.Get() / (1.0 - e);
			}
			case LayerKind::Monolith: {
				double p, a;
				MonolithWidths(g, p, a);
				return a;  // square: 4 a^2 / (4 a); circular: the diameter
			}
			case LayerKind::Laminate:
				return 2.0 * (g.gap_.Get() - 2.0 * g.washcoat_thickness_.Get());
		}
		return 0.0;
	}

	double MixedAdsorbentLayer::GetSpecificSurfaceArea() const {
		const StructuredGeometry& g = structured_;
		switch (kind_) {
			case LayerKind::Packed: {
				const double dp = particle_diameter_.Get();
				return dp > 0.0 ? 6.0 * (1.0 - eb_.Get()) / dp : 0.0;
			}
			case LayerKind::Monolith: {
				double p, a;
				MonolithWidths(g, p, a);
				if (!(p > 0.0)) return 0.0;
				const double perimeter = g.channel_shape_ == ChannelShape::Circular ? kPi * a : 4.0 * a;
				return perimeter / (p * p);
			}
			case LayerKind::Laminate: {
				const double period = g.gap_.Get() + g.sheet_thickness_.Get();
				return period > 0.0 ? 2.0 / period : 0.0;
			}
		}
		return 0.0;
	}

	double MixedAdsorbentLayer::GetFrictionFactorReynolds() const {
		switch (kind_) {
			case LayerKind::Packed: return 0.0;
			case LayerKind::Monolith: return structured_.channel_shape_ == ChannelShape::Circular ? 16.0 : 14.23;
			case LayerKind::Laminate: return 24.0;
		}
		return 0.0;
	}

	double MixedAdsorbentLayer::GetDiffusionLength() const {
		const StructuredGeometry& g = structured_;
		const double t_wc = g.washcoat_thickness_.Get();
		switch (kind_) {
			case LayerKind::Packed: return 0.5 * particle_diameter_.Get();
			case LayerKind::Monolith: return t_wc > 0.0 ? t_wc : 0.5 * g.wall_thickness_.Get();
			case LayerKind::Laminate: return t_wc > 0.0 ? t_wc : 0.5 * g.sheet_thickness_.Get();
		}
		return 0.0;
	}

	double MixedAdsorbentLayer::GetWashcoatFraction() const {
		const StructuredGeometry& g = structured_;
		const double t_wc = g.washcoat_thickness_.Get();
		if (!(t_wc > 0.0)) return 0.0;
		switch (kind_) {
			case LayerKind::Packed: return 0.0;
			case LayerKind::Monolith: {
				double p, a;
				MonolithWidths(g, p, a);
				const double b = p - g.wall_thickness_.Get();  // substrate channel width
				const double k = g.channel_shape_ == ChannelShape::Circular ? 0.25 * kPi : 1.0;
				const double solid = p * p - k * a * a;
				return solid > 0.0 ? k * (b * b - a * a) / solid : 0.0;
			}
			case LayerKind::Laminate: {
				const double solid = g.sheet_thickness_.Get() + 2.0 * t_wc;
				return solid > 0.0 ? 2.0 * t_wc / solid : 0.0;
			}
		}
		return 0.0;
	}

	double MixedAdsorbentLayer::GetChannelNusselt() const {
		switch (kind_) {
			case LayerKind::Packed: return 0.0;
			case LayerKind::Monolith: return structured_.channel_shape_ == ChannelShape::Circular ? 3.657 : 2.976;
			case LayerKind::Laminate: return 7.541;
		}
		return 0.0;
	}

	double MixedAdsorbentLayer::GetTaylorArisKappa() const {
		switch (kind_) {
			case LayerKind::Packed: return 0.0;
			case LayerKind::Monolith: return 192.0;
			case LayerKind::Laminate: return 840.0;
		}
		return 0.0;
	}

	MixedAdsorbentLayer::DiffusionGeometry MixedAdsorbentLayer::GetDiffusionGeometry() const {
		DiffusionGeometry g;
		if (!IsStructured()) {
			g.c = 60.0;
			g.length = particle_diameter_.Get();
			g.volume_to_surface = particle_diameter_.Get() / 6.0;
		} else {
			g.c = 3.0;
			g.length = GetDiffusionLength();
			g.volume_to_surface = g.length;
		}
		return g;
	}

	LayerGeometry LayerGeometrySummary(const MixedAdsorbentLayer& layer, double cross_section_area) {
		LayerGeometry s;
		s.kind = layer.kind_;
		s.open_fraction = layer.GetOpenFraction();
		s.hydraulic_diameter = layer.GetHydraulicDiameter();
		s.specific_surface_area = layer.GetSpecificSurfaceArea();
		s.diffusion_length = layer.GetDiffusionLength();
		s.washcoat_fraction = layer.GetWashcoatFraction();
		s.friction_factor_reynolds = layer.GetFrictionFactorReynolds();
		if (layer.kind_ == LayerKind::Monolith)
			s.channel_count = layer.structured_.cell_density_.Get() * cross_section_area;
		return s;
	}

	std::vector<std::string> CheckLayerGeometry(const MixedAdsorbentLayer& layer) {
		std::vector<std::string> errors;
		if (!layer.IsStructured()) return errors;
		const std::string ln = "Error: Layer \"" + layer.name_ + "\" (" + LayerKindList[static_cast<int>(layer.kind_)] + "): ";
		auto number = [](double v) {
			std::ostringstream ss;
			ss << std::setprecision(6) << v;
			return ss.str();
		};
		const StructuredGeometry& g = layer.structured_;
		auto positive = [&](const char* what, double v) {
			if (!(v > 0.0)) errors.push_back(ln + what + " must be greater than 0, not " + number(v) + ".");
		};
		if (!(g.washcoat_thickness_.Get() >= 0.0))
			errors.push_back(ln + "the washcoat thickness must not be negative, not " + number(g.washcoat_thickness_.Get()) + ".");
		if (layer.kind_ == LayerKind::Monolith) {
			positive("the cell density", g.cell_density_.Get());
			positive("the wall thickness", g.wall_thickness_.Get());
			double p, a;
			MonolithWidths(g, p, a);
			if (p > 0.0 && !(a > 0.0))
				errors.push_back(ln + "the walls (" + number(g.wall_thickness_.Get()) + " m) and twice the washcoat (" +
				                 number(2.0 * g.washcoat_thickness_.Get()) + " m) fill the channel pitch (" + number(p) +
				                 " m), leaving no open channel.");
		} else {
			positive("the sheet thickness", g.sheet_thickness_.Get());
			positive("the gap", g.gap_.Get());
			if (g.gap_.Get() > 0.0 && !(g.gap_.Get() - 2.0 * g.washcoat_thickness_.Get() > 0.0))
				errors.push_back(ln + "twice the washcoat (" + number(2.0 * g.washcoat_thickness_.Get()) +
				                 " m) fills the gap (" + number(g.gap_.Get()) + " m), leaving no open channel.");
		}
		if (errors.empty()) {
			const double ofa = layer.GetOpenFraction();
			if (!(ofa > 0.0 && ofa < 1.0))
				errors.push_back(ln + "the open fraction is " + number(ofa) + "; it must lie between 0 and 1.");
			if (!(layer.GetHydraulicDiameter() > 0.0))
				errors.push_back(ln + "the hydraulic diameter is " + number(layer.GetHydraulicDiameter()) + " m; it must be positive.");
			if (!(layer.GetDiffusionLength() > 0.0))
				errors.push_back(ln + "the diffusion length is " + number(layer.GetDiffusionLength()) + " m; it must be positive.");
		}
		return errors;
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
		const double solid_fraction = 1.0 - GetVoidFraction();
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
    if (model == MomentumModel::Ergun || model == MomentumModel::BlakeKozeny) {
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
