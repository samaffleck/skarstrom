/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "adsim/objects/deck_view.h"

#include <algorithm>
#include <map>
#include <string>
#include <string_view>
#include <vector>

#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/equations/isotherm_kinds.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/manifests.h"
#include "adsim/objects/model_kinds.h"
#include "adsim/objects/user_data.h"
#include "adsim/parameters/parameter_registry.h"
#include "adsim/version.h"
#include "view_json.h"

namespace adsim {

  namespace {

    using namespace view_json;

    std::string KindName(const std::vector<std::string>& kinds, int i) {
      return i >= 0 && i < static_cast<int>(kinds.size()) ? kinds[static_cast<size_t>(i)] : std::string("unknown");
    }

    // The registry's entries, each handed to the first node that claims it.
    class ParameterPool {
    public:
      explicit ParameterPool(const ParameterRegistry& registry) {
        for (const auto& p : registry.All()) entries_.push_back(&p);
        used_.assign(entries_.size(), false);
      }

      // The unclaimed entries whose path starts with `prefix`, in registry order.
      std::vector<const Parameter*> Claim(std::string_view prefix) {
        std::vector<const Parameter*> out;
        for (size_t i = 0; i < entries_.size(); ++i) {
          if (used_[i]) continue;
          const std::string& path = entries_[i]->Path();
          if (path.compare(0, prefix.size(), prefix) == 0) {
            used_[i] = true;
            out.push_back(entries_[i]);
          }
        }
        return out;
      }

      const Parameter* Find(std::string_view path) const {
        for (const auto* p : entries_)
          if (p->Path() == path) return p;
        return nullptr;
      }

    private:
      std::vector<const Parameter*> entries_;
      std::vector<bool> used_;
    };

    // One registry entry: path, display data and value.
    Value ParameterValue(const Parameter& p, std::string_view prefix, Allocator& a) {
      const ParameterInfo& info = p.Info();
      Value v(rapidjson::kObjectType);
      Set(v, "path", info.path, a);
      Set(v, "alias", info.alias.empty() ? info.path : info.alias, a);
      Set(v, "key", std::string_view(info.path).substr(std::min(prefix.size(), info.path.size())), a);
      Set(v, "label", info.name, a);
      Set(v, "category", info.category, a);
      Set(v, "subcategory", info.subcategory, a);
      Set(v, "kind", ParameterKindName(info.kind), a);
      switch (info.kind) {
        case ParameterKind::String:
          Set(v, "value", p.GetText(), a);
          break;
        case ParameterKind::Boolean:
          Set(v, "value", p.GetValue() != 0.0, a);
          break;
        default:
          Set(v, "value", Num(p.GetValue()), a);
          break;
      }
      Set(v, "text", p.GetText(), a);
      const bool formula = info.kind == ParameterKind::Expression && p.ExpressionPointer() &&
                           p.ExpressionPointer()->HasExpression();
      Set(v, "formula", formula, a);
      Set(v, "unit", info.unit, a);
      Set(v, "unit_checked", info.unit_checked, a);
      Set(v, "tunable", info.tunable, a);
      Set(v, "lower", Num(info.lower), a);
      Set(v, "upper", Num(info.upper), a);
      if (!info.choices.empty()) {
        Value choices(rapidjson::kArrayType);
        for (size_t i = 0; i < info.choices.size(); ++i) {
          Value c(rapidjson::kObjectType);
          Set(c, "name", info.choices[i], a);
          SetInt(c, "index", static_cast<long long>(i), a);
          choices.PushBack(c, a);
        }
        Set(v, "choices", std::move(choices), a);
        const double x = p.GetValue();
        const long long i = std::isfinite(x) ? static_cast<long long>(std::llround(x)) : -1;
        if (i >= 0 && i < static_cast<long long>(info.choices.size()))
          Set(v, "choice", info.choices[static_cast<size_t>(i)], a);
        else
          SetNull(v, "choice", a);
      } else {
        Set(v, "choices", Value(rapidjson::kArrayType), a);
        SetNull(v, "choice", a);
      }
      return v;
    }

    class ViewBuilder {
    public:
      ViewBuilder(UserData& data, const ParameterRegistry& registry, Allocator& a)
          : d_(data), pool_(registry), a_(a) {}

      // {"values": [...]} for the entries under `prefix`.
      Value Values(std::string_view prefix) {
        Value list(rapidjson::kArrayType);
        for (const auto* p : pool_.Claim(prefix)) list.PushBack(ParameterValue(*p, prefix, a_), a_);
        return list;
      }
      Value Node(std::string_view prefix) {
        Value n(rapidjson::kObjectType);
        Set(n, "values", Values(prefix), a_);
        return n;
      }

      // The choice name of an enumerated registry entry, or "" when absent.
      std::string ChoiceOf(std::string_view path) const {
        const Parameter* p = pool_.Find(path);
        if (!p || p->Info().choices.empty()) return {};
        const long long i = std::llround(p->GetValue());
        return i >= 0 && i < static_cast<long long>(p->Info().choices.size()) ? p->Info().choices[static_cast<size_t>(i)]
                                                                               : std::string{};
      }

      Value Build() {
        namespace b = build;
        Value root(rapidjson::kObjectType);
        SetInt(root, "view_version", kDeckViewVersion, a_);
        Set(root, "engine_version", kVersionString, a_);

        // The static sections first: their paths are fixed.
        for (const char* section : {"info", "settings", "css", "energy_balance"})
          Set(root, section, Node(std::string(section) + "."), a_);

        Set(root, "components", Components(), a_);
        Set(root, "fluid", Node("fluid."), a_);
        Set(root, "initial_conditions", Node("initial."), a_);
        Set(root, "adsorbents", Adsorbents(), a_);
        Set(root, "layers", Layers(), a_);
        Set(root, "bed", Bed(), a_);

        Value cycle(rapidjson::kObjectType);
        Set(cycle, "step_time_mode", ChoiceOf("cycle.step_time_mode"), a_);
        Set(cycle, "values", Values("cycle."), a_);
        Set(root, "cycle", std::move(cycle), a_);
        Set(root, "steps", Steps(), a_);
        Set(root, "connections", Connections(), a_);
        Set(root, "kpis", Kpis(), a_);
        Set(root, "reaction_kpi", ReactionKpi(), a_);
        Set(root, "variables", Variables(), a_);
        Set(root, "cost", Node("cost."), a_);
        Set(root, "sensitivity", Node("sensitivity."), a_);
        Set(root, "experimental_data", ExperimentalDatasets(), a_);

        // Whatever no node claimed (a registry family the view does not know).
        Set(root, "other", Node(""), a_);
        return root;
      }

    private:
      UserData& d_;
      ParameterPool pool_;
      Allocator& a_;

      Value Components() {
        Value list(rapidjson::kArrayType);
        for (const auto& c : d_.fluid.components_) {
          const ComponentInfo info = DescribeComponent(c);
          Value v(rapidjson::kObjectType);
          Set(v, "name", c, a_);
          Set(v, "formula", info.formula, a_);
          const auto it = d_.fluid.molar_mass_.find(c);
          Set(v, "molar_mass", Num(it != d_.fluid.molar_mass_.end() && it->second > 0 ? it->second : info.molar_mass), a_);
          list.PushBack(v, a_);
        }
        return list;
      }

      // Whether the adsorbent sits in a packed layer and in a structured one
      // (a layer whose composition names it). In no layer: packed.
      // `shape` is the particle shape of the first packed layer that holds it
      // (Sphere when none does).
      // `film` is the film correlation of that layer.
      struct Shapes {
        bool packed = false, structured = false;
        DiffusionShape shape = DiffusionShape::Sphere;
        FilmCorrelation film = FilmCorrelation::WakaoFunazkri;
      };
      Shapes LayerShapes(AdsorbentID id) {
        Shapes s;
        for (const auto& layer : d_.bed.layers_.layers_) {
          if (!layer.composition_.count(id)) continue;
          if (layer.IsStructured()) { s.structured = true; continue; }
          if (!s.packed) {
            s.shape = LayerDiffusionShape(layer);
            s.film = layer.film_correlation_;
          }
          s.packed = true;
        }
        if (!s.packed && !s.structured) s.packed = true;
        return s;
      }

      // The generic coefficient model of one site (prefix "" or "slow_"):
      // <pre>resistances (switch states by tag), <pre>capacity_factor_mode,
      // <pre>preset (the derived preset tag, or custom), <pre>equation (LaTeX;
      // the packed form, sphere or cylinder as the first packed layer that
      // holds the adsorbent, when it is in a packed layer, else the slab
      // form) and <pre>equation_slab (the slab form when it is in layers of
      // both kinds, else null). All null without a generic model on the site.
      void Coefficient(Value& cv, const std::string& pre, IMassTransferCoefficientModel* m, AdsorbentID id) {
        auto* g = dynamic_cast<MassTransferCoefficientSeriesResistance*>(m);
        if (!g) {
          for (const char* f : {"resistances", "capacity_factor_mode", "preset", "equation", "equation_slab"})
            SetNull(cv, pre + f, a_);
          return;
        }
        Value r(rapidjson::kObjectType);
        for (const auto& t : ResistanceTags()) Set(r, t, g->Enabled(t), a_);
        Set(cv, pre + "resistances", std::move(r), a_);
        Set(cv, pre + "capacity_factor_mode", std::string(CapacityFactorModeList[static_cast<int>(g->CapacityMode())]), a_);
        Set(cv, pre + "preset", g->PresetTag(), a_);
        const auto [packed, structured, shape, film] = LayerShapes(id);
        Set(cv, pre + "equation",
            packed ? MassTransferEquationLatex(*g, shape, &film) : MassTransferEquationLatex(*g, DiffusionShape::Slab), a_);
        if (packed && structured) Set(cv, pre + "equation_slab", MassTransferEquationLatex(*g, DiffusionShape::Slab), a_);
        else SetNull(cv, pre + "equation_slab", a_);
      }

      Value Adsorbents() {
        namespace b = build;
        Value list(rapidjson::kArrayType);
        std::vector<AdsorbentID> ids;
        for (const auto& [id, ads] : d_.adsorbent_directory.GetAll()) ids.push_back(id);
        std::sort(ids.begin(), ids.end());
        for (AdsorbentID id : ids) {
          auto& ads = d_.adsorbent_directory.Get(id);
          auto& src = ads.adsorption_source_;
          const std::string ap = "adsorbent[" + std::to_string(id) + "]";
          Value v(rapidjson::kObjectType);
          SetInt(v, "id", static_cast<long long>(id), a_);
          Set(v, "name", ads.name_, a_);
          Set(v, "mixing_rule", KindName(b::MixingRuleKinds(), static_cast<int>(src.mixing_rule_)), a_);

          Value comps(rapidjson::kArrayType);
          for (const auto& c : d_.fluid.components_) {
            Value cv(rapidjson::kObjectType);
            Set(cv, "component", c, a_);
            const auto iso = src.IsothermMap().find(c);
            Set(cv, "isotherm", iso != src.IsothermMap().end() ? std::string(IsothermKindString(iso->second))
                                                                 : std::string("none"), a_);
            const auto mt = src.mass_transfer_models_.find(c);
            const bool has_rate = mt != src.mass_transfer_models_.end() && mt->second && mt->second->mtc_model_;
            if (has_rate) {
              Set(cv, "rate_model", KindName(b::RateModelKinds(), static_cast<int>(mt->second->GetModelType())), a_);
              const KindDescriptor* ck = KindAt(ModelFamily::RateCoefficient,
                                                static_cast<std::size_t>(mt->second->mtc_model_->GetModelType()));
              Set(cv, "rate_coefficient", ck ? ck->tag : std::string("unknown"), a_);
            } else {
              SetNull(cv, "rate_model", a_);
              SetNull(cv, "rate_coefficient", a_);
            }
            Coefficient(cv, "", has_rate ? mt->second->mtc_model_.get() : nullptr, id);
            Coefficient(cv, "slow_", has_rate ? mt->second->SlowCoefficientModel() : nullptr, id);
            Value values(rapidjson::kArrayType);
            for (const std::string head : {".isotherm[", ".kinetics[", ".heat_of_adsorption["}) {
              const std::string prefix = ap + head + c + "]";
              for (const auto* p : pool_.Claim(prefix)) values.PushBack(ParameterValue(*p, ap + ".", a_), a_);
            }
            Set(cv, "values", std::move(values), a_);
            comps.PushBack(cv, a_);
          }
          Set(v, "components", std::move(comps), a_);

          Value reactions(rapidjson::kArrayType);
          for (auto& r : ads.reaction_source_.reactions_) {
            Value rv(rapidjson::kObjectType);
            SetInt(rv, "id", static_cast<long long>(r.id_), a_);
            Set(rv, "name", r.name_, a_);
            Set(rv, "kind", r.rate_ ? std::string(ReactionRateKindString(r.rate_->GetType())) : std::string("none"), a_);
            Set(rv, "basis", KindName(b::ReactionBasisKinds(), static_cast<int>(r.basis_)), a_);
            Value nu(rapidjson::kArrayType);
            for (const auto& [comp, coefficient] : r.nu_) {
              Value e(rapidjson::kObjectType);
              Set(e, "component", comp, a_);
              Set(e, "nu", Num(coefficient), a_);
              nu.PushBack(e, a_);
            }
            Set(rv, "stoichiometry", std::move(nu), a_);
            Set(rv, "effectiveness", KindName(b::EffectivenessKinds(), static_cast<int>(r.effectiveness_.mode_)), a_);
            Set(rv, "key_component", r.effectiveness_.key_component_, a_);
            Set(rv, "values", Values(ap + ".reaction[" + std::to_string(r.id_) + "]."), a_);
            reactions.PushBack(rv, a_);
          }
          Set(v, "reactions", std::move(reactions), a_);
          Set(v, "values", Values(ap + "."), a_);
          list.PushBack(v, a_);
        }
        return list;
      }

      Value Layers() {
        namespace b = build;
        Value list(rapidjson::kArrayType);
        for (auto& layer : d_.bed.layers_.layers_) {
          const std::string lp = "layer[" + std::to_string(layer.id_) + "]";
          Value v(rapidjson::kObjectType);
          SetInt(v, "id", static_cast<long long>(layer.id_), a_);
          Set(v, "name", layer.name_, a_);
          Set(v, "kind", KindName(b::LayerKinds(), static_cast<int>(layer.kind_)), a_);
          if (layer.kind_ == LayerKind::Monolith)
            Set(v, "channel_shape", KindName(b::ChannelShapeKinds(), static_cast<int>(layer.structured_.channel_shape_)), a_);
          else
            SetNull(v, "channel_shape", a_);
          if (layer.IsStructured())
            SetNull(v, "particle_shape", a_);
          else
            Set(v, "particle_shape", KindName(ParticleShapeKinds(), static_cast<int>(layer.particle_shape_)), a_);
          Set(v, "momentum", KindName(b::MomentumKinds(), static_cast<int>(layer.momentum_model_)), a_);
          Set(v, "film_correlation", KindName(b::FilmCorrelationKinds(), static_cast<int>(layer.film_correlation_)), a_);
          auto& le = layer.equations_;
          auto model = [&](const std::vector<std::string>& kinds, auto& ptr) -> Value {
            return ptr ? Str(KindName(kinds, static_cast<int>(ptr->GetModelType())), a_) : Value();
          };
          Set(v, "dispersion", model(b::DispersionKinds(), le.dispersion_model_), a_);
          Set(v, "heat_dispersion", model(b::HeatDispersionKinds(), le.heat_dispersion_model_), a_);
          Set(v, "internal_heat_transfer", model(b::InternalHeatTransferKinds(), le.internal_heat_transfer_model_), a_);
          Set(v, "external_heat_transfer", model(b::ExternalHeatTransferKinds(), le.external_heat_transfer_model_), a_);

          Value comp(rapidjson::kArrayType);
          for (const auto& [ads_id, fraction] : layer.composition_) {
            (void)fraction;
            Value e(rapidjson::kObjectType);
            SetInt(e, "adsorbent", static_cast<long long>(ads_id), a_);
            Set(e, "name", d_.adsorbent_directory.Contains(ads_id) ? d_.adsorbent_directory.Get(ads_id).name_
                                                                    : std::to_string(ads_id), a_);
            const std::string path = lp + ".mass_fraction[" + std::to_string(ads_id) + "]";
            auto claimed = pool_.Claim(path);
            if (!claimed.empty()) Set(e, "mass_fraction", ParameterValue(*claimed.front(), lp + ".", a_), a_);
            else SetNull(e, "mass_fraction", a_);
            comp.PushBack(e, a_);
          }
          Set(v, "composition", std::move(comp), a_);
          Set(v, "values", Values(lp + "."), a_);
          list.PushBack(v, a_);
        }
        return list;
      }

      Value Bed() {
        namespace b = build;
        Value bed(rapidjson::kObjectType);
        auto& eqs = d_.bed.equations_;
        Set(bed, "shape", KindName(b::BedShapeKinds(), static_cast<int>(d_.bed.shape_.index())), a_);
        Set(bed, "thermal_equilibrium", KindName(b::ThermalEquilibriumKinds(), static_cast<int>(eqs.thermal_equilibrium_mode_)), a_);
        Set(bed, "particle_heat_transfer",
            eqs.particle_heat_transfer_model_
                ? KindName(b::ParticleHeatTransferKinds(), static_cast<int>(eqs.particle_heat_transfer_model_->GetModelType()))
                : std::string("constant"), a_);
        Set(bed, "bottom_void", Node("bed.bottom_void."), a_);
        Set(bed, "top_void", Node("bed.top_void."), a_);
        Set(bed, "values", Values("bed."), a_);

        Value wall(rapidjson::kObjectType);
        Set(wall, "material", KindName(b::WallMaterialKinds(), static_cast<int>(d_.bed.wall_.material)), a_);
        Set(wall, "boundary",
            eqs.wall_temperature_boundary_condition_
                ? KindName(b::WallBoundaryKinds(), static_cast<int>(eqs.wall_temperature_boundary_condition_->GetModelType()))
                : std::string("zero_flux"), a_);
        Set(wall, "values", Values("wall."), a_);
        Set(bed, "wall", std::move(wall), a_);
        return bed;
      }

      std::string Profile(IBoundaryValue* v) {
        return v ? KindName(build::ProfileKinds(), static_cast<int>(v->GetType())) : std::string("constant");
      }

      // The pump of an endpoint (claims "endpoint[N].pump."), or null.
      Value Pump(IPump* p, const std::string& ep) {
        namespace b = build;
        if (!p || !p->use_pump_) return Value();
        Value q(rapidjson::kObjectType);
        Set(q, "kind", KindName(b::PumpKinds(), static_cast<int>(p->GetType()) + 1), a_);
        Set(q, "energy", p->energy_eq_ ? KindName(b::PumpEnergyKinds(), static_cast<int>(p->energy_eq_->GetType()))
                                       : std::string("isothermal"), a_);
        Set(q, "efficiency", p->efficiency_model_
                                 ? KindName(b::PumpEfficiencyKinds(), static_cast<int>(p->efficiency_model_->GetType()))
                                 : std::string("constant"), a_);
        Set(q, "values", Values(ep + ".pump."), a_);
        return q;
      }

      // The buffer tank feeding an inlet (claims "endpoint[N].buffer_tank."), or null.
      Value Tank(const BufferTankSource& t, const std::string& ep) {
        if (!t.IsActive() || t.source_endpoint_id_ == 0) return Value();
        Value v(rapidjson::kObjectType);
        SetInt(v, "source_endpoint", static_cast<long long>(t.source_endpoint_id_), a_);
        Set(v, "values", Values(ep + ".buffer_tank."), a_);
        return v;
      }

      Value End(IBoundaryCondition* bc, const char* location) {
        namespace b = build;
        Value e(rapidjson::kObjectType);
        Set(e, "location", location, a_);
        auto* pressure = dynamic_cast<PressureBoundaryCondition*>(bc);
        auto* flow = dynamic_cast<FlowBoundaryCondition*>(bc);
        // The end's tag in the boundary kind table, the one set_boundary_kind takes.
        const KindDescriptor* kind = bc ? KindAt(ModelFamily::Boundary, static_cast<std::size_t>(bc->GetType())) : nullptr;
        Set(e, "kind", kind ? kind->tag : std::string("closed"), a_);
        if (pressure || flow)
          Set(e, "direction", bc->direction_ == BoundaryDirection::INLET ? "inlet" : "outlet", a_);
        else
          SetNull(e, "direction", a_);
        Value streams(rapidjson::kArrayType);
        if (pressure) {
          const std::string ep = "endpoint[" + std::to_string(pressure->endpoint_id_) + "]";
          SetInt(e, "endpoint", static_cast<long long>(pressure->endpoint_id_), a_);
          Set(e, "profile", Profile(pressure->pressure_.get()), a_);
          Set(e, "temperature_profile", Profile(pressure->temperature_profile_.get()), a_);
          Set(e, "pump", Pump(pressure->pump_.get(), ep), a_);
          Set(e, "buffer_tank", Tank(pressure->buffer_tank_, ep), a_);
          Set(e, "values", Values(ep + "."), a_);
        } else {
          SetNull(e, "endpoint", a_);
          SetNull(e, "profile", a_);
          SetNull(e, "temperature_profile", a_);
          SetNull(e, "pump", a_);
          SetNull(e, "buffer_tank", a_);
          Set(e, "values", Value(rapidjson::kArrayType), a_);
        }
        if (flow) {
          for (auto& st : flow->streams_) {
            if (!st) continue;
            const std::string ep = "endpoint[" + std::to_string(st->endpoint_id_) + "]";
            Value s(rapidjson::kObjectType);
            SetInt(s, "endpoint", static_cast<long long>(st->endpoint_id_), a_);
            Set(s, "name", st->name_, a_);
            Set(s, "temperature_profile", Profile(st->temperature_profile_.get()), a_);
            Value flow_kind, profile, power, pump;
            // The device's tag in the stream kind table, the one set_stream_kind
            // takes ("flow", "valve", "pump"), so the view and the edit op agree.
            const KindDescriptor* kind = KindAt(ModelFamily::Stream, static_cast<std::size_t>(st->GetType()));
            const std::string type = kind ? kind->tag : std::string("unsupported");
            if (auto* bs = dynamic_cast<BoundaryStream*>(st.get())) {
              flow_kind = Str(KindName(b::FlowKinds(), static_cast<int>(bs->boundary_flow_->GetType())), a_);
              profile = Str(Profile(bs->boundary_flow_->profile_.get()), a_);
              pump = Pump(bs->pump_.get(), ep);
            } else if (auto* vs = dynamic_cast<ValveBoundaryStream*>(st.get())) {
              profile = Str(Profile(vs->boundary_value_.get()), a_);
            } else if (auto* ps = dynamic_cast<PumpBoundaryStream*>(st.get())) {
              power = Str(ps->power_ ? KindName(b::PumpPowerKinds(), static_cast<int>(ps->power_->GetType()))
                                     : std::string("constant"), a_);
              pump = Pump(ps->pump_.get(), ep);
            }
            Set(s, "type", type, a_);
            Set(s, "flow", std::move(flow_kind), a_);
            Set(s, "profile", std::move(profile), a_);
            Set(s, "power", std::move(power), a_);
            Set(s, "pump", std::move(pump), a_);
            Set(s, "buffer_tank", Tank(st->buffer_tank_, ep), a_);
            Set(s, "values", Values(ep + "."), a_);
            streams.PushBack(s, a_);
          }
        }
        Set(e, "streams", std::move(streams), a_);
        return e;
      }

      Value Steps() {
        Value list(rapidjson::kArrayType);
        const std::string mode = ChoiceOf("cycle.step_time_mode");
        size_t index = 0;
        for (auto& st : d_.cycle.step_list_) {
          const std::string sp = "step[" + std::to_string(st.id_) + "]";
          Value s(rapidjson::kObjectType);
          SetInt(s, "id", static_cast<long long>(st.id_), a_);
          SetInt(s, "index", static_cast<long long>(index++), a_);
          Set(s, "name", st.step_name_, a_);
          Set(s, "time_mode", mode, a_);
          Set(s, "bottom", End(st.bottom_boundary_condition_.get(), "bottom"), a_);
          Set(s, "top", End(st.top_boundary_condition_.get(), "top"), a_);
          Value energy(rapidjson::kArrayType);
          for (const auto& [layer_id, src] : st.energy_sources_) {
            (void)src;
            Value e(rapidjson::kObjectType);
            SetInt(e, "layer", static_cast<long long>(layer_id), a_);
            std::string layer_name = std::to_string(layer_id);
            for (const auto& l : d_.bed.layers_.layers_)
              if (l.id_ == layer_id) layer_name = l.name_;
            Set(e, "layer_name", layer_name, a_);
            Set(e, "values", Values(sp + ".energy[" + std::to_string(layer_id) + "]."), a_);
            energy.PushBack(e, a_);
          }
          Set(s, "energy_sources", std::move(energy), a_);
          Value medium(rapidjson::kObjectType);
          Set(medium, "enabled", st.heating_medium_.enabled_, a_);
          Set(medium, "temperature_profile", Profile(st.heating_medium_.temperature_profile_.get()), a_);
          Set(medium, "values", Values(sp + ".medium."), a_);
          Set(s, "heating_medium", std::move(medium), a_);
          Set(s, "values", Values(sp + "."), a_);
          list.PushBack(s, a_);
        }
        return list;
      }

      Value Connections() {
        std::map<size_t, std::vector<size_t>> by_pair;
        for (auto& h : d_.cycle.Endpoints()) {
          const Connection* c = h.stream ? &h.stream->connection_ : (h.pressure ? &h.pressure->connection_ : nullptr);
          if (c && c->IsConnected()) by_pair[c->GetID()].push_back(h.id);
        }
        Value list(rapidjson::kArrayType);
        for (auto& [pair, eps] : by_pair) {
          if (eps.size() != 2) continue;
          Value c(rapidjson::kObjectType);
          SetInt(c, "id", static_cast<long long>(pair), a_);
          Value e(rapidjson::kArrayType);
          for (size_t id : eps) e.PushBack(Value(static_cast<int64_t>(id)), a_);
          Set(c, "endpoints", std::move(e), a_);
          list.PushBack(c, a_);
        }
        return list;
      }

      Value Kpis() {
        Value list(rapidjson::kArrayType);
        for (const auto& cfg : d_.kpi_config.per_component) {
          auto ids = [&](std::vector<KPIEndpointRef> refs) {
            Value v(rapidjson::kArrayType);
            for (auto& r : refs)
              if (r.Resolve(d_.cycle.step_list_)) v.PushBack(Value(static_cast<int64_t>(r.endpoint_id)), a_);
            return v;
          };
          Value k(rapidjson::kObjectType);
          Set(k, "component", cfg.product_component, a_);
          Set(k, "tracked", cfg.tracked, a_);
          Set(k, "feed", ids(cfg.feed_endpoints), a_);
          Set(k, "product", ids(cfg.product_endpoints), a_);
          Set(k, "recycle", ids(cfg.recycle_endpoints), a_);
          Set(k, "values", Values("kpi[" + cfg.product_component + "]."), a_);
          list.PushBack(k, a_);
        }
        return list;
      }

      // The reaction KPIs: key reactant ("" when unset), products, basis and
      // their registry entries (kpi.reaction.*).
      Value ReactionKpi() {
        const ReactionKPIConfig& rx = d_.kpi_config.reaction;
        Value k(rapidjson::kObjectType);
        Set(k, "key_reactant", rx.key_reactant, a_);
        Value products(rapidjson::kArrayType);
        for (const auto& p : rx.products) products.PushBack(Value(p.c_str(), a_), a_);
        Set(k, "products", std::move(products), a_);
        Set(k, "basis", ChoiceOf("kpi.reaction.basis"), a_);
        Set(k, "values", Values("kpi.reaction."), a_);
        return k;
      }

      Value Variables() {
        Value list(rapidjson::kArrayType);
        for (const auto& var : d_.variable_registry.GetAll()) {
          Value v(rapidjson::kObjectType);
          Set(v, "name", var.name, a_);
          Set(v, "category", build::VariableCategoryKind(var.unit_category), a_);
          Set(v, "formula", var.value_expr, a_);
          Set(v, "description", var.description, a_);
          const std::string path = "variable[" + var.name + "]";
          const Parameter* p = nullptr;
          for (const auto* c : pool_.Claim(path))
            if (c->Path() == path) p = c;
          if (p) Set(v, "value", ParameterValue(*p, "", a_), a_);
          else SetNull(v, "value", a_);
          list.PushBack(v, a_);
        }
        return list;
      }
      static Value Numbers(const std::vector<double>& xs, Allocator& a) {
        Value v(rapidjson::kArrayType);
        for (double x : xs) v.PushBack(Num(x), a);
        return v;
      }

      // The datasets stored in the deck, each with its points and what it
      // measures (experimental_data.h). Not registry values: change them
      // with the experimental-data edits (deck_edit.h).
      Value ExperimentalDatasets() {
        Value list(rapidjson::kArrayType);
        const auto& all = d_.exp_data_.GetAllData();
        for (size_t i = 0; i < all.size(); ++i) {
          const XYData& xy = all[i];
          const ObservableBinding& b = xy.observable;
          Value v(rapidjson::kObjectType);
          SetInt(v, "index", static_cast<long long>(i), a_);
          Set(v, "name", xy.name, a_);
          Set(v, "x_label", xy.x_label, a_);
          Set(v, "y_label", xy.y_label, a_);
          Set(v, "x", Numbers(xy.x, a_), a_);
          Set(v, "y", Numbers(xy.y, a_), a_);
          Value o(rapidjson::kObjectType);
          Set(o, "kind", ObservableKindTags()[static_cast<size_t>(b.kind)], a_);
          Set(o, "variable", b.variable, a_);
          Set(o, "component", b.component, a_);
          Set(o, "location", ObservableLocationTags()[static_cast<size_t>(b.location)], a_);
          Set(o, "z", Num(b.z), a_);
          SetInt(o, "cycle", b.cycle, a_);
          SetInt(o, "cycle_offset", b.cycle_offset, a_);
          Set(o, "x_axis", b.XAxis() == ObservableXAxis::Cycle ? "cycle" : "time", a_);
          Set(o, "x_unit", b.x_unit, a_);
          Set(o, "x_scale", Num(b.x_scale), a_);
          Set(o, "y_unit", b.y_unit, a_);
          Set(o, "si_unit", ObservableSIUnit(b), a_);
          if (b.uncertainty.empty()) SetNull(o, "uncertainty", a_);
          else Set(o, "uncertainty", Numbers(b.uncertainty, a_), a_);
          Set(o, "weight", Num(b.weight), a_);
          Set(v, "observable", std::move(o), a_);
          list.PushBack(v, a_);
        }
        return list;
      }
    };

  }  // namespace

  std::string DeckView(UserData& data, const DeckViewOptions& options) {
    data.cycle.EnsureIds();
    data.cycle.ResolveBufferTankReferences();
    ParameterRegistry registry(data);
    Document doc;
    ViewBuilder builder(data, registry, doc.GetAllocator());
    Value root = builder.Build();
    return Write(root, options.pretty);
  }

}  // namespace adsim

// ---------------------------------------------------------------------------
// The view's JSON Schema
// ---------------------------------------------------------------------------

namespace adsim {

  namespace {

    std::string JsonList(const std::vector<std::string>& items, bool with_null = false) {
      view_json::Document doc;
      auto& a = doc.GetAllocator();
      view_json::Value v = view_json::StrList(items, a);
      if (with_null) v.PushBack(view_json::Value(), a);
      return view_json::Write(v, false);
    }

    void ReplaceAll(std::string& text, const std::string& from, const std::string& to) {
      for (size_t at = text.find(from); at != std::string::npos; at = text.find(from, at + to.size()))
        text.replace(at, from.size(), to);
    }

    // The shape of the view. @NAME@ placeholders are kind lists, filled in
    // from the builder's kind names (deck_builder.h).
    constexpr const char* kSchemaTemplate = R"JSON({
  "$schema": "http://json-schema.org/draft-04/schema#",
  "title": "DeckView",
  "description": "Read-only view of a Skarstrom deck: structure plus every editable value as a parameter-registry entry. Edit values by registry path; edit structure through the builder.",
  "type": "object",
  "additionalProperties": false,
  "required": ["view_version", "engine_version", "info", "settings", "css", "energy_balance", "components", "fluid",
               "initial_conditions", "adsorbents", "layers", "bed", "cycle", "steps", "connections", "kpis",
               "variables", "cost", "sensitivity", "experimental_data", "other"],
  "properties": {
    "view_version": {"type": "integer", "enum": [@VIEW_VERSION@]},
    "engine_version": {"type": "string"},
    "info": {"$ref": "#/definitions/Node"},
    "settings": {"$ref": "#/definitions/Node"},
    "css": {"$ref": "#/definitions/Node"},
    "energy_balance": {"$ref": "#/definitions/Node"},
    "components": {"type": "array", "items": {"$ref": "#/definitions/Component"}},
    "fluid": {"$ref": "#/definitions/Node"},
    "initial_conditions": {"$ref": "#/definitions/Node"},
    "adsorbents": {"type": "array", "items": {"$ref": "#/definitions/Adsorbent"}},
    "layers": {"type": "array", "items": {"$ref": "#/definitions/Layer"}},
    "bed": {"$ref": "#/definitions/Bed"},
    "cycle": {"$ref": "#/definitions/Cycle"},
    "steps": {"type": "array", "items": {"$ref": "#/definitions/Step"}},
    "connections": {"type": "array", "items": {"$ref": "#/definitions/Connection"}},
    "kpis": {"type": "array", "items": {"$ref": "#/definitions/KpiConfig"}},
    "reaction_kpi": {"$ref": "#/definitions/ReactionKpiConfig"},
    "variables": {"type": "array", "items": {"$ref": "#/definitions/Variable"}},
    "cost": {"$ref": "#/definitions/Node"},
    "sensitivity": {"$ref": "#/definitions/Node"},
    "experimental_data": {"type": "array", "items": {"$ref": "#/definitions/ExperimentalDataset"}},
    "other": {"$ref": "#/definitions/Node"}
  },
  "definitions": {
    "Choice": {
      "type": "object", "additionalProperties": false, "required": ["name", "index"],
      "properties": {"name": {"type": "string"}, "index": {"type": "integer", "minimum": 0}}
    },
    "Value": {
      "description": "One parameter-registry entry. Set it with its path (Parameter::SetValue / SetText).",
      "type": "object", "additionalProperties": false,
      "required": ["path", "alias", "key", "label", "category", "subcategory", "kind", "value", "text", "formula",
                   "unit", "unit_checked", "tunable", "lower", "upper", "choices", "choice"],
      "properties": {
        "path": {"type": "string", "description": "stable, id-based registry path"},
        "alias": {"type": "string", "description": "the same address with display names"},
        "key": {"type": "string", "description": "the path relative to its node"},
        "label": {"type": "string"},
        "category": {"type": "string"},
        "subcategory": {"type": "string"},
        "kind": {"type": "string", "enum": ["real", "expression", "integer", "boolean", "string"]},
        "value": {"type": ["number", "boolean", "string", "null"],
                  "description": "number (SI) for real, expression and integer; boolean; string; null when not finite"},
        "text": {"type": "string", "description": "the expression or text form"},
        "formula": {"type": "boolean", "description": "an expression entry written as a formula of user variables"},
        "unit": {"type": "string", "description": "SI unit symbol; '-' dimensionless; '' not numeric"},
        "unit_checked": {"type": "boolean"},
        "tunable": {"type": "boolean"},
        "lower": {"type": ["number", "null"]},
        "upper": {"type": ["number", "null"]},
        "choices": {"type": "array", "items": {"$ref": "#/definitions/Choice"},
                    "description": "the value names of an enumerated integer, in value order; empty otherwise"},
        "choice": {"type": ["string", "null"], "description": "the current choice name of an enumerated integer"}
      }
    },
    "ValueList": {"type": "array", "items": {"$ref": "#/definitions/Value"}},
    "Node": {
      "type": "object", "additionalProperties": false, "required": ["values"],
      "properties": {"values": {"$ref": "#/definitions/ValueList"}}
    },
    "Component": {
      "type": "object", "additionalProperties": false, "required": ["name", "formula", "molar_mass"],
      "properties": {
        "name": {"type": "string", "description": "CoolProp fluid name"},
        "formula": {"type": "string"},
        "molar_mass": {"type": ["number", "null"], "description": "kg/mol"}
      }
    },
    "AdsorbentComponent": {
      "type": "object", "additionalProperties": false,
      "required": ["component", "isotherm", "rate_model", "rate_coefficient", "resistances", "capacity_factor_mode",
                   "preset", "equation", "equation_slab", "slow_resistances", "slow_capacity_factor_mode",
                   "slow_preset", "slow_equation", "slow_equation_slab", "values"],
      "properties": {
        "component": {"type": "string"},
        "isotherm": {"type": "string", "enum": @ISOTHERM@},
        "rate_model": {"type": ["string", "null"], "enum": @RATE_MODEL@},
        "rate_coefficient": {"type": ["string", "null"], "enum": @RATE_COEFFICIENT@},
        "resistances": {"oneOf": [{"type": "null"}, {"$ref": "#/definitions/Resistances"}]},
        "capacity_factor_mode": {"type": ["string", "null"], "enum": @CAPACITY_FACTOR_MODE@},
        "preset": {"type": ["string", "null"], "enum": @RATE_PRESET@,
                   "description": "the preset the switches, mode and values match (derived), or custom"},
        "equation": {"type": ["string", "null"],
                     "description": "LaTeX of the coefficient: sphere form in a packed layer, else slab form"},
        "equation_slab": {"type": ["string", "null"],
                          "description": "the slab form when the adsorbent is also in a structured layer"},
        "slow_resistances": {"oneOf": [{"type": "null"}, {"$ref": "#/definitions/Resistances"}]},
        "slow_capacity_factor_mode": {"type": ["string", "null"], "enum": @CAPACITY_FACTOR_MODE@},
        "slow_preset": {"type": ["string", "null"], "enum": @RATE_PRESET@},
        "slow_equation": {"type": ["string", "null"], "description": "the dual-rate slow site's equation"},
        "slow_equation_slab": {"type": ["string", "null"]},
        "values": {"$ref": "#/definitions/ValueList"}
      }
    },
    "Resistances": {
      "type": "object", "additionalProperties": false,
      "required": ["film", "macropore", "molecular", "knudsen", "micropore", "lumped"],
      "properties": {"film": {"type": "boolean"}, "macropore": {"type": "boolean"}, "molecular": {"type": "boolean"},
                     "knudsen": {"type": "boolean"}, "micropore": {"type": "boolean"}, "lumped": {"type": "boolean"}}
    },
    "StoichiometricCoefficient": {
      "type": "object", "additionalProperties": false, "required": ["component", "nu"],
      "properties": {"component": {"type": "string"}, "nu": {"type": ["number", "null"]}}
    },
)JSON"
        // Split: MSVC limits one string literal to 16 KB.
        R"JSON(    "Reaction": {
      "type": "object", "additionalProperties": false,
      "required": ["id", "name", "kind", "basis", "stoichiometry", "effectiveness", "key_component", "values"],
      "properties": {
        "id": {"type": "integer"},
        "name": {"type": "string"},
        "kind": {"type": "string", "enum": @REACTION_RATE@},
        "basis": {"type": "string", "enum": @REACTION_BASIS@},
        "stoichiometry": {"type": "array", "items": {"$ref": "#/definitions/StoichiometricCoefficient"}},
        "effectiveness": {"type": "string", "enum": @EFFECTIVENESS@},
        "key_component": {"type": "string"},
        "values": {"$ref": "#/definitions/ValueList"}
      }
    },
    "Adsorbent": {
      "type": "object", "additionalProperties": false,
      "required": ["id", "name", "mixing_rule", "components", "reactions", "values"],
      "properties": {
        "id": {"type": "integer"},
        "name": {"type": "string"},
        "mixing_rule": {"type": "string", "enum": @MIXING_RULE@},
        "components": {"type": "array", "items": {"$ref": "#/definitions/AdsorbentComponent"}},
        "reactions": {"type": "array", "items": {"$ref": "#/definitions/Reaction"}},
        "values": {"$ref": "#/definitions/ValueList"}
      }
    },
    "CompositionEntry": {
      "type": "object", "additionalProperties": false, "required": ["adsorbent", "name", "mass_fraction"],
      "properties": {
        "adsorbent": {"type": "integer", "description": "adsorbent id"},
        "name": {"type": "string"},
        "mass_fraction": {"oneOf": [{"$ref": "#/definitions/Value"}, {"type": "null"}]}
      }
    },
    "Layer": {
      "type": "object", "additionalProperties": false,
      "required": ["id", "name", "kind", "channel_shape", "particle_shape", "momentum", "film_correlation", "dispersion",
                   "heat_dispersion", "internal_heat_transfer", "external_heat_transfer", "composition", "values"],
      "properties": {
        "id": {"type": "integer"},
        "name": {"type": "string"},
        "kind": {"type": "string", "enum": @LAYER@},
        "channel_shape": {"type": ["string", "null"], "enum": @CHANNEL_SHAPE@},
        "particle_shape": {"type": ["string", "null"], "enum": @PARTICLE_SHAPE@},
        "momentum": {"type": "string", "enum": @MOMENTUM@},
        "film_correlation": {"type": "string", "enum": @FILM_CORRELATION@},
        "dispersion": {"type": ["string", "null"], "enum": @DISPERSION@},
        "heat_dispersion": {"type": ["string", "null"], "enum": @HEAT_DISPERSION@},
        "internal_heat_transfer": {"type": ["string", "null"], "enum": @INTERNAL_HEAT_TRANSFER@},
        "external_heat_transfer": {"type": ["string", "null"], "enum": @EXTERNAL_HEAT_TRANSFER@},
        "composition": {"type": "array", "items": {"$ref": "#/definitions/CompositionEntry"}},
        "values": {"$ref": "#/definitions/ValueList"}
      }
    },
    "Wall": {
      "type": "object", "additionalProperties": false, "required": ["material", "boundary", "values"],
      "properties": {
        "material": {"type": "string", "enum": @WALL_MATERIAL@},
        "boundary": {"type": "string", "enum": @WALL_BOUNDARY@},
        "values": {"$ref": "#/definitions/ValueList"}
      }
    },
    "Bed": {
      "type": "object", "additionalProperties": false,
      "required": ["shape", "thermal_equilibrium", "particle_heat_transfer", "bottom_void", "top_void", "values", "wall"],
      "properties": {
        "shape": {"type": "string", "enum": @BED_SHAPE@},
        "thermal_equilibrium": {"type": "string", "enum": @THERMAL_EQUILIBRIUM@},
        "particle_heat_transfer": {"type": "string", "enum": @PARTICLE_HEAT_TRANSFER@},
        "bottom_void": {"$ref": "#/definitions/Node"},
        "top_void": {"$ref": "#/definitions/Node"},
        "values": {"$ref": "#/definitions/ValueList"},
        "wall": {"$ref": "#/definitions/Wall"}
      }
    },
    "Cycle": {
      "type": "object", "additionalProperties": false, "required": ["step_time_mode", "values"],
      "properties": {
        "step_time_mode": {"type": "string"},
        "values": {"$ref": "#/definitions/ValueList"}
      }
    },
    "Pump": {
      "type": "object", "additionalProperties": false, "required": ["kind", "energy", "efficiency", "values"],
      "properties": {
        "kind": {"type": "string", "enum": @PUMP@},
        "energy": {"type": "string", "enum": @PUMP_ENERGY@},
        "efficiency": {"type": "string", "enum": @PUMP_EFFICIENCY@},
        "values": {"$ref": "#/definitions/ValueList"}
      }
    },
    "BufferTank": {
      "type": "object", "additionalProperties": false, "required": ["source_endpoint", "values"],
      "properties": {
        "source_endpoint": {"type": "integer", "description": "endpoint id of the source"},
        "values": {"$ref": "#/definitions/ValueList"}
      }
    },
    "Stream": {
      "type": "object", "additionalProperties": false,
      "required": ["endpoint", "name", "temperature_profile", "type", "flow", "profile", "power", "pump",
                   "buffer_tank", "values"],
      "properties": {
        "endpoint": {"type": "integer"},
        "name": {"type": "string"},
        "temperature_profile": {"type": "string", "enum": @PROFILE@},
        "type": {"type": "string", "enum": @STREAM@},
        "flow": {"type": ["string", "null"], "enum": @FLOW@},
        "profile": {"type": ["string", "null"], "enum": @PROFILE_OR_NULL@},
        "power": {"type": ["string", "null"], "enum": @PUMP_POWER@},
        "pump": {"oneOf": [{"$ref": "#/definitions/Pump"}, {"type": "null"}]},
        "buffer_tank": {"oneOf": [{"$ref": "#/definitions/BufferTank"}, {"type": "null"}]},
        "values": {"$ref": "#/definitions/ValueList"}
      }
    },
    "End": {
      "description": "One end of a step. A pressure end has its endpoint, profiles, pump, buffer tank and values; a flow end has streams; a closed end neither.",
      "type": "object", "additionalProperties": false,
      "required": ["location", "kind", "direction", "endpoint", "profile", "temperature_profile", "pump",
                   "buffer_tank", "values", "streams"],
      "properties": {
        "location": {"type": "string", "enum": ["bottom", "top"]},
        "kind": {"type": "string", "enum": @BOUNDARY@},
        "direction": {"type": ["string", "null"], "enum": ["inlet", "outlet", null]},
        "endpoint": {"type": ["integer", "null"]},
        "profile": {"type": ["string", "null"], "enum": @PROFILE_OR_NULL@},
        "temperature_profile": {"type": ["string", "null"], "enum": @PROFILE_OR_NULL@},
        "pump": {"oneOf": [{"$ref": "#/definitions/Pump"}, {"type": "null"}]},
        "buffer_tank": {"oneOf": [{"$ref": "#/definitions/BufferTank"}, {"type": "null"}]},
        "values": {"$ref": "#/definitions/ValueList"},
        "streams": {"type": "array", "items": {"$ref": "#/definitions/Stream"}}
      }
    },
    "EnergySource": {
      "type": "object", "additionalProperties": false, "required": ["layer", "layer_name", "values"],
      "properties": {
        "layer": {"type": "integer", "description": "layer id"},
        "layer_name": {"type": "string"},
        "values": {"$ref": "#/definitions/ValueList"}
      }
    },
    "HeatingMedium": {
      "type": "object", "additionalProperties": false, "required": ["enabled", "temperature_profile", "values"],
      "properties": {
        "enabled": {"type": "boolean"},
        "temperature_profile": {"type": "string", "enum": @PROFILE@},
        "values": {"$ref": "#/definitions/ValueList"}
      }
    },
    "Step": {
      "type": "object", "additionalProperties": false,
      "required": ["id", "index", "name", "time_mode", "bottom", "top", "energy_sources", "heating_medium", "values"],
      "properties": {
        "id": {"type": "integer"},
        "index": {"type": "integer", "minimum": 0},
        "name": {"type": "string"},
        "time_mode": {"type": "string"},
        "bottom": {"$ref": "#/definitions/End"},
        "top": {"$ref": "#/definitions/End"},
        "energy_sources": {"type": "array", "items": {"$ref": "#/definitions/EnergySource"}},
        "heating_medium": {"$ref": "#/definitions/HeatingMedium"},
        "values": {"$ref": "#/definitions/ValueList"}
      }
    },
    "Connection": {
      "type": "object", "additionalProperties": false, "required": ["id", "endpoints"],
      "properties": {
        "id": {"type": "integer"},
        "endpoints": {"type": "array", "items": {"type": "integer"}, "minItems": 2, "maxItems": 2}
      }
    },
    "KpiConfig": {
      "type": "object", "additionalProperties": false, "required": ["component", "tracked", "feed", "product", "recycle", "values"],
      "properties": {
        "component": {"type": "string"},
        "tracked": {"type": "boolean"},
        "feed": {"type": "array", "items": {"type": "integer"}},
        "product": {"type": "array", "items": {"type": "integer"}},
        "recycle": {"type": "array", "items": {"type": "integer"}},
        "values": {"$ref": "#/definitions/ValueList"}
      }
    },
    "ReactionKpiConfig": {
      "type": "object", "additionalProperties": false, "required": ["key_reactant", "products", "basis", "values"],
      "properties": {
        "key_reactant": {"type": "string"},
        "products": {"type": "array", "items": {"type": "string"}},
        "basis": {"type": "string", "enum": ["catalyst_mass", "bed_volume"]},
        "values": {"$ref": "#/definitions/ValueList"}
      }
    },
    "ExperimentalDataset": {
      "description": "A dataset stored in the deck: its points (as measured, in the units of its observable binding) and what it measures. Change it with the experimental-data edits (add_experimental_data, set_experimental_data_binding, remove_experimental_data).",
      "type": "object", "additionalProperties": false, "required": ["index", "name", "x_label", "y_label", "x", "y", "observable"],
      "properties": {
        "index": {"type": "integer", "minimum": 0, "description": "position in the deck: the index the edits take"},
        "name": {"type": "string"},
        "x_label": {"type": "string"},
        "y_label": {"type": "string"},
        "x": {"type": "array", "items": {"type": ["number", "null"]}},
        "y": {"type": "array", "items": {"type": ["number", "null"]}},
        "observable": {"$ref": "#/definitions/Observable"}
      }
    },
    "Observable": {
      "description": "What a dataset measures. kind none: not bound. bulk/component: a bed variable (Series names) at top, bottom or z over one cycle (-1: the last logged); x is time in x_unit (x_scale seconds each). kpi: a per-cycle KPI of a tracked component; x is the cycle number, solver index = x + cycle_offset.",
      "type": "object", "additionalProperties": false,
      "required": ["kind", "variable", "component", "location", "z", "cycle", "cycle_offset", "x_axis", "x_unit", "x_scale",
                   "y_unit", "si_unit", "uncertainty", "weight"],
      "properties": {
        "kind": {"type": "string", "enum": @OBSERVABLE_KIND@},
        "variable": {"type": "string"},
        "component": {"type": "string"},
        "location": {"type": "string", "enum": @OBSERVABLE_LOCATION@},
        "z": {"type": ["number", "null"], "description": "[m] from the bottom, with location position"},
        "cycle": {"type": "integer", "minimum": -1},
        "cycle_offset": {"type": "integer"},
        "x_axis": {"type": "string", "enum": ["time", "cycle"]},
        "x_unit": {"type": "string"},
        "x_scale": {"type": ["number", "null"], "description": "seconds per x unit (1 for cycles)"},
        "y_unit": {"type": "string", "description": "unit of y; empty = si_unit"},
        "si_unit": {"type": "string", "description": "the observable's SI unit; empty when unbound"},
        "uncertainty": {"oneOf": [{"type": "null"}, {"type": "array", "items": {"type": ["number", "null"]}}],
                        "description": "per-point standard deviation in the y unit"},
        "weight": {"type": ["number", "null"]}
      }
    },
    "Variable": {
      "type": "object", "additionalProperties": false, "required": ["name", "category", "formula", "description", "value"],
      "properties": {
        "name": {"type": "string"},
        "category": {"type": "string", "enum": @VARIABLE_CATEGORY@},
        "formula": {"type": "string"},
        "description": {"type": "string"},
        "value": {"oneOf": [{"$ref": "#/definitions/Value"}, {"type": "null"}]}
      }
    }
  }
})JSON";

  }  // namespace

  std::string DeckViewSchema(const DeckViewOptions& options) {
    namespace b = build;
    std::string text = kSchemaTemplate;
    std::vector<std::string> isotherms = b::IsothermKinds();
    isotherms.push_back("none");
    std::vector<std::string> pumps;  // PumpKinds without "none": a view pump is in use
    for (const auto& k : b::PumpKinds())
      if (k != "none") pumps.push_back(k);
    const std::pair<const char*, std::string> lists[] = {
      {"@VIEW_VERSION@", std::to_string(kDeckViewVersion)},
      {"@ISOTHERM@", JsonList(isotherms)},
      {"@RATE_MODEL@", JsonList(b::RateModelKinds(), true)},
      {"@RATE_COEFFICIENT@", JsonList(b::RateCoefficientKinds(), true)},
      {"@CAPACITY_FACTOR_MODE@", JsonList(std::vector<std::string>(std::begin(CapacityFactorModeList), std::end(CapacityFactorModeList)), true)},
      {"@RATE_PRESET@", JsonList([] { auto t = b::RateCoefficientPresets(); t.push_back("custom"); return t; }(), true)},
      {"@REACTION_RATE@", JsonList(b::ReactionRateKinds())},
      {"@REACTION_BASIS@", JsonList(b::ReactionBasisKinds())},
      {"@EFFECTIVENESS@", JsonList(b::EffectivenessKinds())},
      {"@MIXING_RULE@", JsonList(b::MixingRuleKinds())},
      {"@LAYER@", JsonList(b::LayerKinds())},
      {"@CHANNEL_SHAPE@", JsonList(b::ChannelShapeKinds(), true)},
      {"@PARTICLE_SHAPE@", JsonList(ParticleShapeKinds(), true)},
      {"@MOMENTUM@", JsonList(b::MomentumKinds())},
      {"@FILM_CORRELATION@", JsonList(b::FilmCorrelationKinds())},
      {"@DISPERSION@", JsonList(b::DispersionKinds(), true)},
      {"@HEAT_DISPERSION@", JsonList(b::HeatDispersionKinds(), true)},
      {"@INTERNAL_HEAT_TRANSFER@", JsonList(b::InternalHeatTransferKinds(), true)},
      {"@EXTERNAL_HEAT_TRANSFER@", JsonList(b::ExternalHeatTransferKinds(), true)},
      {"@WALL_MATERIAL@", JsonList(b::WallMaterialKinds())},
      {"@WALL_BOUNDARY@", JsonList(b::WallBoundaryKinds())},
      {"@BED_SHAPE@", JsonList(b::BedShapeKinds())},
      {"@THERMAL_EQUILIBRIUM@", JsonList(b::ThermalEquilibriumKinds())},
      {"@PARTICLE_HEAT_TRANSFER@", JsonList(b::ParticleHeatTransferKinds())},
      {"@PUMP@", JsonList(pumps)},
      {"@PUMP_ENERGY@", JsonList(b::PumpEnergyKinds())},
      {"@PUMP_EFFICIENCY@", JsonList(b::PumpEfficiencyKinds())},
      {"@PUMP_POWER@", JsonList(b::PumpPowerKinds(), true)},
      {"@FLOW@", JsonList(b::FlowKinds(), true)},
      {"@BOUNDARY@", JsonList(KindTags(ModelFamily::Boundary))},
      {"@STREAM@", JsonList([] { auto t = KindTags(ModelFamily::Stream); t.push_back("unsupported"); return t; }())},
      {"@PROFILE_OR_NULL@", JsonList(b::ProfileKinds(), true)},
      {"@PROFILE@", JsonList(b::ProfileKinds())},
      {"@VARIABLE_CATEGORY@", JsonList(b::VariableCategoryKinds())},
      {"@OBSERVABLE_KIND@", JsonList(ObservableKindTags())},
      {"@OBSERVABLE_LOCATION@", JsonList(ObservableLocationTags())},
    };
    for (const auto& [placeholder, list] : lists) ReplaceAll(text, placeholder, list);
    view_json::Document doc;
    doc.Parse(text.c_str());
    if (doc.HasParseError()) throw std::logic_error("DeckViewSchema: the schema template is not valid JSON");
    return view_json::Write(doc, options.pretty);
  }

}  // namespace adsim
