/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/objects/reaction.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>
#include <limits>
#include <stdexcept>

#include "adsim/math/constants.h"
#include "adsim/objects/fluid.h"
#include "adsim/objects/porous_media.h"
#include "adsim/objects/reactor.h"

namespace adsim {

  const char* HeatOfReactionModeString(HeatOfReactionMode m) {
    return m == HeatOfReactionMode::Constant ? "constant" : "kirchhoff";
  }

  const char* EffectivenessModeString(EffectivenessMode m) {
    return m == EffectivenessMode::Thiele ? "thiele" : "constant";
  }

  // ---- Reaction ----------------------------------------------------------------

  void Reaction::SetStoichiometry(const std::map<std::string, double>& nu) {
    nu_.clear();
    for (const auto& [c, v] : nu)
      if (v != 0.0) nu_[c] = v;
    if (rate_) rate_->SeedComponents(nu_);
  }

  void Reaction::SetRateType(ReactionRateType type) {
    auto next = CreateReactionRate(type, nu_);
    if (rate_) {
      // LHHW adsorption terms are structural: carry their count over so their
      // leaves exist in the new layout before the values are copied.
      auto* old_lhhw = dynamic_cast<LhhwRate*>(rate_.get());
      auto* new_lhhw = dynamic_cast<LhhwRate*>(next.get());
      if (old_lhhw && new_lhhw) {
        new_lhhw->adsorption_.resize(old_lhhw->adsorption_.size());
        new_lhhw->reversible_ = old_lhhw->reversible_;
      }
      // Orders are per component: make sure the new rate has every order key
      // of the old one before copying.
      std::map<std::string, ExpressionParam>* old_orders = nullptr;
      std::map<std::string, ExpressionParam>* new_orders = nullptr;
      auto orders_of = [](IReactionRate* r) -> std::map<std::string, ExpressionParam>* {
        if (auto* p = dynamic_cast<PowerLawRate*>(r)) return &p->order_;
        if (auto* p = dynamic_cast<ReversiblePowerLawRate*>(r)) return &p->order_;
        if (auto* p = dynamic_cast<LhhwRate*>(r)) return &p->order_;
        return nullptr;
      };
      old_orders = orders_of(rate_.get());
      new_orders = orders_of(next.get());
      if (old_orders && new_orders)
        for (const auto& [c, n] : *old_orders) new_orders->try_emplace(c, 0.0);
      if (old_lhhw && new_lhhw)
        for (size_t j = 0; j < old_lhhw->adsorption_.size(); ++j)
          for (const auto& [c, e] : old_lhhw->adsorption_[j].exponent_)
            new_lhhw->adsorption_[j].exponent_.try_emplace(c, 0.0);

      auto old_fields = ReactionRateFields(rate_.get());
      for (auto& [leaf, param] : ReactionRateFields(next.get())) {
        for (auto& [old_leaf, old_param] : old_fields) {
          if (old_leaf != leaf) continue;
          *param = old_param->Get();
          param->SetExpr(old_param->GetExpr());
          break;
        }
      }
    }
    rate_ = std::move(next);
  }

  std::vector<std::string> Reaction::Components() const {
    std::set<std::string> out;
    for (const auto& [c, v] : nu_) out.insert(c);
    if (rate_) {
      for (auto& [leaf, param] : ReactionRateFields(rate_.get())) {
        (void)param;
        // Per-component leaves end in "[Component]".
        if (leaf.empty() || leaf.back() != ']') continue;
        const size_t open = leaf.rfind('[');
        if (open == std::string::npos) continue;
        const std::string head = leaf.substr(0, open);
        if (head.size() >= 5 && (head.compare(head.size() - 5, 5, "order") == 0 ||
                                 (head.size() >= 8 && head.compare(head.size() - 8, 8, "exponent") == 0)))
          out.insert(leaf.substr(open + 1, leaf.size() - open - 2));
      }
    }
    return {out.begin(), out.end()};
  }

  bool Reaction::UsesComponent(const std::string& component) const {
    auto it = nu_.find(component);
    if (it != nu_.end() && it->second != 0.0) return true;
    return rate_ && rate_->UsesComponent(component);
  }

  bool Reaction::RemoveComponent(const std::string& component) {
    if (UsesComponent(component)) return false;
    nu_.erase(component);
    if (rate_) rate_->RemoveComponent(component);
    return true;
  }

  bool Reaction::IsActive(const std::vector<std::string>& components) const {
    for (const auto& c : Components()) {
      if (!UsesComponent(c)) continue;
      if (std::find(components.begin(), components.end(), c) == components.end()) return false;
    }
    return true;
  }

  bool Reaction::HasReactant() const {
    return std::any_of(nu_.begin(), nu_.end(), [](const auto& e) { return e.second < 0.0; });
  }

  bool Reaction::HasProduct() const {
    return std::any_of(nu_.begin(), nu_.end(), [](const auto& e) { return e.second > 0.0; });
  }

  ElementBalance Reaction::CheckElements(double tolerance) const {
    ElementBalance out;
    std::map<std::string, double> sums;
    for (const auto& [c, v] : nu_) {
      const auto* formula = Components::GetElementFormula(c);
      if (!formula) {
        out.unknown.push_back(c);
        continue;
      }
      for (const auto& [element, atoms] : *formula) sums[element] += v * atoms;
    }
    if (!out.unknown.empty()) return out;
    out.checked = true;
    for (const auto& [element, sum] : sums)
      if (std::abs(sum) > tolerance) out.imbalance[element] = sum;
    return out;
  }

  double Reaction::IntrinsicRate(const ReactionState& state) {
    return rate_ ? rate_->Rate(state, basis_, nu_) : 0.0;
  }

  std::string Reaction::KeyReactant() const {
    auto it = nu_.find(effectiveness_.key_component_);
    if (it != nu_.end() && it->second < 0.0) return it->first;
    for (const auto& [c, v] : nu_)
      if (v < 0.0) return c;
    return {};
  }

  namespace {
    // 8-point Gauss-Legendre nodes and weights on [-1, 1].
    constexpr double kGaussNodes[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
                                       -0.1834346424956498, 0.1834346424956498,  0.5255324099163290,
                                       0.7966664774136267,  0.9602898564975363};
    constexpr double kGaussWeights[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
                                         0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
                                         0.2223810344533745, 0.1012285362903763};

    // eta of a sphere from its generalised modulus phi:
    // (1/phi) (1/tanh(3 phi) - 1/(3 phi)), with its series near 0.
    double SphereEffectiveness(double phi) {
      if (!(phi > 0.0)) return 1.0;
      if (phi < 1e-3) return 1.0 - 0.6 * phi * phi;
      const double x = 3.0 * phi;
      if (x > 40.0) return (1.0 - 1.0 / x) / phi;  // coth(x) = 1 to double precision
      return (1.0 / std::tanh(x) - 1.0 / x) / phi;
    }
  }  // namespace

  double Reaction::Effectiveness(const ReactionState& surface, double particle_radius, double particle_density,
                                 double effective_diffusivity, double* warm_start) {
    if (effectiveness_.mode_ != EffectivenessMode::Thiele) return effectiveness_.effectiveness_.Get();
    const std::string key = KeyReactant();
    if (key.empty() || !rate_ || !surface.y) return 1.0;
    if (!(particle_radius > 0.0) || !(particle_density > 0.0) || !(effective_diffusivity > 0.0)) return 1.0;
    const double nu_A = nu_.at(key);
    const double T = std::max(surface.T, constants::kMinDenomTemperature);
    const double C_total = surface.P / (constants::kIdealGasConstant * T);
    if (!(C_total > 0.0) || !std::isfinite(C_total)) return 1.0;

    // The gas inside the particle: the surface composition, with the
    // stoichiometry's components moving with A.
    std::map<std::string, double> y_line = *surface.y;
    struct Line { double* y; double C_s; double ratio; };
    std::vector<Line> line;
    line.reserve(nu_.size());
    double C_A_s = 0.0;
    for (const auto& [c, v] : nu_) {
      auto it = y_line.find(c);
      const double y_s = it == y_line.end() ? 0.0 : it->second;
      const double C_s = std::max(y_s, 0.0) * C_total;
      line.push_back({&y_line[c], C_s, v / nu_A});
      if (c == key) C_A_s = C_s;
    }
    ReactionState state = surface;
    state.y = &y_line;
    auto rate_at = [&](double C_A) {
      for (const auto& l : line) *l.y = std::max(l.C_s + l.ratio * (C_A - C_A_s), 0.0) / C_total;
      return rate_->Rate(state, basis_, nu_);
    };

    const double r_s = rate_at(C_A_s);
    if (!(r_s != 0.0) || !std::isfinite(r_s)) return 1.0;
    const double dir = r_s > 0.0 ? 1.0 : -1.0;

    // Where every concentration on the line is >= 0.
    double lo = 0.0, hi = std::numeric_limits<double>::infinity();
    for (const auto& l : line) {
      if (l.ratio > 0.0) lo = std::max(lo, C_A_s - l.C_s / l.ratio);
      else if (l.ratio < 0.0) hi = std::min(hi, C_A_s - l.C_s / l.ratio);
    }
    // a: the surface side (rate of the surface's sign); b: the far bound.
    double a = C_A_s;
    double b = dir > 0.0 ? lo : hi;
    if (!std::isfinite(b)) return 1.0;
    double fa = r_s;
    double fb = rate_at(b);
    double C_eq = b;
    if (fb * dir < 0.0) {
      // The rate changes sign on [b, a]: the equilibrium concentration,
      // warm-started from the previous root of this cell when it lies inside.
      if (warm_start && std::isfinite(*warm_start) && (*warm_start - a) * (*warm_start - b) < 0.0) {
        const double fw = rate_at(*warm_start);
        if (fw * dir > 0.0) { a = *warm_start; fa = fw; }
        else if (fw * dir < 0.0) { b = *warm_start; fb = fw; }
        else { a = b = *warm_start; }
      }
      // Illinois regula falsi.
      int side = 0;
      const double tol = 1e-13 * std::max(std::abs(C_A_s), std::abs(b)) + 1e-300;
      for (int it = 0; it < 100 && std::abs(a - b) > tol; ++it) {
        const double c = (a * fb - b * fa) / (fb - fa);
        if (!std::isfinite(c)) break;
        const double fc = rate_at(c);
        if (fc == 0.0) { a = b = c; break; }
        if (fc * dir > 0.0) {
          a = c; fa = fc;
          if (side == 1) fb *= 0.5;
          side = 1;
        } else {
          b = c; fb = fc;
          if (side == -1) fa *= 0.5;
          side = -1;
        }
      }
      C_eq = 0.5 * (a + b);
      if (warm_start) *warm_start = C_eq;
    }

    // integral_{C_eq}^{C_A,s} R_A dC (> 0 either way round).
    const double scale = std::abs(nu_A) * particle_density;
    const double half = 0.5 * (C_A_s - C_eq);
    const double mid = 0.5 * (C_A_s + C_eq);
    double integral = 0.0;
    for (int k = 0; k < 8; ++k) integral += kGaussWeights[k] * rate_at(mid + half * kGaussNodes[k]);
    integral *= half * scale;
    if (!(integral > 0.0) || !std::isfinite(integral)) return 1.0;
    const double phi = (particle_radius / 3.0) * scale * std::abs(r_s) / std::sqrt(2.0 * effective_diffusivity * integral);
    return SphereEffectiveness(phi);
  }

  double Reaction::DeltaHeatCapacity(const std::map<std::string, double>& molar_cp) const {
    double dcp = 0.0;
    for (const auto& [c, v] : nu_) {
      auto it = molar_cp.find(c);
      if (it != molar_cp.end()) dcp += v * it->second;
    }
    return dcp;
  }

  double Reaction::HeatOfReaction(double T, const std::map<std::string, double>& molar_cp) const {
    const double h_ref = heat_of_reaction_.Get();
    if (heat_mode_ == HeatOfReactionMode::Constant) return h_ref;
    return h_ref + DeltaHeatCapacity(molar_cp) * (T - heat_reference_temperature_.Get());
  }

  double Reaction::HeatOfReaction(double T, const Fluid& fluid) const {
    const double h_ref = heat_of_reaction_.Get();
    if (heat_mode_ == HeatOfReactionMode::Constant) return h_ref;
    // Kirchhoff with the gas heat capacities of the energy balance:
    // dH_r(T) = dH_ref + sum_i nu_i integral_{T_ref}^{T} Cp_i dT, exact for the
    // ideal-gas polynomials (Fluid::MolarEnthalpyChange) and dCp (T - T_ref)
    // for held values.
    const double T_ref = heat_reference_temperature_.Get();
    double dh = h_ref;
    for (const auto& [c, v] : nu_) {
      if (fluid.heat_capacity_.count(c) && fluid.molar_mass_.count(c)) dh += v * fluid.MolarEnthalpyChange(c, T_ref, T);
    }
    return dh;
  }

  void Reaction::EvaluateExpressions(const EvalContext& ctx) {
    heat_of_reaction_.Evaluate(ctx);
    heat_reference_temperature_.Evaluate(ctx);
    effectiveness_.EvaluateExpressions(ctx);
    if (!rate_) return;
    if (rate_->GetType() == ReactionRateType::Expression && !ctx.audit) {
      // Outside a run there is no gas state: bind the per-component names of
      // the stoichiometry to an equimolar placeholder at ctx's T and P so the
      // formula parses; its value there is not meaningful.
      std::vector<ComponentStateValue> comps;
      const double n = nu_.empty() ? 1.0 : static_cast<double>(nu_.size());
      const double c_total = ctx.P / (8.31446261815324 * std::max(ctx.T, 1e-6));
      for (const auto& [c, v] : nu_) comps.push_back({c, ctx.P / n, c_total / n, 1.0 / n});
      EvalContext local = ctx;
      local.components = &comps;
      rate_->EvaluateExpressions(local);
      return;
    }
    rate_->EvaluateExpressions(ctx);
  }

  // ---- ReactionSource ------------------------------------------------------------

  ReactionID ReactionSource::Add(const std::string& name, const std::map<std::string, double>& nu,
                                 ReactionRateType type) {
    Reaction r(name, ++next_id_);
    r.SetStoichiometry(nu);
    r.rate_ = CreateReactionRate(type, r.nu_);
    reactions_.push_back(std::move(r));
    return next_id_;
  }

  bool ReactionSource::Remove(ReactionID id) {
    auto it = std::find_if(reactions_.begin(), reactions_.end(), [id](const Reaction& r) { return r.id_ == id; });
    if (it == reactions_.end()) return false;
    reactions_.erase(it);
    return true;
  }

  Reaction* ReactionSource::Find(ReactionID id) {
    for (auto& r : reactions_)
      if (r.id_ == id) return &r;
    return nullptr;
  }

  const Reaction* ReactionSource::Find(ReactionID id) const {
    return const_cast<ReactionSource*>(this)->Find(id);
  }

  Reaction& ReactionSource::Get(ReactionID id) {
    if (Reaction* r = Find(id)) return *r;
    throw std::out_of_range("no reaction with ID " + std::to_string(id));
  }

  const Reaction& ReactionSource::Get(ReactionID id) const {
    return const_cast<ReactionSource*>(this)->Get(id);
  }

  std::vector<std::string> ReactionSource::ReactionsUsing(const std::string& component) const {
    std::vector<std::string> out;
    for (const auto& r : reactions_)
      if (r.UsesComponent(component)) out.push_back(r.name_);
    return out;
  }

  void ReactionSource::RemoveComponent(const std::string& component) {
    for (auto& r : reactions_) r.RemoveComponent(component);
  }

  void ReactionSource::Initialise(const Fluid& /*fluid*/) {
    for (auto& r : reactions_)
      if (r.rate_) r.rate_->SeedComponents(r.nu_);
  }

  void ReactionSource::EvaluateExpressions(const EvalContext& ctx) {
    for (auto& r : reactions_) r.EvaluateExpressions(ctx);
  }

  // ---- checks ----------------------------------------------------------------------

  bool ReactionChecks(const AdsorbentDirectory& dir, const Bed& bed, const Fluid& fluid, std::string& err) {
    std::set<AdsorbentID> used;
    for (const auto& layer : bed.layers_.layers_)
      for (const auto& [id, frac] : layer.composition_)
        if (dir.Contains(id)) used.insert(id);

    bool ok = true;
    auto fail = [&](const std::string& msg) {
      err += "Error: " + msg + "\n";
      ok = false;
    };
    const auto& comps = fluid.components_;
    auto in_fluid = [&](const std::string& c) { return std::find(comps.begin(), comps.end(), c) != comps.end(); };

    for (AdsorbentID id : used) {
      const Adsorbent& ads = dir.Get(id);
      for (const auto& r : ads.reaction_source_.reactions_) {
        const std::string who = "reaction '" + r.name_ + "' (on " + ads.name_ + ")";
        if (!r.rate_) {
          fail(who + " has no rate law.");
          continue;
        }
        for (const auto& c : r.Components()) {
          if (r.UsesComponent(c) && !in_fluid(c))
            fail(who + " uses " + c + ", which is not a component of the deck. Add the component or remove the reaction.");
        }
        if (!r.HasReactant()) fail(who + " has no reactant (no negative stoichiometric coefficient).");
        if (!r.HasProduct()) fail(who + " has no product (no positive stoichiometric coefficient).");
        const ElementBalance eb = r.CheckElements();
        if (eb.checked && !eb.imbalance.empty()) {
          std::ostringstream ss;
          ss << who << " does not balance the elements:";
          for (const auto& [element, sum] : eb.imbalance) ss << " " << element << " " << (sum > 0 ? "+" : "") << sum;
          ss << " (sum of nu times atoms per molecule).";
          fail(ss.str());
        }
        const bool reversible = r.rate_->GetType() == ReactionRateType::ReversiblePowerLaw ||
            (r.rate_->GetType() == ReactionRateType::Lhhw && static_cast<const LhhwRate*>(r.rate_.get())->reversible_);
        for (auto& [leaf, param] : ReactionRateFields(r.rate_.get())) {
          const bool pre_exp = leaf.size() >= 14 && leaf.compare(leaf.size() - 14, 14, "pre_exp_factor") == 0;
          if (!pre_exp) continue;
          const double v = param->Get();
          if (leaf == "equilibrium_pre_exp_factor") {
            if (reversible && !(std::isfinite(v) && v > 0.0))
              fail(who + ": the equilibrium pre-exponential factor must be finite and greater than 0.");
          } else if (!(std::isfinite(v) && v >= 0.0)) {
            fail(who + ": the pre-exponential factor " + leaf + " must be finite and not negative.");
          }
        }
        if (r.effectiveness_.mode_ == EffectivenessMode::Thiele) {
          const auto& key = r.effectiveness_.key_component_;
          auto it = r.nu_.find(key);
          if (key.empty() || it == r.nu_.end() || !(it->second < 0.0))
            fail(who + ": the Thiele effectiveness model needs a key reactant (a component with a negative coefficient).");
        }
      }
    }
    return ok;
  }

}  // namespace adsim
