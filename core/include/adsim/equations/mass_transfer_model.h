/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <string>
#include <memory>
#include <cmath>
#include <vector>

// Cereal includes
#include <cereal/types/memory.hpp>
#include <cereal/types/string.hpp>
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>

// adsim includes
#include "adsim/equations/transport_properties.h"
#include "adsim/math/expression_param.h"
#include "adsim/math/user_variables.h"


namespace adsim {

	// Forward declarations
	class PorousMedia;
	class Adsorbent;
  class MixedAdsorbentLayer;
  class Equations;
  class Fluid;
  class IIsotherm;
  class AdsorptionSource;
  using AdsorbentID = std::size_t; // mirrors porous_media.h; identical typedef redeclaration is legal
	
	enum class MassTransferModelType {
		LDF,
		LDF2,
		DUAL_RATE_LDF  // appended: stored enum values never change
	};

	static const char* RateModelTypeList[] =  {
		"Linear Driving Force (1st Order)",
		"Linear Driving Force (2nd Order)",
		"Dual-Rate Linear Driving Force"
	};

	enum class MassTransferCoefficientType {
		CONSTANT, 
		MACROPORE, 
		MACROPORE_KNUDSEN, 
		HAYNES_1973,
		ARRHENIUS,
		MACROPORE_CAPACITY,
		SERIES_RESISTANCE  // appended: stored enum values never change
	};

	static const char* MassTransferCoefficientTypeList[] = {
		"Constant",		
		"Macropore Diffusion",
		"Macropore + Knudsen Diffusion",
		"Haynes (1973)",
		"Arrhenius",
		"Macropore Diffusion (capacity factor)",
		"Mass-transfer resistances"
	};

	class IMassTransferCoefficientModel {
	public:
		virtual ~IMassTransferCoefficientModel() = default;
		virtual MassTransferCoefficientType GetModelType() = 0;
		virtual unsigned ReadsTransportProperties() const { return kAllTransportProperties; }
		virtual void UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) = 0;
		// The adsorbent's isotherm for this component and its particle density
		// [kg/m3], handed over by AdsorptionSource before every update. Only the
		// capacity-factor model reads them (it evaluates the isotherm itself).
		// The component's isotherm, the adsorbent's solid density and the
		// adsorbent's AdsorptionSource (whose UpdateIsothermAll applies its
		// mixing rule), for the models that evaluate q* themselves.
		virtual void SetEquilibriumContext(IIsotherm* /*isotherm*/, double /*solid_density*/,
		                                   AdsorptionSource* /*source*/) {}
		virtual void EvaluateExpressions(const EvalContext& ctx) {}
		// Compiles the parameters evaluated per cell in the residual (see
		// ExpressionParam::CompileForRun). Called once at the start of a run.
		virtual void CompileForRun(const UserVariableRegistry* registry) {}

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar) {}
	};

	// How the generic coefficient model scales its gas-side resistances (film,
	// macropore) by the capacity factor K, the ratio of the solid's to the
	// gas's capacity for the component. Stored as its integer value.
	enum class CapacityFactorMode {
		NONE = 0,    // K = 1
		SLOPE = 1,   // K = max(rho_p dq*/dc, 1) at max(c, c_floor); the secant where the slope is <= 0
		SECANT = 2,  // K = max(rho_p q*/c, 1) at max(c, c_floor); 1 where rho_p q* <= 0
		FIXED = 3    // K = capacity_factor (not floored)
	};

	static const char* CapacityFactorModeList[] = {"none", "slope", "secant", "fixed"};

	// One switch per resistance of the generic coefficient model. Knudsen
	// diffusion refines the macropore term (Bosanquet: 1/D = 1/D_K + 1/D_m) and
	// is on only with it.
	struct MassTransferResistances {
		bool film = false;
		bool macropore = false;
		bool knudsen = false;
		bool micropore = false;
		bool lumped = true;

		bool operator==(const MassTransferResistances& o) const {
			return film == o.film && macropore == o.macropore && knudsen == o.knudsen && micropore == o.micropore &&
			       lumped == o.lumped;
		}
		bool operator!=(const MassTransferResistances& o) const { return !(*this == o); }

		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const /*version*/) {
			ar(cereal::make_nvp("Film", film), cereal::make_nvp("Macropore", macropore),
			   cereal::make_nvp("Knudsen", knudsen), cereal::make_nvp("Micropore", micropore),
			   cereal::make_nvp("Lumped", lumped));
		}
	};

	// The resistance tags, in the fixed series order of the model.
	inline const std::vector<std::string>& ResistanceTags() {
		static const std::vector<std::string> tags{"film", "macropore", "knudsen", "micropore", "lumped"};
		return tags;
	}

	// The generic LDF coefficient: resistances in series, each switched on or
	// off on its own (Glueckauf 1955; Ruthven 1984):
	//   1/k = (f_film K) R_film + K R_macro + R_micro + R_lumped,
	// a switched-off term contributing nothing. For packed spheres (radius r_p)
	// and a slab of half-thickness l (structured layer):
	//   R_film  = r_p / (3 k_f)  |  l / k_f,   k_f = Sh D_m / d_p  |  Sh D_m / D_h
	//             (the layer's film correlation; the state's D_m);
	//   R_macro = tau r_p^2 / (15 eps_p D)  |  tau l^2 / (3 eps_p D),
	//             D the constant molecular_diffusion_coefficient when > 0, else
	//             the state's D_m; with Knudsen diffusion D -> 1/(1/D_K + 1/D),
	//             D_K = (2/3) r_pore sqrt(8 R T / (pi M));
	//   R_micro = r_c^2 / (15 D_c),  D_c = D_c0 exp(-E_a,c / (R T));
	//   R_lumped = 1 / k_L,  k_L = k exp(-E_a / (R T)) (k when E_a = 0).
	// K multiplies the gas-side terms (film, macropore) only; how it is taken
	// is the capacity-factor mode (CapacityFactorMode). As in
	// Stampi-Bombelli et al. 2024 (Ind. Eng. Chem. Res. 63, 11637,
	// doi:10.1021/acs.iecr.4c01392) for direct-air-capture breakthrough curves.
	//
	// Every older coefficient kind is a preset of this model (PresetType): a
	// switch set, a mode and its parameter values. The single-term presets
	// are evaluated in the operation order of the classes they replace
	// (mass_transfer_legacy.h), so each gives the same double as before:
	//   lumped only:                         k = k_L (constant: k)
	//   macropore (+ knudsen), mode none:    k = c D_eff / l^2 (Glueckauf c = 60 on d_p, 3 on l)
	//   macropore, mode secant:              that times min(c / (rho_p q*), 1)
	//   anything else:                       the series sum above.
	// Every parameter is stored whatever its switch, so switching a term off
	// and on again restores its values.
	//
	// The tag "series_resistance" and the cereal name are those of the first
	// version, in which a term was switched off by a zero parameter
	// (film_resistance_factor, tortuosity, crystal_radius or
	// micropore_diffusivity) and K came from the isotherm slope unless
	// capacity_factor > 0. A version-1 archive loads with the switches and
	// mode derived from those zeros (AdoptVersion1Switches).
	class MassTransferCoefficientSeriesResistance : public IMassTransferCoefficientModel {
	public:
		MassTransferCoefficientType GetModelType() override { return MassTransferCoefficientType::SERIES_RESISTANCE; }
		// The preset (former kind) the switches, mode and values match, or
		// SERIES_RESISTANCE when none does; derived, never stored.
		MassTransferCoefficientType PresetType() const;
		// The preset's tag, or "custom" when no older kind matches.
		std::string PresetTag() const;

		unsigned ReadsTransportProperties() const override;
		void UpdateMassTransferCoefficient(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component, MixedAdsorbentLayer& layer, const UserVariableRegistry* registry) override;
		void SetEquilibriumContext(IIsotherm* isotherm, double solid_density, AdsorptionSource* source) override {
			isotherm_ = isotherm; solid_density_ = solid_density; source_ = source;
		}
		void EvaluateExpressions(const EvalContext& ctx) override;
		void CompileForRun(const UserVariableRegistry* registry) override {
			k_.CompileForRun(registry); activation_energy_.CompileForRun(registry);
		}

		// ---- switches and capacity factor ----
		// Whether a resistance ("film", "macropore", "knudsen", "micropore",
		// "lumped") is on; false for an unknown tag.
		bool Enabled(const std::string& tag) const;
		// The switch set.
		MassTransferResistances Resistances() const { return resistances_; }
		// The stored switch of a resistance; throws std::out_of_range for an
		// unknown tag.
		bool& Switch(const std::string& tag);
		// Sets one switch with the model's rules: at least one of film,
		// macropore, micropore and lumped stays on; knudsen needs macropore,
		// and switching macropore off clears it. Parameter values are kept,
		// except that switching a term on replaces those of its parameters
		// that are out of range (a version-1 zero, say) by the term's defaults.
		// Returns false (err set) when refused.
		bool SetResistance(const std::string& tag, bool on, std::string& err);
		CapacityFactorMode CapacityMode() const;
		// Sets the mode; fixed with capacity_factor <= 0 sets it to 1.
		void SetCapacityMode(CapacityFactorMode mode);
		// K matters only with a gas-side term (film or macropore) on.
		bool UsesCapacityFactor() const { const auto r = Resistances(); return r.film || r.macropore; }
		// The parameters of the current mode: concentration_floor (slope,
		// secant), capacity_factor (fixed), none (none); empty without a
		// gas-side term.
		std::vector<std::string> CapacityParameters() const;
		// The parameters of one resistance, in display order (empty for an
		// unknown tag).
		static const std::vector<std::string>& ResistanceParameters(const std::string& tag);
		// Sets the switches and mode a version-1 model meant: a term is on when
		// its parameter is not zero (film_resistance_factor, tortuosity,
		// crystal_radius and micropore_diffusivity), K from the isotherm slope
		// unless capacity_factor > 0. A switched-off term's zero sentinels then
		// become that term's defaults, which cannot change results: the term
		// is off. Used when a version-1 archive is loaded.
		void AdoptVersion1Switches();
		// Sets the switches, mode and the values the preset fixes (constant:
		// activation_energy = 0; haynes: film_resistance_factor = 1,
		// molecular_diffusion_coefficient = 0). SERIES_RESISTANCE (custom)
		// changes nothing. Other parameters are left as they are.
		void ApplyPreset(MassTransferCoefficientType preset);

		// The capacity factor K at each cell of [start, end] (written into k_out)
		// for the current mode. Public for the unit tests.
		void CapacityFactor(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component,
		                    int start, int end, const UserVariableRegistry* registry, std::vector<double>& k_out);

		// Defaults (also used for a first enable). A new model has the lumped
		// term only, k = 0.5 1/s (the former constant kind's default).
		ExpressionParam film_resistance_factor_{1.0};       // f_film [-], multiplies the film term
		ExpressionParam tortuosity_{3.0};                   // particle tortuosity [-]
		ExpressionParam macropore_void_fraction_{0.35};     // eps_p [-], the macropore term's own
		// Constant molecular diffusivity D [m2/s] of the macropore term when
		// > 0; 0 (the default) keeps the state's binary-diffusion value. The
		// film term always uses the state's D_m, as the film correlation does.
		ExpressionParam molecular_diffusion_coefficient_{0.0};
		ExpressionParam crystal_radius_{1e-6};              // r_c [m]
		ExpressionParam micropore_diffusivity_{1e-13};      // D_c0 [m2/s]
		ExpressionParam micropore_activation_energy_{0.0};  // E_a,c [J/mol]
		ExpressionParam capacity_factor_{0.0};              // K [-] in the fixed mode
		ExpressionParam concentration_floor_{1e-3};         // c_floor [mol/m3] of the slope and secant modes
		ExpressionParam k_{0.5};                            // lumped k [1/s] (LDF2: kg/(mol s))
		ExpressionParam activation_energy_{0.0};            // lumped E_a [J/mol]; 0: k_L = k
		ExpressionParam pore_radius_{1e-9};                 // Knudsen pore radius [m]

		MassTransferResistances resistances_;               // the stored switches
		int capacity_mode_ = static_cast<int>(CapacityFactorMode::NONE);  // CapacityFactorMode

	private:
		// q* of this component at the state's current mole fractions into out;
		// the state's q* (every component) is left as it was.
		void EvaluateLoading(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component,
		                     int start, int end, const UserVariableRegistry* registry, std::vector<double>& out);
		// The evaluation forms (see the class comment).
		void UpdateLumped(Equations& state, AdsorbentID ads_id, const std::string& component, MixedAdsorbentLayer& layer,
		                  const UserVariableRegistry* registry);
		void UpdateMacropore(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component,
		                     MixedAdsorbentLayer& layer, bool knudsen);
		void UpdateMacroporeSecant(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component,
		                           MixedAdsorbentLayer& layer, const UserVariableRegistry* registry);
		void UpdateSeries(Equations& state, AdsorbentID ads_id, const Fluid& fluid, const std::string& component,
		                  MixedAdsorbentLayer& layer, const UserVariableRegistry* registry);
		// True when the lumped term's activation energy is the number 0.
		bool ActivationEnergyIsZero() const;

		IIsotherm* isotherm_ = nullptr;       // transient, set by AdsorptionSource
		AdsorptionSource* source_ = nullptr;  // transient: under IAST, q* comes from its UpdateIsothermAll
		double solid_density_ = 0.0;          // transient [kg/m3]
		std::vector<std::vector<double>> q_all_scratch_;  // transient
		std::vector<double> y_scratch_, q_base_, q_pert_, capacity_;  // transient

		template <class Archive, class Self>
		static void Fields(Archive& ar, Self& s, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IMassTransferCoefficientModel>(&s),
				cereal::make_nvp("Film_Resistance_Factor", s.film_resistance_factor_),
				cereal::make_nvp("Tortuosity", s.tortuosity_),
				cereal::make_nvp("Macropore_Void_Fraction", s.macropore_void_fraction_),
				cereal::make_nvp("Molecular_Diffusion_Coefficient", s.molecular_diffusion_coefficient_),
				cereal::make_nvp("Crystal_Radius", s.crystal_radius_),
				cereal::make_nvp("Micropore_Diffusivity", s.micropore_diffusivity_),
				cereal::make_nvp("Micropore_Activation_Energy", s.micropore_activation_energy_),
				cereal::make_nvp("Capacity_Factor", s.capacity_factor_),
				cereal::make_nvp("Concentration_Floor", s.concentration_floor_)
			);
			// v2: explicit switches and mode, the lumped and Knudsen parameters.
			if (version >= 2)
				ar(cereal::make_nvp("Resistances", s.resistances_),
				   cereal::make_nvp("Capacity_Factor_Mode", s.capacity_mode_),
				   cereal::make_nvp("Mass_Transfer_Coefficient", s.k_),
				   cereal::make_nvp("Activation_Energy", s.activation_energy_),
				   cereal::make_nvp("Pore_Radius", s.pore_radius_));
		}

		friend class cereal::access;
		template <class Archive>
		void save(Archive& ar, std::uint32_t const version) const { Fields(ar, *this, version); }
		template <class Archive>
		void load(Archive& ar, std::uint32_t const version) {
			Fields(ar, *this, version);
			// v1: the switches and mode follow from the zero sentinels.
			if (version < 2) AdoptVersion1Switches();
		}
	};

	// The generic coefficient model equivalent to `model`: the model itself
	// when it already is one (or null), else a MassTransferCoefficientSeriesResistance
	// with the switches, mode and values of the old class's preset, which
	// gives the same coefficient bit for bit (mass_transfer_legacy.h).
	std::unique_ptr<IMassTransferCoefficientModel> ToResistances(std::unique_ptr<IMassTransferCoefficientModel> model);

	class IMassTransferModel {
	public:
		virtual ~IMassTransferModel() = default;
		virtual MassTransferModelType GetModelType() = 0;

		// Kinetic uptake rate dq/dt [mol/kg/s] for equilibrium loading q_star,
		// current loading q (both [mol/kg]) and coefficient k from mtc_model_.
		// This is the single source of truth for the loading residual
		// (UpdateAdsorbedConcentration); the energy and mass source terms all use
		// the resulting dq/dt, so they stay consistent automatically.
		virtual double Rate(double q_star, double q, double k) const = 0;

		virtual void EvaluateExpressions(const EvalContext& ctx) {
			if (mtc_model_) mtc_model_->EvaluateExpressions(ctx);
		}
		virtual void CompileForRun(const UserVariableRegistry* registry) {
			if (mtc_model_) mtc_model_->CompileForRun(registry);
		}

		// Dual-rate models: the fraction of the equilibrium loading held by the
		// fast site (1 for the single-site models) and the slow site's
		// coefficient model (null for the single-site models).
		virtual double FastSiteFraction() const { return 1.0; }
		virtual IMassTransferCoefficientModel* SlowCoefficientModel() { return nullptr; }
		// True when the model carries a second (slow-site) loading in the state:
		// a dual-rate model whose fast-site fraction is below 1. At a fraction
		// of 1 the slow site holds nothing and the model is the plain LDF, with
		// no extra unknown, so its results are those of the LDF exactly.
		bool HasSlowSite() const { return FastSiteFraction() < 1.0; }

    std::unique_ptr<IMassTransferCoefficientModel> mtc_model_;

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar) {
      ar(
				cereal::make_nvp("Mass_Transfer_Coefficient_Model", mtc_model_)
      );
			// Decks that store an old coefficient class load as the generic model.
			if constexpr (Archive::is_loading::value) mtc_model_ = ToResistances(std::move(mtc_model_));
    }
	};

	class LinearDrivingForce : public IMassTransferModel {
	public:
		LinearDrivingForce();

    MassTransferModelType GetModelType() override { return MassTransferModelType::LDF; }

		// dq/dt = k (q* - q), k in [1/s].
		double Rate(double q_star, double q, double k) const override { return k * (q_star - q); }

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IMassTransferModel>(this)
			);
		}
	};

	class LinearDrivingForceSecondOrder : public IMassTransferModel {
	public:
		LinearDrivingForceSecondOrder();

		MassTransferModelType GetModelType() override { return MassTransferModelType::LDF2; }

		// dq/dt = k (q* - q) |q* - q|, k in [1/(s mol/kg)].
		// Sign-preserving form of the second-order driving force: the literal
		// k (q* - q)^2 is never negative and so could never desorb. The derivative
		// d(rate)/dq = -2 k |q* - q| is continuous (zero at q = q*), so plain abs()
		// needs no smoothing for the Newton iteration.
		double Rate(double q_star, double q, double k) const override {
			const double d = q_star - q;
			return k * d * std::abs(d);
		}

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IMassTransferModel>(this)
			);
		}
	};

	// Dual-rate linear driving force: two parallel LDF sites sharing the
	// component's equilibrium loading q* (its isotherm). The fast site holds
	// the fraction phi of it, the slow site the rest:
	//   dq1/dt = k1 (phi q* - q1),   dq2/dt = k2 ((1 - phi) q* - q2),
	// and the component's loading is q = q1 + q2. The state keeps the total q
	// (what the mass and energy balances, the KPIs and the logs read) and the
	// slow-site loading q2, so the total's rate is
	//   dq/dt = k1 (phi q* - (q - q2)) + k2 ((1 - phi) q* - q2).
	// k1 comes from mtc_model_, k2 from slow_mtc_model_, a second generic
	// coefficient model with its own switches and values. At phi = 1 the model is the first-order LDF with k = k1 and adds
	// no unknown (HasSlowSite() is false). Used for the slow uptake tails of
	// amine sorbents at direct-air-capture concentrations (Stampi-Bombelli et
	// al. 2024, Ind. Eng. Chem. Res. 63, 11637, doi:10.1021/acs.iecr.4c01392).
	class DualRateLinearDrivingForce : public IMassTransferModel {
	public:
		DualRateLinearDrivingForce();

		MassTransferModelType GetModelType() override { return MassTransferModelType::DUAL_RATE_LDF; }

		// The single-site rate (phi = 1): k (q* - q). With a slow site the
		// residual uses FastRate and SlowRate instead.
		double Rate(double q_star, double q, double k) const override { return k * (q_star - q); }
		// Rates of the two sites for total loading q and slow-site loading q2.
		static double FastRate(double phi, double q_star, double q, double q2, double k1) {
			return k1 * (phi * q_star - (q - q2));
		}
		static double SlowRate(double phi, double q_star, double q2, double k2) {
			return k2 * ((1.0 - phi) * q_star - q2);
		}

		double FastSiteFraction() const override { return split_.Get(); }
		IMassTransferCoefficientModel* SlowCoefficientModel() override { return slow_mtc_model_.get(); }

		void EvaluateExpressions(const EvalContext& ctx) override {
			IMassTransferModel::EvaluateExpressions(ctx);
			split_.Evaluate(ctx);
			if (slow_mtc_model_) slow_mtc_model_->EvaluateExpressions(ctx);
		}
		void CompileForRun(const UserVariableRegistry* registry) override {
			IMassTransferModel::CompileForRun(registry);
			if (slow_mtc_model_) slow_mtc_model_->CompileForRun(registry);
		}

		// Gives the dual-rate model a slow-site coefficient model when it has
		// none: a copy of the fast site's (switches, mode and values). An
		// existing slow site is left as it is: each site has its own switches.
		// Returns true when it changed anything.
		bool SyncSlowSite();

		ExpressionParam split_{1.0};  // phi, fast-site fraction of q*, in (0, 1]
		std::unique_ptr<IMassTransferCoefficientModel> slow_mtc_model_;

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const /*version*/) {
			ar(
				cereal::virtual_base_class<IMassTransferModel>(this),
				cereal::make_nvp("Split", split_),
				cereal::make_nvp("Slow_Mass_Transfer_Coefficient_Model", slow_mtc_model_)
			);
			if constexpr (Archive::is_loading::value) slow_mtc_model_ = ToResistances(std::move(slow_mtc_model_));
		}
	};

	// Factory functions
	inline std::unique_ptr<IMassTransferModel> CreateMassTransferModel(MassTransferModelType type) {
		switch (type) {
		case MassTransferModelType::LDF: return std::make_unique<LinearDrivingForce>();
		case MassTransferModelType::LDF2: return std::make_unique<LinearDrivingForceSecondOrder>();
		case MassTransferModelType::DUAL_RATE_LDF: return std::make_unique<DualRateLinearDrivingForce>();
		default: return nullptr;
		}
	}

	// The deck checks of one component's mass-transfer model (both sites of a
	// dual-rate model), each message starting with `where`. Errors (a run
	// would be meaningless): no resistance on; an enabled term with a
	// parameter out of range (tortuosity <= 0, macropore void fraction outside
	// (0, 1], pore radius, crystal radius or micropore diffusivity <= 0,
	// lumped k < 0, film_resistance_factor < 0); the fixed capacity-factor
	// mode with capacity_factor <= 0. A lumped k = 0 is allowed (no uptake).
	// Warnings: film on with film_resistance_factor = 0; a second-order LDF
	// with a physical resistance on (its coefficient is in 1/s, while the
	// LDF2 rate constant is in kg/(mol s)).
	void CheckMassTransferModel(const IMassTransferModel& model, const std::string& where,
	                            std::vector<std::string>& errors, std::vector<std::string>& warnings);

	// A new coefficient model of a kind: always the generic model, with the
	// kind's preset and the old kind's default values (SERIES_RESISTANCE: the
	// version-1 defaults, switches derived from them; see the class comment).
	std::unique_ptr<IMassTransferCoefficientModel> CreateMassTransferCoefficientModel(MassTransferCoefficientType type);

} // end adsim namespace

// Class versions live next to the inline serialize templates: every
// translation unit that saves or loads a deck must see the same numbers.
CEREAL_CLASS_VERSION(adsim::MassTransferResistances, 1);
CEREAL_CLASS_VERSION(adsim::MassTransferCoefficientSeriesResistance, 2);
// Its save/load pair, not the serialize it inherits from the interface.
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(adsim::MassTransferCoefficientSeriesResistance,
                                   cereal::specialization::member_load_save);
CEREAL_CLASS_VERSION(adsim::LinearDrivingForce, 1);
CEREAL_CLASS_VERSION(adsim::LinearDrivingForceSecondOrder, 1);
CEREAL_CLASS_VERSION(adsim::DualRateLinearDrivingForce, 1);

// The load-only old coefficient classes.
#include "adsim/equations/mass_transfer_legacy.h"
