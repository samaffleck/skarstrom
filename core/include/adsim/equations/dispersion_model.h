/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <memory>
#include <string>
#include <vector>

// Cereal includes
#include <cereal/types/memory.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>

#include "adsim/equations/transport_properties.h"
#include "adsim/math/expression_param.h"


namespace adsim {

	// forward declarations
	class PorousMedia;
	class Equations;
	class Adsorbent;
	class MixedAdsorbentLayer;

	static const char* AxialDispersionModeList[] = {
		"Constant",						// 0 
		"Ruthven (1984)",				// 1
		"Taylor-Aris (channels)",		// 2
		"Edwards-Richardson (1968)",	// 3
		"Wakao-Funazkri (1978)"			// 4
	};

	// Stored values; never renumber.
	enum class AxialDispersionType {
		CONSTANT = 0,
		RUTHVEN = 1,
		TAYLOR_ARIS = 2,
		EDWARDS_RICHARDSON = 3,
		WAKAO_FUNAZKRI = 4
	};

	static const char* AxialHeatDispersionModeList[] = {
		"Constant",                 // 0
		"Ruthven (1984)",           // 1
		"Yagi and Kunii (1957)",    // 2
		"Channel conduction"        // 3
	};

	// Stored values; never renumber.
	enum class AxialHeatDispersionType {
		CONSTANT = 0,
		RUTHVEN = 1,
		YAGI_AND_KUNII = 2,
		CHANNEL = 3
	};

	class IDispersionModel {
	public:
		virtual ~IDispersionModel() = default;
		virtual void Update(Equations& state, const MixedAdsorbentLayer& layer) = 0;
		virtual AxialDispersionType GetModelType() = 0;
		virtual unsigned ReadsTransportProperties() const { return kAllTransportProperties; }
		virtual void EvaluateExpressions(const EvalContext& ctx) {}
		// True when the model gives each component its own coefficient
		// D_ax,i (UpdateComponent), built on the component's mixture-averaged
		// molecular diffusivity D_m,i (Equations::molecular_diffusion:
		// Chapman-Enskog binary D_ij with Wilke's equation). The component balance
		// then adds the drift correction that keeps the dispersive molar
		// fluxes summing to zero (Ramshaw, J. Non-Equilib. Thermodyn. 15
		// (1990) 295, in the molar-average frame; UpdateDispersionDrift).
		// Update writes the mole-fraction-weighted mean sum_i y_i D_ax,i, for
		// the log only.
		virtual bool PerComponent() const { return false; }
		// D_ax,i of `component` over the layer's cells, into Dxi. Only called
		// when PerComponent().
		virtual void UpdateComponent(Equations& state, const MixedAdsorbentLayer& layer,
		                             const std::string& component, std::vector<double>& Dxi) {}

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {}
	};

	class IHeatDispersionModel {
	public:
		virtual ~IHeatDispersionModel() = default;
    virtual void Update(Equations& state, const MixedAdsorbentLayer& layer) = 0;
		virtual AxialHeatDispersionType GetModelType() = 0;
		virtual unsigned ReadsTransportProperties() const { return kAllTransportProperties; }
		virtual void EvaluateExpressions(const EvalContext& ctx) {}

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {}
	};

	class DispersionModelConstant : public IDispersionModel {
	public:
    void Update(Equations& state, const MixedAdsorbentLayer& layer) override;
		AxialDispersionType GetModelType() override { return AxialDispersionType::CONSTANT; }
		unsigned ReadsTransportProperties() const override { return kNoTransportProperties; }
		void EvaluateExpressions(const EvalContext& ctx) override { dispersion_coefficient_.Evaluate(ctx); }

		ExpressionParam dispersion_coefficient_{0.0}; // [m2/s]

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IDispersionModel>(this),
				cereal::make_nvp("Constant_Longitudinal_Dispersion", dispersion_coefficient_)
			);
		}
	};

	// Ruthven (1984): D_ax,i = 0.7 D_m,i + 0.5 d_p |u| / eb for each
	// component, with D_m,i its mixture-averaged molecular diffusivity.
	class DispersionModelRuthven : public IDispersionModel {
	public:
		void Update(Equations& state, const MixedAdsorbentLayer& layer) override;
		AxialDispersionType GetModelType() override { return AxialDispersionType::RUTHVEN; } // Ruthven (1984, Principles of adsorption and adsorption processes)
		unsigned ReadsTransportProperties() const override { return kMolecularDiffusion; }
		bool PerComponent() const override { return true; }
		void UpdateComponent(Equations& state, const MixedAdsorbentLayer& layer,
		                     const std::string& component, std::vector<double>& Dxi) override;

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(cereal::virtual_base_class<IDispersionModel>(this));
			LoadLegacyDiffusivity(ar, version);
		}
		// Versions 1 and 2 stored one molecular diffusivity for every
		// component (and version 2 a switch to the per-component values).
		// Both are read and dropped: D_m,i now comes from the gas.
		template <class Archive>
		static void LoadLegacyDiffusivity(Archive& ar, std::uint32_t version) {
			if (version >= 3) return;
			ExpressionParam legacy_dm{0.0};
			ar(cereal::make_nvp("Molecular_Diffusion", legacy_dm));
			if (version == 2) {
				int legacy_choice = 0;
				ar(cereal::make_nvp("Component_Diffusivity", legacy_choice));
			}
		}
		friend class DispersionModelTaylorAris;
	};

	// Taylor-Aris dispersion in laminar channels (structured layers), for each
	// component:
	//   D_ax,i = D_m,i + u_ch^2 D_h^2 / (kappa D_m,i),  u_ch = |u| / OFA,
	// with kappa = 192 for a circular tube (Taylor 1953, Aris 1956; also used
	// for square channels, an approximation) and 840 for parallel plates on
	// D_h = 2 x gap (Aris 1956: 210 on the gap width). D_m,i is the
	// component's mixture-averaged molecular diffusivity. A packed layer
	// cannot use it (deck checks).
	class DispersionModelTaylorAris : public IDispersionModel {
	public:
		void Update(Equations& state, const MixedAdsorbentLayer& layer) override;
		AxialDispersionType GetModelType() override { return AxialDispersionType::TAYLOR_ARIS; }
		unsigned ReadsTransportProperties() const override { return kMolecularDiffusion; }
		bool PerComponent() const override { return true; }
		void UpdateComponent(Equations& state, const MixedAdsorbentLayer& layer,
		                     const std::string& component, std::vector<double>& Dxi) override;

	private:
		static double TaylorAris(int n, double dm, const std::vector<double>& u, const std::vector<double>& eb,
		                         const std::vector<double>& dh, const std::vector<double>& kappa);

		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(cereal::virtual_base_class<IDispersionModel>(this));
			DispersionModelRuthven::LoadLegacyDiffusivity(ar, version);
		}
	};

	// Packed-bed correlations on the particle Peclet number Pe = u_i d_p / D_ax,
	// with u_i = |u| / eb the interstitial velocity and, for each component,
	// Re Sc = |u| d_p / D_m,i (u superficial), so eb / (Re Sc) = D_m,i / (u_i d_p).
	// Both tend to Pe = 2 at high flow, as Ruthven's form. A structured layer
	// cannot use them (deck checks).
	//
	// Edwards and Richardson (1968):
	//   1/Pe = 0.73 eb / (Re Sc) + 0.5 / (1 + 13 * 0.73 eb / (Re Sc)), i.e.
	//   D_ax,i = 0.73 D_m,i + 0.5 u_i d_p / (1 + 9.49 D_m,i / (u_i d_p)).
	// At u -> 0 it gives 0.73 D_m,i (the bed's tortuosity factor on molecular
	// diffusion); the second term blends smoothly to 0.5 u_i d_p.
	class DispersionModelEdwardsRichardson : public IDispersionModel {
	public:
		void Update(Equations& state, const MixedAdsorbentLayer& layer) override;
		AxialDispersionType GetModelType() override { return AxialDispersionType::EDWARDS_RICHARDSON; }
		unsigned ReadsTransportProperties() const override { return kMolecularDiffusion; }
		bool PerComponent() const override { return true; }
		void UpdateComponent(Equations& state, const MixedAdsorbentLayer& layer,
		                     const std::string& component, std::vector<double>& Dxi) override;
		// D_ax for a molecular diffusivity dm, interstitial velocity ui and
		// particle diameter dp.
		static double Coefficient(double dm, double ui, double dp);

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			(void)version;
			ar(cereal::virtual_base_class<IDispersionModel>(this));
		}
	};

	// Wakao and Funazkri (1978): eb D_ax,i / D_m,i = 20 + 0.5 Sc Re, i.e.
	//   D_ax,i = 20 D_m,i / eb + 0.5 u_i d_p.
	// The dispersion with which they re-reduced the packed-bed Sherwood data
	// (the default film correlation's), so the consistent pairing. Its u -> 0
	// limit, 20 D_m,i / eb, is about 70 times the 0.73 D_m,i of Edwards and
	// Richardson at eb = 0.4.
	class DispersionModelWakaoFunazkri : public IDispersionModel {
	public:
		void Update(Equations& state, const MixedAdsorbentLayer& layer) override;
		AxialDispersionType GetModelType() override { return AxialDispersionType::WAKAO_FUNAZKRI; }
		unsigned ReadsTransportProperties() const override { return kMolecularDiffusion; }
		bool PerComponent() const override { return true; }
		void UpdateComponent(Equations& state, const MixedAdsorbentLayer& layer,
		                     const std::string& component, std::vector<double>& Dxi) override;
		static double Coefficient(double dm, double eb, double ui, double dp);

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			(void)version;
			ar(cereal::virtual_base_class<IDispersionModel>(this));
		}
	};

	class HeatDispersionModelConstant : public IHeatDispersionModel {
	public:
    void Update(Equations& state, const MixedAdsorbentLayer& layer) override;
		AxialHeatDispersionType GetModelType() override { return AxialHeatDispersionType::CONSTANT; }
		unsigned ReadsTransportProperties() const override { return kNoTransportProperties; }
		void EvaluateExpressions(const EvalContext& ctx) override { dispersion_coefficient_.Evaluate(ctx); }

		ExpressionParam dispersion_coefficient_{0.0};

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IHeatDispersionModel>(this),
				cereal::make_nvp("Constant_Longitudinal_Heat_Dispersion", dispersion_coefficient_)
			);
		}
	};

	class HeatDispersionModelRuthven : public IHeatDispersionModel {
	public:
    void Update(Equations& state, const MixedAdsorbentLayer& layer) override;
		AxialHeatDispersionType GetModelType() override { return AxialHeatDispersionType::RUTHVEN; } // Ruthven (1984, Principles of adsorption and adsorption processes)
		unsigned ReadsTransportProperties() const override { return kNoTransportProperties; }

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(cereal::virtual_base_class<IHeatDispersionModel>(this));

			if (version >= 2) {
				// Future expansion
			}
		}
	};

	class HeatDispersionModelYagiandKunii : public IHeatDispersionModel {
	public:
    void Update(Equations& state, const MixedAdsorbentLayer& layer) override;
		AxialHeatDispersionType GetModelType() override { return AxialHeatDispersionType::YAGI_AND_KUNII; } // (1957, Studies on effective thermal conductivities in packed beds)
		unsigned ReadsTransportProperties() const override { return kReynoldsNumber | kGasThermalConductivity; }

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(cereal::virtual_base_class<IHeatDispersionModel>(this));

			if (version >= 2) {
				// Future expansion
			}
		}
	};

	// Axial conduction in a structured layer (monolith, laminate), on the
	// cross-sections the gas and the solid occupy in parallel (the parallel
	// bound of the effective conductivity of a two-phase medium; for straight
	// channels with walls continuous along the flow it is exact, see e.g.
	// Kaviany, Principles of Heat Transfer in Porous Media, 2nd ed., 1995):
	//   LTNE: lambda = OFA k_g, the gas's alone; the solid conducts in its own
	//         energy balance (DomainVectors::solid_conductivity, (1 - OFA) k_s).
	//   LTE:  lambda = OFA k_g + (1 - OFA) k_s, one temperature for both.
	// k_s is the layer's effective solid conductivity (the composition-weighted
	// value, MixedAdsorbentLayer::GetEffectiveThermalConductivity). Thermal
	// Taylor dispersion in the channels is neglected, as it is small at the
	// channel Peclet numbers of adsorption contactors. A packed layer cannot use
	// it (deck checks).
	class HeatDispersionModelChannel : public IHeatDispersionModel {
	public:
		void Update(Equations& state, const MixedAdsorbentLayer& layer) override;
		AxialHeatDispersionType GetModelType() override { return AxialHeatDispersionType::CHANNEL; }
		unsigned ReadsTransportProperties() const override { return kGasThermalConductivity; }

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			(void)version;
			ar(cereal::virtual_base_class<IHeatDispersionModel>(this));
		}
	};

	// Factory Functions
	inline std::unique_ptr<IDispersionModel> CreateAxialDispersionModel(AxialDispersionType type) {
		switch (type) {
		case AxialDispersionType::CONSTANT: return std::make_unique<DispersionModelConstant>();
		case AxialDispersionType::RUTHVEN: return std::make_unique<DispersionModelRuthven>();
		case AxialDispersionType::TAYLOR_ARIS: return std::make_unique<DispersionModelTaylorAris>();
		case AxialDispersionType::EDWARDS_RICHARDSON: return std::make_unique<DispersionModelEdwardsRichardson>();
		case AxialDispersionType::WAKAO_FUNAZKRI: return std::make_unique<DispersionModelWakaoFunazkri>();
		default: return nullptr;
		}
	}

	inline std::unique_ptr<IHeatDispersionModel> CreateAxialHeatDispersionModel(AxialHeatDispersionType type) {
		switch (type) {
		case AxialHeatDispersionType::CONSTANT: return std::make_unique<HeatDispersionModelConstant>();
		case AxialHeatDispersionType::RUTHVEN: return std::make_unique<HeatDispersionModelRuthven>();
		case AxialHeatDispersionType::YAGI_AND_KUNII: return std::make_unique<HeatDispersionModelYagiandKunii>();
		case AxialHeatDispersionType::CHANNEL: return std::make_unique<HeatDispersionModelChannel>();
		default: return nullptr;
		}
	}

} // end adsim namespace

// Class versions live next to the inline serialize templates: every
// translation unit that saves or loads a deck must see the same numbers.
CEREAL_CLASS_VERSION(adsim::DispersionModelConstant, 1);
CEREAL_CLASS_VERSION(adsim::HeatDispersionModelConstant, 1);
CEREAL_CLASS_VERSION(adsim::DispersionModelRuthven, 3);
CEREAL_CLASS_VERSION(adsim::DispersionModelTaylorAris, 3);
CEREAL_CLASS_VERSION(adsim::DispersionModelEdwardsRichardson, 1);
CEREAL_CLASS_VERSION(adsim::DispersionModelWakaoFunazkri, 1);
CEREAL_CLASS_VERSION(adsim::HeatDispersionModelRuthven, 1);
CEREAL_CLASS_VERSION(adsim::HeatDispersionModelYagiandKunii, 1);
CEREAL_CLASS_VERSION(adsim::HeatDispersionModelChannel, 1);
