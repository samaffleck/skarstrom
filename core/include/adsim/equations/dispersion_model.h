/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <memory>

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
		"Taylor-Aris (channels)"		// 2
	};

	// Stored values; never renumber.
	enum class AxialDispersionType {
		CONSTANT = 0,
		RUTHVEN = 1,
		TAYLOR_ARIS = 2
	};

	static const char* AxialHeatDispersionModeList[] = {
		"Constant", 
		"Ruthven (1984)",	
		"Yagi and Kunii (1957)"	
	};
	
	enum class AxialHeatDispersionType {
		CONSTANT,
		RUTHVEN,
		YAGI_AND_KUNII
	};

	class IDispersionModel {
	public:
		virtual ~IDispersionModel() = default;
		virtual void Update(Equations& state, const MixedAdsorbentLayer& layer) = 0;
		virtual AxialDispersionType GetModelType() = 0;
		virtual unsigned ReadsTransportProperties() const { return kAllTransportProperties; }
		virtual void EvaluateExpressions(const EvalContext& ctx) {}

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

	class DispersionModelRuthven : public IDispersionModel {
	public:
    void Update(Equations& state, const MixedAdsorbentLayer& layer) override;
    AxialDispersionType GetModelType() override { return AxialDispersionType::RUTHVEN; } // Ruthven (1984, Principles of adsorption and adsorption processes)
    unsigned ReadsTransportProperties() const override { return kNoTransportProperties; }
		void EvaluateExpressions(const EvalContext& ctx) override { Dm.Evaluate(ctx); }

    ExpressionParam Dm{1.6e-5}; // [m2/s]

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
        cereal::virtual_base_class<IDispersionModel>(this),
        cereal::make_nvp("Molecular_Diffusion", Dm)
      );
		}
	};

	// Taylor-Aris dispersion in laminar channels (structured layers):
	//   D_ax = D_m + u_ch^2 D_h^2 / (kappa D_m),  u_ch = |u| / OFA,
	// with kappa = 192 for a circular tube (Taylor 1953, Aris 1956; also used
	// for square channels, an approximation) and 840 for parallel plates on
	// D_h = 2 x gap (Aris 1956: 210 on the gap width). D_m is this model's
	// molecular diffusion coefficient. A packed layer cannot use it (deck
	// checks).
	class DispersionModelTaylorAris : public IDispersionModel {
	public:
		void Update(Equations& state, const MixedAdsorbentLayer& layer) override;
		AxialDispersionType GetModelType() override { return AxialDispersionType::TAYLOR_ARIS; }
		unsigned ReadsTransportProperties() const override { return kNoTransportProperties; }
		void EvaluateExpressions(const EvalContext& ctx) override { Dm.Evaluate(ctx); }

		ExpressionParam Dm{1.6e-5}; // [m2/s]

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			(void)version;
			ar(
				cereal::virtual_base_class<IDispersionModel>(this),
				cereal::make_nvp("Molecular_Diffusion", Dm)
			);
		}
	};

	// class DispersionModelWakaoAndFunazkri : public IDispersionModel {
	// public:
  //   void Update(Equations& state, const MixedAdsorbentLayer& layer) override;
  //   void Update(Equations& state, const std::string& component, const MixedAdsorbentLayer& layer) override;
	// 	AxialDispersionType GetModelType() override { return AxialDispersionType::WAKO_AND_FUNAZKRI; } 

	// private:
	// 	friend class cereal::access;
	// 	template <class Archive>
	// 	void serialize(Archive& ar, std::uint32_t const version) {
	// 		ar(cereal::virtual_base_class<IDispersionModel>(this));

	// 		if (version >= 2) {
	// 			// Future expansion
	// 		}
	// 	}
	// };

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

	// Factory Functions
	inline std::unique_ptr<IDispersionModel> CreateAxialDispersionModel(AxialDispersionType type) {
		switch (type) {
		case AxialDispersionType::CONSTANT: return std::make_unique<DispersionModelConstant>();
		case AxialDispersionType::RUTHVEN: return std::make_unique<DispersionModelRuthven>();
		case AxialDispersionType::TAYLOR_ARIS: return std::make_unique<DispersionModelTaylorAris>();
		// case AxialDispersionType::WAKO_AND_FUNAZKRI: return std::make_unique<DispersionModelWakaoAndFunazkri>();
		default: return nullptr;
		}
	}

	inline std::unique_ptr<IHeatDispersionModel> CreateAxialHeatDispersionModel(AxialHeatDispersionType type) {
		switch (type) {
		case AxialHeatDispersionType::CONSTANT: return std::make_unique<HeatDispersionModelConstant>();
		case AxialHeatDispersionType::RUTHVEN: return std::make_unique<HeatDispersionModelRuthven>();
		case AxialHeatDispersionType::YAGI_AND_KUNII: return std::make_unique<HeatDispersionModelYagiandKunii>();
		default: return nullptr;
		}
	}

} // end adsim namespace

// Class versions live next to the inline serialize templates: every
// translation unit that saves or loads a deck must see the same numbers.
CEREAL_CLASS_VERSION(adsim::DispersionModelConstant, 1);
CEREAL_CLASS_VERSION(adsim::HeatDispersionModelConstant, 1);
CEREAL_CLASS_VERSION(adsim::DispersionModelRuthven, 1);
CEREAL_CLASS_VERSION(adsim::DispersionModelTaylorAris, 1);
CEREAL_CLASS_VERSION(adsim::HeatDispersionModelRuthven, 1);
CEREAL_CLASS_VERSION(adsim::HeatDispersionModelYagiandKunii, 1);
