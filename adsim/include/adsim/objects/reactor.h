/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <string>

// adsim includes
#include "adsim/objects/reactor_equations.h"
#include "adsim/objects/porous_media.h"
#include "adsim/objects/bed_shape.h"
#include "adsim/objects/wall.h"
#include "adsim/objects/void_region.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/solver/initial_conditions.h"
#include "adsim/math/expression_param.h"

// Cereal includes
#include <cereal/types/string.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/variant.hpp>
#include <cereal/access.hpp>


namespace adsim {

	// Forward declerations
	struct Fluid;

	class Bed {
	public:
		void Initialise(const AdsorbentDirectory& dir, Fluid& fluid);
		// Also verifies every layer composition id exists in `dir`.
		bool PerformChecks(const AdsorbentDirectory& dir, std::string& error_msg) const;
		int GetNumberOfCells() const;
		void RemoveComponent(const std::string& component, const AdsorbentDirectory& dir);

		double GetArea() const { return adsim::GetArea(shape_); }
		double GetHydraulicDiameter() const { return adsim::GetHydraulicDiameter(shape_); }

		// Cell-index helpers. Indices refer to positions in the global cell
		// array (size = GetNumberOfCells()). Returns -1 when the corresponding
		// void is disabled.
		int GetBottomVoidCellIndex() const { return bottom_void_.enabled_ ? 1 : -1; }
		int GetTopVoidCellIndex() const {
			return top_void_.enabled_ ? GetNumberOfCells() - 2 : -1;
		}
		// First and last bed (porous-media) cell indices. Excludes ghosts and
		// voids.
		int GetFirstBedCellIndex() const { return 1 + (bottom_void_.enabled_ ? 1 : 0); }
		int GetLastBedCellIndex() const {
			return GetNumberOfCells() - 2 - (top_void_.enabled_ ? 1 : 0);
		}

		// Returns the VoidRegion at the given global cell index, or nullptr if
		// the cell is not a void cell.
		const VoidRegion* VoidRegionAt(int cell_idx) const {
			if (bottom_void_.enabled_ && cell_idx == GetBottomVoidCellIndex()) return &bottom_void_;
			if (top_void_.enabled_    && cell_idx == GetTopVoidCellIndex())    return &top_void_;
			return nullptr;
		}

		BedShape shape_ = CylinderBed{ExpressionParam{0.1}};
		ReactorEquations equations_{};
		BedWall wall_{};
		AdsorbentLayers layers_;
		ExpressionParam ambient_temperature_{293.0};
		VoidRegion bottom_void_{};
		VoidRegion top_void_{};

		void EvaluateExpressions(const EvalContext& ctx) {
			ambient_temperature_.Evaluate(ctx);
			EvaluateBedShapeExpressions(shape_, ctx);
			wall_.EvaluateExpressions(ctx);
			equations_.EvaluateExpressions(ctx);
			bottom_void_.EvaluateExpressions(ctx);
			top_void_.EvaluateExpressions(ctx);
			for (auto& layer : layers_.layers_) {
				layer.EvaluateExpressions(ctx);
			}
		}

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::make_nvp("Shape", shape_),
				cereal::make_nvp("Wall", wall_),
				cereal::make_nvp("Layers", layers_),
				cereal::make_nvp("Equations", equations_),
				cereal::make_nvp("Ambient_Temperature", ambient_temperature_)
			);
			if (version >= 4) {
				ar(
					cereal::make_nvp("Bottom_Void", bottom_void_),
					cereal::make_nvp("Top_Void", top_void_)
				);
			}
		}
	};

  // Free functions
  inline double GetAdsorbentMass(const MixedAdsorbentLayer& layer, const AdsorbentDirectory& dir, double area) {
    const double length = layer.domain_.length;
    const double volume = length * area;
    const double bulk_density = layer.GetBulkDensity(dir);
    const double mass = volume * bulk_density;
    return mass;
  }

  inline double GetAdsorbentMass(const Bed& bed, const AdsorbentDirectory& dir) {
    double mass = 0.0;
    const double area = bed.GetArea();

    for (const auto& layer : bed.layers_.layers_) {
      mass += GetAdsorbentMass(layer, dir, area);
    }

    return mass;
  }

} // end adsim namespace
