/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Cereal includes
#include <cereal/access.hpp>

// STL includes
#include <cassert>
#include <map>

// Adsim includes
#include "adsim/math/constants.h"
#include "adsim/math/expression_param.h"


namespace adsim {

  const double PI = adsim::constants::kPI;
  class Bed;
  class Fluid;
  class AdsorbentDirectory;
  class EnergySource;
  using AdsorbentID = std::size_t; // mirrors porous_media.h; identical typedef redeclaration is legal

  /*
  Domain for a 1D cylindrical geometry
  Stores minimal information needed to construct domain
  Diameter is now a Bed-level property (shared across all layers)
  */
  class Domain {
  public:
    int n_cells = 30;
    ExpressionParam length{1.0};
    int start_index{};
    int end_index{};

    double GetCellWidth() const {
      assert(n_cells > 0);
      assert(length.Get() > 0.0);
      return length.Get() / n_cells;
    }

    double GetCellVolume(double area) const {
      return area * GetCellWidth();
    }

    void EvaluateExpressions(const EvalContext& ctx) {
      length.Evaluate(ctx);
    }

  private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::make_nvp("Number_Of_Cells", n_cells),
        cereal::make_nvp("Length", length)
      );
    }
  };


  struct DomainVectors {

    void Initialise(const Bed& bed, const AdsorbentDirectory& dir, const Fluid& fluid);
    void PopulateEnergySources(const std::map<size_t, EnergySource>& energy_sources, const Bed& bed);

    std::vector<double> positions;   // dx[i] - cell width [m]
    std::vector<double> cell_width;   // dx[i] - cell width [m]
    std::vector<double> area;         // A[i] - cross-sectional area [m2]
    std::vector<double> volume;       // V[i] - cell volume [m3]
    std::vector<double> permeability; // viscous permeability K [m2] (A = mu / K)
    std::vector<double> inertial_coefficient; // Ergun inertial beta [1/m] (B = beta * rho); 0 for Darcy
    std::vector<double> Cps;
    std::vector<double> bulk_density;
    std::vector<double> et;
    std::vector<double> eb;
    std::vector<double> ep;
    std::vector<double> dp;
    std::vector<double> ks; // Thermal conductivity of the adsorbent [W/(m-K)]
    std::vector<double> kw; // Thermal conductivity of the wall [W/(m-K)]
    std::vector<double> Cpw; // Specific heat capacity of the wall
    std::vector<double> rho_w; // Density of the wall
    std::vector<double> di; // Internal diameter wall
    std::vector<double> de; // External diameter wall
    std::unordered_map<std::string, std::vector<double>> dH; // Heat of adsorption [J/mol], mass-fraction-weighted aggregate per component
    std::vector<double> energy_source_power; // External energy source [W/m^3]

    // Per-adsorbent decomposition of the packing, used by the source terms so
    // each adsorbent contributes its own (rho_b,k * dq_k/dt) without an
    // effective isotherm. For a single-adsorbent cell, the single entry's
    // bulk_density equals `bulk_density[i]` and dH_ads equals `dH`.
    struct CellAdsorbent {
      AdsorbentID id;
      double mass_fraction;  // w_k in this cell's layer [-]
      double bulk_density;   // (1 - eb) * rho_mix * w_k  [kg adsorbent k / m3 bed]
    };
    std::vector<std::vector<CellAdsorbent>> cell_adsorbents; // per cell; empty in void/ghost cells
    // Heat of adsorption per adsorbent per component, per cell [J/mol].
    std::unordered_map<AdsorbentID, std::unordered_map<std::string, std::vector<double>>> dH_ads;

    // True at cells that represent top/bottom void space (no adsorbent).
    // PorousMediaOnly equations skip these cells; void cells receive lumped
    // 0D well-mixed treatment via the standard interior-cell residuals with
    // adsorbent-related coefficients zeroed.
    std::vector<bool> is_void_cell;

    // Plotting helpers — populated alongside `positions` in Initialise().
    // `bed_cell_indices` lists every cell index that is part of the packed
    // bed (no ghosts, no voids), in axial order.
    //
    // `bed_positions` describes the axial x-axis the plotter should use when
    // showing bed-spatial profiles: cumulative widths inside the packed
    // material, with origin at the bed inlet face (position 0 = entry of
    // porous material). It has the same length as the bed-filtered face
    // value arrays produced by `bed_face_index_begin` .. `bed_face_index_end`
    // (a contiguous subrange of `positions`). When no voids are enabled this
    // matches `positions` exactly (with the origin shift being a no-op since
    // there are no preceding void widths).
    //
    // `bed_face_index_begin` and `bed_face_index_end` give a contiguous slice
    // [begin, end) into the global face-value arrays (e.g. `positions` and
    // `VariableLogEntry::values[t]`) covering only the bed-internal and
    // bed-boundary faces. Used by the plotter to subset variable_data_list.
    std::vector<int> bed_cell_indices;
    std::vector<double> bed_positions;
    int bed_face_index_begin = 0;
    int bed_face_index_end = 0;

    // Global cell index of each enabled void cell, or -1 if disabled.
    // Used by the plotter to extract the void's single-cell time series.
    int bottom_void_cell_index = -1;
    int top_void_cell_index = -1;

    bool IsVoidCell(size_t i) const {
      return i < is_void_cell.size() && is_void_cell[i];
    }
  };

} // end adsim namespace
