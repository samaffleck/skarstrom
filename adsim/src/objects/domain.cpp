/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/objects/domain.h"
#include "adsim/objects/reactor.h"
#include "adsim/objects/fluid.h"
#include "adsim/objects/reactor_equations.h"


namespace adsim {

  void DomainVectors::Initialise(const Bed& bed, const AdsorbentDirectory& dir, const Fluid& fluid) {
    int N = bed.GetNumberOfCells();

    for (const auto& component : fluid.components_) {
      dH[component].resize(N);
    }

    cell_width.resize(N);
    area.resize(N);
    volume.resize(N);
    permeability.resize(N);
    inertial_coefficient.resize(N);
    Cps.resize(N);
    bulk_density.resize(N);
    et.resize(N);
    eb.resize(N);
    ep.resize(N);
    dp.resize(N);
    ks.resize(N);
    kw.resize(N);
    Cpw.resize(N);
    rho_w.resize(N);
    di.resize(N);
    de.resize(N);
    energy_source_power.resize(N, 0.0);
    is_void_cell.assign(N, false);

    // Per-adsorbent packing decomposition (filled in the layer loop below).
    cell_adsorbents.assign(N, {});
    dH_ads.clear();
    for (const auto& layer : bed.layers_.layers_) {
      for (const auto& [ads_id, frac] : layer.composition_) {
        auto& m = dH_ads[ads_id];
        for (const auto& component : fluid.components_) {
          if (m[component].empty()) m[component].assign(N, 0.0);
        }
      }
    }

    double c_kw = bed.wall_.thermal_conductivity_;
    double c_cpw = bed.wall_.heat_capacity_;
    double c_rhow = bed.wall_.density_;

    for (const auto& layer : bed.layers_.layers_) {
      double c_width = layer.domain_.GetCellWidth();
      double c_area = bed.GetArea();
      double c_volume = layer.domain_.GetCellVolume(c_area);
      double c_perm = layer.GetPermeability();
      double c_beta = layer.GetInertialCoefficient();
      double c_heat_cap = layer.GetEffectiveHeatCapacity(dir);
      double c_bulk_density = layer.GetBulkDensity(dir);
      double c_ep = layer.GetEffectiveParticleVoidFraction(dir);
      double c_et = layer.GetTotalVoidFraction(dir);
      double c_eb = layer.eb_;
      double c_dp = layer.particle_diameter_;
      double c_ks = layer.GetEffectiveThermalConductivity(dir);
      const auto c_ads_bulk_density = layer.GetAdsorbentBulkDensities(dir);
      double c_di = bed.GetHydraulicDiameter();
      double c_de = c_di + 2 * bed.wall_.thickness_;

      for (int i = layer.domain_.start_index; i <= layer.domain_.end_index; i++) {
        cell_width[i] = c_width;
        area[i] = c_area;
        volume[i] = c_volume;
        permeability[i] = c_perm;
        inertial_coefficient[i] = c_beta;
        Cps[i] = c_heat_cap;
        bulk_density[i] = c_bulk_density;
        et[i] = c_et;
        eb[i] = c_eb;
        ep[i] = c_ep;
        dp[i] = c_dp;
        ks[i] = c_ks;
        kw[i] = c_kw;
        Cpw[i] = c_cpw;
        rho_w[i] = c_rhow;
        di[i] = c_di;
        de[i] = c_de;

        // Per-adsorbent packing decomposition: bulk density of each adsorbent k
        // in this cell is (1 - eb) * rho_mix * w_k; these sum to bulk_density[i].
        cell_adsorbents[i].clear();
        cell_adsorbents[i].reserve(layer.composition_.size());
        for (const auto& [ads_id, frac] : layer.composition_) {
          cell_adsorbents[i].push_back({ads_id, frac, c_ads_bulk_density.at(ads_id)});
        }

        // Heat of adsorption: per-adsorbent values plus a mass-fraction-weighted
        // aggregate per component (the aggregate is kept for legacy/logging use).
        for (const auto& component : fluid.components_) {
          double weighted_dH = 0.0;
          for (const auto& [ads_id, frac] : layer.composition_) {
            const auto& ads = dir.Get(ads_id);
            auto it = ads.adsorption_source_.heat_of_adsorption_.find(component);
            const double hoa_k = (it != ads.adsorption_source_.heat_of_adsorption_.end())
                                 ? static_cast<double>(it->second) : 0.0;
            weighted_dH += frac * hoa_k;
            dH_ads[ads_id][component][i] = hoa_k;
          }
          dH[component][i] = weighted_dH;
        }
      }
    }

    auto copy_cell_from = [&](int dst, int src) {
      cell_width[dst] = cell_width[src];
      area[dst] = area[src];
      volume[dst] = volume[src];
      permeability[dst] = permeability[src];
      inertial_coefficient[dst] = inertial_coefficient[src];
      Cps[dst] = Cps[src];
      bulk_density[dst] = bulk_density[src];
      et[dst] = et[src];
      eb[dst] = eb[src];
      ep[dst] = ep[src];
      dp[dst] = dp[src];
      ks[dst] = ks[src];
      kw[dst] = kw[src];
      Cpw[dst] = Cpw[src];
      rho_w[dst] = rho_w[src];
      di[dst] = di[src];
      de[dst] = de[src];
      energy_source_power[dst] = energy_source_power[src];
      for (const auto& component : fluid.components_) {
        dH[component][dst] = dH[component][src];
      }
    };

    // Populate void cells (one per end if enabled).
    const double c_area_void = bed.GetArea();
    const double c_di = bed.GetHydraulicDiameter();
    const double c_de = c_di + 2 * bed.wall_.thickness_;

    auto populate_void = [&](int idx, const VoidRegion& vr) {
      const double V_void = vr.volume_.Get();
      const double dx_void = (c_area_void > 0.0) ? V_void / c_area_void : 0.0;
      cell_width[idx] = dx_void;
      area[idx] = c_area_void;
      volume[idx] = V_void;
      // Void cells have no physical permeability — Darcy's law does not apply
      // inside a lumped tank. UpdateVelocity special-cases void-touching faces
      // to use the bed-side cell's resistance only, so this value is never
      // read for momentum.
      permeability[idx] = 0.0;
      inertial_coefficient[idx] = 0.0;
      Cps[idx] = 0.0;
      bulk_density[idx] = 0.0;
      et[idx] = 1.0;
      eb[idx] = 1.0;
      ep[idx] = 0.0;
      dp[idx] = dx_void;  // placeholder; not used by physics in void cells
      ks[idx] = 0.0;
      kw[idx] = c_kw;
      Cpw[idx] = c_cpw;
      rho_w[idx] = c_rhow;
      di[idx] = c_di;
      de[idx] = c_de;
      energy_source_power[idx] = 0.0;
      for (const auto& component : fluid.components_) {
        dH[component][idx] = 0.0;
      }
      is_void_cell[idx] = true;
    };

    const int bottom_void_idx = bed.GetBottomVoidCellIndex();
    const int top_void_idx = bed.GetTopVoidCellIndex();
    bottom_void_cell_index = bottom_void_idx;
    top_void_cell_index = top_void_idx;
    if (bottom_void_idx > 0) populate_void(bottom_void_idx, bed.bottom_void_);
    if (top_void_idx > 0)    populate_void(top_void_idx, bed.top_void_);

    // Ghost cells mirror the first/last *bed* (porous-media) cell so that BC
    // residuals (central-difference across the ghost⇄interior face) reference
    // representative bed properties, not the void cell next door.
    const int first_bed = bed.GetFirstBedCellIndex();
    const int last_bed = bed.GetLastBedCellIndex();
    copy_cell_from(0, first_bed);
    copy_cell_from(N - 1, last_bed);

    // Face positions: cumulative axial coordinate at each cell-cell boundary.
    // positions[0] = 0 sits at the left face of cell 1; subsequent entries
    // add each interior cell's width in turn. For uniform cell widths this
    // reduces to the legacy `i * dx` form; for non-uniform widths (e.g.
    // void cells with different dx_void than the bed) it correctly tracks
    // actual axial position.
    positions.assign(N - 1, 0.0);
    for (int i = 1; i < N - 1; i++) {
      positions[i] = positions[i - 1] + cell_width[i];
    }

    // Bed-only plotting helpers. Cell indices in [first_bed, last_bed] are
    // the packed-bed cells; positions follow the global face-index convention
    // (positions[i] = i * cell_width[i+1]) shifted so face (first_bed - 1)
    // sits at axial 0. The face index slice [bed_face_index_begin,
    // bed_face_index_end) covers exactly the bed inlet face through the bed
    // outlet face, so values from variable_data_list can be subset by the
    // same range.
    bed_cell_indices.clear();
    bed_positions.clear();
    bed_face_index_begin = 0;
    bed_face_index_end = 0;
    if (first_bed >= 0 && last_bed >= first_bed) {
      bed_cell_indices.reserve(last_bed - first_bed + 1);
      for (int i = first_bed; i <= last_bed; ++i) {
        bed_cell_indices.push_back(i);
      }

      // Face indices of interest: faces between bed cells, plus the bed
      // inlet and outlet boundary faces. Face i sits between cell i and
      // cell i+1, so the bed-inlet face has global index (first_bed - 1)
      // and the bed-outlet face has global index last_bed.
      bed_face_index_begin = std::max(0, first_bed - 1);
      bed_face_index_end = std::min(static_cast<int>(positions.size()), last_bed + 1);

      const int n = bed_face_index_end - bed_face_index_begin;
      if (n > 0) {
        bed_positions.reserve(n);
        const double origin = positions[bed_face_index_begin];
        for (int k = 0; k < n; ++k) {
          bed_positions.push_back(positions[bed_face_index_begin + k] - origin);
        }
      }
    }
  }

  void DomainVectors::PopulateEnergySources(
      const std::map<size_t, EnergySource>& energy_sources, const Bed& bed) {
    std::fill(energy_source_power.begin(), energy_source_power.end(), 0.0);

    for (const auto& layer : bed.layers_.layers_) {
      auto it = energy_sources.find(layer.id_);
      if (it == energy_sources.end()) continue;

      const auto& src = it->second;
      if (!src.enabled_ || src.power_ == 0.0) continue;

      double cell_volume = layer.domain_.GetCellVolume(bed.GetArea());
      int n_cells = layer.domain_.n_cells;
      double V_layer = cell_volume * n_cells;
      double q_vol = src.power_ / V_layer;

      for (int i = layer.domain_.start_index; i <= layer.domain_.end_index; ++i) {
        energy_source_power[i] = q_vol;
      }
    }
  }

} // end adsim namespace
