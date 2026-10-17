// Unit tests for top/bottom void space wired into the Bed simulation.
//
// Covers the building blocks (no solver integration): cell counting, cell
// layout, DomainVectors population, ghost mirroring across voids, and the
// EquationDomain abstraction.

#include "../fixtures/user_data_test.h"
#include "adsim/objects/void_region.h"
#include "adsim/equations/isotherm.h"

namespace adsim {

  class VoidRegionTest : public UserDataTest {
  protected:
    // Build a minimal but solver-ready bed setup so DomainVectors::Initialise
    // can run. One layer, one component, simple Henry isotherm.
    void BuildMinimalBed(int n_cells = 10) {
      AddComponents({"Nitrogen"});
      auto [layer_id, ads_id] = AddLayer("Layer1");

      auto& layer = GetLayer(layer_id);
      layer.domain_.length = 1.0;
      layer.domain_.n_cells = n_cells;

      auto& ads = GetAdsorbent(ads_id);
      ads.solid_density_ = 1000.0;
      ads.heat_capacity_ = 900.0;
      ads.thermal_conductivity_ = 0.5;
      ads.ep_ = 0.35;

      auto henry = std::make_unique<HenryIsotherm>("Nitrogen");
      henry->params_["Nitrogen"].k = 0.0;  // inert (no adsorption)
      henry->params_["Nitrogen"].Q = 0.0;
      AddIsotherm({"Nitrogen"}, ads_id, std::move(henry));

      AddLinearDrivingForceModel(ads_id, "Nitrogen", 1.0);
      SetHOA(ads_id, "Nitrogen", 0.0);

      AddConstantInternalHeatTransferCoefficientModel(10.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(10.0, layer_id);

      std::get<CylinderBed>(data.bed.shape_).diameter = 0.05;
      AddInitialConditions(1.0e5, 298.0, 298.0, {{"Nitrogen", 1.0}});
    }
  };

  TEST_F(VoidRegionTest, DefaultsAreDisabled) {
    VoidRegion vr;
    EXPECT_FALSE(vr.enabled_);
    EXPECT_GT(vr.volume_.Get(), 0.0);
  }

  TEST_F(VoidRegionTest, BedNumberOfCells_NoVoids) {
    BuildMinimalBed(/*n_cells=*/10);
    EXPECT_EQ(data.bed.GetNumberOfCells(), 12);  // 10 layer + 2 ghosts
    EXPECT_EQ(data.bed.GetFirstBedCellIndex(), 1);
    EXPECT_EQ(data.bed.GetLastBedCellIndex(), 10);
    EXPECT_EQ(data.bed.GetBottomVoidCellIndex(), -1);
    EXPECT_EQ(data.bed.GetTopVoidCellIndex(), -1);
  }

  TEST_F(VoidRegionTest, BedNumberOfCells_BottomOnly) {
    BuildMinimalBed(10);
    data.bed.bottom_void_.enabled_ = true;
    EXPECT_EQ(data.bed.GetNumberOfCells(), 13);
    EXPECT_EQ(data.bed.GetBottomVoidCellIndex(), 1);
    EXPECT_EQ(data.bed.GetFirstBedCellIndex(), 2);
    EXPECT_EQ(data.bed.GetLastBedCellIndex(), 11);
    EXPECT_EQ(data.bed.GetTopVoidCellIndex(), -1);
  }

  TEST_F(VoidRegionTest, BedNumberOfCells_TopOnly) {
    BuildMinimalBed(10);
    data.bed.top_void_.enabled_ = true;
    EXPECT_EQ(data.bed.GetNumberOfCells(), 13);
    EXPECT_EQ(data.bed.GetBottomVoidCellIndex(), -1);
    EXPECT_EQ(data.bed.GetFirstBedCellIndex(), 1);
    EXPECT_EQ(data.bed.GetLastBedCellIndex(), 10);
    EXPECT_EQ(data.bed.GetTopVoidCellIndex(), 11);
  }

  TEST_F(VoidRegionTest, BedNumberOfCells_BothVoids) {
    BuildMinimalBed(10);
    data.bed.bottom_void_.enabled_ = true;
    data.bed.top_void_.enabled_ = true;
    EXPECT_EQ(data.bed.GetNumberOfCells(), 14);
    EXPECT_EQ(data.bed.GetBottomVoidCellIndex(), 1);
    EXPECT_EQ(data.bed.GetFirstBedCellIndex(), 2);
    EXPECT_EQ(data.bed.GetLastBedCellIndex(), 11);
    EXPECT_EQ(data.bed.GetTopVoidCellIndex(), 12);
  }

  TEST_F(VoidRegionTest, DomainVectors_VoidCellPropertiesPopulated) {
    BuildMinimalBed(10);
    data.bed.bottom_void_.enabled_ = true;
    data.bed.bottom_void_.volume_ = 5e-5;
    data.bed.top_void_.enabled_ = true;
    data.bed.top_void_.volume_ = 2e-5;

    data.Initialise();

    const auto& dv = data.state.domain_vectors;
    const int N = data.bed.GetNumberOfCells();
    ASSERT_EQ((int)dv.is_void_cell.size(), N);

    const int iv_bot = data.bed.GetBottomVoidCellIndex();
    const int iv_top = data.bed.GetTopVoidCellIndex();
    ASSERT_GT(iv_bot, 0);
    ASSERT_GT(iv_top, 0);
    EXPECT_TRUE(dv.IsVoidCell(iv_bot));
    EXPECT_TRUE(dv.IsVoidCell(iv_top));

    const double area = data.bed.GetArea();
    EXPECT_DOUBLE_EQ(dv.eb[iv_bot], 1.0);
    EXPECT_DOUBLE_EQ(dv.ep[iv_bot], 0.0);
    EXPECT_DOUBLE_EQ(dv.et[iv_bot], 1.0);
    EXPECT_DOUBLE_EQ(dv.bulk_density[iv_bot], 0.0);
    EXPECT_DOUBLE_EQ(dv.Cps[iv_bot], 0.0);
    EXPECT_DOUBLE_EQ(dv.ks[iv_bot], 0.0);
    EXPECT_DOUBLE_EQ(dv.volume[iv_bot], 5e-5);
    EXPECT_NEAR(dv.cell_width[iv_bot], 5e-5 / area, 1e-15);
    EXPECT_DOUBLE_EQ(dv.volume[iv_top], 2e-5);
    // Wall props propagate so the wall energy balance remains continuous.
    EXPECT_DOUBLE_EQ(dv.Cpw[iv_bot], data.bed.wall_.heat_capacity_);
    EXPECT_DOUBLE_EQ(dv.rho_w[iv_bot], data.bed.wall_.density_);
    EXPECT_DOUBLE_EQ(dv.kw[iv_bot], data.bed.wall_.thermal_conductivity_);
  }

  TEST_F(VoidRegionTest, DomainVectors_GhostsMirrorBedNotVoid) {
    BuildMinimalBed(10);
    data.bed.bottom_void_.enabled_ = true;
    data.bed.bottom_void_.volume_ = 1e-4;
    data.bed.top_void_.enabled_ = true;
    data.bed.top_void_.volume_ = 1e-4;

    data.Initialise();

    const auto& dv = data.state.domain_vectors;
    const int first_bed = data.bed.GetFirstBedCellIndex();
    const int last_bed = data.bed.GetLastBedCellIndex();
    const int N = (int)dv.eb.size();

    // Ghost bed-related properties should match the first/last *bed* cell,
    // not the adjacent void cell (which has eb=1 and bulk_density=0).
    EXPECT_DOUBLE_EQ(dv.eb[0], dv.eb[first_bed]);
    EXPECT_DOUBLE_EQ(dv.bulk_density[0], dv.bulk_density[first_bed]);
    EXPECT_DOUBLE_EQ(dv.Cps[0], dv.Cps[first_bed]);
    EXPECT_DOUBLE_EQ(dv.eb[N - 1], dv.eb[last_bed]);
    EXPECT_DOUBLE_EQ(dv.bulk_density[N - 1], dv.bulk_density[last_bed]);

    // Sanity: the void cells are different from the ghosts.
    EXPECT_NE(dv.eb[0], dv.eb[1]);  // ghost (eb=bed) vs bottom void (eb=1)
  }

  TEST_F(VoidRegionTest, IsVoidCell_OnlyMarkedAtVoidIndices) {
    BuildMinimalBed(10);
    data.bed.bottom_void_.enabled_ = true;

    data.Initialise();

    const auto& dv = data.state.domain_vectors;
    const int N = data.bed.GetNumberOfCells();
    int void_count = 0;
    for (int i = 0; i < N; ++i) {
      if (dv.is_void_cell[i]) ++void_count;
    }
    EXPECT_EQ(void_count, 1);
    EXPECT_TRUE(dv.is_void_cell[data.bed.GetBottomVoidCellIndex()]);
  }

  // EquationDomain abstraction tests
  TEST_F(VoidRegionTest, EquationDomain_PorousOnlyTagged) {
    BuildMinimalBed(10);
    data.Initialise();

    EXPECT_EQ(data.state.solid_temperature.GetDomain(), EquationDomain::PorousMediaOnly);
    const auto ads_id_dom = data.state.GetAdsorbentIDs().front();
    EXPECT_EQ(data.state.ads_concentration.at(ads_id_dom).at("Nitrogen").GetDomain(), EquationDomain::PorousMediaOnly);

    // Global equations are AllCells by default.
    EXPECT_EQ(data.state.pressure.GetDomain(), EquationDomain::AllCells);
    EXPECT_EQ(data.state.fluid_temperature.GetDomain(), EquationDomain::AllCells);
    EXPECT_EQ(data.state.wall_temperature.GetDomain(), EquationDomain::AllCells);
    EXPECT_EQ(data.state.concentration.GetDomain(), EquationDomain::AllCells);
    EXPECT_EQ(data.state.velocity.GetDomain(), EquationDomain::AllCells);
    EXPECT_EQ(data.state.mole_fraction.at("Nitrogen").GetDomain(), EquationDomain::AllCells);
  }

  TEST_F(VoidRegionTest, IsActiveAt_PorousOnlySkipsVoidCells) {
    BuildMinimalBed(10);
    data.bed.bottom_void_.enabled_ = true;
    data.bed.top_void_.enabled_ = true;
    data.Initialise();

    const auto* mask = &data.state.domain_vectors.is_void_cell;
    const int iv_bot = data.bed.GetBottomVoidCellIndex();
    const int iv_top = data.bed.GetTopVoidCellIndex();
    const int first_bed = data.bed.GetFirstBedCellIndex();

    EXPECT_FALSE(data.state.solid_temperature.IsActiveAt(iv_bot, mask));
    EXPECT_FALSE(data.state.solid_temperature.IsActiveAt(iv_top, mask));
    EXPECT_TRUE(data.state.solid_temperature.IsActiveAt(first_bed, mask));
    EXPECT_TRUE(data.state.solid_temperature.IsActiveAt(0, mask));  // ghost is bed-side

    const auto ads_id_act = data.state.GetAdsorbentIDs().front();
    EXPECT_FALSE(data.state.ads_concentration.at(ads_id_act).at("Nitrogen").IsActiveAt(iv_bot, mask));
    EXPECT_TRUE(data.state.ads_concentration.at(ads_id_act).at("Nitrogen").IsActiveAt(first_bed, mask));

    // Global equations remain active everywhere.
    EXPECT_TRUE(data.state.pressure.IsActiveAt(iv_bot, mask));
    EXPECT_TRUE(data.state.fluid_temperature.IsActiveAt(iv_top, mask));
  }

  TEST_F(VoidRegionTest, BedCellIndices_AreContiguousBetweenFirstAndLastBed) {
    BuildMinimalBed(10);
    data.bed.bottom_void_.enabled_ = true;
    data.bed.top_void_.enabled_ = true;
    data.Initialise();

    const auto& dv = data.state.domain_vectors;
    const int first_bed = data.bed.GetFirstBedCellIndex();
    const int last_bed = data.bed.GetLastBedCellIndex();
    ASSERT_EQ((int)dv.bed_cell_indices.size(), last_bed - first_bed + 1);
    for (int i = 0; i < (int)dv.bed_cell_indices.size(); ++i) {
      EXPECT_EQ(dv.bed_cell_indices[i], first_bed + i);
    }
  }

  TEST_F(VoidRegionTest, BedPositions_StartAtZero_SpansBedLength) {
    BuildMinimalBed(10);
    data.bed.bottom_void_.enabled_ = true;
    data.bed.bottom_void_.volume_ = 1e-4;
    data.bed.top_void_.enabled_ = true;
    data.bed.top_void_.volume_ = 1e-4;
    data.Initialise();

    const auto& dv = data.state.domain_vectors;
    ASSERT_FALSE(dv.bed_positions.empty());
    EXPECT_DOUBLE_EQ(dv.bed_positions.front(), 0.0);
    // The bed-face slice covers bed inlet through bed outlet, so the last
    // entry's axial position equals (face_count - 1) * dx_bed.
    const double layer_length = data.bed.layers_.layers_[0].domain_.length.Get();
    const double dx_bed = layer_length / data.bed.layers_.layers_[0].domain_.n_cells;
    const int n = (int)dv.bed_positions.size();
    EXPECT_NEAR(dv.bed_positions.back(), (n - 1) * dx_bed, 1e-12);
  }

  TEST_F(VoidRegionTest, BedFaceIndexRange_CoversBedFaces) {
    BuildMinimalBed(10);
    data.bed.bottom_void_.enabled_ = true;
    data.bed.bottom_void_.volume_ = 1e-4;
    data.Initialise();

    const auto& dv = data.state.domain_vectors;
    const int first_bed = data.bed.GetFirstBedCellIndex();
    const int last_bed = data.bed.GetLastBedCellIndex();
    EXPECT_EQ(dv.bed_face_index_begin, std::max(0, first_bed - 1));
    EXPECT_EQ(dv.bed_face_index_end, last_bed + 1);
    EXPECT_EQ((int)dv.bed_positions.size(), dv.bed_face_index_end - dv.bed_face_index_begin);
  }

  TEST_F(VoidRegionTest, BedPositions_MatchPositionsSliceWhenVoidsDisabled) {
    BuildMinimalBed(10);
    data.Initialise();

    const auto& dv = data.state.domain_vectors;
    // With no voids, bed_positions equals the full positions vector shifted
    // to start at 0 (the original positions vector already starts at 0).
    ASSERT_EQ((int)dv.bed_positions.size(), (int)dv.positions.size());
    for (int i = 0; i < (int)dv.bed_positions.size(); ++i) {
      EXPECT_DOUBLE_EQ(dv.bed_positions[i], dv.positions[i]);
    }
  }

  TEST_F(VoidRegionTest, SystemSize_PorousOnlySkipsVoidCells) {
    BuildMinimalBed(10);
    data.bed.bottom_void_.enabled_ = true;
    data.bed.top_void_.enabled_ = true;
    data.Initialise();

    // Count expected state-vector entries directly. All-cell equations
    // contribute Size() entries each; porous-only ones contribute
    // Size() - 2 (the two void cells).
    const auto& eqs = data.state.GetEquations();
    int expected = 0;
    for (const auto* eq : eqs) {
      const int n = (int)eq->Size();
      if (eq->GetDomain() == EquationDomain::PorousMediaOnly) {
        expected += std::max(0, n - 2);
      } else {
        expected += n;
      }
    }
    EXPECT_EQ(data.state.GetSystemSize(), expected);
  }

} // end adsim namespace
