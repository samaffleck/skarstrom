// A fluid of three or more components with a cell of one pure component.
// Wilke's equation was 0 / 0 there and the engine returned D_m,i = 0, so every
// diffusion-based mass-transfer coefficient of that component was 0 and its
// loading never moved. With the smooth floor (WilkeMixtureDiffusivity),
// D_m,i is the harmonic mean of the partner D_ij, and the solver integrates
// through the nearly pure region.

// STL includes
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

// Gtest includes
#include <gtest/gtest.h>

// Test helpers
#include "../fixtures/user_data_test.h"

namespace adsim {

  namespace {

    const std::vector<std::string> kComponents = { "Hydrogen", "CarbonDioxide", "Nitrogen" };

    LangmuirIsothermParameters N2Langmuir() {
      LangmuirIsothermParameters p;
      p.m1 = 1.5; p.b1 = 3e-11; p.q1 = 35000;
      return p;
    }

    // H2 / CO2 / N2 on one adsorbent; macropore LDF for CO2 and N2 (D_m from
    // the gas), a constant LDF for H2.
    std::pair<LayerID, AdsorbentID> AddThreeComponentBed(UserDataTest& t) {
      t.AddComponents(kComponents);
      auto [layer_id, ads_id] = t.AddLayer("13X");
      t.GetLayer(layer_id).domain_.n_cells = 20;
      auto iso = std::make_unique<LangmuirIsotherm>();
      iso->params_["Hydrogen"].m1 = 0.05; iso->params_["Hydrogen"].b1 = 3e-11; iso->params_["Hydrogen"].q1 = 10000;
      iso->params_["CarbonDioxide"].m1 = 5.0; iso->params_["CarbonDioxide"].b1 = 3e-11; iso->params_["CarbonDioxide"].q1 = 35000;
      iso->params_["Nitrogen"] = N2Langmuir();
      t.AddIsotherm(kComponents, ads_id, std::move(iso));
      for (const auto& c : kComponents) t.SetHOA(ads_id, c, 0.0);
      t.AddLinearDrivingForceModel(ads_id, "Hydrogen", 1.0);
      t.AddLinearDrivingForceModel_Macropore(ads_id, "CarbonDioxide", 3.0, 0.35);
      t.AddLinearDrivingForceModel_Macropore(ads_id, "Nitrogen", 3.0, 0.35);
      t.AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
      t.AddConstantExternalHeatTransferCoefficientModel(0, layer_id);
      return { layer_id, ads_id };
    }

    const std::vector<double>& PairD(UserDataTest& t, const std::string& a, const std::string& b) {
      return t.data.bed.equations_.binary_diffusion_model_.m_DiffusionCoefficients.at({ a, b }).Dm_ij;
    }

  } // namespace

  // At initialisation every cell is pure N2.
  TEST_F(UserDataTest, WilkePureLimit_DiffusivityOfThePureComponent) {
    auto [layer_id, ads_id] = AddThreeComponentBed(*this);
    AddInitialConditions(1e5, 298, 298, { {"Hydrogen", 0.0}, {"CarbonDioxide", 0.0}, {"Nitrogen", 1.0} });
    data.Initialise();

    const auto& state = data.state;
    const auto& layer = GetLayer(layer_id);
    const auto& Dm_n2 = state.molecular_diffusion.at("Nitrogen").var;
    const auto& Dm_h2 = state.molecular_diffusion.at("Hydrogen").var;
    const auto& Dm_co2 = state.molecular_diffusion.at("CarbonDioxide").var;
    const auto& mtc_n2 = state.mass_transfer_coefficient.at(ads_id).at("Nitrogen").var;
    const auto& D_n2_h2 = PairD(*this, "Nitrogen", "Hydrogen");
    const auto& D_n2_co2 = PairD(*this, "Nitrogen", "CarbonDioxide");

    for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
      // Pure component: the unweighted harmonic mean of its partners' D_ij.
      const double harmonic = 2.0 / (1.0 / D_n2_h2[n] + 1.0 / D_n2_co2[n]);
      EXPECT_NEAR(Dm_n2[n], harmonic, 1e-9 * harmonic) << "cell " << n;
      EXPECT_GT(mtc_n2[n], 0.0) << "cell " << n;
      // Absent components: a trace in pure N2, D_m = D_i,N2 (Blanc's limit).
      EXPECT_NEAR(Dm_h2[n], D_n2_h2[n], 1e-5 * D_n2_h2[n]) << "cell " << n;
      EXPECT_NEAR(Dm_co2[n], D_n2_co2[n], 1e-5 * D_n2_co2[n]) << "cell " << n;
    }
  }

  // Pressurise the pure-N2 bed with N2: the N2 loading must rise to the new
  // equilibrium. With D_m,N2 = 0 the macropore coefficient was 0 and the
  // loading stayed at its initial value.
  TEST_F(UserDataTest, WilkePureLimit_PureComponentAdsorbs) {
    auto [layer_id, ads_id] = AddThreeComponentBed(*this);
    AddInitialConditions(1e5, 298, 298, { {"Hydrogen", 0.0}, {"CarbonDioxide", 0.0}, {"Nitrogen", 1.0} });

    AddStep("Pressurise", 100);
    Step& step = data.cycle.GetStep("Pressurise");
    auto inlet = CreatePressureInletBoundary(CreateLinearBoundaryValue(1e5, 3e5, 10.0),
      { {"Hydrogen", 0.0}, {"CarbonDioxide", 0.0}, {"Nitrogen", 1.0} }, 298);
    inlet->location_ = BoundaryLocation::BOTTOM;
    inlet->direction_ = BoundaryDirection::INLET;
    step.bottom_boundary_condition_ = std::move(inlet);
    AddClosedBC(BoundaryLocation::TOP, "Pressurise");
    ConfigureCycle(1);

    std::string error{};
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    CoupledSolver solver{ data };
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    double q0 = 0.0;
    IsothermFunction::Langmuir(q0, 298.0, 1e5, 1.0, N2Langmuir());
    const auto& layer = GetLayer(layer_id);
    const auto& q = data.state.ads_concentration.at(ads_id).at("Nitrogen").var;
    const auto& q_eq = data.state.eq_ads_concentration.at(ads_id).at("Nitrogen").var;
    for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
      EXPECT_GT(q[n], q0 + 0.05) << "loading did not rise in cell " << n;
      EXPECT_NEAR(q[n], q_eq[n], 0.05 * q_eq[n]) << "cell " << n;
    }
  }

  // Feed H2/CO2 into the pure-N2 bed: the front carries every cell through the
  // nearly pure region where the floor acts. The run must complete, and every
  // D_m,i must stay finite and between its smallest and largest D_ij.
  TEST_F(UserDataTest, WilkePureLimit_SolverIntegratesThroughTheFront) {
    auto [layer_id, ads_id] = AddThreeComponentBed(*this);
    AddInitialConditions(1e5, 298, 298, { {"Hydrogen", 0.0}, {"CarbonDioxide", 0.0}, {"Nitrogen", 1.0} });

    AddStep("Feed", 30.0);
    Step& feed = data.cycle.GetStep("Feed");
    auto feed_bc = std::make_unique<FlowBoundaryCondition>();
    feed_bc->direction_ = BoundaryDirection::INLET;
    feed_bc->location_ = BoundaryLocation::BOTTOM;
    feed_bc->streams_.push_back(CreateNoValveInletStream(
      CreateMolarFlowBC(0.1), { {"Hydrogen", 0.5}, {"CarbonDioxide", 0.5}, {"Nitrogen", 0.0} }, 298));
    feed.bottom_boundary_condition_ = std::move(feed_bc);
    AddPressureOutletBC(CreateConstantBoundaryValue(1e5), BoundaryLocation::TOP, "Feed",
      { {"Hydrogen", 0.0}, {"CarbonDioxide", 0.0}, {"Nitrogen", 1.0} });
    ConfigureCycle(1);

    std::string error{};
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    CoupledSolver solver{ data };
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    const auto& layer = GetLayer(layer_id);
    const auto& y_h2 = data.state.mole_fraction.at("Hydrogen").var;
    bool front_entered = false;
    for (const auto& i : kComponents) {
      const auto& Dm = data.state.molecular_diffusion.at(i).var;
      for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
        double lo = 1e300, hi = 0.0;
        for (const auto& j : kComponents) {
          if (j == i) continue;
          lo = std::min(lo, PairD(*this, i, j)[n]);
          hi = std::max(hi, PairD(*this, i, j)[n]);
        }
        ASSERT_TRUE(std::isfinite(Dm[n])) << i << " cell " << n;
        EXPECT_GE(Dm[n], lo * (1 - 1e-12)) << i << " cell " << n;
        EXPECT_LE(Dm[n], hi * (1 + 1e-12)) << i << " cell " << n;
        if (y_h2[n] > 0.1) front_entered = true;
      }
    }
    EXPECT_TRUE(front_entered) << "the feed front did not enter the bed";
  }

} // namespace adsim
