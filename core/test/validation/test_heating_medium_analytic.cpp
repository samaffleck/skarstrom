// test_heating_medium_analytic.cpp
//
// Analytic heat-up by a heating medium. A closed, non-adsorbing bed with
// insulated wall ends is spatially uniform, so every cell follows the same
// linear ODEs and the response to a medium at fixed T_m has a closed form.
// Three cases, each compared over the whole step with the exact solution:
//
//  1. The jacket heats the wall alone (bed decoupled, h_i = 0):
//       rho_w c_w dTw/dt = K_o (T_m - Tw),  K_o = 4 d_o h_m / (d_o^2 - d_i^2)
//     an exponential approach to T_m with time constant rho_w c_w / K_o.
//
//  2. The jacket heats the wall, which heats the bed (two coupled nodes):
//       C_b dT/dt  = (4 h_i / d_i) (Tw - T)
//       C_w dTw/dt = K_i (T - Tw) + K_o (T_m - Tw),  K_i = 4 d_i h_i / (d_o^2 - d_i^2)
//     solved exactly with the 2x2 matrix exponential.
//
//  3. An internal exchanger heats the bed alone (h_m = 0, h_i = 0):
//       C_b dT/dt = a_x h_x (T_m - T)
//     an exponential with time constant C_b / (a_x h_x).
//
// C_b is the bed's heat capacity per volume at constant volume,
// eb C Cv + rho_b Cps: in a closed bed the solver's P-V work term turns the
// gas Cp into Cv (see ClosedBedInternalEnergyChange). The expected values
// come only from these closed forms and the inputs, never from a solver run.
// The bounds are regression bounds a few times the measured error (relative
// to T_m - T0, at reltol 1e-8).

#include "validation_helpers.h"

#include <array>

namespace adsim {

  using namespace adsim::validation;

  namespace {
    constexpr double kT0 = 300.0;     // [K] initial bed and wall temperature
    constexpr double kTm = 380.0;     // [K] medium temperature
    constexpr double kP0 = 1.0e5;     // [Pa]
    constexpr double kDiameter = 0.03;
    constexpr double kVoidFraction = 0.4;
    constexpr double kSolidDensity = 1000.0;
    constexpr double kSolidCp = 900.0;
    constexpr double kWallThickness = 0.003;
    constexpr double kWallDensity = 7900.0;
    constexpr double kWallCp = 500.0;
  }

  class HeatingMediumAnalyticTest : public UserDataTest {
  protected:
    LayerID layer_id{};

    void Build(double h_internal) {
      AddComponents({"Nitrogen", "Oxygen"});
      auto ids = AddLayer("Inert Packing");
      layer_id = ids.first;
      auto& layer = GetLayer(layer_id);
      layer.domain_.length = 0.3;
      layer.domain_.n_cells = 20;
      layer.eb_ = kVoidFraction;
      layer.particle_diameter_ = 0.002;
      std::get<CylinderBed>(data.bed.shape_).diameter = kDiameter;
      auto& ads = GetAdsorbent(ids.second);
      ads.solid_density_ = kSolidDensity;
      ads.heat_capacity_ = kSolidCp;
      ads.thermal_conductivity_ = 0.0;
      ads.ep_ = 0.0;  // total voidage = eb

      data.bed.wall_.thickness_ = kWallThickness;
      data.bed.wall_.density_ = kWallDensity;
      data.bed.wall_.heat_capacity_ = kWallCp;
      data.bed.wall_.thermal_conductivity_ = 16.2;
      data.bed.equations_.wall_temperature_boundary_condition_ =
          CreateWallTemperatureBoundaryCondition(WallTemperatureBoundaryConditionType::ZERO_FLUX);
      AddConstantInternalHeatTransferCoefficientModel(h_internal, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);
      data.bed.ambient_temperature_ = kT0;
      AddInitialConditions(kP0, kT0, kT0, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}});

      AddStep("Heat", 300.0);
      AddClosedBC(BoundaryLocation::BOTTOM, "Heat");
      AddClosedBC(BoundaryLocation::TOP, "Heat");
      ConfigureCycle(1);
      data.settings.reltol = 1e-8;
      data.settings.max_step = 1.0;
    }

    HeatingMedium& Medium() {
      auto& m = data.cycle.GetStep("Heat").heating_medium_;
      m.enabled_ = true;
      m.temperature_ = kTm;
      return m;
    }

    const EquationLog& Run(CoupledSolver& solver) {
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "solver did not complete";
      return solver.GetLog();
    }

    // Geometry and heat capacities from the inputs (not from the solver).
    static double InnerDiameter() { return kDiameter; }
    static double OuterDiameter() { return kDiameter + 2.0 * kWallThickness; }
    static double AnnulusFactor() {
      return OuterDiameter() * OuterDiameter() - InnerDiameter() * InnerDiameter();
    }
    double BedHeatCapacity() const {  // [J/(m3 K)], constant volume
      const double C0 = kP0 / (kR * kT0);
      const double Cv = MixtureMolarCv(data.fluid, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}});
      return kVoidFraction * C0 * Cv + (1.0 - kVoidFraction) * kSolidDensity * kSolidCp;
    }
    static double WallHeatCapacity() { return kWallDensity * kWallCp; }

    // Largest deviation of a logged temperature (mean over the interior
    // cells) from `exact(t)`, relative to T_m - T0.
    static double MaxRelativeError(const VariableLogEntry& entry, const std::function<double(double)>& exact) {
      double worst = 0.0;
      for (size_t n = 0; n < entry.values.size(); ++n) {
        double sum = 0.0;
        const auto& row = entry.values[n];
        for (size_t i = 1; i + 1 < row.size(); ++i) sum += row[i];
        const double T = sum / static_cast<double>(row.size() - 2);
        worst = std::max(worst, std::abs(T - exact(entry.times[n] - entry.times.front())));
      }
      return worst / (kTm - kT0);
    }
  };

  TEST_F(HeatingMediumAnalyticTest, JacketHeatsTheWallExponentially) {
    data.info.title = "Validation: jacket heats the wall exponentially";
    Build(0.0);  // bed decoupled from the wall
    const double h_m = 200.0;
    Medium().heat_transfer_coefficient_ = h_m;
    CoupledSolver solver{data};
    const auto& log = Run(solver);

    const double K_o = 4.0 * OuterDiameter() * h_m / AnnulusFactor();
    const double tau = WallHeatCapacity() / K_o;
    auto exact = [&](double t) { return kTm + (kT0 - kTm) * std::exp(-t / tau); };
    const auto Tw = log.GetVariable(0, 0, VariableType::WALL_TEMPERATURE);
    ASSERT_FALSE(Tw.values.empty());
    ASSERT_GT(300.0 / tau, 1.0) << "the step should cover more than one time constant";
    EXPECT_LT(MaxRelativeError(Tw, exact), 2e-5);  // measured 5e-6

    // The bed, decoupled, stays at T0.
    const auto T = log.GetVariable(0, 0, VariableType::TEMPERATURE);
    EXPECT_LT(MaxRelativeError(T, [](double) { return kT0; }), 1e-9);
  }

  TEST_F(HeatingMediumAnalyticTest, JacketHeatsWallAndBedAsTwoCoupledNodes) {
    data.info.title = "Validation: jacket heats wall and bed (two nodes)";
    const double h_i = 60.0;
    Build(h_i);
    const double h_m = 300.0;
    Medium().heat_transfer_coefficient_ = h_m;
    CoupledSolver solver{data};
    const auto& log = Run(solver);

    // y = (T, Tw) - T_m obeys y' = A y.
    const double a = 4.0 * h_i / InnerDiameter() / BedHeatCapacity();
    const double b = 4.0 * InnerDiameter() * h_i / AnnulusFactor() / WallHeatCapacity();
    const double c = 4.0 * OuterDiameter() * h_m / AnnulusFactor() / WallHeatCapacity();
    const std::array<std::array<double, 2>, 2> A = {{{-a, a}, {b, -(b + c)}}};
    const double tr = A[0][0] + A[1][1];
    const double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    const double disc = std::sqrt(tr * tr / 4.0 - det);  // real: A is similar to a symmetric matrix
    const double l1 = tr / 2.0 + disc, l2 = tr / 2.0 - disc;
    // exp(A t) = [e^{l1 t} (A - l2 I) - e^{l2 t} (A - l1 I)] / (l1 - l2)
    auto exact = [&](double t, int row) {
      const double e1 = std::exp(l1 * t), e2 = std::exp(l2 * t);
      const double y0 = kT0 - kTm;  // both nodes start at T0
      double y = 0.0;
      for (int col = 0; col < 2; ++col) {
        const double m1 = A[row][col] - (row == col ? l2 : 0.0);
        const double m2 = A[row][col] - (row == col ? l1 : 0.0);
        y += (e1 * m1 - e2 * m2) / (l1 - l2) * y0;
      }
      return kTm + y;
    };
    const auto T = log.GetVariable(0, 0, VariableType::TEMPERATURE);
    const auto Tw = log.GetVariable(0, 0, VariableType::WALL_TEMPERATURE);
    ASSERT_FALSE(T.values.empty());
    // The bed moves appreciably, so the coupling is exercised.
    ASSERT_GT(exact(300.0, 0) - kT0, 0.1 * (kTm - kT0));
    EXPECT_LT(MaxRelativeError(T, [&](double t) { return exact(t, 0); }), 2e-6);   // measured 3e-7
    EXPECT_LT(MaxRelativeError(Tw, [&](double t) { return exact(t, 1); }), 2e-6);  // measured 5e-7
  }

  TEST_F(HeatingMediumAnalyticTest, InternalExchangerHeatsTheBedExponentially) {
    data.info.title = "Validation: internal exchanger heats the bed exponentially";
    Build(0.0);  // no bed-wall exchange
    auto& m = Medium();
    m.heat_transfer_coefficient_ = 0.0;  // no jacket
    m.exchanger_area_per_volume_ = 50.0;
    m.exchanger_heat_transfer_coefficient_ = 40.0;
    CoupledSolver solver{data};
    const auto& log = Run(solver);

    const double tau = BedHeatCapacity() / (50.0 * 40.0);
    auto exact = [&](double t) { return kTm + (kT0 - kTm) * std::exp(-t / tau); };
    const auto T = log.GetVariable(0, 0, VariableType::TEMPERATURE);
    ASSERT_FALSE(T.values.empty());
    ASSERT_GT(300.0 / tau, 0.5);
    EXPECT_LT(MaxRelativeError(T, exact), 5e-7);  // measured 8e-8

    // The wall, with neither jacket nor bed contact, stays at T0.
    const auto Tw = log.GetVariable(0, 0, VariableType::WALL_TEMPERATURE);
    EXPECT_LT(MaxRelativeError(Tw, [](double) { return kT0; }), 1e-9);
  }

} // namespace adsim
