// test_heating_medium_duty.cpp
//
// Energy closure of the logged heating-medium duty. A closed, non-adsorbing
// bed with insulated wall ends exchanges heat only with the step's heating
// medium, so over each step the first law gives, with no approximation,
//
//     medium_heat_in - medium_heat_out = dU_bed + dU_wall
//
// with dU_bed and dU_wall from the independent first-principles helpers. Two
// steps: a medium ramped up from the initial temperature (heat in), then a
// cold medium (heat out). A third step without a medium logs no duty. The
// thermal-energy and cooling-duty sums (CalculateStepThermalDuty) carry the
// medium's figures, alongside the energy sources'.
//
// The same closure with the medium's internal heat exchanger on (alone, and
// with the jacket), and a check that the jacket's and the exchanger's heat
// flows add: at one state, the flow with both equals the sum of the flows
// with each alone.

#include "validation_helpers.h"

#include "adsim/solver/kpis.h"

namespace adsim {

  using namespace adsim::validation;

  // A closed, non-adsorbing bed with insulated wall ends: the step's heating
  // medium (jacket and internal exchanger) is the only way heat gets in or
  // out. Returns the initial temperature [K].
  class HeatingMediumDutyTest : public UserDataTest {
  protected:
    static constexpr double kT0 = 300.0;
    LayerID layer_id{};

    double BuildClosedBed() {
      AddComponents({"Nitrogen", "Oxygen"});
      auto ids = AddLayer("Inert Packing");
      layer_id = ids.first;
      const AdsorbentID ads_id = ids.second;
      auto& layer = GetLayer(layer_id);
      layer.domain_.length = 0.3;
      layer.domain_.n_cells = 20;
      layer.eb_ = 0.4;
      layer.particle_diameter_ = 0.002;
      std::get<CylinderBed>(data.bed.shape_).diameter = 0.03;
      auto& ads = GetAdsorbent(ads_id);
      ads.solid_density_ = 1000.0;
      ads.heat_capacity_ = 900.0;
      ads.thermal_conductivity_ = 0.0;
      ads.ep_ = 0.0;  // total voidage = eb: no extra gas inventory

      data.bed.wall_.thickness_ = 0.003;
      data.bed.wall_.density_ = 7900.0;
      data.bed.wall_.heat_capacity_ = 500.0;
      data.bed.wall_.thermal_conductivity_ = 16.2;
      // Insulated wall ends, so the medium is the only way in or out.
      data.bed.equations_.wall_temperature_boundary_condition_ =
          CreateWallTemperatureBoundaryCondition(WallTemperatureBoundaryConditionType::ZERO_FLUX);
      AddConstantInternalHeatTransferCoefficientModel(80.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);  // no ambient loss without a medium

      data.bed.ambient_temperature_ = kT0;
      AddInitialConditions(1e5, kT0, kT0, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}});
      return kT0;
    }

    void AddClosedStep(const char* name, double time) {
      AddStep(name, time);
      AddClosedBC(BoundaryLocation::BOTTOM, name);
      AddClosedBC(BoundaryLocation::TOP, name);
    }

    void UseTightTolerances() {
      data.settings.reltol = 1e-8;
      data.settings.max_step = 2.0;
    }
  };

  TEST_F(HeatingMediumDutyTest, LoggedDutyClosesTheEnergyBalance) {
    data.info.title = "Validation: heating-medium duty closes the energy balance";
    const double T0 = BuildClosedBed();

    const char* names[] = {"Heat", "Cool", "Hold"};
    for (const char* name : names) {
      AddStep(name, 200.0);
      AddClosedBC(BoundaryLocation::BOTTOM, name);
      AddClosedBC(BoundaryLocation::TOP, name);
    }
    auto& heat = data.cycle.GetStep("Heat").heating_medium_;
    heat.enabled_ = true;
    heat.heat_transfer_coefficient_ = 300.0;
    heat.temperature_profile_ = CreateInletTemperatureProfile(BoundaryValueType::LINEAR, T0, 120.0);
    static_cast<LinearBoundaryValue*>(heat.temperature_profile_.get())->end_value_ = 400.0;
    auto& cool = data.cycle.GetStep("Cool").heating_medium_;
    cool.enabled_ = true;
    cool.heat_transfer_coefficient_ = 300.0;
    cool.temperature_ = 290.0;
    ConfigureCycle(1);
    data.settings.reltol = 1e-8;
    data.settings.max_step = 2.0;

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "solver did not complete";

    const auto& log = solver.GetLog();
    const double Cv = MixtureMolarCv(data.fluid, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}});
    KPIParams params(log, solver.GetBoundaryLog());
    params.cycle_index_ = 0;

    for (int s = 0; s < 2; ++s) {
      SCOPED_TRACE(names[s]);
      const StepHeatDuty& duty = log.GetLog(0, s).heat_duty;
      const double dU = ClosedBedInternalEnergyChange(log, 0, s, Cv) + WallInternalEnergyChange(log, 0, s);
      const double net = duty.medium_heat_in - duty.medium_heat_out;
      ASSERT_GT(std::abs(dU), 1000.0) << "the column barely changed; the test proves nothing";
      EXPECT_NEAR(net, dU, 2e-4 * std::abs(dU))  // measured: 4e-6 heating, 4e-5 cooling
          << "medium in " << duty.medium_heat_in << " J, out " << duty.medium_heat_out
          << " J, dU_bed + dU_wall " << dU << " J";
      // No energy sources: their part of the duty stays zero.
      EXPECT_EQ(duty.heating_power, 0.0);
      EXPECT_EQ(duty.cooling_power, 0.0);

      // The KPI sums carry the medium's duty, and report it separately.
      const StepThermalDuty d = CalculateStepThermalDuty(params, s);
      EXPECT_DOUBLE_EQ(d.heating, duty.medium_heat_in);
      EXPECT_DOUBLE_EQ(d.cooling, duty.medium_heat_out);
      EXPECT_DOUBLE_EQ(d.medium_heating, duty.medium_heat_in);
      EXPECT_DOUBLE_EQ(d.medium_cooling, duty.medium_heat_out);
    }
    // The hot medium heats, the cold one cools.
    EXPECT_GT(log.GetLog(0, 0).heat_duty.medium_heat_in, 0.0);
    EXPECT_GT(log.GetLog(0, 1).heat_duty.medium_heat_out, 0.0);
    EXPECT_GT(log.GetLog(0, 1).heat_duty.medium_heat_out, log.GetLog(0, 1).heat_duty.medium_heat_in);

    // A step without a medium logs none.
    EXPECT_EQ(log.GetLog(0, 2).heat_duty.medium_heat_in, 0.0);
    EXPECT_EQ(log.GetLog(0, 2).heat_duty.medium_heat_out, 0.0);

    const StepThermalDuty cycle = CalculateCycleThermalDuty(params);
    EXPECT_DOUBLE_EQ(cycle.medium_heating,
                     log.GetLog(0, 0).heat_duty.medium_heat_in + log.GetLog(0, 1).heat_duty.medium_heat_in);
    EXPECT_DOUBLE_EQ(cycle.heating, cycle.medium_heating);
  }

  // The internal exchanger heats the bed directly: closure with the exchanger
  // alone (h_m = 0, no jacket) heating and cooling, then with both.
  TEST_F(HeatingMediumDutyTest, ExchangerDutyClosesTheEnergyBalance) {
    data.info.title = "Validation: internal heat exchanger duty closes the energy balance";
    BuildClosedBed();
    const char* names[] = {"ExchangerHeat", "ExchangerCool", "Both"};
    for (const char* name : names) AddClosedStep(name, 200.0);
    auto set = [&](const char* name, double T_m, double h_jacket) {
      auto& m = data.cycle.GetStep(name).heating_medium_;
      m.enabled_ = true;
      m.temperature_ = T_m;
      m.heat_transfer_coefficient_ = h_jacket;
      m.exchanger_area_per_volume_ = 50.0;           // [m2/m3]
      m.exchanger_heat_transfer_coefficient_ = 40.0;  // [W/(m2 K)]
    };
    set("ExchangerHeat", 380.0, 0.0);
    set("ExchangerCool", 290.0, 0.0);
    set("Both", 400.0, 300.0);
    ConfigureCycle(1);
    UseTightTolerances();

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "solver did not complete";

    const auto& log = solver.GetLog();
    const double Cv = MixtureMolarCv(data.fluid, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}});
    for (int s = 0; s < 3; ++s) {
      SCOPED_TRACE(names[s]);
      const StepHeatDuty& duty = log.GetLog(0, s).heat_duty;
      const double dU = ClosedBedInternalEnergyChange(log, 0, s, Cv) + WallInternalEnergyChange(log, 0, s);
      const double net = duty.medium_heat_in - duty.medium_heat_out;
      ASSERT_GT(std::abs(dU), 1000.0) << "the column barely changed; the test proves nothing";
      EXPECT_NEAR(net, dU, 2e-4 * std::abs(dU))  // measured: 1e-6, 5e-6, 3e-5
          << "medium in " << duty.medium_heat_in << " J, out " << duty.medium_heat_out
          << " J, dU_bed + dU_wall " << dU << " J";
    }
    // Exchanger alone: all of the medium's heat went through the exchanger.
    for (int s = 0; s < 2; ++s) {
      const StepHeatDuty& duty = log.GetLog(0, s).heat_duty;
      EXPECT_DOUBLE_EQ(duty.exchanger_heat_in, duty.medium_heat_in);
      EXPECT_DOUBLE_EQ(duty.exchanger_heat_out, duty.medium_heat_out);
    }
    EXPECT_GT(log.GetLog(0, 0).heat_duty.exchanger_heat_in, 0.0);
    EXPECT_GT(log.GetLog(0, 1).heat_duty.exchanger_heat_out, 0.0);
    // Both: the jacket carries part of the heat as well.
    const StepHeatDuty& both = log.GetLog(0, 2).heat_duty;
    EXPECT_GT(both.exchanger_heat_in, 0.0);
    EXPECT_GT(both.medium_heat_in, both.exchanger_heat_in);
  }

  // At one state, the heat flow with jacket and exchanger is the sum of the
  // flows with each alone, and each matches its closed form.
  TEST_F(HeatingMediumDutyTest, JacketAndExchangerFlowsAdd) {
    BuildClosedBed();
    AddClosedStep("Heat", 100.0);
    ConfigureCycle(1);
    data.Initialise();
    auto& step = data.cycle.GetStep("Heat");
    step.Initialise(data.state);
    // A non-uniform state: wall and bed temperatures that vary along the column.
    auto& Tw = data.state.wall_temperature.var;
    auto& T = data.state.fluid_temperature.var;
    for (size_t i = 0; i < T.size(); ++i) {
      T[i] = 300.0 + 2.0 * static_cast<double>(i);
      Tw[i] = 320.0 + 1.0 * static_cast<double>(i);
    }

    auto& m = step.heating_medium_;
    m.enabled_ = true;
    m.temperature_ = 330.0;  // between the coldest and hottest cells: both signs appear
    auto flow = [&](double h_jacket, double a_x) {
      m.heat_transfer_coefficient_ = h_jacket;
      m.exchanger_area_per_volume_ = a_x;
      m.exchanger_heat_transfer_coefficient_ = 40.0;
      return HeatingMediumHeatFlow(data.state, 0.0, false);
    };
    const MediumHeatFlow jacket = flow(250.0, 0.0);
    const MediumHeatFlow exchanger = flow(0.0, 50.0);
    const MediumHeatFlow both = flow(250.0, 50.0);
    EXPECT_GT(jacket.in, 0.0);
    EXPECT_GT(jacket.out, 0.0);
    EXPECT_GT(exchanger.in, 0.0);
    EXPECT_GT(exchanger.out, 0.0);
    EXPECT_EQ(jacket.exchanger_in, 0.0);
    EXPECT_NEAR(both.in, jacket.in + exchanger.in, 1e-12 * both.in);
    EXPECT_NEAR(both.out, jacket.out + exchanger.out, 1e-12 * both.out);
    EXPECT_DOUBLE_EQ(both.exchanger_in, exchanger.in);
    EXPECT_DOUBLE_EQ(both.exchanger_out, exchanger.out);

    // Closed forms: sum over the cells of pi d_o dx h (T_m - Tw) and
    // a_x h_x dV (T_m - T), each split by sign.
    const auto& dom = data.state.domain_vectors;
    double j_net = 0.0, x_net = 0.0;
    for (size_t i = 1; i + 1 < T.size(); ++i) {
      j_net += ExternalWallArea(dom, i) * 250.0 * (330.0 - Tw[i]);
      x_net += 50.0 * 40.0 * dom.volume[i] * (330.0 - T[i]);
    }
    EXPECT_NEAR(jacket.in - jacket.out, j_net, 1e-9 * std::abs(jacket.in));
    EXPECT_NEAR(exchanger.in - exchanger.out, x_net, 1e-9 * std::abs(exchanger.in));
  }

} // namespace adsim
