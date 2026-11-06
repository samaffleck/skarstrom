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

#include "validation_helpers.h"

#include "adsim/solver/kpis.h"

namespace adsim {

  using namespace adsim::validation;

  class HeatingMediumDutyTest : public UserDataTest {};

  TEST_F(HeatingMediumDutyTest, LoggedDutyClosesTheEnergyBalance) {
    data.info.title = "Validation: heating-medium duty closes the energy balance";
    AddComponents({"Nitrogen", "Oxygen"});
    auto [layer_id, ads_id] = AddLayer("Inert Packing");
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

    const double T0 = 300.0;
    data.bed.ambient_temperature_ = T0;
    AddInitialConditions(1e5, T0, T0, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}});

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

} // namespace adsim
