/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Conservation of the dispersive and conductive face fluxes across a layer
// interface where the cell width and the bed voidage both change.
//
// A finite-volume scheme conserves a quantity when the flux through each face
// is one value: what leaves one cell enters its neighbour. The dispersive
// flux of the component balance was eb_n d_f C_f (y_{n+1} - y_n) / dx_n in
// cell n, but eb_{n+1} d_f C_f (y_{n+1} - y_n) / dx_{n+1} in cell n + 1, and
// the conduction of the energy balances and of the wall also divided by the
// cell's own width. Where the width or eb changes (a layer interface), one
// cell lost a different amount from what its neighbour gained: with
// dx_n = 2 dx_{n+1}, twice as much. The single-layer tests (and the method of
// manufactured solutions on one uniform layer) cannot see this, because every
// face there is between two equal cells.
//
// The checks, on a column of a coarse, dense layer (8 cells, eb 0.35) under a
// fine, open one (20 cells, eb 0.48), with a different dispersion coefficient
// and conductivity in each:
//
//  1. Residual level, exact: with no flow, no accumulation and no sources,
//     the residuals hold the face-flux differences alone, so their sum over
//     the cells (in the balance's own units) must vanish to round-off. Every
//     flux path is checked: the species dispersion, the gas conduction in LTE
//     and in LTNE, and the wall conduction; also with void regions at both
//     ends (cells of another width again).
//  2. Solved, species: a closed column with dispersion across the interface
//     (after a feed step leaves an O2 front at the interface) holds the
//     moles of each species, gas and adsorbed, over the closed step, to the
//     integrator's tolerance.
//  3. Solved, energy: a closed column heated in its coarse layer, with gas
//     and wall conduction across the interface (LTE and LTNE), gains exactly
//     the heater's energy in its gas, solid and wall.
//
// The references are written from the physics (ideal gas, cylinder
// geometry); nothing here shares code with core/src.

#include "validation_helpers.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    const std::string kO2 = "Oxygen";
    const std::string kN2 = "Nitrogen";

    constexpr double kT0 = 298.0;      // [K]
    constexpr double kP0 = 1.0e5;      // [Pa]
    constexpr double kCpN2 = 29.1;     // constant molar gas heat capacities [J/(mol K)]
    constexpr double kCpO2 = 29.4;
    constexpr double kDiameter = 0.05; // [m]

    // One layer of the column.
    struct LayerSpec {
      double length;      // [m]
      int cells;
      double eb;          // [-]
      double dispersion;  // constant axial dispersion coefficient [m2/s]
      double lambda;      // constant axial thermal conductivity [W/(m K)]
      double henry_O2;    // Henry constant of O2 [mol/(kg Pa)], 0 = inert solid
      double henry_N2;
    };

    // Bottom: coarse cells (dx = 0.025 m), dense packing, strong dispersion
    // and conduction. Top: fine cells (dx = 0.01 m, 2.5 times smaller), open
    // packing, weaker dispersion and conduction.
    struct Column {
      LayerSpec bottom{0.2, 8, 0.35, 2.0e-3, 2.0, 0.3e-6, 0.1e-6};
      LayerSpec top{0.2, 20, 0.48, 6.0e-4, 0.4, 0.6e-6, 0.2e-6};
      bool adsorbing = true;
      bool ltne = false;
      bool voids = false;      // a void region at each end
      double hi = 0.0;         // bed -> wall [W/(m2 K)]
    };

    void BuildColumn(UserData& data, const Column& c) {
      data.fluid.AddComponent(kN2);
      data.fluid.AddComponent(kO2);
      data.fluid.viscosity_mode_ = ViscosityMode::CONSTANT;
      data.fluid.viscosity_ = 1.8e-5;
      data.fluid.gas_heat_capacity_mode_ = GasHeatCapacityMode::CONSTANT;
      data.fluid.gas_heat_capacity_[kN2] = kCpN2;
      data.fluid.gas_heat_capacity_[kO2] = kCpO2;
      data.fluid.adsorbed_heat_capacity_mode_ = AdsorbedHeatCapacityMode::CONSTANT;
      data.fluid.adsorbed_heat_capacity_[kN2] = 0.0;
      data.fluid.adsorbed_heat_capacity_[kO2] = 0.0;
      data.fluid.heat_of_adsorption_convention_ = HeatOfAdsorptionConvention::ISOSTERIC;

      std::get<CylinderBed>(data.bed.shape_).diameter = kDiameter;
      data.bed.wall_.thickness_ = 0.004;
      data.bed.wall_.density_ = 7900.0;
      data.bed.wall_.heat_capacity_ = 500.0;
      data.bed.wall_.thermal_conductivity_ = 16.0;
      data.bed.ambient_temperature_ = kT0;
      // Insulated wall ends: no heat leaves the wall axially.
      data.bed.equations_.wall_temperature_boundary_condition_ =
        CreateWallTemperatureBoundaryCondition(WallTemperatureBoundaryConditionType::ZERO_FLUX);
      if (c.voids) {
        data.bed.bottom_void_.enabled_ = true;
        data.bed.bottom_void_.volume_ = 6.0e-5;  // [m3]: a cell 3 cm long
        data.bed.top_void_.enabled_ = true;
        data.bed.top_void_.volume_ = 1.0e-5;     // [m3]: 5 mm
      }

      auto add_layer = [&](const std::string& name, const LayerSpec& s) {
        const auto ads_id = data.adsorbent_directory.Add(name + " solid");
        auto& ads = data.adsorbent_directory.Get(ads_id);
        ads.solid_density_ = 1100.0;
        ads.heat_capacity_ = 900.0;
        ads.thermal_conductivity_ = 0.5;
        ads.ep_ = 0.3;
        if (c.adsorbing) {
          for (const auto& [comp, k] : {std::pair{kO2, s.henry_O2}, std::pair{kN2, s.henry_N2}}) {
            auto henry = std::make_unique<HenryIsotherm>();
            henry->params_[comp].k = k;
            henry->params_[comp].Q = 0.0;
            ads.adsorption_source_.AddIsotherm({comp}, std::move(henry));
          }
          for (const auto& comp : {kO2, kN2}) {
            auto ldf = std::make_unique<LinearDrivingForce>();
            auto mtc = std::make_unique<MassTransferCoefficientConstant>();
            mtc->constant_mtc_ = 0.5;
            ldf->mtc_model_ = ToResistances(std::move(mtc));
            ads.adsorption_source_.mass_transfer_models_[comp] = std::move(ldf);
            ads.adsorption_source_.heat_of_adsorption_[comp] = 15000.0;
          }
        }

        const auto lid = data.bed.layers_.Add(name);
        auto& layer = data.bed.layers_.Get(lid);
        layer.momentum_model_ = MomentumModel::KarmenCozeny;
        layer.domain_.length = s.length;
        layer.domain_.n_cells = s.cells;
        layer.eb_ = s.eb;
        layer.particle_diameter_ = 0.002;
        layer.composition_[ads_id] = 1.0;
        auto disp = std::make_unique<DispersionModelConstant>();
        disp->dispersion_coefficient_ = s.dispersion;
        layer.equations_.dispersion_model_ = std::move(disp);
        auto heat_disp = std::make_unique<HeatDispersionModelConstant>();
        heat_disp->dispersion_coefficient_ = s.lambda;
        layer.equations_.heat_dispersion_model_ = std::move(heat_disp);
        auto ihtc = std::make_unique<ConstantInternalHeatTransferModel>();
        ihtc->heat_transfer_coefficient_ = c.hi;
        layer.equations_.internal_heat_transfer_model_ = std::move(ihtc);
        auto ehtc = std::make_unique<ConstantExternalHeatTransferModel>();
        ehtc->heat_transfer_coefficient_ = 0.0;  // no loss to the surroundings
        layer.equations_.external_heat_transfer_model_ = std::move(ehtc);
        return lid;
      };
      add_layer("Coarse", c.bottom);
      add_layer("Fine", c.top);

      if (c.ltne) {
        data.bed.equations_.thermal_equilibrium_mode_ = ThermalEquilibriumMode::NON_LOCAL;
        auto php = std::make_unique<ConstantParticleHeatTransferModel>();
        php->heat_transfer_coefficient_ = 5.0e3;  // [W/(m3 K)]
        data.bed.equations_.particle_heat_transfer_model_ = std::move(php);
      }

      data.ic.pressure = kP0;
      data.ic.temperature = kT0;
      data.ic.wall_temperature = kT0;
      data.ic.velocity = 0.0;
      data.ic.mole_fraction[kN2] = 1.0;
      data.ic.mole_fraction[kO2] = 0.0;
      data.cycle.number_of_cycles_ = 1;
    }

    void SealStep(Step& step) {
      auto bottom = std::make_unique<ClosedBoundaryCondition>();
      bottom->location_ = BoundaryLocation::BOTTOM;
      step.bottom_boundary_condition_ = std::move(bottom);
      auto top = std::make_unique<ClosedBoundaryCondition>();
      top->location_ = BoundaryLocation::TOP;
      step.top_boundary_condition_ = std::move(top);
    }

    // Index of the first cell of the top layer (cell 0 is the ghost).
    int FirstTopCell(const Column& c) { return 1 + (c.voids ? 1 : 0) + c.bottom.cells; }

    // =========================================================================
    // 1. Residual level
    // =========================================================================

    // Sums of a balance's face-flux differences over the interior cells.
    struct FluxSum {
      double net = 0.0;    // sum of the cells' net inflows (0 when conservative)
      double scale = 0.0;  // sum of their magnitudes
      double Relative() const { return std::abs(net) / scale; }
    };

    struct ResidualSums {
      FluxSum species_O2, species_N2, gas_energy, wall_energy;
    };

    // Evaluates the residuals of a still column (u = 0, every time derivative
    // 0, no sources, no bed-wall exchange) whose composition, gas temperature
    // and wall temperature vary along it, so that only the dispersive and
    // conductive face fluxes are left in the residuals.
    ResidualSums StillColumnResiduals(const Column& c) {
      UserData data;
      BuildColumn(data, c);
      data.cycle.AddStep("Hold");
      Step& step = data.cycle.GetStep("Hold");
      step.step_time_ = 1.0;
      SealStep(step);
      std::string error;
      EXPECT_TRUE(data.EvaluateAllExpressions(error) && data.PerformChecks(error)) << error;
      data.Initialise();
      data.cycle.step_list_.front().Initialise(data.state);
      data.state.domain_vectors.PopulateEnergySources(step.energy_sources_, data.bed);

      auto& s = data.state;
      const auto& dom = s.domain_vectors;
      const int N = static_cast<int>(s.pressure.var.size());
      EXPECT_EQ(N, FirstTopCell(c) + c.top.cells + (c.voids ? 1 : 0) + 1);

      // Smooth, non-trivial profiles: a composition front, a temperature
      // hump and a wall gradient, each with a slope at the interface.
      double z = 0.0;
      for (int i = 0; i < N; ++i) {
        const double zc = (i == 0) ? 0.0 : (i == N - 1 ? z : z + 0.5 * dom.cell_width[i]);
        if (i > 0 && i < N - 1) z += dom.cell_width[i];
        const double y = 0.5 + 0.4 * std::tanh((zc - 0.18) / 0.05);
        s.mole_fraction.at(kO2).var[i] = y;
        s.mole_fraction.at(kN2).var[i] = 1.0 - y;
        s.pressure.var[i] = kP0 * (1.0 + 0.1 * zc);
        s.fluid_temperature.var[i] = kT0 + 40.0 * std::exp(-std::pow((zc - 0.15) / 0.08, 2));
        s.solid_temperature.var[i] = s.fluid_temperature.var[i];
        s.wall_temperature.var[i] = kT0 + 30.0 * std::sin(6.0 * zc);
      }
      // Mirror the end cells into the ghosts (no gradient at the bed ends;
      // the dispersion and gas conduction skip those faces anyway).
      for (auto* v : {&s.mole_fraction.at(kO2).var, &s.mole_fraction.at(kN2).var, &s.fluid_temperature.var,
                      &s.solid_temperature.var, &s.wall_temperature.var, &s.pressure.var}) {
        (*v)[0] = (*v)[1];
        (*v)[N - 1] = (*v)[N - 2];
      }
      std::fill(s.velocity.var.begin(), s.velocity.var.end(), 0.0);

      s.UpdateIndependentEquations(&data);
      for (const auto& comp : {kO2, kN2}) {
        UpdateMoleFractionFace(&data, comp);
        UpdateMoleFraction(&data, comp);
      }
      if (c.ltne) UpdateFluidTemperature(&data);
      else UpdateTemperature(&data);
      UpdateWallTemperature(&data);

      // The guards: the coefficients really jump at the interface.
      const int a = FirstTopCell(c) - 1, b = FirstTopCell(c);
      EXPECT_GT(dom.cell_width[a], 2.0 * dom.cell_width[b]);
      EXPECT_GT(dom.eb[b], 1.3 * dom.eb[a]);
      EXPECT_GT(s.axial_dispersion_coefficient.var[a], 3.0 * s.axial_dispersion_coefficient.var[b]);
      EXPECT_GT(s.axial_thermal_conductivity.var[a], 3.0 * s.axial_thermal_conductivity.var[b]);

      // Each residual is the face-flux difference in its own scaling:
      //   species: R T_n / P_n times the molar inflow per unit area [m/s];
      //   gas energy: minus the heat inflow [W] (cell area times flux);
      //   wall energy: minus the heat inflow per unit wall volume [W/m3].
      ResidualSums out;
      for (int i = 1; i < N - 1; ++i) {
        for (const auto& comp : {kO2, kN2}) {
          const double r = s.mole_fraction.at(comp).res_var[i];
          const double moles = r * s.pressure.var[i] / (kR * s.fluid_temperature.var[i]);  // [mol/(m2 s)]
          FluxSum& f = (comp == kO2) ? out.species_O2 : out.species_N2;
          f.net += moles;
          f.scale += std::abs(moles);
        }
        const double q = s.fluid_temperature.res_var[i];
        out.gas_energy.net += q;
        out.gas_energy.scale += std::abs(q);
        const double qw = s.wall_temperature.res_var[i] * WallCellVolume(dom, i);
        out.wall_energy.net += qw;
        out.wall_energy.scale += std::abs(qw);
      }
      return out;
    }

  } // namespace

  // The residual of every cell of a still column is the difference of its two
  // face fluxes; summed over the column the interior faces cancel and the two
  // bed-end faces carry nothing (Danckwerts for the gas, insulated wall
  // ends), so the sum vanishes to round-off in each balance.
  //
  // Measured: at most 1.5e-16 (species), 2.4e-17 (gas energy) and 0 (wall)
  // of the summed magnitudes.
  //
  // GUARD: when each cell divided by its own width, took its own eb and the
  // arithmetic mean conductivity (before the fix), the sums were 0.33 (O2
  // and N2), 0.16 (gas energy, LTE and LTNE) and 0.11 (wall) of the summed
  // magnitudes; with the void regions 0.20, 0.12 and 0.044.
  TEST(ConservativeFaceFluxes, StillColumnResidualsSumToZeroAcrossTheInterface) {
    for (bool voids : {false, true}) {
      for (bool ltne : {false, true}) {
        Column c;
        c.adsorbing = false;
        c.voids = voids;
        c.ltne = ltne;
        const ResidualSums r = StillColumnResiduals(c);
        const std::string what = std::string(ltne ? "LTNE" : "LTE") + (voids ? ", void regions" : "");
        std::printf("[still column, %s] O2 %.2e  N2 %.2e  gas energy %.2e  wall %.2e\n", what.c_str(),
                    r.species_O2.Relative(), r.species_N2.Relative(), r.gas_energy.Relative(),
                    r.wall_energy.Relative());
        for (const auto* f : {&r.species_O2, &r.species_N2, &r.gas_energy, &r.wall_energy})
          EXPECT_GT(f->scale, 0.0) << what << ": a flux path is idle";
        EXPECT_LT(r.species_O2.Relative(), 1e-13) << what << ": O2 dispersion is not conservative";
        EXPECT_LT(r.species_N2.Relative(), 1e-13) << what << ": N2 dispersion is not conservative";
        EXPECT_LT(r.gas_energy.Relative(), 1e-13) << what << ": gas conduction is not conservative";
        EXPECT_LT(r.wall_energy.Relative(), 1e-13) << what << ": wall conduction is not conservative";
      }
    }
  }

  // ===========================================================================
  // 2. Solved: the moles of each species in a closed column
  // ===========================================================================
  //
  // A feed step pushes O2 up the column (bottom inlet, top outlet) until its
  // front sits at the interface; then both ends close and dispersion carries
  // O2 up and N2 down across the interface for a minute. Over the closed step
  // the moles of each species in the gas (between and inside the particles)
  // and on the two solids must not change. The closure is set by the
  // integrator (the inventory P y / (R T) is not linear in the unknowns, so
  // IDA keeps it only to its tolerance): measured 4e-11 to 1.8e-10 of the
  // inventory at reltol 1e-10 in LTE and LTNE.
  //
  // GUARD: dispersion carries O2 across the interface during the closed step
  // (y_O2 of the first fine cell rises from 0.06 to 0.28), and when the face
  // flux was not single valued (before the fix) the closed column gained
  // 7.1 % of its O2 and lost 45 % of its N2.
  TEST(ConservativeFaceFluxes, ClosedColumnKeepsTheMolesOfEachSpecies) {
    for (bool ltne : {false, true}) {
      Column c;
      c.ltne = ltne;
      c.hi = 20.0;
      UserData data;
      BuildColumn(data, c);

      data.cycle.AddStep("Feed");
      Step& feed = data.cycle.GetStep("Feed");
      feed.step_time_ = 30.0;
      auto inlet = std::make_unique<BoundaryStream>();
      inlet->name_ = "Feed inlet";
      auto flow = std::make_unique<MolarFlowBoundaryValue>();
      static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = 5.0e-3;  // [mol/s]
      inlet->boundary_flow_ = std::move(flow);
      inlet->mole_fraction[kN2] = 0.2;
      inlet->mole_fraction[kO2] = 0.8;
      inlet->temperature_ = kT0;
      auto inlet_bc = std::make_unique<FlowBoundaryCondition>();
      inlet_bc->direction_ = BoundaryDirection::INLET;
      inlet_bc->location_ = BoundaryLocation::BOTTOM;
      inlet_bc->streams_.push_back(std::move(inlet));
      feed.bottom_boundary_condition_ = std::move(inlet_bc);
      auto outlet_bc = std::make_unique<PressureBoundaryCondition>();
      outlet_bc->direction_ = BoundaryDirection::OUTLET;
      outlet_bc->location_ = BoundaryLocation::TOP;
      outlet_bc->pressure_ = std::make_unique<ConstantBoundaryValue>(kP0);
      outlet_bc->mole_fraction_[kN2] = 1.0;
      outlet_bc->mole_fraction_[kO2] = 0.0;
      feed.top_boundary_condition_ = std::move(outlet_bc);

      data.cycle.AddStep("Hold");
      Step& hold = data.cycle.GetStep("Hold");
      hold.step_time_ = 60.0;
      SealStep(hold);

      data.settings.reltol = 1e-10;
      data.settings.max_step = 0.5;

      std::string error;
      ASSERT_TRUE(data.EvaluateAllExpressions(error) && data.PerformChecks(error)) << error;
      CoupledSolver solver{data};
      ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
      const auto& log = solver.GetLog();

      // Guard: the O2 front is at the interface when the column closes (y_O2
      // drops by more than 0.2 over the two coarse cells below it), and
      // dispersion carries O2 across it during the closed step (the first
      // fine cell gains more than 0.2).
      const auto y = log.GetComponentVariable(0, 1, ComponentVariableType::MOLE_FRACTION, kO2);
      const int a = FirstTopCell(c) - 1, b = FirstTopCell(c);
      std::printf("[closed column] y_O2 below | above the interface: %.3f | %.3f at the start, %.3f | %.3f at the end\n",
                  y.values.front()[a - 1], y.values.front()[b], y.values.back()[a - 1], y.values.back()[b]);
      ASSERT_GT(y.values.front()[a - 1] - y.values.front()[b], 0.2) << "no O2 front at the interface";
      ASSERT_GT(y.values.back()[b] - y.values.front()[b], 0.2) << "no O2 crossed the interface";

      const char* mode = ltne ? "LTNE" : "LTE";
      for (const auto& comp : {kO2, kN2}) {
        const double before = GasMoles(log, 0, 1, comp, false) + AdsorbedMoles(log, 0, 1, comp, false);
        const double after = GasMoles(log, 0, 1, comp, true) + AdsorbedMoles(log, 0, 1, comp, true);
        const double rel = std::abs(after - before) / before;
        std::printf("[closed column, %s] %s: %.9e -> %.9e mol, relative change %.2e\n", mode, comp.c_str(),
                    before, after, rel);
        EXPECT_LT(rel, 1e-9) << mode << ", " << comp << ": the closed column gained or lost moles";
      }
    }
  }

  // ===========================================================================
  // 3. Solved: the energy of a closed, heated column
  // ===========================================================================
  //
  // A closed column of inert solids, with no loss to the surroundings, is
  // heated in its coarse bottom layer for two minutes. The heat crosses the
  // interface by gas conduction and, through the bed-wall exchange, by
  // conduction along the wall. The first law over gas, solid and wall reads
  //
  //   Q_heater = sum_i V_i [ (Cv / R) et_i dP_i + rho_b,i Cps dT_s,i ]
  //            + sum_i V_w,i rho_w Cp_w dT_w,i,
  //
  // exact for an ideal gas of constant Cv (pure N2 here): its internal energy
  // is Cv n T = (Cv / R) P V whatever the gas moves between cells. Measured
  // 1.0e-8 of the heater's energy at reltol 1e-10, in LTE and LTNE. It is not
  // round-off: besides the integrator's tolerance, the energy balance is
  // written in advective form, rho Cp u_mean (T_e - T_w) with the cell's own
  // rho Cp (numerics.md), which does not conserve the enthalpy the gas
  // carries between cells exactly; the heating moves little gas here.
  //
  // GUARD: heat crosses the interface (the first fine cell warms by 1.6 K
  // and the wall over it by 0.4 K), and when the conductive face fluxes were
  // not single valued (before the fix) the column gained 1.9 % more energy
  // than the heater gave, in LTE and LTNE.
  TEST(ConservativeFaceFluxes, ClosedHeatedColumnGainsTheHeatersEnergy) {
    for (bool ltne : {false, true}) {
      Column c;
      c.adsorbing = false;
      c.ltne = ltne;
      c.hi = 40.0;
      UserData data;
      BuildColumn(data, c);

      constexpr double kPower = 15.0;   // [W]
      constexpr double kTime = 120.0;   // [s]
      data.cycle.AddStep("Heat");
      Step& heat = data.cycle.GetStep("Heat");
      heat.step_time_ = kTime;
      SealStep(heat);
      EnergySource heater;
      heater.enabled_ = true;
      heater.power_ = kPower;
      heat.energy_sources_[data.bed.layers_.layers_.front().id_] = heater;

      data.settings.reltol = 1e-10;
      data.settings.max_step = 0.5;

      std::string error;
      ASSERT_TRUE(data.EvaluateAllExpressions(error) && data.PerformChecks(error)) << error;
      CoupledSolver solver{data};
      ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
      const auto& log = solver.GetLog();
      const auto& dom = log.GetDomain();

      const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
      const auto Ts = log.GetVariable(0, 0, ltne ? VariableType::SOLID_TEMPERATURE : VariableType::TEMPERATURE);
      const auto Tw = log.GetVariable(0, 0, VariableType::WALL_TEMPERATURE);
      ASSERT_FALSE(P.values.empty());
      const double Cv = kCpN2 - kR;
      double dU_gas = 0.0, dU_solid = 0.0, dU_wall = 0.0;
      for (size_t i = 1; i + 1 < P.values.front().size(); ++i) {
        dU_gas += dom.volume[i] * dom.et[i] * (Cv / kR) * (P.values.back()[i] - P.values.front()[i]);
        dU_solid += dom.volume[i] * dom.bulk_density[i] * dom.Cps[i] * (Ts.values.back()[i] - Ts.values.front()[i]);
        dU_wall += WallCellVolume(dom, i) * dom.rho_w[i] * dom.Cpw[i] * (Tw.values.back()[i] - Tw.values.front()[i]);
      }
      const double Q = EnergySourceHeat(log, 0, 0);
      ASSERT_NEAR(Q, kPower * kTime, 1e-9 * kPower * kTime) << "heater duty read back from the log";

      // Guard: heat crossed the interface in the bed and in the wall.
      const int b = FirstTopCell(c);
      ASSERT_GT(Ts.values.back()[b] - kT0, 1.0) << "the top layer did not warm";
      ASSERT_GT(Tw.values.back()[b] - kT0, 0.25) << "the wall over the top layer did not warm";

      const double imbalance = (dU_gas + dU_solid + dU_wall - Q) / Q;
      std::printf("[heated column, %s] Q %.6f J: gas %.6f, solid %.6f, wall %.6f; relative imbalance %.2e\n",
                  ltne ? "LTNE" : "LTE", Q, dU_gas, dU_solid, dU_wall, imbalance);
      EXPECT_LT(std::abs(imbalance), 1e-6) << (ltne ? "LTNE" : "LTE") << ": the column's energy does not close";
    }
  }

} // namespace adsim
