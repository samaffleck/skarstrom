/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Component mass balance with axial dispersion (Danckwerts inlet faces).
//
// The species flux through an inlet face of the bed must be the feed's
// u C y_in, the amount the stream logs (and so the KPIs and the CSS mass-balance
// check) record. When the residual also let a dispersive flux
// eb D C (y_1 - y_in) / dx cross the inlet face, species entered (or left) the
// bed unseen by the logs: the total moles still balanced (the dispersive fluxes
// sum to zero over the components) but each component did not, by an amount
// proportional to D and almost independent of the grid. At cyclic steady state
// (CO2 out) / (CO2 in) of a 13X CO2/N2 cycle was 1.0035 at D = 3.1e-4 m2/s and
// 1.032 at D = 3e-3 m2/s.
//
// Two checks on a dispersive single-bed cycle with a flow-condition inlet at the
// bottom (pressurise, feed) and a pressure-condition inlet at the top (purge),
// so both boundary-condition kinds and both bed ends are covered:
//  * per step and component, the moles the boundary streams log entering the
//    bed equal the change of the bed's gas + adsorbed inventory (holds in every
//    cycle, no CSS needed);
//  * at cyclic steady state, out / in of every component over the cycle is 1.
//
// GUARD: dispersion must matter. Each run checks that the bed's first cell
// differs from the feed composition by a clear margin during the feed and the
// purge, so a dispersive boundary flux, if present, would be large.

#include "psa_cycle_case.h"

#include "adsim/solver/kpis.h"

namespace adsim {

  using namespace adsim::validation;

  namespace {
    constexpr double kDispersion = 3.0e-3;       // [m2/s]: Pe = u L / D ~ 7 in the feed
    constexpr double kHeatDispersion = 1.0;      // [W/(m.K)]
    constexpr double kPurgePressure = 1.001e5;   // [Pa]: ~0.02 m/s through the bed at 1 bar
    constexpr double kPurgeO2 = 0.9;

    // The fixture's pressurise / feed / blowdown cycle with constant axial
    // mass and heat dispersion and a purge step that feeds O2-rich gas into the
    // top through a pressure condition.
    void BuildDispersiveCycle(UserData& data, int n_cells) {
      BuildPsaCycle(data, n_cells);

      for (auto& layer : data.bed.layers_.layers_) {
        auto disp = std::make_unique<DispersionModelConstant>();
        disp->dispersion_coefficient_ = kDispersion;
        layer.equations_.dispersion_model_ = std::move(disp);
        auto heat_disp = std::make_unique<HeatDispersionModelConstant>();
        heat_disp->dispersion_coefficient_ = kHeatDispersion;
        layer.equations_.heat_dispersion_model_ = std::move(heat_disp);
      }

      data.cycle.AddStep("Purge");
      Step& purge = data.cycle.GetStep("Purge");
      purge.step_time_ = 60.0;

      auto inlet = std::make_unique<PressureBoundaryCondition>();
      inlet->direction_ = BoundaryDirection::INLET;
      inlet->location_ = BoundaryLocation::TOP;
      auto p_in = std::make_unique<ConstantBoundaryValue>();
      p_in->constant_value_ = kPurgePressure;
      inlet->pressure_ = std::move(p_in);
      inlet->temperature_ = 298.0;
      inlet->mole_fraction_["Oxygen"] = kPurgeO2;
      inlet->mole_fraction_["Nitrogen"] = 1.0 - kPurgeO2;
      purge.top_boundary_condition_ = std::move(inlet);

      auto outlet = std::make_unique<PressureBoundaryCondition>();
      outlet->direction_ = BoundaryDirection::OUTLET;
      outlet->location_ = BoundaryLocation::BOTTOM;
      auto p_out = std::make_unique<ConstantBoundaryValue>();
      p_out->constant_value_ = 1.0e5;
      outlet->pressure_ = std::move(p_out);
      outlet->mole_fraction_["Oxygen"] = 0.21;
      outlet->mole_fraction_["Nitrogen"] = 0.79;
      purge.bottom_boundary_condition_ = std::move(outlet);

      // The checks integrate the logged stream flows in time (trapezoidal
      // rule over the solver's steps); a capped step keeps that quadrature
      // error (3e-4 of the throughput at max_step = 1 s) well below the
      // tolerance.
      data.settings.reltol = 1e-8;
      data.settings.max_step = 0.1;
    }

    // Moles of `component` the bed gained over one step (gas + adsorbed).
    double InventoryChange(const EquationLog& log, int cycle, int step_index, const std::string& component) {
      return CalculateComponentAccumulation_Volume(log, cycle, step_index,
                                                   ComponentVariableType::COMPONENT_MOLAR_CONCENTRATION, component)
           + CalculateComponentAccumulation_Mass(log, cycle, step_index,
                                                 ComponentVariableType::ADSORBED_CONCENTRATION, component);
    }

    // Moles of `component` entering the bed through each boundary stream of
    // one step, from the streams' logs (what the KPIs integrate): the logged
    // molar flow has the velocity sign (positive bottom -> top), so a top
    // stream is negated. Positive = entering the bed.
    std::vector<double> StreamMolesIn(const BoundaryLog& blog, const Step& step, int cycle, int step_index,
                                      const std::string& component) {
      std::vector<double> moles_in;
      auto add = [&](const IBoundaryCondition* bc, double sign) {
        if (!bc) return;
        std::vector<std::string> names;
        if (auto fbc = dynamic_cast<const FlowBoundaryCondition*>(bc)) {
          for (const auto& stream : fbc->streams_) names.push_back(stream->name_);
        } else if (dynamic_cast<const PressureBoundaryCondition*>(bc)) {
          names.push_back(step.step_name_ + " Pressure Boundary " + BoundaryLocationList[(int)bc->location_]);
        }
        for (const auto& name : names) {
          const auto d = blog.GetComponentVariableData(cycle, step_index, name,
                                                       BoundaryComponentVariableType::MOLAR_FLOW, component);
          moles_in.push_back(sign * IntegrateVariableInTime(d.variable_timeseries, d.time));
        }
      };
      add(step.bottom_boundary_condition_.get(), +1.0);
      add(step.top_boundary_condition_.get(), -1.0);
      return moles_in;
    }

    // Largest |y_first_cell - y_feed| of Oxygen over a step, at the given end
    // (first interior cell 1 or n_cells), from the equation log.
    double MaxInletCellDeparture(const EquationLog& log, int cycle, int step_index, bool top, double y_feed) {
      const auto y = log.GetComponentVariable(cycle, step_index, ComponentVariableType::MOLE_FRACTION, "Oxygen");
      double worst = 0.0;
      for (const auto& snapshot : y.values) {
        const size_t i = top ? snapshot.size() - 2 : 1;
        worst = std::max(worst, std::abs(snapshot[i] - y_feed));
      }
      return worst;
    }
  } // namespace

  TEST(DispersionMassBalance, LoggedStreamsAccountForEverySpeciesInEveryStep) {
    UserData data;
    BuildDispersiveCycle(data, 30);
    data.cycle.number_of_cycles_ = 3;

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    const auto& log = solver.GetLog();
    const auto& blog = solver.GetBoundaryLog();
    const int last = data.cycle.number_of_cycles_ - 1;

    // Guard: the first cell departs from the feed by a clear margin, so a
    // dispersive inlet flux would be a sizeable fraction of the feed.
    EXPECT_GT(MaxInletCellDeparture(log, last, 1, false, 0.21), 0.05) << "feed: first cell sits at the feed";
    EXPECT_GT(MaxInletCellDeparture(log, last, 3, true, kPurgeO2), 0.05) << "purge: first cell sits at the feed";

    for (int cycle = 0; cycle <= last; ++cycle) {
      for (int s = 0; s < static_cast<int>(data.cycle.step_list_.size()); ++s) {
        const Step& step = data.cycle.step_list_[s];
        for (const auto& component : data.fluid.components_) {
          const auto streams = StreamMolesIn(blog, step, cycle, s, component);
          const double logged_in = std::accumulate(streams.begin(), streams.end(), 0.0);
          double throughput = 0.0;
          for (double moles : streams) throughput += std::abs(moles);
          const double accumulated = InventoryChange(log, cycle, s, component);
          // Scale: the moles that crossed the boundaries or the inventory
          // change, whichever is larger. The tolerance allows for the stream
          // log's first point, which holds the pre-switch flow: when Feed
          // opens the top at 3 bar the flow jumps to ~1000 mol/s within 3e-7 s
          // and the quadrature misses ~1.5e-4 of the step's throughput. The
          // unlogged dispersive flux this guards against was 4e-3 to 0.5.
          const double scale = std::max({throughput, std::abs(accumulated), 1e-3});
          EXPECT_LT(std::abs(logged_in - accumulated) / scale, 5e-4)
            << component << ", cycle " << cycle << ", step " << step.step_name_
            << ": streams logged " << logged_in << " mol in, bed inventory changed by " << accumulated << " mol";
        }
      }
    }
  }

  // The reporter's measure: (out) / (in) of each component over a cycle at
  // cyclic steady state, summed over every boundary stream. Isothermal, so the
  // cycle reaches CSS in a few dozen cycles (the wall temperature is the slow
  // mode of the non-isothermal run); CSS is asserted directly as a vanishing
  // per-component inventory change over the last cycle, since the molar-profile
  // criterion looks at the total concentration only.
  TEST(DispersionMassBalance, ComponentBalanceClosesAtCyclicSteadyState) {
    UserData data;
    BuildDispersiveCycle(data, 20);
    data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
    data.cycle.number_of_cycles_ = 60;

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    const auto& log = solver.GetLog();
    const auto& blog = solver.GetBoundaryLog();
    const int last = data.cycle.number_of_cycles_ - 1;
    for (const auto& component : data.fluid.components_) {
      double in = 0.0;
      double out = 0.0;
      double inventory_change = 0.0;
      for (int s = 0; s < static_cast<int>(data.cycle.step_list_.size()); ++s) {
        for (double moles : StreamMolesIn(blog, data.cycle.step_list_[s], last, s, component)) {
          (moles > 0.0 ? in : out) += std::abs(moles);
        }
        inventory_change += InventoryChange(log, last, s, component);
      }
      ASSERT_GT(in, 0.1) << component;
      ASSERT_LT(std::abs(inventory_change) / in, 1e-5)
        << component << ": not at CSS, the bed gained " << inventory_change << " mol over the last cycle";
      // The tolerance is the stream logs' time-quadrature error, as in the
      // per-step test: the first point of a step holds the pre-switch flow,
      // and how much that misses depends on IDA's first step after the switch,
      // which differs between compilers (O2: 1.7e-5 with GCC, 1.0e-4 with
      // MSVC, with the inflow identical to 7 digits). The unlogged dispersive
      // flux this guards against gave O2 1.053 and N2 0.983.
      EXPECT_NEAR(out / in, 1.0, 5e-4)
        << component << " over cycle " << last + 1 << ": in " << in << " mol, out " << out << " mol";
    }
  }

} // namespace adsim
