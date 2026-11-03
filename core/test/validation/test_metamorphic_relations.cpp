/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Metamorphic relations: pairs of runs that pose the same physical problem in
// two different ways and must therefore give the same answer.
//
// None of these needs a reference solution. Each test runs a deck, transforms
// it into an equivalent deck (components listed in another order, the packing
// split into identical layers, a step cut in two, the inputs written in other
// units, the cycle started at another step), runs that, and compares the two
// solutions everywhere: pressure, temperature, every mole fraction and every
// loading in every cell at the end of every step.
//
// THE TOLERANCE comes from the integrator, not from a picked number. Two runs
// of equivalent decks are two adaptive-step integrations of the same
// semi-discrete system; once any floating-point detail differs (a sum taken
// in another order, a unit factor one ulp off, a restart of the integrator)
// they take different steps, and each sits within its global error E of the
// exact semi-discrete solution. So they can differ by up to 2E, and by no
// more if the relation holds. E is measured per field by re-running the base
// deck with reltol tightened 100x: |x(reltol) - x(reltol/100)| estimates the
// error of the reltol run. To that is added one weight of IDA's error test,
// w = reltol |x| + atol (the size of one accepted local error; the absolute
// tolerances are fixed in the code, docs/numerics.md: 1e-8 for mole
// fractions and loadings, 1e-6 Pa, 1e-4 K), which covers fields whose E
// estimate is itself at round-off. A broken relation shows up as a
// discrepancy far above that, because a modelling difference does not shrink
// with the integrator tolerance.
//
// What is NOT tested here, and why: "one bed = two beds in series of half
// length". The engine simulates one column; other beds exist only through
// connections, which replay a partner step's logged stream from that step's
// previous run (docs/connecting-beds.md). Two beds operating in series AT THE
// SAME TIME would need the one simulated column to be in two states at once,
// so the series arrangement cannot be written as a deck that should equal a
// single long bed. The closest exact statement the engine can express is the
// layer split below: a layer interface is a series connection of two packings
// with no dead volume between them.

#include "breakthrough_case.h"
#include "psa_cycle_case.h"
#include "adsim/boundary_conditions/connection.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    // End-of-step snapshots of every field for every logged (cycle, step),
    // keyed "P", "T", "y:<component>", "q:<component>", stored cycle by cycle
    // with the steps in `step_order` (default: deck order), interior cells only.
    using Fields = std::map<std::string, std::vector<std::vector<double>>>;

    Fields CollectFields(const EquationLog& log, const std::vector<std::string>& components,
                         const std::vector<int>& step_order = {}) {
      Fields f;
      auto interior = [](const std::vector<double>& s) {
        return std::vector<double>(s.begin() + 1, s.end() - 1);
      };
      std::vector<int> steps(log.sys_data.step_index_set_.begin(), log.sys_data.step_index_set_.end());
      if (!step_order.empty()) steps = step_order;
      for (int cycle : log.sys_data.cycle_index_set_) {
        for (int step : steps) {
          f["P"].push_back(interior(log.GetVariable(cycle, step, VariableType::PRESSURE).values.back()));
          f["T"].push_back(interior(log.GetVariable(cycle, step, VariableType::TEMPERATURE).values.back()));
          for (const auto& c : components) {
            f["y:" + c].push_back(interior(
              log.GetComponentVariable(cycle, step, ComponentVariableType::MOLE_FRACTION, c).values.back()));
            f["q:" + c].push_back(interior(
              log.GetComponentVariable(cycle, step, ComponentVariableType::ADSORBED_CONCENTRATION, c).values.back()));
          }
        }
      }
      return f;
    }

    // Absolute tolerances of the unknowns (docs/numerics.md).
    double AbsTol(const std::string& key) {
      switch (key[0]) {
        case 'P': return 1e-6;  // [Pa]
        case 'T': return 1e-4;  // [K]
        default: return 1e-8;   // mole fraction [-], loading [mol/kg]
      }
    }

    struct FieldDiff {
      double d = 0.0;      // largest |a - b|, in the field's scale
      double weight = 0.0; // one error-test weight, reltol |x| + atol, in the same scale
    };

    // Largest |a - b| per field, in the field's natural scale: pressure and
    // loading relative to the field's largest magnitude in `a`, temperature in
    // K, mole fractions as they are. A loading that never rises above 100x its
    // absolute tolerance (a component with no isotherm) is integrator noise
    // and is left out.
    std::map<std::string, FieldDiff> FieldDistance(const Fields& a, const Fields& b, double reltol = 0.0) {
      std::map<std::string, FieldDiff> out;
      for (const auto& [key, snaps] : a) {
        const auto& other = b.at(key);
        EXPECT_EQ(snaps.size(), other.size()) << key;
        double magnitude = 0.0;
        for (const auto& s : snaps) for (double v : s) magnitude = std::max(magnitude, std::abs(v));
        if (key[0] == 'q' && magnitude < 100.0 * AbsTol(key)) continue;
        const double scale = (key[0] == 'T' || key[0] == 'y') ? 1.0 : magnitude;
        FieldDiff f;
        for (size_t n = 0; n < snaps.size() && n < other.size(); ++n) {
          EXPECT_EQ(snaps[n].size(), other[n].size()) << key;
          for (size_t i = 0; i < snaps[n].size() && i < other[n].size(); ++i)
            f.d = std::max(f.d, std::abs(snaps[n][i] - other[n][i]) / scale);
        }
        f.weight = (reltol * magnitude + AbsTol(key)) / scale;
        out[key] = f;
      }
      return out;
    }

    double ProfileRange(const Fields& f, const std::string& key) {
      double lo = 1e300, hi = -1e300;
      for (const auto& s : f.at(key)) for (double v : s) { lo = std::min(lo, v); hi = std::max(hi, v); }
      return hi - lo;
    }

    using DeckBuilder = std::function<void(UserData&, double reltol)>;
    // Which logged steps to compare (default: all, in deck order).
    using StepSelector = std::function<std::vector<int>(const UserData&)>;

    // Runs a deck and collects its fields. `ok` is false if it did not run.
    Fields RunDeck(const DeckBuilder& build, double reltol, const std::vector<std::string>& components,
                   bool& ok, const StepSelector& steps = {}) {
      UserData data;
      build(data, reltol);
      std::string error;
      ok = data.EvaluateAllExpressions(error) && data.PerformChecks(error);
      EXPECT_TRUE(ok) << error;
      if (!ok) return {};
      CoupledSolver solver{data};
      ok = solver.IntegrateCycle(data.cycle, false) >= 0;
      EXPECT_TRUE(ok) << "Solver did not complete";
      if (!ok) return {};
      return CollectFields(solver.GetLog(), components, steps ? steps(data) : std::vector<int>{});
    }

    // The relation holds if the transformed deck reproduces the base deck to
    // within 2E + k w per field (see the file comment; k = 1 unless the
    // transformation restarts the integrator). Returns the base deck's fields
    // for the caller's guards.
    Fields ExpectSameSolution(const char* relation, const DeckBuilder& base, const DeckBuilder& transformed,
                              double reltol, const std::vector<std::string>& components,
                              const StepSelector& base_steps = {}, const StepSelector& transformed_steps = {},
                              double k = 1.0) {
      bool ok_a = false, ok_fine = false, ok_b = false;
      const Fields a = RunDeck(base, reltol, components, ok_a, base_steps);
      const Fields fine = RunDeck(base, reltol * 1e-2, components, ok_fine, base_steps);
      const Fields b = RunDeck(transformed, reltol, components, ok_b, transformed_steps);
      if (!(ok_a && ok_fine && ok_b)) {
        ADD_FAILURE() << relation << ": a run failed";
        return a;
      }
      const auto error = FieldDistance(a, fine, reltol);
      const auto discrepancy = FieldDistance(a, b, reltol);
      for (const auto& [key, e] : error) {
        const double d = discrepancy.at(key).d;
        const double tolerance = 2.0 * e.d + k * e.weight;
        std::cout << "[" << relation << "] " << key << ": discrepancy " << d
                  << ", integration error " << e.d << ", tolerance " << tolerance << "\n";
        EXPECT_LE(d, tolerance)
          << relation << ": field " << key << " differs by " << d
          << ", more than the integrator tolerance " << tolerance;
      }
      return a;
    }

  } // namespace

  // ===========================================================================
  // 1. Component order
  // ===========================================================================
  //
  // Nothing in the physics depends on the order in which the components are
  // listed: every per-component quantity is looked up by name. Reversing the
  // fluid's component list on the shared PSA cycle (two adsorbates with
  // different affinities, heat of adsorption, wall heat transfer, flow
  // reversal at blowdown) must give the same solution, component by
  // component. The guard is that the fields move: the cycle swings the
  // pressure by 2 bar and the O2 mole fraction varies across the bed.

  TEST(MetamorphicRelations, ComponentOrderDoesNotChangeTheSolution) {
    auto base = [](UserData& data, double reltol) {
      BuildPsaCycle(data, 30);
      data.cycle.number_of_cycles_ = 2;
      data.settings.reltol = reltol;
      data.settings.max_step = 2.0;
    };
    auto reversed = [&](UserData& data, double reltol) {
      base(data, reltol);
      std::reverse(data.fluid.components_.begin(), data.fluid.components_.end());
    };

    UserData a, b;
    base(a, 1e-6);
    reversed(b, 1e-6);
    ASSERT_EQ(a.fluid.components_.front(), "Nitrogen");
    ASSERT_EQ(b.fluid.components_.front(), "Oxygen");  // the reversal is a real change

    const Fields f = ExpectSameSolution("component order", base, reversed, 1e-6, {"Nitrogen", "Oxygen"});
    ASSERT_FALSE(f.empty());
    EXPECT_GT(ProfileRange(f, "y:Oxygen"), 0.02) << "The cycle barely separates";
    EXPECT_GT(ProfileRange(f, "P"), 1.0e5) << "No pressure swing";
  }

  // ===========================================================================
  // 2. One layer = identical layers
  // ===========================================================================
  //
  // A packing split into two (or four) layers of the same adsorbent, with the
  // same cell width, is the same column: every face between layers must be
  // treated exactly like a face inside a layer. The interface terms are all
  // switched on (axial dispersion, heat of adsorption, bed-wall heat transfer,
  // Darcy pressure drop through the harmonic-mean face permeability), so a
  // layer interface that drops or doubles a flux, or mis-indexes a cell,
  // breaks the relation. The comparison is taken with the front mid-column, so
  // it has crossed the interfaces at z = L/4 and L/2. (Measured: bitwise
  // identical on Linux.)

  TEST(MetamorphicRelations, IdenticalLayersEqualOneLayer) {
    auto make = [](int n_layers) {
      return [n_layers](UserData& data, double reltol) {
        BreakthroughCase c;
        c.adsorbing = true;
        c.n_layers = n_layers;
        c.n_cells = 80;
        c.run_time = 75.0;  // front roughly mid-column
        c.reltol = reltol;
        c.y_feed = 0.2;
        Build(data, c);
        auto& ads = data.adsorbent_directory.Get(data.bed.layers_.layers_.front().composition_.begin()->first);
        ads.adsorption_source_.heat_of_adsorption_["Oxygen"] = 40000.0;
        for (auto& layer : data.bed.layers_.layers_) {
          auto disp = std::make_unique<DispersionModelConstant>();
          disp->dispersion_coefficient_ = 1e-4;
          layer.equations_.dispersion_model_ = std::move(disp);
          auto ihtc = std::make_unique<ConstantInternalHeatTransferModel>();
          ihtc->heat_transfer_coefficient_ = 20.0;
          layer.equations_.internal_heat_transfer_model_ = std::move(ihtc);
        }
        data.bed.ambient_temperature_ = c.bed_temperature;
      };
    };
    const std::vector<std::string> components{"Nitrogen", "Oxygen"};

    {
      UserData four;
      make(4)(four, 1e-6);
      ASSERT_EQ(four.bed.layers_.layers_.size(), 4u);
      ASSERT_EQ(four.bed.layers_.layers_.back().domain_.n_cells, 20);
    }

    const Fields f = ExpectSameSolution("two layers", make(1), make(2), 1e-6, components);
    ExpectSameSolution("four layers", make(1), make(4), 1e-6, components);

    ASSERT_FALSE(f.empty());
    const auto& y = f.at("y:Oxygen").back();
    EXPECT_GT(*std::max_element(y.begin(), y.end()), 0.9 * 0.2) << "Front has not entered the bed";
    EXPECT_LT(*std::min_element(y.begin(), y.end()), 0.1 * 0.2) << "Front has left the bed";
    EXPECT_GT(ProfileRange(f, "T"), 0.5) << "No thermal front";
  }

  // ===========================================================================
  // 3. A step = two half-length steps
  // ===========================================================================
  //
  // Cutting a step in two, with the same boundary conditions on both halves,
  // changes nothing physical. It does change the numerics: every cycle step
  // restarts the integrator (SolverData::RestartForNewStep) from a fresh
  // history at first order and the initial step size, after re-solving the
  // algebraic part of the state. So the relation is approximate, by what the
  // restart may legitimately change: one accepted error test, whose WRMS norm
  // over the N unknowns is at most 1, which allows a single unknown an error
  // of up to sqrt(N) of its weight w. The tolerance is 2E + sqrt(N) w
  // (N ~ 650 here). Measured on y_O2: 3.7e-8, about 0.6 w, with
  // E = 5.6e-9. (On a 150 s run cut at 75 s the discrepancy was 1.5e-7,
  // 3.1e-8 and 6.3e-8 at reltol 1e-6, 1e-8 and 1e-10, against w = 6e-8,
  // 1.05e-8 and 1.0e-8: it follows the error-test weight down to the
  // absolute-tolerance floor, the signature of an integrator effect rather
  // than a modelling one.) The end of the run (the
  // last step of each deck) is compared, with the front inside the bed both
  // at the cut and at the end.

  TEST(MetamorphicRelations, SplitStepEqualsSingleStep) {
    auto make = [](bool split) {
      return [split](UserData& data, double reltol) {
        BreakthroughCase c;
        c.adsorbing = true;
        c.run_time = 100.0;
        c.reltol = reltol;
        Build(data, c);
        if (split) {
          data.DuplicateStep(0);
          data.cycle.step_list_[0].step_time_ = 50.0;
          data.cycle.step_list_[1].step_time_ = 50.0;
        }
      };
    };
    auto last_step = [](const UserData& data) {
      return std::vector<int>{static_cast<int>(data.cycle.step_list_.size()) - 1};
    };

    UserData probe;
    make(true)(probe, 1e-6);
    ASSERT_EQ(probe.cycle.step_list_.size(), 2u);
    std::string error;
    ASSERT_TRUE(probe.PerformChecks(error)) << error;
    CoupledSolver sizing{probe};
    ASSERT_GE(sizing.IntegrateCycle(probe.cycle, false), 0);
    const double n_unknowns = probe.state.GetSystemSize();
    ASSERT_GT(n_unknowns, 0.0);

    const Fields f = ExpectSameSolution("split step", make(false), make(true), 1e-6,
                                        {"Nitrogen", "Oxygen"}, last_step, last_step,
                                        std::sqrt(n_unknowns));
    ASSERT_FALSE(f.empty());
    const auto& y = f.at("y:Oxygen").back();
    EXPECT_GT(*std::max_element(y.begin(), y.end()), 0.9 * 0.05);
    EXPECT_LT(*std::min_element(y.begin(), y.end()), 0.1 * 0.05);
  }

  // ===========================================================================
  // 4. Units
  // ===========================================================================
  //
  // The same deck written with its inputs in other units -- cm, mm, bar,
  // mbar, degC, min, g/cm3, kJ/(kg K), mol/h -- is the same physical deck.
  // The unit parser converts each to SI before the solver sees it, so the
  // converted values must equal the SI ones (to round-off) and the solutions
  // must agree. The guard checks each expression really carried a unit.

  TEST(MetamorphicRelations, InputsInOtherUnitsGiveTheSameSolution) {
    auto make = [](bool units) {
      return [units](UserData& data, double reltol) {
        BreakthroughCase c;
        c.adsorbing = true;
        c.run_time = 75.0;
        c.reltol = reltol;
        c.bed_temperature = 298.15;
        c.feed_temperature = 298.15;
        c.molar_flow = 0.005;  // = 18 mol/h
        Build(data, c);
        if (!units) return;
        auto& layer = data.bed.layers_.layers_.front();
        layer.domain_.length.SetExpr("100 [cm]");
        layer.particle_diameter_.SetExpr("2 [mm]");
        std::get<CylinderBed>(data.bed.shape_).diameter.SetExpr("10 [cm]");
        data.ic.pressure.SetExpr("1 [bar]");
        data.ic.temperature.SetExpr("25 [degC]");
        data.ic.wall_temperature.SetExpr("25 [degC]");
        Step& step = data.cycle.step_list_[0];
        step.step_time_.SetExpr("1.25 [min]");
        auto* inlet_bc = dynamic_cast<FlowBoundaryCondition*>(step.bottom_boundary_condition_.get());
        auto* inlet = dynamic_cast<BoundaryStream*>(inlet_bc->streams_.front().get());
        inlet->temperature_.SetExpr("25 [degC]");
        auto* flow = dynamic_cast<MolarFlowBoundaryValue*>(inlet->boundary_flow_.get());
        static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_.SetExpr("18 [mol/h]");
        auto* outlet = dynamic_cast<PressureBoundaryCondition*>(step.top_boundary_condition_.get());
        static_cast<ConstantBoundaryValue*>(outlet->pressure_.get())->constant_value_.SetExpr("1000 [mbar]");
        auto& ads = data.adsorbent_directory.Get(layer.composition_.begin()->first);
        ads.solid_density_.SetExpr("1 [g/cm3]");
        ads.heat_capacity_.SetExpr("0.9 [kJ/(kg*K)]");
        data.bed.wall_.thickness_.SetExpr("17.5 [mm]");
      };
    };

    // Guard: the expressions carry units and convert to the SI values.
    {
      UserData si, units;
      make(false)(si, 1e-6);
      make(true)(units, 1e-6);
      std::string error;
      ASSERT_TRUE(si.EvaluateAllExpressions(error)) << error;
      ASSERT_TRUE(units.EvaluateAllExpressions(error)) << error;
      auto& lu = units.bed.layers_.layers_.front();
      auto& ls = si.bed.layers_.layers_.front();
      auto& au = units.adsorbent_directory.Get(lu.composition_.begin()->first);
      auto& as = si.adsorbent_directory.Get(ls.composition_.begin()->first);
      const std::vector<std::pair<const ExpressionParam*, const ExpressionParam*>> pairs{
        {&lu.domain_.length, &ls.domain_.length},
        {&lu.particle_diameter_, &ls.particle_diameter_},
        {&std::get<CylinderBed>(units.bed.shape_).diameter, &std::get<CylinderBed>(si.bed.shape_).diameter},
        {&units.ic.pressure, &si.ic.pressure},
        {&units.ic.temperature, &si.ic.temperature},
        {&units.cycle.step_list_[0].step_time_, &si.cycle.step_list_[0].step_time_},
        {&au.solid_density_, &as.solid_density_},
        {&au.heat_capacity_, &as.heat_capacity_},
        {&units.bed.wall_.thickness_, &si.bed.wall_.thickness_},
      };
      for (const auto& [u, s] : pairs) {
        EXPECT_TRUE(u->HasInlineUnit()) << u->GetExpr();
        EXPECT_NE(u->GetExpr(), s->GetExpr());
        EXPECT_NEAR(u->Get(), s->Get(), 1e-14 * std::abs(s->Get())) << u->GetExpr();
      }
    }

    ExpectSameSolution("units", make(false), make(true), 1e-6, {"Nitrogen", "Oxygen"});
  }

  // ===========================================================================
  // 5. Cycle phase: the cyclic steady state does not depend on the first step
  // ===========================================================================
  //
  // The engine solves ONE bed; a multi-bed cycle is that bed's step sequence,
  // and the other beds are the same bed at other phases of the cycle, reached
  // through connections that replay a partner step's last run. "A symmetric
  // N-bed cycle has identical, phase-shifted beds at CSS" is therefore true by
  // construction, and what can be tested is its consequence: the cyclic steady
  // state is a periodic orbit, so it cannot depend on which step the deck
  // lists first. Starting the cycle two steps later (the other bed's phase of
  // a two-bed cycle) must reach the same CSS, step for step.
  //
  // The cycle (isothermal, so the CSS is reached in tens of cycles rather than
  // hundreds): pressurise (pressure ramp to 3 bar at the feed end), feed
  // (product at 3 bar), blowdown to 1 bar, purge with the feed step's whole
  // product replayed into the top. The connection is what makes the rotation
  // non-trivial: listed after the feed, the purge reads the feed of the SAME
  // cycle; listed before it (the rotated deck), it reads the PREVIOUS cycle's
  // feed. The two decks agree only at CSS.
  //
  // Tolerance: each run's own cycle-to-cycle change over its last cycles. A
  // run cannot be expected to agree with another better than with itself one
  // cycle later; both runs are within that change of the CSS (it is at the
  // integrator's noise floor here). Guards: the first cycles of the two runs
  // differ by orders of magnitude more (the transients are different), and
  // the cycle separates strongly.

  namespace {
    void BuildConnectedCycle(UserData& data) {
      BuildPsaCycle(data, 30);
      data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
      data.cycle.GetStep("Feed").step_time_ = 60.0;
      data.cycle.GetStep("Blowdown").step_time_ = 60.0;
      {
        auto bc = std::make_unique<PressureBoundaryCondition>();
        bc->direction_ = BoundaryDirection::INLET;
        bc->location_ = BoundaryLocation::BOTTOM;
        auto ramp = std::make_unique<LinearBoundaryValue>();  // from the bed pressure, over the step
        ramp->end_value_ = 3.0e5;
        bc->pressure_ = std::move(ramp);
        bc->temperature_ = 298.0;
        bc->mole_fraction_["Nitrogen"] = 0.79;
        bc->mole_fraction_["Oxygen"] = 0.21;
        data.cycle.GetStep("Pressurise").bottom_boundary_condition_ = std::move(bc);
      }
      data.cycle.AddStep("Purge");
      Step& purge = data.cycle.GetStep("Purge");
      purge.step_time_ = 60.0;  // a connection pairs steps of equal duration
      auto out = std::make_unique<PressureBoundaryCondition>();
      out->direction_ = BoundaryDirection::OUTLET;
      out->location_ = BoundaryLocation::BOTTOM;
      out->pressure_ = std::make_unique<ConstantBoundaryValue>(1.0e5);
      out->mole_fraction_["Nitrogen"] = 0.79;
      out->mole_fraction_["Oxygen"] = 0.21;
      purge.bottom_boundary_condition_ = std::move(out);
      // First-cycle assumption of the purge inlet, used until the feed has run.
      auto stream = std::make_unique<BoundaryStream>();
      auto flow = std::make_unique<MolarFlowBoundaryValue>();
      static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = 0.005;
      stream->boundary_flow_ = std::move(flow);
      stream->mole_fraction["Nitrogen"] = 0.5;
      stream->mole_fraction["Oxygen"] = 0.5;
      stream->temperature_ = 298.0;
      stream->name_ = "Purge in";
      IBoundaryStream* purge_in = stream.get();
      auto in = std::make_unique<FlowBoundaryCondition>();
      in->direction_ = BoundaryDirection::INLET;
      in->location_ = BoundaryLocation::TOP;
      in->streams_.push_back(std::move(stream));
      purge.top_boundary_condition_ = std::move(in);
      data.cycle.EnsureIds();
      auto* product = dynamic_cast<PressureBoundaryCondition*>(data.cycle.GetStep("Feed").top_boundary_condition_.get());
      auto a = MakeEndpoint(data.cycle.GetStep("Feed"), *product);
      auto b = MakeEndpoint(data.cycle.GetStep("Purge"), *purge_in);
      ASSERT_TRUE(MakeMutualConnection(a, b));
    }
  } // namespace

  TEST(MetamorphicRelations, RotatedCycleReachesTheSameCyclicSteadyState) {
    const int n_cycles = 30;
    const std::vector<std::string> names{"Pressurise", "Feed", "Blowdown", "Purge"};

    // Per cycle, the end-of-step fields with the steps in `names` order.
    auto run = [&](bool rotate) {
      UserData data;
      BuildConnectedCycle(data);
      if (rotate) {  // -> Blowdown, Purge, Pressurise, Feed
        for (int r = 0; r < 2; ++r) { data.SwapSteps(0, 1); data.SwapSteps(1, 2); data.SwapSteps(2, 3); }
        EXPECT_EQ(data.cycle.step_list_.front().step_name_, "Blowdown");
        EXPECT_EQ(data.cycle.step_list_.back().step_name_, "Feed");
      }
      data.cycle.number_of_cycles_ = n_cycles;
      data.settings.reltol = 1e-8;
      data.settings.max_step = 2.0;
      data.settings.max_cycles_to_keep_ = n_cycles;
      std::vector<Fields> per_cycle;
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      CoupledSolver solver{data};
      const int flag = solver.IntegrateCycle(data.cycle, false);
      EXPECT_GE(flag, 0) << "Solver did not complete";
      if (flag < 0) return per_cycle;
      std::vector<int> order;
      for (const auto& n : names)
        for (size_t k = 0; k < data.cycle.step_list_.size(); ++k)
          if (data.cycle.step_list_[k].step_name_ == n) order.push_back(static_cast<int>(k));
      const auto& log = solver.GetLog();
      const Fields all = CollectFields(log, {"Nitrogen", "Oxygen"}, order);
      for (size_t c = 0; c < log.sys_data.cycle_index_set_.size(); ++c) {
        Fields f;
        for (const auto& [key, snaps] : all)
          f[key] = std::vector<std::vector<double>>(snaps.begin() + c * names.size(),
                                                    snaps.begin() + (c + 1) * names.size());
        per_cycle.push_back(std::move(f));
      }
      return per_cycle;
    };

    const auto base = run(false);
    const auto rotated = run(true);
    ASSERT_EQ(base.size(), static_cast<size_t>(n_cycles));
    ASSERT_EQ(rotated.size(), static_cast<size_t>(n_cycles));

    // Guard: strong separation (O2 enriched to > 90 % at the product end).
    double y_max = 0.0;
    for (const auto& s : base.back().at("y:Oxygen")) y_max = std::max(y_max, *std::max_element(s.begin(), s.end()));
    ASSERT_GT(y_max, 0.9) << "The cycle barely separates";

    // Each run's own cycle-to-cycle change over its last five cycles.
    const int window = 5;
    std::map<std::string, double> self_change;
    for (const auto* r : {&base, &rotated})
      for (int c = n_cycles - window; c < n_cycles; ++c)
        for (const auto& [key, d] : FieldDistance((*r)[c - 1], (*r)[c]))
          self_change[key] = std::max(self_change[key], d.d);

    const auto first_cycle = FieldDistance(base.front(), rotated.front());
    const auto css = FieldDistance(base.back(), rotated.back());
    for (const auto& [key, diff] : css) {
      std::cout << "[cycle rotation] " << key << ": CSS discrepancy " << diff.d << ", first-cycle difference "
                << first_cycle.at(key).d << ", cycle-to-cycle change " << self_change.at(key) << "\n";
      EXPECT_LE(diff.d, 2.0 * self_change.at(key))
        << "Field " << key << " of the rotated cycle's CSS differs from the base CSS by " << diff.d;
    }
    // The transients differ: the comparison is not between two copies of one run.
    EXPECT_GT(first_cycle.at("y:Oxygen").d, 100.0 * self_change.at("y:Oxygen"));
  }

} // namespace adsim
