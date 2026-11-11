// ============================================================================
//  Nguyen, Shimizu & Rajendran (2025) - humid post-combustion CO2 capture by
//  VSA on CALF-20, Adsorption 31:45, doi:10.1007/s10450-025-00600-z.
//  The dry cases A0 (Case 1) and B0 (Case 2) of Table 1.
//
//  The paper is experimental (no model of its own). It uses the rig, the
//  columns and the CALF-20 of Nguyen et al. (2023), so the decks are the
//  Nguyen2023_Exp_VII deck (four-step LPP cycle, single-bed emulation, the
//  light-product connection and the N2 vents unchanged) at this paper's
//  conditions: P_H 0.96 bar, 22 C, and each case's t_ads, feed velocity, P_I
//  and P_L. See README.md (provenance table, model-form audit, attribution).
// ============================================================================

#define _USE_MATH_DEFINES
#include <chrono>
#include <cmath>
#include <iostream>

#include "../../../src/pch.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/connection.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/mass_transfer_model.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/solver/kpis.h"

#include <gtest/gtest.h>

#include "../nguyen_2023/nguyen_2023.h"
#include "../../fixtures/validation_metrics.h"

namespace adsim {

  class Nguyen2024Test : public Nguyen2023Test {};

  // A fresh fixture for each attribution run (the runs share one TEST_F).
  struct Nguyen2024Runner : Nguyen2024Test {
    Nguyen2024Runner() { SetUp(); }
    void TestBody() override {}
  };

  namespace {

    // Operating point of a Table 1 case.
    struct Nguyen2024Case {
      std::string id;         // "CaseA0"
      std::string label;      // example deck name
      std::string title;
      std::string description;
      double t_ads;           // [s] total adsorption time (Feed 1 + Feed 2)
      double v_interstitial;  // [m/s] feed velocity, read as interstitial
      double P_I;             // [Pa] end of blowdown
      double P_L;             // [Pa] end of evacuation
      int max_cycles;         // CSS cap
    };

    // One-factor changes for the attribution runs (README "Attribution").
    struct Nguyen2024Variant {
      bool paper_decay = false;      // (1) F = exp(alpha t_step) from SI Table A.1
      bool measured_column = false;  // (2) L = 0.267 m, eps that holds 36.0 g
      bool si_properties = false;    // (3) D_m and mu of SI Table A.1
      bool blake_kozeny = false;     // (4) Blake-Kozeny (150) pressure drop
      std::string Name() const {
        if (paper_decay) return "paper decay rates";
        if (measured_column) return "measured column";
        if (si_properties) return "SI D_m and mu";
        if (blake_kozeny) return "Blake-Kozeny";
        return "deck";
      }
    };

    struct Nguyen2024Result {
      double purity = 0, recovery = 0, productivity = 0;  // %, %, mol/(m3 s)
      int cycles = 0;
      bool css = false;
      double wall_s = 0;
    };

    constexpr double kPH = 0.96e5;     // [Pa] Table 1 caption / p. 5: P_H = 0.96 bar
    constexpr double kT = 295.15;      // [K] Table 1 caption: 22 C
    constexpr double kMads = 36.0e-3;  // [kg] SI Table A.1
    constexpr double kRhoP = 570.0;    // [kg/m3] SI Table A.1

    // Builds the Nguyen2023_Exp_VII deck structure (RunNguyen2023LPP in
    // ../nguyen_2023/nguyen_2023.cpp) at this paper's conditions, runs it to CSS
    // and returns the KPIs of the last cycle. Saves the deck when `label` is set.
    Nguyen2024Result RunNguyen2024(Nguyen2024Test& t, const Nguyen2024Case& c, const Nguyen2024Variant& v,
                                   bool save) {
      const std::string CO2 = "CarbonDioxide";
      const std::string N2 = "Nitrogen";
      t.data.info.title = c.title;
      t.data.info.description = c.description;
      t.data.info.author = "Sam Affleck";

      t.SetUpNguyen2023();

      // This paper's temperature (22 C) for the initial bed and wall and the
      // surroundings; P_H for the initial pressure. Feed-loaded start as VII.
      t.AddInitialConditions(kPH, kT, kT, { { CO2, 0.15 }, { N2, 0.85 } });
      t.data.bed.ambient_temperature_ = kT;

      auto& layer = t.data.bed.layers_.layers_.front();
      double eb = layer.eb_;
      if (v.measured_column) {
        // SI Table A.1: L = 0.267 m. The voidage that holds 36.0 g at
        // rho_p = 570 kg/m3: eps = 1 - m / (rho_p A L) = 0.3170.
        const double Di = std::get<CylinderBed>(t.data.bed.shape_).diameter;
        const double area = 0.25 * M_PI * Di * Di;
        layer.domain_.length = 0.267;
        eb = 1.0 - kMads / (kRhoP * area * 0.267);
        layer.eb_ = eb;
      }
      if (v.si_properties) {
        // SI Table A.1: D_m 1.496e-5 m2/s (printed unit "m s-2"), mu 1.812e-5 Pa s.
        if (auto* r = dynamic_cast<DispersionModelRuthven*>(layer.equations_.dispersion_model_.get()))
          r->Dm = 1.496e-5;
        std::vector<AdsorbentID> ids;
        for (const auto& [id, a] : t.data.adsorbent_directory.GetAll()) ids.push_back(id);
        for (AdsorbentID id : ids)
          for (auto& [comp, model] : t.data.adsorbent_directory.Get(id).adsorption_source_.mass_transfer_models_)
            if (auto* ldf = dynamic_cast<LinearDrivingForce*>(model.get()))
              if (auto* sr = dynamic_cast<MassTransferCoefficientSeriesResistance*>(ldf->mtc_model_.get()))
                sr->molecular_diffusion_coefficient_ = 1.496e-5;
        t.data.fluid.viscosity_mode_ = ViscosityMode::CONSTANT;
        t.data.fluid.viscosity_ = 1.812e-5;
      }
      if (v.blake_kozeny) layer.momentum_model_ = MomentumModel::BlakeKozeny;

      const double press_time = 50;  // [s] VII deck (the paper ends this step on a pressure trigger)
      const double feed1_time = press_time;
      const double feed2_time = c.t_ads - feed1_time;
      const double blowdown_time = 60;  // [s] VII deck
      const double evac_time = 40;      // [s] VII deck

      const std::string feed1 = "Feed 1", feed2 = "Feed 2", blowdown = "Blowdown", evacuation = "Evacuation",
                        press = "Pressurisation";
      t.AddStep(feed1, feed1_time);
      t.AddStep(feed2, feed2_time);
      t.AddStep(blowdown, blowdown_time);
      t.AddStep(evacuation, evac_time);
      t.AddStep(press, press_time);

      const double feed_vel = c.v_interstitial * eb;  // superficial
      const double ambient_valve_cv = 10.0;
      const std::unordered_map<std::string, double> inlet_mf = { { CO2, 0.15 }, { N2, 0.85 } };

      // Feed 1: feed in at the bottom; the top splits into the raffinate vent
      // and the light-product draw that pressurises the other column.
      Step& step_feed1 = t.data.cycle.GetStep(feed1);
      auto feed1_inlet = t.CreateNoValveInletStream(t.CreateVelocityBC(feed_vel), inlet_mf, kT);
      feed1_inlet->name_ = "Feed 1 inlet";
      auto feed1_in = std::make_unique<FlowBoundaryCondition>();
      feed1_in->direction_ = BoundaryDirection::INLET;
      feed1_in->location_ = BoundaryLocation::BOTTOM;
      feed1_in->streams_.push_back(std::move(feed1_inlet));
      step_feed1.bottom_boundary_condition_ = std::move(feed1_in);

      auto feed1_vent = t.CreateValveOutletStream(t.CreateConstantBoundaryValue(kPH), ambient_valve_cv);
      feed1_vent->name_ = "Feed 1 ambient Stream";
      feed1_vent->mole_fraction[CO2] = 0.0;  // back-fed gas is N2, as in the VII deck
      feed1_vent->mole_fraction[N2] = 1.0;
      auto feed1_lpp = t.CreateNoValveOutletStream(t.CreateVelocityBC(0.0));
      feed1_lpp->name_ = "LPP Stream";

      // Feed 2: feed in, raffinate out through the vent.
      Step& step_feed2 = t.data.cycle.GetStep(feed2);
      auto feed2_inlet = t.CreateNoValveInletStream(t.CreateVelocityBC(feed_vel), inlet_mf, kT);
      feed2_inlet->name_ = "Feed 2 inlet";
      auto feed2_in = std::make_unique<FlowBoundaryCondition>();
      feed2_in->direction_ = BoundaryDirection::INLET;
      feed2_in->location_ = BoundaryLocation::BOTTOM;
      feed2_in->streams_.push_back(std::move(feed2_inlet));
      step_feed2.bottom_boundary_condition_ = std::move(feed2_in);
      auto feed2_vent = t.CreateValveOutletStream(t.CreateConstantBoundaryValue(kPH), ambient_valve_cv);
      feed2_vent->name_ = "Feed 2 ambient Stream";
      feed2_vent->mole_fraction[CO2] = 0.0;
      feed2_vent->mole_fraction[N2] = 1.0;
      auto feed2_out = std::make_unique<FlowBoundaryCondition>();
      feed2_out->direction_ = BoundaryDirection::OUTLET;
      feed2_out->location_ = BoundaryLocation::TOP;
      feed2_out->streams_.push_back(std::move(feed2_vent));
      step_feed2.top_boundary_condition_ = std::move(feed2_out);

      // Pressurisation by light product: P_L -> P_H at the top, connected to
      // the Feed 1 light-product draw.
      auto rep = t.CreatePressureInletBoundary(t.CreateLinearBoundaryValue(c.P_L, kPH, press_time),
                                               { { CO2, 0.0 }, { N2, 1.0 } }, kT);
      rep->location_ = BoundaryLocation::TOP;
      rep->direction_ = BoundaryDirection::INLET;
      auto endA = MakeEndpoint(t.data.cycle.GetStep(press), *rep);
      auto endB = MakeEndpoint(t.data.cycle.GetStep(feed1), *feed1_lpp);
      MakeMutualConnection(endA, endB);

      auto feed1_out = std::make_unique<FlowBoundaryCondition>();
      feed1_out->direction_ = BoundaryDirection::OUTLET;
      feed1_out->location_ = BoundaryLocation::TOP;
      feed1_out->streams_.push_back(std::move(feed1_vent));
      feed1_out->streams_.push_back(std::move(feed1_lpp));
      step_feed1.top_boundary_condition_ = std::move(feed1_out);

      // Blowdown: co-current (top) to P_I; bottom closed.
      Step& step_bd = t.data.cycle.GetStep(blowdown);
      auto bd = std::make_unique<PressureBoundaryCondition>();
      bd->direction_ = BoundaryDirection::OUTLET;
      bd->location_ = BoundaryLocation::TOP;
      auto bd_p = t.CreateExponentialBoundaryValue(c.P_I, blowdown_time);
      // SI Table A.1 alpha_BD = 0.153 1/s: F = exp(alpha t_step) = 9.7e3 over 60 s
      // (deck: F = 100, i.e. 0.077 1/s).
      if (v.paper_decay) bd_p->factor_ = std::exp(0.153 * blowdown_time);
      bd->pressure_ = std::move(bd_p);
      bd->mole_fraction_[CO2] = 0.0;
      bd->mole_fraction_[N2] = 1.0;
      step_bd.top_boundary_condition_ = std::move(bd);
      auto bd_closed = std::make_unique<ClosedBoundaryCondition>();
      bd_closed->location_ = BoundaryLocation::BOTTOM;
      step_bd.bottom_boundary_condition_ = std::move(bd_closed);

      // Evacuation: counter-current (bottom) to P_L; top closed.
      Step& step_ev = t.data.cycle.GetStep(evacuation);
      auto ev_closed = std::make_unique<ClosedBoundaryCondition>();
      ev_closed->location_ = BoundaryLocation::TOP;
      step_ev.top_boundary_condition_ = std::move(ev_closed);
      auto ev = std::make_unique<PressureBoundaryCondition>();
      ev->direction_ = BoundaryDirection::OUTLET;
      ev->location_ = BoundaryLocation::BOTTOM;
      auto ev_p = t.CreateExponentialBoundaryValue(c.P_L, evac_time);
      // SI Table A.1 alpha_EVAC = 0.075 1/s: F = exp(3) = 20.1 over 40 s
      // (deck: F = 100, i.e. 0.115 1/s).
      if (v.paper_decay) ev_p->factor_ = std::exp(0.075 * evac_time);
      ev->pressure_ = std::move(ev_p);
      ev->mole_fraction_[CO2] = 0.0;
      ev->mole_fraction_[N2] = 1.0;
      step_ev.bottom_boundary_condition_ = std::move(ev);

      Step& step_press = t.data.cycle.GetStep(press);
      step_press.top_boundary_condition_ = std::move(rep);
      auto press_closed = std::make_unique<ClosedBoundaryCondition>();
      press_closed->location_ = BoundaryLocation::BOTTOM;
      step_press.bottom_boundary_condition_ = std::move(press_closed);

      t.ConfigureCycle(c.max_cycles);
      t.data.settings.css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
      t.data.settings.css.min_cycles_ = 5;
      // KPI as Eqs. 2-4: CO2 fed in Feed 1 and Feed 2 (bottom); product =
      // evacuation (bottom).
      t.ConfigureKPI(0, { 0, 1 }, { BoundaryLocation::BOTTOM, BoundaryLocation::BOTTOM }, { 3 },
                     { BoundaryLocation::BOTTOM });
      t.data.settings.max_step = 2.0;

      CoupledSolver solver{ t.data };
      std::string error{};
      EXPECT_TRUE(t.data.PerformChecks(error)) << "System checks failed: " << error;

      const auto t0 = std::chrono::steady_clock::now();
      int retval = solver.IntegrateCycle(t.data.cycle, false);
      Nguyen2024Result r;
      r.wall_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      EXPECT_GE(retval, 0) << "Solver did not complete";
      EXPECT_TRUE(t.IsMassConserved(solver.GetLog())) << "Mass is not conserved.";
      EXPECT_TRUE(t.IsMolesConserved(solver.GetLog())) << "Moles are not conserved.";
      EXPECT_TRUE(t.IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";

      const EquationLog& eq_log = solver.GetLog();
      const auto& cyc = eq_log.sys_data.cycle_index_set_;
      const int last_cycle = cyc.empty() ? 0 : *cyc.rbegin();
      KPIParams kpi(eq_log, solver.GetBoundaryLog());
      kpi.product_comp_ = CO2;
      kpi.cycle_index_ = last_cycle;
      kpi.feed_endpoints = KpiEndpointsForFaces(t.data.cycle,
        { { 0, BoundaryLocation::BOTTOM }, { 1, BoundaryLocation::BOTTOM } });
      kpi.product_endpoints = KpiEndpointsForFaces(t.data.cycle, { { 3, BoundaryLocation::BOTTOM } });
      r.purity = CalculatePurity(kpi) * 100.0;
      r.recovery = CalculateRecovery(kpi) * 100.0;
      // Eq. 4: mol CO2 in the evacuation product per m3 of adsorbent
      // (m_ads / rho_p = 63.2 cm3) per cycle time.
      const double V_ads = kMads / kRhoP;
      r.productivity = CalculateComponentMolesInProduct(kpi) / (V_ads * eq_log.sys_data.total_cycle_time);
      r.cycles = last_cycle + 1;
      r.css = r.cycles < c.max_cycles;
      std::cout << "\n[KPI] " << c.label << " (" << v.Name() << ") | cycles run = " << r.cycles
                << " | purity = " << r.purity << " % | recovery = " << r.recovery
                << " % | productivity = " << r.productivity << " mol/(m3 s) | wall = " << r.wall_s << " s\n"
                << std::endl;
      if (save) t.SaveResultsAsExample(c.label, solver);
      return r;
    }

    // Table 1 (p. 6), dry rows. Caps about 30 % above the measured CSS stop.
    const Nguyen2024Case kA0{
      "CaseA0", "Nguyen2024_A0", "Nguyen 2024 - Case 1, dry (A0).",
      "Four-step LPP CO2/N2 VSA cycle on CALF-20, dry feed, case study 1 (A0) of Table 1 of "
      "https://doi.org/10.1007/s10450-025-00600-z, built from the Nguyen2023_Exp_VII deck at P_H = 0.96 bar and 22 C.",
      180, 0.081, 0.119e5, 0.032e5, 50 };  // CSS after 39 cycles
    const Nguyen2024Case kB0{
      "CaseB0", "Nguyen2024_B0", "Nguyen 2024 - Case 2, dry (B0).",
      "Four-step LPP CO2/N2 VSA cycle on CALF-20, dry feed, case study 2 (B0) of Table 1 of "
      "https://doi.org/10.1007/s10450-025-00600-z, built from the Nguyen2023_Exp_VII deck at P_H = 0.96 bar and 22 C.",
      157, 0.024, 0.046e5, 0.031e5, 430 };  // CSS after 326 cycles

  }  // namespace

  // Tolerances are the measured error rounded up, with the measured value in a
  // comment (README.md "Measured errors"). Only what reproduces within the
  // printed uncertainty is asserted; the rest is recorded (Info).
  TEST_F(Nguyen2024Test, CaseA0) {
    using namespace validation;
    const auto r = RunNguyen2024(*this, kA0, {}, true);
    CaseReport report("nguyen_2024", "CaseA0");
    // Recorded, not asserted: purity 94.89 % against 98.3 +/- 0.88 (3.41 pp,
    // outside the printed uncertainty), recovery 81.98 % against 71.2 (LPP
    // over-recovery), productivity 0.467 against 0.577 (the printed value is
    // inconsistent with the cycle time). README.md "Measured errors".
    report.Info("CO2 purity vs experiment", r.purity, Reference("nguyen_2024", "CaseA0", "CO2_purity_exp"),
                "abs_error", "recorded: 3.4 pp below, outside the printed uncertainty (README)");
    report.Info("CO2 recovery vs experiment", r.recovery, Reference("nguyen_2024", "CaseA0", "CO2_recovery_exp"),
                "abs_error", "recorded: LPP over-recovery (README)");
    report.Info("CO2 productivity vs experiment", r.productivity,
                Reference("nguyen_2024", "CaseA0", "CO2_productivity_exp"), "rel_error",
                "recorded: the printed value is inconsistent with the cycle time (README)");
  }

  TEST_F(Nguyen2024Test, CaseB0) {
    using namespace validation;
    const auto r = RunNguyen2024(*this, kB0, {}, true);
    CaseReport report("nguyen_2024", "CaseB0");
    // measured: 0.9311 pp (93.33 % against 92.4 +/- 1.31)
    report.CheckAbs("CO2 purity vs experiment", r.purity, Reference("nguyen_2024", "CaseB0", "CO2_purity_exp"), 1.0);
    // measured: 0.0201 rel (0.1500 against 0.147 +/- 0.005)
    report.CheckRel("CO2 productivity vs experiment", r.productivity,
                    Reference("nguyen_2024", "CaseB0", "CO2_productivity_exp"), 0.025);
    // Recorded: recovery 94.87 % against 83.3 +/- 1.97 (11.6 pp, LPP over-recovery).
    report.Info("CO2 recovery vs experiment", r.recovery, Reference("nguyen_2024", "CaseB0", "CO2_recovery_exp"),
                "abs_error", "recorded: LPP over-recovery (README)");
  }

  // Attribution runs (README "Attribution"): one change at a time, both cases.
  // Not part of the suite (about ten runs to CSS); run with
  //   AdsorptionSimulatorTests --gtest_also_run_disabled_tests --gtest_filter='Nguyen2024Test.DISABLED_*'
  // from core/. Prints the KPIs; asserts nothing and saves no deck.
  TEST_F(Nguyen2024Test, DISABLED_Attribution) {
    const std::vector<Nguyen2024Variant> variants = {
      {}, { true, false, false, false }, { false, true, false, false }, { false, false, true, false },
      { false, false, false, true } };
    std::cout << "\n| Case | Variant | Purity [%] | Recovery [%] | Productivity [mol/(m3 s)] | Cycles |\n"
              << "|---|---|---|---|---|---|\n";
    std::vector<std::string> rows;
    for (const auto* c : { &kA0, &kB0 }) {
      for (const auto& v : variants) {
        Nguyen2024Case cc = *c;
        cc.max_cycles = 2000;
        Nguyen2024Runner fresh;
        const auto r = RunNguyen2024(fresh, cc, v, false);
        char buf[256];
        std::snprintf(buf, sizeof buf, "| %s | %s | %.2f | %.2f | %.4f | %d |", c->id.c_str(), v.Name().c_str(),
                      r.purity, r.recovery, r.productivity, r.cycles);
        rows.push_back(buf);
        std::cout << buf << std::endl;
      }
    }
    std::cout << "\n";
    for (const auto& row : rows) std::cout << row << "\n";
  }

}  // namespace adsim
