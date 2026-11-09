/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Analytical validation of a FULL PSA CYCLE: the Skarstrom cycle at
// equilibrium-theory conditions against the Knaebel-Hill purity and recovery.
//
// test_blowdown_composition_path pins one pressure-change step. This test runs
// the whole four-step cycle and checks the two numbers a PSA designer reads
// off it, the light-product purity and recovery, against the closed form of
// the equilibrium theory of Chan, Hill & Wong (1981, Chem. Eng. Sci. 36, 243)
// as generalised to two adsorbing components by Knaebel & Hill (1985, Chem.
// Eng. Sci. 40, 2351); also Ruthven, Farooq & Knaebel, "Pressure Swing
// Adsorption" (1994), ch. 4. The derivation below is written out in full, so
// the reference does not depend on reading the papers.
//
// Conditions: isothermal, linear isotherms q_i = K_i p_i, binary feed of a
// heavy component A (mole fraction y_F) and a light one B, local equilibrium,
// no axial dispersion, no pressure drop. Per unit bed volume a component holds
// a_i P y_i, with
//
//     a_i = eps / (R T) + rho_b K_i,   beta_i = (eps / (R T)) / a_i,
//     beta = beta_A / beta_B = a_B / a_A < 1.
//
// The cycle (feed end = bottom, product end = top), started from the bed full
// of pure B at P_L, which is the state the purge leaves:
//
//   1. Pressurisation with light product, P_L -> P_H, from the top (bottom
//      closed). The bed holds pure B throughout: n_press = V a_B (P_H - P_L).
//   2. High-pressure feed at P_H. The feed meets pure B at a shock whose jump
//      conditions are, for A, w a_A y_F = (eps/(R T)) v_F y_F, i.e. the shock
//      moves at w = beta_A v_F (v_F interstitial), and, for B, fix the product
//      velocity v_2 = v_F (1 - (1 - beta) y_F). The step ends as the shock
//      reaches the top: n_feed = V a_A P_H, and the effluent, pure B, is
//      n_out = V a_A P_H (1 - (1 - beta) y_F) (also the B balance of the step:
//      fed + held at the start - held at the end).
//   3. Countercurrent blowdown P_H -> P_L from the bottom (top closed). The
//      bed is uniform at y_F, so it follows the blowdown path of
//      test_blowdown_composition_path; all of the gas leaves as heavy product.
//   4. Countercurrent purge with light product at P_L from the top. Pure B
//      displaces the A-rich gas in a simple wave (desorption of the heavy
//      component speeds the gas up behind it), whose slowest characteristic,
//      y = 0, moves at beta_A v_P. It sweeps the bed clean once
//      n_purge = gamma V a_A P_L with gamma >= 1 (gamma = 1: the minimum purge;
//      gamma is the purge-to-feed ratio of actual volumes, purge at P_L over
//      feed at P_H). The bed is then pure B at P_L, the state of step 1, so the
//      cycle is at its cyclic steady state from the first cycle.
//
// The light product (feed effluent, less the gas returned for pressurisation
// and purge) is pure B,
//
//     n_prod = V a_A [ (1 - beta) (P_H (1 - y_F) - P_L) - (gamma - 1) P_L ],
//
// so, with the pressure ratio Pi = P_H / P_L,
//
//     purity   = 1,
//     recovery = (1 - beta) (1 - 1 / (Pi (1 - y_F))) - (gamma - 1) / (Pi (1 - y_F)),
//
// the Knaebel-Hill result for pressurisation with product (the first term at
// the minimum purge gamma = 1). All of A leaves in the heavy product
// (blowdown + purge effluent), whose mean A fraction is
//
//     y_W = P_H y_F / ( P_H (beta + (1 - beta) y_F) + P_L (gamma - beta) ).
//
// The simulated cycle is the theory's open single-bed cycle: the gas used for
// pressurisation and purge is supplied as pure B, which is what the theory's
// product is, and the light product counted as the theory counts it (feed
// effluent of B less the B supplied to steps 1 and 4).
//
// What the simulator adds to the theory: a finite LDF rate and the numerical
// dispersion of the finite-volume grid smear the feed shock, so a little A
// leaves with the light product as the shock reaches the top, and the tail of
// the purge wave leaves a little A near the feed end. Both are the width of
// the smeared front, which shrinks as 1/k (the constant-pattern width of a
// self-sharpening front on an LDF bed) and as the cell size. The fast test
// pins the gap on one grid; the slow test refines the grid and the LDF rate
// together and checks that the gap shrinks and extrapolates to the theory.
// Pressure drop: about 20 Pa at the feed velocity with 5 mm particles, 4e-5 of
// P_H and 2e-4 of P_L. Temperature: the isothermal mode holds T exactly.
//
// GUARDS: the cycle separates (recovery about 0.5, well away from 0 and 1;
// heavy product enriched from 0.3 to 0.48 and 0.41), beta < 0.3 so a model without
// the adsorbed phase (beta = 1, recovery 0) fails, and the two purge ratios
// differ in recovery by more than ten times the tolerance, so the test sees
// the purge term.

#include "validation_helpers.h"

#include "adsim/solver/kpis.h"

#include <algorithm>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    constexpr double kT = 298.0;
    constexpr double kEps = 0.4;
    constexpr double kRhoS = 1000.0;
    constexpr double kKA = 1.6e-6;  // [mol/(kg Pa)] heavy: CO2
    constexpr double kKB = 2.0e-7;  // [mol/(kg Pa)] light: N2
    constexpr double kL = 1.0, kD = 0.1;
    constexpr double kPH = 5.0e5, kPL = 1.0e5;  // [Pa], Pi = 5
    constexpr double kYF = 0.3;                  // feed mole fraction of A
    constexpr double kTFeed = 100.0;             // [s] shock transit time
    constexpr double kTPurge = 100.0;            // [s]
    constexpr double kTRamp = 50.0;              // [s] pressure ramps
    constexpr double kTHold = 25.0;              // [s] hold at the end pressure
    const char* const kHeavy = "CarbonDioxide";
    const char* const kLight = "Nitrogen";

    double Volume() { return 0.25 * M_PI * kD * kD * kL; }
    double RhoB() { return (1.0 - kEps) * kRhoS; }
    double Capacity(double K) { return kEps / (kR * kT) + RhoB() * K; }  // a_i
    double Beta() { return Capacity(kKB) / Capacity(kKA); }

    // Knaebel-Hill light-product recovery for pressurisation with product.
    double TheoryRecovery(double gamma) {
      const double pi = kPH / kPL, b = Beta();
      return (1.0 - b) * (1.0 - 1.0 / (pi * (1.0 - kYF))) - (gamma - 1.0) / (pi * (1.0 - kYF));
    }
    // Mean A fraction of the heavy product (blowdown + purge effluent).
    double TheoryHeavyFraction(double gamma) {
      const double b = Beta();
      return kPH * kYF / (kPH * (b + (1.0 - b) * kYF) + kPL * (gamma - b));
    }

    std::unique_ptr<LinearBoundaryValue> RampTo(double to) {
      auto bv = std::make_unique<LinearBoundaryValue>();
      bv->end_value_ = to;  // starts from the bed pressure
      bv->ramp_time_ = kTRamp;
      bv->use_step_time_ = false;
      return bv;
    }
    std::unique_ptr<ConstantBoundaryValue> Constant(double v) {
      auto bv = std::make_unique<ConstantBoundaryValue>();
      bv->constant_value_ = v;
      return bv;
    }
    std::unique_ptr<ClosedBoundaryCondition> Closed(BoundaryLocation where) {
      auto bc = std::make_unique<ClosedBoundaryCondition>();
      bc->location_ = where;
      return bc;
    }
    template <class Map>
    void SetComposition(Map& y, double y_heavy) {
      y[kHeavy] = y_heavy;
      y[kLight] = 1.0 - y_heavy;
    }
    std::unique_ptr<PressureBoundaryCondition> PressureEnd(BoundaryLocation where, BoundaryDirection dir,
                                                           std::unique_ptr<IBoundaryValue> p) {
      auto bc = std::make_unique<PressureBoundaryCondition>();
      bc->location_ = where;
      bc->direction_ = dir;
      bc->pressure_ = std::move(p);
      bc->temperature_ = kT;
      SetComposition(bc->mole_fraction_, 0.0);  // pure B where it is an inlet
      return bc;
    }
    std::unique_ptr<FlowBoundaryCondition> MolarFlowInlet(BoundaryLocation where, const std::string& name,
                                                          double molar_flow, double y_heavy) {
      auto stream = std::make_unique<BoundaryStream>();
      stream->name_ = name;
      auto flow = std::make_unique<MolarFlowBoundaryValue>();
      static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = molar_flow;
      stream->boundary_flow_ = std::move(flow);
      SetComposition(stream->mole_fraction, y_heavy);
      stream->temperature_ = kT;
      auto bc = std::make_unique<FlowBoundaryCondition>();
      bc->direction_ = BoundaryDirection::INLET;
      bc->location_ = where;
      bc->streams_.push_back(std::move(stream));
      return bc;
    }

    struct Config {
      int cells = 50;
      double ldf = 20.0;  // [1/s]
      double gamma = 1.0;
      int cycles = 3;
    };

    void BuildCycle(UserData& data, const Config& cfg) {
      data.bed.wall_.thickness_ = 0.005;
      data.fluid.AddComponent(kLight);
      data.fluid.AddComponent(kHeavy);
      auto ads = data.adsorbent_directory.Add("Packing");
      auto lid = data.bed.layers_.Add("Packing");
      auto& layer = data.bed.layers_.Get(lid);
      layer.momentum_model_ = MomentumModel::KarmenCozeny;
      layer.domain_.length = kL;
      layer.domain_.n_cells = cfg.cells;
      layer.eb_ = kEps;
      layer.particle_diameter_ = 0.005;  // coarse: negligible pressure drop
      layer.composition_[ads] = 1.0;
      std::get<CylinderBed>(data.bed.shape_).diameter = kD;
      auto& adsorbent = data.adsorbent_directory.Get(ads);
      adsorbent.ep_ = 0.0;
      adsorbent.solid_density_ = kRhoS;
      adsorbent.heat_capacity_ = 900.0;

      auto iso = std::make_unique<HenryIsotherm>();
      iso->params_[kHeavy].k = kKA;
      iso->params_[kHeavy].Q = 0.0;
      iso->params_[kLight].k = kKB;
      iso->params_[kLight].Q = 0.0;
      adsorbent.adsorption_source_.AddIsotherm({kHeavy, kLight}, std::move(iso));
      for (const char* c : {kHeavy, kLight}) {
        auto ldf = std::make_unique<LinearDrivingForce>();
        auto mtc = std::make_unique<MassTransferCoefficientConstant>();
        mtc->constant_mtc_ = cfg.ldf;
        ldf->mtc_model_ = ToResistances(std::move(mtc));  // generic model, as the example decks store it
        adsorbent.adsorption_source_.mass_transfer_models_[c] = std::move(ldf);
        adsorbent.adsorption_source_.heat_of_adsorption_[c] = 0.0;
      }
      data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;

      // The purge's end state: pure B at P_L.
      data.ic.pressure = kPL;
      data.ic.temperature = kT;
      data.ic.wall_temperature = kT;
      data.ic.velocity = 0.0;
      SetComposition(data.ic.mole_fraction, 0.0);

      const double feed_flow = Volume() * Capacity(kKA) * kPH / kTFeed;               // n_feed / t
      const double purge_flow = cfg.gamma * Volume() * Capacity(kKA) * kPL / kTPurge;  // n_purge / t

      data.cycle.AddStep("Pressurisation");
      {
        Step& s = data.cycle.GetStep("Pressurisation");
        s.step_time_ = kTRamp + kTHold;
        s.bottom_boundary_condition_ = Closed(BoundaryLocation::BOTTOM);
        s.top_boundary_condition_ = PressureEnd(BoundaryLocation::TOP, BoundaryDirection::INLET, RampTo(kPH));
      }
      data.cycle.AddStep("Feed");
      {
        Step& s = data.cycle.GetStep("Feed");
        s.step_time_ = kTFeed;
        s.bottom_boundary_condition_ = MolarFlowInlet(BoundaryLocation::BOTTOM, "Feed", feed_flow, kYF);
        s.top_boundary_condition_ = PressureEnd(BoundaryLocation::TOP, BoundaryDirection::OUTLET, Constant(kPH));
      }
      data.cycle.AddStep("Blowdown");
      {
        Step& s = data.cycle.GetStep("Blowdown");
        s.step_time_ = kTRamp + kTHold;
        s.bottom_boundary_condition_ = PressureEnd(BoundaryLocation::BOTTOM, BoundaryDirection::OUTLET, RampTo(kPL));
        s.top_boundary_condition_ = Closed(BoundaryLocation::TOP);
      }
      data.cycle.AddStep("Purge");
      {
        Step& s = data.cycle.GetStep("Purge");
        s.step_time_ = kTPurge;
        s.top_boundary_condition_ = MolarFlowInlet(BoundaryLocation::TOP, "Purge", purge_flow, 0.0);
        s.bottom_boundary_condition_ = PressureEnd(BoundaryLocation::BOTTOM, BoundaryDirection::OUTLET, Constant(kPL));
      }

      data.cycle.number_of_cycles_ = cfg.cycles;
      data.settings.css.method_ = CSSMethod::NONE;
      // The KPIs integrate the logged stream flows in time; a capped step keeps
      // the quadrature error well below the tolerances.
      data.settings.reltol = 1e-8;
      data.settings.max_step = 0.5;
    }

    // Moles of `component` crossing one end of one step, positive into the bed.
    double MolesIn(const BoundaryLog& blog, const Step& step, int cycle, int step_index,
                   BoundaryLocation where, const std::string& component) {
      const IBoundaryCondition* bc = where == BoundaryLocation::BOTTOM ? step.bottom_boundary_condition_.get()
                                                                         : step.top_boundary_condition_.get();
      // The logged molar flow has the velocity sign (positive bottom -> top).
      const double sign = where == BoundaryLocation::BOTTOM ? 1.0 : -1.0;
      std::vector<std::string> names;
      if (auto fbc = dynamic_cast<const FlowBoundaryCondition*>(bc)) {
        for (const auto& stream : fbc->streams_) names.push_back(stream->name_);
      } else if (dynamic_cast<const PressureBoundaryCondition*>(bc)) {
        names.push_back(step.step_name_ + " Pressure Boundary " + BoundaryLocationList[(int)where]);
      }
      double moles = 0.0;
      for (const auto& name : names) {
        const auto d = blog.GetComponentVariableData(cycle, step_index, name,
                                                     BoundaryComponentVariableType::MOLAR_FLOW, component);
        moles += sign * IntegrateVariableInTime(d.variable_timeseries, d.time);
      }
      return moles;
    }

    struct CycleKpis {
      double purity = 0.0;          // B in the feed effluent
      double recovery = 0.0;        // net light product B / B fed
      double heavy_fraction = 0.0;  // A in blowdown + purge effluent
      double a_balance = 0.0;       // (A out - A in) / A in over the cycle
    };

    CycleKpis KpisOfCycle(const UserData& data, const BoundaryLog& blog, int cycle) {
      const auto& steps = data.cycle.step_list_;
      auto in = [&](int s, BoundaryLocation where, const char* c) { return MolesIn(blog, steps[s], cycle, s, where, c); };
      using BL = BoundaryLocation;
      const double fed_A = in(1, BL::BOTTOM, kHeavy), fed_B = in(1, BL::BOTTOM, kLight);
      const double eff_A = -in(1, BL::TOP, kHeavy), eff_B = -in(1, BL::TOP, kLight);
      const double press_B = in(0, BL::TOP, kLight) + in(0, BL::TOP, kHeavy);
      const double purge_B = in(3, BL::TOP, kLight) + in(3, BL::TOP, kHeavy);
      const double w_A = -in(2, BL::BOTTOM, kHeavy) - in(3, BL::BOTTOM, kHeavy);
      const double w_B = -in(2, BL::BOTTOM, kLight) - in(3, BL::BOTTOM, kLight);
      CycleKpis k;
      k.purity = eff_B / (eff_A + eff_B);
      k.recovery = (eff_B - press_B - purge_B) / fed_B;
      k.heavy_fraction = w_A / (w_A + w_B);
      k.a_balance = (eff_A + w_A - fed_A) / fed_A;
      return k;
    }

    std::vector<CycleKpis> RunCycle(const Config& cfg) {
      UserData data;
      BuildCycle(data, cfg);
      CoupledSolver solver{data};
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
      std::vector<CycleKpis> out;
      for (int c = 0; c < cfg.cycles; ++c) out.push_back(KpisOfCycle(data, solver.GetBoundaryLog(), c));
      return out;
    }

    void Report(const char* tag, const Config& cfg, const CycleKpis& k) {
      std::cout << "[SkarstromET] " << tag << " cells=" << cfg.cells << " k=" << cfg.ldf
                << " gamma=" << cfg.gamma << ": purity " << k.purity << " (theory 1), recovery "
                << k.recovery << " (theory " << TheoryRecovery(cfg.gamma) << "), heavy y "
                << k.heavy_fraction << " (theory " << TheoryHeavyFraction(cfg.gamma) << "), A balance "
                << k.a_balance << "\n";
    }

  } // namespace

  TEST(SkarstromCycleEquilibriumTheory, PurityAndRecoveryMatchKnaebelHill) {
    // GUARD: a strongly selective pair, so the separation is far from trivial.
    ASSERT_LT(Beta(), 0.3);

    std::vector<CycleKpis> last;
    for (double gamma : {1.0, 1.5}) {
      Config cfg;
      cfg.gamma = gamma;
      const auto kpis = RunCycle(cfg);
      ASSERT_EQ(kpis.size(), static_cast<size_t>(cfg.cycles));
      const CycleKpis& k = kpis.back();
      Report("fast", cfg, k);
      last.push_back(k);

      // GUARD: the cycle separates.
      ASSERT_GT(k.recovery, 0.2);
      ASSERT_GT(k.heavy_fraction, kYF + 0.1);
      // At the theory's cyclic steady state from the first cycle: the KPIs
      // barely move between the last two cycles.
      EXPECT_NEAR(k.recovery, kpis[kpis.size() - 2].recovery, 1e-3);
      // Every mole of A fed leaves the bed over the cycle.
      EXPECT_NEAR(k.a_balance, 0.0, 2e-3);

      // Tolerances: regression bounds over the measured gap on this grid
      // (50 cells, k = 20 1/s): purity 4.6e-3, recovery 1.3e-3, heavy product
      // 3.7e-3. The gap is the smeared feed shock leaking A as it reaches the
      // top; SkarstromCycleEquilibriumTheory.GapVanishesUnderRefinement shows
      // it falling at first order to the theory.
      EXPECT_NEAR(k.purity, 1.0, 6e-3) << "Light product not pure B at gamma = " << gamma;
      EXPECT_NEAR(k.recovery, TheoryRecovery(gamma), 2e-3) << "Recovery at gamma = " << gamma;
      EXPECT_NEAR(k.heavy_fraction, TheoryHeavyFraction(gamma), 5e-3) << "Heavy product at gamma = " << gamma;
    }
    // The purge term: the drop in recovery from gamma = 1 to 1.5 (theory
    // 0.142857, measured 0.142793) is resolved far better than the recovery
    // itself, because the smearing is the same at both purge ratios.
    EXPECT_NEAR(last[0].recovery - last[1].recovery, TheoryRecovery(1.0) - TheoryRecovery(1.5), 5e-4);
  }

  // The strong form: refine the grid and the LDF rate together (cells x2,
  // k x2) at the minimum purge. The gap to the theory halves at each level
  // (first order: the shock is captured at first order, and the LDF part of
  // the gap, 1/k, is a few per cent of it), and the Richardson extrapolation
  // of the two finest levels lands on the theory. Measured (100/40, 200/80,
  // 400/160): purity gap 2.30e-3, 1.14e-3, 5.67e-4; recovery gap 6.35e-4,
  // 3.15e-4, 1.56e-4; extrapolated recovery 0.535030 against 0.535027.
  // About 16 s, so in the slow tier.
  TEST(SkarstromCycleEquilibriumTheory, GapVanishesUnderRefinement) {
    struct Level { int cells; double ldf; CycleKpis k; };
    std::vector<Level> levels{{100, 40.0, {}}, {200, 80.0, {}}, {400, 160.0, {}}};
    for (auto& lv : levels) {
      Config cfg;
      cfg.cells = lv.cells;
      cfg.ldf = lv.ldf;
      cfg.cycles = 2;
      const auto kpis = RunCycle(cfg);
      ASSERT_EQ(kpis.size(), 2u);
      lv.k = kpis.back();
      Report("refine", cfg, lv.k);
    }
    auto purity_gap = [](const CycleKpis& k) { return 1.0 - k.purity; };
    auto recovery_gap = [](const CycleKpis& k) { return TheoryRecovery(1.0) - k.recovery; };
    auto heavy_gap = [](const CycleKpis& k) { return TheoryHeavyFraction(1.0) - k.heavy_fraction; };
    for (size_t i = 0; i + 1 < levels.size(); ++i) {
      const CycleKpis& c = levels[i].k;
      const CycleKpis& f = levels[i + 1].k;
      // GUARD: the gap is there to measure (positive at every level).
      ASSERT_GT(purity_gap(f), 0.0);
      ASSERT_GT(recovery_gap(f), 0.0);
      EXPECT_NEAR(purity_gap(c) / purity_gap(f), 2.0, 0.2) << "purity, level " << i;
      EXPECT_NEAR(recovery_gap(c) / recovery_gap(f), 2.0, 0.2) << "recovery, level " << i;
      EXPECT_NEAR(heavy_gap(c) / heavy_gap(f), 2.0, 0.2) << "heavy product, level " << i;
    }
    // First-order Richardson extrapolation of the two finest levels.
    const CycleKpis& c = levels[1].k;
    const CycleKpis& f = levels[2].k;
    const double purity_x = 2.0 * f.purity - c.purity;
    const double recovery_x = 2.0 * f.recovery - c.recovery;
    const double heavy_x = 2.0 * f.heavy_fraction - c.heavy_fraction;
    std::cout << "[SkarstromET] extrapolated purity " << purity_x << ", recovery " << recovery_x
              << " (theory " << TheoryRecovery(1.0) << "), heavy y " << heavy_x << " (theory "
              << TheoryHeavyFraction(1.0) << ")\n";
    // Tolerances: a fifth of the finest level's gap (measured residuals 9e-6,
    // 3e-6 and 8e-6, against finest gaps of 5.7e-4, 1.6e-4 and 4.5e-4).
    EXPECT_NEAR(purity_x, 1.0, 0.2 * purity_gap(f));
    EXPECT_NEAR(recovery_x, TheoryRecovery(1.0), 0.2 * recovery_gap(f));
    EXPECT_NEAR(heavy_x, TheoryHeavyFraction(1.0), 0.2 * heavy_gap(f));
  }

} // namespace adsim
