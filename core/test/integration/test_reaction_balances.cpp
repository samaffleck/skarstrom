/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// Reactions in the balances (objects/reaction.h, UpdateReactions in
// equations/adsorbent_equations.h): a packed catalyst bed (inert isotherms)
// built from scratch, one feed step through it.
//
//  (a) a zero rate reproduces the deck without the reaction bit for bit, and
//      a deck without reactions registers no reaction fields;
//  (b) isothermal first-order A -> B in plug flow (dispersion off, upwind):
//      the conversion against X = 1 - exp(-Da) and against the upwind
//      scheme's own steady state 1 - (1 + Da/N)^-N, at two residence times;
//  (c) A -> B + C at constant pressure raises the outlet velocity by
//      1 + y_A,in X (sum(nu) != 0 enters continuity);
//  (d) adiabatic temperature rise (-dH) X y_A,in / cp in LTE and LTNE, with
//      the heat on the solid in LTNE, and the energy closure
//      integral(Q dV) = F cp dT at steady state;
//  (e) element and component closure of a reactive step from the logs, with
//      the generation the log integrates;
//  (f) both mass-balance CSS methods converge on a short reactive cycle;
//  (g) the Thiele effectiveness factor against the closed form of a
//      first-order rate, irreversible and reversible, and in a run;
//  (h) inert loadings are not unknowns: the same run with Henry (k = 0)
//      instead of Inert isotherms gives the same results with more unknowns.
//
// Fast tier: 40-100 cells, one or a few short steps.

#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/manifests.h"
#include "adsim/objects/reaction.h"
#include "adsim/simulation.h"
#include "adsim/solver/jacobian.h"
#include "adsim/solver/kpis.h"

namespace adsim {
namespace {

const std::string N2 = "Nitrogen", BUT = "Butane", IBUT = "Isobutane", C2H6 = "Ethane", C2H4 = "Ethylene",
                  H2 = "Hydrogen";

// Full precision text (std::to_string keeps 6 decimals).
std::string Num(double v) {
  char buf[40];
  std::snprintf(buf, sizeof(buf), "%.17g", v);
  return buf;
}

struct Options {
  std::vector<std::string> components{N2, BUT, IBUT};  // the carrier first
  std::map<std::string, double> nu{};                  // empty: no reaction
  std::map<std::string, double> y_feed{{N2, 0.9}, {BUT, 0.1}};
  std::string kind = "power_law";
  std::string basis = "concentration";
  double k = 0.0;                // pre-exponential factor (first order: m3/(kg s))
  double velocity = 0.1;         // [m/s]
  double length = 0.2;           // [m]
  int cells = 100;
  double step_time = 20.0;       // [s]
  bool isothermal = true;
  std::string thermal = "lte";
  double heat_of_reaction = 0.0;  // [J/mol]
  double temperature = 400.0;     // [K]
  double solid_density = 1000.0;  // [kg/m3]
  double solid_heat_capacity = 900.0;  // [J/(kg K)]
  double eb = 0.4;
  double particle_diameter = 5e-3;  // [m]
  std::string isotherm = "inert";   // of every component on the catalyst
};

struct Reactor {
  Simulation sim;
  AdsorbentID cat = 0;
  size_t layer = 0, feed = 0, feed_in = 0, outlet = 0;
  ReactionID rxn = 0;
  Options o;

  // rho_b = (1 - eb) rho_s [kg/m3]
  double BulkDensity() const { return (1.0 - o.eb) * o.solid_density; }
};

void Build(Reactor& r, const Options& o) {
  r.o = o;
  UserData& d = r.sim.Data();
  for (const auto& c : o.components) build::AddComponent(d, c);
  r.cat = build::AddAdsorbent(d, "Catalyst");  // inert isotherms
  if (o.isotherm != "inert") build::SetIsotherm(d, r.cat, o.isotherm, o.components);
  build::SetBedShape(d, "cylinder");
  build::SetWallMaterial(d, "ss316l");
  build::SetWallBoundary(d, "zero_flux");
  build::SetThermalEquilibrium(d, o.thermal);
  if (o.thermal == "ltne") build::SetParticleHeatTransfer(d, "constant");
  r.layer = build::AddLayer(d, "Bed");
  build::SetLayerComposition(d, r.layer, {{r.cat, 1.0}});
  build::SetMomentum(d, r.layer, "kozeny_carman");
  build::SetDispersion(d, r.layer, "constant");
  build::SetHeatDispersion(d, r.layer, "constant");
  build::SetInternalHeatTransfer(d, r.layer, "constant");
  build::SetExternalHeatTransfer(d, r.layer, "constant");

  r.feed = build::AddStep(d, "Feed");
  build::SetFlowBoundary(d, r.feed, "bottom", "inlet");
  r.feed_in = build::AddStream(d, r.feed, "bottom", "Feed inlet", "velocity", "constant");
  r.outlet = build::SetPressureBoundary(d, r.feed, "top", "outlet", "constant");
  if (!o.nu.empty()) {
    r.rxn = build::AddReaction(d, r.cat, "R1", o.nu, o.kind);
    build::SetReactionBasis(d, r.cat, r.rxn, o.basis);
  }
  r.sim.RefreshParameters();

  std::vector<std::pair<std::string, std::string>> values = {
      {"bed.diameter", "0.05"},
      {"adsorbent[Catalyst].solid_density", Num(o.solid_density)},
      {"adsorbent[Catalyst].heat_capacity", Num(o.solid_heat_capacity)},
      {"adsorbent[Catalyst].particle_void_fraction", "0.4"},
      {"layer[Bed].length", Num(o.length)},
      {"layer[Bed].cells", Num(o.cells)},
      {"layer[Bed].bed_void_fraction", Num(o.eb)},
      {"layer[Bed].particle_diameter", Num(o.particle_diameter)},
      {"initial.pressure", "1e5"},
      {"initial.temperature", Num(o.temperature)},
      {"initial.y[" + o.components.front() + "]", "1"},
      {"step[Feed].time", Num(o.step_time)},
      {"step[Feed].bottom[Feed inlet].velocity", Num(o.velocity)},
      {"step[Feed].bottom[Feed inlet].temperature", Num(o.temperature)},
      {"step[Feed].top.pressure", "1e5"},
  };
  for (const auto& [c, y] : o.y_feed) values.push_back({"step[Feed].bottom[Feed inlet].y[" + c + "]", Num(y)});
  if (!o.nu.empty()) {
    values.push_back({"adsorbent[Catalyst].reaction[R1].rate.pre_exp_factor", Num(o.k)});
    values.push_back({"adsorbent[Catalyst].reaction[R1].heat_of_reaction", Num(o.heat_of_reaction)});
  }
  for (const auto& [path, text] : values) r.sim.SetParameterText(path, text);
  if (!o.nu.empty()) r.sim.SetParameterText("adsorbent[Catalyst].reaction[R1].heat_mode", "constant");

  d.settings.advection_scheme_ = AdvectionScheme::Upwind;
  d.settings.css.method_ = CSSMethod::NONE;
  d.cycle.number_of_cycles_ = 1;
  d.bed.equations_.energy_balance_mode_ = o.isothermal ? EnergyBalanceMode::ISOTHERMAL : EnergyBalanceMode::NON_ISOTHERMAL;
}

RunResult RunReactor(Reactor& r, int cycles = 1) {
  RunOptions options;
  options.cycles = cycles;
  options.on_message = [](LogLevel, const std::string&) {};
  return r.sim.Run(options);
}

// Last logged row of a bulk / component variable at (cycle, step).
std::vector<double> LastRow(const Simulation& sim, VariableType t, int cycle = 0, int step = 0) {
  const auto e = sim.Log().GetVariable(cycle, step, t);
  return e.values.empty() ? std::vector<double>{} : e.values.back();
}
std::vector<double> LastRow(const Simulation& sim, ComponentVariableType t, const std::string& c, int cycle = 0,
                            int step = 0) {
  const auto e = sim.Log().GetComponentVariable(cycle, step, t, c);
  return e.values.empty() ? std::vector<double>{} : e.values.back();
}

// Outlet values: the last bed cell (the ghost mirrors it at an outlet).
double Outlet(const std::vector<double>& row) { return row.at(row.size() - 2); }

double MolarCp(const Simulation& sim, size_t cell) {
  const auto cp = LastRow(sim, VariableType::HEAT_CAPACITY);
  const auto mm = LastRow(sim, VariableType::MOLAR_MASS);
  return cp.at(cell) * mm.at(cell);  // [J/(mol K)]
}

// ---- (a) -----------------------------------------------------------------------

TEST(ReactionBalances, ZeroRateIsBitIdenticalAndNoReactionRegistersNothing) {
  Options o;
  o.cells = 40;
  o.step_time = 5.0;
  Reactor plain;
  Build(plain, o);
  ASSERT_TRUE(RunReactor(plain).Ok());
  EXPECT_FALSE(plain.sim.Data().state.HasReactions());
  EXPECT_TRUE(plain.sim.Data().state.component_reaction_source.empty());
  EXPECT_TRUE(plain.sim.Data().state.reaction_fields.empty());
  EXPECT_TRUE(plain.sim.Log().GetVariable(0, 0, VariableType::REACTION_HEAT).values.empty());
  EXPECT_TRUE(plain.sim.Log().GetLog(0, 0).reaction_generation.empty());

  o.nu = {{BUT, -1}, {IBUT, 1}};
  o.k = 0.0;
  o.heat_of_reaction = -5e4;
  Reactor zero;
  Build(zero, o);
  ASSERT_TRUE(RunReactor(zero).Ok());
  EXPECT_TRUE(zero.sim.Data().state.HasReactions());
  EXPECT_EQ(zero.sim.Data().state.reaction_fields.size(), 1u);

  const auto& a = plain.sim.Log().GetLog(0, 0);
  const auto& b = zero.sim.Log().GetLog(0, 0);
  ASSERT_EQ(a.step_data.times, b.step_data.times);
  for (auto t : {VariableType::PRESSURE, VariableType::VELOCITY, VariableType::TEMPERATURE})
    EXPECT_EQ(a.variables.at(t).values, b.variables.at(t).values);
  for (const auto& c : o.components)
    EXPECT_EQ(a.component_variables.at(c).at(ComponentVariableType::MOLE_FRACTION).values,
              b.component_variables.at(c).at(ComponentVariableType::MOLE_FRACTION).values)
        << c;
  // The reaction fields are logged, and zero.
  const auto Q = LastRow(zero.sim, VariableType::REACTION_HEAT);
  ASSERT_FALSE(Q.empty());
  for (double q : Q) EXPECT_EQ(q, 0.0);
  EXPECT_FALSE(zero.sim.Log().GetReactionRate(0, 0, zero.cat, zero.rxn).values.empty());
  EXPECT_EQ(zero.sim.Log().ReactionGeneration(0, 0, BUT), 0.0);
}

// ---- (b) -----------------------------------------------------------------------

// Steady outlet conversion of the plug-flow reactor with Damkoehler number
// Da = rho_b k L / u.
double FirstOrderConversion(double u, double k_scale = 1.0, Reactor* out = nullptr) {
  Options o;
  o.nu = {{BUT, -1}, {IBUT, 1}};
  o.velocity = u;
  // Da = 1 at u = 0.1 m/s.
  o.k = k_scale * 0.1 / ((1.0 - o.eb) * o.solid_density * o.length);
  Reactor local;
  Reactor& r = out ? *out : local;
  Build(r, o);
  const RunResult run = RunReactor(r);
  EXPECT_TRUE(run.Ok()) << run.message;
  const double y_out = Outlet(LastRow(r.sim, ComponentVariableType::MOLE_FRACTION, BUT));
  return 1.0 - y_out / o.y_feed.at(BUT);  // sum(nu) = 0: y_A / y_A,in = 1 - X
}

TEST(ReactionBalances, FirstOrderPlugFlowConversion) {
  // Kozeny-Carman with 5 mm particles: a pressure drop of a few Pa over the
  // bed, so u, C and the residence time are uniform to ~1e-4.
  const int N = 100;
  for (double u : {0.1, 0.2}) {
    SCOPED_TRACE(u);
    Reactor r;
    const double X = FirstOrderConversion(u, 1.0, &r);
    const double Da = 0.1 / u;
    const double analytic = 1.0 - std::exp(-Da);
    const double upwind = 1.0 - std::pow(1.0 + Da / N, -N);
    // The upwind scheme's steady state to the pressure drop and the end of
    // the transient (20 s, 15 to 30 residence times).
    EXPECT_NEAR(X, upwind, 2e-4 * upwind);
    // The analytic plug flow to the scheme's numerical dispersion (O(Da/N)).
    EXPECT_NEAR(X, analytic, 0.01 * analytic);
    // The logged rate is eta k C_A at the outlet cell.
    const auto rate = r.sim.Log().GetReactionRate(0, 0, r.cat, r.rxn).values.back();
    const auto C = LastRow(r.sim, VariableType::MOLAR_CONCENTRATION);
    const auto y = LastRow(r.sim, ComponentVariableType::MOLE_FRACTION, BUT);
    const double k = 0.1 / (r.BulkDensity() * r.o.length);
    EXPECT_NEAR(Outlet(rate), k * Outlet(C) * Outlet(y), 1e-9 * Outlet(rate));
    // The component sources: S_A = -rho_b r, S_B = +rho_b r, and they sum to 0.
    const auto SA = LastRow(r.sim, ComponentVariableType::REACTION_SOURCE, BUT);
    const auto SB = LastRow(r.sim, ComponentVariableType::REACTION_SOURCE, IBUT);
    EXPECT_NEAR(Outlet(SA), -r.BulkDensity() * Outlet(rate), 1e-12);
    EXPECT_DOUBLE_EQ(Outlet(SA), -Outlet(SB));
    EXPECT_EQ(SA.front(), 0.0);  // ghost cell
  }
}

// ---- (c) -----------------------------------------------------------------------

TEST(ReactionBalances, MoleChangingReactionRaisesTheOutletVelocity) {
  Options o;
  o.components = {N2, C2H6, C2H4, H2};
  o.nu = {{C2H6, -1}, {C2H4, 1}, {H2, 1}};
  o.y_feed = {{N2, 0.8}, {C2H6, 0.2}};
  o.k = 5.0 * 0.1 / ((1.0 - o.eb) * o.solid_density * o.length);  // Da = 5
  Reactor r;
  Build(r, o);
  const RunResult run = RunReactor(r);
  ASSERT_TRUE(run.Ok()) << run.message;

  const double yA_in = 0.2;
  const double yA_out = Outlet(LastRow(r.sim, ComponentVariableType::MOLE_FRACTION, C2H6));
  // y_A,out = y_A,in (1 - X) / (1 + y_A,in X)
  const double X = (yA_in - yA_out) / (yA_in * (1.0 + yA_out));
  EXPECT_GT(X, 0.9);
  const auto u = LastRow(r.sim, VariableType::VELOCITY);
  const auto P = LastRow(r.sim, VariableType::PRESSURE);
  // Constant T; the molar flux u P / (R T) grows by 1 + y_A,in X. Corrected
  // for the small pressure drop between the inlet and outlet faces.
  const double P_in = 0.5 * (P[0] + P[1]);
  const double P_out = 0.5 * (P[P.size() - 2] + P.back());
  const double ratio = (u.back() * P_out) / (u.front() * P_in);
  EXPECT_NEAR(ratio, 1.0 + yA_in * X, 2e-4);
  // The uncorrected velocity ratio too, to the pressure drop.
  EXPECT_NEAR(u.back() / u.front(), 1.0 + yA_in * X, 1e-3);
  // Generation over the step: the net moles made equal the moles of C2H4 made.
  const double G_tot = r.sim.Log().ReactionGeneration(0, 0, C2H6) + r.sim.Log().ReactionGeneration(0, 0, C2H4) +
                       r.sim.Log().ReactionGeneration(0, 0, H2);
  EXPECT_GT(G_tot, 0.0);
  EXPECT_NEAR(G_tot, r.sim.Log().ReactionGeneration(0, 0, C2H4), 1e-12 * G_tot);
}

// ---- (d) -----------------------------------------------------------------------

struct AdiabaticResult {
  double dT = 0.0, predicted = 0.0, closure = 0.0;
  double max_solid_minus_gas = 0.0;
};

AdiabaticResult Adiabatic(const std::string& thermal) {
  Options o;
  o.nu = {{BUT, -1}, {IBUT, 1}};
  o.y_feed = {{N2, 0.95}, {BUT, 0.05}};
  o.k = 2.0 * 0.1 / ((1.0 - o.eb) * o.solid_density * o.length);  // Da = 2 at T_in
  o.isothermal = false;
  o.thermal = thermal;
  o.heat_of_reaction = -2e4;
  o.cells = 50;
  o.step_time = 60.0;
  // A light solid, so the thermal front leaves the bed within the step
  // (front speed u (rho cp)_gas / (rho cp)_bed, about 0.02 m/s).
  o.solid_heat_capacity = 5.0;
  Reactor r;
  Build(r, o);
  // Adiabatic column: the constant models' defaults, no bed-to-wall
  // exchange (h_i = 0) and no axial conduction.
  const RunResult run = RunReactor(r);
  EXPECT_TRUE(run.Ok()) << run.message;

  AdiabaticResult out;
  const auto T = LastRow(r.sim, VariableType::TEMPERATURE);
  const double T_in = o.temperature;
  out.dT = Outlet(T) - T_in;
  const double yA_in = 0.05;
  const double X = 1.0 - Outlet(LastRow(r.sim, ComponentVariableType::MOLE_FRACTION, BUT)) / yA_in;
  const double cp = 0.5 * (MolarCp(r.sim, 1) + MolarCp(r.sim, T.size() - 2));
  out.predicted = -o.heat_of_reaction * X * yA_in / cp;

  // Energy closure: the heat released in the bed equals the gas enthalpy
  // gained, F cp (T_out - T_in), F the molar flow (constant: sum(nu) = 0).
  const auto Q = LastRow(r.sim, VariableType::REACTION_HEAT);
  const auto& dV = r.sim.Log().GetDomain().volume;
  double heat = 0.0;
  for (size_t i = 0; i < Q.size(); ++i) heat += Q[i] * dV[i];
  const auto u = LastRow(r.sim, VariableType::VELOCITY);
  const auto P = LastRow(r.sim, VariableType::PRESSURE);
  const double A = r.sim.Log().GetDomain().area[1];
  const double F = u.front() * A * 0.5 * (P[0] + P[1]) / (R * T_in);
  out.closure = heat / (F * cp * out.dT) - 1.0;

  if (thermal == "ltne") {
    const auto Ts = LastRow(r.sim, VariableType::SOLID_TEMPERATURE);
    for (size_t i = 1; i + 1 < T.size(); ++i) out.max_solid_minus_gas = std::max(out.max_solid_minus_gas, Ts[i] - T[i]);
  }
  return out;
}

TEST(ReactionBalances, AdiabaticTemperatureRiseLteAndLtne) {
  const AdiabaticResult lte = Adiabatic("lte");
  EXPECT_GT(lte.dT, 5.0);
  EXPECT_NEAR(lte.dT, lte.predicted, 0.01 * lte.predicted);
  EXPECT_NEAR(lte.closure, 0.0, 0.01);

  const AdiabaticResult ltne = Adiabatic("ltne");
  EXPECT_NEAR(ltne.dT, ltne.predicted, 0.01 * ltne.predicted);
  EXPECT_NEAR(ltne.closure, 0.0, 0.01);
  // The heat is released on the catalyst: the solid runs hotter than the gas.
  EXPECT_GT(ltne.max_solid_minus_gas, 0.01);
}

// ---- (e) -----------------------------------------------------------------------

TEST(ReactionBalances, ElementAndComponentBalancesCloseOverAReactiveStep) {
  Options o;
  o.components = {N2, C2H6, C2H4, H2};
  o.nu = {{C2H6, -1}, {C2H4, 1}, {H2, 1}};
  o.y_feed = {{N2, 0.8}, {C2H6, 0.2}};
  o.k = 1.0 * 0.1 / ((1.0 - o.eb) * o.solid_density * o.length);
  o.cells = 50;
  o.step_time = 4.0;  // a transient step: the inventory changes
  Reactor r;
  Build(r, o);
  ASSERT_TRUE(RunReactor(r).Ok());
  const EquationLog& log = r.sim.Log();
  const auto& dom = log.GetDomain();

  // Gas inventory sum_n et dV C y_i at the first and last logged point.
  auto inventory = [&](const std::string& c, bool end) {
    const auto Ce = log.GetVariable(0, 0, VariableType::MOLAR_CONCENTRATION);
    const auto ye = log.GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION, c);
    const auto& C = end ? Ce.values.back() : Ce.values.front();
    const auto& y = end ? ye.values.back() : ye.values.front();
    double n = 0.0;
    for (size_t i = 1; i + 1 < C.size(); ++i) n += dom.et[i] * dom.volume[i] * C[i] * y[i];
    return n;
  };
  std::map<std::string, double> closure;  // in + G - dN [mol]
  double throughput = 0.0;
  for (const auto& c : o.components) {
    const double in = GetTotalComponentMolesIn(log, 0, 0, BoundaryLocation::BOTTOM, c) +
                      GetTotalComponentMolesIn(log, 0, 0, BoundaryLocation::TOP, c);
    const double G = log.ReactionGeneration(0, 0, c);
    closure[c] = in + G - (inventory(c, true) - inventory(c, false));
    throughput += std::abs(GetTotalComponentMolesIn(log, 0, 0, BoundaryLocation::BOTTOM, c));
  }
  EXPECT_GT(log.ReactionGeneration(0, 0, C2H4), 1e-3 * throughput);
  for (const auto& [c, e] : closure) EXPECT_NEAR(e, 0.0, 2e-3 * throughput) << c;
  // Elements: the reactions make and destroy none (sum_i a_ei G_i = 0).
  for (const auto& element : {"C", "H", "N"}) {
    double g = 0.0, net = 0.0;
    for (const auto& c : o.components) {
      const auto* f = Components::GetElementFormula(c);
      ASSERT_NE(f, nullptr);
      auto it = f->find(element);
      const double a = it == f->end() ? 0.0 : it->second;
      g += a * log.ReactionGeneration(0, 0, c);
      net += a * closure[c];
    }
    EXPECT_NEAR(g, 0.0, 1e-12 * throughput) << element;
    EXPECT_NEAR(net, 0.0, 4e-3 * throughput) << element;
  }
}

// ---- (f) -----------------------------------------------------------------------

// Feed A -> B + C through the bed, then purge it back with the carrier.
void BuildCycle(Reactor& r) {
  Options o;
  o.components = {N2, C2H6, C2H4, H2};
  o.nu = {{C2H6, -1}, {C2H4, 1}, {H2, 1}};
  o.y_feed = {{N2, 0.8}, {C2H6, 0.2}};
  o.k = 1.0 * 0.1 / ((1.0 - o.eb) * o.solid_density * o.length);
  o.cells = 30;
  o.step_time = 4.0;
  Build(r, o);
  UserData& d = r.sim.Data();
  const size_t purge = build::AddStep(d, "Purge");
  build::SetFlowBoundary(d, purge, "top", "inlet");
  build::AddStream(d, purge, "top", "Purge inlet", "velocity", "constant");
  build::SetPressureBoundary(d, purge, "bottom", "outlet", "constant");
  r.sim.RefreshParameters();
  for (const auto& [path, text] : std::vector<std::pair<std::string, std::string>>{
           {"step[Purge].time", "4"},
           {"step[Purge].top[Purge inlet].velocity", "0.1"},
           {"step[Purge].top[Purge inlet].temperature", "400"},
           {"step[Purge].top[Purge inlet].y[Nitrogen]", "1"},
           {"step[Purge].bottom.pressure", "1e5"}})
    r.sim.SetParameterText(path, text);
}

TEST(ReactionBalances, MassBalanceCssMethodsCountGeneration) {
  for (CSSMethod method : {CSSMethod::MASS_BALANCE_TOTAL, CSSMethod::MASS_BALANCE_COMPONENT}) {
    SCOPED_TRACE(static_cast<int>(method));
    Reactor r;
    BuildCycle(r);
    UserData& d = r.sim.Data();
    d.settings.css.method_ = method;
    d.settings.css.mass_balance_tolerance_ = 1e-3;
    const RunResult run = RunReactor(r, 15);
    ASSERT_TRUE(run.Ok()) << run.message;
    EXPECT_TRUE(run.css_reached) << run.message;
    EXPECT_LT(run.cycles_run, 15);

    // The test without generation would read |in - out| / in = the net moles
    // made over the moles fed, far above the tolerance.
    const int c = run.cycles_run - 1;
    double G = 0.0, in = 0.0;
    for (int s = 0; s < 2; ++s) {
      for (const auto& comp : d.fluid.components_) G += r.sim.Log().ReactionGeneration(c, s, comp);
      for (auto loc : {BoundaryLocation::BOTTOM, BoundaryLocation::TOP})
        in += std::max(GetTotalMolesIn(r.sim.Log(), c, s, loc), 0.0);
    }
    EXPECT_GT(G / in, 0.05);
  }
}

// ---- (g) -----------------------------------------------------------------------

double ClosedFormEta(double phi1) { return 3.0 / (phi1 * phi1) * (phi1 / std::tanh(phi1) - 1.0); }

TEST(ReactionBalances, ThieleEffectivenessMatchesTheFirstOrderClosedForm) {
  Reaction rx("A to B", 1);
  rx.nu_ = {{BUT, -1.0}, {IBUT, 1.0}};
  rx.rate_ = CreateReactionRate(ReactionRateType::PowerLaw, rx.nu_);
  rx.basis_ = RateBasis::Concentration;
  rx.effectiveness_.mode_ = EffectivenessMode::Thiele;
  rx.effectiveness_.key_component_ = BUT;
  const double k = 2e-3, rho_p = 1200.0, D = 1e-6;
  dynamic_cast<PowerLawRate&>(*rx.rate_).pre_exp_factor_ = k;
  std::map<std::string, double> y{{N2, 0.8}, {BUT, 0.15}, {IBUT, 0.05}};
  ReactionState s;
  s.T = 500.0;
  s.P = 2e5;
  s.y = &y;
  for (double Rp : {1e-5, 1e-4, 1e-3, 5e-3, 2e-2}) {
    const double phi1 = Rp * std::sqrt(rho_p * k / D);
    EXPECT_NEAR(rx.Effectiveness(s, Rp, rho_p, D), ClosedFormEta(phi1), 1e-10) << phi1;
  }

  // Reversible first order A <-> B, r = k (C_A - C_B / K): linear along the
  // line, so the same closed form with k (1 + 1/K) toward equilibrium.
  rx.SetRateType(ReactionRateType::ReversiblePowerLaw);
  auto& rev = dynamic_cast<ReversiblePowerLawRate&>(*rx.rate_);
  rev.pre_exp_factor_ = k;
  const double K = 0.5;
  rev.equilibrium_pre_exp_factor_ = K;
  double warm = std::numeric_limits<double>::quiet_NaN();
  for (double Rp : {1e-4, 1e-3, 5e-3}) {
    const double phi1 = Rp * std::sqrt(rho_p * k * (1.0 + 1.0 / K) / D);
    EXPECT_NEAR(rx.Effectiveness(s, Rp, rho_p, D, &warm), ClosedFormEta(phi1), 1e-9) << phi1;
    // The equilibrium concentration on the line: C_B = K C_A, C_A + C_B fixed.
    const double C = s.P / (R * s.T);
    EXPECT_NEAR(warm, C * (0.15 + 0.05) / (1.0 + K), 1e-9 * C);
  }
  // Backwards at the surface (C_B > K C_A): the same eta.
  std::map<std::string, double> y_back{{N2, 0.8}, {BUT, 0.02}, {IBUT, 0.18}};
  s.y = &y_back;
  const double Rp = 1e-3;
  EXPECT_NEAR(rx.Effectiveness(s, Rp, rho_p, D), ClosedFormEta(Rp * std::sqrt(rho_p * k * (1.0 + 1.0 / K) / D)), 1e-9);

  // Constant mode returns the user value.
  rx.effectiveness_.mode_ = EffectivenessMode::Constant;
  rx.effectiveness_.effectiveness_ = 0.7;
  EXPECT_EQ(rx.Effectiveness(s, Rp, rho_p, D), 0.7);
}

TEST(ReactionBalances, ThieleEffectivenessInARun) {
  // First order with a constant eta: the plug-flow conversion is
  // 1 - exp(-eta Da) with eta from the cell's D_m.
  Options o;
  o.nu = {{BUT, -1}, {IBUT, 1}};
  o.k = 4.0 * 0.1 / ((1.0 - o.eb) * o.solid_density * o.length);  // Da = 4 at eta = 1
  o.step_time = 10.0;
  Reactor r;
  Build(r, o);
  build::SetEffectiveness(r.sim.Data(), r.cat, r.rxn, "thiele", BUT);
  r.sim.RefreshParameters();
  r.sim.SetParameterText("adsorbent[Catalyst].reaction[R1].tortuosity", "4");
  ASSERT_TRUE(RunReactor(r).Ok());
  const auto rate = r.sim.Log().GetReactionRate(0, 0, r.cat, r.rxn).values.back();
  const auto C = LastRow(r.sim, VariableType::MOLAR_CONCENTRATION);
  const auto y = LastRow(r.sim, ComponentVariableType::MOLE_FRACTION, BUT);
  const auto Dm = LastRow(r.sim, ComponentVariableType::MOLECULAR_DIFFUSION, BUT);
  const size_t n = 10;
  const double D_eff = 0.4 / 4.0 * Dm.at(n);
  const double phi1 = 0.5 * o.particle_diameter * std::sqrt(o.solid_density * o.k / D_eff);
  const double eta = ClosedFormEta(phi1);
  EXPECT_LT(eta, 0.9);
  EXPECT_NEAR(rate.at(n), eta * o.k * C.at(n) * y.at(n), 1e-8 * rate.at(n));
  // The upwind steady state with each cell's own eta (D_m follows the
  // composition): y_n = y_{n-1} / (1 + 4 eta_n / N).
  const int N = o.cells;
  double ratio = 1.0;
  for (int i = 1; i <= N; ++i) {
    const double phi = 0.5 * o.particle_diameter * std::sqrt(o.solid_density * o.k / (0.4 / 4.0 * Dm.at(i)));
    ratio /= 1.0 + 4.0 * ClosedFormEta(phi) / N;
  }
  const double X = 1.0 - Outlet(y) / 0.1;
  EXPECT_NEAR(X, 1.0 - ratio, 2e-4);
}

// ---- (h) -----------------------------------------------------------------------

TEST(ReactionBalances, InertLoadingsAreNotUnknowns) {
  auto run = [](const std::string& isotherm, int& size) {
    Options o;
    o.nu = {{BUT, -1}, {IBUT, 1}};
    o.k = 0.1 / ((1.0 - o.eb) * o.solid_density * o.length);
    o.cells = 40;
    o.step_time = 5.0;
    o.isotherm = isotherm;
    Reactor r;
    Build(r, o);
    if (isotherm == "henry")
      for (const auto& c : o.components) r.sim.SetParameterText("adsorbent[Catalyst].isotherm[" + c + "].k", "0");
    const RunResult run = RunReactor(r);
    EXPECT_TRUE(run.Ok()) << run.message;
    size = r.sim.Data().state.GetSystemSize();
    // The inert loadings are logged, and exactly zero (a Henry k = 0 loading
    // carries round-off).
    for (const auto& c : o.components)
      if (isotherm == "inert")
        for (double q : LastRow(r.sim, ComponentVariableType::ADSORBED_CONCENTRATION, c)) EXPECT_EQ(q, 0.0);
    return Outlet(LastRow(r.sim, ComponentVariableType::MOLE_FRACTION, BUT));
  };
  int n_inert = 0, n_henry = 0;
  const double y_inert = run("inert", n_inert);
  const double y_henry = run("henry", n_henry);
  // 3 components x (40 bed cells + 2 ghosts) fewer unknowns.
  EXPECT_EQ(n_henry - n_inert, 3 * 42);
  EXPECT_NEAR(y_inert, y_henry, 1e-6 * y_henry);
}

TEST(ReactionBalances, SeriesManifestListsTheReactionFields) {
  const std::string m = SeriesManifest();
  EXPECT_NE(m.find("\"reaction_heat\""), std::string::npos);
  EXPECT_NE(m.find("\"reaction_source\""), std::string::npos);
  EXPECT_NE(m.find("{reaction} Rate ({adsorbent})"), std::string::npos);
}

}  // namespace
}  // namespace adsim

namespace adsim {
namespace {

// Jacobian cost on a reactive deck (Coloured mode, solver/jacobian.h): the
// rate couples the unknowns of a cell (y_j, P, T), all inside the existing
// band, so the learned pattern only fills in the cell's block. Reports the
// colours and the residual calls per Jacobian with and without the reaction
// (a zero rate learns the pattern of the deck without it).
TEST(ReactionBalances, JacobianColoursOnAReactiveDeck) {
  struct Cost { int colours = 0; long evaluations = 0; long residual_calls = 0; };
  auto measure = [](double da) {
    Options o;
    o.components = {N2, C2H6, C2H4, H2};
    o.nu = {{C2H6, -1}, {C2H4, 1}, {H2, 1}};
    o.y_feed = {{N2, 0.8}, {C2H6, 0.2}};
    o.k = da * 0.1 / ((1.0 - o.eb) * o.solid_density * o.length);
    o.isothermal = false;
    o.heat_of_reaction = -2e4;
    o.solid_heat_capacity = 50.0;
    o.cells = 40;
    o.step_time = 10.0;
    Reactor r;
    Build(r, o);
    Cost c;
    BandedJacobian::SetEvaluationObserverForTesting([&](const BandedJacobian::EvaluationView& v) {
      ++c.evaluations;
      if (v.udata && v.udata->jacobian_) {
        c.colours = v.udata->jacobian_->NumColours();
        c.residual_calls = v.udata->jacobian_->NumResidualEvals();
      }
    });
    const RunResult run = RunReactor(r);
    BandedJacobian::SetEvaluationObserverForTesting({});
    EXPECT_TRUE(run.Ok()) << run.message;
    return c;
  };
  const Cost off = measure(0.0);
  const Cost on = measure(2.0);
  std::cout << "[reaction jacobian] colours " << off.colours << " -> " << on.colours << "; residual calls per Jacobian "
            << static_cast<double>(off.residual_calls) / std::max(1L, off.evaluations) << " -> "
            << static_cast<double>(on.residual_calls) / std::max(1L, on.evaluations) << " (" << off.evaluations
            << " / " << on.evaluations << " Jacobians)\n";
  ASSERT_GT(off.colours, 0);
  EXPECT_GE(on.colours, off.colours);
}

}  // namespace
}  // namespace adsim
