/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// Analytic acceptance tests of reactions in the bed (objects/reaction.h,
// UpdateReactions in equations/adsorbent_equations.h), complementing
// test_reaction_balances.cpp (plug-flow first order, velocity rise, adiabatic
// rise, single-step closure, CSS methods, Thiele, zero rate):
//
//  (a) first order A -> B with axial dispersion between Danckwerts closed-closed
//      boundaries, against the Wehner-Wilhelm closed form at Pe = 5 and 50,
//      with a grid-refinement check of the scheme's error;
//  (b) reversible A <-> B: equilibrium conversion y_B / y_A = K(T) at a long
//      residence time, at two temperatures (van't Hoff), and below it at a
//      short one, against the scheme's own steady state;
//  (c) Kirchhoff heat of reaction: -Q / (rho_b eta r) is dH_ref +
//      sum_i nu_i integral Cp_i dT of the engine's own Cp(T) in every cell, the
//      constant mode is dH_ref, and the adiabatic Kirchhoff run conserves the
//      enthalpy of the stream;
//  (d) A -> 2B at a constant outlet pressure: the stream logs' molar outflow /
//      inflow is 1 + y_A,in X, and X follows the variable-density plug flow;
//  (e) element, component and energy closure on a cyclic reactive deck over
//      several cycles (per step and per cycle);
//  (f) a product starting at y = 0 stays above -atol and leaves at the
//      conversion-consistent fraction.
//
// Every tolerance is set from a measured discretisation or integration error
// (quoted next to it), with a margin of 2-5x; none is fitted to the run. The
// whole file runs in about 2 s (20-160 cells, steps of 5-60 s), so the grid
// refinement of (a) stays in the fast tier.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/objects/deck_builder.h"
#include "adsim/objects/reaction.h"
#include "adsim/simulation.h"
#include "adsim/solver/kpis.h"

namespace adsim {
namespace {

const std::string N2 = "Nitrogen", BUT = "Butane", IBUT = "Isobutane", C2H6 = "Ethane", C2H4 = "Ethylene",
                  H2 = "Hydrogen", C6H12 = "Cyclohexane", C3H6 = "Propylene";

std::string Num(double v) {
  char buf[40];
  std::snprintf(buf, sizeof(buf), "%.17g", v);
  return buf;
}

// A packed catalyst bed (inert isotherms) with one feed step; the minimal
// setup of test_reaction_balances.cpp plus the switches these cases need.
struct Options {
  std::vector<std::string> components{N2, BUT, IBUT};  // the carrier first
  std::map<std::string, double> nu{{BUT, -1}, {IBUT, 1}};
  std::map<std::string, double> y_feed{{N2, 0.9}, {BUT, 0.1}};
  std::string kind = "power_law";
  double k = 0.0;                      // pre-exponential factor (first order: m3/(kg s))
  double K_pre = 1.0, K_energy = 0.0;  // reversible: K = K_pre exp(-K_energy / (R T))
  double velocity = 0.1;               // [m/s]
  double length = 0.2;                 // [m]
  int cells = 30;
  double step_time = 20.0;  // [s]
  bool isothermal = true;
  double heat_of_reaction = 0.0;  // [J/mol]
  std::string heat_mode = "constant";
  double temperature = 400.0;          // [K]
  double solid_density = 1000.0;       // [kg/m3]
  double solid_heat_capacity = 900.0;  // [J/(kg K)]
  double eb = 0.4;
  double particle_diameter = 5e-3;  // [m]
  double dispersion = 0.0;          // constant D_ax [m2/s]
  AdvectionScheme scheme = AdvectionScheme::Upwind;
};

struct Reactor {
  Simulation sim;
  AdsorbentID cat = 0;
  size_t layer = 0, feed = 0, feed_in = 0, outlet = 0;
  ReactionID rxn = 0;
  Options o;
  double BulkDensity() const { return (1.0 - o.eb) * o.solid_density; }
};

void Build(Reactor& r, const Options& o) {
  r.o = o;
  UserData& d = r.sim.Data();
  for (const auto& c : o.components) build::AddComponent(d, c);
  r.cat = build::AddAdsorbent(d, "Catalyst");
  build::SetBedShape(d, "cylinder");
  build::SetWallMaterial(d, "ss316l");
  build::SetWallBoundary(d, "zero_flux");
  build::SetThermalEquilibrium(d, "lte");
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
  r.rxn = build::AddReaction(d, r.cat, "R1", o.nu, o.kind);
  build::SetReactionBasis(d, r.cat, r.rxn, "concentration");
  r.sim.RefreshParameters();

  const std::string rx = "adsorbent[Catalyst].reaction[R1].";
  std::vector<std::pair<std::string, std::string>> values = {
      {"bed.diameter", "0.05"},
      {"adsorbent[Catalyst].solid_density", Num(o.solid_density)},
      {"adsorbent[Catalyst].heat_capacity", Num(o.solid_heat_capacity)},
      {"adsorbent[Catalyst].particle_void_fraction", "0.4"},
      {"layer[Bed].length", Num(o.length)},
      {"layer[Bed].cells", Num(o.cells)},
      {"layer[Bed].bed_void_fraction", Num(o.eb)},
      {"layer[Bed].particle_diameter", Num(o.particle_diameter)},
      {"layer[Bed].dispersion_coefficient", Num(o.dispersion)},
      {"initial.pressure", "1e5"},
      {"initial.temperature", Num(o.temperature)},
      {"initial.y[" + o.components.front() + "]", "1"},
      {"step[Feed].time", Num(o.step_time)},
      {"step[Feed].bottom[Feed inlet].velocity", Num(o.velocity)},
      {"step[Feed].bottom[Feed inlet].temperature", Num(o.temperature)},
      {"step[Feed].top.pressure", "1e5"},
      {rx + "rate.pre_exp_factor", Num(o.k)},
      {rx + "heat_of_reaction", Num(o.heat_of_reaction)},
      {rx + "heat_mode", o.heat_mode},
  };
  if (o.kind == "reversible_power_law") {
    values.push_back({rx + "rate.equilibrium_pre_exp_factor", Num(o.K_pre)});
    values.push_back({rx + "rate.equilibrium_energy", Num(o.K_energy)});
  }
  for (const auto& [c, y] : o.y_feed) values.push_back({"step[Feed].bottom[Feed inlet].y[" + c + "]", Num(y)});
  for (const auto& [path, text] : values) r.sim.SetParameterText(path, text);

  d.settings.advection_scheme_ = o.scheme;
  d.settings.css.method_ = CSSMethod::NONE;
  d.cycle.number_of_cycles_ = 1;
  d.bed.equations_.energy_balance_mode_ =
      o.isothermal ? EnergyBalanceMode::ISOTHERMAL : EnergyBalanceMode::NON_ISOTHERMAL;
}

RunResult RunReactor(Reactor& r, int cycles = 1) {
  RunOptions options;
  options.cycles = cycles;
  options.on_message = [](LogLevel, const std::string&) {};
  return r.sim.Run(options);
}

std::vector<double> LastRow(const Simulation& sim, VariableType t, int cycle = 0, int step = 0) {
  const auto e = sim.Log().GetVariable(cycle, step, t);
  return e.values.empty() ? std::vector<double>{} : e.values.back();
}
std::vector<double> LastRow(const Simulation& sim, ComponentVariableType t, const std::string& c, int cycle = 0,
                            int step = 0) {
  const auto e = sim.Log().GetComponentVariable(cycle, step, t, c);
  return e.values.empty() ? std::vector<double>{} : e.values.back();
}
// The last bed cell (an outlet ghost mirrors it: zero gradient).
double Outlet(const std::vector<double>& row) { return row.at(row.size() - 2); }

// ---- (a) Danckwerts closed-closed reactor ----------------------------------
//
// Steady first order A -> B with constant u and C (sum(nu) = 0), in the
// engine's variables (docs/mass-balance.md: the dispersive flux is
// eb D C dy/dz):
//   u dy/dz = eb D d2y/dz2 - rho_b k y,   Pe = u L / (eb D),  Da = rho_b k L / u.
// How the engine's ends map to the closed-closed assumptions (docs/
// numerics.md, "Boundary faces"):
//   * inlet: no dispersive flux crosses the boundary face, so the feed's
//     u C y_in is the whole flux into the first cell, which balances the
//     convective-plus-dispersive flux inside the bed: the finite-volume form of
//     u y_in = u y(0+) - eb D dy/dz(0+) (Danckwerts; y(0+) differs from y_in);
//   * outlet: the ghost cell mirrors the last cell, dy/dz(L) = 0, and no
//     dispersive flux leaves through the face (no dispersion downstream).
// The closed form for these conditions (Wehner-Wilhelm):
//   1 - X = 4 a exp(Pe/2) / [(1+a)^2 exp(a Pe/2) - (1-a)^2 exp(-a Pe/2)],
//   a = sqrt(1 + 4 Da / Pe).
double DanckwertsConversion(double Pe, double Da) {
  const double a = std::sqrt(1.0 + 4.0 * Da / Pe);
  // Divided through by exp(a Pe / 2) to keep the exponents bounded.
  const double num = 4.0 * a * std::exp(0.5 * Pe * (1.0 - a));
  const double den = (1.0 + a) * (1.0 + a) - (1.0 - a) * (1.0 - a) * std::exp(-a * Pe);
  return 1.0 - num / den;
}

double DanckwertsRun(double Pe, double Da, int cells, AdvectionScheme scheme) {
  Options o;
  o.cells = cells;
  o.scheme = scheme;
  // 2 cm particles: a Kozeny-Carman drop of ~0.2 Pa, so C and u are uniform to
  // ~2e-6. (With 5 mm particles the ~3 Pa drop leaves a 2e-5 floor on X that
  // no refinement removes; measured.)
  o.particle_diameter = 2e-2;
  o.k = Da * o.velocity / ((1.0 - o.eb) * o.solid_density * o.length);
  o.dispersion = o.velocity * o.length / (o.eb * Pe);
  Reactor r;
  Build(r, o);
  const RunResult run = RunReactor(r);
  EXPECT_TRUE(run.Ok()) << run.message;
  return 1.0 - Outlet(LastRow(r.sim, ComponentVariableType::MOLE_FRACTION, BUT)) / o.y_feed.at(BUT);
}

TEST(ReactionAnalytic, DanckwertsClosedClosedFirstOrder) {
  const double Da = 2.0;
  for (double Pe : {5.0, 50.0}) {
    SCOPED_TRACE(Pe);
    const double X = DanckwertsConversion(Pe, Da);
    // Dispersion matters here: plug flow would give 1 - exp(-Da) = 0.865
    // (Pe = 5: X = 0.796; Pe = 50: X = 0.854).
    EXPECT_GT(1.0 - std::exp(-Da) - X, 0.01);

    // Upwind convection adds a numerical dispersion u dz / 2, so its error is
    // first order: measured 5.76e-3, 2.84e-3, 1.41e-3, 7.03e-4 (Pe = 5) and
    // 1.20e-2 ... 1.51e-3 (Pe = 50) on 20-160 cells, halving per refinement to
    // within 1.5%.
    const double X40 = DanckwertsRun(Pe, Da, 40, AdvectionScheme::Upwind);
    const double X80 = DanckwertsRun(Pe, Da, 80, AdvectionScheme::Upwind);
    const double X160 = DanckwertsRun(Pe, Da, 160, AdvectionScheme::Upwind);
    EXPECT_LT(X40, X80);  // numerical dispersion lowers X; refinement raises it
    EXPECT_LT(X80, X160);
    EXPECT_LT(X160, X);
    const double ratio1 = (X - X40) / (X - X80), ratio2 = (X - X80) / (X - X160);
    EXPECT_NEAR(ratio1, 2.0, 0.1);
    EXPECT_NEAR(ratio2, 2.0, 0.1);
    // Richardson extrapolation removes the first-order term; what is left is
    // O(dz^2), estimated by the change of the extrapolant between the two
    // pairs (measured 1.6e-5 and 1.1e-5 for an actual error of 4.4e-6 and
    // 4.9e-6), plus the 2e-6 pressure-drop floor.
    const double R1 = 2.0 * X80 - X40, R2 = 2.0 * X160 - X80;
    EXPECT_NEAR(R2, X, std::abs(R2 - R1) + 2e-6);

    // The high-resolution (van Leer) scheme on 40 cells. Its error over 20-160
    // cells was measured at most 3.3e-5 (the limiter makes the convergence
    // non-monotone: -1.0e-5, 3.3e-5, 1.2e-5, 1.9e-6 at Pe = 50); bound 1e-4.
    EXPECT_NEAR(DanckwertsRun(Pe, Da, 40, AdvectionScheme::VanLeer), X, 1e-4);
  }
}

// ---- (b) reversible A <-> B --------------------------------------------------

// r = k (C_A - C_B / K), K(T) = K_pre exp(-E_K / (R T)) (van't Hoff). With
// sum(nu) = 0 and a uniform C, each upwind cell relaxes y_A - y_A,eq by
// 1 / (1 + Da' / N), Da' = Da (1 + 1 / K), toward y_A,eq = y_A,in / (1 + K).
struct ReversibleResult {
  double K = 0, ratio = 0, yA = 0, yA_scheme = 0;
};

ReversibleResult ReversibleRun(double T, double Da, double K_pre, double E_K) {
  Options o;
  o.kind = "reversible_power_law";
  o.temperature = T;
  o.K_pre = K_pre;
  o.K_energy = E_K;
  o.particle_diameter = 2e-2;  // pressure drop ~0.2 Pa (see (a))
  o.k = Da * o.velocity / ((1.0 - o.eb) * o.solid_density * o.length);
  Reactor r;
  Build(r, o);
  const RunResult run = RunReactor(r);
  EXPECT_TRUE(run.Ok()) << run.message;
  ReversibleResult out;
  out.K = K_pre * std::exp(-E_K / (R * T));
  out.yA = Outlet(LastRow(r.sim, ComponentVariableType::MOLE_FRACTION, BUT));
  out.ratio = Outlet(LastRow(r.sim, ComponentVariableType::MOLE_FRACTION, IBUT)) / out.yA;
  const double yA_in = o.y_feed.at(BUT), y_eq = yA_in / (1.0 + out.K);
  out.yA_scheme = y_eq + (yA_in - y_eq) * std::pow(1.0 + Da * (1.0 + 1.0 / out.K) / o.cells, -o.cells);
  return out;
}

TEST(ReactionAnalytic, ReversibleReactionReachesEquilibriumConversion) {
  // Exothermic: K(400 K) = 3, K(500 K) = 1.5, E_K = -R ln 2 / (1/400 - 1/500).
  const double E_K = -R * std::log(2.0) / (1.0 / 400.0 - 1.0 / 500.0);
  const double K_pre = 3.0 * std::exp(E_K / (R * 400.0));
  for (double T : {400.0, 500.0}) {
    SCOPED_TRACE(T);
    // Long residence time, Da = 20 (Da' = 27-33): the scheme leaves
    // y_A - y_A,eq below 1e-8 of y_A,in at the outlet, so the outlet ratio is K
    // to the integrator's accuracy. Measured |ratio / K - 1| = 2.1e-8 (400 K)
    // and 4.6e-10 (500 K); bound 1e-6, the integrator's reltol.
    const ReversibleResult eq = ReversibleRun(T, 20.0, K_pre, E_K);
    EXPECT_NEAR(eq.ratio / eq.K, 1.0, 1e-6) << eq.ratio << " vs K " << eq.K;
    // Short residence time, Da = 0.3: well short of equilibrium, and on the
    // scheme's own steady state. Measured |y_A / y_A,scheme - 1| = 1.4e-6 and
    // 1.5e-6 (the remaining pressure drop); bound 1e-5.
    const ReversibleResult kin = ReversibleRun(T, 0.3, K_pre, E_K);
    EXPECT_LT(kin.ratio, 0.25 * kin.K);
    EXPECT_NEAR(kin.yA / kin.yA_scheme, 1.0, 1e-5);
  }
  // The van't Hoff dependence moved the equilibrium: K halves from 400 to 500 K.
  EXPECT_NEAR(K_pre * std::exp(-E_K / (R * 500.0)), 1.5, 1e-12);
}

// ---- (c) Kirchhoff heat of reaction --------------------------------------------

// integral_{T1}^{T2} Cp_i dT of the engine's own Cp_i(T) (Fluid::
// MolarHeatCapacity) by composite Simpson, independent of the closed-form
// enthalpy the engine integrates (Fluid::MolarEnthalpyChange). The NASA
// polynomials are quartic in T: 400 intervals are exact to ~1e-14.
double SimpsonEnthalpy(const Fluid& f, const std::string& c, double T1, double T2, int n = 400) {
  const double h = (T2 - T1) / n;
  double s = f.MolarHeatCapacity(c, T1) + f.MolarHeatCapacity(c, T2);
  for (int i = 1; i < n; ++i) s += (i % 2 ? 4.0 : 2.0) * f.MolarHeatCapacity(c, T1 + i * h);
  return s * h / 3.0;
}

struct KirchhoffResult {
  double max_cell_error = 0;  // max over cells |dH_engine(T) - dH_expected(T)| / |dH_ref|
  double term_in = 0, term_out = 0;  // sum_i nu_i int_{T_ref}^{T} Cp_i dT at T_in and T_out [J/mol]
  double enthalpy_closure = 0;  // (H_in - H_out - xi dH_ref) / (xi |dH_ref|)
  double dH_ref = 0;
};

// Ethane -> ethylene + H2 (dCp ~ +20 J/(mol K)) in an adiabatic bed fed at
// 700 K, 400 K above T_ref = 298.15 K: the Kirchhoff term is ~6 kJ/mol, 20%
// of dH_ref.
KirchhoffResult KirchhoffRun(const std::string& mode, int cells) {
  Options o;
  o.components = {N2, C2H6, C2H4, H2};
  o.nu = {{C2H6, -1}, {C2H4, 1}, {H2, 1}};
  o.y_feed = {{N2, 0.9}, {C2H6, 0.1}};
  o.isothermal = false;
  o.heat_mode = mode;
  o.heat_of_reaction = -3e4;
  o.temperature = 700.0;
  o.k = 2.0 * 0.1 / ((1.0 - o.eb) * o.solid_density * o.length);
  o.cells = cells;
  o.step_time = 60.0;
  o.solid_heat_capacity = 5.0;  // a light solid: the thermal front leaves within the step
  Reactor r;
  Build(r, o);
  const RunResult run = RunReactor(r);
  EXPECT_TRUE(run.Ok()) << run.message;
  KirchhoffResult out;
  out.dH_ref = o.heat_of_reaction;
  const Fluid& fluid = r.sim.Data().fluid;
  const double T_ref = 298.15;  // the default reference temperature
  auto term = [&](double T) {
    double sum = 0.0;
    for (const auto& [c, v] : o.nu) sum += v * SimpsonEnthalpy(fluid, c, T_ref, T);
    return sum;
  };
  // Cell by cell: Q = rho_b (-dH_r(T)) eta r (LTE: at the gas temperature).
  const auto T = LastRow(r.sim, VariableType::TEMPERATURE);
  const auto Q = LastRow(r.sim, VariableType::REACTION_HEAT);
  const auto rate = r.sim.Log().GetReactionRate(0, 0, r.cat, r.rxn).values.back();
  int cells_checked = 0;
  for (size_t n = 1; n + 1 < T.size(); ++n) {
    if (rate[n] <= 0.0) continue;
    const double dH_engine = -Q[n] / (r.BulkDensity() * rate[n]);
    const double expected = mode == "kirchhoff" ? out.dH_ref + term(T[n]) : out.dH_ref;
    out.max_cell_error = std::max(out.max_cell_error, std::abs(dH_engine - expected) / std::abs(out.dH_ref));
    ++cells_checked;
  }
  EXPECT_EQ(cells_checked, cells);
  // Steady stream enthalpy with h_i(T) = int_{T_ref}^{T} Cp_i dT and dH_ref
  // per mole of extent: sum F_i,in h_i(T_in) - xi dH_ref = sum F_i,out h_i(T_out),
  // xi = F_C2H4,out (none fed). Exact for the Kirchhoff heat; the constant
  // heat misses xi sum nu_i int Cp_i dT.
  const StreamData* in = r.sim.Boundaries().GetStreamData(0, 0, "Feed inlet");
  const StreamData* outlet = r.sim.Boundaries().GetStreamData(0, 0, "Feed Pressure Boundary Top");
  EXPECT_TRUE(in && outlet);
  if (!in || !outlet) return out;
  const double T_in = in->temperature_boundary_.back(), T_out = outlet->temperature_boundary_.back();
  double H_in = 0.0, H_out = 0.0;
  for (const auto& c : o.components) {
    H_in += std::abs(in->component_molar_flow_.at(c).back()) * SimpsonEnthalpy(fluid, c, T_ref, T_in);
    H_out += std::abs(outlet->component_molar_flow_.at(c).back()) * SimpsonEnthalpy(fluid, c, T_ref, T_out);
  }
  const double xi = std::abs(outlet->component_molar_flow_.at(C2H4).back());
  EXPECT_GT(T_out - T_in, 30.0);
  out.enthalpy_closure = (H_in - H_out - xi * out.dH_ref) / (xi * std::abs(out.dH_ref));
  out.term_in = term(T_in);
  out.term_out = term(T_out);
  return out;
}

TEST(ReactionAnalytic, KirchhoffHeatOfReactionFollowsTheEnginesHeatCapacity) {
  const KirchhoffResult k30 = KirchhoffRun("kirchhoff", 30);
  const KirchhoffResult c30 = KirchhoffRun("constant", 30);
  // dCp != 0 matters here: the Kirchhoff term is ~20% of dH_ref.
  EXPECT_GT(k30.term_in, 0.15 * std::abs(k30.dH_ref));
  // In every cell -Q / (rho_b eta r) is dH_ref + sum nu_i int Cp_i dT of the
  // engine's own Cp(T) (kirchhoff) or dH_ref (constant), to round-off:
  // measured 2.0e-13 and 1.2e-16 of dH_ref.
  EXPECT_LT(k30.max_cell_error, 1e-10);
  EXPECT_LT(c30.max_cell_error, 1e-12);

  // The adiabatic Kirchhoff run conserves the stream enthalpy to the energy
  // balance's discretisation error (convection of T in non-conservative form,
  // first order): measured -6.1e-4 on 30 cells, -2.8e-4 on 60.
  const KirchhoffResult k60 = KirchhoffRun("kirchhoff", 60);
  EXPECT_LT(std::abs(k30.enthalpy_closure), 1.5e-3);
  EXPECT_LT(std::abs(k60.enthalpy_closure), 0.6 * std::abs(k30.enthalpy_closure));
  // The constant heat misses the Kirchhoff term, between its values at the
  // inlet and outlet temperatures (measured -0.194 in [-0.200, -0.185]), up
  // to the same discretisation error.
  const double lo = -c30.term_out / std::abs(c30.dH_ref) - 1.5e-3;
  const double hi = -c30.term_in / std::abs(c30.dH_ref) + 1.5e-3;
  EXPECT_GT(c30.enthalpy_closure, lo);
  EXPECT_LT(c30.enthalpy_closure, hi);
}

// ---- (d) A -> 2B at a constant outlet pressure ---------------------------------

// Isothermal, isobaric plug flow, first order in C_A, eps = y_A,in
// (Levenspiel): (1 + eps) ln(1 / (1 - X)) - eps X = Da = rho_b k L / u_in.
double VariableDensityPlugFlow(double Da, double eps) {
  double X = 0.5;
  for (int i = 0; i < 50; ++i) {
    const double f = -(1.0 + eps) * std::log(1.0 - X) - eps * X - Da;
    X -= f / ((1.0 + eps) / (1.0 - X) - eps);
  }
  return X;
}
// The upwind scheme's steady state: per cell, with f = F_A / F_0 and the total
// flow F / F_0 = 1 + y0 - f, f_n - f_{n-1} = -(Da / N) f_n / (1 + y0 - f_n).
double VariableDensityUpwind(double Da, double y0, int N) {
  double f = y0;
  const double a = Da / N;
  for (int n = 0; n < N; ++n) {
    // f^2 - f (1 + y0 + a + f_prev) + f_prev (1 + y0) = 0, the smaller root.
    const double b = 1.0 + y0 + a + f, c = f * (1.0 + y0);
    f = 0.5 * (b - std::sqrt(b * b - 4.0 * c));
  }
  return 1.0 - f / y0;
}

struct MoleChangingResult {
  double X = 0, flow_ratio = 0;
};

// Cyclohexane -> 2 propylene (element balanced), 30% in nitrogen, at 500 K.
MoleChangingResult MoleChangingRun(double Da, int cells) {
  Options o;
  o.components = {N2, C6H12, C3H6};
  o.nu = {{C6H12, -1}, {C3H6, 2}};
  o.y_feed = {{N2, 0.7}, {C6H12, 0.3}};
  o.temperature = 500.0;
  o.particle_diameter = 2e-2;  // pressure drop ~0.2 Pa (see (a))
  o.cells = cells;
  o.k = Da * o.velocity / ((1.0 - o.eb) * o.solid_density * o.length);
  Reactor r;
  Build(r, o);
  const RunResult run = RunReactor(r);
  EXPECT_TRUE(run.Ok()) << run.message;
  MoleChangingResult out;
  const double y0 = o.y_feed.at(C6H12);
  const double yA = Outlet(LastRow(r.sim, ComponentVariableType::MOLE_FRACTION, C6H12));
  out.X = (y0 - yA) / (y0 * (1.0 + yA));  // y_A,out = y0 (1 - X) / (1 + y0 X)
  const StreamData* in = r.sim.Boundaries().GetStreamData(0, 0, "Feed inlet");
  const StreamData* outlet = r.sim.Boundaries().GetStreamData(0, 0, "Feed Pressure Boundary Top");
  EXPECT_TRUE(in && outlet);
  if (in && outlet) out.flow_ratio = std::abs(outlet->molar_flow_.back()) / std::abs(in->molar_flow_.back());
  return out;
}

TEST(ReactionAnalytic, MoleChangingReactionMolarFlowRatio) {
  const double Da = 1.5, y0 = 0.3;
  const MoleChangingResult r30 = MoleChangingRun(Da, 30);
  const MoleChangingResult r60 = MoleChangingRun(Da, 60);
  EXPECT_GT(r30.X, 0.5);
  // The stream logs: molar outflow / inflow = 1 + y_A,in X (one extra mole per
  // mole of A converted). Measured to 7e-14; bound 1e-9.
  for (const auto& r : {r30, r60}) EXPECT_NEAR(r.flow_ratio, 1.0 + y0 * r.X, 1e-9);
  // X on the scheme's own steady state, to the remaining pressure drop:
  // measured 2.6e-6 relative; bound 2e-5.
  EXPECT_NEAR(r30.X, VariableDensityUpwind(Da, y0, 30), 2e-5 * r30.X);
  EXPECT_NEAR(r60.X, VariableDensityUpwind(Da, y0, 60), 2e-5 * r60.X);
  // ... and the variable-density plug flow, first order in 1/N: measured
  // errors 8.2e-3 and 4.1e-3 (ratio 1.98); the Richardson extrapolant leaves
  // 6.7e-5 (its O(dz^2) remainder); bound 3e-4.
  const double X = VariableDensityPlugFlow(Da, y0);
  EXPECT_NEAR((X - r30.X) / (X - r60.X), 2.0, 0.1);
  EXPECT_NEAR(2.0 * r60.X - r30.X, X, 3e-4);
}

// ---- (e) cyclic closure ---------------------------------------------------------

double Trapezoid(const std::vector<double>& t, const std::vector<double>& f) {
  double s = 0.0;
  for (size_t i = 1; i < t.size(); ++i) s += 0.5 * (f[i] + f[i - 1]) * (t[i] - t[i - 1]);
  return s;
}

// Butane -> isobutane, exothermic, adiabatic (no bed-to-wall exchange): a feed
// step, then a hold with both ends closed in which the reaction goes on and
// heats the bed, repeated. Equal constant gas heat capacities, so dCp = 0 and
// the constant heat of reaction is exact thermodynamics.
TEST(ReactionAnalytic, ElementAndEnergyClosureOverSeveralReactiveCycles) {
  Options o;
  o.isothermal = false;
  o.heat_of_reaction = -4e4;
  o.k = 1.0 * 0.1 / ((1.0 - o.eb) * o.solid_density * o.length);
  o.step_time = 5.0;
  o.solid_heat_capacity = 200.0;
  Reactor r;
  Build(r, o);
  UserData& d = r.sim.Data();
  const size_t hold = build::AddStep(d, "Hold");
  build::SetClosed(d, hold, "bottom");
  build::SetClosed(d, hold, "top");
  r.sim.RefreshParameters();
  r.sim.SetParameterText("step[Hold].time", "5");
  r.sim.SetParameterText("fluid.gas_heat_capacity_mode", "constant");
  r.sim.RefreshParameters();
  const double Cp = 40.0, T_ref = 298.15;  // [J/(mol K)]
  for (const auto& c : o.components) r.sim.SetParameterText("fluid.gas_heat_capacity[" + c + "]", Num(Cp));
  const int cycles = 4;
  const RunResult run = RunReactor(r, cycles);
  ASSERT_TRUE(run.Ok()) << run.message;
  const EquationLog& log = r.sim.Log();
  const auto& dom = log.GetDomain();

  for (int c = 0; c < cycles; ++c) {
    SCOPED_TRACE("cycle " + std::to_string(c));
    std::map<std::string, double> in_cycle, gen_cycle;
    double feed_throughput = 0.0;
    std::map<std::string, double> inv0, inv1;
    for (int s = 0; s < 2; ++s) {
      SCOPED_TRACE("step " + std::to_string(s));
      const auto Ce = log.GetVariable(c, s, VariableType::MOLAR_CONCENTRATION);
      const auto Te = log.GetVariable(c, s, VariableType::TEMPERATURE);
      const auto Pe = log.GetVariable(c, s, VariableType::PRESSURE);
      const auto Qe = log.GetVariable(c, s, VariableType::REACTION_HEAT);
      ASSERT_FALSE(Ce.values.empty());
      const size_t last = Ce.values.size() - 1;

      // Component balance of the step: moles in through the faces +
      // generation - change of the gas inventory (inert solid: no loadings).
      std::map<std::string, double> closure;
      for (const auto& comp : o.components) {
        const auto ye = log.GetComponentVariable(c, s, ComponentVariableType::MOLE_FRACTION, comp);
        auto inventory = [&](size_t k) {
          double n = 0.0;
          for (size_t i = 1; i + 1 < Ce.values[k].size(); ++i)
            n += dom.et[i] * dom.volume[i] * Ce.values[k][i] * ye.values[k][i];
          return n;
        };
        const double in = GetTotalComponentMolesIn(log, c, s, BoundaryLocation::BOTTOM, comp) +
                          GetTotalComponentMolesIn(log, c, s, BoundaryLocation::TOP, comp);
        const double G = log.ReactionGeneration(c, s, comp);
        closure[comp] = in + G - (inventory(last) - inventory(0));
        in_cycle[comp] += in;
        gen_cycle[comp] += G;
        if (s == 0) inv0[comp] = inventory(0);
        if (s == 1) inv1[comp] = inventory(last);
        if (s == 0) feed_throughput += std::abs(GetTotalComponentMolesIn(log, c, s, BoundaryLocation::BOTTOM, comp));
      }
      ASSERT_GT(feed_throughput, 0.0);
      // Measured |closure| <= 4.4e-8 mol = 1.5e-6 of the feed (trapezoids over
      // the logged points, worst after the pressure jump at the start of a
      // feed); bound 1e-5.
      for (const auto& [comp, e] : closure) EXPECT_NEAR(e, 0.0, 1e-5 * feed_throughput) << comp;
      // Elements: the generation makes and destroys none, to round-off.
      for (const char* element : {"C", "H", "N"}) {
        double g = 0.0;
        for (const auto& comp : o.components) {
          const auto* f = Components::GetElementFormula(comp);
          ASSERT_NE(f, nullptr);
          auto it = f->find(element);
          if (it != f->end()) g += it->second * log.ReactionGeneration(c, s, comp);
        }
        EXPECT_NEAR(g, 0.0, 1e-12 * feed_throughput) << element;
      }

      // The heat logged is -dH times the extent the generation logged
      // (both trapezoids over the same points of proportional fields).
      std::vector<double> heat(Qe.values.size(), 0.0);
      for (size_t k = 0; k < heat.size(); ++k)
        for (size_t i = 1; i + 1 < Qe.values[k].size(); ++i) heat[k] += Qe.values[k][i] * dom.volume[i];
      const double Q_rxn = Trapezoid(Qe.times, heat);
      const double extent = log.ReactionGeneration(c, s, IBUT);
      EXPECT_GT(extent, 0.0);
      EXPECT_NEAR(Q_rxn, -o.heat_of_reaction * extent, 1e-10 * Q_rxn);

      // First law of the bed: H_in - H_out + Q_rxn = dU, with
      // h = Cp (T - T_ref) and u = h - R T per mole of gas, so
      // U = sum dV [et (C Cp (T - T_ref) - P) + rho_b cs (T - T_ref)].
      auto U = [&](size_t k) {
        double u = 0.0;
        for (size_t i = 1; i + 1 < Te.values[k].size(); ++i)
          u += dom.volume[i] * (dom.et[i] * (Ce.values[k][i] * Cp * (Te.values[k][i] - T_ref) - Pe.values[k][i]) +
                                dom.bulk_density[i] * dom.Cps[i] * (Te.values[k][i] - T_ref));
        return u;
      };
      double H_in = 0.0, H_out = 0.0;  // into the bed, out of it [J]
      for (const auto& info : r.sim.Streams()) {
        if (info.step != s) continue;
        const StreamData* st = r.sim.Boundaries().GetStreamData(c, s, info.name);
        if (!st) continue;
        std::vector<double> h(st->time_series_.size());
        for (size_t k = 0; k < h.size(); ++k) {
          const double sign = info.location == BoundaryLocation::BOTTOM ? 1.0 : -1.0;  // + into the bed
          h[k] = sign * st->molar_flow_[k] * Cp * (st->temperature_boundary_[k] - T_ref);
        }
        const double H = Trapezoid(st->time_series_, h);
        (H > 0.0 ? H_in : H_out) += std::abs(H);
      }
      const double dU = U(last) - U(0);
      const double imbalance = H_in - H_out + Q_rxn - dU;
      if (s == 1) {
        // Sealed: exact first law up to the trapezoid of Q over the logged
        // points; measured 2.2e-5 of Q_rxn; bound 1e-4.
        EXPECT_EQ(H_in + H_out, 0.0);
        EXPECT_NEAR(imbalance, 0.0, 1e-4 * Q_rxn);
      } else {
        // Open: also the energy balance's first-order convection of T in
        // non-conservative form; measured 1.5e-5 to 1.4e-4 of H_in, growing
        // as the cycles steepen the temperature profile; bound 5e-4.
        EXPECT_GT(Q_rxn, 0.3 * H_in);  // the reaction heat is a large term
        EXPECT_NEAR(imbalance, 0.0, 5e-4 * H_in);
      }
    }
    // The per-step generation integrals sum to the cycle's component balance:
    // in + G = change of inventory over the cycle (same bound as per step,
    // two steps).
    for (const auto& comp : o.components)
      EXPECT_NEAR(in_cycle[comp] + gen_cycle[comp], inv1[comp] - inv0[comp], 2e-5 * feed_throughput) << comp;
    // What the reaction made of isobutane it took from butane, over the cycle.
    EXPECT_NEAR(gen_cycle[IBUT], -gen_cycle[BUT], 1e-12 * gen_cycle[IBUT]);
  }
}

// ---- (f) trace product ----------------------------------------------------------

// The product enters at y = 0 (bed full of carrier, none fed) behind a sharp
// reaction front (Da = 8): no logged value of y_P may undershoot -atol, the
// absolute tolerance of the mole-fraction unknowns, and the outlet carries the
// product of the conversion.
TEST(ReactionAnalytic, TraceProductStaysAboveMinusAtol) {
  for (auto scheme : {AdvectionScheme::WENO, AdvectionScheme::VanLeer, AdvectionScheme::Upwind}) {
    SCOPED_TRACE(static_cast<int>(scheme));
    Options o;
    o.scheme = scheme;
    o.k = 8.0 * 0.1 / ((1.0 - o.eb) * o.solid_density * o.length);
    o.step_time = 10.0;
    Reactor r;
    Build(r, o);
    const RunResult run = RunReactor(r);
    ASSERT_TRUE(run.Ok()) << run.message;
    const double atol = r.sim.Data().state.mole_fraction.at(IBUT).abstol;
    ASSERT_GT(atol, 0.0);
    const auto e = r.sim.Log().GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION, IBUT);
    ASSERT_GT(e.values.size(), 10u);
    double y_min = 1.0;
    for (const auto& row : e.values)
      for (size_t n = 1; n + 1 < row.size(); ++n) y_min = std::min(y_min, row[n]);
    // Measured minima: -5.1e-11 (WENO), -1.4e-11 (van Leer), ~0 (upwind),
    // against atol = 1e-8.
    EXPECT_GE(y_min, -atol);
    // sum(nu) = 0: y_P,out = y_A,in - y_A,out = y_A,in X. Each fraction is
    // its own unknown, so the identity holds to the integrator's global
    // error, of order its reltol (1e-6): measured differences with WENO
    // 4e-8 (Linux GCC) and 4.5e-7 (Windows MSVC, a different trajectory
    // through the same error control), < 1e-15 with van Leer and upwind;
    // bound reltol.
    const double reltol = r.sim.Data().settings.reltol.Get();
    ASSERT_GT(reltol, 0.0);
    const double yA = Outlet(LastRow(r.sim, ComponentVariableType::MOLE_FRACTION, BUT));
    const double yP = Outlet(LastRow(r.sim, ComponentVariableType::MOLE_FRACTION, IBUT));
    EXPECT_GT(yP, 0.099);  // X > 0.99 (1 - exp(-8) = 0.9997)
    EXPECT_NEAR(yP, o.y_feed.at(BUT) - yA, reltol);
  }
}

}  // namespace
}  // namespace adsim
