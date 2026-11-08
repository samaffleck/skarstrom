/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// Reaction KPIs (KPIConfig::reaction, KPIRecord::conversion ...) and the
// per-reaction series (Simulation::Series / Profile with a reaction): a packed
// catalyst bed (inert isotherms) built from scratch, fed continuously.
//
//  (a) first-order A -> B in plug flow: the KPI conversion equals the upwind
//      scheme's steady state 1 - (1 + Da/N)^-N and 1 - exp(-Da), equals
//      (in - out) / in from the streams, and is the same with A's KPI feed
//      endpoints and without (every endpoint's inflow);
//  (b) A -> B + C with a side reaction A -> D: selectivity and yield from the
//      logged generations, Y_B + Y_D = X, the space-time yield by hand on both
//      bases, recovery of a product not fed is NaN;
//  (c) the peak temperature is the maximum of the logged temperature;
//  (d) without reactions (or without a key reactant) there are no reaction
//      KPIs and the records are unchanged;
//  (e) registry round trip through a saved deck, and the cascade on
//      component removal;
//  (f) the reaction-rate series and profile, by name and by ID.
//
// Fast tier: 30-50 cells, short steps.

#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/objects/deck_builder.h"
#include "adsim/objects/deck_edit.h"
#include "adsim/objects/manifests.h"
#include "adsim/simulation.h"
#include "adsim/solver/kpis.h"

namespace adsim {
namespace {

const std::string N2 = "Nitrogen", BUT = "Butane", IBUT = "Isobutane", C2H6 = "Ethane", C2H4 = "Ethylene";

std::string Num(double v) {
  char buf[40];
  std::snprintf(buf, sizeof(buf), "%.17g", v);
  return buf;
}

struct Options {
  std::vector<std::string> components{N2, BUT, IBUT};
  // name -> stoichiometry; empty: no reaction
  std::vector<std::pair<std::string, std::map<std::string, double>>> reactions{};
  std::map<std::string, double> k{};  // name -> first-order pre-exponential factor [m3/(kg s)]
  std::map<std::string, double> y_feed{{N2, 0.9}, {BUT, 0.1}};
  double velocity = 0.1;  // [m/s]
  double length = 0.2;    // [m]
  double diameter = 0.05;  // [m]
  int cells = 50;
  double step_time = 10.0;  // [s]
  bool isothermal = true;
  double heat_of_reaction = 0.0;  // [J/mol], every reaction
  double temperature = 400.0;     // [K]
  double solid_density = 1000.0;  // [kg/m3]
  double solid_heat_capacity = 900.0;
  double eb = 0.4;
};

struct Reactor {
  Simulation sim;
  AdsorbentID cat = 0;
  size_t layer = 0, feed = 0, feed_in = 0, outlet = 0;
  std::map<std::string, ReactionID> rxn;
  Options o;

  double CatalystMass() const {
    return (1.0 - o.eb) * o.solid_density * M_PI * o.diameter * o.diameter / 4.0 * o.length;
  }
  double BedVolume() const { return M_PI * o.diameter * o.diameter / 4.0 * o.length; }
};

void Build(Reactor& r, const Options& o) {
  r.o = o;
  UserData& d = r.sim.Data();
  for (const auto& c : o.components) build::AddComponent(d, c);
  r.cat = build::AddAdsorbent(d, "Catalyst");
  build::SetBedShape(d, "cylinder");
  build::SetWallMaterial(d, "ss316l");
  build::SetWallBoundary(d, "zero_flux");
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
  for (const auto& [name, nu] : o.reactions) {
    r.rxn[name] = build::AddReaction(d, r.cat, name, nu, "power_law");
    build::SetReactionBasis(d, r.cat, r.rxn[name], "concentration");
  }
  r.sim.RefreshParameters();

  std::vector<std::pair<std::string, std::string>> values = {
      {"bed.diameter", Num(o.diameter)},
      {"adsorbent[Catalyst].solid_density", Num(o.solid_density)},
      {"adsorbent[Catalyst].heat_capacity", Num(o.solid_heat_capacity)},
      {"adsorbent[Catalyst].particle_void_fraction", "0.4"},
      {"layer[Bed].length", Num(o.length)},
      {"layer[Bed].cells", Num(o.cells)},
      {"layer[Bed].bed_void_fraction", Num(o.eb)},
      {"layer[Bed].particle_diameter", "5e-3"},
      {"initial.pressure", "1e5"},
      {"initial.temperature", Num(o.temperature)},
      {"initial.y[" + o.components.front() + "]", "1"},
      {"step[Feed].time", Num(o.step_time)},
      {"step[Feed].bottom[Feed inlet].velocity", Num(o.velocity)},
      {"step[Feed].bottom[Feed inlet].temperature", Num(o.temperature)},
      {"step[Feed].top.pressure", "1e5"},
  };
  for (const auto& [c, y] : o.y_feed) values.push_back({"step[Feed].bottom[Feed inlet].y[" + c + "]", Num(y)});
  for (const auto& [name, nu] : o.reactions) {
    const std::string head = "adsorbent[Catalyst].reaction[" + name + "].";
    values.push_back({head + "rate.pre_exp_factor", Num(o.k.at(name))});
    values.push_back({head + "heat_of_reaction", Num(o.heat_of_reaction)});
  }
  for (const auto& [path, text] : values) r.sim.SetParameterText(path, text);
  for (const auto& [name, nu] : o.reactions)
    r.sim.SetParameterText("adsorbent[Catalyst].reaction[" + name + "].heat_mode", "constant");

  d.settings.advection_scheme_ = AdvectionScheme::Upwind;
  d.settings.css.method_ = CSSMethod::NONE;
  d.bed.equations_.energy_balance_mode_ =
      o.isothermal ? EnergyBalanceMode::ISOTHERMAL : EnergyBalanceMode::NON_ISOTHERMAL;
}

RunResult RunReactor(Reactor& r, int cycles) {
  RunOptions options;
  options.cycles = cycles;
  options.on_message = [](LogLevel, const std::string&) {};
  return r.sim.Run(options);
}

const KPIRecord* Find(const std::vector<KPIRecord>& rows, const std::string& component) {
  for (const auto& k : rows)
    if (k.component == component) return &k;
  return nullptr;
}
// A pointer into a temporary would dangle.
const KPIRecord* Find(std::vector<KPIRecord>&& rows, const std::string& component) = delete;

double CycleGeneration(const Simulation& sim, int cycle, const std::string& c) {
  double g = 0.0;
  for (int s = 0; s < static_cast<int>(sim.Data().cycle.step_list_.size()); ++s)
    g += sim.Log().ReactionGeneration(cycle, s, c);
  return g;
}

// Da = rho_b k L / u = `da` at velocity u.
double FirstOrderK(const Options& o, double da) {
  return da * o.velocity / ((1.0 - o.eb) * o.solid_density * o.length);
}

// ---- (a) -----------------------------------------------------------------------

TEST(ReactionKPIs, FirstOrderPlugFlowConversion) {
  Options o;
  o.reactions = {{"R1", {{BUT, -1}, {IBUT, 1}}}};
  o.k = {{"R1", FirstOrderK(o, 1.0)}};
  Reactor r;
  Build(r, o);
  UserData& d = r.sim.Data();
  build::SetKPI(d, BUT, {r.feed_in}, {r.outlet});
  build::SetReactionKPI(d, BUT, {IBUT});
  r.sim.RefreshParameters();
  // Two cycles of one continuous feed step: the second has no start-up.
  const RunResult run = RunReactor(r, 2);
  ASSERT_TRUE(run.Ok()) << run.message;

  const auto rows = r.sim.KPIs(1);
  const KPIRecord* a = Find(rows, BUT);
  ASSERT_NE(a, nullptr);
  ASSERT_TRUE(a->reaction_kpis);
  const double N = o.cells, Da = 1.0;
  const double upwind = 1.0 - std::pow(1.0 + Da / N, -N);
  const double analytic = 1.0 - std::exp(-Da);
  EXPECT_NEAR(a->conversion, upwind, 1e-3 * upwind);
  EXPECT_NEAR(a->conversion, analytic, 0.02 * analytic);
  EXPECT_TRUE(std::isnan(a->selectivity));  // the key reactant's record carries the conversion only
  EXPECT_EQ(a->reaction_generation, CycleGeneration(r.sim, 1, BUT));

  // (in - out) / in from the streams.
  KPIEndpointRef in, out;
  in.endpoint_id = r.feed_in;
  out.endpoint_id = r.outlet;
  ASSERT_TRUE(in.Resolve(d.cycle.step_list_));
  ASSERT_TRUE(out.Resolve(d.cycle.step_list_));
  const double F_in = EndpointComponentMolesIn(r.sim.Boundaries(), 1, in, BUT);
  const double F_out = -EndpointComponentMolesIn(r.sim.Boundaries(), 1, out, BUT);
  EXPECT_NEAR(a->conversion, (F_in - F_out) / F_in, 2e-3 * a->conversion);
  EXPECT_NEAR(a->conversion, -a->reaction_generation / F_in, 1e-12);

  // The product's record: one mole of B per mole of A, so S = 1 and Y = X.
  const KPIRecord* b = Find(rows, IBUT);
  ASSERT_NE(b, nullptr);
  EXPECT_TRUE(b->reaction_kpis);
  EXPECT_FALSE(b->valid);  // its separation KPIs are not configured
  EXPECT_NEAR(b->selectivity, 1.0, 1e-9);
  EXPECT_NEAR(b->reaction_yield, a->conversion, 1e-9);
  EXPECT_TRUE(std::isnan(b->conversion));
  EXPECT_TRUE(std::isnan(b->purity));

  // Without A's KPI feed endpoints the feed is every endpoint's inflow of A:
  // here the one feed stream, so the same conversion.
  d.kpi_config.FindByComponent(BUT)->feed_endpoints.clear();
  d.kpi_config.FindByComponent(BUT)->tracked = false;
  const auto rows2 = r.sim.KPIs(1);
  const KPIRecord* a2 = Find(rows2, BUT);
  ASSERT_NE(a2, nullptr);
  EXPECT_NEAR(a2->conversion, a->conversion, 1e-12);
}

// ---- (b) -----------------------------------------------------------------------

Reactor& SideReactionReactor(Reactor& r, const std::string& basis) {
  Options o;
  o.components = {N2, BUT, C2H4, C2H6, IBUT};
  o.reactions = {{"Cracking", {{BUT, -1}, {C2H4, 1}, {C2H6, 1}}}, {"Isomerisation", {{BUT, -1}, {IBUT, 1}}}};
  o.k = {{"Cracking", FirstOrderK(o, 1.0)}, {"Isomerisation", FirstOrderK(o, 0.3)}};
  o.cells = 30;
  Build(r, o);
  UserData& d = r.sim.Data();
  // A's separation KPIs: fed at the bottom, the outlet is its "product".
  build::SetKPI(d, BUT, {r.feed_in}, {r.outlet});
  build::SetKPI(d, C2H4, {r.feed_in}, {r.outlet});  // tracked, not fed
  build::SetReactionKPI(d, BUT, {C2H4, IBUT}, basis);
  r.sim.RefreshParameters();
  const RunResult run = RunReactor(r, 2);
  EXPECT_TRUE(run.Ok()) << run.message;
  return r;
}

TEST(ReactionKPIs, SelectivityYieldAndSpaceTimeYield) {
  Reactor r;
  SideReactionReactor(r, "catalyst_mass");
  const auto rows = r.sim.KPIs(1);
  const KPIRecord* a = Find(rows, BUT);
  const KPIRecord* b = Find(rows, C2H4);
  const KPIRecord* dd = Find(rows, IBUT);
  ASSERT_TRUE(a && b && dd);
  EXPECT_EQ(Find(rows, C2H6), nullptr);  // neither tracked nor a reaction KPI component

  const double GA = CycleGeneration(r.sim, 1, BUT), GB = CycleGeneration(r.sim, 1, C2H4),
               GD = CycleGeneration(r.sim, 1, IBUT);
  ASSERT_LT(GA, 0.0);
  KPIEndpointRef in;
  in.endpoint_id = r.feed_in;
  ASSERT_TRUE(in.Resolve(r.sim.Data().cycle.step_list_));
  const double F = EndpointComponentMolesIn(r.sim.Boundaries(), 1, in, BUT);

  EXPECT_NEAR(a->conversion, -GA / F, 1e-12);
  EXPECT_NEAR(b->selectivity, GB / -GA, 1e-12);
  EXPECT_NEAR(dd->selectivity, GD / -GA, 1e-12);
  EXPECT_NEAR(b->reaction_yield, GB / F, 1e-12);
  EXPECT_NEAR(dd->reaction_yield, GD / F, 1e-12);
  // The rates are in the ratio of their constants (both first order in A).
  EXPECT_NEAR(b->selectivity, 1.0 / 1.3, 1e-6);
  // Every mole of A consumed makes one of B or one of D.
  EXPECT_LE(b->reaction_yield + dd->reaction_yield, a->conversion * (1.0 + 1e-9));
  EXPECT_NEAR(b->reaction_yield + dd->reaction_yield, a->conversion, 1e-9);
  EXPECT_NEAR(b->selectivity + dd->selectivity, 1.0, 1e-9);

  // Space-time yield by hand: G_P M_P / (t_cycle m_cat).
  const auto& mm = r.sim.Data().fluid.molar_mass_;
  const double t_cycle = r.o.step_time;
  EXPECT_NEAR(b->space_time_yield, GB * mm.at(C2H4) / (t_cycle * r.CatalystMass()), 1e-6 * b->space_time_yield);
  EXPECT_EQ(b->space_time_yield_unit, "kg/(kg s)");
  EXPECT_TRUE(std::isnan(a->space_time_yield));

  // Purity and recovery still apply: the product is not fed, so its recovery
  // does not apply (NaN), while its purity does.
  EXPECT_TRUE(b->valid);
  EXPECT_TRUE(std::isnan(b->recovery));
  EXPECT_GT(b->purity, 0.0);
  EXPECT_LT(b->purity, 1.0);
  EXPECT_TRUE(a->valid);
  EXPECT_NEAR(a->recovery, 1.0 - a->conversion, 2e-3);

  // The bed-volume basis.
  Reactor v;
  SideReactionReactor(v, "bed_volume");
  const auto rows_v = v.sim.KPIs(1);  // keep the records alive: Find returns a pointer into them
  const KPIRecord* bv = Find(rows_v, C2H4);
  ASSERT_NE(bv, nullptr);
  const double GBv = CycleGeneration(v.sim, 1, C2H4);
  EXPECT_NEAR(bv->space_time_yield, GBv * mm.at(C2H4) / (t_cycle * v.BedVolume()), 1e-6 * bv->space_time_yield);
  EXPECT_EQ(bv->space_time_yield_unit, "kg/(m3 s)");

  // The completed-cycle record (KPIHistory) is the one computed from the log.
  const auto& recorded = r.sim.KPIHistoryRecords().Cycle(1);
  const KPIRecord* br = Find(recorded, C2H4);
  ASSERT_NE(br, nullptr);
  EXPECT_DOUBLE_EQ(br->selectivity, b->selectivity);
  EXPECT_DOUBLE_EQ(br->space_time_yield, b->space_time_yield);
}

// ---- (c) -----------------------------------------------------------------------

TEST(ReactionKPIs, PeakTemperatureIsTheLoggedMaximum) {
  Options o;
  o.reactions = {{"R1", {{BUT, -1}, {IBUT, 1}}}};
  o.k = {{"R1", FirstOrderK(o, 2.0)}};
  o.y_feed = {{N2, 0.95}, {BUT, 0.05}};
  o.isothermal = false;
  o.heat_of_reaction = -1e5;
  o.cells = 30;
  o.step_time = 20.0;
  Reactor r;
  Build(r, o);
  build::SetReactionKPI(r.sim.Data(), BUT, {IBUT});
  r.sim.RefreshParameters();
  const RunResult run = RunReactor(r, 1);
  ASSERT_TRUE(run.Ok()) << run.message;

  const auto rows = r.sim.KPIs(0);  // keep the records alive: Find returns a pointer into them
  const KPIRecord* a = Find(rows, BUT);
  ASSERT_NE(a, nullptr);
  // Brute force over the logged temperature of the bed cells.
  const auto e = r.sim.Log().GetVariable(0, 0, VariableType::TEMPERATURE);
  const auto centres = r.sim.Positions();
  double peak = -1.0, z = 0.0, t = 0.0;
  for (size_t i = 0; i < e.values.size(); ++i)
    for (size_t c = 1; c + 1 < e.values[i].size(); ++c)
      if (e.values[i][c] > peak) {
        peak = e.values[i][c];
        z = centres[c - 1];
        t = e.times[i] - e.times.front();
      }
  EXPECT_GT(peak, o.temperature + 1.0);  // the exotherm heats the bed
  EXPECT_DOUBLE_EQ(a->peak_temperature, peak);
  EXPECT_DOUBLE_EQ(a->peak_temperature_position, z);
  EXPECT_DOUBLE_EQ(a->peak_temperature_time, t);
  const KPIRecord* b = Find(rows, IBUT);
  ASSERT_NE(b, nullptr);
  EXPECT_DOUBLE_EQ(b->peak_temperature, peak);
}

// ---- (d) -----------------------------------------------------------------------

void ExpectSameKPIs(const std::vector<KPIRecord>& x, const std::vector<KPIRecord>& y) {
  ASSERT_EQ(x.size(), y.size());
  auto same = [](double p, double q) { return (std::isnan(p) && std::isnan(q)) || p == q; };
  for (size_t i = 0; i < x.size(); ++i) {
    EXPECT_EQ(x[i].component, y[i].component);
    EXPECT_EQ(x[i].valid, y[i].valid);
    EXPECT_TRUE(same(x[i].purity, y[i].purity));
    EXPECT_TRUE(same(x[i].recovery, y[i].recovery));
    EXPECT_TRUE(same(x[i].productivity_mass, y[i].productivity_mass));
    EXPECT_TRUE(same(x[i].energy_consumption, y[i].energy_consumption));
    EXPECT_TRUE(same(x[i].mass_balance_error, y[i].mass_balance_error));
    EXPECT_FALSE(y[i].reaction_kpis);
    EXPECT_TRUE(std::isnan(y[i].conversion));
    EXPECT_TRUE(std::isnan(y[i].selectivity));
    EXPECT_TRUE(std::isnan(y[i].peak_temperature));
    EXPECT_TRUE(y[i].space_time_yield_unit.empty());
  }
}

TEST(ReactionKPIs, NoReactionsNoReactionKPIs) {
  // A deck without reactions: a configured key reactant changes nothing.
  Options o;
  o.cells = 20;
  o.step_time = 3.0;
  Reactor r;
  Build(r, o);
  build::SetKPI(r.sim.Data(), BUT, {r.feed_in}, {r.outlet});
  r.sim.RefreshParameters();
  ASSERT_TRUE(RunReactor(r, 1).Ok());
  EXPECT_TRUE(r.sim.Log().Reactions().empty());
  const auto plain = r.sim.KPIs(0);
  ASSERT_EQ(plain.size(), 1u);
  build::SetReactionKPI(r.sim.Data(), BUT, {IBUT});
  ExpectSameKPIs(plain, r.sim.KPIs(0));

  // A deck with reactions but no key reactant: no reaction KPIs either.
  Options ro = o;
  ro.reactions = {{"R1", {{BUT, -1}, {IBUT, 1}}}};
  ro.k = {{"R1", FirstOrderK(ro, 1.0)}};
  Reactor rr;
  Build(rr, ro);
  build::SetKPI(rr.sim.Data(), BUT, {rr.feed_in}, {rr.outlet});
  rr.sim.RefreshParameters();
  ASSERT_TRUE(RunReactor(rr, 1).Ok());
  const auto rows = rr.sim.KPIs(0);
  ASSERT_EQ(rows.size(), 1u);
  EXPECT_FALSE(rows[0].reaction_kpis);
  EXPECT_TRUE(std::isnan(rows[0].conversion));
}

// ---- (e) -----------------------------------------------------------------------

TEST(ReactionKPIs, RegistryRoundTripAndComponentRemoval) {
  Options o;
  o.components = {N2, BUT, IBUT, C2H4};
  Reactor r;
  Build(r, o);
  Simulation& sim = r.sim;
  EXPECT_EQ(sim.GetParameter("kpi.reaction.key_reactant"), 0.0);  // "none"
  sim.SetParameterText("kpi.reaction.key_reactant", BUT);
  sim.SetParameterText("kpi.reaction.product[" + IBUT + "]", "true");
  sim.SetParameterText("kpi.reaction.product[" + C2H4 + "]", "true");
  sim.SetParameterText("kpi.reaction.basis", "bed_volume");
  EXPECT_THROW(sim.SetParameterText("kpi.reaction.key_reactant", "Unobtainium"), std::invalid_argument);
  EXPECT_EQ(sim.Data().kpi_config.reaction.key_reactant, BUT);
  EXPECT_THROW(build::SetReactionKPI(sim.Data(), BUT, {IBUT}, "per_hour"), std::invalid_argument);
  EXPECT_THROW(build::SetReactionKPI(sim.Data(), BUT, {BUT}), std::invalid_argument);

  // Saved and reloaded.
  Simulation copy = Simulation::FromJson(sim.ToJson());
  const ReactionKPIConfig& rx = copy.Data().kpi_config.reaction;
  EXPECT_EQ(rx.key_reactant, BUT);
  EXPECT_EQ(rx.products, (std::vector<std::string>{IBUT, C2H4}));
  EXPECT_EQ(rx.basis, SpaceTimeYieldBasis::BED_VOLUME);
  EXPECT_EQ(copy.GetParameter("kpi.reaction.key_reactant"), 2.0);  // none, Nitrogen, Butane
  EXPECT_EQ(copy.GetParameter("kpi.reaction.product[" + IBUT + "]"), 1.0);
  EXPECT_EQ(copy.GetParameter("kpi.reaction.product[" + N2 + "]"), 0.0);
  EXPECT_EQ(copy.GetParameter("kpi.reaction.basis"), 1.0);

  // A removed product leaves the list; a removed key reactant clears the KPIs.
  std::string err;
  ASSERT_TRUE(RemoveComponentFromDeck(copy.Data(), C2H4, err)) << err;
  EXPECT_EQ(copy.Data().kpi_config.reaction.products, (std::vector<std::string>{IBUT}));
  ASSERT_TRUE(RemoveComponentFromDeck(copy.Data(), BUT, err)) << err;
  EXPECT_FALSE(copy.Data().kpi_config.reaction.IsSet());
  copy.RefreshParameters();
  EXPECT_EQ(copy.GetParameter("kpi.reaction.key_reactant"), 0.0);

  // Clearing.
  build::SetReactionKPI(sim.Data(), "", {});
  EXPECT_FALSE(sim.Data().kpi_config.reaction.IsSet());
  EXPECT_TRUE(sim.Data().kpi_config.reaction.products.empty());

  // A deck written before the reaction KPIs loads without them.
  Simulation old = Simulation::FromFile("../examples/decks/Perez2019_A.json");
  EXPECT_FALSE(old.Data().kpi_config.reaction.IsSet());
}

TEST(ReactionKPIs, ManifestListsTheReactionKPIs) {
  std::vector<std::string> keys;
  for (const auto& k : KpiDefinitions())
    if (k.reaction) keys.push_back(k.key);
  EXPECT_EQ(keys, (std::vector<std::string>{"conversion", "selectivity", "reaction_yield", "space_time_yield",
                                            "peak_temperature", "peak_temperature_position",
                                            "peak_temperature_time", "reaction_generation"}));
}

// ---- (f) -----------------------------------------------------------------------

TEST(ReactionKPIs, ReactionRateSeriesAndProfile) {
  Options o;
  o.reactions = {{"R1", {{BUT, -1}, {IBUT, 1}}}};
  o.k = {{"R1", FirstOrderK(o, 1.0)}};
  o.cells = 20;
  o.step_time = 5.0;
  Reactor r;
  Build(r, o);
  ASSERT_TRUE(RunReactor(r, 1).Ok());
  const Simulation& sim = r.sim;

  const auto reactions = sim.SeriesReactions();
  ASSERT_EQ(reactions.size(), 1u);
  EXPECT_EQ(reactions[0].name, "R1");
  EXPECT_EQ(reactions[0].adsorbent, "Catalyst");
  EXPECT_EQ(reactions[0].variable_name, "R1 Rate (Catalyst)");
  EXPECT_EQ(ReactionVariableNames(), std::vector<std::string>{"reaction_rate"});

  const auto log = sim.Log().GetReactionRate(0, 0, r.cat, r.rxn.at("R1"));
  const AxialProfile p = sim.Profile("reaction_rate", "", 0, 0, "", "R1");
  EXPECT_EQ(p.reaction, "R1");
  EXPECT_EQ(p.adsorbent, "Catalyst");
  EXPECT_EQ(p.unit, "mol/(kg.s)");
  ASSERT_EQ(p.rows, log.values.size());
  ASSERT_EQ(p.cols, static_cast<size_t>(o.cells));
  EXPECT_EQ(p.values[(p.rows - 1) * p.cols + 3], log.values.back()[4]);
  // By ID, and with the adsorbent named.
  const AxialProfile by_id = sim.Profile("reaction_rate", "", 0, 0, "Catalyst", std::to_string(r.rxn.at("R1")));
  EXPECT_EQ(by_id.values, p.values);

  const TimeSeries top = sim.Series("reaction_rate", "", SeriesLocation::Top, 0.0, 0, "", "R1");
  const TimeSeries bottom = sim.Series("reaction_rate", "", SeriesLocation::Bottom, 0.0, 0, "", "R1");
  ASSERT_EQ(top.values.size(), log.values.size());
  EXPECT_EQ(top.values.back(), log.values.back()[o.cells]);
  EXPECT_EQ(bottom.values.back(), log.values.back()[1]);
  EXPECT_GT(bottom.values.back(), top.values.back());  // A is consumed along the bed
  EXPECT_DOUBLE_EQ(top.time.front(), 0.0);
  // Midway between two cell centres: their mean.
  const auto z = sim.Positions();
  const TimeSeries mid = sim.Series("reaction_rate", "", SeriesLocation::Position, 0.5 * (z[4] + z[5]), 0, "", "R1");
  EXPECT_NEAR(mid.values.back(), 0.5 * (log.values.back()[5] + log.values.back()[6]), 1e-12 * mid.values.back());

  // The reactions' totals stay reachable the usual way.
  EXPECT_FALSE(sim.Series("reaction_heat").values.empty());
  EXPECT_FALSE(sim.Profile("reaction_source", BUT, 0, 0).values.empty());

  // Misuse.
  EXPECT_THROW(sim.Series("reaction_rate"), std::invalid_argument);                        // no reaction
  EXPECT_THROW(sim.Series("pressure", "", SeriesLocation::Top, 0.0, -1, "", "R1"), std::invalid_argument);
  EXPECT_THROW(sim.Profile("reaction_rate", "", 0, 0, "", "R9"), std::invalid_argument);   // unknown
  EXPECT_THROW(sim.Profile("reaction_rate", BUT, 0, 0, "", "R1"), std::invalid_argument);  // component
}

}  // namespace
}  // namespace adsim
