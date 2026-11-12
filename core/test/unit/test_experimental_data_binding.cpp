// test_experimental_data_binding.cpp
//
// What a stored dataset measures (ObservableBinding, experimental_data.h):
// saving and loading it with the deck (older decks load unbound), checking it
// (NormaliseObservableBinding and the deck_edit.h edits), the deck view, and
// the simulated counterpart and comparison of parameterisation.h
// (SimulatedSeries, CompareWithSimulation) on a breakthrough deck (Ward 2022,
// outlet CO2 and a thermocouple) and on a cyclic deck with per-cycle KPIs.
// The runs are cut to 20 cells and 3 cycles: they check the plumbing, not the
// decks' results.

#include <cmath>
#include <filesystem>
#include <sstream>
#include <string>
#include <vector>

#include <cereal/archives/json.hpp>
#include <cereal/external/rapidjson/document.h>
#include <gtest/gtest.h>

#include "adsim/objects/deck_edit.h"
#include "adsim/objects/deck_view.h"
#include "adsim/simulation.h"
#include "adsim/solver/parameterisation.h"

namespace adsim {
namespace {

const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";
std::string Deck(const std::string& name) { return (kExamples / name).string(); }

constexpr int kShortCells = 20;

XYData Points(std::vector<double> x, std::vector<double> y) {
  XYData d;
  d.name = "data";
  d.x = std::move(x);
  d.y = std::move(y);
  return d;
}

ObservableBinding Bind(ObservableKind kind, std::string variable, std::string component = {}) {
  ObservableBinding b;
  b.kind = kind;
  b.variable = std::move(variable);
  b.component = std::move(component);
  return b;
}

std::string SaveJson(const ExperimentalData& data) {
  std::ostringstream oss;
  {
    cereal::JSONOutputArchive ar(oss);
    ar(cereal::make_nvp("Experimental_Data", data));
  }
  return oss.str();
}

ExperimentalData LoadJson(const std::string& text) {
  ExperimentalData out;
  std::istringstream iss(text);
  cereal::JSONInputArchive ar(iss);
  ar(cereal::make_nvp("Experimental_Data", out));
  return out;
}

TEST(ExperimentalDataBinding, RoundTripsThroughTheDeck) {
  ExperimentalData data;
  XYData xy = Points({0, 60, 120}, {20, 25, 30});
  xy.observable = Bind(ObservableKind::Bulk, "temperature");
  xy.observable.location = ObservableLocation::Position;
  xy.observable.z = 0.0715;
  xy.observable.cycle = 2;
  xy.observable.x_unit = "min";
  xy.observable.x_scale = 60.0;
  xy.observable.y_unit = "degC";
  xy.observable.uncertainty = {0.5, 0.5, 1.0};
  xy.observable.weight = 2.0;
  data.GetAllData().push_back(xy);
  XYData kpi = Points({1, 2}, {0.9, 0.95});
  kpi.observable = Bind(ObservableKind::Kpi, "purity", "CarbonDioxide");
  kpi.observable.cycle_offset = -1;
  data.GetAllData().push_back(kpi);

  const ExperimentalData back = LoadJson(SaveJson(data));
  ASSERT_EQ(back.size(), 2u);
  const ObservableBinding& b = back.GetAllData()[0].observable;
  EXPECT_EQ(b.kind, ObservableKind::Bulk);
  EXPECT_EQ(b.variable, "temperature");
  EXPECT_EQ(b.location, ObservableLocation::Position);
  EXPECT_DOUBLE_EQ(b.z, 0.0715);
  EXPECT_EQ(b.cycle, 2);
  EXPECT_EQ(b.x_unit, "min");
  EXPECT_DOUBLE_EQ(b.x_scale, 60.0);
  EXPECT_EQ(b.y_unit, "degC");
  EXPECT_EQ(b.uncertainty, (std::vector<double>{0.5, 0.5, 1.0}));
  EXPECT_DOUBLE_EQ(b.weight, 2.0);
  const ObservableBinding& k = back.GetAllData()[1].observable;
  EXPECT_EQ(k.kind, ObservableKind::Kpi);
  EXPECT_EQ(k.component, "CarbonDioxide");
  EXPECT_EQ(k.cycle_offset, -1);
  EXPECT_EQ(k.XAxis(), ObservableXAxis::Cycle);

  // The tags, not the enum values, are stored.
  rapidjson::Document doc;
  doc.Parse(SaveJson(data).c_str());
  const auto& first = doc["Experimental_Data"]["Data"][0];
  EXPECT_EQ(first["cereal_class_version"].GetInt(), 1);
  EXPECT_STREQ(first["Observable"]["Kind"].GetString(), "bulk");
  EXPECT_STREQ(first["Observable"]["Location"].GetString(), "position");
}

TEST(ExperimentalDataBinding, DatasetsSavedBeforeBindingsLoadUnbound) {
  // An XYData of version 0: no Observable member.
  const std::string old = R"({"Experimental_Data": {"cereal_class_version": 1, "Data": [
      {"cereal_class_version": 0, "Name": "T", "x_data": [0.0, 1.0], "y_data": [300.0, 301.0],
       "x_label": "Time [s]", "y_label": "T [K]"}]}})";
  const ExperimentalData data = LoadJson(old);
  ASSERT_EQ(data.size(), 1u);
  const XYData& xy = data.GetAllData()[0];
  EXPECT_EQ(xy.name, "T");
  EXPECT_FALSE(xy.observable.Bound());
  EXPECT_EQ(xy.observable.x_unit, "s");
  EXPECT_DOUBLE_EQ(xy.observable.weight, 1.0);
}

TEST(ExperimentalDataBinding, NormaliseChecksTheBinding) {
  const std::vector<std::string> comps = {"CarbonDioxide", "Nitrogen"};
  const XYData xy = Points({0, 1, 2}, {0.1, 0.2, 0.3});
  std::string err;
  auto ok = [&](ObservableBinding b) { return NormaliseObservableBinding(b, xy, comps, err); };

  EXPECT_TRUE(ok(ObservableBinding{}));  // unbound is fine
  EXPECT_TRUE(ok(Bind(ObservableKind::Component, "mole_fraction", "CarbonDioxide")));
  EXPECT_FALSE(ok(Bind(ObservableKind::Component, "mole_fraction", "Helium")));
  EXPECT_NE(err.find("no component"), std::string::npos);
  EXPECT_FALSE(ok(Bind(ObservableKind::Component, "mole_fractions", "CarbonDioxide")));
  EXPECT_FALSE(ok(Bind(ObservableKind::Bulk, "temperature", "CarbonDioxide")));  // bulk takes no component
  EXPECT_FALSE(ok(Bind(ObservableKind::Bulk, "mole_fraction")));
  EXPECT_FALSE(ok(Bind(ObservableKind::Kpi, "purity")));                          // needs a component
  EXPECT_FALSE(ok(Bind(ObservableKind::Kpi, "purityy", "CarbonDioxide")));

  // x unit: a time unit, its multiplier stored; a KPI's x is the cycle.
  ObservableBinding t = Bind(ObservableKind::Bulk, "temperature");
  t.x_unit = "min";
  ASSERT_TRUE(NormaliseObservableBinding(t, xy, comps, err)) << err;
  EXPECT_DOUBLE_EQ(t.x_scale, 60.0);
  t.x_unit = "h";
  ASSERT_TRUE(NormaliseObservableBinding(t, xy, comps, err)) << err;
  EXPECT_DOUBLE_EQ(t.x_scale, 3600.0);
  t.x_unit = "ms";
  ASSERT_TRUE(NormaliseObservableBinding(t, xy, comps, err)) << err;
  EXPECT_DOUBLE_EQ(t.x_scale, 1e-3);
  t.x_unit = "m";
  EXPECT_FALSE(NormaliseObservableBinding(t, xy, comps, err));
  ObservableBinding k = Bind(ObservableKind::Kpi, "recovery", "CarbonDioxide");
  k.x_unit = "min";
  k.x_scale = 7;
  EXPECT_FALSE(NormaliseObservableBinding(k, xy, comps, err));
  k.x_unit = "";
  ASSERT_TRUE(NormaliseObservableBinding(k, xy, comps, err)) << err;
  EXPECT_EQ(k.x_unit, "cycle");
  EXPECT_DOUBLE_EQ(k.x_scale, 1.0);

  // y unit: of the observable's dimension.
  ObservableBinding y = Bind(ObservableKind::Bulk, "temperature");
  for (const char* u : {"K", "degC", "degF"}) {
    y.y_unit = u;
    EXPECT_TRUE(ok(y)) << u << ": " << err;
  }
  y.y_unit = "kPa";
  EXPECT_FALSE(ok(y));
  y.y_unit = "furlong";
  EXPECT_FALSE(ok(y));
  ObservableBinding m = Bind(ObservableKind::Component, "mole_fraction", "Nitrogen");
  m.y_unit = "%";
  EXPECT_TRUE(ok(m)) << err;
  ObservableBinding e = Bind(ObservableKind::Kpi, "energy_consumption", "CarbonDioxide");
  e.y_unit = "kWh/kg";
  EXPECT_TRUE(ok(e)) << err;

  // Position, uncertainty, weight.
  ObservableBinding p = Bind(ObservableKind::Bulk, "pressure");
  p.location = ObservableLocation::Position;
  p.z = -0.1;
  EXPECT_FALSE(ok(p));
  p.z = 0.1;
  p.uncertainty = {1, 1};
  EXPECT_FALSE(ok(p));
  p.uncertainty = {1, 0, 1};
  EXPECT_FALSE(ok(p));
  p.uncertainty = {1, 2, 3};
  EXPECT_TRUE(ok(p)) << err;
  p.weight = 0;
  EXPECT_FALSE(ok(p));
}

TEST(ExperimentalDataBinding, EditsKeepTheDeckValidAndTheViewCarriesThem) {
  UserData data;
  LoadFromJson(data, Deck("Ward2022_Q100ccm.json"));
  const size_t n0 = data.exp_data_.size();
  ASSERT_GT(n0, 0u);

  XYData xy = Points({0, 1}, {0, 1});
  xy.name = "outlet";
  xy.observable = Bind(ObservableKind::Component, "mole_fraction", "Nitrogen");  // not in the deck
  size_t index = 0;
  std::string err;
  EXPECT_FALSE(AddExperimentalData(data, xy, index, err));
  EXPECT_EQ(data.exp_data_.size(), n0);
  xy.observable.component = "CarbonDioxide";
  xy.observable.x_unit = "min";
  ASSERT_TRUE(AddExperimentalData(data, xy, index, err)) << err;
  EXPECT_EQ(index, n0);
  EXPECT_DOUBLE_EQ(data.exp_data_.GetAllData()[index].observable.x_scale, 60.0);
  EXPECT_FALSE(AddExperimentalData(data, Points({0, 1}, {0}), index, err));

  ObservableBinding bad = Bind(ObservableKind::Bulk, "temperature");
  bad.y_unit = "Pa";
  EXPECT_FALSE(SetExperimentalDataBinding(data, n0, bad, err));
  EXPECT_EQ(data.exp_data_.GetAllData()[n0].observable.variable, "mole_fraction");  // unchanged
  EXPECT_FALSE(SetExperimentalDataBinding(data, n0 + 5, ObservableBinding{}, err));
  ObservableBinding kpi = Bind(ObservableKind::Kpi, "purity", "CarbonDioxide");
  kpi.cycle_offset = -1;
  ASSERT_TRUE(SetExperimentalDataBinding(data, n0, kpi, err)) << err;

  rapidjson::Document view;
  view.Parse(DeckView(data).c_str());
  ASSERT_TRUE(view.HasMember("experimental_data"));
  const auto& sets = view["experimental_data"];
  ASSERT_EQ(sets.Size(), n0 + 1);
  const auto& last = sets[static_cast<rapidjson::SizeType>(n0)];
  EXPECT_EQ(last["index"].GetInt(), static_cast<int>(n0));
  EXPECT_STREQ(last["name"].GetString(), "outlet");
  EXPECT_EQ(last["x"].Size(), 2u);
  EXPECT_STREQ(last["observable"]["kind"].GetString(), "kpi");
  EXPECT_STREQ(last["observable"]["x_axis"].GetString(), "cycle");
  EXPECT_STREQ(last["observable"]["x_unit"].GetString(), "cycle");
  EXPECT_STREQ(last["observable"]["si_unit"].GetString(), "mol/mol");
  EXPECT_EQ(last["observable"]["cycle_offset"].GetInt(), -1);
  EXPECT_TRUE(last["observable"]["uncertainty"].IsNull());

  ASSERT_TRUE(RemoveExperimentalData(data, n0, err)) << err;
  EXPECT_EQ(data.exp_data_.size(), n0);
  EXPECT_FALSE(RemoveExperimentalData(data, n0, err));
}

// The Ward 2022 breakthrough (CO2 into a He-filled activated-carbon bed) on
// 20 cells: the outlet CO2 and the thermocouple at z = 0.0715 m.
class WardBinding : public ::testing::Test {
 protected:
  static void SetUpTestSuite() {
    sim_ = new Simulation(Simulation::FromFile(Deck("Ward2022_Q100ccm.json")));
    for (auto& layer : sim_->Data().bed.layers_.layers_) layer.domain_.n_cells = kShortCells;
    ASSERT_TRUE(sim_->Run().Ok());
  }
  static void TearDownTestSuite() {
    delete sim_;
    sim_ = nullptr;
  }
  static SimulatedCounterpart Simulated(const XYData& xy) {
    return SimulatedSeries(xy, sim_->Log(), sim_->KPIs(),
                           BedBoundarySources{&sim_->Boundaries(), &sim_->Data().cycle.step_list_});
  }
  static XYData Named(const std::string& prefix) {
    for (const auto& xy : sim_->Data().exp_data_.GetAllData())
      if (xy.name.rfind(prefix, 0) == 0) return xy;
    ADD_FAILURE() << "no dataset " << prefix;
    return {};
  }
  static Simulation* sim_;
};
Simulation* WardBinding::sim_ = nullptr;

TEST_F(WardBinding, OutletMoleFractionMatchesTheSeries) {
  XYData xy = Named("outlet_CO2_mole_fraction");
  xy.observable = Bind(ObservableKind::Component, "mole_fraction", "CarbonDioxide");
  const SimulatedCounterpart s = Simulated(xy);
  const TimeSeries ref = sim_->Series("mole_fraction", "CarbonDioxide", SeriesLocation::Top);
  EXPECT_EQ(s.x, ref.time);
  EXPECT_EQ(s.y, ref.values);
  EXPECT_EQ(s.si_unit, "mol/mol");

  const DatasetComparison c = CompareWithSimulation(xy, s);
  EXPECT_GT(c.points, xy.x.size() / 2);
  EXPECT_LE(c.points, xy.x.size());
  EXPECT_TRUE(std::isfinite(c.rmse));
  EXPECT_LT(c.nrmse, 0.3);  // a coarse grid, but the same breakthrough
  EXPECT_GE(c.max_abs_error, c.rmse);
  EXPECT_DOUBLE_EQ(c.weighted_sse, c.rmse * c.rmse * static_cast<double>(c.points));  // sigma = 1, weight 1

  // The data in % and minutes: the same comparison, scaled.
  XYData pct = xy;
  for (double& x : pct.x) x /= 60.0;
  for (double& y : pct.y) y *= 100.0;
  pct.observable.x_unit = "min";
  pct.observable.y_unit = "%";
  std::string err;
  ASSERT_TRUE(NormaliseObservableBinding(pct.observable, pct, sim_->Data().fluid.components_, err)) << err;
  const SimulatedCounterpart sp = Simulated(pct);
  ASSERT_EQ(sp.x.size(), s.x.size());
  EXPECT_NEAR(sp.x.back(), s.x.back() / 60.0, 1e-9);
  EXPECT_NEAR(sp.y.back(), s.y.back() * 100.0, 1e-9);
  const DatasetComparison cp = CompareWithSimulation(pct, sp);
  EXPECT_EQ(cp.points, c.points);
  EXPECT_NEAR(cp.rmse, 100.0 * c.rmse, 1e-9 * 100.0);
  EXPECT_NEAR(cp.nrmse, c.nrmse, 1e-9);

  // Uncertainty and weight enter only the weighted sum.
  pct.observable.uncertainty.assign(pct.y.size(), 2.0);
  pct.observable.weight = 3.0;
  const DatasetComparison cw = CompareWithSimulation(pct, sp);
  EXPECT_NEAR(cw.rmse, cp.rmse, 1e-12);
  EXPECT_NEAR(cw.weighted_sse, 3.0 * cp.weighted_sse / 4.0, 1e-9 * cp.weighted_sse);
}

TEST_F(WardBinding, ThermocoupleInDegreesCelsius) {
  XYData xy = Named("temperature[K]_z=0.0715");
  xy.observable = Bind(ObservableKind::Bulk, "temperature");
  xy.observable.location = ObservableLocation::Position;
  xy.observable.z = 0.0715;
  const SimulatedCounterpart k = Simulated(xy);
  const TimeSeries ref = sim_->Series("temperature", "", SeriesLocation::Position, 0.0715);
  EXPECT_EQ(k.y, ref.values);
  const DatasetComparison ck = CompareWithSimulation(xy, k);
  EXPECT_GT(ck.points, 0u);
  EXPECT_LT(ck.rmse, 10.0);  // [K]

  XYData c = xy;
  for (double& y : c.y) y -= 273.15;
  c.observable.y_unit = "degC";
  const SimulatedCounterpart sc = Simulated(c);
  ASSERT_EQ(sc.y.size(), k.y.size());
  EXPECT_NEAR(sc.y.front(), k.y.front() - 273.15, 1e-9);
  const DatasetComparison cc = CompareWithSimulation(c, sc);
  EXPECT_EQ(cc.points, ck.points);
  EXPECT_NEAR(cc.rmse, ck.rmse, 1e-8);
}

TEST_F(WardBinding, UnboundOrUnloggedIsAnError) {
  XYData xy = Named("outlet_CO2_mole_fraction");
  xy.observable = ObservableBinding{};
  EXPECT_THROW(Simulated(xy), std::invalid_argument);
  xy.observable = Bind(ObservableKind::Component, "mole_fraction", "CarbonDioxide");
  xy.observable.cycle = 7;  // one cycle was run
  EXPECT_THROW(Simulated(xy), std::out_of_range);
  xy.observable = Bind(ObservableKind::Kpi, "purity", "CarbonDioxide");  // Ward tracks no KPIs
  EXPECT_THROW(Simulated(xy), std::out_of_range);
}

// Perez 2019 A (CO2 VSA, CO2 tracked): three cycles on 20 cells. KPIs come
// from every completed cycle, also the ones evicted from the log.
TEST(ExperimentalDataBinding, PerCycleKpisMatchedByCycleNumber) {
  Simulation sim = Simulation::FromFile(Deck("Perez2019_A.json"));
  sim.SetParameter("cycle.number_of_cycles", 3);
  for (auto& layer : sim.Data().bed.layers_.layers_) layer.domain_.n_cells = kShortCells;
  ASSERT_TRUE(sim.Run().Ok());
  const std::vector<KPIRecord> kpis = sim.KPIs();

  // 1-based cycle numbers 1..3 (and a cycle the run does not have), in %.
  XYData xy = Points({1, 2, 3, 9}, {0, 0, 0, 50});
  for (const auto& k : kpis)
    if (k.component == "CarbonDioxide") xy.y[static_cast<size_t>(k.cycle)] = 100.0 * k.purity + 1.0;
  xy.observable = Bind(ObservableKind::Kpi, "purity", "CarbonDioxide");
  xy.observable.cycle_offset = -1;
  xy.observable.y_unit = "%";
  std::string err;
  ASSERT_TRUE(NormaliseObservableBinding(xy.observable, xy, sim.Data().fluid.components_, err)) << err;

  const SimulatedCounterpart s =
      SimulatedSeries(xy, sim.Log(), kpis, BedBoundarySources{&sim.Boundaries(), &sim.Data().cycle.step_list_});
  EXPECT_EQ(s.x, (std::vector<double>{1, 2, 3}));
  ASSERT_EQ(s.y.size(), 3u);
  EXPECT_NEAR(s.y[0], xy.y[0] - 1.0, 1e-9);
  EXPECT_EQ(s.si_unit, "mol/mol");

  const DatasetComparison c = CompareWithSimulation(xy, s);
  EXPECT_EQ(c.points, 3u);  // cycle 9 not run
  EXPECT_NEAR(c.rmse, 1.0, 1e-9);
  EXPECT_NEAR(c.max_abs_error, 1.0, 1e-9);
}

}  // namespace
}  // namespace adsim
