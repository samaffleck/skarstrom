"""Tests for building a simulation from scratch (skarstrom.build and the builder methods)."""

from pathlib import Path

import numpy as np
import pytest

import skarstrom as sk
from deck_compare import connections, diff, normalise, parameter_diff

RUN_CHECKED = ["Perez2019_A", "Perez2019_E", "Ward2022_Q100ccm", "Illustrative_TSA_13X", "Arora2021_LiLSX_PSA"]


def assert_same_deck(a, b):
    """The same structure and values: deck JSON, registry values by alias, connections."""
    assert diff(normalise(a.to_json()), normalise(b.to_json())) == []
    assert parameter_diff(a, b) == []
    assert connections(a) == connections(b)


def assert_same_run(a, b, cycles):
    ra, rb = a.run(cycles=cycles), b.run(cycles=cycles)
    assert ra.ok and rb.ok, (ra.message, rb.message)
    assert ra.cycles_run == rb.cycles_run
    ka, kb = a.kpis(), b.kpis()
    assert len(ka) == len(kb)
    for x, y in zip(ka, kb):
        dx, dy = x.as_dict(), y.as_dict()
        assert dx.keys() == dy.keys()
        for k in dx:
            if isinstance(dx[k], float):
                np.testing.assert_allclose(dy[k], dx[k], rtol=1e-9, atol=0, err_msg=k)
            else:
                assert dx[k] == dy[k], k
    for var, comps in (("pressure", [None]), ("temperature", [None]), ("mole_fraction", a.components)):
        for c in comps:
            for loc in ("top", "bottom"):
                sa, sb = a.series(var, c, location=loc), b.series(var, c, location=loc)
                np.testing.assert_allclose(sb.time, sa.time, rtol=1e-12)
                np.testing.assert_allclose(sb.values, sa.values, rtol=1e-9, atol=1e-12 * np.max(np.abs(sa.values)))
    assert [(s["step"], s["name"]) for s in a.streams()] == [(s["step"], s["name"]) for s in b.streams()]
    for s in a.streams():
        fa = a.stream_series(s["step"], s["name"], "molar_flow_rate")
        fb = b.stream_series(s["step"], s["name"], "molar_flow_rate")
        np.testing.assert_allclose(fb.values, fa.values, rtol=1e-9, atol=1e-15)


# ---- kinds ------------------------------------------------------------------------------------


def test_build_kinds():
    kinds = sk.build_kinds()
    for family in ("isotherm", "rate_model", "rate_coefficient", "momentum", "dispersion", "heat_dispersion",
                   "internal_heat_transfer", "external_heat_transfer", "particle_heat_transfer",
                   "thermal_equilibrium", "bed_shape", "wall_material", "wall_boundary", "profile", "flow",
                   "pump", "pump_energy", "pump_efficiency"):
        assert kinds[family], family
    assert set(sk.isotherm_kinds()) <= set(kinds["isotherm"]) and "inert" in kinds["isotherm"]
    assert kinds["momentum"] == ["kozeny_carman", "ergun", "blake_kozeny"]
    assert kinds["flow"] == ["velocity", "volumetric_flow", "molar_flow", "mass_flow"]


def test_unknown_kinds_list_the_valid_ones():
    sim = sk.Simulation.new()
    sim.add_component("Nitrogen")
    ads = sim.add_adsorbent("5A")
    with pytest.raises(ValueError, match="langmuir"):
        ads.set_isotherm("not_an_isotherm", components=["Nitrogen"])
    with pytest.raises(ValueError, match="ldf2"):
        ads.set_rate_model("Nitrogen", "fickian")
    with pytest.raises(ValueError, match="cylinder"):
        sim.bed.set_shape("sphere", diameter=0.1)
    with pytest.raises(ValueError, match="ergun"):
        sim.add_layer("bed", composition={"5A": 1.0}, momentum="darcy")
    layer = sim.add_layer("bed", composition={"5A": 1.0})
    with pytest.raises(ValueError, match="yagi_kunii"):
        layer.set_heat_dispersion("magic")
    step = sim.add_step("Feed", time=10)
    with pytest.raises(ValueError, match="molar_flow"):
        step.bottom.set_streams("inlet", [sk.Stream("feed", flow=sk.Flow("litres", sk.constant(1.0)))])
    with pytest.raises(ValueError, match="inlet"):
        step.top.set_pressure("sideways", 1e5)
    with pytest.raises(ValueError, match="bed pressure"):
        step.top.set_pressure("outlet", sk.linear(1e5, start=2e5))
    with pytest.raises(ValueError, match="pressure"):
        sim.add_variable("p_feed", 1e5, category="pressures")
    with pytest.raises(ValueError, match="built-in"):
        sim.add_variable("T", 300.0, category="temperature")
    with pytest.raises(KeyError):
        ads.set_rate_model("Argon")
    with pytest.raises(TypeError):
        step.top.set_streams("outlet", ["not a stream"])


# ---- a deck from scratch ----------------------------------------------------------------------


def build_vsa():
    """A small two-component VSA that touches every builder operation."""
    sim = sk.Simulation.new(title="Builder test", author="tests")
    sim.add_variable("t_feed", "1 min", category="time", description="feed and evacuation time")
    sim.add_variable("t_blowdown", category="time", formula="t_feed / 3")
    sim.add_component("CarbonDioxide")
    sim.add_component("Nitrogen")
    sim.fluid.adsorbed_heat_capacity["CO2"] = "30.7 J/(mol*K)"

    ads = sim.add_adsorbent("13X", solid_density="1138 kg/m3", heat_capacity=900.0,
                            thermal_conductivity=0.63, particle_void_fraction=0.35)
    ads.set_isotherm("extended_dual_site_langmuir_c", {
        "CO2": {"m1": 2.768, "b1": 3.934e-7, "q1": "40.98 kJ/mol", "m2": 2.758, "b2": 1.256e-7, "q2": 34190.0},
        "N2": {"m1": 2.768, "b1": 3.134e-6, "q1": 15710.0, "m2": 2.758, "b2": 3.134e-6, "q2": 15720.0},
    })
    ads.set_rate_model("CO2", "ldf", "macropore", tortuosity=3.0, macropore_void_fraction=0.35)
    ads.set_rate_model("N2", "ldf", "constant", k=0.5)
    ads.heat_of_adsorption["CO2"] = "30.85 kJ/mol"
    ads.heat_of_adsorption["N2"] = 15930.0
    inert = sim.add_adsorbent("Beads", solid_density=2500.0, heat_capacity=800.0, thermal_conductivity=1.0,
                              particle_void_fraction=0.3)
    inert.set_isotherm("inert", components=["CO2", "N2"])

    sim.bed.set_shape("cylinder", diameter="3.2 cm")
    sim.bed.set_thermal_equilibrium("lte")
    sim.wall.set_material("ss316l")
    sim.wall.thickness = "3.1 mm"
    sim.wall.set_boundary("zero_flux")

    layer = sim.add_layer("13X", composition={"13X": 0.9, "Beads": 0.1}, momentum="kozeny_carman",
                          length="0.3 m", cells=20, bed_void_fraction=0.41, particle_diameter="1.5 mm")
    layer.set_dispersion("ruthven", molecular_diffusion_coefficient=1.6e-5)
    layer.set_heat_dispersion("constant", heat_dispersion_coefficient=0.0903)
    layer.set_internal_heat_transfer("constant", internal_htc=8.6)
    layer.set_external_heat_transfer("constant", external_htc=10.0)

    sim.initial_conditions.pressure = "1 bar"
    sim.initial_conditions.temperature = 298.0
    sim.initial_conditions.wall_temperature = 298.0
    sim.initial_conditions.y["N2"] = 1.0
    sim.initial_conditions.y["CO2"] = 0.0

    feed = sim.add_step("Feed", time="t_feed")
    feed.bottom.set_streams("inlet", [
        sk.Stream("Feed inlet", flow=sk.velocity(sk.linear("8 cm/s", start="2 cm/s", time="5 s")),
                  y={"CO2": 0.15, "N2": 0.85}, temperature=sk.exponential(310.0, start=298.0, time=10.0)),
    ])
    product, purge = feed.top.set_streams("outlet", [
        sk.ValveStream("Product", pressure="1 bar", cv=10.0, xt=0.53),
        sk.Stream("Purge out", flow=sk.molar_flow(1e-4)),
    ])
    feed.add_energy_source(layer, power="-1 W")

    blowdown = sim.add_step("Blowdown", time="t_blowdown")
    blowdown.top.set_pressure("outlet", sk.exponential("0.2 bar"), y={"N2": 1.0})

    evacuation = sim.add_step("Evacuation", time="t_feed")  # connected to Feed: same duration
    evacuation.bottom.set_pressure("outlet", sk.exponential("0.05 bar"),
                                   pump=sk.Pump("vacuum_pump", "adiabatic", "vacuum_pump", eta_max=0.7, alpha=2.0,
                                                discharge_pressure="1 bar"))
    evacuation.top.set_streams("inlet", [
        sk.Stream("Purge in", flow=sk.molar_flow(1e-4), y={"N2": 1.0}, temperature=298.0),
    ])

    pressurisation = sim.add_step("Pressurisation", time=20.0)
    pressurisation.bottom.set_pressure("inlet", sk.linear("1 bar"), y={"CO2": 0.15, "N2": 0.85}, temperature=298.0)

    sim.connect(feed.top["Purge out"], evacuation.top["Purge in"])
    sim.kpi.track("CO2", feed=[feed.bottom, pressurisation.bottom], product=[evacuation.bottom])
    return sim, dict(feed=feed, blowdown=blowdown, evacuation=evacuation, product=product, purge=purge, layer=layer)


def test_from_scratch_structure():
    sim, parts = build_vsa()
    d = sim.describe()
    assert d["components"] == ["CarbonDioxide", "Nitrogen"]
    ads = {a["name"]: a for a in d["adsorbents"]}
    assert ads["13X"]["isotherms"] == [{"kind": "extended_dual_site_langmuir_c",
                                        "components": ["CarbonDioxide", "Nitrogen"]}]
    assert ads["13X"]["rate_models"]["Nitrogen"] == {"model": "ldf", "coefficient": "constant"}
    assert ads["Beads"]["isotherms"][0]["kind"] == "inert"
    assert d["bed"]["shape"] == "cylinder" and d["bed"]["wall_material"] == "ss316l"
    (layer,) = d["layers"]
    assert layer["momentum"] == "kozeny_carman" and layer["dispersion"] == "ruthven"
    assert sorted(f for _, f in layer["composition"]) == [0.1, 0.9]
    assert [s["name"] for s in d["steps"]] == ["Feed", "Blowdown", "Evacuation", "Pressurisation"]
    feed, blowdown, evac, press = d["steps"]
    assert feed["bottom"]["kind"] == "flow" and feed["bottom"]["direction"] == "inlet"
    assert [s["type"] for s in feed["top"]["streams"]] == ["valve", "stream"]
    assert feed["bottom"]["streams"][0]["flow"] == "velocity"
    assert feed["bottom"]["streams"][0]["profile"] == "linear"
    assert feed["bottom"]["streams"][0]["temperature_profile"] == "exponential"
    assert blowdown["bottom"]["kind"] == "closed"
    assert blowdown["top"]["profile"] == "exponential"
    assert evac["bottom"]["pump"] == {"kind": "vacuum_pump", "energy": "adiabatic", "efficiency": "vacuum_pump"}
    assert press["top"]["kind"] == "closed"
    assert len(d["connections"]) == 1
    assert d["kpis"] == [{"component": "CarbonDioxide",
                          "feed": [feed["bottom"]["streams"][0]["endpoint"], press["bottom"]["endpoint"]],
                          "product": [evac["bottom"]["endpoint"]], "recycle": []}]
    assert connections(sim) == {frozenset({"Feed/top/1", "Evacuation/top/0"})}

    assert d["variables"] == [
        {"name": "t_feed", "category": "time", "formula": "", "description": "feed and evacuation time"},
        {"name": "t_blowdown", "category": "time", "formula": "t_feed / 3", "description": ""},
    ]
    assert sim.variables["t_feed"] == 60.0
    assert sim.steps["Feed"].time == 60.0 and sim.get_text("step[Feed].time") == "t_feed"
    assert sim.steps["Blowdown"].time == pytest.approx(20.0)

    # values went through the registry, units converted
    assert sim.get("bed.diameter") == pytest.approx(0.032)
    assert sim.layers["13X"].particle_diameter == pytest.approx(1.5e-3)
    assert sim.layers["13X"].cells == 20
    assert sim.adsorbents["13X"].isotherm["CO2"].q1 == pytest.approx(40980.0)
    assert sim.adsorbents["13X"].kinetics["N2"].k == 0.5
    assert sim.fluid.adsorbed_heat_capacity["CO2"] == 30.7
    inlet = parts["feed"].bottom["Feed inlet"]
    assert inlet.velocity_linear_start == pytest.approx(0.02)
    assert inlet.velocity_linear_end == pytest.approx(0.08)
    assert inlet.temperature_exponential_end == 310.0
    assert parts["product"].valve.cv == 10.0
    assert parts["purge"].molar_flow == 1e-4
    pump = parts["evacuation"].bottom.pump
    assert (pump.eta_max, pump.alpha, pump.discharge_pressure) == (0.7, 2.0, 1e5)
    assert parts["feed"].energy_sources["13X"].power == -1.0
    assert sim.check()[0], sim.check()


def test_from_scratch_runs():
    sim, _ = build_vsa()
    r = sim.run(cycles=2)
    assert r.ok, r.message
    (kpi,) = sim.kpis(cycle=-1)
    assert kpi.component == "CarbonDioxide" and kpi.valid and 0 < kpi.recovery <= 1


def test_recycle_counts_against_recovery_not_purity():
    # N2 product leaves at the Feed top; the purge stream carries part of it back
    # into the bed during Evacuation. As recycle it lowers the net recovery and
    # leaves the purity alone.
    sim, parts = build_vsa()
    feed, evacuation = parts["feed"], parts["evacuation"]
    sim.kpi.track("N2", feed=[feed.bottom], product=[feed.top])
    assert sim.run(cycles=2).ok
    gross = {k.component: k for k in sim.kpis(cycle=-1)}["Nitrogen"]
    sim.kpi.track("N2", feed=[feed.bottom], product=[feed.top], recycle=[evacuation.top["Purge in"]])
    assert sim.describe()["kpis"][1]["recycle"] == [evacuation.top["Purge in"].id]
    net = {k.component: k for k in sim.kpis(cycle=-1)}["Nitrogen"]
    assert net.purity == pytest.approx(gross.purity, rel=1e-12)
    assert 0 < net.recovery < gross.recovery
    assert net.production_rate < gross.production_rate
    script = sk.to_script(sim)
    assert "recycle=[" in script


def test_from_scratch_round_trips_through_a_script():
    sim, _ = build_vsa()
    ns = {}
    exec(compile(sk.to_script(sim), "<vsa>", "exec"), ns)
    assert_same_deck(sim, ns["sim"])


def test_blake_kozeny_momentum_builds_and_round_trips():
    sim, parts = build_vsa()
    layer = parts["layer"].set_momentum("blake_kozeny")
    assert sim.describe()["layers"][0]["momentum"] == "blake_kozeny"
    assert layer.momentum_model == "blake_kozeny"
    assert sim.get("layer[13X].momentum_model") == 2
    assert '"Momentum_Model": 2' in sim.to_json()
    copy = sk.Simulation.from_json(sim.to_json())
    assert copy.describe()["layers"][0]["momentum"] == "blake_kozeny"
    assert_same_deck(sim, copy)
    script = sk.to_script(sim)
    assert "momentum='blake_kozeny'" in script
    ns = {}
    exec(compile(script, "<blake_kozeny>", "exec"), ns)
    assert_same_deck(sim, ns["sim"])
    # By registry name too, and back.
    sim["layer[13X].momentum_model"] = "ergun"
    assert sim.describe()["layers"][0]["momentum"] == "ergun"
    sim["layer[13X].momentum_model"] = "blake_kozeny"
    assert sim.describe()["layers"][0]["momentum"] == "blake_kozeny"
    assert sim.run(cycles=1).ok


def test_constant_macropore_diffusivity_builds_and_round_trips():
    sim, _ = build_vsa()
    ads = sim.adsorbents["13X"]
    for kind in ("macropore", "macropore_capacity", "macropore_knudsen"):
        ads = ads.set_rate_model("CO2", "ldf", kind, tortuosity=3.0, macropore_void_fraction=0.35)
        assert ads.kinetics["CO2"].molecular_diffusion_coefficient == 0.0  # default: from T and P
        ads = ads.set_rate_model("CO2", "ldf", kind, tortuosity=3.0, macropore_void_fraction=0.35,
                                 molecular_diffusion_coefficient=1.6e-5)
        kn = ads.kinetics["CO2"]
        assert kn.molecular_diffusion_coefficient == pytest.approx(1.6e-5)
        path = kn.path("molecular_diffusion_coefficient")
        assert path.endswith("kinetics[CarbonDioxide].molecular_diffusion_coefficient")
        sim[path] = "0.16 [cm**2/s]"
        assert sim.get(path) == pytest.approx(1.6e-5)
        with pytest.raises(ValueError):
            kn.molecular_diffusion_coefficient = 2 * sk.ureg.bar
        assert '"Molecular_Diffusion_Coefficient"' in sim.to_json()
        copy = sk.Simulation.from_json(sim.to_json())
        assert copy.adsorbents["13X"].kinetics["CO2"].molecular_diffusion_coefficient == pytest.approx(1.6e-5)
        assert_same_deck(sim, copy)
        script = sk.to_script(sim)
        assert "molecular_diffusion_coefficient=" in script
        ns = {}
        exec(compile(script, f"<{kind}>", "exec"), ns)
        assert_same_deck(sim, ns["sim"])
    assert sim.run(cycles=1).ok


def test_log_linear_profiles_build_and_round_trip():
    sim, parts = build_vsa()
    press = sim.steps["Pressurisation"]
    press.bottom.set_pressure("inlet", sk.log_linear("1 bar"), y={"CO2": 0.15, "N2": 0.85}, temperature=298.0)
    parts["feed"].bottom.set_streams("inlet", [
        sk.Stream("Feed inlet", flow=sk.velocity(sk.log_linear("8 cm/s", start="2 cm/s", time="5 s")),
                  y={"CO2": 0.15, "N2": 0.85}, temperature=sk.log_linear(310.0, start=298.0, time=10.0)),
    ])
    d = {s["name"]: s for s in sim.describe()["steps"]}
    assert d["Pressurisation"]["bottom"]["profile"] == "log_linear"
    assert d["Feed"]["bottom"]["streams"][0]["profile"] == "log_linear"
    assert d["Feed"]["bottom"]["streams"][0]["temperature_profile"] == "log_linear"
    assert press.bottom.pressure_log_linear_end == pytest.approx(1e5)
    inlet = parts["feed"].bottom["Feed inlet"]
    assert inlet.velocity_log_linear_start == pytest.approx(0.02)
    assert inlet.velocity_log_linear_time == pytest.approx(5.0)
    assert inlet.temperature_log_linear_end == 310.0
    script = sk.to_script(sim)
    assert "sk.log_linear(" in script
    ns = {}
    exec(compile(script, "<log_linear>", "exec"), ns)
    assert_same_deck(sim, ns["sim"])
    assert sim.run(cycles=1).ok


def test_structure_edits_replace_what_was_there():
    sim, parts = build_vsa()
    feed = parts["feed"]
    feed.top.set_closed()
    d = sim.describe()
    assert d["steps"][0]["top"]["kind"] == "closed"
    assert d["connections"] == []  # the connected stream went with the end
    ads = sim.adsorbents["13X"]
    ads.set_isotherm("henry", {"N2": {"k": 1e-4, "Q": 15000.0}})
    kinds = {i["kind"]: i["components"] for i in sim.describe()["adsorbents"][0]["isotherms"]}
    assert kinds == {"extended_dual_site_langmuir_c": ["CarbonDioxide"], "henry": ["Nitrogen"]}
    assert sim.adsorbents["13X"].isotherm["N2"].k == 1e-4
    parts["layer"].set_momentum("ergun")
    assert sim.describe()["layers"][0]["momentum"] == "ergun"
    parts["evacuation"].bottom.set_pump(None)
    assert sim.describe()["steps"][2]["bottom"]["pump"] is None


DOCS = Path(__file__).resolve().parents[2] / "docs" / "python.md"


@pytest.mark.skipif(not DOCS.is_file(), reason="docs/ not found (not a repository checkout)")
def test_the_documented_example_runs():
    section = DOCS.read_text(encoding="utf-8").split("## Building a simulation in Python", 1)[1]
    code = section.split("```python\n", 1)[1].split("```", 1)[0]
    ns = {}
    exec(code.replace("number_of_cycles = 50", "number_of_cycles = 1"), ns)
    assert ns["result"].ok, ns["result"].message
    (kpi,) = ns["sim"].kpis(cycle=-1)
    assert kpi.valid and kpi.energy_consumption > 0  # the vacuum pump counts


def test_objects_held_in_variables_follow_structural_edits():
    sim, parts = build_vsa()
    bed, feed, layer = sim.bed, parts["feed"], parts["layer"]
    sim.bed.set_shape("square", width="5 cm")
    assert bed.width == pytest.approx(0.05) and "diameter" not in bed.fields
    blowdown_top = parts["blowdown"].top
    blowdown_top.set_streams("outlet", [sk.Stream("Vent", flow=sk.velocity(0.01))])
    assert blowdown_top["Vent"].velocity == 0.01
    layer.set_dispersion("constant", dispersion_coefficient=1e-4)
    assert layer.dispersion_coefficient == 1e-4
    assert feed.energy_sources["13X"].power == -1.0


# ---- the example decks, rebuilt -------------------------------------------------------------


@pytest.mark.parametrize("example", sk.examples.ALL, ids=lambda e: e.name)
def test_every_example_round_trips_through_to_script(example):
    deck = example.load()
    ns = {}
    exec(compile(sk.to_script(deck), example.name, "exec"), ns)
    assert_same_deck(deck, ns["sim"])


@pytest.mark.parametrize("name", RUN_CHECKED)
def test_generated_scripts_run_like_their_decks(name):
    deck = sk.examples.get(name).load()
    ns = {}
    exec(compile(sk.to_script(deck), name, "exec"), ns)
    assert_same_run(deck, ns["sim"], cycles=min(2, deck.cycle.number_of_cycles))


# ---- pump-driven streams, buffer tanks, experimental data ----------------------------


def test_pump_stream_sets_the_flow_from_the_power():
    from skarstrom import ureg as u

    sim = sk.examples.Perez2019_A.load()
    feed, evac = sim.steps["Feed"], sim.steps["Evacuation"]
    (vac,) = evac.bottom.set_streams("outlet", [
        sk.PumpStream("Vacuum pump", power=40 * u.W,
                      pump=sk.Pump(efficiency_model="vacuum_pump", eta_max=0.7, alpha=19.55,
                                   discharge_pressure=1 * u.bar)),
    ])
    assert vac.power == 40.0 and vac.pump.eta_max == 0.7 and vac.pump.discharge_pressure == 1e5
    st = sim.describe()["steps"][2]["bottom"]["streams"][0]
    assert st["type"] == "pump" and st["power"] == "constant"
    assert st["pump"] == {"kind": "vacuum_pump", "energy": "isothermal", "efficiency": "vacuum_pump"}
    sim.kpi.track("CO2", feed=[feed.bottom], product=[vac])
    r = sim.run(cycles=2)
    assert r.ok, r.message
    out = sim.endpoint_series(vac.id, "molar_flow_rate")
    assert np.all(out.values < 0) and out.values.min() < -1e-5  # gas leaves through the bottom
    (kpi,) = sim.kpis(cycle=-1)
    assert kpi.valid and kpi.energy_consumption > 0

    rebuilt = {}
    exec(compile(sk.to_script(sim), "<pump>", "exec"), rebuilt)
    assert_same_deck(sim, rebuilt["sim"])


def test_pump_stream_with_a_power_polynomial():
    sim = sk.examples.Perez2019_A.load()
    (vac,) = sim.steps["Evacuation"].bottom.set_streams("outlet", [
        sk.PumpStream("Vacuum pump", power=sk.polynomial_power(a=20.0, b=1e-3), pump=sk.Pump(efficiency=0.6)),
    ])
    assert (vac.power_poly_a, vac.power_poly_b, vac.power_poly_f) == (20.0, 1e-3, 0.0)
    assert vac.pump.efficiency == 0.6
    assert sim.describe()["steps"][2]["bottom"]["streams"][0]["power"] == "polynomial"
    sim.kpi.track("CO2", feed=[sim.steps["Feed"].bottom], product=[vac])
    assert sim.check()[0], sim.check()
    with pytest.raises(ValueError, match="always has a pump"):
        vac.set_pump(None)
    rebuilt = {}
    exec(compile(sk.to_script(sim), "<poly>", "exec"), rebuilt)
    assert_same_deck(sim, rebuilt["sim"])


def test_buffer_tank_feeds_an_inlet_from_another_steps_product():
    from skarstrom import ureg as u

    sim = sk.examples.Perez2019_A.load()
    product = sim.steps["Feed"].top                       # the N2 product outlet
    press = sim.steps["Pressurisation"].bottom
    assert press.buffer_tank is None
    press.set_buffer_tank(product, percentage=40 * u.percent)
    tank = sim.steps["Pressurisation"].bottom.buffer_tank
    assert tank.percentage == pytest.approx(40.0) and tank.source.id == product.id
    assert sim.describe()["steps"][3]["bottom"]["buffer_tank"] == product.id
    with pytest.raises(ValueError, match="only an inlet"):
        sim.steps["Blowdown"].top.set_buffer_tank(product)
    assert sim.run(cycles=2).ok

    rebuilt = {}
    script = sk.to_script(sim)
    assert ".set_buffer_tank(" in script
    exec(compile(script, "<tank>", "exec"), rebuilt)
    assert_same_deck(sim, rebuilt["sim"])

    sim.steps["Pressurisation"].bottom.set_buffer_tank(None)
    assert sim.steps["Pressurisation"].bottom.buffer_tank is None


def test_experimental_data_is_stored_in_the_deck(tmp_path):
    sim = sk.examples.Ward2022_Q100ccm.load()
    sim.experimental_data.clear()
    assert len(sim.experimental_data) == 0
    csv = tmp_path / "outlet.csv"
    csv.write_text("time [s],y_CO2 [-]\n0,0\n1000,0.2\n2000,0.9\n3000,1.0\n")
    s = sim.experimental_data.load_csv(csv)
    assert s.name == "outlet" and s.x_label == "time [s]" and list(s.y) == [0, 0.2, 0.9, 1.0]
    sim.experimental_data.add("T at 7 cm", [0, 60], [291.4, 300.0], "time [s]", "T [K]")
    assert sim.experimental_data.names() == ["outlet", "T at 7 cm"]

    copy = sk.Simulation.from_json(sim.to_json())          # saved with the deck
    assert copy.experimental_data["T at 7 cm"].y.tolist() == [291.4, 300.0]

    assert sim.run().ok
    y = sim.series("mole_fraction", "CarbonDioxide", location="top", cycle=0)
    fit = sim.experimental_data["outlet"].compare(y.time, y.values)
    assert fit["points"] == 4 and 0 <= fit["nrmse"] < 1
    sim.experimental_data.remove("outlet")
    assert "outlet" not in sim.experimental_data and len(sim.experimental_data) == 1


def test_double_exponential_pressure_and_capacity_factor_build_and_round_trip():
    """The Perez 2019 forms: a two-term exponential evacuation pressure and the
    capacity-factor macropore coefficient (SI Eqs. S4 and the evacuation profile)."""
    sim = sk.examples.Perez2019_A.load()
    evac = sim.steps["Evacuation"]
    evac.bottom.set_pressure("outlet", sk.double_exponential("3.2 kPa", rate_1=0.019, rate_2=0.101, weight=0.3),
                             y={"CarbonDioxide": 0.0, "Nitrogen": 1.0}, temperature=298.0)
    assert evac.bottom.pressure_double_exponential_end == pytest.approx(3200.0)
    assert evac.bottom.pressure_double_exponential_rate_1 == pytest.approx(0.019)
    assert evac.bottom.pressure_double_exponential_weight == pytest.approx(0.3)
    ads = sim.adsorbents[0]
    ads.set_rate_model("CarbonDioxide", "ldf", "macropore_capacity", tortuosity=3.0, macropore_void_fraction=0.35)
    assert ads.kinetics["CarbonDioxide"].concentration_floor == pytest.approx(1e-3)
    assert "macropore_capacity" in sk.build_kinds()["rate_coefficient"]
    d = {s["name"]: s for s in sim.describe()["steps"]}
    assert d["Evacuation"]["bottom"]["profile"] == "double_exponential"
    script = sk.to_script(sim)
    assert "sk.double_exponential(" in script and "macropore_capacity" in script
    ns = {}
    exec(compile(script, "<double_exponential>", "exec"), ns)
    assert_same_deck(sim, ns["sim"])
    result = sim.run(cycles=1)
    assert result.ok, result.message
    k = sim.series("mass_transfer_coefficient", "CarbonDioxide", location="bottom")
    import numpy as np
    assert np.all(np.isfinite(k.values)) and np.all(k.values >= 0.0)
