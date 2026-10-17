"""Tests for building a simulation from scratch (skarstrom.build and the builder methods)."""

import runpy
from pathlib import Path

import numpy as np
import pytest

import skarstrom as sk
from deck_compare import connections, diff, normalise, parameter_diff

SEEDS = Path(__file__).resolve().parents[2] / "examples" / "_build"
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
    assert "ergun" in kinds["momentum"] and "kozeny_carman" in kinds["momentum"]
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
                          "product": [evac["bottom"]["endpoint"]]}]
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


def test_from_scratch_round_trips_through_a_script():
    sim, _ = build_vsa()
    ns = {}
    exec(compile(sk.to_script(sim), "<vsa>", "exec"), ns)
    assert_same_deck(sim, ns["sim"])


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


@pytest.mark.skipif(not SEEDS.is_dir(), reason="examples/_build not found (not a repository checkout)")
@pytest.mark.parametrize("example", sk.examples.ALL, ids=lambda e: e.name)
def test_seed_scripts_reproduce_their_decks(example):
    seed = SEEDS / f"{example.name}.py"
    assert seed.is_file(), "run scripts/gen_build_seeds.py"
    assert_same_deck(example.load(), runpy.run_path(str(seed))["sim"])


@pytest.mark.skipif(not SEEDS.is_dir(), reason="examples/_build not found (not a repository checkout)")
@pytest.mark.parametrize("name", RUN_CHECKED)
def test_seed_scripts_run_like_their_decks(name):
    deck = sk.examples.get(name).load()
    rebuilt = runpy.run_path(str(SEEDS / f"{name}.py"))["sim"]
    assert_same_run(deck, rebuilt, cycles=min(2, deck.cycle.number_of_cycles))
