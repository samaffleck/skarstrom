"""The reaction builder: ``Adsorbent.add_reaction`` and the ``Reaction`` setters, ``describe()`` and
``to_script`` for reactive decks, and agreement with the low-level ``build`` operations that wrap
``adsim::build`` (the path the C++ reaction tests use)."""

import pytest

import skarstrom as sk

from deck_compare import diff, normalise

CH4, H2O, CO, H2, CO2 = "Methane", "Water", "CarbonMonoxide", "Hydrogen", "CarbonDioxide"
SMR = {CH4: -1, H2O: -1, CO: 1, H2: 3}
WGS = {CO: -1, H2O: -1, CO2: 1, H2: 1}


def _base():
    """Five components, a catalyst and a CO2 sorbent mixed in one layer, one feed step."""
    sim = sk.Simulation.new(title="reformer")
    for c in (CH4, H2O, CO, H2, CO2):
        sim.add_component(c)
    ni = sim.add_adsorbent("Ni catalyst", solid_density=2000, heat_capacity=900, thermal_conductivity=0.5,
                           particle_void_fraction=0.4)
    sorbent = sim.add_adsorbent("Sorbent", solid_density=1500, heat_capacity=900, thermal_conductivity=0.5,
                                particle_void_fraction=0.4)
    # a real isotherm on every component of both solids (the sorbent takes up CO2)
    for ads in (ni, sorbent):
        ads.set_isotherm("henry", {c: {"k": 1e-5 if (ads is sorbent and c == CO2) else 0.0}
                                   for c in (CH4, H2O, CO, H2, CO2)})
    sim.bed.set_shape("cylinder", diameter=0.05)
    sim.wall.set_material("ss316l")
    sim.wall.set_boundary("zero_flux")
    sim.add_layer("Layer 1", composition={"Ni catalyst": 0.3, "Sorbent": 0.7}, momentum="kozeny_carman",
                  length=0.2, cells=10, bed_void_fraction=0.4, particle_diameter=3e-3)
    sim.initial_conditions.pressure = 1e5
    sim.initial_conditions.temperature = 800
    sim.initial_conditions.y[H2O] = 1.0
    feed = sim.add_step("Feed", time=2.0)
    feed.bottom.set_streams("inlet", [sk.Stream("Feed inlet", flow=sk.velocity(0.05),
                                                y={CH4: 0.2, H2O: 0.8}, temperature=800)])
    feed.top.set_pressure("outlet", 1e5)
    return sim, ni


def build_reactor():
    sim, ni = _base()
    smr = ni.add_reaction("SMR", SMR, kind="lhhw", adsorption_terms=2, heat_of_reaction=206100.0)
    smr.rate.pre_exp_factor = 1e3
    smr.rate.activation_energy = 8e4
    smr.rate.order[H2] = 0.5
    smr.rate.adsorption[2].exponent[H2] = -1.0
    wgs = ni.add_reaction("WGS", WGS, kind="reversible_power_law", basis="concentration",
                          effectiveness="thiele", key_component=CO)
    wgs.rate.pre_exp_factor = 5.0
    sim.kpi.track_reaction(CH4, products=[H2])
    return sim


def build_reactor_low_level():
    """The same deck through the raw build operations."""
    sim, ni = _base()
    a = ni.id
    smr = sim._build("add_reaction", a, "SMR", SMR, "lhhw")
    sim._build("add_adsorption_term", a, smr)
    sim._build("add_adsorption_term", a, smr)
    wgs = sim._build("add_reaction", a, "WGS", WGS, "reversible_power_law")
    sim._build("set_reaction_basis", a, wgs, "concentration")
    sim._build("set_effectiveness", a, wgs, "thiele", CO)
    sim["adsorbent[%d].reaction[%d].heat_of_reaction" % (a, smr)] = 206100.0
    sim["adsorbent[%d].reaction[%d].rate.pre_exp_factor" % (a, smr)] = 1e3
    sim["adsorbent[%d].reaction[%d].rate.activation_energy" % (a, smr)] = 8e4
    sim["adsorbent[%d].reaction[%d].rate.order[%s]" % (a, smr, H2)] = 0.5
    sim["adsorbent[%d].reaction[%d].rate.adsorption[2].exponent[%s]" % (a, smr, H2)] = -1.0
    sim["adsorbent[%d].reaction[%d].rate.pre_exp_factor" % (a, wgs)] = 5.0
    sim.kpi.track_reaction(CH4, products=[H2])
    return sim


def assert_same_deck(a, b):
    d = diff(normalise(a.to_json()), normalise(b.to_json()))
    assert not d, d[:5]


def test_add_reaction_creates_structure_and_fields():
    sim, ni = _base()
    assert len(ni.reactions) == 0
    smr = ni.add_reaction("SMR", SMR, kind="lhhw", adsorption_terms=2)
    assert (smr.id, smr.name) == (1, "SMR")
    assert smr.nu[H2] == 3.0
    assert smr.basis == "partial_pressure"
    assert "adsorption" in smr.rate.keys() if hasattr(smr.rate, "keys") else True
    assert sim.describe()["adsorbents"][0]["reactions"] == [{
        "id": 1, "name": "SMR", "nu": SMR, "kind": "lhhw", "basis": "partial_pressure",
        "effectiveness": "constant", "key_component": "", "adsorption_terms": 2}]
    assert ni.reactions["SMR"].id == 1
    assert sim.adsorbents["Sorbent"].reactions.names == []


def test_setters_and_errors():
    sim, ni = _base()
    r = ni.add_reaction("R", {CH4: -1, H2: 2})
    assert r.set_kind("lhhw").add_adsorption_term() == 1
    assert sim.describe()["adsorbents"][0]["reactions"][0]["adsorption_terms"] == 1
    r.remove_adsorption_term(1)
    assert sim.describe()["adsorbents"][0]["reactions"][0]["adsorption_terms"] == 0
    r.set_stoichiometry({CH4: -1, CO: 1, H2: 3}).set_basis("concentration")
    d = sim.describe()["adsorbents"][0]["reactions"][0]
    assert d["nu"] == {CH4: -1, CO: 1, H2: 3} and d["basis"] == "concentration"
    with pytest.raises(ValueError, match="key reactant"):
        r.set_effectiveness("thiele", H2)
    r.set_effectiveness("thiele", CH4)
    assert sim.describe()["adsorbents"][0]["reactions"][0]["key_component"] == CH4
    r.set_effectiveness("constant")
    with pytest.raises(ValueError):
        ni.add_reaction("bad", {CH4: -1}, kind="michaelis")
    with pytest.raises(ValueError):
        r.set_kind("power_law").add_adsorption_term()
    ni.remove_reaction("R")
    assert len(ni.reactions) == 0


def test_python_build_equals_low_level_build():
    assert_same_deck(build_reactor(), build_reactor_low_level())


def test_reactive_deck_round_trips_through_a_script():
    sim = build_reactor()
    script = sk.to_script(sim)
    assert "add_reaction('SMR'" in script and "adsorption_terms=2" in script
    assert "effectiveness='thiele', key_component='CarbonMonoxide'" in script
    ns = {}
    exec(compile(script, "<reactor>", "exec"), ns)
    assert_same_deck(sim, ns["sim"])
    # the script of the rebuilt deck is the same script
    assert sk.to_script(ns["sim"]) == script
    assert_same_deck(sim, sk.Simulation.from_json(sim.to_json()))


def test_builds_run_to_the_same_kpis():
    a, b = build_reactor(), build_reactor_low_level()
    assert a.run(cycles=1).ok and b.run(cycles=1).ok
    ka = {k.component: k for k in a.kpis(cycle=-1)}
    kb = {k.component: k for k in b.kpis(cycle=-1)}
    assert ka.keys() == kb.keys()
    for c in ka:
        for f in ("conversion", "reaction_generation", "peak_temperature"):
            x, y = getattr(ka[c], f), getattr(kb[c], f)
            assert x == pytest.approx(y, rel=1e-12, nan_ok=True)
    assert ka[CH4].conversion > 0


def test_build_kinds_list_reaction_kinds():
    k = sk.build_kinds()
    assert k["reaction_rate"] == ["power_law", "reversible_power_law", "lhhw", "expression"]
    assert k["reaction_basis"] == ["partial_pressure", "concentration"]
    assert k["effectiveness"] == ["constant", "thiele"]
