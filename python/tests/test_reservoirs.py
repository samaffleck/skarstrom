"""Named sources and sinks (docs/sources-and-sinks.md): the object layer, the
string paths, the structural edits and the script writer."""

import pytest

import skarstrom as sk


@pytest.fixture
def sim():
    # Feed: a flow inlet at the bottom, a pressure outlet at the top.
    # Pressurisation: a pressure inlet at the bottom.
    return sk.examples.Perez2019_A.load()


def test_add_attach_and_shared_values(sim):
    assert len(sim.reservoirs) == 0
    feed = sim.add_source("Feed gas", temperature=300.0, pressure=1.2e5, y={"CO2": 0.2, "N2": 0.8})
    ambient = sim.add_sink("Ambient", pressure=1.0e5)
    assert [r.name for r in sim.reservoirs] == ["Feed gas", "Ambient"]
    assert (feed.kind, ambient.kind) == ("source", "sink")
    assert feed.temperature == 300.0 and feed.pressure == 1.2e5
    assert feed.y["CO2"] == pytest.approx(0.2)

    inlet = sim.steps["Feed"].bottom.attach(feed)
    press = sim.steps["Pressurisation"].bottom.attach(feed)
    sim.steps["Feed"].top.attach("Ambient")
    assert inlet.reservoir.name == "Feed gas"
    assert sorted(e.step for e in feed.endpoints) == ["Feed", "Pressurisation"]

    # The values are the reservoir's: listed once, and one set reaches both inlets.
    paths = {p.path for p in sim.parameters()}
    rp = f"reservoir[{feed.id}]"
    assert {f"{rp}.temperature", f"{rp}.pressure"} <= paths
    assert f"endpoint[{inlet.id}].temperature" not in paths
    assert f"endpoint[{press.id}].temperature" not in paths
    assert sim.get("reservoir[Feed gas].temperature") == 300.0   # the alias uses the name
    # A sink is a pressure only.
    assert f"reservoir[{ambient.id}].pressure" in paths
    assert f"reservoir[{ambient.id}].temperature" not in paths

    feed.temperature = 310.0
    press.detach()
    assert sim.get(f"endpoint[{press.id}].temperature") == 310.0   # kept as its own value
    assert sim.steps["Pressurisation"].bottom.reservoir is None


def test_kinds_must_agree(sim):
    feed = sim.add_source("Feed gas")
    ambient = sim.add_sink("Ambient")
    with pytest.raises(ValueError, match="source"):
        sim.steps["Feed"].bottom.attach(ambient)      # an inlet on a sink
    with pytest.raises(ValueError, match="sink"):
        sim.steps["Feed"].top.attach(feed)            # an outlet on a source
    with pytest.raises(KeyError):
        sim.steps["Feed"].bottom.attach("No such source")


def test_edits(sim):
    inlet = sim.steps["Feed"].bottom
    t = inlet.temperature
    rid = sim.edit.reservoir_from_endpoint(inlet, "Feed gas")
    r = sim.reservoirs.by_id(rid)
    assert (r.name, r.kind, r.temperature) == ("Feed gas", "source", t)
    assert sim.steps["Feed"].bottom.reservoir.id == rid

    sim.edit.rename_reservoir(rid, "Header")
    assert sim.reservoirs["Header"].id == rid
    other = sim.edit.add_reservoir("source")
    with pytest.raises(ValueError):
        sim.edit.rename_reservoir(other, "Header")    # names are unique
    with pytest.raises(ValueError):
        sim.edit.add_reservoir("tank")

    sim.edit.remove_reservoir("Header")
    assert sim.steps["Feed"].bottom.reservoir is None
    assert sim.steps["Feed"].bottom.temperature == t  # the endpoint keeps the values
    assert "reservoir" in sk.build_kinds()


def test_runs_like_the_same_values_on_the_endpoints(sim):
    plain = sk.examples.Perez2019_A.load()
    for step in sim.steps:
        for end in (step.bottom, step.top):
            for e in end:
                sim.edit.reservoir_from_endpoint(e)
    a = plain.run(cycles=2)
    b = sim.run(cycles=2)
    assert a.ok and b.ok
    for ka, kb in zip(plain.kpis(), sim.kpis()):
        assert ka.purity == kb.purity and ka.recovery == kb.recovery


def test_script_round_trip(sim):
    feed = sim.add_source("Feed gas", temperature=300.0, pressure=1.2e5, y={"CO2": 0.2, "N2": 0.8})
    sim.steps["Feed"].bottom.attach(feed)
    sim.steps["Feed"].top.attach(sim.add_sink("Ambient", pressure=1.0e5))
    scope = {}
    exec(sk.to_script(sim), scope)   # the script builds `sim`
    built = scope["sim"]
    assert [(r.name, r.kind) for r in built.reservoirs] == [("Feed gas", "source"), ("Ambient", "sink")]
    assert built.reservoirs["Feed gas"].temperature == 300.0
    assert built.steps["Feed"].bottom.reservoir.name == "Feed gas"
    assert built.steps["Feed"].top.reservoir.name == "Ambient"
