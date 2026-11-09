"""The check valve of a valve stream and the outlet backflow warning."""

import skarstrom as sk

VENT = "endpoint[2].valve.check_valve"  # Perez2019_D1: the Feed 1 vent


def _backflow(sim):
    return [w for w in sim.operating_limits().warnings if w.limit == "outlet_backflow"]


def test_registry_model_and_defaults():
    sim = sk.examples.Perez2019_D1.load()
    paths = {p.path: p for p in sim.parameters(search="check_valve")}
    assert set(paths) == {"endpoint[2].valve.check_valve", "endpoint[5].valve.check_valve"}
    assert paths[VENT].kind == "boolean" and paths[VENT].value is False
    assert sim.get("settings.warn_outlet_backflow") is True
    vent = sim.steps["Feed 1"].top["Feed 1 ambient Stream"].valve
    assert vent.check_valve is False
    vent.check_valve = True
    assert sim.get(VENT) is True
    sim.settings.warn_outlet_backflow = False
    assert sim.get("settings.warn_outlet_backflow") is False


def test_builder_and_script():
    sim = sk.Simulation.new(title="check valve")
    sim.add_component("CarbonDioxide")
    sim.add_component("Nitrogen")
    step = sim.add_step("Feed", time=10)
    step.bottom.set_streams("inlet", [sk.Stream("Feed in", flow=sk.velocity(0.01),
                                                y={"CO2": 0.15, "N2": 0.85}, temperature=298)])
    plain, checked = step.top.set_streams("outlet", [
        sk.ValveStream("Vent", pressure=1e5, cv=1.0),
        sk.ValveStream("Checked vent", pressure=1e5, cv=1.0, check_valve=True),
    ])
    assert plain.valve.check_valve is False
    assert checked.valve.check_valve is True
    script = sk.to_script(sim)
    assert script.count("check_valve=True") == 1


def test_d1_vent_backflow_is_reported_and_stopped_by_the_check_valve():
    sim = sk.examples.Perez2019_D1.load()
    r = sim.run(cycles=3)
    assert r.status in ("completed", "reached_css")
    (w,) = _backflow(sim)
    assert w.step == 0 and w.value > w.limit_value
    assert "Feed 1 ambient Stream" in w.message and VENT in w.message
    assert r.healthy  # an operating-limit warning, not a run-health one

    sim = sk.examples.Perez2019_D1.load()
    sim[VENT] = True
    sim["endpoint[5].valve.check_valve"] = True
    r = sim.run(cycles=3)
    assert r.status in ("completed", "reached_css")
    assert not _backflow(sim)
