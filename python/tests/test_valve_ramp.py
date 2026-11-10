"""The opening ramp of a valve stream (endpoint[ID].valve.opening_start / .opening_ramp_time)."""

import skarstrom as sk


def _deck():
    sim = sk.Simulation.new(title="valve ramp")
    sim.add_component("CarbonDioxide")
    sim.add_component("Nitrogen")
    step = sim.add_step("Pressurisation", time=30)
    step.top.set_closed()
    return sim, step


def test_registry_model_and_defaults():
    sim = sk.examples.Perez2019_D1.load()   # valves without a ramp: fully open
    starts = [p for p in sim.parameters(search="opening_start")]
    ramps = [p for p in sim.parameters(search="opening_ramp_time")]
    assert starts and len(starts) == len(ramps)
    assert all(p.value == 1.0 and p.unit == "-" for p in starts)
    assert all(p.value == 0.0 and p.unit == "s" for p in ramps)
    valve = sim.steps["Feed 1"].top["Feed 1 ambient Stream"].valve
    assert valve.opening_start == 1.0 and valve.opening_ramp_time == 0.0
    valve.opening_start = 0.25
    valve.opening_ramp_time = 12
    head = starts[0].path.rsplit(".", 1)[0]
    assert sim.get(head + ".opening_start") == 0.25
    assert sim.get(head + ".opening_ramp_time") == 12.0


def test_builder_and_script():
    sim, step = _deck()
    plain, ramped = step.bottom.set_streams("inlet", [
        sk.ValveStream("Header", pressure=1e5, cv=0.01, y={"CO2": 0.15, "N2": 0.85}, temperature=298),
        sk.ValveStream("Ramped", pressure=1e5, cv=0.01, y={"CO2": 0.15, "N2": 0.85}, temperature=298,
                       opening_start=0.1, opening_ramp_time=20),
    ])
    assert plain.valve.opening_start == 1.0 and plain.valve.opening_ramp_time == 0.0
    assert ramped.valve.opening_start == 0.1 and ramped.valve.opening_ramp_time == 20.0
    script = sk.to_script(sim)
    # Only the ramped valve writes the keywords; the defaults are left out.
    assert script.count("opening_start=") == 1
    assert "opening_start=0.1, opening_ramp_time=20.0" in script


def test_expression_ramp_round_trips_through_the_script():
    sim, step = _deck()
    sim.add_variable("ramp_start", 0.2)
    sim.add_variable("ramp_time", 8, category="time")
    (valve,) = step.bottom.set_streams("inlet", [
        sk.ValveStream("Header", pressure=1e5, cv=0.01, y={"CO2": 0.15, "N2": 0.85}, temperature=298,
                       opening_start="ramp_start", opening_ramp_time="ramp_time"),
    ])
    assert valve.valve.opening_start == 0.2
    script = sk.to_script(sim)
    assert "opening_start='ramp_start', opening_ramp_time='ramp_time'" in script


def test_case_study_decks_ramp_their_pressurisation_valves():
    for example, n in ((sk.examples.Perez2019_CaseStudy_Base, 1), (sk.examples.Perez2019_CaseStudy_Base_LPP, 2)):
        sim = example.load()
        ramped = [p for p in sim.parameters(search="opening_start") if sim.get_text(p.path) == "valve_ramp_start"]
        assert len(ramped) == n, example.name
        assert 0 < sim.get("variable[valve_ramp_start]") < 1
        assert sim.get("variable[valve_ramp_time]") > 0
