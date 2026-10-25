"""Wall pressure design (sim.wall_design) and the operating-limit report of a
run (sim.operating_limits)."""

import math

import pytest

import skarstrom

PEREZ = skarstrom.examples.Perez2019_A  # four-step CO2 VSA (evacuates below 1 atm)


def test_wall_design_ug27_and_vacuum():
    sim = skarstrom.Simulation(PEREZ)
    sim["wall.design_pressure"] = 6e5
    sim["wall.design_temperature"] = 293.15
    sim["wall.corrosion_allowance"] = 0.001
    d = sim.wall_design()
    assert d.valid, d.message
    S, E, c = d.allowable_stress, sim.get("wall.weld_joint_efficiency"), 0.001
    R = sim.get("bed.diameter") / 2 + c
    assert d.required_thickness == pytest.approx(6e5 * R / (S * E - 0.6 * 6e5) + c, rel=1e-12)
    assert not d.external_governs

    sim["wall.vacuum_service"] = True
    v = sim.wall_design()
    assert v.valid and v.external_thickness > 0
    assert v.required_thickness == max(v.internal_thickness, v.external_thickness)

    sim["wall.design_code"] = "asme_b31_3"
    assert sim.wall_design().valid


def test_wall_design_needs_a_design():
    d = skarstrom.Simulation(PEREZ).wall_design()
    assert not d.valid and "design pressure" in d.message


def test_operating_limits_of_a_vsa_run():
    sim = skarstrom.Simulation(PEREZ)
    sim["cycle.number_of_cycles"] = 2
    assert sim.operating_limits().cycle == -1  # no results yet
    assert sim.run().ok

    names = [s.name for s in sim.steps]
    r = sim.operating_limits()
    assert r.cycle == 1 and not r.design_checked
    assert [s.step_name for s in r.steps] == names
    assert all(math.isfinite(s.max_pressure) and math.isfinite(s.max_approach_to_fluidisation) for s in r.steps)
    assert min(s.min_pressure for s in r.steps) < 101325  # vacuum steps
    assert not [w for w in r.warnings if w.limit != "fluidisation"]

    # A design in use, but no vacuum service: the evacuation is reported.
    sim["wall.design_pressure"] = 10e5
    sim["wall.design_temperature"] = 400.0
    r = sim.operating_limits()
    assert r.design_checked
    vacuum = [w for w in r.warnings if w.limit == "vacuum"]
    assert vacuum and all(w.limit_value == 101325 for w in vacuum)
    assert all(names[w.step] in w.message for w in vacuum)

    sim["wall.vacuum_service"] = True
    assert not [w for w in sim.operating_limits().warnings if w.limit != "fluidisation"]


def test_warnings_can_be_switched_off():
    sim = skarstrom.Simulation(PEREZ)
    sim["cycle.number_of_cycles"] = 1
    sim["wall.design_pressure"] = 10e5
    sim["wall.design_temperature"] = 400.0
    assert sim.settings.warn_fluidisation and sim.get("settings.warn_vacuum")
    assert sim.run().ok
    limits = {w.limit for w in sim.operating_limits().warnings}
    assert {"fluidisation", "vacuum"} <= limits  # Perez 2019 A's pressurisation and evacuation

    sim.settings.warn_fluidisation = False
    sim["settings.warn_vacuum"] = False
    r = sim.operating_limits()
    assert not r.warnings
    assert max(s.max_approach_to_fluidisation for s in r.steps) > 1  # the summary is complete
