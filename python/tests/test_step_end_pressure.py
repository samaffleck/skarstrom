"""A vacuum pump stream's ultimate pressure and the step-end pressure outputs.

The engine tests (core/test/unit/test_pump_ultimate_pressure.cpp) check the
model; these cover the Python surface on the Perez 2019 column, coarse, one
cycle, with the evacuation driven by a constant-power vacuum pump.
"""

import json
import math

import pytest

import skarstrom as sk
from skarstrom import optimisation as opt
from skarstrom import sensitivity as sa
from skarstrom import ureg as u

COARSE = {"layer[1].cells": 10}


def _pumped(ultimate=None):
    sim = sk.examples.Perez2019_A.load()
    sim.update(COARSE)
    values = {"discharge_pressure": 1e5}
    if ultimate is not None:
        values["ultimate_pressure"] = ultimate
    sim.steps["Evacuation"].bottom.set_streams("outlet", [sk.PumpStream(
        "Vacuum pump", power=10, pump=sk.Pump("vacuum_pump", efficiency_model="vacuum_pump", **values))])
    return sim


def _evacuation(sim):
    return next(r for r in sim.step_end_pressures() if r.step_name == "Evacuation")


def test_ultimate_pressure_is_a_pump_leaf_that_defaults_to_off():
    sim = _pumped()
    eid = sim.steps["Evacuation"].bottom[0].id
    assert sim.get(f"endpoint[{eid}].pump.ultimate_pressure") == 0.0
    sim.set(f"endpoint[{eid}].pump.ultimate_pressure", 10 * u.kPa)
    assert sim.get("step[Evacuation].bottom[Vacuum pump].pump.ultimate_pressure") == pytest.approx(1e4)
    assert "ultimate_pressure=" in sk.to_script(sim)


def test_the_pump_holds_its_ultimate_pressure():
    free = _pumped()
    assert free.run(cycles=1).ok
    floor = _pumped(ultimate=2e4)
    assert floor.run(cycles=1).ok
    assert _evacuation(free).bottom < 1e4
    end = _evacuation(floor)
    assert 2e4 <= end.bottom < 1.02 * 2e4


def test_step_end_pressures_match_the_series_and_feed_the_outputs():
    sim = _pumped()
    assert sim.run(cycles=1).ok
    ends = sim.step_end_pressures()
    assert [r.step_name for r in ends] == [s.name for s in sim.steps]
    last = ends[-1]
    bottom = sim.series("pressure", location="bottom")
    assert last.bottom == bottom.values[-1]
    assert last.time == bottom.time[-1]
    assert last.bed_min <= last.bed_mean <= last.bed_max
    assert set(last.as_dict()) >= {"step_name", "bottom", "top", "bed_min", "bed_mean", "bed_max"}
    assert sa.step_output_names() == ["end_pressure_bottom", "end_pressure_top", "end_pressure_min",
                                      "end_pressure_mean", "end_pressure_max"]


def test_optimiser_terms_take_a_step_reference():
    with pytest.raises(ValueError, match="unknown step KPI"):
        opt.optimise(_pumped(), {"step[Evacuation].time": (60, 120)},
                     objectives=[opt.maximise("step[Evacuation]", "purity")], evaluations=1)
    sim = _pumped(ultimate=2e4)
    r = opt.optimise(sim, {"step[Evacuation].bottom[Vacuum pump].power": (5, 15)},
                     objectives=[opt.minimise("step[Evacuation]", "end_pressure_bottom")],
                     constraints=[opt.at_least("step[Evacuation]", "end_pressure_bottom", 1e4, penalty=1e-3)],
                     evaluations=1, init_samples=1, cycles=1, rerun_cold=False)
    value = r.kpis["step[Evacuation]"]["end_pressure_bottom"]
    assert 2e4 <= value < 1.02 * 2e4
    assert r.cost == pytest.approx(value)


def test_cli_json_lists_the_step_end_pressures(tmp_path, capsys):
    from skarstrom.cli import main

    deck = tmp_path / "pumped.json"
    _pumped().save(str(deck))
    assert main(["run", str(deck), "--cycles", "1", "--json", "-q"]) == 0
    payload = json.loads(capsys.readouterr().out)
    names = [r["step_name"] for r in payload["step_end_pressures"]]
    assert "Evacuation" in names
    assert all(r["bottom"] is None or math.isfinite(r["bottom"]) for r in payload["step_end_pressures"])
