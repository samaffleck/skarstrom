"""Tests for skarstrom.optimisation and warm-started studies.

The loop itself is checked in the engine's test suite
(core/test/integration/test_warm_start_cache.cpp); these tests cover the
Python surface on the TSA example with a coarse grid, run to CSS.
"""

import math

import pytest

import skarstrom
from skarstrom import optimisation as opt
from skarstrom import sensitivity as sa

DECK = "Illustrative_TSA_13X"
COARSE = {"layer[1].cells": 6}


@pytest.fixture(scope="module")
def sim():
    return skarstrom.examples.get(DECK).load()


def test_optimise_warm_starts_and_reports_the_cold_optimum(sim):
    t = sim.get("step[1].time")
    seen = []
    r = opt.optimise(sim, {"step[1].time": (0.8 * t, 1.2 * t)},
                     objectives=[opt.maximise("CarbonDioxide", "recovery")],
                     constraints=[opt.at_least("CarbonDioxide", "purity", 0.5, penalty=10)],
                     evaluations=3, init_samples=2, settings=COARSE, progress=seen.append)
    assert len(r.evaluations) == 3 and len(seen) == 4  # and the cold re-run
    assert not r.evaluations[0]["warm"]
    assert all(e["warm"] and e["css_reached"] for e in r.evaluations[1:])
    assert r.cold is not None and not r.cold["warm"] and r.cold["index"] == -1
    assert r.cost == r.cold["cost"]
    assert r.cost == pytest.approx(r.search_cost, abs=100 * sim.get("css.state_tolerance"))
    assert 0.8 * t <= r.x["step[1].time"] <= 1.2 * t
    assert 0 < r.kpis["CarbonDioxide"]["recovery"] <= 1
    table = r.table()
    assert list(table["index"]) == [0, 1, 2, -1]
    assert {"step[1].time", "cost", "warm", "cycles_run", "CarbonDioxide.recovery"} <= set(table.columns)
    copy = sim.copy()
    r.apply(copy)
    assert copy.get("step[1].time") == pytest.approx(r.x["step[1].time"])
    assert sim.get("step[1].time") == t  # optimise leaves the deck alone

    cold = opt.optimise(sim, {"step[1].time": (0.8 * t, 1.2 * t)},
                        objectives=[opt.maximise("CarbonDioxide", "recovery")],
                        evaluations=2, init_samples=2, settings=COARSE, warm_start=False, rerun_cold=False)
    assert not any(e["warm"] for e in cold.evaluations) and cold.cold is None
    assert cold.cost == cold.search_cost


def test_optimise_checks_its_set_up(sim):
    obj = [opt.maximise("CarbonDioxide", "recovery")]
    with pytest.raises(ValueError):
        opt.optimise(sim, {}, obj)
    with pytest.raises(ValueError):
        opt.optimise(sim, {"step[1].time": (1, 2)}, [])
    with pytest.raises(ValueError, match="unknown KPI"):
        opt.optimise(sim, {"step[1].time": (1, 2)}, [opt.maximise("CarbonDioxide", "happiness")])
    with pytest.raises(ValueError, match="no such parameter"):
        opt.optimise(sim, {"step[999].time": (1, 2)}, obj)
    with pytest.raises(ValueError, match="lower < upper"):
        opt.optimise(sim, {"step[1].time": (2, 1)}, obj)
    with pytest.raises(ValueError, match="not tracked"):
        opt.optimise(sim, {"step[1].time": (1, 2)}, [opt.maximise("Nitrogen", "purity")])


def test_state_cache_and_run_warm(sim):
    t = sim.get("step[1].time")
    cache = skarstrom.StateCache([0.8 * t], [1.2 * t])
    assert len(cache) == 0 and cache.dimension == 1 and cache.nearest([t]) is None
    a = sim.copy()
    a.set("layer[1].cells", 6)
    run, info = opt.run_warm(a, cache, [t])
    assert run.css_reached and not info["warm"] and info["cached"] and len(cache) == 1
    state, distance = cache.nearest([1.1 * t])
    assert isinstance(state, skarstrom.CycleState) and distance == pytest.approx(0.25)
    b = sim.copy()
    b.set("layer[1].cells", 6)
    b.set("step[1].time", 1.1 * t)
    run2, info2 = opt.run_warm(b, cache, [1.1 * t])
    assert info2["warm"] and info2["note"] == "" and run2.cycles_run < run.cycles_run
    assert info2["start_cycles"] == run.cycles_run and len(cache) == 2
    assert cache.normalise([t]) == pytest.approx([0.5])


def test_study_warm_starts_its_runs(sim):
    inputs = [sa.relative("step[1].time", 0.1)]
    outputs = ["CarbonDioxide.purity", "CarbonDioxide.recovery"]
    warm = sa.Study(sim, inputs, outputs=outputs, settings=COARSE)
    cold = sa.Study(sim, inputs, outputs=outputs, settings=COARSE, warm_start=False)
    assert warm.warm_start and warm.describe()["warm_start"] and not cold.warm_start
    assert "warm_start" not in warm.setup()  # a run option, not saved with the deck
    ws = warm.screen(trajectories=1)
    cs = cold.screen(trajectories=1)
    assert [s["warm"] for s in ws.samples] == [False, True]
    assert not any(s["warm"] for s in cs.samples)
    assert ws.samples[1]["cycles_run"] < cs.samples[1]["cycles_run"]
    tol = 100 * sim.get("css.state_tolerance")
    for a, b in zip(ws.samples, cs.samples):
        for x, y in zip(a["outputs"], b["outputs"]):
            assert math.isfinite(x) and abs(x - y) <= tol
    assert {"warm", "cycles_run", "note"} <= set(ws.samples_table().columns)
