"""The cycle-state API: final_state, run_cycle, warm start (set_initial_state)."""

import copy

import numpy as np
import pytest

import skarstrom as sk

# A deck whose pressure-equalisation step replays a partner step that runs later
# in the cycle (the previous cycle's logged outlet), on a coarse grid with the
# CSS check off so every run does its cycles.
DECK = sk.examples.Vysyaraju2024_PSA_exp1


def make(cells=20):
    sim = DECK.load()
    for layer in sim.layers:
        layer.cells = cells
    sim.css.method = 0  # None
    return sim


def difference(a, b):
    """Largest relative difference over the differential entries of each field."""
    worst = 0.0
    ya, yb = a.y, b.y
    for f in b.fields:
        if not f["differential"]:
            continue
        idx = np.asarray(f["index"])
        scale = max(np.sqrt(np.mean(yb[idx] ** 2)), 1e-8)
        worst = max(worst, float(np.max(np.abs(ya[idx] - yb[idx])) / scale))
    return worst


@pytest.fixture(scope="module")
def states():
    straight = make()
    assert straight.run(cycles=3).ok
    first = make()
    assert first.run(cycles=2).ok
    return first.final_state(), straight.final_state()


def test_round_trip_matches_the_straight_run(states):
    x, reference = states
    assert x.cycles == 2 and reference.cycles == 3
    assert x.has_history
    sim = make()
    n_cycles = sim.cycle.number_of_cycles
    r = sim.run_cycle(x)
    assert r.ok and r.cycles_run == 1
    assert sim.cycle.number_of_cycles == n_cycles  # run_cycle leaves the deck alone
    gx = sim.final_state()
    assert gx.cycles == 3
    assert difference(gx, reference) < 1e-4


def test_repeated_run_cycle_steps_the_cycle_map(states):
    _, reference = states
    sim = make()
    for _ in range(3):
        assert sim.run_cycle().ok
    assert sim.final_state().cycles == 3
    assert difference(sim.final_state(), reference) < 1e-4


def test_warm_start(states):
    x, reference = states
    sim = make()
    assert sim.initial_state is None
    sim.set_initial_state(x)
    assert sim.initial_state is not None and sim.initial_state.cycles == 2
    assert sim.run(cycles=1).ok
    assert difference(sim.final_state(), reference) < 1e-4
    # warm_start() defaults to this run's final state
    s = sim.warm_start()
    assert s.cycles == 3 and sim.initial_state.cycles == 3
    sim.clear_initial_state()
    assert sim.initial_state is None


def test_state_contents(states):
    x, _ = states
    assert x.components == make().components
    assert [l["cells"] for l in x.layers] == [20] * len(x.layers)
    assert len(x.steps) == len(make().steps)
    assert len(x) == x.size == len(x.y) == len(x.yp)
    assert sum(len(f["index"]) for f in x.fields) == x.size
    assert any(e["history_points"] > 0 for e in x.endpoints)

    p = x.field("Pressure")
    assert p.shape[0] > 0 and np.all(p > 0)
    comp = x.components[0]
    assert x.field("Mole Fraction", comp).shape[0] > 0
    with pytest.raises(IndexError):
        x.field("No Such Field")

    d = x.differential
    assert 0 < d.size < x.size
    c = x.copy()
    c.differential = d * 1.0
    np.testing.assert_array_equal(c.y, x.y)
    c2 = copy.deepcopy(x)
    c2.differential = d * 0.5
    assert not np.array_equal(c2.y, x.y)
    assert np.array_equal(x.differential, d)  # copies are independent
    with pytest.raises(ValueError):
        c.differential = np.zeros(d.size + 1)
    with pytest.raises(ValueError):
        c.y = np.zeros(x.size + 1)


def test_incompatible_state_is_refused(states):
    x, _ = states
    sim = make()
    sim.add_component("Argon")
    why = sim.check_state(x)
    assert "components" in why
    with pytest.raises(ValueError, match="components"):
        sim.set_initial_state(x)
    assert make().check_state(x) == ""
    # A different cell count is not a refusal: the run remaps the state onto the deck's grid.
    assert make(cells=24).check_state(x) == ""


def test_no_final_state_before_a_run():
    sim = make()
    assert not sim.has_final_state
    with pytest.raises(RuntimeError):
        sim.final_state()
