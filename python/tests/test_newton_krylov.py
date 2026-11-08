"""Newton-Krylov CSS acceleration from Python: the registry choice and the
diagnostics in the CSS history."""

import math

import skarstrom


def coarse_perez(cells=10):
    sim = skarstrom.examples.Perez2019_A.load()
    for i in range(len(sim.layers)):
        sim.layers[i].cells = cells
    return sim


def test_newton_krylov_is_a_choice():
    sim = coarse_perez()
    sim.settings.acceleration_method = "newton_krylov"
    assert sim.settings.acceleration_method == "newton_krylov"
    assert sim.get("settings.acceleration_method") == 3  # appended after Anderson


def test_newton_krylov_reports_its_solve():
    sim = coarse_perez()
    sim.settings.acceleration_method = "newton_krylov"
    result = sim.run(cycles=6)
    assert result.ok, result.message
    css = sim.css_history()
    solves = [c for c in css if c.newton_iterations >= 0]
    assert solves, "a solve runs from the fourth cycle"
    first = solves[0]
    assert first.cycle == 3
    assert first.newton_cycles > 0 and first.newton_outcome
    assert len(first.newton_residuals) == first.newton_iterations + 1
    assert len(first.newton_krylov_vectors) >= first.newton_iterations
    assert math.isnan(first.floquet_multiplier) or 0 <= first.floquet_multiplier < 1.5
    # Cycles without a solve keep the defaults.
    quiet = [c for c in css if c.newton_iterations < 0]
    assert all(c.newton_cycles == 0 and not c.newton_residuals for c in quiet)
