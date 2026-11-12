"""Tests for skarstrom.calibration (the engine's calibration object from Python).

The solver itself is checked in core/test/unit/test_calibration.cpp; these
tests cover the Python surface. The deck is the diluted-bed breakthrough
example on 12 cells (a run takes a fraction of a second) and the data are
synthetic: the deck is run at known values, its outlet CO2 is read at the
measured times and seeded noise is added, so a fit must recover the values.
"""

import json
import math

import numpy as np
import pytest

import skarstrom
from skarstrom import calibration as cal
from skarstrom import cli

DECK = "Illustrative_Diluted_Bed"
OUTLET = "outlet_CO2_mole_fraction[-]_Qin=100[ccm]"
CAPACITY = "adsorbent[1].isotherm[CarbonDioxide].m2"  # 6.04 mol/kg in the deck
AFFINITY = "adsorbent[1].isotherm[CarbonDioxide].b2"  # 3.32e-6 m3/mol
FLOW = "endpoint[1].volumetric_flow"
NOISE = 0.01


def synthetic(values, seed, flow_factor=1.0):
    """The deck on 12 cells, its outlet CO2 dataset replaced by the simulated
    curve at `values` plus seeded noise (with that noise as its uncertainty).
    The returned deck keeps its own values (6.04 and 3.32e-6)."""
    sim = skarstrom.examples.get(DECK).load()
    sim["layer[1].cells"] = 12
    sim[FLOW] = sim[FLOW] * flow_factor
    truth = sim.copy()
    for path, v in values.items():
        truth[path] = v
    assert truth.run().ok
    data = sim.experimental_data[OUTLET]
    xs, ys = data.simulated(truth)
    rng = np.random.default_rng(seed)
    y = np.interp(data.x, xs, ys) + NOISE * rng.standard_normal(len(data.x))
    observable = data.observable
    sim.experimental_data.remove(OUTLET)
    sim.experimental_data.add(OUTLET, data.x, y, observable=observable, uncertainty=[NOISE] * len(y))
    return sim


@pytest.fixture(scope="module")
def one_deck():
    sim = synthetic({CAPACITY: 6.5}, seed=3)
    result = sim.calibrate([cal.Parameter(CAPACITY, 4.0, 8.0, start=5.5)], datasets=[OUTLET],
                           max_simulations=30, threads=2)
    return sim, result


def test_one_deck_fit_recovers_the_parameter(one_deck):
    sim, r = one_deck
    assert r.ok and not r.cancelled, r.message
    (p,) = r.parameters
    assert p.path == CAPACITY and p.case is None and p.unit == "mol/kg"
    assert math.isfinite(p.std_error) and p.std_error > 0
    assert abs(p.value - 6.5) < max(4 * p.std_error, 0.05)
    assert p.start == 5.5 and (p.lower, p.upper) == (4.0, 8.0) and not p.at_bound
    assert r[CAPACITY] is p and r.values == {CAPACITY: p.value}
    assert r.correlation.shape == (1, 1)
    # The trace goes down; the counters are filled.
    trace = r.objective_trace
    assert len(trace) >= 2 and trace[-1] < trace[0]
    assert 0 < r.simulations <= 31 and r.iterations > 0 and r.wall_time_s > 0
    assert r.method == "least_squares" and r.stop_reason in ("converged", "small_step", "budget", "no_improvement")
    (d,) = r.datasets
    assert d is r.dataset(OUTLET)
    assert d.name == OUTLET and d.points == len(d.x) == len(d.residuals) > 10
    np.testing.assert_allclose(d.residuals, d.simulated - d.measured)
    assert d.rmse < 3 * NOISE and len(d.series_x) == len(d.series_y) > len(d.x)
    text = r.table()
    assert CAPACITY in text and OUTLET in text and "std error" in text
    # The deck passed is not changed by the fit.
    assert sim[CAPACITY] == pytest.approx(6.04)


def test_apply_and_calibrated(one_deck):
    sim, r = one_deck
    value = r.parameters[0].value
    calibrated = r.calibrated()
    assert isinstance(calibrated, skarstrom.Simulation)
    assert calibrated[CAPACITY] == pytest.approx(value)
    other = sim.copy()
    assert r.apply(other) is other and other[CAPACITY] == pytest.approx(value)
    assert sim[CAPACITY] == pytest.approx(6.04)
    target = sim.copy()
    r2 = cal.CalibrationResult(r._data, [cal.Case(target)])
    r2.apply()  # the case's own deck
    assert target[CAPACITY] == pytest.approx(value)


def test_to_json_round_trip(one_deck, tmp_path):
    _, r = one_deck
    path = tmp_path / "result.json"
    text = r.to_json(str(path))
    json.loads(text)  # strict JSON (NaN written as null)
    back = cal.CalibrationResult.from_json(str(path))
    assert back.to_dict() == r.to_dict()
    assert back.parameters == r.parameters
    np.testing.assert_array_equal(back.datasets[0].residuals, r.datasets[0].residuals)
    assert back.calibrated()[CAPACITY] == pytest.approx(r.parameters[0].value)
    with pytest.raises(ValueError, match="from JSON"):
        back.apply()


def test_two_cases_with_a_shared_and_a_per_case_parameter():
    a = synthetic({CAPACITY: 6.5, AFFINITY: 3.0e-6}, seed=11)
    b = synthetic({CAPACITY: 6.5, AFFINITY: 4.0e-6}, seed=12, flow_factor=2.0)
    r = skarstrom.calibrate(
        [cal.Case(a, [OUTLET], name="Q100"), (b, [OUTLET], 1.0)],
        [cal.Parameter(CAPACITY, 4.0, 8.0),
         cal.Parameter(AFFINITY, 1e-6, 2e-5, scale="log", shared=False)],
        max_simulations=80, threads=3, differences="central")
    assert r.ok, r.message
    assert [p.name for p in r.parameters] == [CAPACITY, f"{AFFINITY} [Q100]", f"{AFFINITY} [case 2]"]
    assert [p.case for p in r.parameters] == [None, 0, 1]
    assert r.case_names == ["Q100", "case 2"]
    assert r.correlation.shape == (3, 3)
    assert len(r.datasets) == 2 and {d.case for d in r.datasets} == {0, 1}
    assert r[CAPACITY].value == pytest.approx(6.5, rel=0.1)
    with pytest.raises(KeyError, match="per case"):
        r[AFFINITY]
    # Each case's deck gets the shared value and its own per-case value.
    decks = r.calibrated()
    assert len(decks) == 2
    for i, deck in enumerate(decks):
        assert deck[CAPACITY] == pytest.approx(r.parameters[0].value)
        assert deck[AFFINITY] == pytest.approx(r.parameters[1 + i].value)
    r.apply(1)
    assert b[AFFINITY] == pytest.approx(r.parameters[2].value)


def test_global_with_a_tiny_budget():
    sim = synthetic({CAPACITY: 6.5}, seed=5)
    seen = []
    r = sim.calibrate([(CAPACITY, 4.0, 8.0)], method="global", global_evaluations=4, threads=2,
                      progress=seen.append)
    assert r.ok and r.method == "global" and r.stop_reason == "global_done"
    assert len(seen) == 4 and all(s["solver"] == "global" for s in seen)
    assert 4 <= r.simulations <= 6  # and the Jacobian for the standard error
    assert 4.0 <= r.parameters[0].value <= 8.0


def test_progress_false_and_cancel_token_stop_the_fit():
    sim = synthetic({CAPACITY: 6.5}, seed=7)
    r = sim.calibrate([cal.Parameter(CAPACITY, 4.0, 8.0, start=5.0)], progress=lambda p: False, threads=1)
    assert r.cancelled and r.stop_reason == "cancelled"
    token = cal.CancelToken()
    r = sim.calibrate([cal.Parameter(CAPACITY, 4.0, 8.0, start=5.0)], cancel=token, threads=1,
                      progress=lambda p: token.cancel())
    assert token.cancelled and r.cancelled


def test_progress_errors_are_raised():
    sim = synthetic({CAPACITY: 6.5}, seed=7)

    def boom(p):
        raise RuntimeError("from progress")

    with pytest.raises(RuntimeError, match="from progress"):
        sim.calibrate([(CAPACITY, 4.0, 8.0)], progress=boom, threads=1)


def test_invalid_problems_raise_before_any_run():
    sim = skarstrom.examples.get(DECK).load()
    with pytest.raises(ValueError, match="bounds"):
        sim.calibrate([(CAPACITY, 8.0, 4.0)])
    with pytest.raises(ValueError, match="no dataset"):
        sim.calibrate([(CAPACITY, 4.0, 8.0)], datasets=["no such dataset"])
    with pytest.raises(KeyError, match="unknown parameter"):
        sim.calibrate([("layer[1].no_such_thing", 0, 1)])
    with pytest.raises(ValueError, match="method"):
        sim.calibrate([CAPACITY], method="newton")
    unbound = sim.copy()
    for s in unbound.experimental_data:
        s.observable = None
    with pytest.raises(ValueError, match="bound"):
        unbound.calibrate([(CAPACITY, 4.0, 8.0)])


def test_bounds_and_start_take_units():
    sim = skarstrom.examples.get(DECK).load()
    t = sim.get("step[1].time")
    d = cal._parameter_to_engine(cal.Parameter("step[1].time", "10 min", 3600, start=40, unit=None), sim)
    assert (d["lower"], d["upper"], d["initial"]) == (600.0, 3600.0, 40.0)
    d = cal._parameter_to_engine(cal.Parameter("step[1].time", 10, 60, unit="min"), sim)
    assert (d["lower"], d["upper"]) == (600.0, 3600.0)
    d = cal._parameter_to_engine(
        cal.Parameter("step[1].time", 10 * skarstrom.ureg.min, "1 [h]", start="50 min"), sim)
    assert (d["lower"], d["upper"], d["initial"]) == (600.0, 3600.0, 3000.0)
    with pytest.raises(ValueError):
        cal._parameter_to_engine(cal.Parameter("step[1].time", "2 bar", 60), sim)
    # A bare path takes the registry's suggested bounds.
    d = cal._parameter_to_engine(cal.Parameter(CAPACITY), sim)
    bounds = sim._index().find(CAPACITY).bounds
    assert (d["lower"], d["upper"]) == bounds and t > 0


def test_cli_calibrate(tmp_path, capsys):
    deck = tmp_path / "synthetic.json"
    synthetic({CAPACITY: 6.5}, seed=9).save(str(deck))
    out = tmp_path / "calibrated.json"
    js = tmp_path / "result.json"
    code = cli.main(["calibrate", str(deck), "--param", f"{CAPACITY}:4:8", "--dataset", OUTLET,
                     "--max-simulations", "20", "--threads", "2", "--out", str(out), "--json", str(js), "-q"])
    assert code == 0
    printed = capsys.readouterr().out
    assert CAPACITY in printed and "std error" in printed
    r = cal.CalibrationResult.from_json(str(js))
    assert r.ok and r.parameters[0].value == pytest.approx(skarstrom.Simulation(str(out))[CAPACITY])
    assert skarstrom.Simulation(str(deck))[CAPACITY] == pytest.approx(6.04)


def test_cli_parsing():
    p = cli._parse_param(f"{CAPACITY}:1e-6:2e-5:log:per_case")
    assert (p.path, p.lower, p.upper, p.scale, p.shared) == (CAPACITY, "1e-6", "2e-5", "log", False)
    p = cli._parse_param(CAPACITY)
    assert (p.path, p.lower, p.upper, p.scale, p.shared) == (CAPACITY, None, None, "linear", True)
    assert cli._parse_dataset("a_b[-]:2") == ("a_b[-]", 2.0)
    assert cli._parse_dataset("a_b[-]") == ("a_b[-]", None)


def test_cli_failed_fit_exits_non_zero(tmp_path, capsys):
    deck = tmp_path / "unbound.json"
    sim = skarstrom.examples.get(DECK).load()
    for s in sim.experimental_data:
        s.observable = None
    sim.save(str(deck))
    assert cli.main(["calibrate", str(deck), "--param", f"{CAPACITY}:4:8", "-q"]) == 1
    assert "bound" in capsys.readouterr().err
