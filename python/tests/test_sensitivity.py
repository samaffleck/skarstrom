"""Tests for skarstrom.sensitivity and ``skarstrom sensitivity``.

The statistics themselves are checked against closed-form answers in the
engine's test suite (core/test/unit/test_sensitivity.cpp); these tests cover
the Python surface on a real deck run for one cycle.
"""

import json
import math

import pytest

import skarstrom
from skarstrom import sensitivity as sa
from skarstrom import ureg as u
from skarstrom.cli import main

DECK = "Arora2021_LiLSX_PSA"


def first(sim, suffix, prefix):
    return next(p.path for p in sim.parameters(tunable_only=True) if p.path.startswith(prefix) and p.path.endswith(suffix))


@pytest.fixture(scope="module")
def sim():
    return skarstrom.examples.get(DECK).load()


@pytest.fixture(scope="module")
def paths(sim):
    return first(sim, ".time", "step["), first(sim, ".bed_void_fraction", "layer[")


def quick(sim, inputs, **kw):
    """One cycle per run, so no run reaches CSS: the requirement is off."""
    return sa.Study(sim, inputs, outputs=["Oxygen.purity", "Oxygen.recovery"], cycles=1, require_css=False, **kw)


def test_inputs():
    assert sa.uniform("a", 1, 2).as_dict() == {"path": "a", "distribution": "uniform", "lower": 1.0, "upper": 2.0,
                                                "mode": None, "mean": None, "sd": None}
    assert list(sa.uniform("a", 0, 10).quantile([0.25, 0.5])) == [2.5, 5.0]
    assert sa.normal("n", 5, 2).quantile(0.5)[0] == pytest.approx(5.0)
    assert sa.triangular("t", 0, 2, 10).quantile(0.2)[0] == pytest.approx(2.0)
    with pytest.raises(ValueError):
        sa.uniform("a", 2, 1).quantile(0.5)
    assert "purity" in sa.kpi_names() and "energy_consumption" in sa.kpi_names()


def test_study_resolves_and_checks_inputs(sim, paths):
    time, void = paths
    t = sim.get(time)
    study = quick(sim, [sa.relative(time, 0.1), sa.uniform(void, 0.30, 0.40)])
    assert study.inputs[0].lower == pytest.approx(0.9 * t) and study.inputs[0].upper == pytest.approx(1.1 * t)
    assert study.inputs[1].distribution == "uniform"
    # A quantity is converted to the parameter's SI unit.
    minutes = quick(sim, [sa.uniform(time, t / 60 * 0.9 * u.min, t / 60 * 1.1 * u.min)])
    assert minutes.inputs[0].lower == pytest.approx(0.9 * t)
    # Defaults: every tracked component's four headline KPIs.
    default = sa.Study(sim, [sa.relative(time, 0.1)], require_css=False)
    assert default.outputs == ["Oxygen.purity", "Oxygen.recovery", "Oxygen.productivity_mass",
                               "Oxygen.energy_consumption"]
    with pytest.raises(KeyError):
        quick(sim, [sa.uniform("step[999].time", 1, 2)])
    with pytest.raises(ValueError, match="cannot be varied"):
        quick(sim, [sa.uniform("settings.reltol", 1e-6, 1e-4)])
    with pytest.raises(ValueError, match="not tracked"):
        sa.Study(sim, [sa.relative(time, 0.1)], outputs=["Nitrogen.purity"], require_css=False)
    with pytest.raises(ValueError, match="CSS"):
        sa.Study(sim, [sa.relative(time, 0.1)], settings={"css.method": 0})
    with pytest.raises(ValueError):
        quick(sim, [sa.relative(time, 0.0)])


def test_screen(sim, paths, tmp_path):
    time, void = paths
    study = quick(sim, [sa.relative(time, 0.2), sa.uniform(void, 0.30, 0.40)])
    s = study.screen(trajectories=2)
    assert s.runs == 2 * 3 and s.failed == 0 and s.warnings == []
    table = s.table()
    assert len(table) == 2 * 2 and set(table["output"]) == {"Oxygen.purity", "Oxygen.recovery"}
    assert (table["mu_star"] >= 0).all() and (table["effects"] == 2).all()
    ranking = s.ranking("Oxygen.purity")
    assert {r[0] for r in ranking} == {time, void} and ranking[0][1] >= ranking[1][1]
    assert set(s.important(0.0)) == {time, void}
    assert s.important(1.0) and len(s.important(1.0)) <= 2
    assert len(s.samples_table()) == 6

    s.save(tmp_path / "screen.json")
    saved = json.loads((tmp_path / "screen.json").read_text())
    assert saved["method"] == "morris" and saved["study"]["policy"]["cycles"] == 1
    assert len(saved["samples"]) == 6

    # Parallel workers give the same numbers.
    parallel = quick(sim, [sa.relative(time, 0.2), sa.uniform(void, 0.30, 0.40)], workers=2).screen(trajectories=2)
    assert parallel.table()["mu_star"].tolist() == table["mu_star"].tolist()


def test_propagate(sim, paths):
    time, void = paths
    study = quick(sim, [sa.relative(time, 0.1), sa.uniform(void, 0.36, 0.40)])
    calls = []
    p = study.propagate(runs=4, batch=2, max_runs=6, progress=lambda done, total: calls.append((done, total)))
    assert p.runs == 6 and p.batches == 2 and not p.converged
    assert calls[-1] == (6, 6)
    table = p.table(nominal=True)
    assert list(table["output"]) == ["Oxygen.purity", "Oxygen.recovery"]
    row = table.iloc[0]
    assert row["p05"] <= row["p50"] <= row["p95"] and row["valid"] == 6
    assert 0 < row["nominal"] <= 1
    assert len(p.shares()) == 4
    bars = p.error_bars()
    assert bars["Oxygen.purity"][0] == row["p05"]
    # Only the inputs given are varied.
    one = study.propagate(runs=3, batch=1, max_runs=3, inputs=[time])
    assert one.runs == 3 and [i.path for i in one.inputs] == [time]
    with pytest.raises(KeyError):
        study.propagate(runs=3, max_runs=3, inputs=["layer[999].length"])


def test_cli(capsys, tmp_path, sim, paths):
    time, void = paths
    argv = ["sensitivity", DECK, "--input", f"{time}=10%", "--input", f"{void}=0.36:0.40", "--output",
            "Oxygen.purity", "--cycles", "1", "--no-require-css", "--trajectories", "1", "--runs", "3",
            "--max-runs", "3", "-q"]
    code = main(argv + ["--json"])
    out = json.loads(capsys.readouterr().out)
    assert code == 0
    assert out["screening"]["plan"]["trajectories"] == 1 and len(out["screening"]["samples"]) == 3
    assert len(out["propagation"]["samples"]) == 3
    code = main(argv + ["--method", "screen", "--save", str(tmp_path / "s.json")])
    text = capsys.readouterr().out
    assert code == 0 and "Screening (Morris)" in text and "mu*" in text
    assert json.loads((tmp_path / "s.json").read_text())["screening"]["method"] == "morris"
    with pytest.raises(SystemExit, match="expected 10%"):
        main(["sensitivity", DECK, "--input", f"{time}=lognormal:1:2", "--cycles", "1"])
    assert math.isfinite(out["propagation"]["result"]["outputs"][0]["p50"])
