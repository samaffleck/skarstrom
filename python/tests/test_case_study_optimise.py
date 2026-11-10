"""The case-study optimisation script (examples/perez_2019_case_study_optimise.py).

Fast checks only: one evaluation of the base deck for a few cycles (every
constraint measured and written to results.csv), and a two-candidate search
that is resumed from its results file. The full search runs each candidate to
cyclic steady state (minutes per generation) and is not run here.
"""

import csv
import importlib
import math
import sys
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parents[2] / "examples"

pytestmark = pytest.mark.skipif(not EXAMPLES.is_dir(), reason="examples/ not found (not a repository checkout)")


@pytest.fixture
def script(monkeypatch, tmp_path):
    monkeypatch.syspath_prepend(str(EXAMPLES))
    monkeypatch.setenv("SKARSTROM_EXAMPLES_OUTPUT", str(tmp_path))
    for name in ("_common", "perez_2019_case_study_optimise"):
        sys.modules.pop(name, None)
    return importlib.import_module("perez_2019_case_study_optimise")


def _rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def test_base_point_is_evaluated_with_every_constraint(script, tmp_path):
    assert script.main(["--deck", "Base", "--base-only", "--cycles", "2", "--out", str(tmp_path)]) == 0
    rows = _rows(tmp_path / "results.csv")
    assert len(rows) == 1
    row = rows[0]
    assert row["deck"] == "Base" and row["study"] == "base" and row["status"] == "completed"
    # The decision variables are the deck's own values.
    assert float(row["variable[t_feed]"]) == pytest.approx(99.0)
    assert float(row["variable[vacuum_pump_power]"]) == pytest.approx(10.2)
    assert row["variable[cv_lpp]"] == ""      # a variable of the other deck
    for key in ("purity", "recovery", "energy_kwh_per_t", "productivity_kg_per_t_h", "max_fluidisation",
                "max_dp_per_m", "min_evacuation_pressure", "backflow_fraction", "wall_time_s"):
        assert math.isfinite(float(row[key])), key
    for name in script.CONSTRAINTS:
        assert row[f"margin_{name}"] != "", name
    # Two cycles do not reach CSS, so the run is infeasible whatever its KPIs.
    assert row["css_reached"] == "False" and row["margin_converged"] == "-1.0" and row["feasible"] == "False"
    # Margins are normalised by the bound: >= 0 when satisfied.
    floor = 10e3
    assert float(row["margin_evacuation_floor"]) == pytest.approx(
        (float(row["min_evacuation_pressure"]) - floor) / floor)
    assert float(row["margin_fluidisation"]) == pytest.approx(1 - float(row["max_fluidisation"]))


def test_constraint_sets():
    sys.path.insert(0, str(EXAMPLES))
    try:
        mod = importlib.import_module("perez_2019_case_study_optimise")
    finally:
        sys.path.remove(str(EXAMPLES))
    assert mod.active_constraints("energy", False) == (
        "purity", "recovery", "fluidisation", "pressure_drop", "evacuation_floor", "product", "converged")
    assert mod.active_constraints("purity", False, backflow=True) == (
        "fluidisation", "pressure_drop", "evacuation_floor", "backflow", "product", "converged")
    assert mod.active_constraints("purity", True) == ("product", "converged")
    # The registry paths of the decision variables exist in their decks.
    import skarstrom as sk

    for deck, variables in mod.VARIABLES.items():
        sim = sk.Simulation(str(EXAMPLES / "decks" / f"{mod.DECKS[deck]}.json"))
        x = mod.base_point(sim, deck)
        lower = [v[1] for v in variables]
        upper = [v[2] for v in variables]
        assert all(lo <= v <= hi for v, lo, hi in zip(x, lower, upper)), deck


def test_search_resumes_from_results(script, tmp_path, capsys):
    pytest.importorskip("pymoo")
    args = ["--deck", "Base", "--pop", "2", "--gen", "1", "--workers", "1", "--cycles", "2", "--no-plot", "-q",
            "--out", str(tmp_path)]
    assert script.main(args) == 0
    first = _rows(tmp_path / "results.csv")
    assert len(first) == 2 and {r["study"] for r in first} == {"energy"}
    # The same command again reads both candidates back instead of running them.
    assert script.main(args) == 0
    assert len(_rows(tmp_path / "results.csv")) == 2
    # One more generation runs only the new candidates.
    assert script.main(args[:5] + ["2"] + args[6:]) == 0
    assert len(_rows(tmp_path / "results.csv")) == 4
    assert "completed" in capsys.readouterr().out


def test_summary_tables_are_written_from_the_results_file(script, tmp_path):
    assert script.main(["--deck", "Base", "--base-only", "--cycles", "2", "--out", str(tmp_path)]) == 0
    summary = tmp_path / "summary"
    assert script.main(["--out", str(tmp_path), "--summary", str(summary)]) == 0
    fronts = _rows(summary / "fronts.csv")
    # The paper's points are always there; a run that did not reach CSS is on no front.
    assert {r["point"] for r in fronts if r["front"] == "paper"} == {"A", "B", "C1", "D1", "E", "F"}
    assert not [r for r in fronts if r["front"] != "paper"]
    cost = _rows(summary / "constraint_cost.csv")
    assert {r["constraints"] for r in cost} == {"none", "fluidisation", "pressure_drop", "evacuation_floor", "all"}
    counts = {r["deck"]: r for r in _rows(summary / "counts.csv")}
    assert counts["Base"]["evaluations"] == "1" and counts["Base"]["converged"] == "0"
