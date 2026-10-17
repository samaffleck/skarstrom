"""Smoke test for the example scripts in examples/ (a repository checkout).

Every script is imported, and each runs once with fast options: one
breakthrough case, or two cycles of one cyclic case. The full runs (to
cyclic steady state) are timed in examples/README.md.
"""

import importlib
import os
import sys
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parents[2] / "examples"

FAST = {
    "ward_2022": ["--case", "100"],
    "casas_2012": ["--case", "Fig4_50_50"],
    "nguyen_2022": ["--case", "Fig5"],
    "streb_2021": ["--case", "A1"],
    "arora_2021": ["--case", "LiLSX_PSA", "--cycles", "2"],
    "ferre_2026": ["--cycles", "2"],
    "nguyen_2023": ["--case", "I", "--cycles", "2"],
    "perez_2019": ["--case", "A", "--cycles", "2"],
    "vysyaraju_2024": ["--cycles", "2"],
    "ward_pini_2024": ["--cycles", "2"],
    "tsa_illustrative": ["--cycles", "2"],
    "five_step_product_purge": ["--cycles", "2"],
    "optimise_nsga2": ["--pop", "4", "--gen", "1", "--workers", "2", "--cycles", "2"],
}
NEEDS = {"optimise_nsga2": "pymoo"}

pytestmark = pytest.mark.skipif(not EXAMPLES.is_dir(), reason="examples/ not found (not a repository checkout)")


@pytest.fixture(autouse=True)
def _examples_on_path(monkeypatch, tmp_path):
    monkeypatch.syspath_prepend(str(EXAMPLES))
    monkeypatch.setenv("SKARSTROM_EXAMPLES_OUTPUT", str(tmp_path))


def test_every_script_is_covered():
    scripts = {p.stem for p in EXAMPLES.glob("*.py") if not p.stem.startswith("_")}
    assert scripts == set(FAST), "add new example scripts to FAST"


@pytest.mark.parametrize("name", sorted(FAST))
def test_example_runs(name, capsys):
    if name in NEEDS:
        pytest.importorskip(NEEDS[name])
    sys.modules.pop("_common", None)
    mod = importlib.import_module(name)
    assert mod.main(FAST[name] + ["--no-plot", "-q"]) == 0
    out = capsys.readouterr().out
    assert "completed" in out or "reached_css" in out
