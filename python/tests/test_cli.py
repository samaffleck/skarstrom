"""Tests for the ``skarstrom`` command-line tool."""

import csv
import json
import os
import subprocess
import sys

import pytest

import skarstrom
from skarstrom.cli import main

PEREZ = "Perez2019_A"


def run_cli(capsys, *argv):
    code = main(list(argv))
    out, err = capsys.readouterr()
    return code, out, err


def test_version(capsys):
    with pytest.raises(SystemExit) as e:
        main(["--version"])
    assert e.value.code == 0
    assert skarstrom.__version__ in capsys.readouterr().out


def test_examples(capsys, tmp_path):
    code, out, _ = run_cli(capsys, "examples")
    assert code == 0 and PEREZ in out.split()
    code, out2, _ = run_cli(capsys, "examples", "list")
    assert code == 0 and out2 == out
    code, out, _ = run_cli(capsys, "examples", "--json")
    rows = json.loads(out)
    assert code == 0 and {"name", "family", "path"} <= set(rows[0])
    assert any(r["name"] == PEREZ and r["family"] == "Perez2019" for r in rows)
    code, out, _ = run_cli(capsys, "examples", "path", "Perez2019_A")
    assert code == 0 and os.path.isfile(out.strip())
    code, out, _ = run_cli(capsys, "examples", "copy", PEREZ, str(tmp_path))
    assert code == 0 and (tmp_path / (PEREZ + ".json")).is_file()


def test_examples_run(capsys):
    code, out, err = run_cli(capsys, "examples", "run", "Perez2019_A", "--cycles", "1", "--json", "-q")
    assert code == 0, err
    payload = json.loads(out)
    assert payload["cycles_run"] == 1 and payload["status"] in ("completed", "reached_css")


def test_params_lists_tunable_and_all(capsys):
    deck = skarstrom.examples.get(PEREZ).path
    code, out, _ = run_cli(capsys, "params", deck)
    assert code == 0
    assert "step[1].time" in out and "settings.reltol" not in out
    code, out, _ = run_cli(capsys, "params", deck, "--all")
    assert "settings.reltol" in out and "info.title" in out
    code, out, _ = run_cli(capsys, "params", deck, "--json", "--search", "time")
    rows = json.loads(out)
    assert rows and all("time" in (r["path"] + r["alias"] + r["name"]).lower() for r in rows)


def test_run_prints_kpis_and_css(capsys, tmp_path):
    code, out, err = run_cli(
        capsys, "run", PEREZ, "--cycles", "2", "--set", "step[1].time=40", "--set", "info.title=a=b",
        "--csv", str(tmp_path), "-q",
    )
    assert code == 0, out + err
    assert "Status:" in out and "CSS:" in out
    assert "KPIs at cycle 2" in out and "CarbonDioxide" in out
    with open(tmp_path / "kpis.csv") as f:
        rows = list(csv.DictReader(f))
    assert [r["cycle"] for r in rows] == ["1", "2"]
    assert float(rows[-1]["purity [mol/mol]"]) > 0
    run = json.loads((tmp_path / "run.json").read_text())
    assert run["cycles_run"] == 2
    assert (tmp_path / "css.csv").exists() and (tmp_path / "logs").is_dir()


def test_run_json(capsys):
    code, out, _ = run_cli(capsys, "run", PEREZ, "--cycles", "1", "--json", "-q")
    assert code == 0
    payload = json.loads(out)
    assert payload["cycles_run"] == 1 and payload["kpis"][0]["component"] == "CarbonDioxide"
    assert payload["warnings"] == []  # no wall design, no fluidisation
    code, out, _ = run_cli(capsys, "run", PEREZ, "--cycles", "1", "--set", "wall.design_pressure=1e5", "--json", "-q")
    assert any("vacuum service" in w for w in json.loads(out)["warnings"])


def test_run_rejects_unknown_parameter_and_bad_value(capsys):
    with pytest.raises(SystemExit) as e:
        main(["run", PEREZ, "--set", "no.such=1", "-q"])
    assert "unknown parameter" in str(e.value.code)
    with pytest.raises(SystemExit) as e:
        main(["run", PEREZ, "--set", "cycle.number_of_cycles=many", "-q"])
    assert "cannot set" in str(e.value.code)
    with pytest.raises(SystemExit):
        main(["run", "missing-deck.json"])


def test_console_script_entry_point(tmp_path):
    exe = [sys.executable, "-m", "skarstrom"]
    r = subprocess.run(exe + ["run", PEREZ, "--cycles", "1", "-q"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "KPIs at cycle 1" in r.stdout
