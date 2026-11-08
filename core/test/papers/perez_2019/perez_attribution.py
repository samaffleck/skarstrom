#!/usr/bin/env python3
"""Evacuation-profile split attribution for the Perez 2019 decks.

Runs the six validation cases (A, B, C1, D1, E, F) of the test binary at each
split of the two-term evacuation profile and prints a Markdown table of the
differences to the experiment and to the paper's model (reference_values.csv).
The split is passed through ADSIM_PEREZ_EVAC_WEIGHT, which the fixture reads
(perez_2019.cpp); the decks and the fixture default (0.5) are not touched, and a
run with the variable set neither asserts tolerances nor re-saves a deck.

    python core/test/papers/perez_2019/perez_attribution.py \
        --binary build/linux-release/core/test/AdsorptionSimulatorTests \
        --splits 0.3 0.4 0.5 0.6 0.7 --jobs 4

Run it from `core/` (the tests' working directory). Productivity is the model's
CO2 mass productivity [kg/(kg adsorbent h)]; Table 3 values for it are not in
reference_values.csv, so it is shown without a difference. Energy is n/a for
decks without a pump or compressor model (these decks: all of them).
"""
import argparse
import concurrent.futures as cf
import csv
import math
import os
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = ["A", "B", "C1", "D1", "E", "F"]
GTEST = {"A": "ConfigA", "B": "ConfigB", "C1": "ConfigC1", "D1": "ConfigD1", "E": "ConfigE", "F": "ConfigF"}


def references():
    ref = {}
    with open(HERE / "reference_values.csv") as f:
        for row in csv.DictReader(f):
            ref[(row["case"], row["quantity"])] = float(row["value"])
    return ref


def run(binary, split, case, cwd):
    env = dict(os.environ, ADSIM_PEREZ_EVAC_WEIGHT=str(split))
    t0 = time.time()
    p = subprocess.run([binary, f"--gtest_filter=Perez2019Test.{GTEST[case]}"], cwd=cwd, env=env,
                       capture_output=True, text=True)
    wall = time.time() - t0
    m = re.search(r"\[ATTRIB\] split=(\S+) case=(\S+) cycles=(\d+) purity=(\S+) recovery=(\S+) "
                  r"productivity=(\S+) energy=(\S+)", p.stdout)
    if not m or p.returncode != 0:
        return dict(split=split, case=case, failed=True, tail=p.stdout[-400:], wall=wall)
    return dict(split=split, case=case, failed=False, cycles=int(m[3]), purity=float(m[4]),
                recovery=float(m[5]), productivity=float(m[6]), energy=float(m[7]), wall=wall)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", required=True)
    ap.add_argument("--splits", type=float, nargs="+", default=[0.3, 0.4, 0.5, 0.6, 0.7])
    ap.add_argument("--cases", nargs="+", default=CASES, choices=CASES)
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--cwd", default=".")
    a = ap.parse_args()
    ref = references()
    jobs = [(s, c) for c in a.cases for s in a.splits]
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        res = list(ex.map(lambda j: run(a.binary, j[0], j[1], a.cwd), jobs))
    print("| Case | split | CSS cycles | wall [s] | purity | recovery | productivity [kg/(kg h)] | energy [MJ/kg] "
          "| purity - exp | recovery - exp | purity - model | recovery - model |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in res:
        c = f"Config{r['case']}"
        if r["failed"]:
            print(f"| {r['case']} | {r['split']} | FAILED | {r['wall']:.0f} |" + " |" * 8)
            print(r["tail"], file=sys.stderr)
            continue
        d = lambda q, v: f"{v - ref[(c, q)]:+.2f}"
        en = "n/a" if math.isnan(r["energy"]) or r["energy"] == 0 else f"{r['energy']:.3f}"
        pr = "n/a" if math.isnan(r["productivity"]) else f"{r['productivity']:.4f}"
        print(f"| {r['case']} | {r['split']} | {r['cycles']} | {r['wall']:.0f} | {r['purity']:.2f} | "
              f"{r['recovery']:.2f} | {pr} | {en} | {d('CO2_purity_exp', r['purity'])} | "
              f"{d('CO2_recovery_exp', r['recovery'])} | {d('CO2_purity_model', r['purity'])} | "
              f"{d('CO2_recovery_model', r['recovery'])} |")


if __name__ == "__main__":
    main()
