#!/usr/bin/env python3
"""Calibration of the Streb & Mazzotti (2021) decks on experiment A3 (README, "Calibration on A3").

Builds every Table 2 case like streb_2021.cpp, with the calibrated values as
parameters, and runs the three stages of the calibration:

    flows      flow factor of each case from the speed of its CO2 temperature
               front between the thermocouples at 35 and 110 cm
    fit        weighted least squares on A3 alone: LDF coefficients of CO2 and
               CH4, wall-to-thermofluid coefficient and sampling-line volume,
               with A3's flow factor held at 1.02 (its thermocouples give 1.033;
               +2 % of the set point is the bound)
    predict    every case with the calibrated values and its flow factor;
               --clip 0.02 holds the flow factors to +-2 % of the set point

    python core/test/papers/streb_2021/streb_calibration.py flows
    python core/test/papers/streb_2021/streb_calibration.py fit
    python core/test/papers/streb_2021/streb_calibration.py predict [--clip 0.02]

Needs the skarstrom package (pip install .) and scipy. The values printed by
`fit` and `flows` are the ones written into streb_2021.cpp and
examples/streb_2021.py. A fit takes about 30 A3 runs (10 s each, --jobs 4).
"""
from __future__ import annotations

import argparse
import json
import math
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

import skarstrom as sk
from skarstrom import ureg as u

REPO = Path(__file__).resolve().parents[4]
DATA = REPO / "core" / "data" / "streb_2021"
Z_CM = [10, 35, 60, 85, 110]
LENGTH = 1.2
AREA = math.pi * 0.025 ** 2 / 4

# name: (P [Pa], T [K], Q [m3/s at P, T; Table 2], ternary, duration [s])
CASES = {
    "A3": (10e5, 298.0, 20e-6, False, 4100.0),
    "A1": (1e5, 298.0, 20e-6, False, 4100.0),
    "A2": (10e5, 298.0, 10e-6, False, 1500.0),
    "A4": (10e5, 298.0, 40e-6, False, 1000.0),
    "A5": (25e5, 298.0, 20e-6, False, 750.0),
    "A6": (10e5, 318.15, 20e-6, False, 1050.0),
    "A7": (10e5, 338.15, 20e-6, False, 1100.0),
    "B1": (1e5, 298.0, 20e-6, True, 10400.0),
    "B2": (10e5, 298.0, 10e-6, True, 3350.0),
}
# The decks' flow factors (streb_2021.cpp), from `flows` with the calibrated values.
FLOW_FACTORS = {"A3": 1.020, "A1": 1.051, "A2": 0.939, "A4": 0.906, "A5": 0.978,
                "A6": 0.965, "A7": 0.995, "B1": 1.057, "B2": 0.887}

SKELETAL, PARTICLE, BED = 2359.0, 1085.0, 595.0
EP = 1 - PARTICLE / SKELETAL
EB = 1 - BED / PARTICLE
SIGMA_Y, SIGMA_T = 0.03, 3.0          # residual weights: mole fraction, thermocouple [K]


@dataclass
class Params:
    k_co2: float = 0.0592              # [1/s] (paper 0.06)
    k_ch4: float = 0.137               # [1/s] (paper 0.2)
    h_out: float = 49.7                # [W/(m2 K)] (paper 220)
    line_volume: float = 1.42e-4       # [m3] sampling line, gas volume at column pressure
    flow_factor: float = 1.0
    line_cells: int = 40


def build(case: str, p: Params) -> sk.Simulation:
    P, T, Q, ternary, duration = CASES[case]
    comps = ["CarbonDioxide", "Methane"] + (["Hydrogen"] if ternary else []) + ["Helium"]
    sim = sk.Simulation.new(title=f"Streb 2021 {case}", description="calibration", author="streb_calibration.py")
    for c in comps:
        sim.add_component(c)
    sim.fluid.adsorbed_heat_capacity_mode = "constant"
    sim.fluid.heat_of_adsorption_convention = "isosteric"
    sim.fluid.gas_heat_capacity_mode = "constant"
    for c in comps:
        sim.fluid.adsorbed_heat_capacity[c] = 0 * u.J / (u.mol * u.K)
        sim.fluid.gas_heat_capacity[c] = sk.ideal_gas_heat_capacity(c, T)

    z = sim.add_adsorbent("13X", solid_density=PARTICLE * u.kg / u.m**3, heat_capacity=1100 * u.J / (u.kg * u.K),
                          thermal_conductivity=0.15 * u.W / (u.m * u.K), particle_void_fraction=EP)
    iso = {
        "CarbonDioxide": {"A1": 7.268 * u.mol / u.kg, "A2": -0.6168, "B1": 1.129e-4 / u.bar, "B2": 28.389 * u.kJ / u.mol,
                          "C1": 0.425, "C2": 0.7238, "T_ref": 298 * u.K},
        "Methane": {"A1": 4.473 * u.mol / u.kg, "A2": -0.6569, "B1": 3.605e-4 / u.bar, "B2": 15.561 * u.kJ / u.mol,
                    "C1": 1.073, "C2": 0.0, "T_ref": 298 * u.K},
    }
    if ternary:
        iso["Hydrogen"] = {"A1": 5.013 * u.mol / u.kg, "A2": 0.0, "B1": 1.034e-4 / u.bar, "B2": 9.453 * u.kJ / u.mol,
                           "C1": 1.006, "C2": 0.0, "T_ref": 298 * u.K}
    z.set_isotherm("extended_sips_2", iso, components=list(iso))
    z.set_isotherm("inert", components=["Helium"])
    z.set_rate_model("CarbonDioxide", "ldf", "constant", k=p.k_co2 / u.s)
    z.set_rate_model("Methane", "ldf", "constant", k=p.k_ch4 / u.s)
    if ternary:
        z.set_rate_model("Hydrogen", "ldf", "constant", k=1.0 / u.s)
    z.set_rate_model("Helium", "ldf", "constant", k=0.5 / u.s)
    for c, hoa in (("CarbonDioxide", 36.3), ("Methane", 18.8), ("Hydrogen", 10.7), ("Helium", 0.0)):
        if c in comps:
            z.heat_of_adsorption[c] = hoa * u.kJ / u.mol

    sim.bed.set_shape("cylinder", diameter=2.5 * u.cm)
    sim.bed.ambient_temperature = T * u.K
    sim.bed.set_thermal_equilibrium("lte")
    sim.bed.set_particle_heat_transfer("constant", particle_htc=100 * u.W / (u.m**2 * u.K))
    for void in (sim.bed.bottom_void, sim.bed.top_void):
        void.enabled = False
    sim.wall.set_material("ss316l", density=8000 * u.kg / u.m**3, heat_capacity=500 * u.J / (u.kg * u.K),
                          thermal_conductivity=0 * u.W / (u.m * u.K))
    sim.wall.thickness = 2.5 * u.mm
    sim.wall.set_boundary("zero_flux")

    layer = sim.add_layer("13X", composition={"13X": 1.0}, momentum="kozeny_carman", length=LENGTH * u.m,
                          cells=160, bed_void_fraction=EB, particle_diameter=1.8 * u.mm)
    layer.set_dispersion("constant", dispersion_coefficient=0 * u.m**2 / u.s)
    layer.set_heat_dispersion("constant", heat_dispersion_coefficient=0 * u.W / (u.m * u.K))
    layer.set_internal_heat_transfer("dixon")
    layer.set_external_heat_transfer("constant", external_htc=p.h_out * u.W / (u.m**2 * u.K))

    # Sampling line: inert, nearly empty, adiabatic layer above the bed (plug flow).
    line_solid = sim.add_adsorbent("Sampling line", solid_density=1000 * u.kg / u.m**3,
                                   heat_capacity=500 * u.J / (u.kg * u.K),
                                   thermal_conductivity=0.15 * u.W / (u.m * u.K), particle_void_fraction=0.0)
    line_solid.set_isotherm("inert", components=comps)
    for c in comps:
        line_solid.set_rate_model(c, "ldf", "constant", k=0.5 / u.s)
        line_solid.heat_of_adsorption[c] = 0 * u.J / u.mol
    line = sim.add_layer("Sampling line", composition={"Sampling line": 1.0}, momentum="kozeny_carman",
                         length=p.line_volume / (0.99 * AREA) * u.m, cells=p.line_cells,
                         bed_void_fraction=0.99, particle_diameter=10 * u.mm)
    line.set_dispersion("constant", dispersion_coefficient=0 * u.m**2 / u.s)
    line.set_heat_dispersion("constant", heat_dispersion_coefficient=0 * u.W / (u.m * u.K))
    line.set_internal_heat_transfer("constant", internal_htc=0 * u.W / (u.m**2 * u.K))
    line.set_external_heat_transfer("constant", external_htc=0 * u.W / (u.m**2 * u.K))

    ic = sim.initial_conditions
    ic.pressure = P * u.Pa
    ic.temperature = T * u.K
    ic.wall_temperature = T * u.K
    ic.velocity = 0 * u.m / u.s
    for c in comps:
        ic.y[c] = 1.0 if c == "Helium" else 0.0

    sim.cycle.number_of_cycles = 1
    sim.cycle.step_time_mode = "time"
    feed = sim.add_step("Feed", time=duration * u.s)
    y_feed = ({"CarbonDioxide": 0.20, "Methane": 0.05, "Hydrogen": 0.75, "Helium": 0.0} if ternary
              else {"CarbonDioxide": 0.5, "Methane": 0.5, "Helium": 0.0})
    feed.bottom.set_streams("inlet", [sk.Stream("Feed", flow=sk.volumetric_flow(Q * p.flow_factor * u.m**3 / u.s),
                                                y=y_feed, temperature=T * u.K)])
    feed.top.set_pressure("outlet", P * u.Pa, y={c: (1.0 if c == "Helium" else 0.0) for c in comps},
                          temperature=T * u.K)

    s = sim.settings
    s.init_step, s.min_step, s.max_step = 1e-10 * u.s, 1e-15 * u.s, 0 * u.s
    s.reltol, s.nonlin_conv_coef, s.max_order = 1e-6, 0.1, 5
    s.max_nonlin_iters, s.max_conv_fails, s.max_err_test_fails = 15, 200, 40
    s.max_num_steps, s.max_iters_per_step = 100000, 50000
    s.advection_scheme = "weno"
    s.max_cycles_to_keep = 1
    s.max_log_points_per_step = 400
    sim.css.method = "none"
    sim.energy_balance.energy_balance_mode = "non_isothermal"
    return sim


# ---- metrics: the definitions of core/test/fixtures/validation_metrics.h -----------------

def load_csv(f: Path):
    rows = sorted(tuple(float(v) for v in line.split(",")[:2])
                  for line in f.read_text(encoding="utf-8-sig").splitlines()[1:] if line.strip())
    return np.array([r[0] for r in rows]), np.array([r[1] for r in rows])


def increasing(t, y):
    t, y = np.asarray(t), np.asarray(y)
    keep = np.concatenate([[True], np.diff(t) > 0])
    return t[keep], y[keep]


def crossing(t, y, level):
    idx = np.nonzero(y >= level)[0]
    if len(idx) == 0:
        return float("nan")
    i = idx[0]
    if i == 0:
        return float(t[0])
    return float(t[i - 1] + (level - y[i - 1]) / (y[i] - y[i - 1]) * (t[i] - t[i - 1]))


def score(sim, case):
    ternary = CASES[case][3]
    y_feed = ({"CarbonDioxide": 0.2, "Methane": 0.05, "Hydrogen": 0.75} if ternary
              else {"CarbonDioxide": 0.5, "Methane": 0.5})
    label = {"CarbonDioxide": "CO2", "Methane": "CH4", "Hydrogen": "H2"}
    metrics, residuals = {}, []
    for c, yf in y_feed.items():
        t, y = increasing(*(lambda s: (s.time, s.values))(sim.series("mole_fraction", c, location="top", cycle=0)))
        te, ye = load_csv(DATA / case / f"{case}_mole_fraction_{label[c]}.csv")
        m = (te >= t[0]) & (te <= t[-1])
        e = np.interp(te[m], t, y) - ye[m]
        metrics[f"nrmse_{label[c]}"] = float(np.sqrt(np.mean(e * e))) / yf
        if c != "Hydrogen":
            metrics[f"front_{label[c]}"] = crossing(t, y, yf / 2) / crossing(te, ye, yf / 2) - 1
        residuals.append(e / SIGMA_Y)
    sq, n = 0.0, 0
    for zc in Z_CM:
        s = sim.series("temperature", z=zc / 100, cycle=0)
        t, y = increasing(s.time, s.values)
        te, ye = load_csv(DATA / case / f"{case}_temperature_z={zc}[cm].csv")
        m = (te >= t[0]) & (te <= t[-1])
        e = np.interp(te[m], t, y) - ye[m]
        sq, n = sq + float(e @ e), n + int(m.sum())
        residuals.append(e / SIGMA_T)
    metrics["T_rmse"] = math.sqrt(sq / n)
    return metrics, np.concatenate(residuals)


def run(case, p):
    sim = build(case, p)
    r = sim.run()
    if not r.ok:
        raise RuntimeError(f"{case}: {r.message}")
    return sim, *score(sim, case)


def fmt(m):
    return "  ".join(f"{k} {v:+.4f}" if k.startswith("front") else f"{k} {v:.4f}" for k, v in m.items())


# ---- stage 1: flow factor from the thermal-front speed -------------------------------------

ZS = [35, 60, 85, 110]


def front_arrival(t, y, t0):
    return crossing(t, y, t0 + 0.5 * (y.max() - t0))


def front_slope(times):
    return np.polyfit(ZS, times, 1)[0]


def flow_factor(args):
    case, p = args
    exp = []
    for zc in ZS:
        te, ye = load_csv(DATA / case / f"{case}_temperature_z={zc}[cm].csv")
        exp.append(front_arrival(te, ye, ye[0]))
    f = 1.0
    for _ in range(4):
        sim, m, _ = run(case, replace(p, flow_factor=f))
        T0 = CASES[case][1]
        sim_t = [front_arrival(*increasing(*(lambda s: (s.time, s.values))(sim.series("temperature", z=zc / 100, cycle=0))), T0)
                 for zc in ZS]
        ratio = front_slope(sim_t) / front_slope(exp)
        if abs(ratio - 1) < 2e-3:
            break
        f *= ratio
    return case, f, m


# ---- stage 2: least squares on A3 -----------------------------------------------------------

NAMES = ["k_co2", "k_ch4", "h_out", "line_volume"]


def to_params(x):
    return Params(k_co2=math.exp(x[0]), k_ch4=math.exp(x[1]), h_out=math.exp(x[2]),
                  line_volume=x[3] * 1e-3, flow_factor=FLOW_FACTORS["A3"])


def residual(x):
    return run("A3", to_params(x))[2]


def fit(jobs):
    from scipy.optimize import least_squares
    x0 = np.array([math.log(0.084), math.log(0.39), math.log(108.0), 0.18])   # the first fit's values
    lo = [math.log(0.01), math.log(0.02), math.log(10.0), 0.0]
    hi = [math.log(1.0), math.log(5.0), math.log(2000.0), 0.6]
    steps = np.array([0.05, 0.05, 0.05, 0.01])
    with ProcessPoolExecutor(jobs) as pool:
        def jac(x):
            xs = [x] + [x + np.eye(len(x))[i] * steps[i] for i in range(len(x))]
            r = list(pool.map(residual, xs))
            return np.column_stack([(r[i + 1] - r[0]) / steps[i] for i in range(len(x))])
        sol = least_squares(residual, x0, jac=jac, bounds=(lo, hi), max_nfev=30, xtol=1e-3, ftol=1e-4)
    cov = np.linalg.pinv(sol.jac.T @ sol.jac) * (2 * sol.cost / max(1, len(sol.fun) - len(sol.x)))
    sd = np.sqrt(np.diag(cov))
    p = to_params(sol.x)
    print(json.dumps({"params": {n: getattr(p, n) for n in NAMES}, "sd (log for k, h; L for the line)": dict(zip(NAMES, sd.tolist())),
                      "cost": sol.cost}, indent=1))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("stage", choices=["flows", "fit", "predict"])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--clip", type=float, default=0.0, help="hold flow factors to +-clip of the set point")
    args = ap.parse_args()
    p = Params()
    if args.stage == "flows":
        with ProcessPoolExecutor(args.jobs) as pool:
            for case, f, m in pool.map(flow_factor, [(c, p) for c in CASES]):
                print(f"{case} flow factor {f:.4f}  {fmt(m)}", flush=True)
    elif args.stage == "fit":
        fit(args.jobs)
    else:
        def clipped(f):
            return min(max(f, 1 - args.clip), 1 + args.clip) if args.clip > 0 else f
        with ProcessPoolExecutor(args.jobs) as pool:
            futures = {c: pool.submit(run, c, replace(p, flow_factor=clipped(FLOW_FACTORS[c]))) for c in CASES}
            for c, fut in futures.items():
                print(f"{c} flow factor {clipped(FLOW_FACTORS[c]):.3f}  {fmt(fut.result()[1])}", flush=True)


if __name__ == "__main__":
    main()
