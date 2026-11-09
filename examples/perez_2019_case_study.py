"""Perez 2019 case study: the base-case decks, written the way a plant runs them.

The validated Perez 2019 decks (examples/perez_2019.py) impose the paper's
step times and pressure profiles directly. The case-study decks keep the
paper's column, adsorbent, feed and operating points, but drive every step by
a device and time the cycle as an n-bed schedule:

* **Schedule (fraction mode).** ``cycle.total_time`` and a fraction per step,
  written as expressions of the user variable ``n_beds``. The feed is
  continuous (the feed steps take 1/n of the cycle) and so is the vacuum pump
  (one bed is always in blowdown or evacuation: together 1/n of the cycle).
  The light-product pressurisation (LPP) step is as long as the first feed
  step of the partner bed that supplies it.
* **Devices.** The feed is delivered by a positive-displacement blower (a set
  flow with the compressor power calculated against the bed and vent-valve
  pressure drop); the light product leaves through a valve to the atmosphere;
  the four-step cycle pressurises from the feed header through a valve, the
  LPP cycle through a valve from the partner bed; blowdown and evacuation are
  driven by one vacuum pump of constant shaft power with the pressure-
  dependent vacuum-pump efficiency curve. The Ergun momentum balance is on.
* **Constraints as parameters.** Every device setting is a user variable, so
  each practical constraint is one registry path (``variable[...]``).

The pump power, the blowdown share of the pump's time and the valve
coefficients were chosen so that the device-driven steps reach the pressures
the validated decks reach at cyclic steady state (``REFERENCE`` below); the
docs page (docs/case-study-vsa.md) gives the comparison.

    python examples/perez_2019_case_study.py                # both decks, to CSS
    python examples/perez_2019_case_study.py --case Base_LPP
    python examples/perez_2019_case_study.py --save-decks   # rewrite examples/decks/

``build("Base")`` returns the configured simulation without running it
(python/tests/test_examples_build.py checks it against examples/decks/).

Reference: L. Estupinan Perez, P. Sarkar, A. Rajendran, Separation and
Purification Technology 224 (2019) 553-563. doi:10.1016/j.seppur.2019.05.039
"""

from __future__ import annotations

import sys

import numpy as np

import _common as c
import perez_2019 as paper
import skarstrom as sk
from skarstrom import ureg as u

FEED, BACKFLOW, NO_GAS = paper.FEED, paper.BACKFLOW, paper.NO_GAS

# The validated configuration each base case starts from, its feed (kept as in
# the paper) and the device settings chosen to reach its operating point.
CASES = {
    # Four-step cycle of configuration A: feed, blowdown, evacuation, pressurisation.
    "Base": dict(lpp=False, configuration="A", feed_time=99 * u.s, feed_velocity=0.08486999999999999 * u.m / u.s,
                 vacuum_pump_power=10.4 * u.W, blowdown_share=0.197, cv_vent=0.5, cv_pressurisation=0.006),
    # LPP cycle of configuration D1: feed 1 (supplies the LPP), feed 2, blowdown, evacuation, LPP.
    "Base_LPP": dict(lpp=True, configuration="D1", feed_times=(27.5 * u.s, 117 * u.s),
                     feed_velocity=0.032799999999999996 * u.m / u.s,
                     vacuum_pump_power=5.8 * u.W, blowdown_share=0.237, cv_vent=0.5, cv_lpp=0.008),
}
DECKS = {"Base": "Perez2019_CaseStudy_Base", "Base_LPP": "Perez2019_CaseStudy_Base_LPP"}

# What the validated decks reach at cyclic steady state (this engine, 80 cells; vacuum-pump
# energy with the same pump curve as the case-study decks): the operating points matched.
REFERENCE = {
    "Base": dict(deck="Perez2019_A", cycle_time=205.0, blowdown=(19.0, 7205.4), evacuation=(57.0, 3995.5),
                 pressurisation=(30.0, 99980.4), purity=0.99551, recovery=0.60144, productivity=3.678e-05,
                 energy=0.7517e6),
    "Base_LPP": dict(deck="Perez2019_D1", cycle_time=360.0, blowdown=(21.0, 8161.2), evacuation=(167.0, 3114.1),
                     pressurisation=(27.5, 101301.5), purity=0.94084, recovery=0.89480, productivity=1.559e-05,
                     energy=0.8211e6),
}

# Registry paths of the settings the case study varies (all user variables).
CONSTRAINT_PATHS = (
    "variable[n_beds]", "variable[t_feed]", "variable[feed_velocity]", "variable[blowdown_share]",
    "variable[vacuum_pump_power]", "variable[vp_eta_max]", "variable[vp_alpha]", "variable[blower_efficiency]",
    "variable[cv_vent]", "variable[p_atm]", "variable[p_evac_floor]", "variable[purity_target]",
    "variable[recovery_target]",
)


def build(case: str) -> sk.Simulation:
    """The base-case simulation ("Base" or "Base_LPP"), ready to run."""
    cfg = CASES[case]
    # Column, adsorbent, heat transfer, initial state, solver and CSS settings of the
    # validated configuration; its imposed-profile cycle is replaced below.
    sim = paper.build(cfg["configuration"])
    for step in reversed(sim.steps.names):
        sim.edit.remove_step(step)
    sim.info.title = f"Perez 2019 case study - base case{' with LPP' if cfg['lpp'] else ''}"
    sim.info.description = description(case)
    sim.layers[0].set_momentum("ergun")   # the blower and pump see the full Ergun pressure drop

    # ---- the settings the case study varies: one registry path each ---------------------
    sim.add_variable("n_beds", 3, description="beds in the schedule: feed 1/n, vacuum pump 1/n of the cycle")
    sim.add_variable("blowdown_share", cfg["blowdown_share"],
                     description="share of the vacuum pump's 1/n of the cycle spent in blowdown")
    sim.add_variable("feed_velocity", cfg["feed_velocity"], category="velocity",
                     description="superficial feed velocity delivered by the blower")
    sim.add_variable("p_atm", 1 * u.bar, category="pressure",
                     description="atmosphere: vent, blower suction and vacuum-pump discharge pressure")
    sim.add_variable("blower_efficiency", 0.7, description="feed blower efficiency")
    sim.add_variable("vacuum_pump_power", cfg["vacuum_pump_power"], category="custom",
                     description="vacuum pump shaft power at full vacuum [W] (W = W_max (1 - P / p_atm))")
    sim.add_variable("vp_eta_max", 0.8, description="vacuum pump efficiency curve: maximum efficiency")
    sim.add_variable("vp_alpha", 19.55, description="vacuum pump efficiency curve: alpha [1/bar]")
    sim.add_variable("cv_vent", cfg["cv_vent"], category="custom",
                     description="light-product vent valve coefficient [mol/(s bar)]")
    # Constraint values the case study checks (they do not change the run).
    sim.add_variable("p_evac_floor", 10 * u.kPa, category="pressure",
                     description="constraint: lowest allowed evacuation pressure")
    sim.add_variable("purity_target", 0.95, description="constraint: CO2 purity target")
    sim.add_variable("recovery_target", 0.90, description="constraint: CO2 recovery target")

    sim.cycle.step_time_mode = "fraction"
    if cfg["lpp"]:
        add_lpp_cycle(sim, cfg)
    else:
        add_basic_cycle(sim, cfg)
    return sim


def blower() -> sk.Pump:
    """A positive-displacement feed blower: it delivers the set flow; the power follows from the
    pressure it works against (bed inlet pressure over the suction pressure)."""
    return sk.Pump("compressor", "isothermal", "constant", efficiency="blower_efficiency", suction_pressure="p_atm",
                   polytropic_index=1.4)


def vacuum_pump() -> sk.Pump:
    """The vacuum pump, discharging to the atmosphere, with the pressure-dependent efficiency
    eta = eta_max alpha p / (1 + alpha p), p the suction pressure in bar."""
    return sk.Pump("vacuum_pump", "isothermal", "vacuum_pump", eta_max="vp_eta_max", alpha="vp_alpha",
                   discharge_pressure="p_atm", polytropic_index=1.4)


def vacuum_pump_power() -> sk.PolynomialPower:
    """The vacuum pump's shaft power against its suction pressure P: W = W_max (1 - P / p_atm).

    The pump unloads as the pressure ratio goes to 1, so near the atmosphere it moves a
    finite volume flow (eta W_max / p_atm per unit of the relative pressure rise, like a
    constant-speed pump) instead of the unbounded flow a constant power gives there; the
    power rises to W_max at full vacuum.
    """
    return sk.polynomial_power(a="vacuum_pump_power", b="-vacuum_pump_power / p_atm")


def vent(name: str) -> sk.ValveStream:
    return sk.ValveStream(name, pressure="p_atm", cv="cv_vent", xt=0.53, y=NO_GAS, temperature=293 * u.K)


def feed_inlet(name: str) -> sk.Stream:
    return sk.Stream(name, flow=sk.velocity("feed_velocity"), y=FEED, temperature=298 * u.K, pump=blower())


def add_vacuum_steps(sim: sk.Simulation):
    """Blowdown (top) then evacuation (bottom), both on the one vacuum pump, together 1/n of the cycle."""
    blowdown = sim.add_step("Blowdown", fraction="blowdown_share / n_beds")
    blowdown.bottom.set_closed()
    blowdown.top.set_streams("outlet", [
        sk.PumpStream("Blowdown pump", power=vacuum_pump_power(), pump=vacuum_pump(), y=BACKFLOW,
                      temperature=293 * u.K),
    ])
    evacuation = sim.add_step("Evacuation", fraction="(1 - blowdown_share) / n_beds")
    evacuation.bottom.set_streams("outlet", [
        sk.PumpStream("Evacuation pump", power=vacuum_pump_power(), pump=vacuum_pump(), y=BACKFLOW,
                      temperature=293 * u.K),
    ])
    evacuation.top.set_closed()
    return blowdown, evacuation


def add_basic_cycle(sim: sk.Simulation, cfg: dict) -> None:
    """Feed 1/n, blowdown + evacuation 1/n, pressurisation from the feed header for the rest."""
    # The feed time is kept: the cycle is n times the feed step.
    sim.add_variable("t_feed", cfg["feed_time"], category="time", description="feed step time (1/n of the cycle)")
    sim.cycle.total_time = "n_beds * t_feed"
    sim.add_variable("cv_pressurisation", cfg["cv_pressurisation"], category="custom",
                     description="pressurisation valve coefficient [mol/(s bar)], from the feed header")
    sim.add_variable("p_feed_header", 1 * u.bar, category="pressure",
                     description="feed header pressure, the source of the pressurisation gas")

    feed = sim.add_step("Feed", fraction="1 / n_beds")
    feed.bottom.set_streams("inlet", [feed_inlet("Feed inlet")])
    feed.top.set_streams("outlet", [vent("Vent")])

    _, evacuation = add_vacuum_steps(sim)

    pressurisation = sim.add_step("Pressurisation", fraction="1 - 2 / n_beds")
    pressurisation.bottom.set_streams("inlet", [
        sk.ValveStream("Feed header", pressure="p_feed_header", cv="cv_pressurisation", xt=0.53, y=FEED,
                       temperature=298 * u.K),
    ])
    pressurisation.top.set_closed()

    sim.kpi.track("CarbonDioxide", feed=[feed.bottom["Feed inlet"], pressurisation.bottom],
                  product=[evacuation.bottom])


def add_lpp_cycle(sim: sk.Simulation, cfg: dict) -> None:
    """Feed 1 + feed 2 1/n, blowdown + evacuation 1/n, LPP during the partner's feed 1, then idle."""
    feed_1_time, feed_2_time = (t.to("s").magnitude for t in cfg["feed_times"])
    sim.add_variable("t_feed", (feed_1_time + feed_2_time) * u.s, category="time",
                     description="feed time, feed 1 + feed 2 (1/n of the cycle)")
    sim.cycle.total_time = "n_beds * t_feed"
    sim.add_variable("feed1_share", round(feed_1_time / (feed_1_time + feed_2_time), 6),
                     description="share of the feed 1/n in feed 1 (the step that supplies the LPP)")
    sim.add_variable("cv_lpp", cfg["cv_lpp"], category="custom",
                     description="light-product pressurisation valve coefficient [mol/(s bar)]")

    feed_1 = sim.add_step("Feed 1", fraction="feed1_share / n_beds")
    feed_1.bottom.set_streams("inlet", [feed_inlet("Feed 1 inlet")])
    feed_1.top.set_streams("outlet", [
        vent("Feed 1 vent"),
        # The LPP valve between this bed and the bed being pressurised; its driving
        # pressure is that bed's pressure (through the connection).
        sk.ValveStream("LPP out", pressure="p_atm", cv="cv_lpp", xt=0.53, y=NO_GAS, temperature=293 * u.K),
    ])
    feed_2 = sim.add_step("Feed 2", fraction="(1 - feed1_share) / n_beds")
    feed_2.bottom.set_streams("inlet", [feed_inlet("Feed 2 inlet")])
    feed_2.top.set_streams("outlet", [vent("Feed 2 vent")])

    _, evacuation = add_vacuum_steps(sim)

    pressurisation = sim.add_step("Pressurisation", fraction="feed1_share / n_beds")
    pressurisation.bottom.set_closed()
    pressurisation.top.set_streams("inlet", [
        sk.ValveStream("LPP in", pressure="p_atm", cv="cv_lpp", xt=0.53, y=BACKFLOW, temperature=298 * u.K),
    ])
    # The rest of the bed's third of the cycle: closed, waiting for its feed slot.
    sim.add_step("Idle", fraction="1 - (2 + feed1_share) / n_beds")

    # The two ends of the LPP valve: each reads the other bed's pressure.
    sim.connect(feed_1.top["LPP out"], pressurisation.top["LPP in"])
    sim.kpi.track("CarbonDioxide", feed=[feed_1.bottom["Feed 1 inlet"], feed_2.bottom["Feed 2 inlet"]],
                  product=[evacuation.bottom])


def description(case: str) -> str:
    cfg = CASES[case]
    cycle = ("feed (1/n), blowdown and evacuation on one vacuum pump (1/n), LPP through a valve from the "
             "partner bed in its feed 1, idle" if cfg["lpp"] else
             "feed (1/n), blowdown and evacuation on one vacuum pump (1/n), pressurisation from the feed "
             "header through a valve")
    return (f"Case-study base deck: configuration {cfg['configuration']} of Estupinan Perez et al. (2019), "
            f"CO2/N2 VSA on 13X, run by devices on a {'6' if cfg['lpp'] else '4'}-step n-bed schedule "
            f"(fraction mode, n_beds = 3): {cycle}. Feed by a blower, vent valve to the atmosphere. Ergun momentum "
            f"(the validated Perez2019 decks use Kozeny-Carman). The pump and valves are sized to reach the "
            f"step-end pressures the validated deck {REFERENCE[case]['deck']} reaches at cyclic steady state "
            f"(blowdown {REFERENCE[case]['blowdown'][1] / 1e3:.1f} kPa, evacuation "
            f"{REFERENCE[case]['evacuation'][1] / 1e3:.1f} kPa), not the nominal end values of its imposed "
            f"profiles, which it does not reach in the step time. Settings the case study varies: " + ", ".join(CONSTRAINT_PATHS)
            + (", variable[cv_lpp], variable[feed1_share]" if cfg["lpp"] else
               ", variable[cv_pressurisation], variable[p_feed_header]")
            + ". p_evac_floor, purity_target and recovery_target are constraint values, not used by the run.")


# ---- running and comparing ---------------------------------------------------------------

TOLERANCE = 0.02   # step-end pressures are matched to 2 %

def step_report(sim: sk.Simulation, name: str, target: float) -> tuple:
    """(step time, end pressure, time to come within TOLERANCE of target) of a step in the last cycle.

    A pressurisation counts as arrived once it is within TOLERANCE of the lower of the
    target and its own end pressure (it cannot rise above its source).
    """
    i = sim.steps.names.index(name)
    p = sim.profile("pressure", cycle=-1, step=i)
    mean = p.values.mean(axis=1)    # bed-averaged pressure
    t = p.time - p.time[0]
    rising = mean[-1] > mean[0]
    reached = (mean >= (1 - TOLERANCE) * min(target, mean[-1])) if rising else (mean <= (1 + TOLERANCE) * target)
    return float(t[-1]), float(mean[-1]), float(t[np.argmax(reached)]) if reached.any() else float("nan")


def report(sim: sk.Simulation, case: str) -> None:
    ref = REFERENCE[case]
    rows = [["cycle time [s]", f"{ref['cycle_time']:.1f}", f"{sim.cycle.total_time:.1f}", ""]]
    for label, step, key in (("blowdown", "Blowdown", "blowdown"), ("evacuation", "Evacuation", "evacuation"),
                             ("pressurisation", "Pressurisation", "pressurisation")):
        t_ref, p_ref = ref[key]
        t, p_end, t_hit = step_report(sim, step, p_ref)
        rows.append([f"{label}: time [s] / end pressure [kPa]", f"{t_ref:.1f} / {p_ref / 1e3:.2f}",
                     f"{t:.1f} / {p_end / 1e3:.2f}", f"{100 * (p_end / p_ref - 1):+.1f} %, within {100 * TOLERANCE:.0f} % after {t_hit:.1f} s"])
    k = c.kpi(sim, "CarbonDioxide")
    for label, key, value, scale, unit in (("CO2 purity", "purity", k.purity, 100, "%"),
                                           ("CO2 recovery", "recovery", k.recovery, 100, "%"),
                                           ("productivity", "productivity", k.productivity_mass, 1e5, "1e-5 kg/(kg s)"),
                                           ("specific energy", "energy", k.energy_consumption, 1e-6, "MJ/kg")):
        rows.append([f"{label} [{unit}]", f"{scale * ref[key]:.3g}", f"{scale * value:.3g}", ""])
    c.table(rows, ["quantity", ref["deck"], DECKS[case], "match"])


def main(argv=None) -> int:
    ap = c.parser(__doc__.splitlines()[0], list(CASES))
    ap.add_argument("--save-decks", action="store_true", help="write the decks to examples/decks/ and stop")
    args = ap.parse_args(argv)
    for case in c.selected(args):
        sim = build(case)
        if args.save_decks:
            path = c.DECKS / f"{DECKS[case]}.json"
            sim.save(str(path))
            print(f"saved {path}")
            continue
        if args.cycles is not None:
            sim.cycle.number_of_cycles = int(args.cycles)
        c.run(sim, f"case study {case}", args.quiet)
        print(f"  {c.css_status(sim)}")
        report(sim, case)
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
