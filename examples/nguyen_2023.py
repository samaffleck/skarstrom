"""Nguyen et al. (2023): CO2/N2 VPSA on structured CALF-20, experiments I-X.

Builds each experiment from scratch in code: the fluid, CALF-20 with a
competitive dual-site Langmuir isotherm, the column (sized to hold the 36 g of
adsorbent the paper reports), the cycle steps with their boundary conditions,
the KPI definition and the solver settings. Every value is set explicitly,
with its unit.

Experiments I-V use the basic 4-step cycle, VI-X the 4-step cycle with
light-product pressurisation. Each experiment runs to cyclic steady state; the
CO2 purity and recovery of the evacuation product (as the validation test
defines them) are compared with the measured values and with the paper's
detailed model (Table 2).

    python examples/nguyen_2023.py                 # all ten (several minutes)
    python examples/nguyen_2023.py --case I --case VI

``build("VI")`` returns the configured simulation without running it
(python/tests/test_examples_build.py checks it against examples/decks/).

Reference: T.T.T. Nguyen, G.K.H. Shimizu, A. Rajendran, Chemical Engineering
Journal 452 (2023) 139550. doi:10.1016/j.cej.2022.139550
"""

from __future__ import annotations

import math
import sys

import _common as c
import skarstrom as sk
from skarstrom import ureg as u  # pint units: 3.2 * u.cm, 40.98 * u.kJ / u.mol

# Operating conditions of the ten experiments (core/test/papers/nguyen_2023/
# nguyen_2023.cpp): total feed time, interstitial feed velocity, the pressure
# P_I that ends blowdown and, for the LPP cycle, the evacuation pressure.
EXPERIMENTS = {
    #       cycle         feed time [s]  interstitial velocity [m/s]   P_I [Pa]   evacuation [Pa]
    "I":    dict(lpp=False, feed_time=176, interstitial_velocity=0.093, blowdown_pressure=0.119e5),
    "II":   dict(lpp=False, feed_time=174, interstitial_velocity=0.065, blowdown_pressure=0.148e5),
    "III":  dict(lpp=False, feed_time=174, interstitial_velocity=0.058, blowdown_pressure=0.189e5),
    "IV":   dict(lpp=False, feed_time=174, interstitial_velocity=0.055, blowdown_pressure=0.233e5),
    "V":    dict(lpp=False, feed_time=174, interstitial_velocity=0.046, blowdown_pressure=0.263e5),
    "VI":   dict(lpp=True, feed_time=193, interstitial_velocity=0.041, blowdown_pressure=0.063e5,
                 evacuation_pressure=0.031e5),
    "VII":  dict(lpp=True, feed_time=157, interstitial_velocity=0.024, blowdown_pressure=0.046e5,
                 evacuation_pressure=0.031e5),
    "VIII": dict(lpp=True, feed_time=162, interstitial_velocity=0.023, blowdown_pressure=0.060e5,
                 evacuation_pressure=0.034e5),
    "IX":   dict(lpp=True, feed_time=156, interstitial_velocity=0.021, blowdown_pressure=0.066e5,
                 evacuation_pressure=0.031e5),
    "X":    dict(lpp=True, feed_time=156, interstitial_velocity=0.022, blowdown_pressure=0.080e5,
                 evacuation_pressure=0.031e5),
}

HIGH_PRESSURE = 1.0e5                                  # [Pa]
BED_VOID_FRACTION = 0.4
FEED = {"Nitrogen": 0.85, "CarbonDioxide": 0.15}       # mole fractions
BACKFLOW = {"CarbonDioxide": 0.0, "Nitrogen": 1.0}     # gas entering an outlet if the flow reverses
NO_GAS = {"CarbonDioxide": 0.0, "Nitrogen": 0.0}       # outlet streams: no backflow composition

# Column (Table 1). Table 1 is internally inconsistent (36.0 g at 570 kg/m3 and
# eps = 0.4 does not fit a 0.021 m x 0.267 m column), so the length is derived
# from the measured mass, keeping the adsorbent inventory: L ~ 0.304 m.
ADSORBENT_MASS = 36.0e-3        # [kg]
INNER_DIAMETER = 0.021          # [m]
OUTER_DIAMETER = 0.025          # [m]
PARTICLE_DENSITY = 570.0        # [kg/m3]
BED_LENGTH = ADSORBENT_MASS / (PARTICLE_DENSITY * (1 - BED_VOID_FRACTION)) / (math.pi / 4 * INNER_DIAMETER**2)
WALL_THICKNESS = (OUTER_DIAMETER - INNER_DIAMETER) / 2


def build(case: str) -> sk.Simulation:
    """The simulation of one experiment ("I" ... "X"), ready to run."""
    exp = EXPERIMENTS[case]
    cycle = "Four-step LPP cycle" if exp["lpp"] else "Basic four-step cycle"
    sim = sk.Simulation.new(
        title=f"Nguyen 2023 - Experiment {case}.",
        description=f"{cycle} for a binary CO2/N2 VPSA cycle using CALF-20 (Nguyen et al. 2023), experiment {case}.",
        author="Skarstrom examples",
    )

    # ---- fluid ---------------------------------------------------------------------
    sim.add_component("CarbonDioxide")
    sim.add_component("Nitrogen")
    sim.fluid.adsorbed_heat_capacity_mode = "constant"
    sim.fluid.heat_of_adsorption_convention = "isosteric"
    sim.fluid.adsorbed_heat_capacity["CarbonDioxide"] = 37.435 * u.J / (u.mol * u.K)
    sim.fluid.adsorbed_heat_capacity["Nitrogen"] = 29.171 * u.J / (u.mol * u.K)

    # ---- adsorbent: structured CALF-20, competitive dual-site Langmuir (SI Table S1) ----
    calf20 = sim.add_adsorbent("CALF-20", solid_density=PARTICLE_DENSITY * u.kg / u.m**3,
                               heat_capacity=1350 * u.J / (u.kg * u.K), thermal_conductivity=0.63 * u.W / (u.m * u.K),
                               particle_void_fraction=0.35)
    calf20.set_isotherm("extended_dual_site_langmuir_c", {
        "CarbonDioxide": {"m1": 2.687 * u.mol / u.kg, "b1": 5.742e-7 * u.m**3 / u.mol, "q1": 35.18 * u.kJ / u.mol,
                          "m2": 4.597 * u.mol / u.kg, "b2": 3.684e-9 * u.m**3 / u.mol, "q2": 35.03 * u.kJ / u.mol},
        "Nitrogen":      {"m1": 2.687 * u.mol / u.kg, "b1": 7.466e-7 * u.m**3 / u.mol, "q1": 17.91 * u.kJ / u.mol,
                          "m2": 4.597 * u.mol / u.kg, "b2": 7.466e-7 * u.m**3 / u.mol, "q2": 17.91 * u.kJ / u.mol},
    }, components=["CarbonDioxide", "Nitrogen"])
    for component in ("CarbonDioxide", "Nitrogen"):
        calf20.set_rate_model(component, "ldf", "macropore", tortuosity=3.0, macropore_void_fraction=0.35)
    calf20.heat_of_adsorption["CarbonDioxide"] = 38 * u.kJ / u.mol
    calf20.heat_of_adsorption["Nitrogen"] = 15 * u.kJ / u.mol

    # ---- column ------------------------------------------------------------------------
    sim.bed.set_shape("cylinder", diameter=INNER_DIAMETER * u.m)
    sim.bed.ambient_temperature = 293 * u.K                      # not stated in the paper: an assumption
    sim.bed.set_thermal_equilibrium("lte")
    sim.bed.set_particle_heat_transfer("constant", particle_htc=100 * u.W / (u.m**2 * u.K))
    for void in (sim.bed.bottom_void, sim.bed.top_void):       # no dead volume at the column ends
        void.enabled = False
        void.volume = 1e-4 * u.m**3
        void.resistance_coefficient = 0 * u.Pa * u.s / u.m
    sim.wall.set_material("ss316l", density=7800 * u.kg / u.m**3, heat_capacity=502 * u.J / (u.kg * u.K),
                          thermal_conductivity=10 * u.W / (u.m * u.K))
    sim.wall.thickness = WALL_THICKNESS * u.m
    sim.wall.set_boundary("zero_flux")

    # ~1 mm pellets (paper Section 2).
    layer = sim.add_layer("CALF-20", composition={"CALF-20": 1.0}, momentum="kozeny_carman",
                          length=BED_LENGTH * u.m, cells=100, bed_void_fraction=BED_VOID_FRACTION,
                          particle_diameter=1 * u.mm)
    layer.set_dispersion("ruthven", molecular_diffusion_coefficient=1.6e-5 * u.m**2 / u.s)
    layer.set_heat_dispersion("constant", heat_dispersion_coefficient=0.4 * u.W / (u.m * u.K))
    layer.set_internal_heat_transfer("constant", internal_htc=8 * u.W / (u.m**2 * u.K))
    layer.set_external_heat_transfer("constant", external_htc=10 * u.W / (u.m**2 * u.K))

    # ---- initial state: the bed full of N2 at 1 bar and 293 K --------------------------
    sim.initial_conditions.pressure = 1 * u.bar
    sim.initial_conditions.temperature = 293 * u.K
    sim.initial_conditions.wall_temperature = 293 * u.K
    sim.initial_conditions.velocity = 0 * u.m / u.s
    sim.initial_conditions.y["Nitrogen"] = 1.0
    sim.initial_conditions.y["CarbonDioxide"] = 0.0

    # ---- cycle -------------------------------------------------------------------------
    # The light-product pressurisation cycles need up to about 500 cycles to close the
    # CO2 balance (CSS below).
    sim.cycle.number_of_cycles = 1000 if exp["lpp"] else 250
    sim.cycle.step_time_mode = "time"
    if exp["lpp"]:
        add_lpp_cycle(sim, exp)
    else:
        add_basic_cycle(sim, exp)

    # ---- solver (SUNDIALS IDA), cyclic steady state (the paper's mass-balance
    # criterion, at least 100 cycles), energy balance ------------------------------------
    sim.settings.init_step = 1e-10 * u.s
    sim.settings.min_step = 1e-15 * u.s
    sim.settings.max_step = 2 * u.s
    sim.settings.reltol = 1e-6
    sim.settings.nonlin_conv_coef = 0.1
    sim.settings.max_order = 5
    sim.settings.max_nonlin_iters = 15
    sim.settings.max_conv_fails = 200
    sim.settings.max_err_test_fails = 40
    sim.settings.max_num_steps = 100000
    sim.settings.max_iters_per_step = 50000
    sim.settings.advection_scheme = "upwind"
    sim.settings.max_cycles_to_keep = 10
    sim.settings.max_log_points_per_step = 200
    sim.settings.acceleration_method = "none"
    sim.css.method = "mass_balance_component"
    sim.css.min_cycles = 100
    sim.css.molar_absolute_tolerance = 1e-4
    sim.css.molar_relative_tolerance = 1e-4
    sim.css.mass_balance_tolerance = 0.005
    sim.energy_balance.energy_balance_mode = "non_isothermal"
    sim.energy_balance.isothermal_use_initial_temperature = True
    sim.energy_balance.isothermal_temperature = 298.15 * u.K
    return sim


def feed_velocity(exp: dict) -> str:
    """Superficial feed velocity: the paper's interstitial velocity times the bed voidage."""
    return (exp['interstitial_velocity'] * BED_VOID_FRACTION) * u.m / u.s


def add_basic_cycle(sim: sk.Simulation, exp: dict) -> None:
    """Feed, blowdown (top) to P_I, evacuation (bottom), pressurisation with feed gas (bottom)."""
    # The paper ends blowdown with a pressure trigger at P_I on an exponential decay
    # P = P_H exp(-lambda_BD t), lambda_BD = 0.153 1/s (Table 1), so it lasts
    # t = ln(P_H / P_I) / lambda_BD.
    blowdown_time = math.log(HIGH_PRESSURE / exp["blowdown_pressure"]) / 0.153

    feed = sim.add_step("Feed", time=exp['feed_time'] * u.s)
    feed.bottom.set_streams("inlet", [
        sk.Stream("Feed inlet", flow=sk.velocity(feed_velocity(exp)), y=FEED, temperature=293 * u.K),
    ])
    feed.top.set_pressure("outlet", HIGH_PRESSURE * u.Pa, y=BACKFLOW, temperature=293 * u.K)

    blowdown = sim.add_step("Blowdown", time=blowdown_time * u.s)
    blowdown.bottom.set_closed()
    blowdown.top.set_pressure("outlet", sk.exponential(exp['blowdown_pressure'] * u.Pa), y=BACKFLOW,
                              temperature=293 * u.K)

    evacuation = sim.add_step("Evacuation", time=70 * u.s)
    evacuation.bottom.set_pressure("outlet", sk.exponential(3100 * u.Pa), y=BACKFLOW, temperature=293 * u.K)
    evacuation.top.set_closed()

    pressurisation = sim.add_step("Pressurisation", time=29 * u.s)
    pressurisation.bottom.set_pressure("inlet", sk.linear(HIGH_PRESSURE * u.Pa), y=FEED, temperature=293 * u.K)
    pressurisation.top.set_closed()

    # CO2 is fed in the feed and pressurisation steps and recovered by evacuation.
    sim.kpi.track("CarbonDioxide", feed=[feed.bottom["Feed inlet"], pressurisation.bottom],
                  product=[evacuation.bottom])


def add_lpp_cycle(sim: sk.Simulation, exp: dict) -> None:
    """Feed split in two steps; the first one's light product pressurises the bed from the top."""
    pressurisation_time = 50.0                       # [s]; feed 1 runs alongside pressurisation
    feed_1_time = pressurisation_time
    feed_2_time = exp["feed_time"] - feed_1_time
    # Light product leaves through a valve (flow coefficient Cv, pressure-drop ratio factor xT) to 1 bar.
    ambient = dict(pressure=1 * u.bar, cv=10.0, xt=0.53, y=NO_GAS, temperature=293 * u.K)

    feed_1 = sim.add_step("Feed 1", time=feed_1_time * u.s)
    feed_1.bottom.set_streams("inlet", [
        sk.Stream("Feed 1 inlet", flow=sk.velocity(feed_velocity(exp)), y=FEED, temperature=293 * u.K),
    ])
    feed_1.top.set_streams("outlet", [
        sk.ValveStream("Feed 1 ambient Stream", **ambient),           # light product to the atmosphere
        # The LPP draw: its flow is whatever the connected pressurisation step takes.
        sk.Stream("LPP Stream", flow=sk.velocity(0 * u.m / u.s), y=NO_GAS, temperature=293 * u.K),
    ])

    feed_2 = sim.add_step("Feed 2", time=feed_2_time * u.s)
    feed_2.bottom.set_streams("inlet", [
        sk.Stream("Feed 2 inlet", flow=sk.velocity(feed_velocity(exp)), y=FEED, temperature=293 * u.K),
    ])
    feed_2.top.set_streams("outlet", [sk.ValveStream("Feed 2 ambient Stream", **ambient)])

    blowdown = sim.add_step("Blowdown", time=60 * u.s)
    blowdown.bottom.set_closed()
    blowdown.top.set_pressure("outlet", sk.exponential(exp['blowdown_pressure'] * u.Pa), y=BACKFLOW,
                              temperature=293 * u.K)

    evacuation = sim.add_step("Evacuation", time=40 * u.s)
    evacuation.bottom.set_pressure("outlet", sk.exponential(exp['evacuation_pressure'] * u.Pa), y=BACKFLOW,
                                   temperature=293 * u.K)
    evacuation.top.set_closed()

    pressurisation = sim.add_step("Pressurisation", time=pressurisation_time * u.s)
    pressurisation.bottom.set_closed()
    pressurisation.top.set_pressure("inlet", sk.linear(HIGH_PRESSURE * u.Pa), y=BACKFLOW, temperature=293 * u.K)

    # Light-product pressurisation: the LPP stream of Feed 1 supplies the pressurisation inlet.
    sim.connect(feed_1.top["LPP Stream"], pressurisation.top)

    # CO2 is fed in both feed steps and recovered by evacuation. The LPP stream is N2-rich
    # gas moving within the cycle, not CO2 product, so it is not a KPI endpoint.
    sim.kpi.track("CarbonDioxide", feed=[feed_1.bottom["Feed 1 inlet"], feed_2.bottom["Feed 2 inlet"]],
                  product=[evacuation.bottom])


NUMERALS = list(EXPERIMENTS)
CASES = {n: f"Experiment{n}" for n in NUMERALS}                  # case -> reference_values.csv case
DECKS = {n: f"Nguyen2023_Exp_{n}" for n in NUMERALS}             # the shipped deck of each case
QUANTITIES = [
    ("CO2 purity vs experiment", "CarbonDioxide.purity", "CO2_purity_exp"),
    ("CO2 recovery vs experiment", "CarbonDioxide.recovery", "CO2_recovery_exp"),
    ("CO2 purity vs paper model", "CarbonDioxide.purity", "CO2_purity_model"),
    ("CO2 recovery vs paper model", "CarbonDioxide.recovery", "CO2_recovery_model"),
]


def main(argv=None) -> int:
    return c.cyclic_family(argv, __doc__, "nguyen_2023", CASES, QUANTITIES, "CarbonDioxide", "nguyen_2023.png",
                           values=c.evacuation_product_kpis, build=build)


if __name__ == "__main__":
    sys.exit(main())
