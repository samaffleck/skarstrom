"""Estupinan Perez et al. (2019): CO2/N2 VSA on zeolite 13X, configurations A-F.

Builds each configuration from scratch in code: the fluid, the adsorbent and
its isotherm and rate model, the column and its heat transfer, the cycle
steps with their boundary conditions, the KPI definition and the solver
settings. Every value is set explicitly, with its unit, so this file is also a
template for a VSA cycle of your own.

A, B and C1 use the basic four-step cycle (feed, blowdown, evacuation,
pressurisation with feed gas). D1, E and F split the feed step and use part of
the light product of the first feed step to pressurise the bed from the top
(light-product pressurisation, LPP). Each configuration runs to cyclic steady
state; the CO2 purity and recovery of the evacuation product (as the
validation test defines them) are compared with the measured values and the
paper's detailed model (Table 3).

    python examples/perez_2019.py                  # all six
    python examples/perez_2019.py --case A
    python examples/perez_2019.py --case E --cycles 20

``build("A")`` returns the configured simulation without running it
(python/tests/test_examples_build.py checks it against examples/decks/).

Reference: L. Estupinan Perez, P. Sarkar, A. Rajendran, Separation and
Purification Technology 224 (2019) 553-563. doi:10.1016/j.seppur.2019.05.039
"""

from __future__ import annotations

import sys

import _common as c
import skarstrom as sk
from skarstrom import ureg as u  # pint units: 3.2 * u.cm, 40.98 * u.kJ / u.mol

# Operating conditions of the six configurations. Pressures are the end
# values of the blowdown and evacuation ramps; velocities are superficial at
# the feed end.
CONFIGURATIONS = {
    #       cycle   feed step(s)        feed velocity                    blowdown          evacuation         pressurisation
    "A":  dict(lpp=False, feed_time=99 * u.s, feed_velocity=0.08486999999999999 * u.m / u.s,
               blowdown_time=19 * u.s, blowdown_pressure=6000 * u.Pa, evacuation_time=57 * u.s,
               evacuation_pressure=3200 * u.Pa, pressurisation_time=30 * u.s, high_pressure=1 * u.bar),
    "B":  dict(lpp=False, feed_time=99 * u.s, feed_velocity=0.05985999999999999 * u.m / u.s,
               blowdown_time=19 * u.s, blowdown_pressure=7200 * u.Pa, evacuation_time=57 * u.s,
               evacuation_pressure=3200 * u.Pa, pressurisation_time=30 * u.s, high_pressure=1 * u.bar),
    "C1": dict(lpp=False, feed_time=114 * u.s, feed_velocity=0.03198 * u.m / u.s,
               blowdown_time=19 * u.s, blowdown_pressure=10700 * u.Pa, evacuation_time=57 * u.s,
               evacuation_pressure=3200 * u.Pa, pressurisation_time=30 * u.s, high_pressure=1 * u.bar),
    "D1": dict(lpp=True, feed_times=(27.5 * u.s, 117 * u.s), feed_velocity=0.032799999999999996 * u.m / u.s,
               blowdown_time=21 * u.s, blowdown_pressure=7200 * u.Pa, evacuation_time=167 * u.s,
               evacuation_pressure=3000 * u.Pa, pressurisation_time=27.5 * u.s, high_pressure=1 * u.atm),
    "E":  dict(lpp=True, feed_times=(27.5 * u.s, 112.5 * u.s), feed_velocity=0.024189999999999996 * u.m / u.s,
               blowdown_time=21 * u.s, blowdown_pressure=9100 * u.Pa, evacuation_time=167 * u.s,
               evacuation_pressure=3100 * u.Pa, pressurisation_time=27.5 * u.s, high_pressure=1 * u.bar),
    "F":  dict(lpp=True, feed_times=(27.5 * u.s, 89.5 * u.s), feed_velocity=0.025009999999999998 * u.m / u.s,
               blowdown_time=21 * u.s, blowdown_pressure=13400 * u.Pa, evacuation_time=167 * u.s,
               evacuation_pressure=3200 * u.Pa, pressurisation_time=27.5 * u.s, high_pressure=1 * u.bar),
}

FEED = {"CarbonDioxide": 0.15, "Nitrogen": 0.85}       # flue gas, mole fractions
BACKFLOW = {"CarbonDioxide": 0.0, "Nitrogen": 1.0}     # gas entering an outlet if the flow reverses
NO_GAS = {"CarbonDioxide": 0.0, "Nitrogen": 0.0}       # outlet streams: no backflow composition


def add_zeolite_13x(sim: sk.Simulation) -> sk.Adsorbent:
    """Zeolite 13X with the paper's competitive dual-site Langmuir isotherm (Table 2)."""
    ads = sim.add_adsorbent(
        "13X",
        solid_density=1138.0252816422599 * u.kg / u.m**3,      # particle density
        heat_capacity=900 * u.J / (u.kg * u.K),
        thermal_conductivity=0.63 * u.W / (u.m * u.K),
        particle_void_fraction=0.35,
    )
    # Extended dual-site Langmuir in concentration form: b_i(T) = b0_i exp(q_i / RT).
    ads.set_isotherm("extended_dual_site_langmuir_c", {
        "CarbonDioxide": {"m1": 2.768 * u.mol / u.kg, "b1": 3.934e-7 * u.m**3 / u.mol, "q1": 40.98 * u.kJ / u.mol,
                          "m2": 2.758 * u.mol / u.kg, "b2": 1.256e-7 * u.m**3 / u.mol, "q2": 34.19 * u.kJ / u.mol},
        "Nitrogen":      {"m1": 2.768 * u.mol / u.kg, "b1": 3.134e-6 * u.m**3 / u.mol, "q1": 15.71 * u.kJ / u.mol,
                          "m2": 2.758 * u.mol / u.kg, "b2": 3.134e-6 * u.m**3 / u.mol, "q2": 15.72 * u.kJ / u.mol},
    }, components=["CarbonDioxide", "Nitrogen"])
    # Linear driving force with the macropore coefficient times the capacity factor
    # c / (rho_p q*) (paper SI Eq. S4).
    for component in ("CarbonDioxide", "Nitrogen"):
        ads.set_rate_model(component, "ldf", "macropore_capacity", tortuosity=3.0, macropore_void_fraction=0.35)
    ads.heat_of_adsorption["CarbonDioxide"] = 30.85 * u.kJ / u.mol
    ads.heat_of_adsorption["Nitrogen"] = 15.93 * u.kJ / u.mol
    return ads


def build(case: str) -> sk.Simulation:
    """The simulation of one configuration ("A", "B", "C1", "D1", "E" or "F"), ready to run."""
    cfg = CONFIGURATIONS[case]
    sim = sk.Simulation.new(
        title=f"Perez 2019 - Configuration {case}",
        description=f"Configuration {case} of Estupinan Perez et al. (2019), CO2/N2 VSA on 13X, "
                    f"{'with light-product pressurisation' if cfg['lpp'] else 'basic four-step cycle'}.",
        author="Skarstrom examples",
    )

    # ---- fluid ---------------------------------------------------------------------
    sim.add_component("CarbonDioxide")
    sim.add_component("Nitrogen")
    sim.fluid.adsorbed_heat_capacity_mode = "constant"
    sim.fluid.heat_of_adsorption_convention = "isosteric"
    sim.fluid.adsorbed_heat_capacity["CarbonDioxide"] = 30.7 * u.J / (u.mol * u.K)
    sim.fluid.adsorbed_heat_capacity["Nitrogen"] = 30.7 * u.J / (u.mol * u.K)

    # ---- adsorbent -------------------------------------------------------------------
    add_zeolite_13x(sim)

    # ---- column: 0.30 m x 3.2 cm stainless-steel tube (Table 1) ------------------------
    sim.bed.set_shape("cylinder", diameter=3.2 * u.cm)
    sim.bed.ambient_temperature = 293 * u.K
    sim.bed.set_thermal_equilibrium("lte")                     # gas and solid at one temperature
    sim.bed.set_particle_heat_transfer("constant", particle_htc=100 * u.W / (u.m**2 * u.K))
    for void in (sim.bed.bottom_void, sim.bed.top_void):       # no dead volume at the column ends
        void.enabled = False
        void.volume = 1e-4 * u.m**3
        void.resistance_coefficient = 0 * u.Pa * u.s / u.m
    sim.wall.set_material("ss316l", density=7800 * u.kg / u.m**3, heat_capacity=502 * u.J / (u.kg * u.K),
                          thermal_conductivity=16 * u.W / (u.m * u.K))
    sim.wall.thickness = 3.1 * u.mm
    sim.wall.set_boundary("fixed", fixed_temperature=293 * u.K)     # SI Table S2: wall ends at T_a

    layer = sim.add_layer("13X", composition={"13X": 1.0}, momentum="kozeny_carman",
                          length=0.3 * u.m, cells=30, bed_void_fraction=0.41, particle_diameter=1.5 * u.mm)
    layer.set_dispersion("ruthven", molecular_diffusion_coefficient=1.6e-5 * u.m**2 / u.s)
    layer.set_heat_dispersion("constant", heat_dispersion_coefficient=0.0903 * u.W / (u.m * u.K))
    layer.set_internal_heat_transfer("constant", internal_htc=8.6 * u.W / (u.m**2 * u.K))
    layer.set_external_heat_transfer("constant", external_htc=10 * u.W / (u.m**2 * u.K))

    # ---- initial state: the bed full of N2 at 1 bar and 298 K --------------------------
    sim.initial_conditions.pressure = 1 * u.bar
    sim.initial_conditions.temperature = 298 * u.K
    sim.initial_conditions.wall_temperature = 298 * u.K
    sim.initial_conditions.velocity = 0 * u.m / u.s
    sim.initial_conditions.y["Nitrogen"] = 1.0
    sim.initial_conditions.y["CarbonDioxide"] = 0.0

    # ---- cycle -------------------------------------------------------------------------
    # E and F need about 450 and 1000 cycles to close the CO2 balance (CSS below).
    sim.cycle.number_of_cycles = 1500 if case in ("E", "F") else 250
    sim.cycle.step_time_mode = "time"
    if cfg["lpp"]:
        add_lpp_cycle(sim, cfg)
    else:
        add_basic_cycle(sim, cfg)

    # ---- solver (SUNDIALS IDA), cyclic steady state, energy balance ----------------------
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
    sim.css.min_cycles = 5
    sim.css.molar_absolute_tolerance = 1e-4
    sim.css.molar_relative_tolerance = 1e-4
    sim.css.mass_balance_tolerance = 0.005
    sim.energy_balance.energy_balance_mode = "non_isothermal"
    sim.energy_balance.isothermal_use_initial_temperature = True
    sim.energy_balance.isothermal_temperature = 298.15 * u.K
    return sim


def add_basic_cycle(sim: sk.Simulation, cfg: dict) -> None:
    """Feed, blowdown (top), evacuation (bottom), pressurisation with feed gas (bottom)."""
    feed = sim.add_step("Feed", time=cfg["feed_time"])
    feed.bottom.set_streams("inlet", [
        sk.Stream("Feed inlet", flow=sk.velocity(cfg["feed_velocity"]), y=FEED, temperature=298 * u.K),
    ])
    feed.top.set_pressure("outlet", cfg["high_pressure"], y=BACKFLOW, temperature=293 * u.K)

    blowdown = sim.add_step("Blowdown", time=cfg["blowdown_time"])
    blowdown.bottom.set_closed()
    blowdown.top.set_pressure("outlet", sk.exponential(cfg["blowdown_pressure"]), y=BACKFLOW, temperature=293 * u.K)

    evacuation = sim.add_step("Evacuation", time=cfg["evacuation_time"])
    # Sect. 4.3: two-term exponential decay (0.019 and 0.101 1/s, Table 1); the split is not
    # printed and an equal weight is used (core/test/papers/perez_2019/README.md).
    evacuation.bottom.set_pressure("outlet", sk.double_exponential(cfg["evacuation_pressure"], rate_1=0.019, rate_2=0.101, weight=0.5),
                                   y=BACKFLOW, temperature=293 * u.K)
    evacuation.top.set_closed()

    pressurisation = sim.add_step("Pressurisation", time=cfg["pressurisation_time"])
    pressurisation.bottom.set_pressure("inlet", sk.linear(cfg["high_pressure"]), y=FEED, temperature=298 * u.K)
    pressurisation.top.set_closed()

    # CO2 is fed in the feed and pressurisation steps and recovered by evacuation.
    sim.kpi.track("CarbonDioxide", feed=[feed.bottom["Feed inlet"], pressurisation.bottom],
                  product=[evacuation.bottom])


def add_lpp_cycle(sim: sk.Simulation, cfg: dict) -> None:
    """Feed split in two steps; part of the first one's light product pressurises the bed from the top."""
    feed_times = cfg["feed_times"]
    # Light product leaves through a valve (flow coefficient Cv, pressure-drop ratio factor xT) to 1 bar.
    ambient = dict(pressure=1 * u.bar, cv=10.0, xt=0.53, y=NO_GAS, temperature=293 * u.K)

    feed_1 = sim.add_step("Feed 1", time=feed_times[0])
    feed_1.bottom.set_streams("inlet", [
        sk.Stream("Feed 1 inlet", flow=sk.velocity(cfg["feed_velocity"]), y=FEED, temperature=298 * u.K),
    ])
    feed_1.top.set_streams("outlet", [
        sk.ValveStream("Feed 1 ambient Stream", **ambient),           # light product to the atmosphere
        # The LPP draw: its flow is whatever the connected pressurisation step takes.
        sk.Stream("LPP Stream", flow=sk.velocity(0 * u.m / u.s), y=NO_GAS, temperature=293 * u.K),
    ])

    feed_2 = sim.add_step("Feed 2", time=feed_times[1])
    feed_2.bottom.set_streams("inlet", [
        sk.Stream("Feed 2 inlet", flow=sk.velocity(cfg["feed_velocity"]), y=FEED, temperature=298 * u.K),
    ])
    feed_2.top.set_streams("outlet", [sk.ValveStream("Feed 2 ambient Stream", **ambient)])

    blowdown = sim.add_step("Blowdown", time=cfg["blowdown_time"])
    blowdown.bottom.set_closed()
    blowdown.top.set_pressure("outlet", sk.exponential(cfg["blowdown_pressure"]), y=BACKFLOW, temperature=293 * u.K)

    evacuation = sim.add_step("Evacuation", time=cfg["evacuation_time"])
    # Sect. 4.3: two-term exponential decay (0.019 and 0.101 1/s, Table 1); the split is not
    # printed and an equal weight is used (core/test/papers/perez_2019/README.md).
    evacuation.bottom.set_pressure("outlet", sk.double_exponential(cfg["evacuation_pressure"], rate_1=0.019, rate_2=0.101, weight=0.5),
                                   y=BACKFLOW, temperature=293 * u.K)
    evacuation.top.set_closed()

    pressurisation = sim.add_step("Pressurisation", time=cfg["pressurisation_time"])
    pressurisation.bottom.set_closed()
    pressurisation.top.set_pressure("inlet", sk.linear(cfg["high_pressure"]), y=BACKFLOW, temperature=298 * u.K)

    # Light-product pressurisation: the LPP stream of Feed 1 supplies the pressurisation inlet.
    sim.connect(feed_1.top["LPP Stream"], pressurisation.top)

    # CO2 is fed in both feed steps and recovered by evacuation. The LPP stream is N2-rich
    # gas moving within the cycle, not CO2 product, so it is not a KPI endpoint.
    sim.kpi.track("CarbonDioxide", feed=[feed_1.bottom["Feed 1 inlet"], feed_2.bottom["Feed 2 inlet"]],
                  product=[evacuation.bottom])


CASES = {name: f"Config{name}" for name in CONFIGURATIONS}   # case -> reference_values.csv case
QUANTITIES = [
    ("CO2 purity vs experiment", "CarbonDioxide.purity", "CO2_purity_exp"),
    ("CO2 recovery vs experiment", "CarbonDioxide.recovery", "CO2_recovery_exp"),
    ("CO2 purity vs paper model", "CarbonDioxide.purity", "CO2_purity_model"),
    ("CO2 recovery vs paper model", "CarbonDioxide.recovery", "CO2_recovery_model"),
]
DECKS = {name: f"Perez2019_{name}" for name in CONFIGURATIONS}   # the shipped deck of each case


def main(argv=None) -> int:
    return c.cyclic_family(argv, __doc__, "perez_2019", CASES, QUANTITIES, "CarbonDioxide", "perez_2019.png",
                           values=c.evacuation_product_kpis, build=build)


if __name__ == "__main__":
    sys.exit(main())
