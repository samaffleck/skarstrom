"""Ward & Pini (2024): CO2/N2 4-step VSA on Mg-MOF-74.

Builds the cycle from scratch in code: the fluid, the Mg-MOF-74 adsorbent and
its competitive isotherm, the column and its heat transfer, the four steps
(feed, blowdown, evacuation, pressurisation) with their boundary conditions,
the KPI definition and the solver settings. Every value is set explicitly,
with its unit, so this file is also a template for a VSA cycle of your own.

The cycle runs to cyclic steady state; CO2 purity and recovery are compared
with the paper's optimised detailed-model result (Table 4, maximum separation
effectiveness, single adsorbent).

    python examples/ward_pini_2024.py

``build("FourStepPSA")`` returns the configured simulation without running it
(python/tests/test_examples_build.py checks it against examples/decks/).

Reference: A. Ward, R. Pini, Industrial & Engineering Chemistry Research 63
(2024) 13787-13800. doi:10.1021/acs.iecr.4c01345
"""

from __future__ import annotations

import sys

import _common as c
import skarstrom as sk
from skarstrom import ureg as u  # pint units: 3.2 * u.cm, 40.98 * u.kJ / u.mol

FEED = {"CarbonDioxide": 0.15, "Nitrogen": 0.85}       # flue gas, mole fractions
BACKFLOW = {"CarbonDioxide": 0.0, "Nitrogen": 1.0}     # gas entering an outlet if the flow reverses


def add_mg_mof_74(sim: sk.Simulation) -> sk.Adsorbent:
    """Mg-MOF-74 with the extended (competitive) dual-site Langmuir isotherm of Table 1."""
    ads = sim.add_adsorbent(
        "Mg-MOF-74",
        solid_density=910.4166666666667 * u.kg / u.m**3,       # particle density
        heat_capacity=900 * u.J / (u.kg * u.K),
        thermal_conductivity=0.5 * u.W / (u.m * u.K),
        particle_void_fraction=0.35,
    )
    ads.set_isotherm("extended_dual_site_langmuir_c", {
        "CarbonDioxide": {"m1": 6.8 * u.mol / u.kg, "b1": 1.81e-7 * u.m**3 / u.mol, "q1": 39.3 * u.kJ / u.mol,
                          "m2": 9.9 * u.mol / u.kg, "b2": 1.06e-6 * u.m**3 / u.mol, "q2": 21.2 * u.kJ / u.mol},
        # N2 adsorbs on one site only.
        "Nitrogen":      {"m1": 14 * u.mol / u.kg,  "b1": 3.45e-6 * u.m**3 / u.mol, "q1": 15.5 * u.kJ / u.mol,
                          "m2": 0 * u.mol / u.kg,   "b2": 0 * u.m**3 / u.mol,       "q2": 0 * u.J / u.mol},
    }, components=["CarbonDioxide", "Nitrogen"])
    # Linear driving force with a constant coefficient.
    for component in ("CarbonDioxide", "Nitrogen"):
        ads.set_rate_model(component, "ldf", "constant", k=1 / u.s)
    ads.heat_of_adsorption["CarbonDioxide"] = 30.85 * u.kJ / u.mol
    ads.heat_of_adsorption["Nitrogen"] = 15.93 * u.kJ / u.mol
    return ads


def add_activated_carbon(sim: sk.Simulation) -> sk.Adsorbent:
    """The paper's second material; defined in the deck but not packed in this single-adsorbent case."""
    return sim.add_adsorbent(
        "Activated Carbon",
        solid_density=1200 * u.kg / u.m**3,
        heat_capacity=500 * u.J / (u.kg * u.K),
        thermal_conductivity=0.5 * u.W / (u.m * u.K),
        particle_void_fraction=0.35,
    )


def build(case: str = "FourStepPSA") -> sk.Simulation:
    """The four-step VSA cycle of the single-adsorbent case, ready to run."""
    if case != "FourStepPSA":
        raise KeyError(f"no case {case!r}; the only case is 'FourStepPSA'")
    sim = sk.Simulation.new(
        title="Ward 2024",
        description="Four-step VSA for CO2/N2 separation on Mg-MOF-74, based on Ward & Pini (2024).",
        author="Skarstrom examples",
    )

    # ---- fluid ---------------------------------------------------------------------
    sim.add_component("CarbonDioxide")
    sim.add_component("Nitrogen")
    sim.fluid.adsorbed_heat_capacity_mode = "constant"
    sim.fluid.heat_of_adsorption_convention = "internal_energy"
    sim.fluid.adsorbed_heat_capacity["CarbonDioxide"] = 30.7 * u.J / (u.mol * u.K)
    sim.fluid.adsorbed_heat_capacity["Nitrogen"] = 30.7 * u.J / (u.mol * u.K)

    # ---- adsorbents ------------------------------------------------------------------
    add_mg_mof_74(sim)
    add_activated_carbon(sim)

    # ---- column: 1 m x 0.29 m, stainless steel ----------------------------------------------
    sim.bed.set_shape("cylinder", diameter=0.29 * u.m)
    sim.bed.ambient_temperature = 293 * u.K
    sim.bed.set_thermal_equilibrium("lte")                     # gas and solid at one temperature
    sim.bed.set_particle_heat_transfer("constant", particle_htc=100 * u.W / (u.m**2 * u.K))
    for void in (sim.bed.bottom_void, sim.bed.top_void):       # no dead volume at the column ends
        void.enabled = False
        void.volume = 1e-4 * u.m**3
        void.resistance_coefficient = 0 * u.Pa * u.s / u.m
    sim.wall.set_material("ss316l", density=7800 * u.kg / u.m**3, heat_capacity=502 * u.J / (u.kg * u.K),
                          thermal_conductivity=16 * u.W / (u.m * u.K))
    sim.wall.thickness = 17 * u.mm
    sim.wall.set_boundary("zero_flux")

    layer = sim.add_layer("Adsorbent", composition={"Mg-MOF-74": 1.0}, momentum="kozeny_carman",
                          length=1 * u.m, cells=30, bed_void_fraction=0.52, particle_diameter=2 * u.mm)
    layer.set_dispersion("ruthven", molecular_diffusion_coefficient=1.6e-5 * u.m**2 / u.s)
    layer.set_heat_dispersion("constant", heat_dispersion_coefficient=0 * u.W / (u.m * u.K))
    layer.set_internal_heat_transfer("constant", internal_htc=8.6 * u.W / (u.m**2 * u.K))
    layer.set_external_heat_transfer("constant", external_htc=2.5 * u.W / (u.m**2 * u.K))

    # ---- initial state: the bed full of N2 at 1 bar and 298 K --------------------------
    sim.initial_conditions.pressure = 1 * u.bar
    sim.initial_conditions.temperature = 298 * u.K
    sim.initial_conditions.wall_temperature = 298 * u.K
    sim.initial_conditions.velocity = 0 * u.m / u.s
    sim.initial_conditions.y["Nitrogen"] = 1.0
    sim.initial_conditions.y["CarbonDioxide"] = 0.0

    # ---- cycle: feed, blowdown (top), evacuation (bottom), pressurisation with feed gas ----
    sim.cycle.number_of_cycles = 250
    sim.cycle.step_time_mode = "time"
    feed = sim.add_step("Feed", time=57.6 * u.s)
    feed.bottom.set_streams("inlet", [
        sk.Stream("Feed inlet", flow=sk.velocity(0.7956 * u.m / u.s), y=FEED, temperature=298.15 * u.K),
    ])
    feed.top.set_pressure("outlet", 1 * u.bar, y=BACKFLOW, temperature=293 * u.K)

    blowdown = sim.add_step("Blowdown", time=60.7 * u.s)
    blowdown.bottom.set_closed()
    blowdown.top.set_pressure("outlet", sk.exponential(9050 * u.Pa), y=BACKFLOW, temperature=293 * u.K)

    evacuation = sim.add_step("Evacuation", time=200 * u.s)
    evacuation.bottom.set_pressure("outlet", sk.exponential(1000 * u.Pa), y=BACKFLOW, temperature=293 * u.K)
    evacuation.top.set_closed()

    pressurisation = sim.add_step("Pressurisation", time=30 * u.s)
    pressurisation.bottom.set_pressure("inlet", sk.linear(1 * u.bar), y=FEED, temperature=298.15 * u.K)
    pressurisation.top.set_closed()

    # CO2 is fed in the feed and pressurisation steps and recovered by evacuation.
    sim.kpi.track("CarbonDioxide", feed=[feed.bottom["Feed inlet"], pressurisation.bottom],
                  product=[evacuation.bottom])

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


CASES = {"FourStepPSA": "FourStepPSA"}                  # case -> reference_values.csv case
DECKS = {"FourStepPSA": "Ward2024"}                     # the shipped deck of each case
QUANTITIES = [
    ("CO2 purity", "CarbonDioxide.purity", "CO2_purity"),
    ("CO2 recovery", "CarbonDioxide.recovery", "CO2_recovery"),
]


def main(argv=None) -> int:
    return c.cyclic_family(argv, __doc__, "ward_pini_2024", CASES, QUANTITIES, "CarbonDioxide", "ward_pini_2024.png",
                           build=build)


if __name__ == "__main__":
    sys.exit(main())
