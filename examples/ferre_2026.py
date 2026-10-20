"""Ferre et al. (2026): biogas upgrading (CO2/CH4) by 4-step PVSA on zeolite NaY.

Builds the NaY (Feng isotherm data) cycle from scratch in code: the biogas
components, the adsorbent with its competitive dual-site Langmuir isotherm, the
column, the four steps at the paper's optimum operating point and the solver
settings. Every value is set explicitly, with its unit.

Runs the cycle to cyclic steady state and compares CO2 and CH4 purity and
recovery with the paper's detailed-model result (Table 5, recovery-recovery
case).

    python examples/ferre_2026.py

``build("NaY_F")`` returns the configured simulation without running it
(python/tests/test_examples_build.py checks it against examples/decks/).

Reference: A. Ferre, C. Wedler, A. Ward, J. Burger, R. Pini, C. Petit,
Separation and Purification Technology 406 (2026) 138821.
doi:10.1016/j.seppur.2026.138821
"""

from __future__ import annotations

import sys

import _common as c
import skarstrom as sk
from skarstrom import ureg as u  # pint units: 3.2 * u.cm, 40.98 * u.kJ / u.mol

BIOGAS = {"CarbonDioxide": 0.4, "Methane": 0.6}   # feed (and initial) gas, mole fractions

# Table 5 (recovery-recovery, maximum CO2 recovery) operating point of NaY(F)
# (core/test/papers/ferre_2026/ferre_2026.cpp). The paper gives the feed as an
# interstitial velocity; the superficial velocity is that times the bed voidage.
OPERATING_POINT = dict(
    high_pressure=4.6 * u.bar, intermediate_pressure=0.7 * u.bar, low_pressure=0.03 * u.bar,
    feed_time=40.6 * u.s, blowdown_time=87.5 * u.s, evacuation_time=183.2 * u.s, pressurisation_time=20 * u.s,
    interstitial_velocity=0.5,      # [m/s]
)
BED_VOID_FRACTION = 0.37


def add_nay(sim: sk.Simulation) -> None:
    """Zeolite NaY(F): SI Tables S1-S3 of the paper."""
    ads = sim.add_adsorbent("NaY(F)", solid_density=1330 * u.kg / u.m**3, heat_capacity=1070 * u.J / (u.kg * u.K),
                            thermal_conductivity=0.5 * u.W / (u.m * u.K), particle_void_fraction=0.35)
    # Extended (competitive) dual-site Langmuir, concentration form (Sect. 3.1).
    ads.set_isotherm("extended_dual_site_langmuir_c", {
        "CarbonDioxide": {"m1": 6.39 * u.mol / u.kg, "b1": 1.23e-6 * u.m**3 / u.mol, "q1": 27028 * u.J / u.mol,
                          "m2": 0 * u.mol / u.kg, "b2": 2.1e-6 * u.m**3 / u.mol, "q2": 21169 * u.J / u.mol},
        "Methane":       {"m1": 4.09 * u.mol / u.kg, "b1": 8.16e-6 * u.m**3 / u.mol, "q1": 12963 * u.J / u.mol,
                          "m2": 2.85 * u.mol / u.kg, "b2": 5.39e-6 * u.m**3 / u.mol, "q2": 12756 * u.J / u.mol},
    }, components=["CarbonDioxide", "Methane"])
    # LDF with the macropore-diffusion constant 15 ep Dm / (tau rp^2) ~ 46.7 1/s. The paper
    # scales it by the capacity factor c/q*; the separation is equilibrium-controlled, so
    # the constant reproduces the same KPIs (see the validation test).
    ads.set_rate_model("CarbonDioxide", "ldf", "constant", k=46.7 / u.s)
    ads.set_rate_model("Methane", "ldf", "constant", k=46.7 / u.s)
    # Constant heats of adsorption: the DSL isosteric heat at the feed loading.
    ads.heat_of_adsorption["CarbonDioxide"] = 27028 * u.J / u.mol
    ads.heat_of_adsorption["Methane"] = 12896 * u.J / u.mol


def build(case: str = "NaY_F") -> sk.Simulation:
    """The NaY(F) cycle (the only case), ready to run."""
    if case != "NaY_F":
        raise KeyError(f"no case {case!r}; the only case is 'NaY_F'")
    op = OPERATING_POINT
    sim = sk.Simulation.new(title="Ferre 2026 - NaY(F)", description="", author="Ferre et al. (2026)")

    # ---- fluid ----------------------------------------------------------------------------
    sim.add_component("CarbonDioxide")
    sim.add_component("Methane")
    sim.fluid.adsorbed_heat_capacity_mode = "constant"
    # The paper's detailed model uses internal-energy heats of adsorption (Ward & Pini 2022).
    sim.fluid.heat_of_adsorption_convention = "internal_energy"
    sim.fluid.adsorbed_heat_capacity["CarbonDioxide"] = 30.7 * u.J / (u.mol * u.K)
    sim.fluid.adsorbed_heat_capacity["Methane"] = 30.7 * u.J / (u.mol * u.K)

    # ---- adsorbent ----------------------------------------------------------------------------
    add_nay(sim)

    # ---- column: 1 m long, 28.9 cm inner diameter, steel wall to r = 0.162 m (Table 2) ------
    r_in, r_out = 0.1445, 0.162                                # [m]
    sim.bed.set_shape("cylinder", diameter=(2 * r_in) * u.m)
    sim.bed.ambient_temperature = 298.15 * u.K
    sim.bed.set_thermal_equilibrium("lte")                     # gas and solid at one temperature
    sim.bed.set_particle_heat_transfer("constant", particle_htc=100 * u.W / (u.m**2 * u.K))
    for void in (sim.bed.bottom_void, sim.bed.top_void):       # no dead volume at the column ends
        void.enabled = False
        void.volume = 1e-4 * u.m**3
        void.resistance_coefficient = 0 * u.Pa * u.s / u.m
    sim.wall.set_material("ss316l", density=7800 * u.kg / u.m**3, heat_capacity=502 * u.J / (u.kg * u.K),
                          thermal_conductivity=16 * u.W / (u.m * u.K))
    sim.wall.thickness = (r_out - r_in) * u.m
    sim.wall.set_boundary("zero_flux")

    layer = sim.add_layer("Adsorbent", composition={"NaY(F)": 1.0}, momentum="kozeny_carman",
                          length=1 * u.m, cells=30, bed_void_fraction=BED_VOID_FRACTION, particle_diameter=1.5 * u.mm)
    layer.set_dispersion("ruthven", molecular_diffusion_coefficient=1.6e-5 * u.m**2 / u.s)
    layer.set_heat_dispersion("constant", heat_dispersion_coefficient=0 * u.W / (u.m * u.K))
    layer.set_internal_heat_transfer("constant", internal_htc=8.6 * u.W / (u.m**2 * u.K))
    layer.set_external_heat_transfer("constant", external_htc=2.5 * u.W / (u.m**2 * u.K))

    # ---- initial state: the bed full of biogas at the high pressure ----------------------------
    sim.initial_conditions.pressure = op["high_pressure"]
    sim.initial_conditions.temperature = 298.15 * u.K
    sim.initial_conditions.wall_temperature = 298.15 * u.K
    sim.initial_conditions.velocity = 0 * u.m / u.s
    sim.initial_conditions.y["Methane"] = BIOGAS["Methane"]
    sim.initial_conditions.y["CarbonDioxide"] = BIOGAS["CarbonDioxide"]

    # ---- cycle: feed, blowdown (top), evacuation (bottom), feed pressurisation -----------------
    sim.cycle.number_of_cycles = 250
    sim.cycle.step_time_mode = "time"
    feed_velocity = op["interstitial_velocity"] * BED_VOID_FRACTION   # superficial [m/s]

    # 1. Feed: biogas in at the bottom, CH4-rich raffinate (light product) out of the top.
    feed = sim.add_step("Feed", time=op["feed_time"])
    feed.bottom.set_streams("inlet", [
        sk.Stream("", flow=sk.velocity(feed_velocity * u.m / u.s), y=BIOGAS, temperature=298.15 * u.K),
    ])
    feed.top.set_pressure("outlet", op["high_pressure"], y=BIOGAS, temperature=293 * u.K)

    # 2. Blowdown to the intermediate pressure through the top.
    blowdown = sim.add_step("Blowdown", time=op["blowdown_time"])
    blowdown.bottom.set_closed()
    blowdown.top.set_pressure("outlet", sk.exponential(op["intermediate_pressure"]), y=BIOGAS, temperature=293 * u.K)

    # 3. Evacuation through the bottom: the CO2 heavy product.
    evacuation = sim.add_step("Evacuation", time=op["evacuation_time"])
    evacuation.bottom.set_pressure("outlet", sk.exponential(op["low_pressure"]), y=BIOGAS, temperature=293 * u.K)
    evacuation.top.set_closed()

    # 4. Pressurisation with feed gas from the bottom.
    pressurisation = sim.add_step("Pressurisation", time=op["pressurisation_time"])
    pressurisation.bottom.set_pressure("inlet", sk.linear(op["high_pressure"]), y=BIOGAS, temperature=298.15 * u.K)
    pressurisation.top.set_closed()
    # (No KPI configuration: the script integrates the stream logs, see kpis().)

    # ---- solver (SUNDIALS IDA), cyclic steady state, energy balance -----------------------------
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


CASES = {"NaY_F": "NaY_F"}                 # case -> reference_values.csv case
DECKS = {"NaY_F": "Ferre2026_NaYF"}        # the shipped deck
QUANTITIES = [
    ("CO2 purity", "CarbonDioxide.purity", "CO2_purity"),
    ("CO2 recovery", "CarbonDioxide.recovery", "CO2_recovery"),
    ("CH4 purity", "Methane.purity", "CH4_purity"),
    ("CH4 recovery", "Methane.recovery", "CH4_recovery"),
]


def kpis(sim):
    """The cycle has no KPI configuration; integrate the stream logs of the last cycle.

    Steps: 0 Feed, 1 Blowdown, 2 Evacuation, 3 Pressurisation. CO2 is the heavy
    product (evacuation, bottom), CH4 the light product (feed raffinate, top);
    the feed is the feed and pressurisation inlets (bottom), as in the test.
    """
    m = c.endpoint_moles
    out = {}
    for comp, (step, loc) in (("CarbonDioxide", (2, "bottom")), ("Methane", (0, "top"))):
        fed = m(sim, 0, "bottom", comp) + m(sim, 3, "bottom", comp)
        out[f"{comp}.purity"] = m(sim, step, loc, comp) / m(sim, step, loc)
        out[f"{comp}.recovery"] = m(sim, step, loc, comp) / fed
    return out


def main(argv=None) -> int:
    return c.cyclic_family(argv, __doc__, "ferre_2026", CASES, QUANTITIES, "Methane", "ferre_2026.png",
                           values=kpis, build=build)


if __name__ == "__main__":
    sys.exit(main())
