"""Illustrative TSA: CO2 capture on zeolite 13X with indirect heating and cooling.

Not a reproduction of a paper: a four-step temperature swing cycle
(adsorption, electric heating with the CO2 product leaving the feed end,
cooling with both ends closed, feed pressurisation) driven by heater and
cooler energy sources, built for the TSA tutorial in the documentation.

The script builds the cycle from scratch in code: the fluid, the adsorbent
(the values of the adsorbent library entry "Zeolite 13X (Perez 2019)"), the
column, the four steps with their boundary conditions and energy sources, the
KPI definition and the solver settings. Every value is set explicitly, with
its unit. It runs the cycle to cyclic steady state, prints the CSS check and
the KPIs, and plots the bed temperature at both ends over the last cycle.

    python examples/tsa_illustrative.py
    python examples/tsa_illustrative.py --cycles 5

``build()`` returns the configured simulation without running it
(python/tests/test_examples_build.py checks it against examples/decks/).
"""

from __future__ import annotations

import sys

import _common as c
import skarstrom as sk
from skarstrom import ureg as u  # pint units: 3.2 * u.cm, 40.98 * u.kJ / u.mol

FEED = {"Nitrogen": 0.85, "CarbonDioxide": 0.15}     # flue gas, mole fractions


def build(case: str = "TSA") -> sk.Simulation:
    """The illustrative TSA cycle, ready to run (there is one case, "TSA")."""
    sim = sk.Simulation.new(
        title="Illustrative TSA - CO2/N2 on 13X",
        description="ILLUSTRATIVE, not a validated reproduction. Four-step indirectly heated TSA (adsorption, "
                    "electric heating with CO2 product at the feed end, cooling, feed pressurisation) for a "
                    "15 % CO2 / 85 % N2 flue gas.",
        author="Skarstrom examples",
    )

    # ---- fluid ---------------------------------------------------------------------
    sim.add_component("CarbonDioxide")
    sim.add_component("Nitrogen")
    sim.fluid.adsorbed_heat_capacity_mode = "constant"
    sim.fluid.heat_of_adsorption_convention = "isosteric"
    sim.fluid.adsorbed_heat_capacity["CarbonDioxide"] = 30.7 * u.J / (u.mol * u.K)
    sim.fluid.adsorbed_heat_capacity["Nitrogen"] = 30.7 * u.J / (u.mol * u.K)

    # ---- adsorbent: the library entry "Zeolite 13X (Perez 2019)", written out ------------
    # (sim.add_library_adsorbent("zeolite_13x_perez2019") adds the same values.)
    # The isotherm was fitted at 298-348 K; the heating step goes above that range.
    ads = sim.add_adsorbent("Zeolite 13X (Perez 2019)", solid_density=1138.0252816422599 * u.kg / u.m**3,
                            heat_capacity=900 * u.J / (u.kg * u.K), thermal_conductivity=0.63 * u.W / (u.m * u.K),
                            particle_void_fraction=0.35)
    ads.set_isotherm("extended_dual_site_langmuir_c", {
        "CarbonDioxide": {"m1": 2.768 * u.mol / u.kg, "b1": 3.934e-7 * u.m**3 / u.mol, "q1": 40.98 * u.kJ / u.mol,
                          "m2": 2.758 * u.mol / u.kg, "b2": 1.256e-7 * u.m**3 / u.mol, "q2": 34.19 * u.kJ / u.mol},
        "Nitrogen":      {"m1": 2.768 * u.mol / u.kg, "b1": 3.134e-6 * u.m**3 / u.mol, "q1": 15.71 * u.kJ / u.mol,
                          "m2": 2.758 * u.mol / u.kg, "b2": 3.134e-6 * u.m**3 / u.mol, "q2": 15.72 * u.kJ / u.mol},
    }, components=["CarbonDioxide", "Nitrogen"])
    for component in ("CarbonDioxide", "Nitrogen"):
        ads.set_rate_model(component, "ldf", "macropore", tortuosity=3.0, macropore_void_fraction=0.35)
    ads.heat_of_adsorption["CarbonDioxide"] = 30.85 * u.kJ / u.mol
    ads.heat_of_adsorption["Nitrogen"] = 15.93 * u.kJ / u.mol

    # ---- column: the Perez 2019 geometry, wall and heat transfer ---------------------------
    sim.bed.set_shape("cylinder", diameter=3.2 * u.cm)
    sim.bed.ambient_temperature = 298 * u.K
    sim.bed.set_thermal_equilibrium("lte")
    sim.bed.set_particle_heat_transfer("constant", particle_htc=100 * u.W / (u.m**2 * u.K))
    for void in (sim.bed.bottom_void, sim.bed.top_void):       # no dead volume at the column ends
        void.enabled = False
        void.volume = 1e-4 * u.m**3
        void.resistance_coefficient = 0 * u.Pa * u.s / u.m
    sim.wall.set_material("ss316l", density=7990 * u.kg / u.m**3, heat_capacity=500 * u.J / (u.kg * u.K),
                          thermal_conductivity=16.3 * u.W / (u.m * u.K))
    sim.wall.thickness = 3.1 * u.mm
    sim.wall.set_boundary("zero_flux")

    # Ergun pressure drop (the default for a new layer).
    layer = sim.add_layer("13X", composition={"Zeolite 13X (Perez 2019)": 1.0}, momentum="ergun",
                          length=0.3 * u.m, cells=30, bed_void_fraction=0.41, particle_diameter=1.5 * u.mm)
    layer.set_dispersion("ruthven", molecular_diffusion_coefficient=1.6e-5 * u.m**2 / u.s)
    layer.set_heat_dispersion("constant", heat_dispersion_coefficient=0.0903 * u.W / (u.m * u.K))
    layer.set_internal_heat_transfer("constant", internal_htc=8.6 * u.W / (u.m**2 * u.K))
    layer.set_external_heat_transfer("constant", external_htc=10 * u.W / (u.m**2 * u.K))

    # ---- initial state: the bed full of N2 at 1 bar and 298 K ------------------------------
    sim.initial_conditions.pressure = 1 * u.bar
    sim.initial_conditions.temperature = 298 * u.K
    sim.initial_conditions.wall_temperature = 298 * u.K
    sim.initial_conditions.velocity = 0 * u.m / u.s
    sim.initial_conditions.y["Nitrogen"] = 1.0
    sim.initial_conditions.y["CarbonDioxide"] = 0.0

    # ---- cycle (round numbers chosen for the tutorial) --------------------------------------
    sim.cycle.number_of_cycles = 60
    sim.cycle.step_time_mode = "time"

    # 1. Adsorption: flue gas in at the bottom, N2-rich raffinate out of the top.
    adsorption = sim.add_step("Adsorption", time=600 * u.s)
    adsorption.bottom.set_streams("inlet", [
        sk.Stream("Flue gas", flow=sk.velocity(0.05 * u.m / u.s), y=FEED, temperature=298 * u.K),
    ])
    adsorption.top.set_pressure("outlet", 1 * u.bar, y={"Nitrogen": 1.0, "CarbonDioxide": 0.0}, temperature=293 * u.K)

    # 2. Heating: top closed; CO2 leaves the feed end at 1 bar while a heater puts 30 W into the layer.
    heating = sim.add_step("Heating", time=900 * u.s)
    heating.bottom.set_pressure("outlet", 1 * u.bar, y={"Nitrogen": 0.0, "CarbonDioxide": 1.0}, temperature=293 * u.K)
    heating.top.set_closed()
    heating.add_energy_source(layer, power=30 * u.W, enabled=True)

    # 3. Cooling: both ends closed while 12 W is removed from the layer.
    cooling = sim.add_step("Cooling", time=900 * u.s)
    cooling.bottom.set_closed()
    cooling.top.set_closed()
    cooling.add_energy_source(layer, power=-12 * u.W, enabled=True)

    # 4. Pressurisation: feed gas back to 1 bar through the bottom, top closed.
    pressurisation = sim.add_step("Pressurisation", time=30 * u.s)
    pressurisation.bottom.set_pressure("inlet", sk.linear(1 * u.bar), y=FEED, temperature=298 * u.K)
    pressurisation.top.set_closed()

    # CO2 is fed in adsorption and pressurisation (bottom); the product is the heating-step outflow.
    sim.kpi.track("CarbonDioxide", feed=[adsorption.bottom["Flue gas"], pressurisation.bottom],
                  product=[heating.bottom])

    # ---- solver (SUNDIALS IDA), cyclic steady state, energy balance ----------------------------
    sim.settings.init_step = 1e-10 * u.s
    sim.settings.min_step = 1e-15 * u.s
    sim.settings.max_step = 0 * u.s                 # 0: no limit on the integrator step
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


DECKS = {"TSA": "Illustrative_TSA_13X"}       # the shipped deck


def main(argv=None) -> int:
    args = c.parser(__doc__.splitlines()[0], ["TSA"]).parse_args(argv)
    sim = build()
    if args.cycles is not None:
        sim.cycle.number_of_cycles = int(args.cycles)
    c.run(sim, "Illustrative TSA", args.quiet)
    print(f"  {c.css_status(sim)}")
    c.print_kpis(sim)
    plt = c.pyplot(args)
    if plt:
        fig, ax = plt.subplots(figsize=(6, 4))
        for loc in ("bottom", "top"):
            T = sim.series("temperature", location=loc, cycle=-1)
            ax.plot(T.time, T.values, label=f"bed temperature, {loc}")
        ax.set(xlabel="time in the last cycle [s]", ylabel="temperature [K]")
        ax.legend(fontsize=7)
        c.save(fig, "tsa_illustrative.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
