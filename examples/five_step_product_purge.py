"""CO2/N2 VSA with product purge (template cycle, no literature source).

Builds the cycle from scratch in code: a 13X-type adsorbent, the column and
six steps (feed 1, feed 2, blowdown, evacuation, pressurisation and product
purge), two of them fed from another step's light product through endpoint
connections. Every value is set explicitly, with its unit. The cycle comes
from an integration test (core/test/integration/test_fluid_solver.cpp) and has
no KPI configuration, so the script reports the run statistics, the CSS check
and the CO2 mole fraction at both ends over the last cycle. Use it as a
starting point for your own cycles, not as data.

    python examples/five_step_product_purge.py

``build()`` returns the configured simulation without running it
(python/tests/test_examples_build.py checks it against examples/decks/).
"""

from __future__ import annotations

import sys

import _common as c
import skarstrom as sk

FEED = {"CarbonDioxide": 0.15, "Nitrogen": 0.85}   # flue gas, mole fractions
N2 = {"CarbonDioxide": 0.0, "Nitrogen": 1.0}        # pressurisation/purge gas and backflow composition


def add_13x(sim: sk.Simulation) -> None:
    """A 13X-type adsorbent: dual-site Langmuir in partial pressure, constant LDF coefficients."""
    ads = sim.add_adsorbent("13X", solid_density="1200 kg/m3", heat_capacity="500 J/(kg*K)",
                            thermal_conductivity="0.5 W/(m*K)", particle_void_fraction=0.37)
    ads.set_isotherm("dual_site_langmuir_p", {
        "CarbonDioxide": {"m1": "3.16 mol/kg", "b1": "9.23e-13 1/Pa", "q1": "47612.86 J/mol",
                          "m2": "3.16 mol/kg", "b2": "4.88e-13 1/Pa", "q2": "38545.63 J/mol"},
        "Nitrogen":      {"m1": "1.2 mol/kg", "b1": "5.85e-8 1/Pa", "q1": "16059.73 J/mol",
                          "m2": "1.2 mol/kg", "b2": "5.85e-8 1/Pa", "q2": "16059.73 J/mol"},
    }, components=["CarbonDioxide", "Nitrogen"])
    ads.set_rate_model("CarbonDioxide", "ldf", "constant", k="0.1 1/s")
    ads.set_rate_model("Nitrogen", "ldf", "constant", k="0.1 1/s")
    ads.heat_of_adsorption["CarbonDioxide"] = "40 kJ/mol"
    ads.heat_of_adsorption["Nitrogen"] = "0 J/mol"


def build(case: str = "FiveStepProductPurge") -> sk.Simulation:
    """The template cycle (the only case), ready to run."""
    if case != "FiveStepProductPurge":
        raise KeyError(f"no case {case!r}; the only case is 'FiveStepProductPurge'")
    sim = sk.Simulation.new(title="", description="", author="")

    # ---- fluid ------------------------------------------------------------------------------
    sim.add_component("CarbonDioxide")
    sim.add_component("Nitrogen")
    sim.fluid.adsorbed_heat_capacity_mode = "constant"
    sim.fluid.heat_of_adsorption_convention = "isosteric"
    sim.fluid.adsorbed_heat_capacity["CarbonDioxide"] = "0 J/(mol*K)"
    sim.fluid.adsorbed_heat_capacity["Nitrogen"] = "0 J/(mol*K)"

    # ---- adsorbent ----------------------------------------------------------------------------
    add_13x(sim)

    # ---- column: 1 m x 10 cm, the wall held at 293 K -------------------------------------------
    sim.bed.set_shape("cylinder", diameter="0.1 m")
    sim.bed.ambient_temperature = "293 K"
    sim.bed.set_thermal_equilibrium("lte")                     # gas and solid at one temperature
    sim.bed.set_particle_heat_transfer("constant", particle_htc="100 W/(m^2*K)")
    for void in (sim.bed.bottom_void, sim.bed.top_void):       # no dead volume at the column ends
        void.enabled = False
        void.volume = "1e-4 m3"
        void.resistance_coefficient = "0 Pa*s/m"
    sim.wall.set_material("ss316l", density="7800 kg/m3", heat_capacity="537 J/(kg*K)",
                          thermal_conductivity="16.2 W/(m*K)")
    sim.wall.thickness = "17.5 mm"
    sim.wall.set_boundary("fixed", fixed_temperature="293 K")

    layer = sim.add_layer("13X", composition={"13X": 1.0}, momentum="kozeny_carman",
                          length="1 m", cells=20, bed_void_fraction=0.37, particle_diameter="1.1 mm")
    layer.set_dispersion("ruthven", molecular_diffusion_coefficient="1.6e-5 m2/s")
    layer.set_heat_dispersion("constant", heat_dispersion_coefficient="0 W/(m*K)")
    layer.set_internal_heat_transfer("constant", internal_htc="59.1 W/(m^2*K)")
    layer.set_external_heat_transfer("constant", external_htc="0 W/(m^2*K)")

    # ---- initial state: the bed full of N2 at 1 atm ---------------------------------------------
    sim.initial_conditions.pressure = "1 atm"
    sim.initial_conditions.temperature = "298 K"
    sim.initial_conditions.wall_temperature = "298 K"
    sim.initial_conditions.velocity = "0 m/s"
    sim.initial_conditions.y["Nitrogen"] = 1.0
    sim.initial_conditions.y["CarbonDioxide"] = 0.0

    # ---- cycle (100 s per step) ----------------------------------------------------------------
    sim.cycle.number_of_cycles = 5
    sim.cycle.step_time_mode = "time"
    step_time = "100 s"

    def feed_step(n: int, draw: str) -> sk.Step:
        """Feed at the bottom; light product to the atmosphere through a valve, plus a draw stream."""
        step = sim.add_step(f"Feed{n}", time=step_time)
        step.bottom.set_streams("inlet", [
            sk.Stream(f"Feed {n} inlet", flow=sk.velocity("0.5 m/s"), y=FEED, temperature="298 K"),
        ])
        step.top.set_streams("outlet", [
            sk.ValveStream("Ambient Stream", pressure="1 atm", cv=100.0, xt=0.53, y={}, temperature="293 K"),
            # The draw's flow is whatever the connected step takes from it.
            sk.Stream(draw, flow=sk.velocity("0 m/s"), y={}, temperature="293 K"),
        ])
        return step

    # 1-2. Feed: the first feed step supplies the pressurisation gas, the second the purge gas.
    feed1 = feed_step(1, "Feed 1 Press Stream")
    feed2 = feed_step(2, "Feed 2 Purge Stream")

    # 3. Blowdown through a valve at the top to 0.4 bar.
    blowdown = sim.add_step("Blowdown", time=step_time)
    blowdown.bottom.set_closed()
    blowdown.top.set_streams("outlet", [
        sk.ValveStream("Blowdown waste", pressure="0.4 bar", cv=1.0, xt=0.53, y={}, temperature="293 K"),
    ])

    # 4. Evacuation through a valve at the bottom to 0.1 bar: the CO2 product.
    evacuation = sim.add_step("Evacuation", time=step_time)
    evacuation.bottom.set_streams("outlet", [
        sk.ValveStream("Product stream", pressure="0.1 bar", cv=1.0, xt=0.53, y={}, temperature="293 K"),
    ])
    evacuation.top.set_closed()

    # 5. Pressurisation from the top with light product from Feed1.
    press = sim.add_step("Press", time=step_time)
    press.bottom.set_closed()
    press.top.set_pressure("inlet", sk.exponential("1 atm"), y=N2, temperature="298 K")

    # 6. Product purge: light product from Feed2 in at the top, out of the bottom at 1 atm.
    purge = sim.add_step("Purge", time=step_time)
    purge.bottom.set_pressure("outlet", "1 atm", y=N2, temperature="293 K")
    purge.top.set_streams("inlet", [
        sk.Stream("Purge stream", flow=sk.molar_flow("0 mol/s"), y=N2, temperature="298 K"),
    ])

    # Connections: each draw stream supplies the matching inlet.
    sim.connect(feed1.top["Feed 1 Press Stream"], press.top)
    sim.connect(feed2.top["Feed 2 Purge Stream"], purge.top["Purge stream"])

    # ---- solver (SUNDIALS IDA), cyclic steady state, energy balance --------------------------------
    sim.settings.init_step = "1e-10 s"
    sim.settings.min_step = "1e-15 s"
    sim.settings.max_step = "0 s"                              # no limit
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
    sim.css.method = "none"                                    # run all cycles
    sim.css.min_cycles = 5
    sim.css.molar_absolute_tolerance = 1e-4
    sim.css.molar_relative_tolerance = 1e-4
    sim.css.mass_balance_tolerance = 0.005
    sim.energy_balance.energy_balance_mode = "non_isothermal"
    sim.energy_balance.isothermal_use_initial_temperature = True
    sim.energy_balance.isothermal_temperature = "298.15 K"
    return sim


DECKS = {"FiveStepProductPurge": "FiveStepProductPurge"}   # the shipped deck


def main(argv=None) -> int:
    args = c.parser(__doc__.splitlines()[0], ["FiveStepProductPurge"]).parse_args(argv)
    sim = build()
    if args.cycles is not None:
        sim.cycle.number_of_cycles = int(args.cycles)
    result = c.run(sim, "FiveStepProductPurge", args.quiet)
    print(f"  steps: {', '.join(sim.steps.names)}")
    print(f"  {c.css_status(sim)}")
    print(f"  integrator: {result.stats['steps']} steps, {result.stats['residual_evaluations']} residual evaluations")
    c.print_kpis(sim)
    plt = c.pyplot(args)
    if plt:
        fig, ax = plt.subplots(figsize=(6, 4))
        for loc in ("bottom", "top"):
            y = sim.series("mole_fraction", "CarbonDioxide", location=loc, cycle=-1)
            ax.plot(y.time, y.values, label=f"y_CO2 at the {loc}")
        ax.set(xlabel="time in the last cycle [s]", ylabel="CO2 mole fraction [-]")
        ax.legend(fontsize=7)
        c.save(fig, "five_step_product_purge.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
