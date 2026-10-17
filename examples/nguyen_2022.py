"""Nguyen et al. (2022): CO2 breakthrough on CALF-20 (Figs. 4 and 5).

Builds each case from scratch in code: the CO2/N2/He fluid, CALF-20 with the
paper's dual-site Langmuir isotherm, the column (sized to hold 16.7 g of
adsorbent), the single feed step and the solver settings. Every value is set
explicitly, with its unit.

Fig. 4 is a 70 % CO2 feed and Fig. 5 a 15 % CO2 feed (balance helium) into a
bed initially full of helium. The outlet CO2 mole fraction is compared with the
digitised experiments (core/data/nguyen_2022/): NRMSE normalised by the feed
mole fraction, and the time at which the outlet reaches half the feed value.

    python examples/nguyen_2022.py
    python examples/nguyen_2022.py --case Fig4

``build("Fig4")`` returns the configured simulation without running it
(python/tests/test_examples_build.py checks it against examples/decks/).

Reference: T.T.T. Nguyen, J.-B. Lin, G.K.H. Shimizu, A. Rajendran, Chemical
Engineering Journal 442 (2022) 136263. doi:10.1016/j.cej.2022.136263
"""

from __future__ import annotations

import math
import sys

import _common as c
import skarstrom as sk
from skarstrom import ureg as u  # pint units: 3.2 * u.cm, 40.98 * u.kJ / u.mol

# The two figures differ only in the feed: composition, flow rate and how long it runs.
FIGURES = {
    "Fig4": dict(title="Nguyen 2022 - Fig. 4", feed_co2=0.70, feed_flow=50 * u.mL / u.min, feed_time=6000 * u.s,
                 description="Binary CO2/He breakthrough on CALF-20 (Fig. 4), inlet CO2 mole fraction 0.7."),
    "Fig5": dict(title="Nguyen 2022 - Fig. 5", feed_co2=0.15, feed_flow=100 * u.mL / u.min, feed_time=12000 * u.s,
                 description="Ternary CO2/N2/He breakthrough on CALF-20 (Fig. 5), inlet CO2 mole fraction 0.15."),
}

# Column geometry (core/test/papers/nguyen_2022/nguyen_2022.h): the bed length
# is chosen so that it holds exactly 16.7 g of CALF-20, with rho_b = rho_p (1 - eb).
ADSORBENT_MASS = 16.7e-3        # [kg]
INNER_DIAMETER = 0.0282         # [m]
OUTER_DIAMETER = 0.0318         # [m]
PARTICLE_DENSITY = 570.0        # [kg/m3]
BED_VOID_FRACTION = 0.4
BED_LENGTH = ADSORBENT_MASS / (PARTICLE_DENSITY * (1 - BED_VOID_FRACTION)) / (math.pi / 4 * INNER_DIAMETER**2)
WALL_THICKNESS = (OUTER_DIAMETER - INNER_DIAMETER) / 2


def build(case: str) -> sk.Simulation:
    """The simulation of one figure ("Fig4" or "Fig5"), ready to run."""
    fig = FIGURES[case]
    sim = sk.Simulation.new(title=fig["title"], description=fig["description"], author="Skarstrom examples")

    # ---- fluid: CO2 and N2 adsorb, helium is the inert carrier ---------------------
    for component in ("CarbonDioxide", "Nitrogen", "Helium"):
        sim.add_component(component)
    sim.fluid.adsorbed_heat_capacity_mode = "constant"
    sim.fluid.heat_of_adsorption_convention = "isosteric"
    sim.fluid.adsorbed_heat_capacity["CarbonDioxide"] = 37.435 * u.J / (u.mol * u.K)
    sim.fluid.adsorbed_heat_capacity["Nitrogen"] = 29.171 * u.J / (u.mol * u.K)
    sim.fluid.adsorbed_heat_capacity["Helium"] = 0 * u.J / (u.mol * u.K)

    # ---- adsorbent: CALF-20, dual-site Langmuir from the ESS model of Table 3 ---------
    calf20 = sim.add_adsorbent("CALF-20", solid_density=PARTICLE_DENSITY * u.kg / u.m**3,
                               heat_capacity=1371 * u.J / (u.kg * u.K), thermal_conductivity=0.63 * u.W / (u.m * u.K),
                               particle_void_fraction=0.35)
    calf20.set_isotherm("extended_dual_site_langmuir_c", {
        "CarbonDioxide": {"m1": 2.387 * u.mol / u.kg, "b1": 5.519e-7 * u.m**3 / u.mol, "q1": 35.06 * u.kJ / u.mol,
                          "m2": 3.271 * u.mol / u.kg, "b2": 5.187e-8 * u.m**3 / u.mol, "q2": 28.95 * u.kJ / u.mol},
        "Nitrogen":      {"m1": 2.387 * u.mol / u.kg, "b1": 8.139e-7 * u.m**3 / u.mol, "q1": 17.96 * u.kJ / u.mol,
                          "m2": 3.271 * u.mol / u.kg, "b2": 8.139e-7 * u.m**3 / u.mol, "q2": 17.96 * u.kJ / u.mol},
    }, components=["CarbonDioxide", "Nitrogen"])
    calf20.set_isotherm("inert", components=["Helium"])
    calf20.set_rate_model("CarbonDioxide", "ldf", "macropore", tortuosity=3.0, macropore_void_fraction=0.35)
    calf20.set_rate_model("Nitrogen", "ldf", "macropore", tortuosity=3.0, macropore_void_fraction=0.35)
    calf20.set_rate_model("Helium", "ldf", "constant", k=0.5 / u.s)
    calf20.heat_of_adsorption["CarbonDioxide"] = 38 * u.kJ / u.mol
    calf20.heat_of_adsorption["Nitrogen"] = 15 * u.kJ / u.mol
    calf20.heat_of_adsorption["Helium"] = 0 * u.kJ / u.mol

    # ---- column ------------------------------------------------------------------------
    sim.bed.set_shape("cylinder", diameter=INNER_DIAMETER * u.m)
    # Not stated in the material available: 295 K is the initial thermocouple
    # reading of the digitised experiments (core/data/nguyen_2022/Fig*/T_CO2_*.csv).
    sim.bed.ambient_temperature = 295 * u.K
    sim.bed.set_thermal_equilibrium("lte")
    sim.bed.set_particle_heat_transfer("constant", particle_htc=100 * u.W / (u.m**2 * u.K))
    for void in (sim.bed.bottom_void, sim.bed.top_void):       # no dead volume at the column ends
        void.enabled = False
        void.volume = 1e-4 * u.m**3
        void.resistance_coefficient = 0 * u.Pa * u.s / u.m
    sim.wall.set_material("ss316l", density=7800 * u.kg / u.m**3, heat_capacity=502 * u.J / (u.kg * u.K),
                          thermal_conductivity=16 * u.W / (u.m * u.K))
    sim.wall.thickness = WALL_THICKNESS * u.m
    sim.wall.set_boundary("zero_flux")

    # Particle diameter: not stated in the material available; 1.5 mm is the deck author's assumption.
    layer = sim.add_layer("CALF-20", composition={"CALF-20": 1.0}, momentum="kozeny_carman",
                          length=BED_LENGTH * u.m, cells=30, bed_void_fraction=BED_VOID_FRACTION,
                          particle_diameter=1.5 * u.mm)
    layer.set_dispersion("ruthven", molecular_diffusion_coefficient=1.6e-5 * u.m**2 / u.s)
    layer.set_heat_dispersion("constant", heat_dispersion_coefficient=0.4 * u.W / (u.m * u.K))
    layer.set_internal_heat_transfer("constant", internal_htc=6.963 * u.W / (u.m**2 * u.K))
    layer.set_external_heat_transfer("constant", external_htc=14 * u.W / (u.m**2 * u.K))

    # ---- initial state: the bed full of helium at 0.97 bar ---------------------------------
    sim.initial_conditions.pressure = 0.97 * u.bar
    sim.initial_conditions.temperature = 295 * u.K
    sim.initial_conditions.wall_temperature = 295 * u.K
    sim.initial_conditions.velocity = 0 * u.m / u.s
    sim.initial_conditions.y["Helium"] = 1.0
    sim.initial_conditions.y["Nitrogen"] = 0.0
    sim.initial_conditions.y["CarbonDioxide"] = 0.0

    # ---- one feed step: breakthrough --------------------------------------------------------
    sim.cycle.number_of_cycles = 1
    sim.cycle.step_time_mode = "time"
    feed = sim.add_step("Feed", time=fig["feed_time"])
    feed_gas = {"Helium": 1.0 - fig["feed_co2"], "Nitrogen": 0.0, "CarbonDioxide": fig["feed_co2"]}
    feed.bottom.set_streams("inlet", [
        sk.Stream("CO2 Inlet", flow=sk.volumetric_flow(fig["feed_flow"]), y=feed_gas, temperature=295 * u.K),
    ])
    feed.top.set_pressure("outlet", 0.97 * u.bar, y={"CarbonDioxide": 0.0, "Nitrogen": 0.0, "Helium": 1.0},
                          temperature=293 * u.K)

    # ---- solver (SUNDIALS IDA); a single step, so no cyclic-steady-state check ----------------
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
    sim.css.method = "none"
    sim.css.min_cycles = 5
    sim.css.molar_absolute_tolerance = 1e-4
    sim.css.molar_relative_tolerance = 1e-4
    sim.css.mass_balance_tolerance = 0.005
    sim.energy_balance.energy_balance_mode = "non_isothermal"
    sim.energy_balance.isothermal_use_initial_temperature = True
    sim.energy_balance.isothermal_temperature = 298.15 * u.K
    return sim


CASES = {"Fig4": 70, "Fig5": 15}                                 # case -> feed CO2 [%] (digitised data file)
DECKS = {"Fig4": "Nguyen2022_Fig4", "Fig5": "Nguyen2022_Fig5"}   # the shipped deck of each case


def main(argv=None) -> int:
    args = c.parser(__doc__.splitlines()[0], list(CASES)).parse_args(argv)
    plt = c.pyplot(args)
    fig, ax = (plt.subplots(figsize=(6, 4)) if plt else (None, None))
    rows = []
    for name in c.selected(args):
        feed_pct = CASES[name]
        y_feed = feed_pct / 100.0
        sim = build(name)
        if args.cycles is not None:
            sim.cycle.number_of_cycles = int(args.cycles)
        c.run(sim, f"Nguyen 2022 {name}", args.quiet)
        y = sim.series("mole_fraction", "CarbonDioxide", location="top", cycle=0)
        t50 = c.crossing_time(y.time, y.values, 0.5 * y_feed)
        row = [name, f"{feed_pct} %", c.fmt(t50), "-", "-"]
        data = c.DATA / "nguyen_2022" / name / f"Y_CO2_{feed_pct}%.csv"
        if data.is_file():
            t_e, y_e = c.load_csv(data)
            row[3] = c.fmt(c.crossing_time(t_e, y_e, 0.5 * y_feed))
            row[4] = c.fmt(c.nrmse(y.time, y.values, t_e, y_e, y_feed))
            if ax is not None:
                ax.plot(t_e, y_e, "o", ms=3, label=f"{name} experiment")
        if ax is not None:
            ax.plot(y.time, y.values, "-", label=f"{name} Skarstrom")
        rows.append(row)
    print()
    c.table(rows, ["case", "feed CO2", "t50 simulated [s]", "t50 experiment [s]", "NRMSE y_CO2"])
    if fig is not None:
        ax.set(xlabel="time [s]", ylabel="outlet CO2 mole fraction [-]")
        ax.legend(fontsize=7)
        c.save(fig, "nguyen_2022.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
