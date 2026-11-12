"""Streb & Mazzotti (2021): CO2/CH4 breakthrough on zeolite 13X with thermocouples.

Builds each case from scratch in code: the CO2/CH4(/H2)/He fluid, zeolite 13X with
the paper's extended-Sips competitive isotherm, the 1.2 m column, the sampling
line to the mass spectrometer, the single feed step and the solver settings.
Every value is set explicitly, with its unit.

The decks are calibrated on A3 (core/test/papers/streb_2021/README.md,
"Calibration on A3"): the LDF coefficients of CO2 and CH4, the wall-to-
thermofluid coefficient and the gas volume of the sampling line. Each case's
feed flow is its Table 2 set point times a flow factor fitted to the speed of
its own CO2 temperature front between the thermocouples at 35 and 110 cm.

The nine experiments of Table 2 with digitised data feed a 50:50 CO2/CH4 mixture
(A1-A7: 1-25 bar, 25-65 C, 10-40 cm3/s) or a 75:20:5 H2/CO2/CH4 mixture (B1, B2)
into a helium-filled column. The outlet CO2 and CH4 mole fractions and the five
centre-line thermocouples (z = 10, 35, 60, 85, 110 cm) are compared with the
digitised experiments (core/data/streb_2021/).

    python examples/streb_2021.py               # A1 and A3
    python examples/streb_2021.py --case B2     # any of the nine (--case repeatable)

``build("A3")`` returns the configured simulation without running it
(python/tests/test_examples_build.py checks it against examples/decks/).

Reference: A. Streb, M. Mazzotti, Adsorption 27 (2021) 541-558.
doi:10.1007/s10450-021-00306-y
"""

from __future__ import annotations

import math
import sys

import numpy as np

import _common as c
import skarstrom as sk
from skarstrom import ureg as u  # pint units: 3.2 * u.cm, 40.98 * u.kJ / u.mol

# The nine experiments of Table 2 with digitised data: column pressure, temperature
# (feed = initial = thermostat), flow set point (at column conditions), feed time
# (covers the digitised record) and feed (binary CO2:CH4 = 50:50, ternary
# H2:CO2:CH4 = 75:20:5). The flow factor is the actual flow over the set point, from
# the speed of the CO2 temperature front (A3 calibrated within 2 % of the set point;
# README "Calibration on A3").
EXPERIMENTS = {
    "A1": dict(pressure=1, temperature=298.0, flow=20, flow_factor=1.051, time=4100.0, ternary=False),
    "A2": dict(pressure=10, temperature=298.0, flow=10, flow_factor=0.939, time=1500.0, ternary=False),
    "A3": dict(pressure=10, temperature=298.0, flow=20, flow_factor=1.020, time=4100.0, ternary=False),
    "A4": dict(pressure=10, temperature=298.0, flow=40, flow_factor=0.906, time=1000.0, ternary=False),
    "A5": dict(pressure=25, temperature=298.0, flow=20, flow_factor=0.978, time=750.0, ternary=False),
    "A6": dict(pressure=10, temperature=318.15, flow=20, flow_factor=0.965, time=1050.0, ternary=False),
    "A7": dict(pressure=10, temperature=338.15, flow=20, flow_factor=0.995, time=1100.0, ternary=False),
    "B1": dict(pressure=1, temperature=298.0, flow=20, flow_factor=1.057, time=10400.0, ternary=True),
    "B2": dict(pressure=10, temperature=298.0, flow=10, flow_factor=0.887, time=3350.0, ternary=True),
}

# Densities of the 13X beads and the packed bed (Table 1) give both voidages.
SKELETAL_DENSITY = 2359.0       # [kg/m3]
PARTICLE_DENSITY = 1085.0       # [kg/m3]
BED_DENSITY = 595.0             # [kg/m3]
PARTICLE_VOID_FRACTION = 1 - PARTICLE_DENSITY / SKELETAL_DENSITY
BED_VOID_FRACTION = 1 - BED_DENSITY / PARTICLE_DENSITY

# Calibrated on A3 (README "Calibration on A3"; the paper's values in brackets).
LDF_CO2 = 0.0592 / u.s                          # (Table 4: 0.06 1/s)
LDF_CH4 = 0.137 / u.s                           # (Table 4: 0.2 1/s)
EXTERNAL_HTC = 49.7 * u.W / (u.m**2 * u.K)      # wall to thermofluid (Sect. 3.2.1: 220)
LINE_VOLUME = 0.142 * u.L                       # gas volume of the sampling line
LINE_VOID_FRACTION = 0.99
COLUMN_AREA = math.pi * (2.5 * u.cm) ** 2 / 4


def build(case: str) -> sk.Simulation:
    """The simulation of one experiment (a key of EXPERIMENTS), ready to run."""
    exp = EXPERIMENTS[case]
    ternary = exp["ternary"]
    T = exp["temperature"] * u.K
    P = exp["pressure"] * u.bar
    feed_name = "ternary H2/CO2/CH4 75:20:5" if ternary else "binary CO2/CH4 50:50"
    sim = sk.Simulation.new(title=f"Streb 2021 - Experiment {case}",
                            description=f"Breakthrough on zeolite 13X (Streb & Mazzotti 2021, experiment {case}, {feed_name}) "
                                        "into a column full of helium.",
                            author="Skarstrom examples")
    adsorbing = ["CarbonDioxide", "Methane"] + (["Hydrogen"] if ternary else [])
    components = adsorbing + ["Helium"]

    # ---- fluid: CO2, CH4 (and H2) adsorb, helium is the inert initial gas --------------
    for component in components:
        sim.add_component(component)
    sim.fluid.adsorbed_heat_capacity_mode = "constant"
    sim.fluid.heat_of_adsorption_convention = "isosteric"
    for component in components:
        sim.fluid.adsorbed_heat_capacity[component] = 0 * u.J / (u.mol * u.K)
    # The paper uses constant, average gas heat capacities and prints none; ours is
    # each component's ideal-gas Cp0 at the feed temperature.
    sim.fluid.gas_heat_capacity_mode = "constant"
    for component in components:
        sim.fluid.gas_heat_capacity[component] = sk.ideal_gas_heat_capacity(component, exp["temperature"])

    # ---- adsorbent: zeolite 13X, extended-Sips competitive isotherm (Table 3) ---------
    zeolite = sim.add_adsorbent("13X", solid_density=PARTICLE_DENSITY * u.kg / u.m**3,
                                heat_capacity=1100 * u.J / (u.kg * u.K), thermal_conductivity=0.15 * u.W / (u.m * u.K),
                                particle_void_fraction=PARTICLE_VOID_FRACTION)
    # Extended Sips with temperature-dependent capacity and exponent (Table 3).
    sips = {
        "CarbonDioxide": {"A1": 7.268 * u.mol / u.kg, "A2": -0.6168, "B1": 1.129e-4 / u.bar, "B2": 28.389 * u.kJ / u.mol,
                          "C1": 0.425, "C2": 0.7238, "T_ref": 298 * u.K},
        "Methane":       {"A1": 4.473 * u.mol / u.kg, "A2": -0.6569, "B1": 3.605e-4 / u.bar, "B2": 15.561 * u.kJ / u.mol,
                          "C1": 1.073, "C2": 0.0, "T_ref": 298 * u.K},
    }
    if ternary:
        sips["Hydrogen"] = {"A1": 5.013 * u.mol / u.kg, "A2": 0.0, "B1": 1.034e-4 / u.bar, "B2": 9.453 * u.kJ / u.mol,
                            "C1": 1.006, "C2": 0.0, "T_ref": 298 * u.K}
    zeolite.set_isotherm("extended_sips_2", sips, components=adsorbing)
    zeolite.set_isotherm("inert", components=["Helium"])
    # LDF coefficients calibrated on A3 (Table 4: 0.06 and 0.2 1/s).
    zeolite.set_rate_model("CarbonDioxide", "ldf", "constant", k=LDF_CO2)
    zeolite.set_rate_model("Methane", "ldf", "constant", k=LDF_CH4)
    if ternary:
        zeolite.set_rate_model("Hydrogen", "ldf", "constant", k=1.0 / u.s)     # Table 4 (b): assumed
    zeolite.set_rate_model("Helium", "ldf", "constant", k=0.5 / u.s)
    # Loading-averaged isosteric heats (Sect. 4.1).
    zeolite.heat_of_adsorption["CarbonDioxide"] = 36.3 * u.kJ / u.mol
    zeolite.heat_of_adsorption["Methane"] = 18.8 * u.kJ / u.mol
    if ternary:
        zeolite.heat_of_adsorption["Hydrogen"] = 10.7 * u.kJ / u.mol
    zeolite.heat_of_adsorption["Helium"] = 0 * u.kJ / u.mol

    # ---- column: 1.2 m x 2.5 cm, thermostatted at the experiment's temperature -----------
    sim.bed.set_shape("cylinder", diameter=2.5 * u.cm)
    sim.bed.ambient_temperature = T
    sim.bed.set_thermal_equilibrium("lte")
    sim.bed.set_particle_heat_transfer("constant", particle_htc=100 * u.W / (u.m**2 * u.K))
    for void in (sim.bed.bottom_void, sim.bed.top_void):       # no dead volume at the column ends
        void.enabled = False
        void.volume = 1e-4 * u.m**3
        void.resistance_coefficient = 0 * u.Pa * u.s / u.m
    # The paper does not model wall conduction, so the wall conductivity is 0.
    sim.wall.set_material("ss316l", density=8000 * u.kg / u.m**3, heat_capacity=500 * u.J / (u.kg * u.K),
                          thermal_conductivity=0 * u.W / (u.m * u.K))
    sim.wall.thickness = 2.5 * u.mm
    sim.wall.set_boundary("zero_flux")

    layer = sim.add_layer("13X", composition={"13X": 1.0}, momentum="kozeny_carman",
                          length=1.2 * u.m, cells=160, bed_void_fraction=BED_VOID_FRACTION, particle_diameter=1.8 * u.mm)
    layer.set_dispersion("constant", dispersion_coefficient=0 * u.m**2 / u.s)         # no axial dispersion
    layer.set_heat_dispersion("constant", heat_dispersion_coefficient=0 * u.W / (u.m * u.K))
    layer.set_internal_heat_transfer("dixon")                                  # Dixon (1996) correlation
    layer.set_external_heat_transfer("constant", external_htc=EXTERNAL_HTC)      # calibrated on A3

    # ---- sampling line to the MS: an inert, nearly empty, adiabatic layer above the bed ----
    # It carries the outlet gas as plug flow over its 40 cells (the engine's void region is
    # one well-mixed cell, which smears the CH4 front). Gas volume calibrated on A3.
    line_solid = sim.add_adsorbent("Sampling line", solid_density=1000 * u.kg / u.m**3,
                                   heat_capacity=500 * u.J / (u.kg * u.K), thermal_conductivity=0.15 * u.W / (u.m * u.K),
                                   particle_void_fraction=0.0)
    line_solid.set_isotherm("inert", components=components)
    for component in components:
        line_solid.set_rate_model(component, "ldf", "constant", k=0.5 / u.s)
        line_solid.heat_of_adsorption[component] = 0 * u.kJ / u.mol
    line = sim.add_layer("Sampling line", composition={"Sampling line": 1.0}, momentum="kozeny_carman",
                         length=LINE_VOLUME / (LINE_VOID_FRACTION * COLUMN_AREA), cells=40,
                         bed_void_fraction=LINE_VOID_FRACTION, particle_diameter=10 * u.mm)
    line.set_dispersion("constant", dispersion_coefficient=0 * u.m**2 / u.s)
    line.set_heat_dispersion("constant", heat_dispersion_coefficient=0 * u.W / (u.m * u.K))
    line.set_internal_heat_transfer("constant", internal_htc=0 * u.W / (u.m**2 * u.K))
    line.set_external_heat_transfer("constant", external_htc=0 * u.W / (u.m**2 * u.K))

    # ---- initial state: the column full of helium at the experiment's pressure -----------
    sim.initial_conditions.pressure = P
    sim.initial_conditions.temperature = T
    sim.initial_conditions.wall_temperature = T
    sim.initial_conditions.velocity = 0 * u.m / u.s
    for component in components:
        sim.initial_conditions.y[component] = 1.0 if component == "Helium" else 0.0

    # ---- one feed step ----------------------------------------------------------------------
    feed_y = ({"CarbonDioxide": 0.20, "Methane": 0.05, "Hydrogen": 0.75, "Helium": 0.0} if ternary
              else {"CarbonDioxide": 0.5, "Methane": 0.5, "Helium": 0.0})
    sim.cycle.number_of_cycles = 1
    sim.cycle.step_time_mode = "time"
    feed = sim.add_step("Feed", time=exp["time"] * u.s)
    feed.bottom.set_streams("inlet", [
        sk.Stream("Ternary Feed" if ternary else "Binary Feed",
                  flow=sk.volumetric_flow(exp["flow"] * 1e-6 * exp["flow_factor"] * u.m**3 / u.s), y=feed_y, temperature=T),
    ])
    feed.top.set_pressure("outlet", P, y={component: (1.0 if component == "Helium" else 0.0) for component in components},
                          temperature=293 * u.K)

    # ---- solver (SUNDIALS IDA); a single step, so no cyclic-steady-state check --------------
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
    sim.settings.advection_scheme = "weno"
    sim.settings.max_cycles_to_keep = 1
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


CASES = list(EXPERIMENTS)
DECKS = {case: f"Streb2021_{case}" for case in EXPERIMENTS}   # the shipped deck of each case
Z_CM = [10, 35, 60, 85, 110]                                # thermocouple positions [cm]


def main(argv=None) -> int:
    args = c.parser(__doc__.splitlines()[0], list(CASES), default=["A1", "A3"]).parse_args(argv)
    plt = c.pyplot(args)
    rows = []
    for name in c.selected(args):
        sim = build(name)
        if args.cycles is not None:
            sim.cycle.number_of_cycles = int(args.cycles)
        c.run(sim, f"Streb 2021 {name}", args.quiet)
        fig, axes = (plt.subplots(1, 2, figsize=(10, 4)) if plt else (None, None))
        d = c.DATA / "streb_2021" / name
        row = [name]
        for comp, label in (("CarbonDioxide", "CO2"), ("Methane", "CH4")):
            y = sim.series("mole_fraction", comp, location="top", cycle=0)
            t50 = c.crossing_time(y.time, y.values, 0.25)
            err = "-"
            f = d / f"{name}_mole_fraction_{label}.csv"
            if f.is_file():
                t_e, y_e = c.load_csv(f)
                err = c.fmt(c.nrmse(y.time, y.values, t_e, y_e, 0.5))
                if axes is not None:
                    axes[0].plot(t_e, y_e, "o", ms=3, label=f"{label} experiment")
            if axes is not None:
                axes[0].plot(y.time, y.values, "-", label=f"{label} Skarstrom")
            row += [c.fmt(t50), err]
        sq, n = 0.0, 0
        for z in Z_CM:
            T = sim.series("temperature", z=z / 100.0, cycle=0)
            f = d / f"{name}_temperature_z={z}[cm].csv"
            if f.is_file():
                t_e, T_e = c.load_csv(f)
                sq += float(np.sum((np.interp(t_e, T.time, T.values) - T_e) ** 2))
                n += len(t_e)
                if axes is not None:
                    axes[1].plot(t_e, T_e, "o", ms=2)
            if axes is not None:
                axes[1].plot(T.time, T.values, "-", label=f"z = {z} cm")
        row.append(c.fmt(math.sqrt(sq / n)) if n else "-")
        rows.append(row)
        if fig is not None:
            axes[0].set(xlabel="time [s]", ylabel="outlet mole fraction [-]")
            axes[1].set(xlabel="time [s]", ylabel="temperature [K] (lines: Skarstrom, markers: experiment)")
            for a in axes:
                a.legend(fontsize=7)
            c.save(fig, f"streb_2021_{name}.png")
    print()
    c.table(rows, ["case", "t50 CO2 [s]", "NRMSE y_CO2", "t50 CH4 [s]", "NRMSE y_CH4", "RMSE T, 5 TCs [K]"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
