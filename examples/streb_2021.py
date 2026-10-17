"""Streb & Mazzotti (2021): CO2/CH4 breakthrough on zeolite 13X with thermocouples.

Builds each case from scratch in code: the CO2/CH4/He fluid, zeolite 13X with
the paper's extended-Sips competitive isotherm and LDF coefficients, the 1.2 m
column, the single feed step and the solver settings. Every value is set
explicitly, with its unit.

Experiments A1 (1 bar) and A3 (10 bar, the reference experiment) feed a
50:50 CO2/CH4 mixture into a helium-filled column. The outlet mole fractions
and the five centre-line thermocouples (z = 10, 35, 60, 85, 110 cm) are
compared with the digitised experiments (core/data/streb_2021/).

    python examples/streb_2021.py
    python examples/streb_2021.py --case A3

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

# The two experiments differ only in the column pressure (Table 2); both run at
# 25 C with the same 20 cm3/s feed (at column conditions) for 4100 s.
EXPERIMENTS = {
    "A1": dict(pressure="1 bar", title="Streb 2021 - Experiment A1",
               description="Binary CO2/CH4 breakthrough on zeolite 13X at 1 bar (Streb & Mazzotti 2021, "
                           "experiment A1); the bed is initially full of helium."),
    "A3": dict(pressure="10 bar", title="Streb 2021 - Experiment A3",
               description="Binary CO2/CH4 breakthrough on zeolite 13X at 10 bar, the reference experiment "
                           "(Streb & Mazzotti 2021, experiment A3); the bed is initially full of helium."),
}

# Densities of the 13X beads and the packed bed (Table 1) give both voidages.
SKELETAL_DENSITY = 2359.0       # [kg/m3]
PARTICLE_DENSITY = 1085.0       # [kg/m3]
BED_DENSITY = 595.0             # [kg/m3]
PARTICLE_VOID_FRACTION = 1 - PARTICLE_DENSITY / SKELETAL_DENSITY
BED_VOID_FRACTION = 1 - BED_DENSITY / PARTICLE_DENSITY


def build(case: str) -> sk.Simulation:
    """The simulation of one experiment ("A1" or "A3"), ready to run."""
    exp = EXPERIMENTS[case]
    sim = sk.Simulation.new(title=exp["title"], description=exp["description"], author="Skarstrom examples")

    # ---- fluid: CO2 and CH4 adsorb, helium is the inert initial gas ----------------
    for component in ("CarbonDioxide", "Methane", "Helium"):
        sim.add_component(component)
    sim.fluid.adsorbed_heat_capacity_mode = "constant"
    sim.fluid.heat_of_adsorption_convention = "isosteric"
    for component in ("CarbonDioxide", "Methane", "Helium"):
        sim.fluid.adsorbed_heat_capacity[component] = "0 J/(mol*K)"

    # ---- adsorbent: zeolite 13X, extended-Sips competitive isotherm (Table 3) ---------
    zeolite = sim.add_adsorbent("13X", solid_density=f"{PARTICLE_DENSITY} kg/m3",
                                heat_capacity="1100 J/(kg*K)", thermal_conductivity="0.15 W/(m*K)",
                                particle_void_fraction=PARTICLE_VOID_FRACTION)
    # Extended Sips with temperature-dependent capacity and exponent (Table 3).
    zeolite.set_isotherm("extended_sips_2", {
        "CarbonDioxide": {"A1": "7.268 mol/kg", "A2": -0.6168, "B1": "1.129e-4 1/bar", "B2": "28.389 kJ/mol",
                          "C1": 0.425, "C2": 0.7238, "T_ref": "298 K"},
        "Methane":       {"A1": "4.473 mol/kg", "A2": -0.6569, "B1": "3.605e-4 1/bar", "B2": "15.561 kJ/mol",
                          "C1": 1.073, "C2": 0.0, "T_ref": "298 K"},
    }, components=["CarbonDioxide", "Methane"])
    zeolite.set_isotherm("inert", components=["Helium"])
    # LDF coefficients estimated on A3 (Table 4).
    zeolite.set_rate_model("CarbonDioxide", "ldf", "constant", k="0.06 1/s")
    zeolite.set_rate_model("Methane", "ldf", "constant", k="0.2 1/s")
    zeolite.set_rate_model("Helium", "ldf", "constant", k="0.5 1/s")
    # Loading-averaged isosteric heats (Sect. 4.1).
    zeolite.heat_of_adsorption["CarbonDioxide"] = "36.3 kJ/mol"
    zeolite.heat_of_adsorption["Methane"] = "18.8 kJ/mol"
    zeolite.heat_of_adsorption["Helium"] = "0 kJ/mol"

    # ---- column: 1.2 m x 2.5 cm, thermostatted at 25 C --------------------------------
    sim.bed.set_shape("cylinder", diameter="2.5 cm")
    sim.bed.ambient_temperature = "298 K"
    sim.bed.set_thermal_equilibrium("lte")
    sim.bed.set_particle_heat_transfer("constant", particle_htc="100 W/(m^2*K)")
    for void in (sim.bed.bottom_void, sim.bed.top_void):       # no dead volume at the column ends
        void.enabled = False
        void.volume = "1e-4 m3"
        void.resistance_coefficient = "0 Pa*s/m"
    # The paper does not model wall conduction, so the wall conductivity is 0.
    sim.wall.set_material("ss316l", density="8000 kg/m3", heat_capacity="500 J/(kg*K)",
                          thermal_conductivity="0 W/(m*K)")
    sim.wall.thickness = "2.5 mm"
    sim.wall.set_boundary("zero_flux")

    layer = sim.add_layer("13X", composition={"13X": 1.0}, momentum="kozeny_carman",
                          length="1.2 m", cells=60, bed_void_fraction=BED_VOID_FRACTION, particle_diameter="1.8 mm")
    layer.set_dispersion("constant", dispersion_coefficient="0 m2/s")         # no axial dispersion
    layer.set_heat_dispersion("constant", heat_dispersion_coefficient="0 W/(m*K)")
    layer.set_internal_heat_transfer("dixon")                                  # Dixon (1996) correlation
    layer.set_external_heat_transfer("constant", external_htc="220 W/(m^2*K)")

    # ---- initial state: the column full of helium at the experiment's pressure -----------
    sim.initial_conditions.pressure = exp["pressure"]
    sim.initial_conditions.temperature = "298 K"
    sim.initial_conditions.wall_temperature = "298 K"
    sim.initial_conditions.velocity = "0 m/s"
    sim.initial_conditions.y["Helium"] = 1.0
    sim.initial_conditions.y["Methane"] = 0.0
    sim.initial_conditions.y["CarbonDioxide"] = 0.0

    # ---- one feed step: 50:50 CO2/CH4 breakthrough ---------------------------------------
    sim.cycle.number_of_cycles = 1
    sim.cycle.step_time_mode = "time"
    feed = sim.add_step("Feed", time="4100 s")
    feed.bottom.set_streams("inlet", [
        sk.Stream("Binary Feed", flow=sk.volumetric_flow("20 cm^3/s"),
                  y={"Helium": 0.0, "Methane": 0.5, "CarbonDioxide": 0.5}, temperature="298 K"),
    ])
    feed.top.set_pressure("outlet", exp["pressure"], y={"CarbonDioxide": 0.0, "Methane": 0.0, "Helium": 1.0},
                          temperature="293 K")

    # ---- solver (SUNDIALS IDA); a single step, so no cyclic-steady-state check --------------
    sim.settings.init_step = "1e-10 s"
    sim.settings.min_step = "1e-15 s"
    sim.settings.max_step = "0 s"                 # 0: no limit on the integrator step
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
    sim.energy_balance.isothermal_temperature = "298.15 K"
    return sim


CASES = list(EXPERIMENTS)
DECKS = {"A1": "Streb2021_A1", "A3": "Streb2021_A3"}       # the shipped deck of each case
Z_CM = [10, 35, 60, 85, 110]                                # thermocouple positions [cm]


def main(argv=None) -> int:
    args = c.parser(__doc__.splitlines()[0], list(CASES)).parse_args(argv)
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
