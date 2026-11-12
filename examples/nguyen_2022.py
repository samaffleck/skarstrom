"""Nguyen et al. (2022): CO2 breakthrough on CALF-20 (Figs. 4 and 5), all six experiments.

Builds each case from scratch in code: the CO2/N2/He fluid, CALF-20 with the
paper's extended dual-site Langmuir isotherm (equal-energy-site N2), the
column (sized to hold 16.7 g of adsorbent), the single feed step and the solver
settings. Every value is set explicitly, with its unit.

Fig. 4 feeds 50, 70 or 100 % CO2 in helium at 50 cm3/min; Fig. 5 feeds 5, 15
or 30 % CO2 in nitrogen at 100 cm3/min; the column starts full of helium. The
flows, temperatures and pressures are those of Table A2.1 of the author's PhD
thesis. Each case is compared with the measured markers re-digitised from the
paper's figures (core/data/nguyen_2022/): the outlet CO2 mole fraction (NRMSE
normalised by the feed fraction, and the time at which it reaches half the feed
value), the thermocouple 1.82 cm from the outlet and, in Fig. 5, the outlet N2.

    python examples/nguyen_2022.py
    python examples/nguyen_2022.py --case Fig4_70pct

``build("Fig4_70pct")`` returns the configured simulation without running it
(python/tests/test_examples_build.py checks it against examples/decks/).

Reference: T.T.T. Nguyen, J.-B. Lin, G.K.H. Shimizu, A. Rajendran, Chemical
Engineering Journal 442 (2022) 136263. doi:10.1016/j.cej.2022.136263
"""

from __future__ import annotations

import math
import sys

import numpy as np

import _common as c
import skarstrom as sk
from skarstrom import ureg as u  # pint units: 3.2 * u.cm, 40.98 * u.kJ / u.mol

# The experiments (thesis Table A2.1): figure, feed CO2 [%], flow [cm3/min at feed
# conditions], temperature [C] (feed, initial and ambient), pressure [bar], deck name.
CASES = {
    "Fig4_100pct": dict(fig="Fig4", percent=100, flow=50, temperature=22.4, pressure=0.97, deck="Nguyen2022_Fig4_100pct",
                        what="pure CO2 (balance He)"),
    "Fig4_70pct": dict(fig="Fig4", percent=70, flow=50, temperature=22.3, pressure=0.97, deck="Nguyen2022_Fig4",
                       what="70 % CO2 in He"),
    "Fig4_50pct": dict(fig="Fig4", percent=50, flow=50, temperature=22.4, pressure=0.97, deck="Nguyen2022_Fig4_50pct",
                       what="50 % CO2 in He"),
    "Fig5_30pct": dict(fig="Fig5", percent=30, flow=100, temperature=21.8, pressure=0.98, deck="Nguyen2022_Fig5_30pct",
                       what="30 % CO2 in N2"),
    "Fig5_15pct": dict(fig="Fig5", percent=15, flow=100, temperature=20.3, pressure=0.97, deck="Nguyen2022_Fig5",
                       what="15 % CO2 in N2"),
    "Fig5_5pct": dict(fig="Fig5", percent=5, flow=100, temperature=23.4, pressure=0.97, deck="Nguyen2022_Fig5_5pct",
                      what="5 % CO2 in N2"),
}
DECKS = {name: case["deck"] for name, case in CASES.items()}
FEED_TIME = {"Fig4": 7100.0, "Fig5": 11800.0}   # [s] covers the figure

# Column and CALF-20 (Table 2). The bed length is chosen so that it holds exactly
# 16.7 g of CALF-20, with rho_b = rho_p (1 - eb) (0.0782 m; Table 2: 0.0786 m).
ADSORBENT_MASS = 16.7e-3        # [kg]
INNER_DIAMETER = 0.0282         # [m]
OUTER_DIAMETER = 0.0318         # [m]
PARTICLE_DENSITY = 570.0        # [kg/m3]
BED_VOID_FRACTION = 0.4
PARTICLE_DIAMETER = 1.5e-3      # [m] assumption: the pellets are 1.0-3.0 mm (Sect. 2)
MOLECULAR_DIFFUSIVITY = 1.496e-5  # [m2/s] Table 2
BED_LENGTH = ADSORBENT_MASS / (PARTICLE_DENSITY * (1 - BED_VOID_FRACTION) * (0.25 * math.pi * INNER_DIAMETER * INNER_DIAMETER))
THERMOCOUPLE_Z = BED_LENGTH - 0.0182   # [m] 1.82 cm from the outlet (Sect. 3.2)

# Constant gas heat capacities of Table 2 (per kg: CO2 841, N2 1040, He 5200 J/(kg K))
# per mole; the adsorbed phase takes the gas's ("assumed", Table 2).
CP_GAS = {"CarbonDioxide": 841.0 * 0.0440095, "Nitrogen": 1040.0 * 0.0280134, "Helium": 5200.0 * 0.004002602}


def build(case: str) -> sk.Simulation:
    """The simulation of one experiment (a key of CASES), ready to run."""
    k = CASES[case]
    fig_label = "Fig. 4" if k["fig"] == "Fig4" else "Fig. 5"
    T = (k["temperature"] + 273.15) * u.K
    P = k["pressure"] * 1e5 * u.Pa
    sim = sk.Simulation.new(
        title=f"Nguyen 2022 - {fig_label}, {k['what']}",
        description=f"CO2 breakthrough on CALF-20 (Nguyen et al. 2022, {fig_label}): {k['what']} fed at {k['flow']} "
                    "cm3/min into a column full of helium. The measured outlet composition and the thermocouple 1.82 cm "
                    "from the outlet, re-digitised from the paper's figure, are bound to the deck (Experimental Data).",
        author="Skarstrom examples")

    # ---- fluid: CO2 and N2 adsorb, helium is the inert initial gas --------------------
    for component in ("CarbonDioxide", "Nitrogen", "Helium"):
        sim.add_component(component)
    sim.fluid.adsorbed_heat_capacity_mode = "constant"
    sim.fluid.heat_of_adsorption_convention = "isosteric"
    sim.fluid.gas_heat_capacity_mode = "constant"
    for component, cp in CP_GAS.items():
        sim.fluid.gas_heat_capacity[component] = cp * u.J / (u.mol * u.K)
        sim.fluid.adsorbed_heat_capacity[component] = (0.0 if component == "Helium" else cp) * u.J / (u.mol * u.K)

    # ---- adsorbent: CALF-20 with 20 wt% polysulfone binder (Table 2) ----------------------
    # The paper's balances hold up gas in the bed voidage only: the particle voidage
    # (0.35) enters the LDF coefficient alone.
    calf20 = sim.add_adsorbent("CALF-20", solid_density=PARTICLE_DENSITY * u.kg / u.m**3,
                               heat_capacity=1371 * u.J / (u.kg * u.K), thermal_conductivity=0.63 * u.W / (u.m * u.K),
                               particle_void_fraction=0.0)
    # Extended dual-site Langmuir, concentration form; N2 on equal-energy sites (Table 3).
    calf20.set_isotherm("extended_dual_site_langmuir_c", {
        "CarbonDioxide": {"m1": 2.387 * u.mol / u.kg, "b1": 5.519e-7 * u.m**3 / u.mol, "q1": 35.06 * u.kJ / u.mol,
                          "m2": 3.271 * u.mol / u.kg, "b2": 5.187e-8 * u.m**3 / u.mol, "q2": 28.95 * u.kJ / u.mol},
        "Nitrogen": {"m1": 2.387 * u.mol / u.kg, "b1": 8.139e-7 * u.m**3 / u.mol, "q1": 17.96 * u.kJ / u.mol,
                     "m2": 3.271 * u.mol / u.kg, "b2": 8.139e-7 * u.m**3 / u.mol, "q2": 17.96 * u.kJ / u.mol},
    }, components=["CarbonDioxide", "Nitrogen"])
    calf20.set_isotherm("inert", components=["Helium"])
    # LDF with the macropore coefficient and the capacity factor (Eq. 5).
    for component in ("CarbonDioxide", "Nitrogen"):
        calf20.set_rate_model(component, "ldf", "macropore_capacity", tortuosity=3.0, macropore_void_fraction=0.35,
                              molecular_diffusion_coefficient=MOLECULAR_DIFFUSIVITY * u.m**2 / u.s)
    calf20.set_rate_model("Helium", "ldf", "constant", k=0.5 / u.s)
    # Isosteric heats of Fig. 2d (about 37.5 and 15 kJ/mol).
    calf20.heat_of_adsorption["CarbonDioxide"] = 38.0 * u.kJ / u.mol
    calf20.heat_of_adsorption["Nitrogen"] = 15.0 * u.kJ / u.mol
    calf20.heat_of_adsorption["Helium"] = 0 * u.kJ / u.mol

    # ---- column: stainless steel, 2.82 cm bore -----------------------------------------
    sim.bed.set_shape("cylinder", diameter=INNER_DIAMETER * u.m)
    sim.bed.ambient_temperature = T
    sim.bed.set_thermal_equilibrium("lte")
    sim.bed.set_particle_heat_transfer("constant", particle_htc=100 * u.W / (u.m**2 * u.K))
    for void in (sim.bed.bottom_void, sim.bed.top_void):       # no dead volume at the column ends
        void.enabled = False
        void.volume = 1e-4 * u.m**3
        void.resistance_coefficient = 0 * u.Pa * u.s / u.m
    sim.wall.set_material("ss316l", density=7800 * u.kg / u.m**3, heat_capacity=502 * u.J / (u.kg * u.K),
                          thermal_conductivity=16 * u.W / (u.m * u.K))
    sim.wall.thickness = 0.5 * (OUTER_DIAMETER - INNER_DIAMETER) * u.m
    sim.wall.set_boundary("zero_flux")                       # the paper does not state the wall ends

    layer = sim.add_layer("CALF-20", composition={"CALF-20": 1.0}, momentum="blake_kozeny",   # Eq. 7
                          length=BED_LENGTH * u.m, cells=30, bed_void_fraction=BED_VOID_FRACTION,
                          particle_diameter=PARTICLE_DIAMETER * u.m)
    layer.set_dispersion("ruthven", molecular_diffusion_coefficient=MOLECULAR_DIFFUSIVITY * u.m**2 / u.s)   # Eq. 6
    layer.set_heat_dispersion("constant", heat_dispersion_coefficient=0.4 * u.W / (u.m * u.K))
    layer.set_internal_heat_transfer("constant", internal_htc=6.963 * u.W / (u.m**2 * u.K))
    layer.set_external_heat_transfer("constant", external_htc=14 * u.W / (u.m**2 * u.K))

    # ---- initial state: the column full of helium -----------------------------------------
    sim.initial_conditions.pressure = P
    sim.initial_conditions.temperature = T
    sim.initial_conditions.wall_temperature = T
    sim.initial_conditions.velocity = 0 * u.m / u.s
    sim.initial_conditions.y["Helium"] = 1.0
    sim.initial_conditions.y["Nitrogen"] = 0.0
    sim.initial_conditions.y["CarbonDioxide"] = 0.0

    # ---- one feed step: CO2 in He (Fig. 4) or in N2 (Fig. 5) -------------------------------
    y = k["percent"] / 100.0
    feed_y = {"CarbonDioxide": y, "Nitrogen": 0.0, "Helium": 0.0}
    feed_y["Helium" if k["fig"] == "Fig4" else "Nitrogen"] = 1.0 - y
    sim.cycle.number_of_cycles = 1
    sim.cycle.step_time_mode = "time"
    feed = sim.add_step("Feed", time=FEED_TIME[k["fig"]] * u.s)
    feed.bottom.set_streams("inlet", [
        sk.Stream("Feed", flow=sk.volumetric_flow(k["flow"] * 1e-6 / 60.0 * u.m**3 / u.s), y=feed_y, temperature=T),
    ])
    feed.top.set_pressure("outlet", P, y={"CarbonDioxide": 0.0, "Nitrogen": 0.0, "Helium": 1.0}, temperature=293 * u.K)

    # ---- solver (SUNDIALS IDA); a single step, so no cyclic-steady-state check --------------
    sim.settings.init_step = 1e-10 * u.s
    sim.settings.min_step = 1e-15 * u.s
    sim.settings.max_step = 0 * u.s
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
    sim.energy_balance.energy_balance_mode = "non_isothermal"
    return sim


def main(argv=None) -> int:
    args = c.parser(__doc__.splitlines()[0], list(CASES)).parse_args(argv)
    plt = c.pyplot(args)
    rows = []
    for name in c.selected(args):
        k = CASES[name]
        y_feed = k["percent"] / 100.0
        sim = build(name)
        if args.cycles is not None:
            sim.cycle.number_of_cycles = int(args.cycles)
        c.run(sim, f"Nguyen 2022 {name}", args.quiet)
        y = sim.series("mole_fraction", "CarbonDioxide", location="top", cycle=0)
        T = sim.series("temperature", z=THERMOCOUPLE_Z, cycle=0)
        row = [name, c.fmt(c.crossing_time(y.time, y.values, 0.5 * y_feed)), "-", "-", "-"]
        d = c.DATA / "nguyen_2022" / k["fig"]
        fig, axes = (plt.subplots(1, 2, figsize=(10, 4)) if plt else (None, None))
        f_y, f_T = d / f"Y_CO2_{k['percent']}%.csv", d / f"T_CO2_{k['percent']}%.csv"
        if f_y.is_file():
            t_e, y_e = c.load_csv(f_y)
            row[2] = c.fmt(c.crossing_time(t_e, y_e, 0.5 * y_feed))
            row[3] = c.fmt(c.nrmse(y.time, y.values, t_e, y_e, y_feed))
            if axes is not None:
                axes[0].plot(t_e, y_e, "o", ms=3, label="CO2 measured")
        if f_T.is_file():
            tT_e, T_e = c.load_csv(f_T)
            row[4] = c.fmt(math.sqrt(float(np.mean((np.interp(tT_e, T.time, T.values) - T_e) ** 2))))
            if axes is not None:
                axes[1].plot(tT_e, T_e, "o", ms=3, label="measured")
        rows.append(row)
        if fig is not None:
            axes[0].plot(y.time, y.values, "-", label="CO2 Skarstrom")
            axes[1].plot(T.time, T.values, "-", label="Skarstrom")
            axes[0].set(xlabel="time [s]", ylabel="outlet CO2 mole fraction [-]", title=f"Nguyen 2022 {name}")
            axes[1].set(xlabel="time [s]", ylabel="temperature 1.82 cm from the outlet [K]")
            for a in axes:
                a.legend(fontsize=7)
            c.save(fig, f"nguyen_2022_{name}.png")
    print()
    c.table(rows, ["case", "t50 simulated [s]", "t50 measured [s]", "NRMSE y_CO2", "RMSE T [K]"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
