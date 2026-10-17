"""Casas et al. (2012): CO2/H2 breakthrough on activated carbon AP3-60 at 15 bar, 25 C.

Builds each breakthrough experiment from scratch in code: CO2, H2 and the
helium the column is first filled with, the AP3-60 carbon with the competitive
Sips isotherm of Schell et al. (2012), the column with its heating shells and
fitted heat-transfer correlations, and the single feed step. Every value is
set explicitly, with its unit.

Runs the three feed compositions (50/50 of Fig. 4, 25/75 and 75/25 of Fig. 6)
and compares the H2 and CO2 breakthrough times (outlet mole fraction above
0.002, the paper's Sect. 4.4 definition) with the paper's simulated values
(Table 5).

    python examples/casas_2012.py
    python examples/casas_2012.py --case Fig4_50_50

``build("Fig4_50_50")`` returns the configured simulation without running it
(python/tests/test_examples_build.py checks it against examples/decks/).

Reference: N. Casas, J. Schell, R. Pini, M. Mazzotti, Adsorption 18 (2012)
143-161. doi:10.1007/s10450-012-9389-z
"""

from __future__ import annotations

import sys

import _common as c
import skarstrom as sk

# Feed CO2 mole fraction (the rest is H2) and the reference case name of each experiment.
EXPERIMENTS = {
    "Fig4_50_50": dict(y_co2=0.5, ref="Fig4_Reference_50_50_15bar_25C", figure="Fig. 4 (reference, 50/50, 15 bar, 25 C)"),
    "Fig6_25_75": dict(y_co2=0.25, ref="Fig6_25_75_15bar_25C", figure="Fig. 6 (25/75, 15 bar, 25 C)"),
    "Fig6_75_25": dict(y_co2=0.75, ref="Fig6_75_25_15bar_25C", figure="Fig. 6 (75/25, 15 bar, 25 C)"),
}
PRESSURE = "15 bar"          # total pressure, held at the outlet
TEMPERATURE = "298.15 K"     # feed, initial and ambient temperature (25 C)
COMPONENTS = ("CarbonDioxide", "Hydrogen", "Helium")


def add_ap3_60(sim: sk.Simulation) -> None:
    """Activated carbon AP3-60 (Casas 2012 Table 2) with the Schell (2012) competitive Sips isotherm."""
    particle_density = 850.0     # [kg/m3] (manufacturer)
    skeletal_density = 1970.0    # [kg/m3] (He pycnometry, 1.97 g/cm3)
    ads = sim.add_adsorbent(
        "AP3-60",
        solid_density=f"{particle_density} kg/m3",
        heat_capacity="1000 J/(kg*K)",
        thermal_conductivity="0.15 W/(m*K)",      # not tabulated; nominal (core/test/papers/casas_2012/README.md)
        particle_void_fraction=1.0 - particle_density / skeletal_density,
    )
    # Sips with temperature-dependent saturation capacity, affinity and exponent:
    # n_inf = a exp(b / RT), K = A exp(B / RT), c = c_ref + c_amp tanh(c_beta (T - c_Tref)).
    ads.set_isotherm("arrhenius_sips", {
        "CarbonDioxide": {"n_inf_a": "1.38286 mol/kg", "n_inf_b": "5625.93 J/mol", "k_A": "1.6814e-8 1/Pa",
                          "k_B": "9160.22 J/mol", "c_ref": 0.832365, "c_amp": 0.116571, "c_beta": "0.0928714 1/K",
                          "c_Tref": "330.129 K"},
        "Hydrogen":      {"n_inf_a": "6.66 mol/kg", "n_inf_b": "0 J/mol", "k_A": "6.6886e-10 1/Pa",
                          "k_B": "9954.73 J/mol", "c_ref": 0.96, "c_amp": 0.0, "c_beta": "0 1/K", "c_Tref": "298.15 K"},
    }, components=["CarbonDioxide", "Hydrogen"])
    ads.set_isotherm("inert", components=["Helium"])        # helium does not adsorb
    # Lumped LDF coefficients fitted in the paper (Table 3); helium keeps the default.
    ads.set_rate_model("CarbonDioxide", "ldf", "constant", k="0.15 1/s")
    ads.set_rate_model("Hydrogen", "ldf", "constant", k="1 1/s")
    ads.set_rate_model("Helium", "ldf", "constant", k="0.5 1/s")
    # Isosteric heats (Table 2, Clausius-Clapeyron).
    ads.heat_of_adsorption["CarbonDioxide"] = "26 kJ/mol"
    ads.heat_of_adsorption["Hydrogen"] = "9.8 kJ/mol"
    ads.heat_of_adsorption["Helium"] = "0 J/mol"


def build(case: str) -> sk.Simulation:
    """The simulation of one experiment ("Fig4_50_50", "Fig6_25_75" or "Fig6_75_25"), ready to run."""
    exp = EXPERIMENTS[case]
    y_co2 = exp["y_co2"]
    sim = sk.Simulation.new(
        title=f"Casas 2012 - {exp['figure']}",
        description=f"Binary CO2/H2 breakthrough on AP3-60 activated carbon (Casas et al. 2012): a "
                    f"{round(100 * y_co2)}/{round(100 * (1 - y_co2))} CO2/H2 feed at 15 bar and 25 C into a "
                    "column filled with helium.",
        author="Skarstrom examples",
    )

    # ---- fluid ------------------------------------------------------------------------------
    for component in COMPONENTS:
        sim.add_component(component)
    sim.fluid.adsorbed_heat_capacity_mode = "constant"
    sim.fluid.heat_of_adsorption_convention = "isosteric"
    for component in COMPONENTS:
        sim.fluid.adsorbed_heat_capacity[component] = "0 J/(mol*K)"

    # ---- adsorbent ---------------------------------------------------------------------------
    add_ap3_60(sim)

    # ---- column: 1.20 m x 2.5 cm (Table 2) -----------------------------------------------------
    inner_diameter = 0.025                                     # [m]
    sim.bed.set_shape("cylinder", diameter=f"{inner_diameter} m")
    sim.bed.ambient_temperature = TEMPERATURE
    sim.bed.set_thermal_equilibrium("lte")                     # gas and solid at one temperature
    sim.bed.set_particle_heat_transfer("constant", particle_htc="100 W/(m^2*K)")
    for void in (sim.bed.bottom_void, sim.bed.top_void):       # no dead volume at the column ends
        void.enabled = False
        void.volume = "1e-4 m3"
        void.resistance_coefficient = "0 Pa*s/m"
    # The wall and its copper/brass heating shells, lumped (Table 2): volumetric heat
    # capacity 8000 kg/m3 x 500 J/(kg K) = 4e6 J/(m3 K) out to R_o = 0.02 m; no axial conduction.
    sim.wall.set_material("ss316l", density="8000 kg/m3", heat_capacity="500 J/(kg*K)",
                          thermal_conductivity="0 W/(m*K)")
    sim.wall.thickness = f"{0.02 - 0.5 * inner_diameter!r} m"
    sim.wall.set_boundary("zero_flux")

    layer = sim.add_layer("AP3-60", composition={"AP3-60": 1.0}, momentum="kozeny_carman",
                          length="1.2 m", cells=60, bed_void_fraction=0.403, particle_diameter="3 mm")
    layer.set_dispersion("ruthven", molecular_diffusion_coefficient="1.6e-5 m2/s")
    layer.set_heat_dispersion("constant", heat_dispersion_coefficient="0 W/(m*K)")
    # Gas-to-wall coefficient: Leva (1947) with the coefficient and exponent fitted by
    # Casas (Table 3; Leva's own values are 0.813 and 0.90). Wall-to-ambient: fitted constant.
    layer.set_internal_heat_transfer("leva", leva_eta1=41.13, leva_eta2=0.32, leva_exp_coeff=6.0)
    layer.set_external_heat_transfer("constant", external_htc="5 W/(m^2*K)")

    # ---- initial state: the column pressurised with helium -----------------------------------
    sim.initial_conditions.pressure = PRESSURE
    sim.initial_conditions.temperature = TEMPERATURE
    sim.initial_conditions.wall_temperature = TEMPERATURE
    sim.initial_conditions.velocity = "0 m/s"
    sim.initial_conditions.y["Helium"] = 1.0
    sim.initial_conditions.y["Hydrogen"] = 0.0
    sim.initial_conditions.y["CarbonDioxide"] = 0.0

    # ---- one feed step: CO2/H2 at a constant 10 cm3/s, the outlet held at 15 bar ---------------
    sim.cycle.number_of_cycles = 1
    sim.cycle.step_time_mode = "time"
    feed = sim.add_step("Feed", time="1800 s")
    feed.bottom.set_streams("inlet", [
        sk.Stream("CO2/H2 Feed", flow=sk.volumetric_flow("10 cm^3/s"),
                  y={"Helium": 0.0, "Hydrogen": 1.0 - y_co2, "CarbonDioxide": y_co2}, temperature=TEMPERATURE),
    ])
    feed.top.set_pressure("outlet", PRESSURE, y={"CarbonDioxide": 0.0, "Hydrogen": 0.0, "Helium": 1.0},
                          temperature="293 K")

    # ---- solver (SUNDIALS IDA), cyclic steady state, energy balance ------------------------------
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
    sim.css.method = "none"                                    # a single breakthrough, no cycle
    sim.css.min_cycles = 5
    sim.css.molar_absolute_tolerance = 1e-4
    sim.css.molar_relative_tolerance = 1e-4
    sim.css.mass_balance_tolerance = 0.005
    sim.energy_balance.energy_balance_mode = "non_isothermal"
    sim.energy_balance.isothermal_use_initial_temperature = True
    sim.energy_balance.isothermal_temperature = "298.15 K"
    return sim


DECKS = {
    "Fig4_50_50": "Casas2012_Fig4_50_50_15bar_25C",
    "Fig6_25_75": "Casas2012_Fig6_25_75_15bar_25C",
    "Fig6_75_25": "Casas2012_Fig6_75_25_15bar_25C",
}
THRESHOLD = 0.002


def main(argv=None) -> int:
    args = c.parser(__doc__.splitlines()[0], list(EXPERIMENTS)).parse_args(argv)
    plt = c.pyplot(args)
    fig, ax = (plt.subplots(figsize=(6, 4)) if plt else (None, None))
    rows = []
    for name in c.selected(args):
        sim = build(name)
        if args.cycles is not None:
            sim.cycle.number_of_cycles = int(args.cycles)
        c.run(sim, f"Casas 2012 {name}", args.quiet)
        for comp, label in (("Hydrogen", "H2"), ("CarbonDioxide", "CO2")):
            s = sim.series("mole_fraction", comp, location="top", cycle=0)
            t = c.crossing_time(s.time, s.values, THRESHOLD)
            ref = c.reference("casas_2012", EXPERIMENTS[name]["ref"], f"{label}_breakthrough_time")
            err = "-" if ref is None else f"{100 * (t - ref) / ref:+.1f} %"
            rows.append([name, label, c.fmt(t), c.fmt(ref), err])
            if ax is not None:
                ax.plot(s.time, s.values, label=f"{name} {label}")
    print()
    c.table(rows, ["case", "component", "t_br simulated [s]", "t_br paper model [s]", "difference"])
    if fig is not None:
        ax.set(xlabel="time [s]", ylabel="outlet mole fraction [-]")
        ax.legend(fontsize=7)
        c.save(fig, "casas_2012.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
