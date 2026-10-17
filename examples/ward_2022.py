"""Ward & Pini (2022): CO2 breakthrough on activated carbon at 50, 100 and 200 ccm.

Builds each breakthrough experiment from scratch in code (pure CO2 fed into a
bed initially full of air): the fluid, the activated carbon and its isotherm
and rate model, the column and its heat transfer, the single feed step and the
solver settings. Every value is set explicitly, with its unit, so this file is
also a template for a breakthrough simulation of your own.

Each case is run and the outlet CO2 mole fraction and the thermocouple at
z = 0.0715 m are compared with the digitised experiments (core/data/ward_2022/).
The breakthrough time (outlet y_CO2 = 0.5) is set against the 95 % credibility
interval of the paper's own model (Table 10).

    python examples/ward_2022.py                 # all three flow rates
    python examples/ward_2022.py --case 100      # the fastest single case

``build("100")`` returns the configured simulation without running it
(python/tests/test_examples_build.py checks it against examples/decks/).

Reference: A. Ward, R. Pini, Adsorption 28 (2022) 161-183.
doi:10.1007/s10450-022-00361-z
"""

from __future__ import annotations

import sys

import _common as c
import skarstrom as sk

# The three experiments differ only in the CO2 feed flow and in the bed-to-wall
# heat transfer coefficient fitted for each flow (core/test/papers/ward_2022/ward_2022.cpp).
FLOW_RATES = {
    "50":  dict(feed_flow="50 ccm",  internal_htc="77.5 W/(m^2*K)"),
    "100": dict(feed_flow="100 ccm", internal_htc="57.9 W/(m^2*K)"),
    "200": dict(feed_flow="200 ccm", internal_htc="45.2 W/(m^2*K)"),
}
DECKS = {q: f"Ward2022_Q{q}ccm" for q in FLOW_RATES}   # the shipped deck of each case
Z_TC = 0.0715  # thermocouple position [m]

SPEED = {"50": "low", "100": "medium", "200": "high"}


def add_activated_carbon(sim: sk.Simulation) -> sk.Adsorbent:
    """The activated carbon of the lab-scale rig (Tables 1, 2 and 6 of the paper)."""
    ads = sim.add_adsorbent(
        "Activated Carbon",
        solid_density="799 kg/m3",                     # particle density
        heat_capacity="1040 J/(kg*K)",
        thermal_conductivity="0.63 W/(m*K)",
        particle_void_fraction=0.0,
    )
    # Dual-site Langmuir in concentration form for CO2; air does not adsorb.
    ads.set_isotherm("dual_site_langmuir_c", {
        "CarbonDioxide": {"m1": "0.418 mol/kg", "b1": "2.05e-6 m3/mol", "q1": "29.98 kJ/mol",
                          "m2": "6.04 mol/kg",  "b2": "3.32e-6 m3/mol", "q2": "20.367 kJ/mol"},
    }, components=["CarbonDioxide"])
    ads.set_isotherm("inert", components=["Air"])
    # Linear driving force; the CO2 coefficient from the Haynes (1973) macropore correlation.
    ads.set_rate_model("CarbonDioxide", "ldf", "haynes", tortuosity=3.0, macropore_void_fraction=0.59)
    ads.set_rate_model("Air", "ldf", "constant", k="0 1/s")
    ads.heat_of_adsorption["CarbonDioxide"] = "29.98 kJ/mol"
    ads.heat_of_adsorption["Air"] = "0 J/mol"
    return ads


def build(case: str) -> sk.Simulation:
    """The breakthrough experiment at one feed flow ("50", "100" or "200" ccm), ready to run."""
    cfg = FLOW_RATES[case]
    sim = sk.Simulation.new(
        title=f"Ward 2022 - {case} ccm Flow Rate",
        description=f"CO2 breakthrough on activated carbon (Ward & Pini 2022), {SPEED[case]} flow rate "
                    f"({case} ccm): pure CO2 fed into a bed initially full of air.",
        author="Skarstrom examples",
    )

    # ---- fluid ---------------------------------------------------------------------
    sim.add_component("CarbonDioxide")
    sim.add_component("Air")
    sim.fluid.adsorbed_heat_capacity_mode = "constant"
    sim.fluid.heat_of_adsorption_convention = "isosteric"
    sim.fluid.adsorbed_heat_capacity["CarbonDioxide"] = "36.96 J/(mol*K)"
    sim.fluid.adsorbed_heat_capacity["Air"] = "0 J/(mol*K)"

    # ---- adsorbent -------------------------------------------------------------------
    add_activated_carbon(sim)

    # ---- column: 0.286 m x 18.9 mm tube ---------------------------------------------------
    sim.bed.set_shape("cylinder", diameter="18.9 mm")
    sim.bed.ambient_temperature = "291.4 K"                    # Table 1
    sim.bed.set_thermal_equilibrium("lte")                     # gas and solid at one temperature
    sim.bed.set_particle_heat_transfer("constant", particle_htc="100 W/(m^2*K)")
    for void in (sim.bed.bottom_void, sim.bed.top_void):       # no dead volume at the column ends
        void.enabled = False
        void.volume = "1e-4 m3"
        void.resistance_coefficient = "0 Pa*s/m"
    # The tube's properties are set explicitly (the deck keeps the SS316L material label).
    sim.wall.set_material("ss316l", density="2700 kg/m3", heat_capacity="902 J/(kg*K)",
                          thermal_conductivity="205 W/(m*K)")
    sim.wall.thickness = "3.25 mm"
    sim.wall.set_boundary("zero_flux")

    layer = sim.add_layer("Activated Carbon", composition={"Activated Carbon": 1.0}, momentum="kozeny_carman",
                          length="0.286 m", cells=60, bed_void_fraction=0.456, particle_diameter="3 mm")
    layer.set_dispersion("ruthven", molecular_diffusion_coefficient="1.6e-5 m2/s")
    layer.set_heat_dispersion("yagi_kunii")                    # axial conductivity from Yagi & Kunii
    layer.set_internal_heat_transfer("constant", internal_htc=cfg["internal_htc"])
    # Wall to ambient by natural convection, Holman correlation Nu = c Re^m Pr^n.
    layer.set_external_heat_transfer("holman", holman_c=0.027, holman_m=0.8, holman_n=1 / 3,
                                     external_velocity="0.5 m/s")

    # ---- initial state: the bed full of air at 1 bar and 291.4 K --------------------------
    sim.initial_conditions.pressure = "1 bar"
    sim.initial_conditions.temperature = "291.4 K"
    sim.initial_conditions.wall_temperature = "291.4 K"
    sim.initial_conditions.velocity = "0 m/s"
    sim.initial_conditions.y["Air"] = 1.0
    sim.initial_conditions.y["CarbonDioxide"] = 0.0

    # ---- one feed step: pure CO2 in at the bottom, the top open at 1 bar -------------------
    sim.cycle.number_of_cycles = 1
    sim.cycle.step_time_mode = "time"
    feed = sim.add_step("Feed", time="3000 s")
    feed.bottom.set_streams("inlet", [
        sk.Stream("CO2 Inlet", flow=sk.volumetric_flow(cfg["feed_flow"]), y={"Air": 0.0, "CarbonDioxide": 1.0},
                  temperature="293 K"),
    ])
    feed.top.set_pressure("outlet", "1 bar", y={"CarbonDioxide": 0.0, "Air": 1.0}, temperature="293 K")

    # ---- solver (SUNDIALS IDA), cyclic steady state, energy balance ----------------------
    sim.settings.init_step = "1e-10 s"
    sim.settings.min_step = "1e-15 s"
    sim.settings.max_step = "0 s"                              # no limit on the integrator step
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
    sim.css.method = "none"                                    # a single breakthrough, no cycles
    sim.css.min_cycles = 5
    sim.css.molar_absolute_tolerance = 1e-4
    sim.css.molar_relative_tolerance = 1e-4
    sim.css.mass_balance_tolerance = 0.005
    sim.energy_balance.energy_balance_mode = "non_isothermal"
    sim.energy_balance.isothermal_use_initial_temperature = True
    sim.energy_balance.isothermal_temperature = "298.15 K"
    return sim


def main(argv=None) -> int:
    args = c.parser(__doc__.splitlines()[0], list(FLOW_RATES)).parse_args(argv)
    plt = c.pyplot(args)
    fig, axes = (plt.subplots(1, 2, figsize=(10, 4)) if plt else (None, None))
    rows = []
    for q in c.selected(args):
        sim = build(q)
        if args.cycles is not None:
            sim.cycle.number_of_cycles = int(args.cycles)
        c.run(sim, f"Ward 2022, {q} ccm", args.quiet)
        y = sim.series("mole_fraction", "CarbonDioxide", location="top", cycle=0)
        T = sim.series("temperature", z=Z_TC, cycle=0)
        t_br = c.crossing_time(y.time, y.values, 0.5)
        lo = c.reference("ward_2022", f"Q{q}ccm", "breakthrough_time_CI_low")
        hi = c.reference("ward_2022", f"Q{q}ccm", "breakthrough_time_CI_high")
        row = [f"{q} ccm", c.fmt(t_br), "-" if lo is None else f"[{lo:g}, {hi:g}]", "-", "-"]
        data_dir = c.DATA / "ward_2022" / f"Qin={q}"
        if data_dir.is_dir():
            t_e, y_e = c.load_csv(data_dir / f"outlet_CO2_mole_fraction[-]_Qin={q}[ccm].csv")
            tT_e, T_e = c.load_csv(data_dir / f"temperature[K]_z=0.0715[m]_Qin={q}[ccm].csv")
            row[3] = c.fmt(c.nrmse(y.time, y.values, t_e, y_e))  # feed is pure CO2: already normalised
            row[4] = c.fmt(c.nrmse(T.time, T.values, tT_e, T_e))
            if axes is not None:
                axes[0].plot(t_e, y_e, "o", ms=3, label=f"{q} ccm, experiment")
                axes[1].plot(tT_e, T_e, "o", ms=3, label=f"{q} ccm, experiment")
        if axes is not None:
            axes[0].plot(y.time, y.values, "-", label=f"{q} ccm, Skarstrom")
            axes[1].plot(T.time, T.values, "-", label=f"{q} ccm, Skarstrom")
        rows.append(row)
    print()
    c.table(rows, ["case", "t_br (y=0.5) [s]", "paper model 95% CI [s]", "NRMSE y_CO2", "RMSE T [K]"])
    if not c.have_reference_data():
        print("(run from a repository checkout to compare with the digitised data)")
    if fig is not None:
        axes[0].set(xlabel="time [s]", ylabel="outlet CO2 mole fraction [-]")
        axes[1].set(xlabel="time [s]", ylabel=f"temperature at z = {Z_TC} m [K]")
        for a in axes:
            a.legend(fontsize=7)
        c.save(fig, "ward_2022.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
