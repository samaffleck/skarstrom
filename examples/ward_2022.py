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
from skarstrom import ureg as u  # pint units: 3.2 * u.cm, 40.98 * u.kJ / u.mol

# The three experiments differ only in the CO2 feed flow, in the bed-to-wall
# heat transfer coefficient fitted for each flow and in the paper's dispersion
# coefficient (core/test/papers/ward_2022/ward_2022.cpp).
FLOW_RATES = {
    "50":  dict(feed_flow=50 * u.cm**3 / u.min,  internal_htc=77.5 * u.W / (u.m**2 * u.K)),
    "100": dict(feed_flow=100 * u.cm**3 / u.min, internal_htc=57.9 * u.W / (u.m**2 * u.K)),
    "200": dict(feed_flow=200 * u.cm**3 / u.min, internal_htc=45.2 * u.W / (u.m**2 * u.K)),
}
DECKS = {q: f"Ward2022_Q{q}ccm" for q in FLOW_RATES}   # the shipped deck of each case
Z_TC = 0.0715  # thermocouple position [m]

SPEED = {"50": "low", "100": "medium", "200": "high"}


def paper_dispersion(feed_flow) -> float:
    """Ward & Pini's dispersion coefficient [m2/s] for one feed flow: 4 x (0.7 D_m + v_F r_p)."""
    import math
    Dm, rp, Di, eb = 0.605e-5, 1.5e-3, 0.0189, 0.456
    q = feed_flow.to("m**3/s").magnitude
    vF = q / (math.pi / 4 * Di**2) / eb
    return 4.0 * (0.7 * Dm + vF * rp)


def add_activated_carbon(sim: sk.Simulation) -> sk.Adsorbent:
    """The activated carbon of the lab-scale rig (Tables 1, 2 and 6 of the paper)."""
    ads = sim.add_adsorbent(
        "Activated Carbon",
        solid_density=799 * u.kg / u.m**3,                     # particle density
        heat_capacity=1040 * u.J / (u.kg * u.K),
        thermal_conductivity=0.63 * u.W / (u.m * u.K),
        particle_void_fraction=0.0,
    )
    # Dual-site Langmuir in concentration form for CO2; helium does not adsorb.
    ads.set_isotherm("dual_site_langmuir_c", {
        "CarbonDioxide": {"m1": 0.418 * u.mol / u.kg, "b1": 2.05e-6 * u.m**3 / u.mol, "q1": 29.98 * u.kJ / u.mol,
                          "m2": 6.04 * u.mol / u.kg,  "b2": 3.32e-6 * u.m**3 / u.mol, "q2": 20.367 * u.kJ / u.mol},
    }, components=["CarbonDioxide"])
    ads.set_isotherm("inert", components=["Helium"])
    # Linear driving force; the CO2 coefficient from the Haynes (1973) macropore correlation.
    ads.set_rate_model("CarbonDioxide", "ldf", "haynes", tortuosity=3.0, macropore_void_fraction=0.59)
    ads.set_rate_model("Helium", "ldf", "constant", k=0 / u.s)
    # The paper's energy balance releases -dU_1 + R T per mole adsorbed (Eq. A4), at T_F = 293 K.
    ads.heat_of_adsorption["CarbonDioxide"] = (29980.0 + 8.314462618 * 293.0) * u.J / u.mol
    ads.heat_of_adsorption["Helium"] = 0 * u.J / u.mol
    return ads


def build(case: str) -> sk.Simulation:
    """The breakthrough experiment at one feed flow ("50", "100" or "200" ccm), ready to run."""
    cfg = FLOW_RATES[case]
    sim = sk.Simulation.new(
        title=f"Ward 2022 - {case} ccm Flow Rate",
        description=f"CO2 breakthrough on activated carbon (Ward & Pini 2022), {SPEED[case]} flow rate "
                    f"({case} ccm): pure CO2 fed into a bed initially full of helium.",
        author="Skarstrom examples",
    )

    # ---- fluid ---------------------------------------------------------------------
    sim.add_component("CarbonDioxide")
    sim.add_component("Helium")                                # the CO2/He experiment (Sect. 2)
    sim.fluid.adsorbed_heat_capacity_mode = "constant"
    sim.fluid.heat_of_adsorption_convention = "isosteric"
    sim.fluid.adsorbed_heat_capacity["CarbonDioxide"] = 36.96 * u.J / (u.mol * u.K)
    sim.fluid.adsorbed_heat_capacity["Helium"] = 0 * u.J / (u.mol * u.K)

    # ---- adsorbent -------------------------------------------------------------------
    add_activated_carbon(sim)

    # ---- column: 0.286 m x 18.9 mm tube ---------------------------------------------------
    sim.bed.set_shape("cylinder", diameter=18.9 * u.mm)
    sim.bed.ambient_temperature = 291.4 * u.K                    # Table 1
    sim.bed.set_thermal_equilibrium("lte")                     # gas and solid at one temperature
    sim.bed.set_particle_heat_transfer("constant", particle_htc=100 * u.W / (u.m**2 * u.K))
    for void in (sim.bed.bottom_void, sim.bed.top_void):       # no dead volume at the column ends
        void.enabled = False
        void.volume = 1e-4 * u.m**3
        void.resistance_coefficient = 0 * u.Pa * u.s / u.m
    # Aluminium tube: pure aluminium's properties (not a preset alloy), so a custom material.
    sim.wall.set_material("custom", density=2700 * u.kg / u.m**3, heat_capacity=902 * u.J / (u.kg * u.K),
                          thermal_conductivity=205 * u.W / (u.m * u.K))
    sim.wall.thickness = 3.25 * u.mm
    sim.wall.set_boundary("fixed", fixed_temperature=291.4 * u.K)   # Table 4: wall ends at T_a

    layer = sim.add_layer("Activated Carbon", composition={"Activated Carbon": 1.0}, momentum="ergun",
                          length=0.286 * u.m, cells=60, bed_void_fraction=0.456, particle_diameter=3 * u.mm)
    # The paper's dispersion: Eq. 4, D_L = 0.7 D_m + v_F r_p (D_m = 0.605e-5 m2/s, Table 1), times the
    # factor 4 fitted to the outlet composition (Sect. 3.2); a constant per flow rate.
    layer.set_dispersion("constant", dispersion_coefficient=paper_dispersion(cfg["feed_flow"]) * u.m**2 / u.s)
    layer.set_heat_dispersion("yagi_kunii")                    # axial conductivity from Yagi & Kunii
    layer.set_internal_heat_transfer("constant", internal_htc=cfg["internal_htc"])
    # Wall to ambient by natural convection, Holman correlation Nu = c Re^m Pr^n.
    layer.set_external_heat_transfer("holman", holman_c=0.027, holman_m=0.8, holman_n=1 / 3,
                                     external_velocity=0.5 * u.m / u.s)

    # ---- initial state: the bed full of helium at 1 bar and 291.4 K ----------------------
    sim.initial_conditions.pressure = 1 * u.bar
    sim.initial_conditions.temperature = 291.4 * u.K
    sim.initial_conditions.wall_temperature = 291.4 * u.K
    sim.initial_conditions.velocity = 0 * u.m / u.s
    sim.initial_conditions.y["Helium"] = 1.0
    sim.initial_conditions.y["CarbonDioxide"] = 0.0

    # ---- one feed step: pure CO2 in at the bottom, the top open at 1 bar -------------------
    sim.cycle.number_of_cycles = 1
    sim.cycle.step_time_mode = "time"
    feed = sim.add_step("Feed", time=3000 * u.s)
    feed.bottom.set_streams("inlet", [
        sk.Stream("CO2 Inlet", flow=sk.volumetric_flow(cfg["feed_flow"]), y={"Helium": 0.0, "CarbonDioxide": 1.0},
                  temperature=293 * u.K),
    ])
    feed.top.set_pressure("outlet", 1 * u.bar, y={"CarbonDioxide": 0.0, "Helium": 1.0}, temperature=293 * u.K)

    # ---- solver (SUNDIALS IDA), cyclic steady state, energy balance ----------------------
    sim.settings.init_step = 1e-10 * u.s
    sim.settings.min_step = 1e-15 * u.s
    sim.settings.max_step = 0 * u.s                              # no limit on the integrator step
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
    sim.energy_balance.isothermal_temperature = 298.15 * u.K
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
