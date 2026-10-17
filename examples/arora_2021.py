"""Arora & Hasan (2021): medical oxygen from air on LiLSX, 4-step PSA and PVSA.

Builds both cycles from scratch in code: the air components, the three
adsorbents of the paper's Table 2 (LiX, LiLSX, 5A; the column is packed with
LiLSX, the other two are there to swap in), the column, the four steps and the
KPI definition. Every value is set explicitly, with its unit.

Each cycle runs 50 cycles and the O2 purity of the collected product and the
net O2 recovery are compared with the paper's optimised model results
(Table 4). The purge and pressurisation inlets are KPI recycle endpoints: the
O2 they return is subtracted from the product for the recovery (net, as in
the paper), while the purity is that of the collected product.

    python examples/arora_2021.py
    python examples/arora_2021.py --case LiLSX_PVSA

``build("LiLSX_PSA")`` returns the configured simulation without running it
(python/tests/test_examples_build.py checks it against examples/decks/).

Reference: A. Arora, M.M.F. Hasan, Flexible oxygen concentrators for medical
applications, Scientific Reports 11 (2021) 14317. doi:10.1038/s41598-021-93796-3
"""

from __future__ import annotations

import sys

import _common as c
import skarstrom as sk
from skarstrom import ureg as u  # pint units: 3.2 * u.cm, 40.98 * u.kJ / u.mol

# Table 4 optimal operating points (core/test/papers/arora_hasan_2021/). The
# purge flow is the paper's purge velocity factor gamma = v_purge / v_feed
# turned into a molar flow at the purge pressure: n_purge = gamma n_feed P3/P1.
CYCLES = {
    #               product    high      depress. / purge   step times [s]               air feed     purge
    #               O2 spec    pressure  pressure           prod, depr, purge, press     [mol/s]      factor
    "LiLSX_PSA":  dict(o2_spec=0.85, p_high=3.5e5, p_low=1.0e5, times=(3.0, 1.0, 1.4, 1.0), feed=0.23, gamma=1.97,
                       title="Arora 2021 - LiLSX PSA", cycle="PSA"),
    "LiLSX_PVSA": dict(o2_spec=0.9, p_high=3.5e5, p_low=0.5e5, times=(3.2, 1.0, 1.3, 1.0), feed=0.25, gamma=2.32,
                       title="Arora 2021 - LiLSX PVSA", cycle="PVSA"),
}

AIR = {"Oxygen": 0.21, "Nitrogen": 0.79}   # feed and initial gas, mole fractions


def add_adsorbents(sim: sk.Simulation) -> None:
    """LiX, LiLSX and 5A (Table 2): dual-site Langmuir in partial pressure, constant LDF coefficients."""
    adsorbents = {
        #  name             particle    heat         particle  isotherm sites (m [mol/kg], b0 [1/Pa], q [J/mol])            LDF k [1/s]     heat of adsorption
        #                   density     capacity     voidage   O2: (m1, q1, m2, q2), N2: (m1, q1, m2, q2)                   O2, N2          O2, N2 [J/mol]
        "LiX Zeolite":   (1200 * u.kg / u.m**3, 1172 * u.J / (u.kg * u.K), 0.33, (1.5, 21774.6, 0.6075, 24624.1), (1.5, 29195.0, 1.5, 26260.5),
                          (104.07, 42.84), (13221.44, 23430.4)),
        "LiLSX Zeolite": (1035 * u.kg / u.m**3, 1172 * u.J / (u.kg * u.K), 0.33, (1.5, 19606.9, 1.5, 22906.2), (1.5, 27923.9, 0.6613, 27923.9),
                          (104.07, 42.84), (14661.44, 25451.65)),
        "Zeolite 5A":    (1160 * u.kg / u.m**3, 1338.9 * u.J / (u.kg * u.K), 0.65, (1.5, 22312.1, 1.5, 22312.1), (1.5, 26140.8, 1.5, 23284.3),
                          (12.07, 4.02), (13221.44, 22886.48)),
    }
    for name, (density, cp, ep, o2, n2, k, hoa) in adsorbents.items():
        ads = sim.add_adsorbent(name, solid_density=density, heat_capacity=cp,
                                thermal_conductivity=0.15 * u.W / (u.m * u.K), particle_void_fraction=ep)

        def sites(m1, q1, m2, q2):
            return {"m1": m1 * u.mol / u.kg, "b1": 1e-10 / u.Pa, "q1": q1 * u.J / u.mol,
                    "m2": m2 * u.mol / u.kg, "b2": 1e-10 / u.Pa, "q2": q2 * u.J / u.mol}

        ads.set_isotherm("dual_site_langmuir_p", {"Oxygen": sites(*o2), "Nitrogen": sites(*n2)},
                         components=["Oxygen", "Nitrogen"])
        ads.set_rate_model("Oxygen", "ldf", "constant", k=k[0] / u.s)
        ads.set_rate_model("Nitrogen", "ldf", "constant", k=k[1] / u.s)
        ads.heat_of_adsorption["Oxygen"] = hoa[0] * u.J / u.mol
        ads.heat_of_adsorption["Nitrogen"] = hoa[1] * u.J / u.mol


def build(case: str) -> sk.Simulation:
    """The simulation of one cycle ("LiLSX_PSA" or "LiLSX_PVSA"), ready to run."""
    cfg = CYCLES[case]
    sim = sk.Simulation.new(
        title=cfg["title"],
        description=f"4-step Skarstrom {cfg['cycle']}, air separation to O2, LiLSX adsorbent (Table 4)",
        author="Arora & Hasan (2021)",
    )

    # ---- fluid: air as O2/N2 -------------------------------------------------------------
    sim.add_component("Oxygen")
    sim.add_component("Nitrogen")
    sim.fluid.adsorbed_heat_capacity_mode = "constant"
    sim.fluid.heat_of_adsorption_convention = "isosteric"

    # ---- adsorbents ------------------------------------------------------------------------
    add_adsorbents(sim)

    # ---- column: 0.127 m x 10 cm packed with LiLSX ---------------------------------------
    sim.bed.set_shape("cylinder", diameter=0.1 * u.m)
    sim.bed.ambient_temperature = 298 * u.K
    sim.bed.set_thermal_equilibrium("lte")                     # gas and solid at one temperature
    sim.bed.set_particle_heat_transfer("constant", particle_htc=100 * u.W / (u.m**2 * u.K))
    for void in (sim.bed.bottom_void, sim.bed.top_void):       # no dead volume at the column ends
        void.enabled = False
        void.volume = 1e-4 * u.m**3
        void.resistance_coefficient = 0 * u.Pa * u.s / u.m
    sim.wall.set_material("ss316l", density=7800 * u.kg / u.m**3, heat_capacity=537 * u.J / (u.kg * u.K),
                          thermal_conductivity=16.2 * u.W / (u.m * u.K))
    sim.wall.thickness = 17.5 * u.mm
    sim.wall.set_boundary("zero_flux")

    layer = sim.add_layer("LiLSX Zeolite", composition={"LiLSX Zeolite": 1.0}, momentum="kozeny_carman",
                          length=0.127 * u.m, cells=30, bed_void_fraction=0.35, particle_diameter=35 * u.mm)
    layer.set_dispersion("constant", dispersion_coefficient=0 * u.m**2 / u.s)   # no axial dispersion
    layer.set_heat_dispersion("constant", heat_dispersion_coefficient=0.29 * u.W / (u.m * u.K))
    layer.set_internal_heat_transfer("constant", internal_htc=70 * u.W / (u.m**2 * u.K))
    layer.set_external_heat_transfer("constant", external_htc=70 * u.W / (u.m**2 * u.K))

    # ---- initial state: the bed full of air at the high pressure (paper assumption 7) -----
    sim.initial_conditions.pressure = cfg['p_high'] * u.Pa
    sim.initial_conditions.temperature = 298 * u.K
    sim.initial_conditions.wall_temperature = 298 * u.K
    sim.initial_conditions.velocity = 0 * u.m / u.s
    sim.initial_conditions.y["Nitrogen"] = AIR["Nitrogen"]
    sim.initial_conditions.y["Oxygen"] = AIR["Oxygen"]

    # ---- cycle: production, reverse depressurisation, reverse purge, pressurisation -----
    sim.cycle.number_of_cycles = 50
    sim.cycle.step_time_mode = "time"
    add_cycle(sim, cfg)

    # ---- solver (SUNDIALS IDA), cyclic steady state, energy balance ----------------------
    sim.settings.init_step = 1e-10 * u.s
    sim.settings.min_step = 1e-15 * u.s
    sim.settings.max_step = 0.01 * u.s
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
    sim.css.method = "none"                                    # run all 50 cycles
    sim.css.min_cycles = 5
    sim.css.molar_absolute_tolerance = 1e-4
    sim.css.molar_relative_tolerance = 1e-4
    sim.css.mass_balance_tolerance = 0.005
    sim.energy_balance.energy_balance_mode = "non_isothermal"
    sim.energy_balance.isothermal_use_initial_temperature = True
    sim.energy_balance.isothermal_temperature = 298.15 * u.K
    return sim


def add_cycle(sim: sk.Simulation, cfg: dict) -> None:
    """The four steps of the Skarstrom cycle; the purge and pressurisation gas is O2 product."""
    p_high, p_low = cfg["p_high"], cfg["p_low"]
    t_prod, t_depr, t_purge, t_press = cfg["times"]
    product = {"Oxygen": cfg["o2_spec"], "Nitrogen": 1.0 - cfg["o2_spec"]}   # recycled O2 at the product spec
    n_purge = cfg["gamma"] * cfg["feed"] * (p_low / p_high)                  # [mol/s]

    # 1. Production: air in at the bottom, O2 product out at the top at the high pressure.
    production = sim.add_step("Production", time=t_prod * u.s)
    production.bottom.set_streams("inlet", [
        sk.Stream("", flow=sk.molar_flow(cfg['feed'] * u.mol / u.s), y=AIR, temperature=298 * u.K),
    ])
    production.top.set_pressure("outlet", p_high * u.Pa, y=product, temperature=293 * u.K)

    # 2. Reverse depressurisation: top closed, the feed end vents to the low pressure.
    depressurization = sim.add_step("Depressurization", time=t_depr * u.s)
    depressurization.bottom.set_pressure("outlet", sk.linear(p_low * u.Pa), y=product, temperature=293 * u.K)
    depressurization.top.set_closed()

    # 3. Reverse purge: O2 product in at the top, waste out of the bottom.
    purge = sim.add_step("Purge", time=t_purge * u.s)
    purge.bottom.set_pressure("outlet", p_low * u.Pa, y=AIR, temperature=293 * u.K)
    purge.top.set_streams("inlet", [
        sk.Stream("", flow=sk.molar_flow(n_purge * u.mol / u.s), y=product, temperature=298 * u.K),
    ])

    # 4. Pressurisation with O2 product from the top, bottom closed.
    pressurization = sim.add_step("Pressurization", time=t_press * u.s)
    pressurization.bottom.set_closed()
    pressurization.top.set_pressure("inlet", sk.linear(p_high * u.Pa), y=product, temperature=298 * u.K)

    # O2 KPIs: fed with the air, collected at the production outlet; the O2
    # recycled to purge and pressurisation is subtracted for the net recovery.
    sim.kpi.track("Oxygen", feed=[production.bottom[0]], product=[production.top],
                  recycle=[purge.top[0], pressurization.top])


CASES = {name: name for name in CYCLES}                      # case -> reference_values.csv case
DECKS = {name: f"Arora2021_{name}" for name in CYCLES}       # the shipped deck of each case
QUANTITIES = [("O2 purity", "Oxygen.purity", "O2_purity"), ("O2 recovery", "Oxygen.recovery", "O2_recovery")]


def main(argv=None) -> int:
    return c.cyclic_family(argv, __doc__, "arora_hasan_2021", CASES, QUANTITIES, "Oxygen", "arora_2021.png",
                           build=build)


if __name__ == "__main__":
    sys.exit(main())
