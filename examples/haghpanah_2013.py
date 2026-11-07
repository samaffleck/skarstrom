"""Haghpanah et al. (2013): four-step CO2/N2 PVSA on zeolite 13X.

Builds each case from scratch in code: the fluid, the adsorbent and its
isotherm and rate model, the column and its heat transfer, the four steps with
their pressure profiles, pumps and compressors, the KPI definition and the
solver settings. Every value is set explicitly, with its unit.

The cycle is pressurisation with feed gas, high-pressure adsorption, forward
blowdown and reverse evacuation (the CO2 product), with exponential pressure
profiles of time constant 2 s. Cases:

    Sec61          the reference design of Sect. 6.1 (Table 4)
    PL002 ... PL010  the minimum-energy designs of Table 6 at evacuation pressures of
                   0.02-0.1 bar (also verified by Glover et al. 2026, Table 4)
    Fig5Benchmark  the CSS benchmark design of Glover et al. (2026), Fig. 5, on 10 cells

Each runs to cyclic steady state; the CO2 purity and recovery of the
evacuation product are compared with the published values.

    python examples/haghpanah_2013.py                  # Sec61 and the five Table 6 designs
    python examples/haghpanah_2013.py --case PL002
    python examples/haghpanah_2013.py --case Sec61 --cycles 20

``build("Sec61")`` returns the configured simulation without running it
(python/tests/test_examples_build.py checks it against examples/decks/).

References: R. Haghpanah, A. Majumder, R. Nilam, A. Rajendran, S. Farooq,
I. A. Karimi, M. Amanullah, Ind. Eng. Chem. Res. 52 (2013) 4249-4265,
doi:10.1021/ie302658y (isotherm: Krishnamurthy, Haghpanah, Rajendran, Farooq,
APCChE 2012, as cited there as ref. 38); T. Glover, T. Papathanasiou, R. Pini,
arXiv:2606.21523 (2026).
"""

from __future__ import annotations

import sys

import _common as c
import skarstrom as sk
from skarstrom import ureg as u  # pint units

# Operating points: step times, the pressure at the end of each pressure step,
# and the interstitial feed velocity. Table 6 does not list the pressurisation
# time: 20 s, as fixed in both papers' optimisations.
DESIGNS = {
    "Sec61": dict(t_press=15, t_ads=15, t_bd=30, t_evac=40, p_h=1.0, p_i=0.2, p_l=0.1, v_feed=1.0, cells=30),
    "PL002": dict(t_press=20, t_ads=94.89, t_bd=122.84, t_evac=189.46, p_h=1.00, p_i=0.18, p_l=0.02, v_feed=0.58, cells=30),
    "PL003": dict(t_press=20, t_ads=69.36, t_bd=30.95, t_evac=124.27, p_h=1.76, p_i=0.22, p_l=0.03, v_feed=0.41, cells=30),
    "PL005": dict(t_press=20, t_ads=91.63, t_bd=161.32, t_evac=195.24, p_h=3.54, p_i=0.43, p_l=0.05, v_feed=0.19, cells=30),
    "PL007": dict(t_press=20, t_ads=70.93, t_bd=54.13, t_evac=88.45, p_h=5.73, p_i=0.53, p_l=0.07, v_feed=0.15, cells=30),
    "PL010": dict(t_press=20, t_ads=97.36, t_bd=94.52, t_evac=128.45, p_h=9.04, p_i=1.07, p_l=0.10, v_feed=0.11, cells=30),
    "Fig5Benchmark": dict(t_press=20, t_ads=50, t_bd=100, t_evac=100, p_h=8.0, p_i=1.5, p_l=0.05, v_feed=0.8, cells=10),
}

FEED = {"CarbonDioxide": 0.15, "Nitrogen": 0.85}       # flue gas, mole fractions
BACKFLOW = {"CarbonDioxide": 0.0, "Nitrogen": 1.0}     # gas entering an outlet if the flow reverses
LAMBDA = 0.5 / u.s                                     # pressure-profile time constant (Table 4)
BED_VOID = 0.37
PARTICLE_RADIUS = 1.0 * u.mm
D_M = 1.60e-5 * u.m**2 / u.s                           # molecular diffusivity (Table 3)


def ramp(end):
    """exp(-lambda t) approach from the bed pressure at the start of the step to ``end`` (Eqs. 57-59)."""
    return sk.double_exponential(end, rate_1=LAMBDA, rate_2=LAMBDA, weight=1.0)


def pump(kind: str) -> sk.Pump:
    """Isentropic compression (gamma 1.4) at 72 % efficiency from or to 1 bar (Eqs. 66-69).
    The gas enters the bed at the feed temperature (isothermal outlet model)."""
    side = {"compressor": dict(suction_pressure=1 * u.bar), "vacuum_pump": dict(discharge_pressure=1 * u.bar)}[kind]
    return sk.Pump(kind, energy="isothermal", efficiency_model="constant", efficiency=0.72,
                   polytropic_index=1.4, **side)


def build(case: str) -> sk.Simulation:
    """The simulation of one case (see DESIGNS), ready to run."""
    d = DESIGNS[case]
    p_h, p_i, p_l = d["p_h"] * u.bar, d["p_i"] * u.bar, d["p_l"] * u.bar
    v_feed = d["v_feed"] * u.m / u.s
    sim = sk.Simulation.new(title=f"Haghpanah 2013 - {case}", author="Skarstrom examples",
                            description=f"Case {case} of the four-step CO2/N2 PVSA on 13X of Haghpanah et al. (2013).")

    # ---- fluid: constant properties (Table 3) ---------------------------------------
    sim.add_component("CarbonDioxide")
    sim.add_component("Nitrogen")
    sim.fluid.gas_heat_capacity_mode = "constant"
    sim.fluid.adsorbed_heat_capacity_mode = "constant"
    for component in ("CarbonDioxide", "Nitrogen"):
        sim.fluid.gas_heat_capacity[component] = 30.7 * u.J / (u.mol * u.K)
        sim.fluid.adsorbed_heat_capacity[component] = 30.7 * u.J / (u.mol * u.K)
    sim.fluid.viscosity_mode = "constant"
    sim.fluid.viscosity = 1.72e-5 * u.Pa * u.s
    # The energy balance (Eq. A.8) carries -C_pa T dq/dt next to the heat of adsorption.
    sim.fluid.heat_of_adsorption_convention = "internal_energy"

    # ---- adsorbent: zeolite 13X ----------------------------------------------------------
    # The paper's balances hold gas in the bed voidage only (Eqs. A.2-A.3): no macropore
    # holdup, so particle_void_fraction is 0; the particle voidage 0.35 enters the rate (below).
    ads = sim.add_adsorbent("13X", solid_density=1130 * u.kg / u.m**3, heat_capacity=1070 * u.J / (u.kg * u.K),
                            thermal_conductivity=0.63 * u.W / (u.m * u.K), particle_void_fraction=0.0)
    # Table 2 (Krishnamurthy et al. 2012 fit): competitive dual-site Langmuir in concentration
    # form, b(T) = b0 exp(-dU / RT); CO2 on two sites, N2 on one.
    ads.set_isotherm("extended_dual_site_langmuir_c", {
        "CarbonDioxide": {"m1": 3.09 * u.mol / u.kg, "b1": 8.65e-7 * u.m**3 / u.mol, "q1": 36641.21 * u.J / u.mol,
                          "m2": 2.54 * u.mol / u.kg, "b2": 2.63e-8 * u.m**3 / u.mol, "q2": 35690.66 * u.J / u.mol},
        "Nitrogen":      {"m1": 5.84 * u.mol / u.kg, "b1": 2.50e-6 * u.m**3 / u.mol, "q1": 15800 * u.J / u.mol,
                          "m2": 0 * u.mol / u.kg, "b2": 0 * u.m**3 / u.mol, "q2": 0 * u.J / u.mol},
    }, components=["CarbonDioxide", "Nitrogen"])
    # Eq. 40: LDF coefficient (c / q*) 15 eps_p D_m / (tau r_p^2), with the constant D_m of Table 3.
    for component in ("CarbonDioxide", "Nitrogen"):
        ads.set_rate_model(component, "ldf", "macropore_capacity", tortuosity=3.0, macropore_void_fraction=0.35,
                           molecular_diffusion_coefficient=D_M)
    # Constant heats (the paper prints none): CO2 the isosteric weighting of its two site
    # energies at the feed partial pressure, N2 its single site energy.
    ads.heat_of_adsorption["CarbonDioxide"] = 35983.8 * u.J / u.mol
    ads.heat_of_adsorption["Nitrogen"] = 15800 * u.J / u.mol

    # ---- column: 1 m, inner radius 0.1445 m, outer 0.162 m (Table 3) ------------------
    sim.bed.set_shape("cylinder", diameter=0.289 * u.m)
    sim.bed.ambient_temperature = 298.15 * u.K
    sim.wall.set_material("ss316l", density=7800 * u.kg / u.m**3, heat_capacity=502 * u.J / (u.kg * u.K),
                          thermal_conductivity=16 * u.W / (u.m * u.K))
    sim.wall.thickness = 0.0175 * u.m
    sim.wall.set_boundary("fixed", fixed_temperature=298.15 * u.K)    # T_w = T_a at both ends (Eq. 12)
    layer = sim.add_layer("13X", composition={"13X": 1.0}, momentum="blake_kozeny", length=1 * u.m,
                          cells=d["cells"], bed_void_fraction=BED_VOID, particle_diameter=2 * PARTICLE_RADIUS)
    # Eq. 42: D_L = 0.7 D_m + v_feed r_p, constant for the design.
    layer.set_dispersion("constant", dispersion_coefficient=0.7 * D_M + v_feed * PARTICLE_RADIUS)
    layer.set_heat_dispersion("constant", heat_dispersion_coefficient=0.09 * u.W / (u.m * u.K))
    layer.set_internal_heat_transfer("constant", internal_htc=8.6 * u.W / (u.m**2 * u.K))
    layer.set_external_heat_transfer("constant", external_htc=2.5 * u.W / (u.m**2 * u.K))

    # ---- initial state: N2 at P_L (Sect. 6.1) ------------------------------------------
    sim.initial_conditions.pressure = p_l
    sim.initial_conditions.temperature = 298.15 * u.K
    sim.initial_conditions.wall_temperature = 298.15 * u.K
    sim.initial_conditions.y["CarbonDioxide"] = 0.0
    sim.initial_conditions.y["Nitrogen"] = 1.0

    # ---- cycle -------------------------------------------------------------------------
    sim.cycle.number_of_cycles = 500
    press = sim.add_step("Pressurisation", time=d["t_press"] * u.s)
    press.bottom.set_pressure("inlet", ramp(p_h), y=FEED, temperature=298.15 * u.K, pump=pump("compressor"))
    press.top.set_closed()

    ads_step = sim.add_step("Adsorption", time=d["t_ads"] * u.s)
    ads_step.bottom.set_streams("inlet", [sk.Stream("Feed inlet", flow=sk.velocity(BED_VOID * v_feed), y=FEED,
                                                    temperature=298.15 * u.K, pump=pump("compressor"))])
    ads_step.top.set_pressure("outlet", p_h, y=BACKFLOW)

    blowdown = sim.add_step("Blowdown", time=d["t_bd"] * u.s)
    blowdown.bottom.set_closed()
    blowdown.top.set_pressure("outlet", ramp(p_i), y=BACKFLOW, pump=pump("vacuum_pump"))

    evacuation = sim.add_step("Evacuation", time=d["t_evac"] * u.s)
    evacuation.top.set_closed()
    evacuation.bottom.set_pressure("outlet", ramp(p_l), y=BACKFLOW, pump=pump("vacuum_pump"))

    # CO2 is fed in pressurisation and adsorption and recovered by evacuation (Eqs. 60-61).
    sim.kpi.track("CarbonDioxide", feed=[press.bottom, ads_step.bottom["Feed inlet"]], product=[evacuation.bottom])

    # ---- solver: run to cyclic steady state ------------------------------------------------
    sim.settings.max_step = 2 * u.s
    sim.css.method = "mass_balance_component"
    sim.css.min_cycles = 5
    return sim


CASES = {name: name for name in DESIGNS}   # case -> reference_values.csv case
QUANTITIES = [
    ("CO2 purity vs Haghpanah (Sect. 6.1)", "CarbonDioxide.purity", "CO2_purity"),
    ("CO2 recovery vs Haghpanah (Sect. 6.1)", "CarbonDioxide.recovery", "CO2_recovery"),
    ("CO2 purity vs Glover (Table 4)", "CarbonDioxide.purity", "purity_G"),
    ("CO2 recovery vs Glover (Table 4)", "CarbonDioxide.recovery", "recovery_G"),
]
DECKS = {name: f"Haghpanah2013_{name}" for name in DESIGNS}   # the shipped deck of each case


def main(argv=None) -> int:
    return c.cyclic_family(argv, __doc__, "haghpanah_2013", CASES, QUANTITIES, "CarbonDioxide",
                           "haghpanah_2013.png", default=[n for n in DESIGNS if n != "Fig5Benchmark"],
                           values=c.evacuation_product_kpis, build=build)


if __name__ == "__main__":
    sys.exit(main())
