"""Vysyaraju et al. (2024): helium recovery from He/N2 by 5-step PSA on binderless 5A.

Builds the experiment-1 cycle from scratch in code: the fluid, the zeolite 5A
adsorbent and its N2 isotherm, the column and its heat transfer, the five steps
(pressurisation, adsorption, pressure-equalisation depressurisation, blowdown,
pressure-equalisation repressurisation) with their boundary conditions, the
connection between the two equalisation steps and the solver settings. Every
value is set explicitly, with its unit, so this file is also a template for a
PSA cycle with pressure equalisation.

The cycle runs to cyclic steady state and the He purity and recovery are
compared with the measured values and the paper's model (Table 3).

    python examples/vysyaraju_2024.py

``build("PSA_Exp1")`` returns the configured simulation without running it
(python/tests/test_examples_build.py checks it against examples/decks/).

Reference: R. Vysyaraju, J.A. Sawada, L. Estupinan Perez, A. Rajendran,
Chemical Engineering Journal 499 (2024) 155929. doi:10.1016/j.cej.2024.155929
"""

from __future__ import annotations

import sys

import _common as c
import skarstrom as sk
from skarstrom import ureg as u  # pint units: 3.2 * u.cm, 40.98 * u.kJ / u.mol

FEED = {"Nitrogen": 0.9, "Helium": 0.1}                # He/N2 feed, mole fractions
NO_GAS = {"Nitrogen": 0.0, "Helium": 0.0}              # outlet stream: no backflow composition


def add_zeolite_5a(sim: sk.Simulation) -> sk.Adsorbent:
    """Binderless zeolite 5A: N2 adsorbs (dual-site Langmuir, Table 2), He is inert."""
    ads = sim.add_adsorbent(
        "Zeolite5A",
        solid_density=1183.6 * u.kg / u.m**3,                  # particle density
        heat_capacity=960 * u.J / (u.kg * u.K),
        thermal_conductivity=0.15 * u.W / (u.m * u.K),
        particle_void_fraction=0.33,
    )
    ads.set_isotherm("dual_site_langmuir_c", {
        "Nitrogen": {"m1": 3.072 * u.mol / u.kg, "b1": 1.576e-6 * u.m**3 / u.mol, "q1": 18.09 * u.kJ / u.mol,
                     "m2": 0.646 * u.mol / u.kg, "b2": 9.195e-6 * u.m**3 / u.mol, "q2": 19.36 * u.kJ / u.mol},
    }, components=["Nitrogen"])
    ads.set_isotherm("inert", components=["Helium"])
    # Macropore-diffusion-controlled constant LDF coefficient, k = 15 ep Dm / (tau rp^2)
    # (core/test/papers/vysyarajua_2024/README.md); He's coefficient is irrelevant (inert).
    ads.set_rate_model("Nitrogen", "ldf", "constant", k=104.98161937977187 / u.s)
    ads.set_rate_model("Helium", "ldf", "constant", k=0.5 / u.s)
    # The full isosteric heat is released (isosteric convention, see the README).
    ads.heat_of_adsorption["Nitrogen"] = 18.7 * u.kJ / u.mol
    ads.heat_of_adsorption["Helium"] = 0 * u.J / u.mol
    return ads


def build(case: str = "PSA_Exp1") -> sk.Simulation:
    """The five-step PSA cycle of experiment 1, ready to run."""
    if case != "PSA_Exp1":
        raise KeyError(f"no case {case!r}; the only case is 'PSA_Exp1'")
    sim = sk.Simulation.new(
        title="Vysyaraju 2024 - PSA exp1",
        description="Five-step He/N2 PSA with pressure equalisation on binderless zeolite 5A, "
                    "experiment 1 of Vysyaraju et al. (2024).",
        author="Skarstrom examples",
    )

    # ---- fluid ---------------------------------------------------------------------
    sim.add_component("Nitrogen")
    sim.add_component("Helium")
    sim.fluid.adsorbed_heat_capacity_mode = "constant"
    sim.fluid.heat_of_adsorption_convention = "isosteric"
    sim.fluid.adsorbed_heat_capacity["Nitrogen"] = 29.1 * u.J / (u.mol * u.K)
    sim.fluid.adsorbed_heat_capacity["Helium"] = 0 * u.J / (u.mol * u.K)

    # ---- adsorbent -------------------------------------------------------------------
    add_zeolite_5a(sim)

    # ---- column: 0.305 m x 3.4 cm stainless-steel tube --------------------------------------
    sim.bed.set_shape("cylinder", diameter=3.4 * u.cm)
    sim.bed.ambient_temperature = 298 * u.K
    sim.bed.set_thermal_equilibrium("lte")                     # gas and solid at one temperature
    sim.bed.set_particle_heat_transfer("constant", particle_htc=100 * u.W / (u.m**2 * u.K))
    for void in (sim.bed.bottom_void, sim.bed.top_void):       # no dead volume at the column ends
        void.enabled = False
        void.volume = 1e-4 * u.m**3
        void.resistance_coefficient = 0 * u.Pa * u.s / u.m
    sim.wall.set_material("ss316l", density=7800 * u.kg / u.m**3, heat_capacity=502 * u.J / (u.kg * u.K),
                          thermal_conductivity=16 * u.W / (u.m * u.K))
    sim.wall.thickness = 2 * u.mm
    sim.wall.set_boundary("zero_flux")

    layer = sim.add_layer("Zeolite5A", composition={"Zeolite5A": 1.0}, momentum="kozeny_carman",
                          length=0.305 * u.m, cells=30, bed_void_fraction=0.35, particle_diameter=2.06 * u.mm)
    layer.set_dispersion("ruthven", molecular_diffusion_coefficient=1.6e-5 * u.m**2 / u.s)
    layer.set_heat_dispersion("constant", heat_dispersion_coefficient=0 * u.W / (u.m * u.K))
    layer.set_internal_heat_transfer("constant", internal_htc=15 * u.W / (u.m**2 * u.K))
    layer.set_external_heat_transfer("constant", external_htc=15 * u.W / (u.m**2 * u.K))

    # ---- initial state: the bed full of He at 14.043 bar and 298 K --------------------------
    sim.initial_conditions.pressure = 14.043 * u.bar
    sim.initial_conditions.temperature = 298 * u.K
    sim.initial_conditions.wall_temperature = 298 * u.K
    sim.initial_conditions.velocity = 0 * u.m / u.s
    sim.initial_conditions.y["Nitrogen"] = 0.0
    sim.initial_conditions.y["Helium"] = 1.0

    # ---- cycle -------------------------------------------------------------------------
    sim.cycle.number_of_cycles = 60
    sim.cycle.step_time_mode = "time"

    # 1. PRES: feed gas through a valve from the 15.1 bar supply, top closed.
    pres = sim.add_step("PRES", time=30 * u.s)
    pres.bottom.set_streams("inlet", [
        sk.ValveStream("", pressure=15.1 * u.bar, cv=0.012, xt=0.53, y=FEED, temperature=298 * u.K),
    ])
    pres.top.set_closed()

    # 2. ADS: feed at a set velocity, He-rich product leaves the top at 15.1 bar.
    adsorption = sim.add_step("ADS", time=15 * u.s)
    adsorption.bottom.set_streams("inlet", [
        sk.Stream("", flow=sk.velocity(0.005 * u.m / u.s), y=FEED, temperature=298 * u.K),
    ])
    adsorption.top.set_pressure("outlet", 15.1 * u.bar, y=FEED, temperature=293 * u.K)

    # 3. PED: pressure-equalisation depressurisation; gas leaves the top through a valve.
    ped = sim.add_step("PED", time=30 * u.s)
    ped.bottom.set_closed()
    ped.top.set_streams("outlet", [
        sk.ValveStream("PE stream", pressure=7 * u.bar, cv=0.2, xt=0.53, y=NO_GAS, temperature=293 * u.K),
    ])

    # 4. BLO: blowdown of the N2-rich gas from the bottom down to 1 bar.
    blo = sim.add_step("BLO", time=30 * u.s)
    blo.bottom.set_pressure("outlet", sk.exponential(1 * u.bar), y=FEED, temperature=293 * u.K)
    blo.top.set_closed()

    # 5. PER: pressure-equalisation repressurisation; its inflow is the PED stream.
    per = sim.add_step("PER", time=30 * u.s)
    per.bottom.set_closed()
    per.top.set_streams("inlet", [
        sk.Stream("PE stream", flow=sk.velocity(0 * u.m / u.s), y=FEED, temperature=298 * u.K),
    ])
    sim.connect(ped.top["PE stream"], per.top["PE stream"])

    # ---- solver (SUNDIALS IDA), cyclic steady state, energy balance ----------------------
    sim.settings.init_step = 1e-10 * u.s
    sim.settings.min_step = 1e-15 * u.s
    sim.settings.max_step = 1 * u.s
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
    sim.css.method = "mass_balance_component"
    sim.css.min_cycles = 60
    sim.css.molar_absolute_tolerance = 1e-4
    sim.css.molar_relative_tolerance = 1e-4
    sim.css.mass_balance_tolerance = 0.005
    sim.energy_balance.energy_balance_mode = "non_isothermal"
    sim.energy_balance.isothermal_use_initial_temperature = True
    sim.energy_balance.isothermal_temperature = 298.15 * u.K
    return sim


CASES = {"PSA_Exp1": "PSA_Exp1"}                        # case -> reference_values.csv case
DECKS = {"PSA_Exp1": "Vysyaraju2024_PSA_exp1"}          # the shipped deck of each case
QUANTITIES = [
    ("He purity vs experiment", "Helium.purity", "He_purity_exp"),
    ("He recovery vs experiment", "Helium.recovery", "He_recovery_exp"),
    ("He purity vs paper model", "Helium.purity", "He_purity_model"),
    ("He recovery vs paper model", "Helium.recovery", "He_recovery_model"),
]


def kpis(sim):
    """The cycle has no KPI configuration; integrate the stream logs of the last cycle.

    He purity is the He fraction of the adsorption-step top product (step 1);
    He recovery is 1 - (He lost in the blowdown, step 3 bottom) / (He fed in
    steps 0 and 1, bottom), the conservation form the validation test uses.
    """
    m = c.endpoint_moles
    fed = m(sim, 0, "bottom", "Helium") + m(sim, 1, "bottom", "Helium")
    return {
        "Helium.purity": m(sim, 1, "top", "Helium") / m(sim, 1, "top"),
        "Helium.recovery": 1.0 - m(sim, 3, "bottom", "Helium") / fed,
    }


def main(argv=None) -> int:
    return c.cyclic_family(argv, __doc__, "vysyarajua_2024", CASES, QUANTITIES, "Helium", "vysyaraju_2024.png",
                           values=kpis, build=build)


if __name__ == "__main__":
    sys.exit(main())
