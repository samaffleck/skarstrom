"""Illustrative humid DAC TVSA: a cycle on APDES-NFC that crosses water saturation.

Not a reproduction of a paper: a five-step temperature-vacuum swing cycle for
direct air capture on the library adsorbent "APDES-NFC" (amine-functionalised
nanofibrillated cellulose; GAB water isotherm, humid Toth CO2 isotherm), fed
with humid air at 20 C and 49 % RH:

1. Adsorption (4 h): humid air in at the bottom, CO2- and water-lean air out of the top.
2. Evacuation: the top closed, the bottom pulled down to 0.2 bar.
3. Heating: a 95 C medium in a jacket and an internal exchanger heats the bed at
   0.2 bar; water and CO2 leave the bottom (the product).
4. Cooling: both ends closed, 20 C water in the same jacket and exchanger.
5. Pressurisation: humid air back to 1 bar through the bottom.

The heating step leaves water-rich gas in the bed. In the cooling step the jacket cools
the steel wall faster than the packing gives up its heat, and the water partial pressure
of the gas passes the saturation pressure at the wall temperature: the run-health check
SUPERSATURATION warns (the cold-wall measure, up to 2.6 times saturation). The gas itself
comes close but stays just below saturation (0.96 to 1.00 in the first cycles, highest
when the feed compresses the remaining vapour in the pressurisation step), because the
sorbent takes the water back up. The engine does not model condensation, so the warning
says that the run has left the model's range; the deck exists to show that warning on a
realistic cycle.

The feed, the column diameter and the bed-to-wall and wall-to-ambient
coefficients are those of the Illustrative_Humid_DAC_Breakthrough deck. The step
times, the vacuum level, the medium temperatures and the exchanger are round
numbers chosen for the example, not values from a paper. The heating step takes
the CO2 isotherm above its 296-343 K fit range.

    python examples/illustrative_humid_dac_tvsa.py

``build()`` returns the configured simulation without running it
(python/tests/test_examples_build.py checks it against examples/decks/).
"""

from __future__ import annotations

import sys

import _common as c
import skarstrom as sk
from skarstrom import ureg as u

T_AMBIENT = 293.15 * u.K                                   # 20 C
FEED = {"Nitrogen": 1 - 400e-6 - 0.0115, "CarbonDioxide": 400e-6, "Water": 0.0115}   # 49 % RH at 20 C, 1 bar
P_VACUUM = 0.2 * u.bar


def build(case: str = "TVSA") -> sk.Simulation:
    """The illustrative humid DAC TVSA cycle, ready to run (there is one case, "TVSA")."""
    sim = sk.Simulation.new(
        title="Illustrative humid DAC TVSA - APDES-NFC, crossing water saturation",
        description="ILLUSTRATIVE, not a validated reproduction. Five-step temperature-vacuum swing cycle for direct "
                    "air capture on the library adsorbent APDES-NFC (GAB water, humid Toth CO2): humid air "
                    "(400 ppm CO2, 49 % RH at 20 C) adsorption, evacuation to 0.2 bar, heating with a 95 C medium "
                    "at 0.2 bar (water and CO2 product at the feed end), closed cooling with 20 C water, feed "
                    "pressurisation. Step times, vacuum and medium are round numbers for the example. On cooling, "
                    "the water-rich gas left by the heating step passes its saturation pressure at the cold wall: "
                    "condensation is not modelled, and the SUPERSATURATION run-health warning says so.",
        author="Skarstrom examples",
    )
    for component in ("Nitrogen", "CarbonDioxide", "Water"):
        sim.add_component(component)

    # ---- adsorbent: the library entry, unchanged (N2 is inert on it) ----------------------
    ads_id = sim.add_library_adsorbent("apdes_nfc_stampi2020")
    ads = sim.adsorbents.by_id(ads_id)
    ads.set_rate_model("Nitrogen", "ldf", "constant", k=0.1 / u.s)

    # ---- column: the humid breakthrough deck's, with a shorter bed -------------------------
    sim.bed.set_shape("cylinder", diameter=8 * u.cm)
    sim.bed.ambient_temperature = T_AMBIENT
    sim.wall.set_material("ss316l")
    sim.wall.thickness = 2 * u.mm
    layer = sim.add_layer("APDES-NFC", composition={ads.name: 1.0}, momentum="ergun", length=0.1 * u.m, cells=20,
                          bed_void_fraction=0.4, particle_diameter=7.5 * u.mm)
    layer.set_internal_heat_transfer("constant", internal_htc=3 * u.W / (u.m**2 * u.K))
    layer.set_external_heat_transfer("constant", external_htc=26 * u.W / (u.m**2 * u.K))
    sim.energy_balance.energy_balance_mode = "non_isothermal"

    sim.initial_conditions.pressure = 1 * u.bar
    sim.initial_conditions.temperature = T_AMBIENT
    sim.initial_conditions.wall_temperature = T_AMBIENT
    sim.initial_conditions.velocity = 0 * u.m / u.s
    sim.initial_conditions.y["Nitrogen"] = 1.0
    sim.initial_conditions.y["CarbonDioxide"] = 0.0
    sim.initial_conditions.y["Water"] = 0.0

    # ---- cycle -------------------------------------------------------------------------------
    sim.cycle.number_of_cycles = 60

    adsorption = sim.add_step("Adsorption", time=4 * u.h)
    adsorption.bottom.set_streams("inlet", [
        sk.Stream("Humid air", flow=sk.velocity(0.01 * u.m / u.s), y=FEED, temperature=T_AMBIENT)])
    adsorption.top.set_pressure("outlet", 1 * u.bar)

    evacuation = sim.add_step("Evacuation", time=2 * u.min)
    evacuation.bottom.set_pressure("outlet", sk.linear(P_VACUUM))
    evacuation.top.set_closed()

    heating = sim.add_step("Heating", time=30 * u.min)
    heating.bottom.set_pressure("outlet", P_VACUUM)
    heating.top.set_closed()

    cooling = sim.add_step("Cooling", time=20 * u.min)
    cooling.bottom.set_closed()
    cooling.top.set_closed()

    # Heating and cooling media in a jacket on the wall and an internal exchanger in the packing
    # (50 m2 of tube per m3 of bed, 40 W/(m2 K) on the bed side), as in the Illustrative_TSA_13X deck.
    for step, medium_temperature in ((heating, u.Quantity(95, "degC")), (cooling, u.Quantity(20, "degC"))):
        step.medium.enabled = True
        step.medium.temperature = medium_temperature
        step.medium.heat_transfer_coefficient = 500 * u.W / (u.m**2 * u.K)
        step.medium.exchanger_area_per_volume = 50 * u.m**2 / u.m**3
        step.medium.exchanger_heat_transfer_coefficient = 40 * u.W / (u.m**2 * u.K)

    pressurisation = sim.add_step("Pressurisation", time=2 * u.min)
    pressurisation.bottom.set_pressure("inlet", sk.linear(1 * u.bar), y=FEED, temperature=T_AMBIENT)
    pressurisation.top.set_closed()

    sim.kpi.track("CarbonDioxide", feed=[adsorption.bottom, pressurisation.bottom],
                  product=[evacuation.bottom, heating.bottom])

    sim.settings.max_cycles_to_keep = 1
    sim.css.method = "mass_balance_component"
    sim.css.min_cycles = 5
    sim.css.molar_absolute_tolerance = 1e-4
    sim.css.molar_relative_tolerance = 1e-4
    sim.css.mass_balance_tolerance = 0.005
    return sim


DECKS = {"TVSA": "Illustrative_Humid_DAC_TVSA"}       # the shipped deck


def main(argv=None) -> int:
    args = c.parser(__doc__.splitlines()[0], ["TVSA"]).parse_args(argv)
    sim = build()
    if args.cycles is not None:
        sim.cycle.number_of_cycles = int(args.cycles)
    result = c.run(sim, "Illustrative humid DAC TVSA", args.quiet)
    print(f"  {c.css_status(sim)}")
    c.print_kpis(sim)
    for w in result.warnings:
        if w.code == "SUPERSATURATION":
            print(f"  {w.code}: {w.message}")
    cycle = sim.run_health().cycles[-1]
    print(f"  last cycle: max relative saturation {cycle.supersaturation_gas:.2f} in the gas, "
          f"{cycle.supersaturation_wall:.2f} at the wall")
    return 0


if __name__ == "__main__":
    sys.exit(main())
