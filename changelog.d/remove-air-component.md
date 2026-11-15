### Removed — Air
- Air is no longer a component. It was a lumped pseudo-component left over from CoolProp: one isotherm and one set of kinetics for N₂ and O₂ together, no element formula and no CAS number, and a heat-capacity polynomial mixed from the N₂, O₂ and Ar entries of TM-4513. A gas separation needs the real species of air, and an inert carrier is better modelled as Nitrogen. No example deck, adsorbent library entry or cycle template used it.
- With it go its rows in the generated component table, the ideal-gas heat-capacity, viscosity and conductivity tables, and its CoolProp reference values in the tests.

### Changed — decks with Air load with Nitrogen
- A deck, adsorbent or cycle file that has Air loads with Nitrogen in its place: the fluid's component list, every component map (mole fractions, isotherms, kinetics, heats of adsorption, heat capacities, KPI and reaction entries), the component name fields and the parameter paths (`adsorbent[1].isotherm[Air].b` becomes `adsorbent[1].isotherm[Nitrogen].b`). Stream and endpoint names are not changed. The physics follow the Nitrogen properties (molar mass, transport, heat capacity) from then on, so results move slightly.
- If the deck already has Nitrogen, Air's mole fractions are added to Nitrogen's and Air's other entries are dropped, so Nitrogen keeps its own isotherm and kinetics.
- New: `Components::renamedComponents` (`fluid.h`) and `MigrateRenamedComponents` (`adsim/objects/component_migration.h`), which every load function runs on the JSON text before it reads it. New tests: `ComponentMigration.*`.
