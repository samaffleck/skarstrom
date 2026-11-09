# Sources and Sinks

!!! note "New in 0.9"
    Named sources and sinks: `sim.add_source`, `sim.add_sink`, `endpoint.attach`, the `reservoir[ID].*` parameter paths and the `reservoirs` list of the deck view.

A cycle usually takes the same gas in more than one place. Two feed steps and a feed pressurisation step all draw the feed gas. A blowdown and a purge both vent to ambient. Without sources and sinks, you type that temperature, composition and pressure on every [boundary](boundary-conditions.md), and you change it in every place when the feed changes.

A **source** is a named gas state that inlets attach to. A **sink** is a named pressure that outlets attach to. Together the engine calls them **reservoirs**.

| Value | Path | A source | A sink |
|-------|------|----------|--------|
| Pressure [Pa] | `reservoir[ID].pressure` | Supply pressure | Back pressure |
| Temperature [K] | `reservoir[ID].temperature` | Temperature of the gas that enters | Not used |
| Mole fractions [-] | `reservoir[ID].y[Component]` | Composition of the gas that enters | Not used |

A sink has a pressure only. An outlet is zero-gradient in temperature and composition, also when the flow through it reverses (refer to [Backflow through an outlet](boundary-conditions.md#backflow-through-an-outlet)). Thus the engine does not use a temperature or a composition of a sink, and the registry does not list them.

The reservoir holds the **state** of the gas only. **How much** gas passes, and **through what**, stays on the endpoint: its set flow, its valve, its pump or its pressure ramp. So two feed steps can share one source and run at different flows.

## What an attached endpoint takes

An endpoint is a flow stream or a pressure boundary. While it is attached to a reservoir:

- Its temperature is the reservoir's temperature. A reservoir has one temperature, so an [inlet temperature profile](boundary-conditions.md#inlet-temperature-profile) is removed.
- Its mole fractions are the reservoir's mole fractions. (For an outlet on a sink, the engine does not use the temperature and the mole fractions.)
- Its pressure target is the reservoir's pressure, if it has one:
    - For a **pressure boundary**, this is the constant pressure or the end value of the ramp. The shape and the time of the ramp stay on the boundary.
    - For a **valve stream**, this is the pressure on the far side of the valve.
    - A **set-flow stream** and a **pump stream** have no pressure value. They do not use the reservoir's pressure.

The parameter registry lists these values one time, on the reservoir. It does not list them on the attached endpoints. A sweep or an optimiser that changes `reservoir[1].temperature` changes every endpoint attached to it.

An inlet attaches to a source. An outlet attaches to a sink. The engine refuses the other combination, and the deck checks report it before a run.

A [connection to another step](connecting-beds.md) or a [buffer tank](boundary-conditions.md#drawing-gas-from-another-step-buffer-tank) still decides the gas that an inlet receives. For such an inlet, the reservoir gives the first-cycle assumption, as the endpoint's own values do.

## In Python

<!-- doc-setup
import skarstrom
sim = skarstrom.examples.Perez2019_A.load()
-->

```python
feed = sim.add_source("Feed", temperature=298.15, pressure=1.0e5,
                      y={"CO2": 0.15, "N2": 0.85})
ambient = sim.add_sink("Ambient", pressure=1.0e5)       # a sink is a pressure

sim.steps["Feed"].bottom.attach(feed)             # the feed stream keeps its own flow
sim.steps["Pressurisation"].bottom.attach(feed)   # the ramp now ends at the source's pressure
sim.steps["Feed"].top.attach(ambient)

feed.temperature = 303.15                          # both inlets change
print(sim.reservoirs["Feed"].endpoints)
print(sim.steps["Feed"].bottom.reservoir.name)
```

`endpoint.detach()` (or `endpoint.attach(None)`) detaches an endpoint. The endpoint keeps the values that it had from the reservoir as its own values.

The structural edits are on `sim.edit`:

```python
rid = sim.edit.reservoir_from_endpoint(sim.steps["Evacuation"].bottom, "Vacuum")
sim.edit.rename_reservoir(rid, "Vacuum header")
sim.edit.remove_reservoir(rid)                    # attached endpoints keep the values
```

`reservoir_from_endpoint` makes a source (for an inlet) or a sink (for an outlet) from the values that the endpoint has now, and attaches the endpoint to it. Use it to turn the values of an existing deck into a shared source or sink. A set-flow or pump stream has no pressure. For such a stream, the reservoir takes the pressure target of the other end of the step. If that end has none, it takes the initial pressure.

## In C++

`adsim/objects/deck_edit.h` has the cascades: `AddReservoir`, `AddReservoirFromEndpoint`, `RemoveReservoir`, `RenameReservoir`, `AttachReservoir` and `DetachMismatchedReservoirs`. The builder (`adsim/objects/deck_builder.h`) has `build::AddReservoir` and `build::AttachReservoir`. The class is `adsim::Reservoir` (`adsim/boundary_conditions/reservoir.h`), and a cycle keeps its reservoirs in `Cycle::reservoirs_`.

## In the deck

The cycle stores the reservoirs (`Reservoirs`, each with an `ID`, a `Name`, a `Kind`, a `Temperature`, a `Pressure` and a `Mole_Fraction` map). Each endpoint stores the ID of its reservoir (`Reservoir_ID`, 0 for none). A deck from an earlier version has no reservoirs and loads as before. The [deck view](deck-view.md) gives the list as `reservoirs`, and each stream and pressure end gives its `reservoir` (an ID or null).

## Limits

- A reservoir is a fixed state. It has no holdup and no dynamics. A vessel that fills and empties is not a reservoir.
- The KPI configuration still names the feed, product and recycle endpoints. A reservoir does not set these roles.
