# Dynamic Tanks

A **dynamic tank** is a well-mixed gas vessel of fixed volume that is solved together with the bed: its species moles and its temperature are unknowns of the same equation system as the bed's, in every step. A bed end can draw gas from a tank or send gas into it, and the tank can have its own streams to the outside (a make-up feed, a vent valve, a back-pressure regulator or a pressure hold). Use it for a product or surge tank on a single-bed rig, a feed tank, or any vessel whose pressure and composition matter to the cycle.

A tank is not the [buffer-tank source](light-reflux-recycle.md#buffer-tank-averaging), which re-injects a share of another step's time-averaged output and has no volume or pressure of its own. Both remain available.

> **One bed per tank.** A tank shared by several staggered beds (the `beds` field) is not supported yet: `beds` must be 1. The tank is then exact at all times, not only at the cyclic steady state.

## Equations

With $n_i$ the moles of component $i$, $N = \sum_i n_i$ and $V$ the volume, the pressure follows from the ideal gas law and the composition from the moles:

$$
P = \frac{N R T}{V}, \qquad y_i = \frac{n_i}{N}.
$$

The balances are

$$
\frac{dn_i}{dt} = \sum_p F_{p,i},
\qquad
\frac{d}{dt}\Big(\sum_i n_i c_{v,i} T\Big) = \sum_p \sum_i F_{p,i}\, c_{p,i}\, T_p + Q,
$$

where $F_{p,i}$ is the molar flow of component $i$ *into* the tank through port or stream $p$ (negative when gas leaves), $T_p$ the temperature of the gas that crosses $p$ (the tank's own for gas leaving it), $c_{p,i}$ the component's molar heat capacity (the same value the bed uses) and $c_{v,i} = c_{p,i} - R$. Written this way the energy balance is exactly conservative; it is the same as

$$
N c_v \frac{dT}{dt} = \sum_\text{in} F c_p (T_\text{in} - T) + R T \frac{dN}{dt} + Q .
$$

A port of the bed takes exactly the flow the bed's boundary passes (the stream's own flow function, or the pressure boundary's face flux), so whatever leaves the bed enters the tank and the two together conserve moles, including when a flow reverses.

### Thermal modes

| Mode | Heat to the gas $Q$ | Unknowns |
|---|---|---|
| `isothermal` (default) | the temperature is held at the initial temperature | $n_i$ |
| `adiabatic` | $0$ | $n_i$, $T$ |
| `heat_loss` | $UA\,(T_\text{amb} - T)$ | $n_i$, $T$ |
| `wall` | $h_\text{in} A (T_w - T)$, with a wall node $C_w \, dT_w/dt = h_\text{in} A (T - T_w) - h_\text{out} A (T_w - T_\text{amb})$ | $n_i$, $T$, $T_w$ |

In an [isothermal run](energy-balance.md#isothermal-mode) every tank is held at the run's isothermal temperature, as the bed is.

## Ports: bed ends that use a tank

Any bed endpoint (a set-flow stream, a valve stream, a pump stream or a pressure boundary) can name a tank. Its own direction and bed end say which way the gas goes; it then sees the tank as follows, in the steps it belongs to:

| Endpoint | Inlet (tank → bed) | Outlet (bed → tank) |
|---|---|---|
| Valve stream | upstream pressure, composition and temperature are the tank's | downstream pressure is the tank's |
| Set-flow stream | the set flow at the tank's composition and temperature | the set flow goes into the tank |
| Pump stream | the compressor draws from the tank pressure | the vacuum pump discharges against the tank pressure |
| Pressure boundary | bed end held at its set pressure (an ideal regulator) with the tank's gas | the outflow goes into the tank |

A tank port excludes a connection and a buffer-tank source on the same endpoint. Gas that flows back through an *outlet* (for example a valve to a tank whose pressure has risen above the bed's) carries the bed-end composition and temperature, as for every outlet. Several steps of the cycle can use the same tank, and a step can feed and drain it at the same time.

A tank is placed in the solver's equation layout next to the bed end its ports use. A tank used at **both** ends (a recycle around the bed) is correct but makes the solver's banded matrix span the whole system, which is slower; the solver prints an advisory.

## External streams

| Kind | Flow out of the tank | Fields |
|---|---|---|
| `make_up` | $-F$ (inflow) with $F$ set (`mode = "set_flow"`), or $F = C_v \max(P_h - P, 0)/10^5$ from a header (`mode = "valve"`) | `molar_flow`, `cv`, `header_pressure`, `y`, `temperature` |
| `valve` | $C_v (P - P_d)/10^5$; with `check_valve` (default) never negative | `cv`, `downstream_pressure`, `check_valve`; `y`, `temperature` of the downstream gas if it flows back |
| `regulator` | $C_v \,\mathrm{clamp}\!\big(\tfrac{P - P_\text{set}}{P_\text{full} - P_\text{set}}, 0, 1\big) \max(P - P_d, 0)/10^5$ | `cv`, `set_pressure`, `full_open_pressure`, `downstream_pressure` |
| `hold` | a regulator with $P_\text{full} = P_\text{set}(1 + \delta)$, $\delta$ = `band` (default $10^{-4}$) | `set_pressure`, `band`, `cv`, `downstream_pressure` |

$C_v$ is in mol/(s·bar). The `max` and `clamp` are smoothed over 1 Pa and 1 % of the regulator's band. The **hold** keeps the tank at its set pressure (to within $\delta P_\text{set}$) and draws off whatever net inflow arrives while that is positive and below $C_v (P_\text{full} - P_d)/10^5$: an overflow product tank.

Every external stream is logged like a bed stream (molar flow positive *out of the tank*) under the name `"<tank> / <stream>"`, in every step, and its ID is an endpoint ID (`endpoint[<id>]` paths). The KPIs are still configured on bed endpoints: at the cyclic steady state a tank's inventory repeats, so over a cycle what the bed sends into a product tank equals what the tank's draw-off delivers.

## Cyclic steady state

A large tank approaches the cyclic steady state slowly, and the bed alone can pass the [CSS test](cyclic-steady-state.md) while the tank still drifts. So the tanks join every method:

- the **state test** adds each tank's end-of-cycle relative change of total moles and temperature and the absolute change of its mole fractions;
- the **mass-balance methods** require every tank's cycle balance, $|\text{in} - \text{out}|/\text{in}$ over its ports and streams, to pass as well as the bed's;
- the **molar-profile method** compares the tanks' concentrations $n_i/V$ along with the bed's profile;
- **Aitken acceleration** extrapolates the tank moles and temperatures with the bed's slow fields.

## Results

The equation log holds each tank's pressure, temperature, wall temperature, total moles, heat and, per component, its mole fraction, moles and in- and outflow, at every logged time. `save_logs` writes them to `tanks.csv`.

## Building a tank

In Python:

```python
import skarstrom as sk
from skarstrom import ureg as u

tank = sim.add_tank("Product", volume=0.0314, pressure=11.12e5, temperature=303.15,
                    y={"H2": 0.999, "CO": 0.001})
tank.set_thermal("heat_loss", ua=2.0, ambient_temperature=298.15)
tank.add_stream("hold", "Product", set_pressure=11.12e5, cv=10, downstream_pressure=1e5)
sim.steps["Feed"].top["To tank"].set_tank(tank)          # the feed step's product goes into the tank
sim.steps["Purge"].top["Purge in"].set_tank(tank)        # and the purge draws from it
result = sim.run()
p = sim.tank_series(tank, "pressure")                    # .time from the cycle start, .values [Pa]
y = sim.tank_series(tank, "mole_fraction", "CO")
```

In C++ (`adsim/objects/deck_builder.h`): `AddTank(data, name)` returns the tank ID; `SetTankThermal(data, tank, kind)`; `AddTankStream(data, tank, kind, name)` returns the stream's endpoint ID; `SetTankMakeUpMode(data, endpoint, mode)`; `SetTankPort(data, endpoint, tank)` (0 removes the port); `RemoveTank`, `RemoveTankStream`. The values are set through the [parameter registry](parameterisation.md): `tank[<id>].volume`, `initial_pressure`, `initial_temperature`, `y[<component>]`, `beds`, `thermal_mode`, `ambient_temperature`, `ua`, `wall_area`, `h_in`, `h_out`, `wall_heat_capacity`, `initial_wall_temperature`; a stream's fields are `endpoint[<id>].<field>` (aliases `tank[<name>]` and `tank[<name>].stream[<name>]`).

Decks with tanks do not load in engines older than this feature; older decks load unchanged, with no tanks.
