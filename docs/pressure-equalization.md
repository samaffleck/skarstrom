# Pressure Equalisation

Pressure equalisation recovers work and gas by connecting a bed that is finishing a high-pressure step to a bed that is at low pressure, so gas flows from one to the other until their pressures approach each other. It improves recovery and reduces compression duty, and is a standard feature of multi-bed PSA cycles.

There is no dedicated "equalisation" object — you build it from a **valve on each bed** [connected together](connecting-beds.md). The valve turns the flow into a pressure-driven exchange: gas moves from the higher-pressure endpoint to the lower-pressure one at a rate set by the valve coefficient.

## Setting it up

You need two steps, running for the **same duration** (they represent the two beds equalising at the same time):

**Equalisation-down step** (the bed releasing gas):

1. On the end that connects to the other bed, set the boundary to **Flow**, direction **Outlet**, and add a stream with **Boundary Device Type = Valve**.
2. Enter a **Valve Flow Coefficient (Cv)** and a target pressure.
3. Set the opposite end to **Closed**.

**Equalisation-up step** (the bed receiving gas):

1. Mirror the setup: a **Valve** stream, direction **Inlet**, on the connecting end.
2. Set the opposite end to **Closed**.

**Connect them:** connect the two valve endpoints. Once connected, each valve uses the *other* bed's pressure as its driving pressure, so gas flows down the pressure gradient between them.

## Tuning

- The **valve coefficient (Cv)** controls how quickly the beds equalise. A small Cv throttles the exchange so the pressures approach gently over the step; a large Cv equalises almost immediately. Choose Cv together with the step time so the beds reach the intended intermediate pressure by the end of the step.
- For multiple equalisation stages, repeat the pattern with additional step pairs, each equalising to a successively lower pressure level.

> Remember the coupling is [lagged by one cycle](connecting-beds.md): the two beds converge on the true shared equalisation pressure as the simulation approaches [cyclic steady state](cyclic-steady-state.md), not on the very first cycle.
