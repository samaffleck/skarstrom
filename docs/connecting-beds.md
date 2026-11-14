# Connecting Beds

Real PSA/TSA processes use several beds so that one can produce while another regenerates, and beds exchange gas with each other (equalisation, product pressurisation, reflux). Skarstrom captures these multi-bed interactions through **connections** between steps.

## How multi-bed processes are represented

Skarstrom solves **one bed** in detail, marching a single column through its sequence of [steps](cycles.md). A multi-bed process is modelled by laying all the relevant steps out on **one timeline** and **connecting** the steps that exchange gas. Each connection lets one step's boundary read the gas the *other* step produces.

Because the interacting steps are connected as a pair, a step that receives gas (say a pressurisation step) is linked to the step that supplies it (a producing step), and the simulator reproduces the coupled behaviour of the two beds without solving both columns simultaneously.

> **The coupling is lagged.** A connected boundary reads its partner's log from the partner's **most recent run**, linearly interpolated at the same time since the start of the step. Each step clears its own log when it starts. So the lag depends on the order of the two steps in the cycle: when the partner step comes **later** in the cycle, its most recent run is from the previous cycle (a lag of one cycle); when it comes **earlier**, its most recent run is from the current cycle. Until the partner has run once, there is no logged partner data, so that side falls back to the values you entered (the *first-cycle assumptions*); this happens only in the first cycle, and only for the side whose partner comes later. At cyclic steady state each run repeats the previous one, so the lagged coupling gives the coupled solution. Expect the first cycle or two to be approximate; the solution becomes self-consistent as it approaches [cyclic steady state](cyclic-steady-state.md).
>
> **Moles are conserved across a connection only at cyclic steady state.** Each side computes its own flow from the other side's log, not from the other bed's current state. A connected valve's flow uses the partner's logged pressure, and the two valves of a pair give equal and opposite flows only when the logs agree with the current runs. A Pressure inlet takes only the temperature and composition from its partner: its flow follows its own pressure profile, not the partner's outflow. A *None* stream copies the partner's logged flow, which equals the partner's current flow only when the runs repeat. Before cyclic steady state, the gas that enters one bed through a connection is in general not the gas that left the other.

!!! note "New in 0.9"
    The warning `settings.warn_connection_pressure` reports gas that passes through a connection from the lower to the higher pressure.

### Pressure order through a connection

A connected endpoint replays its partner's previous run: a flow stream replays the partner's flow, and a pressure boundary follows its own profile. Neither checks the pressures of the two beds, so gas can pass from the lower to the higher pressure, for example when a producing bed falls below the bed that its light product pressurises. The engine reports this but does not prevent it.

At the end of a run, the warning **Connection pressure** (`settings.warn_connection_pressure`) reports gas that enters a bed through a connection while that bed's end pressure is more than 1 kPa above the partner's logged pressure at the same instant, by more than 0.1 % of the cycle feed in one step, in the last cycle. Pump streams are not checked. A valve stream on a connection sets its flow from the pressure difference, so it passes gas from the higher to the lower pressure by its own law.

## Making a connection

Connections are stored in the deck with the cycle (each connected endpoint carries a pair id).

1. Give each of the two steps a boundary of the right kind (for example a valve stream on each, or a producing stream on one and a matching inlet on the other).
2. Connect the **From** step and endpoint to the **To** step and endpoint (in C++, `MakeMutualConnection`; front ends offer a connection editor; `adsim/objects/deck_edit.h` has the helpers that disconnect a pair cleanly).
3. Press **Connect**. The pair now appears under **Connected Endpoints**, where you can **Disconnect** it.

What passes through the connection depends on the endpoint types:

| Connected endpoint | Reads from its partner |
|--------------------|------------------------|
| Flow stream (no device) | Molar flow (with its sign), temperature and composition |
| Valve stream | The partner's pressure, used as the driving pressure (the bed-end face pressure, or the end-cell pressure with `settings.valve_bed_pressure` = `end_cell`); temperature and composition at an inlet |
| Pump stream | The partner's pressure (suction/discharge); temperature and composition at an inlet |
| Pressure boundary | Inlet temperature and composition |

### Rules

- **Interacting steps must have equal step time.** Two connected steps represent two beds acting at the same time, so their durations must match.
- A pressure boundary can connect to a flow/valve/pump stream, but **two pressure boundaries cannot be connected to each other** — one side must set the flow.
- An endpoint cannot connect to itself or to another endpoint of the same step.
- An endpoint can either be part of a live connection **or** use a [buffer-tank source](light-reflux-recycle.md#buffer-tank-averaging) — not both.

## Two ways to couple

Skarstrom offers two coupling mechanisms, suited to different purposes:

- **Interacting-step connections** (this page) reproduce the *instantaneous* profile of the partner stream through the step, using the time series of the partner's most recent run. Use these when the timing of the exchange matters — equalisation, direct product pressurisation, direct reflux.
- **Buffer-tank sources** feed an inlet with the *time-averaged* gas (total moles, average temperature and composition) produced by another step's endpoint over its last run. This models a well-mixed surge tank between the two beds. See [Light Reflux & Recycle](light-reflux-recycle.md).

---

The following pages show these tools applied to specific process steps:

- [Pressure Equalisation](pressure-equalization.md)
- [Light Product Pressurisation](light-product-pressurization.md)
- [Light Reflux & Recycle](light-reflux-recycle.md)
