# Connecting Beds

Real PSA/TSA processes use several beds so that one can produce while another regenerates, and beds exchange gas with each other (equalisation, product pressurisation, reflux). Skarstrom captures these multi-bed interactions through **connections** between steps.

## How multi-bed processes are represented

Skarstrom solves **one bed** in detail, marching a single column through its sequence of [steps](cycles.md). A multi-bed process is modelled by laying all the relevant steps out on **one timeline** and **connecting** the steps that exchange gas. Each connection lets one step's boundary read the gas the *other* step produces.

Because the interacting steps are connected as a pair, a step that receives gas (say a pressurisation step) is linked to the step that supplies it (a producing step), and the simulator reproduces the coupled behaviour of the two beds without solving both columns simultaneously.

> **The coupling is lagged by one cycle.** A connected boundary uses the gas logged by its partner on the partner's **previous** run. This is exactly how a real multi-bed process behaves at cyclic steady state — each bed sees what the other bed did a cycle ago — and it converges to the correct coupled solution as the cycles repeat. On the **first** cycle there is no logged partner data yet, so each side falls back to the values you entered (the *first-cycle assumptions*). Expect the first cycle or two to be approximate; the solution becomes self-consistent as it approaches [cyclic steady state](cyclic-steady-state.md).

## Making a connection

Connections are stored in the deck with the cycle (each connected endpoint carries a pair id).

1. Give each of the two steps a boundary of the right kind (for example a valve stream on each, or a producing stream on one and a matching inlet on the other).
2. Connect the **From** step and endpoint to the **To** step and endpoint (in C++, `MakeMutualConnection`; front ends offer a connection editor; `adsim/objects/deck_edit.h` has the helpers that disconnect a pair cleanly).
3. Press **Connect**. The pair now appears under **Connected Endpoints**, where you can **Disconnect** it.

What passes through the connection depends on the endpoint types:

| Connected endpoint | Reads from its partner |
|--------------------|------------------------|
| Flow stream (no device) | Molar flow, temperature and composition |
| Valve stream | The partner's pressure, used as the driving pressure |
| Pump stream | The partner's pressure (suction/discharge) |
| Pressure boundary | Inlet temperature and composition |

### Rules

- **Interacting steps must have equal step time.** Two connected steps represent two beds acting at the same time, so their durations must match.
- A pressure boundary can connect to a flow/valve/pump stream, but **two pressure boundaries cannot be connected to each other** — one side must set the flow.
- An endpoint cannot connect to itself or to another endpoint of the same step.
- An endpoint can either be part of a live connection **or** use a [buffer-tank source](light-reflux-recycle.md#buffer-tank-averaging) — not both.

## Two ways to couple

Skarstrom offers two coupling mechanisms, suited to different purposes:

- **Interacting-step connections** (this page) reproduce the *instantaneous* profile of the partner stream through the step, using the partner's previous-cycle time series. Use these when the timing of the exchange matters — equalisation, direct product pressurisation, direct reflux.
- **Buffer-tank sources** feed an inlet with the *time-averaged* gas (total moles, average temperature and composition) produced by another step's endpoint over its last run. This models a well-mixed surge tank between the two beds. See [Light Reflux & Recycle](light-reflux-recycle.md).

---

The following pages show these tools applied to specific process steps:

- [Pressure Equalisation](pressure-equalization.md)
- [Light Product Pressurisation](light-product-pressurization.md)
- [Light Reflux & Recycle](light-reflux-recycle.md)
