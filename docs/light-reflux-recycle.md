# Light Reflux & Recycle

Recycling gas within a cycle is used to sharpen separations and lift purity or recovery. The most common form is **light reflux** (a light-product purge): part of the purified light product is sent back through a bed during its low-pressure purge step to displace the heavy component. Skarstrom supports two ways to model any recycle, depending on whether the *timing* of the recycled stream matters.

## Direct light reflux (connection)

Model reflux as a [connection](connecting-beds.md) between a producing step and a purge step. The purge inlet then follows the producing stream instant by instant (using the previous cycle's profile).

The two steps run for the **same duration**.

**Producing step** — on the light-product end, add a **Flow** boundary, direction **Outlet**, with a stream carrying the product to be refluxed (for example a set molar flow, or the full product stream).

**Purge / reflux step** — on the corresponding end of the bed being purged, add a **Flow** boundary, direction **Inlet**, with a stream. Enter the light-product composition as the first-cycle assumption.

**Connect them** (see [Connecting Beds](connecting-beds.md)). The connected purge inlet then copies the producing stream's **molar flow, temperature and composition**, so the bed is refluxed with real product gas.

Use this when the reflux stream's variation through the step is important — as it usually is for the sharp fronts in a purge step. This is how the built-in **6-Step VPSA with Rinse** example is wired.

## Buffer-tank averaging

Alternatively, an inlet can draw from a **buffer tank**: a model of a well-mixed surge tank that collects the gas produced by another step and feeds it back **time-averaged**. Use it when a physical buffer vessel sits between the beds, or when you only care about the average recycled amount rather than its instantaneous profile.

Configure it on an inlet under **Inlet Source**:

1. Tick **Use buffer tank from another step**.
2. Choose the **Interacting step** (the source) and its **Source endpoint**.
3. Set **Take % of last-cycle moles** — the fraction of the source endpoint's total produced moles to recycle.

The inlet then receives that fraction of the source's total moles, spread uniformly over the consumer step's duration, at the source's molar-flow-weighted **average temperature and composition** from its last completed run.

In Python, `inlet.set_buffer_tank(source, percentage=40 * u.percent)` sets the same thing (`source` is the producing endpoint, such as `sim.steps["Feed"].top`); see [the Python guide](python.md#pump-driven-streams-and-buffer-tanks).

> The buffer tank is a **time-averaging re-injection model**, not a vessel with its own volume or pressure dynamics. It and a live [interacting-step connection](connecting-beds.md) are mutually exclusive on the same endpoint — use one or the other. Until the source step has completed a run, the inlet falls back to the temperature, composition and flow you entered (the first-cycle assumption).

## Accounting for recycle in KPIs

Recycled gas should not be counted as net product. In the [KPI configuration](performance-metrics.md#configuring-kpis), give the inlets that receive recycled product the **recycle** role (`sim.kpi.track(..., recycle=[...])` in Python) so that recovery is computed as *(product − recycled product) / feed* while purity stays the purity of the collected product. See [Performance Metrics](performance-metrics.md) for the exact definitions.
