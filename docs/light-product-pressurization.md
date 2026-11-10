# Light Product Pressurisation

Light product pressurisation (LPP) re-pressurises a bed using purified **light product** (the weakly-adsorbed gas) taken from a bed that is producing. Because the pressurising gas is already clean, LPP raises product purity compared with feed-end pressurisation, and it recovers product that would otherwise be spent on pressurisation.

LPP is built as a [connection](connecting-beds.md) between a **producing step** and a **pressurisation step**.

## Setting it up

The two steps run for the **same duration**.

**Producing step** — on the light-product end (usually the Top), add a **Flow** boundary, direction **Outlet**, with a stream that carries the product. This is the gas that will be diverted to pressurise the other bed.

**Pressurisation step (LPP)** — on the same (light-product) end of the receiving bed:

1. Set the boundary to **Pressure**, direction **Inlet**.
2. Give it a **Linear** (or Exponential) pressure profile that ramps from the low pressure up to the high pressure over the step. This drives the pressurisation.
3. Enter a temperature and composition — these are used only as the [first-cycle assumption](connecting-beds.md).

**Connect them:** connect the producing outlet to the LPP inlet. Once connected, the LPP inlet takes its **temperature and composition from the producing stream**, so the bed is genuinely pressurised with light product rather than with the placeholder composition you typed.

## Notes

- The receiving bed follows the **pressure ramp** you set; the connection supplies the *quality* (composition and temperature) of the gas, while the pressure profile supplies the *driving force*. Match the ramp to the step time so the bed reaches the target pressure by the end of the step.
- LPP is included in the built-in **4-Step VPSA with LPP** example cycle — load it from **Cycle Design → Load Example Cycle** to see a complete working setup.
- To account for the light product consumed by pressurisation, classify the LPP step appropriately in the [KPI Configuration](performance-metrics.md) so recovery is computed correctly.
