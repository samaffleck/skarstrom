# Run health

The engine's tests check the engine's own decks. Nothing in them checks *your* deck while it runs: a deck can
be wrong in ways no test covers (too few cells, steps too short, a tolerance it cannot carry) and still finish
with status `completed` and plausible KPIs. Run-health checks are cheap tests the solver makes after every
cycle for the symptoms of such a run. They never stop a run, change its status, exit code or results: they add
warnings and set `healthy` to false.

```python
result = sim.run()
result.ok, result.healthy          # finished, and no run-health check fired
for w in result.warnings:          # RunHealthWarning
    print(w.code, w.value, w.threshold, w.cycle, w.step_name, w.suspects, w.hint)
sim.run_health().cycles            # the measures of every completed cycle
```

The CLI's `--json` summary has `healthy` and `run_health` (the same structured warnings; their messages are
also appended to `warnings`), the text summary has a `Health:` line, and run-health warnings are printed to
stderr even with `-q`. The solver logs each warning once at the end of the run, at log level *low*.

Each warning is `{code, value, threshold, cycle, step, suspects, hint}` (plus `step_name`, `component` and a
one-line `message`). Each code fires at most once per run, at the cycle where its value was largest. Each has
a switch, `settings.warn_<code in lower case>` (on by default); a switched-off check is not reported and does
not make the run unhealthy. The warning manifest (`skarstrom.warning_manifest()`) lists them with their
`code`.

## The checks

### `MOLE_FRACTION_SUM`

The largest $|\sum_i y_i - 1|$ over the cycle's logged times and the column's cells (packed and void cells).
Every component balance is solved for its own mole fraction, so nothing forces the sum to 1: it drifts when a
convective scheme or a source does not conserve the gas composition. A past advection defect let the sum
reach about 3 while the run still reported `completed`.

Suspects: the advection scheme (`settings.advection_scheme`), the grid, the step and outlet where the
largest error is. Try: re-run with `upwind` (or `weno`) and with more cells to compare; a sum that moves with
them is numerical, one that does not is in the model.

### `MASS_BALANCE_STUCK`

The cycle mole balance of each component (the `mass_balance_error` diagnostic of the KPIs, from the KPI
quadratures when they are on, else the trapezoid): the net moles that entered through every endpoint minus the
change of the column's gas and adsorbed inventory, over the moles that entered. It compares flows with the
inventory, so a conservative run closes it on *every* cycle up to the integration error, whether or not it is
near the cyclic steady state. One that stays open as the run settles is a model error (moles gained or lost),
not slow convergence. The rule, at cycle $k$, with the window $W$ = 5 cycles:

- the end-of-cycle state change at $k$ (the CSS state test's measure,
  [Cyclic steady state](cyclic-steady-state.md)) is at most $10^{-3}$: the run has settled;
- every $|e|$ in the window $k-W+1 \dots k$ is above the threshold; and
- the balance stopped falling: $(|e_k| / |e_{k-W+1}|)^{1/(W-1)} \ge 0.9$, it fell by less than 10 % per
  cycle over the window.

Needs the KPI quadratures (on by default, [Performance metrics](performance-metrics.md)): without them the column inventory is
not recorded and the check is skipped.

### `STEP_RETRIES`

The integrator's rejected attempts in a cycle (error-test failures, nonlinear convergence failures and the
solver's own step recoveries) per accepted internal step, over a cycle of at least 50 accepted steps. Many
rejections mean the integrator is fighting the deck: a tolerance it cannot meet, an abrupt boundary condition,
a grid too coarse for a front.

### `HISTORY_TRUNCATED`

The run outlived its log (`settings.max_cycles_to_keep`) and the log kept fewer cycles than the diagnostics
read back (6: the window of the CSS remaining-distance estimate). The CSS history, profiles and series read
after the run are then incomplete. KPIs are kept for every cycle, and the checks above run as each cycle
completes, so neither is affected.

## Thresholds

| Code | Value | Threshold |
|---|---|---|
| `MOLE_FRACTION_SUM` | max $\lvert\sum y - 1\rvert$ | $10^{-3}$ |
| `MASS_BALANCE_STUCK` | $\lvert e\rvert$ of the component | $2 \times 10^{-3}$ (with the state change $\le 10^{-3}$ and the rate $\ge 0.9$) |
| `STEP_RETRIES` | rejected attempts per accepted step | 0.5 |
| `HISTORY_TRUNCATED` | cycles logged | cycles needed (6, or the cycles run if fewer) |

Every example deck and paper case runs without a run-health warning; the test suite checks it for each deck it
regenerates.
