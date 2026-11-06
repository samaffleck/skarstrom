### Fixed — memory grows without bound when a step takes very many solver steps
- The result log records a full snapshot of every variable at every accepted
  solver step (about 100 kB per solver step at 320 cells), so a stalled or
  very stiff step could exhaust memory before it failed. Each (cycle, step)
  now has a memory budget of 2 GiB, sized from the row size: when a step
  reaches it, every other interior row of that step is dropped (the first
  row and the newest row are always kept, so the end-of-step state used by
  the CSS check is exact), and the solver log says so once at the LOW level.
  The budget is above every bundled deck (the largest, Streb2021_A1, logs
  about 1.07 GB in its one step), so no result changes.
- The `max_log_points_per_step` setting is now documented as what it does:
  it decimates the profile CSV files and adds fill-in points in long solver
  steps; it does not limit the in-memory log.
