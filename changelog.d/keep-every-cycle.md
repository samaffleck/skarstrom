### Changed — every cycle's KPIs and CSS record are kept
- `settings.max_cycles_to_keep` now limits only the logged series and profiles (the memory-heavy data).
  `sim.kpis()`, `sim.kpi_table()`, `sim.css_history()` and the CLI's `kpis.csv` / `css.csv` cover every
  cycle of the run: the solver records each cycle's CSS check (`CSSHistory`, `GetCSSHistory()`) when the
  cycle completes, as it already did for the KPIs. One record is about 200 bytes, so a 2000-cycle run adds
  well under 1 MB. Results are unchanged.
- The `HISTORY_TRUNCATED` run-health warning now says only that series and profiles cover the logged cycles.
