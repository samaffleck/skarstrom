### Changed — the log keeps one cycle by default
- `settings.max_cycles_to_keep` defaults to 1 (was 10), and every example deck keeps 1. The log stores
  every solver step of each kept cycle, and with 10 cycles some example decks (Arora2021_LiLSX_PSA,
  Ferre2026_NaYF) ran out of the 2 GB that the web solver has. The KPIs and the CSS history still cover
  every cycle, and the log still keeps the cycles the CSS check reads. Raise the setting only to plot the
  series or profiles of earlier cycles.
- There is no run-health warning when the log keeps fewer cycles than the run: the `HISTORY_TRUNCATED`
  check and its switch were removed before release. `WarningSettings` is class version 4; a deck saved
  with that switch (version 2 or 3) still loads, and the value is ignored.
