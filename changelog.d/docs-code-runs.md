### Added — documentation code blocks run in CI
- `scripts/docs/run_doc_blocks.py` executes every Python and shell code block of `docs/*.md` and
  `AGENTS.md` against the installed package and reports each block's file, line, result and time; CI runs
  it in the Python wheel job (about 35 s, each run capped at 3 cycles). `--full` runs the uncapped set,
  including blocks marked `<!-- long -->`, for a nightly run. `<!-- no-run -->` excludes a block and a
  hidden `<!-- doc-setup -->` block defines what a fragment assumes; the convention is in AGENTS.md.

### Fixed — documentation examples the engine rejected
- A call to `sim.edit.set_pump_kind` on a valve stream, `sim.stream_series` with a wrong argument list and
  an `evaluate_iast` call with a pressure array (AGENTS.md); a `fit_isotherm` example whose data were
  overwritten by the example before it, and one that wrote a fit of one isotherm kind into a deck of
  another (isotherm fitting); `enumerate(sim.steps)` formatted as a name (TSA regeneration); a variable
  and pressure path the example deck does not have (global variables); the evacuation pressure path of the
  NSGA-II example (parameterisation, Perez 2019); and several fragments that used
  names no earlier block defined.
