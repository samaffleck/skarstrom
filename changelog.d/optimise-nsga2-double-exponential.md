### Fixed — NSGA-II example uses the double-exponential evacuation pressure
- `examples/optimise_nsga2.py` and the `skarstrom.optimisation` docstring varied `pressure_exponential_end`, which the evacuation step of `Perez2019_A` no longer has. They now vary `pressure_double_exponential_end`. The Python wheel CI job installs `pymoo`, so this example test runs there too.
