# Getting Started

Skarstrom has two layers: the C++ engine (`adsim`), and the `skarstrom` Python package and command-line tool that wrap it. Most users only need the Python package.

---

## 1. Install the Python package

Install a pre-built wheel from PyPI (Linux x86_64/aarch64, macOS 13.3+ arm64/x86_64, Windows x64; Python 3.10–3.13):

```bash
pip install skarstrom
pip install "skarstrom[examples]"     # also pandas, matplotlib, scipy and Jupyter
```

### From source

For the development version, or on a platform without a wheel, build it from a checkout. You need Python 3.9 or newer, CMake 3.15 or newer, a C++20 compiler (GCC 11+, Clang 14+, Apple clang from Xcode 15+, or MSVC 2022) and, recommended, Ninja:

```bash
git clone --recursive https://github.com/samaffleck/skarstrom.git
cd skarstrom
pip install .                 # compiles the engine: several minutes the first time
```

Forgot `--recursive`? Run `git submodule update --init --recursive` before building. A source distribution (`skarstrom-<version>.tar.gz`, attached to every GitHub release) contains the submodule sources too, so `pip install skarstrom-<version>.tar.gz` builds offline.

Check the install:

```bash
skarstrom --version
skarstrom examples            # the bundled example decks
```

## 2. Run a deck from the command line

```bash
skarstrom run Perez2019_A.json --cycles 3
```

This runs three cycles of a bundled deck — the four-step CO₂/N₂ vacuum swing cycle of Estupiñan Perez et al. (2019), configuration A — and prints the run status, the cyclic-steady-state check and the KPIs of the last cycle. Drop `--cycles` to run the deck to cyclic steady state as saved. Other useful forms:

```bash
skarstrom examples run "Ward2022_Q100ccm" --json -q     # JSON summary on stdout
skarstrom run my_deck.json --set "step[Feed].time=120" --csv out/
skarstrom params Perez2019_A.json                           # the parameters you can set
```

## 3. Run a deck from Python

```python
import skarstrom
from skarstrom import ureg as u

sim = skarstrom.examples.Ward2022_Q100ccm.load()
result = sim.run()
print(result.status, f"{result.wall_time_s:.1f} s")

y = sim.series("mole_fraction", "CarbonDioxide", location="top")   # outlet breakthrough curve
print(y.time[-1], y.values[-1])

sim = skarstrom.examples.Perez2019_A.load()
sim.steps["Feed"].time = 2 * u.min        # a bare number would be seconds (SI)
result = sim.run()                      # to cyclic steady state
for k in sim.kpis(cycle=-1):
    print(k.component, f"purity {k.purity:.3f}", f"recovery {k.recovery:.3f}")
```

The Ward (2022) deck is a single breakthrough experiment: pure CO₂ fed at 100 cm³/min into a column of activated carbon that starts full of air. It runs in about a second. The outlet reaches 50 % CO₂ at about 1180 s, close to the measured curve (see [Validation](validation.md#ward-2022)).

[Python and the Command Line](python.md) covers the deck objects, units, parameter paths, results, the isotherm functions, the CLI and the bundled data.

## 4. Next steps

- Work through a tutorial: the [O₂ PSA (Arora 2021)](arora-2021-o2-psa.md), the [CO₂ VPSA with light-product pressurisation (Perez 2019)](perez-2019-vpsa.md) or the [TSA regeneration](tsa-regeneration.md).
- Browse the [Examples](examples.md): one runnable script per deck family in `examples/`, each compared with its paper.
- Read how the model works, starting from the [Governing Equations](equations.md).

---

## Building the C++ library and tests

The engine, its GoogleTest suite and the `bench_solver` timing tool build with CMake. The presets in `CMakePresets.json` cover Linux (GCC or Clang), macOS (Apple clang) and Windows (Visual Studio 2022):

```bash
git clone --recursive https://github.com/samaffleck/skarstrom.git
cd skarstrom
cmake --preset linux-release          # or macos-release, windows
cmake --build --preset linux-release  # (windows: --preset windows-release)
ctest --preset linux-release -j3
```

- The full suite (about 670 tests) includes the paper replications. It takes around 10 minutes on four cores and writes about 4 GB of results to `core/Results/`; delete that folder afterwards. `ctest -R` runs a subset, e.g. `ctest --preset linux-release -R "Isotherm|AdsorbentLibrary|Version"`.
- `build/linux-release/core/bench/bench_solver deck.json` solves any deck and prints the wall time and solver statistics. Add `--log medium` to see the per-cycle log, including the cycle at which cyclic steady state was reached.
- Without presets: `cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DADSIM_BUILD_TESTS=ON`.

[Using the C++ library](cpp-library.md) explains how to use the engine from your own C++ project, and [Numerics](numerics.md) what the solver does with a deck.
