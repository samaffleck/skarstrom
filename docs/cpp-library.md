# Using the C++ Library

The engine is the static library `adsim` (CMake target `Adsim::adsim`), written in C++20. Its public headers are in `adsim/include/adsim/`; the most useful entry points are:

| Header | What it gives you |
|---|---|
| `adsim/simulation.h` | `adsim::Simulation`: load a deck, get/set parameters, run (blocking or in slices), read KPIs, the CSS history, time series and axial profiles. The Python package wraps exactly this class. |
| `adsim/objects/user_data.h` | `adsim::UserData`, the whole deck (fluid, adsorbents, bed, cycle, settings), and `LoadFromJson` / `SaveToJson` / `SaveToString` / `LoadFromString`. |
| `adsim/objects/deck_io.h` | `LoadDeckFlexible`: load deck JSON text with either root shape (`"simulation"` or cereal's default `"value0"`). |
| `adsim/objects/deck_edit.h` | Structural edits that keep a deck consistent: add/remove a component (with every mole-fraction map and adsorbent), add/remove adsorbents, add a library adsorbent, disconnect connection pairs. |
| `adsim/parameters/parameter_registry.h` | The parameter registry: every settable deck value by stable path (`step[2].time`, `adsorbent[1].isotherm[CarbonDioxide].b1`, ...). |
| `adsim/objects/adsorbent_library.h` | Read the curated adsorbent catalogue (`data/adsorbents/library.json`). |
| `adsim/equations/isotherm_kinds.h` | Isotherm kind tags (`"sips"`, `"dual_site_langmuir_p"`, ...), generic parameter access and single-component evaluation. |
| `adsim/solver/isotherm_fit.h` | The isotherm fitter (global MaxLIPO search plus bounded Levenberg–Marquardt). |
| `adsim/solver/coupled_solver.h` | `CoupledSolver` (blocking cycle loop) and `CoupledSolverChuncked` (time-sliced, for single-threaded hosts). |
| `adsim/solver/kpis.h`, `adsim/solver/kpi_config.h` | KPI calculation and the deck's KPI configuration. |
| `adsim/solver/optimisation.h`, `adsim/solver/parameterisation.h` | dlib-based optimiser and the parameter-fitting helpers. |
| `adsim/version.h` | Generated version macros (`ADSIM_VERSION_STRING`, `adsim::kVersionString`). |

```cpp
#include "adsim/simulation.h"

int main() {
  adsim::Simulation sim = adsim::Simulation::FromFile("examples/decks/Perez2019_A.json");
  sim.SetParameter("step[Feed].time", 120.0);
  adsim::RunResult r = sim.Run();
  for (const auto& k : sim.KPIs(-1)) std::printf("%s purity %.3f\n", k.component.c_str(), k.purity);
  return r.Ok() ? 0 : 1;
}
```

## As an installed package (`find_package`)

```bash
cmake --preset linux-release -DCMAKE_INSTALL_PREFIX=$HOME/.local
cmake --build --preset linux-release
cmake --install build/linux-release
```

```cmake
find_package(Adsim 0.1 CONFIG REQUIRED)
target_link_libraries(my_app PRIVATE Adsim::adsim)
```

The install tree contains the static libraries (adsim, SUNDIALS, muparser, dlib, CoolProp), the public headers with the third-party headers they include (under `include/adsim/third_party`), the CMake package files and the licence texts. `scripts/test_install_consumer.sh` builds `adsim/install_consumer/` against a fresh install to prove it works.

## As a subproject (`add_subdirectory`)

A front end that pins the engine as a git submodule (or a release tarball) builds it in-tree:

```cmake
add_subdirectory(external/skarstrom)          # this repository
target_link_libraries(my_frontend PRIVATE Adsim::adsim)
```

- The options default to what a front end wants: `ADSIM_BUILD_TESTS`, `ADSIM_BUILD_BENCH` and `ADSIM_INSTALL` are `OFF` when the engine is not the top-level project. Set them before `add_subdirectory` to override.
- The engine's settings are directory-scoped (C++20, the `FMT_CONSTEVAL=` and `CMAKE_POLICY_VERSION_MINIMUM` workarounds for CoolProp) and do not change the parent's. The dependency options it sets (`BUILD_SHARED_LIBS`, the SUNDIALS and dlib switches) are cache entries without `FORCE`, so a value the parent already set wins.
- `Adsim::adsim` carries the include directories, `EIGEN_MPL2_ONLY` and `FMT_CONSTEVAL=` to your targets.
- The engine's targets (`adsim`, `CoolProp`, `dlib`, `muparser`, the `sundials_*` libraries) must not clash with targets of the same name in the parent.
- The bundled data are at `examples/decks/`, `data/adsorbents/` and `data/cycles/` of the checkout, for front ends that ship them as assets. File names there are shell-safe (checked by the `DataFileNamesAreShellSafe` test).

`adsim/subproject_consumer/` is a minimal front end built this way, and `scripts/test_subproject_consumer.sh` builds and runs it.

## WebAssembly (Emscripten)

The engine compiles with Emscripten, as a subproject of a front end configured with `emcmake`. The engine's own code is compiled with `-fexceptions` (it throws and catches C++ exceptions); the front end must link its module with `-fexceptions` (or `-sNO_DISABLE_EXCEPTION_CATCHING=1`) and chooses everything else (threads, memory growth, the file system, exported functions). For a browser, drive the solver with `CoupledSolverChuncked` (or `Simulation::Start` / `Advance(max_time_ms)`) so the host stays responsive; see [Numerics](numerics.md#the-chunked-solver-for-cooperative-front-ends).

```bash
source /path/to/emsdk/emsdk_env.sh
bash scripts/test_subproject_consumer.sh --emscripten   # builds adsim + a tiny module, runs it under Node.js
```

CI builds the engine this way on every push.

## Deck files and versioning

A deck is the JSON written by cereal from `UserData`. Every serialised class carries a `cereal_class_version`; loading an older deck runs the class's migration code for its version, so old decks keep loading. A deck written by a newer engine than the reader's may not load. Bump the class version (`CEREAL_CLASS_VERSION`) whenever you add or change a serialised field, and handle the old versions in `serialize` / `load`.
