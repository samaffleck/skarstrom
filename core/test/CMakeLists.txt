# CONFIGURE_DEPENDS so a newly added test file is picked up by the next build
# rather than being silently left out until someone re-runs cmake by hand.
file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS "*.cpp")

# Add the test executable
add_executable(AdsorptionSimulatorTests ${TEST_SOURCES})

# Link the adsim library (includes are inherited from the target)
target_link_libraries(AdsorptionSimulatorTests PRIVATE adsim gtest_main)

# Precompiled header: the test sources include ../src/pch.h (cereal, Eigen,
# CoolProp, SUNDIALS, the STL) and GoogleTest, which every translation unit
# otherwise parses from scratch. The PCH is the test target's own rather than
# REUSE_FROM adsim, because adsim compiles with private definitions the tests
# do not share. A PCH is included before anything in a source file, so a
# #define placed above <cmath> no longer reaches it; _USE_MATH_DEFINES (M_PI on
# MSVC) is therefore a target definition. It is defined empty, exactly as the
# test files' own `#define _USE_MATH_DEFINES` lines define it, so those lines
# are identical redefinitions and raise no warning.
target_compile_definitions(AdsorptionSimulatorTests PRIVATE _USE_MATH_DEFINES=)
target_precompile_headers(AdsorptionSimulatorTests PRIVATE
    ${ADSIM_SOURCE_DIR}/core/src/pch.h
    <gtest/gtest.h>)
# As for adsim (core/CMakeLists.txt): no timestamp in a Clang PCH, so a
# compiler cache can reuse it on a fresh checkout.
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  target_compile_options(AdsorptionSimulatorTests PRIVATE "SHELL:-Xclang -fno-pch-timestamp")
endif()

# Auto-register tests with CTest.
# The working directory is core/, so relative paths resolve as follows:
#   data/...                digitised experimental curves (core/data/)
#   Results/...             test output, ~4 GB on a full run (git-ignored)
#   ../examples/decks/...   the example decks
#   ../data/adsorbents/...  the curated adsorbent library
#   ../data/cycles/...      cycle templates (written by test_cycle_designs.cpp)
include(GoogleTest)
# Paper-validation tests write one metrics JSON per case into
# <build>/validation (see test/fixtures/validation_metrics.h);
# scripts/docs/gen_validation_page.py merges them into
# core/test/papers/validation_summary.json and renders docs/validation.md.
#
# Test tiers (ctest labels). The per-PR CI jobs run everything except the
# labels `paper` and `slow` (`ctest -LE "paper|slow"`); main, the nightly run
# and manual runs run everything.
#   paper  the paper replications (core/test/papers/): agreement with
#          published results. A new paper suite goes in this list.
#   slow   verification that cannot be made cheap without losing its point
#          (a convergence order on fine grids, a long approach to CSS, a
#          regression that only shows on a fine grid or after many cycles).
#          Everything else should answer "does it run / does the feature
#          work" on a small grid and a few cycles (AGENTS.md, "Keep the
#          fast tier fast").
# Entries are GoogleTest filter patterns (Suite.Name, '*' wildcards).
set(ADSIM_PAPER_TESTS
    "Arora2021Test.*" "Casas2012Test.*" "Ferre2026Test.*" "Nguyen2022Test.*"
    "Nguyen2023Test.*" "Perez2019Test.*" "Streb2021Test.*" "Vysyaraju2024Test.*"
    "Ward2022Test.*" "WardPini2024Test.*")
set(ADSIM_SLOW_TESTS
    # Cycles to CSS with and without Aitken acceleration on a slow cycle, and
    # on the thermal four-step LPP deck (87 vs 39 cycles, ~10 s).
    "CSSAcceleration.ReducesCyclesToCSS*"
    "FourStepLPPTest.FourStepLPP_ThermalCSS_AccelerationReducesCycles"
    "CSSAcceleration.SlowGeometricModeReachesCssSooner"
    # Anderson acceleration on Perez 2019 A and E as shipped (~50 s).
    "CSSAndersonDecks.*"
    # The CSS remaining-distance estimate at the plain stop of Perez 2019 E
    # (205 cycles, ~90 s).
    "CSSDistanceDecks.Perez2019EReportsTheRemainingDistance"
    # CSS stops (plain, Aitken, Anderson) and a perturbed warm start against a
    # long reference run, full state and KPIs, on eight example decks
    # (1-8 min each).
    "*CSSReference.*"
    # The same decks from a clean, a feed-saturated, a hot and a cold bed and
    # from a neighbouring deck's CSS: all must reach the same CSS (2-10 min each).
    "*CSSStartIndependence.*"
    # Grid sequencing: coarse-to-fine warm starts against a cold fine grid on
    # six decks, three paths each (2-8 min per deck).
    "*GridSequencing.*"
    # CSS behaviour on a slow (rate 0.97) cycle: needs 100+ cycles.
    "CSSConvergence.SlowCycleFollowsTheGeometricModel"
    "CSSConvergence.AitkenAccelerationReachesTheSameKpis"
    # The stopping rule on Perez2019_A against its 1000-cycle KPIs (stops after
    # cycle 50; the old rule's early stop is at cycle 43).
    "CSSRobustDeck.*"
    # Four-grid KPI studies of three schemes on 10-80 cells (the upwind study,
    # CycleKpiGridConvergence.UpwindConvergesMonotonicallyAtFirstOrder, stays
    # in the fast tier; WENO's order is also ManufacturedSolutions.WenoOrder).
    "CycleKpiGridConvergence.HigherOrderSchemesConvergeFasterToTheSameKpis"
    # Regressions that only show on a fine grid or after many cycles: the
    # 320-cell stall, the stall at cycle 112, CSS without a max-step cap (83
    # cycles), a pressurisation pulse stepped over at cycle 47.
    "FineGridBreakthrough.Ward2022Q50ccmConvergesTo320Cells"
    "FourStepLPPTest.FourStepLPP_LongRunDoesNotStall"
    "FourStepLPPTest.FourStepLPP_NoMaxStepReachesCSS"
    "StepBoundaryRestart.VysyarajuPressurisationApproachesCSSSmoothly"
    # Cross-code regressions: frozen outlet curves of two other simulators,
    # replayed at the reference's grid with bounds set from the measured
    # difference (like the paper suites, an agreement check, not a smoke test).
    "CrossCodePeese.*"
    "CrossCodeToPSAil.*"
    # Cost and accuracy of the cycle-map Jacobian by forward sensitivities
    # (the measurement docs/cyclic-steady-state.md reports): 1-16 directions,
    # both correctors, finite differences at reltol 1e-10 (~1-2 min).
    "CycleSensitivity.MeasureCostAndAccuracy"
    # Knaebel-Hill Skarstrom cycle refined to 400 cells (the gap halves per
    # level and extrapolates to the theory): ~16 s.
    "SkarstromCycleEquilibriumTheory.GapVanishesUnderRefinement"
)
set(_adsim_test_props ENVIRONMENT "ADSIM_VALIDATION_DIR=${ADSIM_BINARY_DIR}/validation")
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.22)  # TEST_FILTER
  list(JOIN ADSIM_PAPER_TESTS ":" _adsim_paper_filter)
  list(JOIN ADSIM_SLOW_TESTS ":" _adsim_slow_filter)
  gtest_discover_tests(AdsorptionSimulatorTests
      WORKING_DIRECTORY ${ADSIM_SOURCE_DIR}/core
      TEST_FILTER "-${_adsim_paper_filter}:${_adsim_slow_filter}"
      PROPERTIES ${_adsim_test_props}
      DISCOVERY_TIMEOUT 120)
  gtest_discover_tests(AdsorptionSimulatorTests
      WORKING_DIRECTORY ${ADSIM_SOURCE_DIR}/core
      TEST_FILTER "${_adsim_paper_filter}"
      PROPERTIES ${_adsim_test_props} LABELS paper
      DISCOVERY_TIMEOUT 120)
  gtest_discover_tests(AdsorptionSimulatorTests
      WORKING_DIRECTORY ${ADSIM_SOURCE_DIR}/core
      TEST_FILTER "${_adsim_slow_filter}"
      PROPERTIES ${_adsim_test_props} LABELS slow
      DISCOVERY_TIMEOUT 120)
else()  # no labels: every test in one tier
  gtest_discover_tests(AdsorptionSimulatorTests
      WORKING_DIRECTORY ${ADSIM_SOURCE_DIR}/core
      PROPERTIES ${_adsim_test_props}
      DISCOVERY_TIMEOUT 120)
endif()

# Data file names (examples/decks, data/adsorbents, data/cycles) must be shell-safe (letters, digits and . _ - , = + @):
# downstream front ends copy the decks and adsorbent files as build assets with
# unquoted shell commands, and the Python wheel ships them as package data.
add_test(NAME DataFileNamesAreShellSafe
    COMMAND ${CMAKE_COMMAND}
            "-DASSET_DIRS=${ADSIM_SOURCE_DIR}/examples/decks;${ADSIM_SOURCE_DIR}/data/adsorbents;${ADSIM_SOURCE_DIR}/data/cycles"
            -P ${ADSIM_SOURCE_DIR}/cmake/check_data_names.cmake)
