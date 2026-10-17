# CONFIGURE_DEPENDS so a newly added test file is picked up by the next build
# rather than being silently left out until someone re-runs cmake by hand.
file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS "*.cpp")

# Add the test executable
add_executable(AdsorptionSimulatorTests ${TEST_SOURCES})

# Link the adsim library (includes are inherited from the target)
target_link_libraries(AdsorptionSimulatorTests PRIVATE adsim gtest_main)

# Auto-register tests with CTest.
# The working directory is adsim/, so relative paths resolve as follows:
#   data/...                digitised experimental curves (adsim/data/)
#   Results/...             test output, ~4 GB on a full run (git-ignored)
#   ../examples/decks/...   the example decks
#   ../data/adsorbents/...  the curated adsorbent library
#   ../data/cycles/...      cycle templates (written by test_cycle_designs.cpp)
include(GoogleTest)
# Paper-validation tests write one metrics JSON per case into
# <build>/validation (see test/fixtures/validation_metrics.h);
# scripts/docs/gen_validation_page.py merges them into
# adsim/test/papers/validation_summary.json and renders docs/validation.md.
gtest_discover_tests(AdsorptionSimulatorTests
    WORKING_DIRECTORY ${ADSIM_SOURCE_DIR}/adsim
    PROPERTIES ENVIRONMENT "ADSIM_VALIDATION_DIR=${ADSIM_BINARY_DIR}/validation"
    DISCOVERY_TIMEOUT 120
)

# Data file names (examples/decks, data/adsorbents, data/cycles) must be shell-safe (letters, digits and . _ - , = + @):
# downstream front ends copy the decks and adsorbent files as build assets with
# unquoted shell commands, and the Python wheel ships them as package data.
add_test(NAME DataFileNamesAreShellSafe
    COMMAND ${CMAKE_COMMAND}
            "-DASSET_DIRS=${ADSIM_SOURCE_DIR}/examples/decks;${ADSIM_SOURCE_DIR}/data/adsorbents;${ADSIM_SOURCE_DIR}/data/cycles"
            -P ${ADSIM_SOURCE_DIR}/cmake/check_data_names.cmake)
