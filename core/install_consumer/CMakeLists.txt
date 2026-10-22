# Minimal downstream project that consumes an INSTALLED adsim through its CMake
# package — no add_subdirectory, no submodules. Used by
# scripts/test_install_consumer.sh to prove the install/export works:
#
#   cmake -S core/install_consumer -B <dir> -DCMAKE_PREFIX_PATH=<prefix> \
#         -DADSIM_REQUIRED_VERSION=<major.minor>
#   cmake --build <dir>
#   <dir>/adsim_consumer examples/decks/<deck>.json
#
# It is NOT part of the main build (and lives outside core/test/, whose
# sources are globbed into the gtest binary).
cmake_minimum_required(VERSION 3.15)
project(AdsimInstallConsumer LANGUAGES C CXX)

# The version to ask for comes from the caller (the script reads it from the
# root CMakeLists.txt), so a version bump needs no edit here. The package is
# exported with SameMinorVersion, so this also checks that compatibility rule.
set(ADSIM_REQUIRED_VERSION "" CACHE STRING
    "Adsim version to request in find_package (empty: any)")
find_package(Adsim ${ADSIM_REQUIRED_VERSION} CONFIG REQUIRED)
message(STATUS "Found Adsim ${Adsim_VERSION} in ${Adsim_DIR}")

add_executable(adsim_consumer main.cpp)
target_link_libraries(adsim_consumer PRIVATE Adsim::adsim)
