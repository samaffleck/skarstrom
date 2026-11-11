cmake_minimum_required(VERSION 3.15)
# Single source of truth for the Skarstrom version (SemVer MAJOR.MINOR.PATCH).
# It flows into the generated adsim/version.h (ADSIM_VERSION_STRING), the CMake
# package version (find_package(Adsim)) and the Python package version (read
# from these two lines by pyproject.toml; keep them adjacent).
# ADSIM_VERSION_SUFFIX marks a pre-release ("rc1", "a1", "b2"; empty for a
# final release) and/or a development version (".dev0"). Between releases
# main is a development version (PEP 440): after releasing 0.9.0 it is 0.9.1
# with ".dev0", so a build of main reports 0.9.1.dev0, never a released
# version. The release workflow sets and strips the marker
# (scripts/release.py); do not edit these lines by hand. The suffix is
# appended to ADSIM_VERSION_STRING and the Python version, and a released one
# must match the tag (v0.2.0rc1). See RELEASING.md.
set(ADSIM_PROJECT_VERSION 0.9.2)
set(ADSIM_VERSION_SUFFIX ".dev0")

# Minimum macOS version. It must be set before project(). 13.3 is the lowest
# release whose system libc++ ships floating-point std::to_chars, which
# core/src/solver/log.cpp uses; std::filesystem only needs 10.15. Without a
# value the compiler targets the build machine's own macOS version, so a binary
# built on macOS 15 would refuse to start on 14.
# Pass -DCMAKE_OSX_DEPLOYMENT_TARGET=... or set MACOSX_DEPLOYMENT_TARGET to
# override (the Python wheels set the same 13.3 in pyproject.toml).
# APPLE is not defined until project() runs, hence the host check.
set(ADSIM_MIN_MACOS_VERSION 13.3)
if(CMAKE_HOST_APPLE AND NOT CMAKE_SYSTEM_NAME
   AND NOT DEFINED CMAKE_OSX_DEPLOYMENT_TARGET
   AND NOT DEFINED ENV{MACOSX_DEPLOYMENT_TARGET})
    set(CMAKE_OSX_DEPLOYMENT_TARGET ${ADSIM_MIN_MACOS_VERSION} CACHE STRING
        "Minimum macOS version to target")
endif()

project(AdSim VERSION ${ADSIM_PROJECT_VERSION} LANGUAGES C CXX)

# Development builds name their commit: ADSIM_VERSION_LOCAL is a PEP 440 local
# version ("+g0a19084", from `git rev-parse --short HEAD` at configure time)
# appended to ADSIM_VERSION_STRING, so adsim::kVersionString, skarstrom
# --version and skarstrom.__version__ read 0.9.1.dev0+g0a19084. Only a
# development version (".devN" suffix) built from a git checkout gets one: a
# released version (the release build at its tag) is always clean, and a
# development build without git (a source tarball) has no local part. The
# commit is read when CMake configures, so re-run CMake to refresh it. Pass
# -DADSIM_VERSION_LOCAL=... to override (empty for none).
if(NOT DEFINED ADSIM_VERSION_LOCAL)
    set(ADSIM_VERSION_LOCAL "")
    if(ADSIM_VERSION_SUFFIX MATCHES "\\.dev[0-9]+$" AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
        find_package(Git QUIET)
        if(GIT_FOUND)
            execute_process(
                COMMAND "${GIT_EXECUTABLE}" rev-parse --short HEAD
                WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                OUTPUT_VARIABLE _adsim_git_sha
                RESULT_VARIABLE _adsim_git_rc
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET)
            if(_adsim_git_rc EQUAL 0 AND _adsim_git_sha MATCHES "^[0-9a-f]+$")
                set(ADSIM_VERSION_LOCAL "+g${_adsim_git_sha}")
            endif()
        endif()
    endif()
endif()
set(ADSIM_VERSION_STRING "${ADSIM_PROJECT_VERSION}${ADSIM_VERSION_SUFFIX}${ADSIM_VERSION_LOCAL}")
message(STATUS "Skarstrom (adsim) version ${ADSIM_VERSION_STRING}")

# Paths of this project, valid whether it is the top-level project or added to
# a front end with add_subdirectory(). Never use CMAKE_SOURCE_DIR /
# CMAKE_BINARY_DIR in this tree: they point at the parent project.
set(ADSIM_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(ADSIM_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(ADSIM_IS_TOP_LEVEL ON)
else()
    set(ADSIM_IS_TOP_LEVEL OFF)
endif()

if(APPLE AND CMAKE_OSX_DEPLOYMENT_TARGET
   AND CMAKE_OSX_DEPLOYMENT_TARGET VERSION_LESS ADSIM_MIN_MACOS_VERSION)
    message(WARNING
        "CMAKE_OSX_DEPLOYMENT_TARGET is ${CMAKE_OSX_DEPLOYMENT_TARGET}, but Skarstrom needs "
        "macOS ${ADSIM_MIN_MACOS_VERSION} or later (floating-point std::to_chars). "
        "Expect 'is unavailable: introduced in macOS 13.3' compile errors.")
endif()

# These are directory-scoped: when this project is a subdirectory they apply to
# the engine and its dependencies only, not to the parent project.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(ADSIM_IS_TOP_LEVEL)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

# Allow absolute install paths (needed for subdirectories with bad install rules)
set(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION OFF)
set(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION ON)

# ============================================================
# Options
# ============================================================
# Every option defaults to the right value for both uses: a top-level build of
# the engine, and a front end that does add_subdirectory(<this repo>). A front
# end normally needs none of them; it links the `adsim` target (alias
# Adsim::adsim).
option(ADSIM_BUILD_TESTS "Build the GoogleTest suite" OFF)
option(ADSIM_BUILD_BENCH "Build the bench_solver timing tool" ${ADSIM_BUILD_TESTS})
# The Python extension module (skarstrom._core). Normally driven by
# scikit-build-core from pyproject.toml (pip install .), which turns this on
# and the tests and install rules off.
option(ADSIM_BUILD_PYTHON "Build the skarstrom Python extension module (pybind11)" OFF)
# Install rules + CMake package config for the adsim library, so another project
# can use find_package(Adsim) instead of add_subdirectory().
# Off by default when the engine is a subproject (the parent decides what it
# installs) and for the Python build.
if(ADSIM_BUILD_PYTHON OR NOT ADSIM_IS_TOP_LEVEL)
    set(_adsim_install_default OFF)
else()
    set(_adsim_install_default ON)
endif()
option(ADSIM_INSTALL "Install the adsim library, its headers and the Adsim CMake package" ${_adsim_install_default})
# Release builds are portable by default: no -march=native, so a binary or wheel
# built on one machine runs on any CPU of the same architecture. Opt in with
# -DADSIM_NATIVE_ARCH=ON for a local build tuned to the build machine's CPU
# (GCC/Clang only).
option(ADSIM_NATIVE_ARCH "Compile Release builds with -march=native (non-portable binaries)" OFF)

# ============================================================
# Emscripten exceptions
# ============================================================
# Emscripten (a front end compiling the engine to WebAssembly with emcmake):
# the engine and cereal throw and catch C++ exceptions, which
# Emscripten disables by default. ADSIM_EMSCRIPTEN_EXCEPTIONS picks the scheme:
#   wasm (default)  native WebAssembly exceptions (-fwasm-exceptions): throw and
#                   catch stay inside the module. Measured on the web solver:
#                   10 % faster and a 12 % smaller .wasm than js. Needs
#                   Chrome 95, Firefox 100, Safari 15.2 or Node 17.
#   js              Emscripten's JavaScript-based exceptions (-fexceptions):
#                   every call that may throw goes through a JS trampoline.
# The whole tree (engine and dependencies) is compiled with the chosen flag, and
# Adsim::adsim passes it on to the consumer's compile and link, because every
# object in a module must use the same scheme. A consumer must not add
# -sNO_DISABLE_EXCEPTION_CATCHING=1 with wasm (Emscripten rejects the mix);
# threading and the other link flags are its choice.
set(ADSIM_EMSCRIPTEN_EXCEPTIONS "wasm" CACHE STRING
    "C++ exception scheme for Emscripten builds: wasm (native) or js (JavaScript-based)")
set_property(CACHE ADSIM_EMSCRIPTEN_EXCEPTIONS PROPERTY STRINGS wasm js)
set(ADSIM_EMSCRIPTEN_EXCEPTIONS_FLAG "")
if(EMSCRIPTEN)
    if(ADSIM_EMSCRIPTEN_EXCEPTIONS STREQUAL "wasm")
        set(ADSIM_EMSCRIPTEN_EXCEPTIONS_FLAG "-fwasm-exceptions")
    elseif(ADSIM_EMSCRIPTEN_EXCEPTIONS STREQUAL "js")
        set(ADSIM_EMSCRIPTEN_EXCEPTIONS_FLAG "-fexceptions")
    else()
        message(FATAL_ERROR "ADSIM_EMSCRIPTEN_EXCEPTIONS must be wasm or js, not '${ADSIM_EMSCRIPTEN_EXCEPTIONS}'")
    endif()
    add_compile_options(${ADSIM_EMSCRIPTEN_EXCEPTIONS_FLAG})
    message(STATUS "Adsim: Emscripten C++ exceptions: ${ADSIM_EMSCRIPTEN_EXCEPTIONS} (${ADSIM_EMSCRIPTEN_EXCEPTIONS_FLAG})")
endif()

# The Python extension is a shared module that links every static library
# below, so they all need position-independent code.
if(ADSIM_BUILD_PYTHON)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Older configures FORCE-wrote "-march=native" into the cached Release flags, so
# an existing build directory would silently keep it. Strip it once (unless the
# opt-in above is on, which adds it back below in a non-cached way).
foreach(_adsim_lang C CXX)
    if(CMAKE_${_adsim_lang}_FLAGS_RELEASE MATCHES "-march=native")
        string(REPLACE "-march=native" "" _adsim_flags "${CMAKE_${_adsim_lang}_FLAGS_RELEASE}")
        string(REGEX REPLACE "  +" " " _adsim_flags "${_adsim_flags}")
        string(STRIP "${_adsim_flags}" _adsim_flags)
        set(CMAKE_${_adsim_lang}_FLAGS_RELEASE "${_adsim_flags}" CACHE STRING
            "Flags used by the ${_adsim_lang} compiler during RELEASE builds." FORCE)
        message(STATUS "Removed -march=native from the cached CMAKE_${_adsim_lang}_FLAGS_RELEASE "
                       "(use -DADSIM_NATIVE_ARCH=ON to opt in)")
    endif()
endforeach()

if(ADSIM_NATIVE_ARCH)
    if(MSVC OR EMSCRIPTEN)
        message(WARNING "ADSIM_NATIVE_ARCH is ignored for MSVC and Emscripten builds")
    else()
        add_compile_options($<$<CONFIG:Release>:-march=native>)
    endif()
endif()

# ============================================================
# Third-party dependencies (git submodules in external/, plus the vendored
# external/eigen-3.4.0)
# ============================================================
foreach(_adsim_dep sundials cereal muparser dlib units)
    if(NOT EXISTS "${ADSIM_SOURCE_DIR}/external/${_adsim_dep}/CMakeLists.txt"
       AND NOT EXISTS "${ADSIM_SOURCE_DIR}/external/${_adsim_dep}/include")
        message(FATAL_ERROR
            "external/${_adsim_dep} is empty. Fetch the submodules first:\n"
            "  git submodule update --init --recursive")
    endif()
endforeach()

# The dependency settings below are cache entries without FORCE, so a parent
# project that already set one (e.g. BUILD_SHARED_LIBS) keeps its value.

# -- SUNDIALS (IDA/IDAS DAE solvers) --
set(SUNDIALS_INDEX_SIZE "64" CACHE STRING "SUNDIALS index type size (32 or 64)")
set(BUILD_STATIC_LIBS ON CACHE BOOL "Build static libraries")
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries")
set(EXAMPLES_ENABLE OFF CACHE BOOL "Build SUNDIALS examples")
set(EXAMPLES_INSTALL OFF CACHE BOOL "Install SUNDIALS examples")
set(EXAMPLES_ENABLE_C OFF CACHE BOOL "Build SUNDIALS C examples")
set(EXAMPLES_ENABLE_CXX OFF CACHE BOOL "Build SUNDIALS C++ examples")
set(BUILD_IDA ON CACHE BOOL "Build IDA DAE solver")
set(BUILD_IDAS ON CACHE BOOL "Build IDAS DAE solver with sensitivity")
set(BUILD_CVODE OFF CACHE BOOL "Build CVODE")
set(BUILD_CVODES OFF CACHE BOOL "Build CVODES")
set(BUILD_ARKODE OFF CACHE BOOL "Build ARKODE")
set(BUILD_KINSOL OFF CACHE BOOL "Build KINSOL")
add_subdirectory(external/sundials)

# -- muparser (expression evaluator) --
set(ENABLE_SAMPLES OFF CACHE BOOL "Don't build muparser samples")
set(ENABLE_OPENMP  OFF CACHE BOOL "Disable muparser OpenMP")
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build muparser as static")
add_subdirectory(external/muparser)

# -- dlib (global optimisation, used by adsim/solver/optimisation.h and the
# isotherm fitter) --
# Only find_min_global is used, so build just the six dlib sources it needs
# (see cmake/dlib_minimal.cmake) rather than all of dlib.
include(${ADSIM_SOURCE_DIR}/cmake/dlib_minimal.cmake)

# ============================================================
# Core solver library
# ============================================================
add_subdirectory(core)

# ============================================================
# Python extension module (optional; see python/ and pyproject.toml)
# ============================================================
if(ADSIM_BUILD_PYTHON)
    add_subdirectory(python)
endif()

# ============================================================
# Tests and benchmark (optional)
# ============================================================
if(ADSIM_BUILD_TESTS)
    enable_testing()
    set(gtest_force_shared_crt ON CACHE BOOL "")
    set(INSTALL_GTEST OFF CACHE BOOL "")
    if(NOT TARGET gtest_main)
        add_subdirectory(external/googletest)
    endif()
    include(CTest)
    add_subdirectory(core/test)
endif()
if(ADSIM_BUILD_BENCH)
    add_subdirectory(core/bench)
endif()
