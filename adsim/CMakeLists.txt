# adsim - Adsorption Simulator Library
# Dependencies (sundials, CoolProp, cereal, Eigen) are managed by the root CMakeLists.txt

include(GNUInstallDirs)

# Add the static library
add_library(adsim STATIC)
add_library(Adsim::adsim ALIAS adsim)

set(EIGEN_DIR ${ADSIM_SOURCE_DIR}/external/eigen-3.4.0)

# Third-party include paths are build-tree only; the installed package ships
# the headers it needs under include/adsim/third_party (see install rules below).
target_include_directories(adsim PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${EIGEN_DIR}>
    $<BUILD_INTERFACE:${EIGEN_DIR}/Eigen>
    $<BUILD_INTERFACE:${ADSIM_SOURCE_DIR}/external/muparser/include>
)

# Generated version header (adsim/version.h) — the version comes from the root
# project(AdSim VERSION ...) call. Written into the build tree, never committed.
set(ADSIM_GENERATED_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated/include)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/adsim/version.h.in
    ${ADSIM_GENERATED_INCLUDE_DIR}/adsim/version.h
    @ONLY
)
target_include_directories(adsim PUBLIC
    $<BUILD_INTERFACE:${ADSIM_GENERATED_INCLUDE_DIR}>
)
# Installed layout: <prefix>/include/adsim/... (public + generated headers) and
# <prefix>/include/adsim/third_party/{Eigen,cereal,CoolProp}.
target_include_directories(adsim PUBLIC
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/adsim/third_party>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/adsim/third_party/CoolProp>
)

# Cereal (header-only)
target_include_directories(adsim PUBLIC
    $<BUILD_INTERFACE:${ADSIM_SOURCE_DIR}/external/cereal/include>
)

# Guarantee only MPL2-licensed Eigen code is compiled in (excludes the
# LGPL-licensed Eigen files), so binaries carry no LGPL obligations.
target_compile_definitions(adsim PUBLIC EIGEN_MPL2_ONLY)

# The public headers include CoolProp.h, which pulls in CoolProp's bundled fmt;
# every consumer needs the same consteval workaround the engine is built with
# (see the root CMakeLists.txt). The installed package gets it through
# Adsim::CoolProp (cmake/AdsimConfig.cmake.in).
target_compile_definitions(adsim PUBLIC $<BUILD_INTERFACE:FMT_CONSTEVAL=>)

# Add all .cpp files in the src directory
file(GLOB_RECURSE ADS_SRC CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

# Add source files
target_sources(adsim PRIVATE ${ADS_SRC})

# Enable precompiled headers
target_precompile_headers(adsim PRIVATE src/pch.h)

# Link dependencies. CoolProp has no CMake export of its own (and its include
# paths point into the source tree), so the installed package re-creates it as
# the imported target Adsim::CoolProp (see cmake/AdsimConfig.cmake.in).
target_link_libraries(adsim PUBLIC
    $<BUILD_INTERFACE:CoolProp>
    $<INSTALL_INTERFACE:Adsim::CoolProp>
    sundials_ida sundials_idas dlib::dlib muparser)

# Set C++ standard
target_compile_features(adsim PUBLIC cxx_std_20)

# MSVC only, C++ sources only (never the resource compiler):
# - /utf-8: the sources are UTF-8 and contain UTF-8 string literals (unit
#   symbols such as "m²" and "µs"). Without it MSVC reads files in the system
#   code page.
# - /bigobj: cereal + Eigen translation units can pass the default object-file
#   section limit (error C1128).
# Passed on to the in-tree consumers (tests, bench, Python module), not to
# installed packages.
target_compile_options(adsim PUBLIC
    "$<BUILD_INTERFACE:$<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/utf-8>>"
    "$<BUILD_INTERFACE:$<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/bigobj>>")

# MSVC debugger visualization
if(MSVC)
    set_property(TARGET adsim PROPERTY
                VS_DEBUGGER_VISUALIZERS
                ${CMAKE_CURRENT_SOURCE_DIR}/Eigen.natvis)
endif()

# ============================================================
# Install + CMake package (find_package(Adsim)) — opt-in via ADSIM_INSTALL
# ============================================================
# Installed tree (relative to CMAKE_INSTALL_PREFIX):
#   lib/libadsim.a, lib/libdlib.a, lib/libCoolProp.a      (static libraries)
#   lib/libsundials_*.a, lib/libmuparser.a                (installed by SUNDIALS /
#   lib/cmake/sundials, lib/cmake/muparser                 muparser's own rules)
#   include/adsim/...                                      (public API + version.h)
#   include/adsim/third_party/{Eigen,cereal,CoolProp,CoolProp/fmt} (headers the API includes)
#   include/dlib, include/sundials, include/ida, ...       (dlib / SUNDIALS headers)
#   lib/cmake/Adsim/AdsimConfig.cmake (+ Version, Targets)
# Consumers: find_package(Adsim CONFIG REQUIRED) and link Adsim::adsim.
if(ADSIM_INSTALL)
    include(CMakePackageConfigHelpers)
    set(ADSIM_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/Adsim)

    # dlib is added in "project build" mode, which skips dlib's own install
    # rules, so adsim's export set carries it (exported as Adsim::dlib, with
    # the same compile options the library was built with).
    install(TARGETS adsim dlib
        EXPORT AdsimTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

    # Public headers + the generated version header.
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/adsim
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
    install(FILES ${ADSIM_GENERATED_INCLUDE_DIR}/adsim/version.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/adsim)

    # Header-only / private third-party headers reachable from the public API.
    set(_adsim_tp ${CMAKE_INSTALL_INCLUDEDIR}/adsim/third_party)
    install(DIRECTORY ${EIGEN_DIR}/Eigen DESTINATION ${_adsim_tp})
    install(DIRECTORY ${ADSIM_SOURCE_DIR}/external/cereal/include/cereal DESTINATION ${_adsim_tp})
    # Only the CoolProp.h include closure is part of adsim's API (objects/fluid.h).
    set(_adsim_coolprop_inc ${ADSIM_SOURCE_DIR}/external/CoolProp/include)
    install(FILES
        ${_adsim_coolprop_inc}/CoolProp.h
        ${_adsim_coolprop_inc}/DataStructures.h
        ${_adsim_coolprop_inc}/CPnumerics.h
        ${_adsim_coolprop_inc}/CPstrings.h
        ${_adsim_coolprop_inc}/Exceptions.h
        ${_adsim_coolprop_inc}/PlatformDetermination.h
        DESTINATION ${_adsim_tp}/CoolProp)
    # CPstrings.h pulls in CoolProp's bundled fmt (header-only mode), which is
    # also why consumers get FMT_CONSTEVAL= (see AdsimConfig.cmake.in).
    install(DIRECTORY ${ADSIM_SOURCE_DIR}/external/CoolProp/externals/fmtlib/include/fmt
        DESTINATION ${_adsim_tp}/CoolProp)
    install(FILES $<TARGET_FILE:CoolProp> DESTINATION ${CMAKE_INSTALL_LIBDIR})

    # dlib headers (dlib's own install rule is skipped in project-build mode).
    install(DIRECTORY ${ADSIM_SOURCE_DIR}/external/dlib/dlib/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dlib
        FILES_MATCHING PATTERN "*.h"
        PATTERN "test" EXCLUDE
        PATTERN "external" EXCLUDE
        PATTERN "cmake_utils" EXCLUDE)

    # Licence texts for everything the package redistributes.
    install(FILES
        ${ADSIM_SOURCE_DIR}/LICENSE
        ${ADSIM_SOURCE_DIR}/NOTICE
        ${ADSIM_SOURCE_DIR}/THIRD_PARTY_NOTICES.md
        DESTINATION ${CMAKE_INSTALL_DATADIR}/doc/Adsim)
    install(FILES ${EIGEN_DIR}/COPYING.MPL2 DESTINATION ${CMAKE_INSTALL_DATADIR}/doc/Adsim/licenses/Eigen)
    install(FILES ${ADSIM_SOURCE_DIR}/external/cereal/LICENSE DESTINATION ${CMAKE_INSTALL_DATADIR}/doc/Adsim/licenses/cereal)
    install(FILES ${ADSIM_SOURCE_DIR}/external/CoolProp/LICENSE DESTINATION ${CMAKE_INSTALL_DATADIR}/doc/Adsim/licenses/CoolProp)
    install(FILES ${ADSIM_SOURCE_DIR}/external/CoolProp/externals/fmtlib/LICENSE.rst DESTINATION ${CMAKE_INSTALL_DATADIR}/doc/Adsim/licenses/fmt)
    install(FILES ${ADSIM_SOURCE_DIR}/external/dlib/LICENSE.txt DESTINATION ${CMAKE_INSTALL_DATADIR}/doc/Adsim/licenses/dlib)

    install(EXPORT AdsimTargets
        FILE AdsimTargets.cmake
        NAMESPACE Adsim::
        DESTINATION ${ADSIM_INSTALL_CMAKEDIR})

    configure_package_config_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/AdsimConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/AdsimConfig.cmake
        INSTALL_DESTINATION ${ADSIM_INSTALL_CMAKEDIR}
        PATH_VARS CMAKE_INSTALL_LIBDIR CMAKE_INSTALL_INCLUDEDIR)
    # SemVer, pre-1.0: patch releases are drop-in; a new MINOR may change the
    # API (as a new MAJOR will after 1.0).
    write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/AdsimConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMinorVersion)
    install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/AdsimConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/AdsimConfigVersion.cmake
        DESTINATION ${ADSIM_INSTALL_CMAKEDIR})
endif()
