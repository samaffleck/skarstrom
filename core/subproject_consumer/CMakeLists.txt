# Minimal front-end project that consumes this repository as a SUBPROJECT:
#
#   add_subdirectory(<path to the skarstrom checkout> skarstrom)
#   target_link_libraries(app PRIVATE Adsim::adsim)
#
# This is how a front end (desktop GUI, WebAssembly module, ...) that pins the
# engine as a git submodule builds it. Used by scripts/test_subproject_consumer.sh
# natively and with emcmake (Emscripten). Not part of the main build.
#
#   cmake -S core/subproject_consumer -B <dir> [-DSKARSTROM_DIR=<checkout>]
#   cmake --build <dir>
#   <dir>/frontend examples/decks/<deck>.json      (node <dir>/frontend.js ... under Emscripten)
cmake_minimum_required(VERSION 3.20)
project(SkarstromFrontendExample LANGUAGES C CXX)

# The parent's own settings must survive add_subdirectory().
set(CMAKE_CXX_STANDARD 20)
set(FRONTEND_MARKER "parent" CACHE STRING "")

get_filename_component(_default_dir "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(SKARSTROM_DIR "${_default_dir}" CACHE PATH "Skarstrom checkout")
add_subdirectory(${SKARSTROM_DIR} skarstrom)

foreach(_v ADSIM_INSTALL ADSIM_BUILD_TESTS ADSIM_BUILD_BENCH)
    if(${_v})
        message(FATAL_ERROR "${_v} should default to OFF in a subproject build")
    endif()
endforeach()
if(NOT FRONTEND_MARKER STREQUAL "parent")
    message(FATAL_ERROR "the engine clobbered a parent cache variable")
endif()

add_executable(frontend main.cpp)
target_link_libraries(frontend PRIVATE Adsim::adsim)

if(EMSCRIPTEN)
    # The C++ exception flag comes with Adsim::adsim (ADSIM_EMSCRIPTEN_EXCEPTIONS);
    # NODERAWFS lets the test read the deck from the host file system under Node.js.
    target_link_options(frontend PRIVATE -sALLOW_MEMORY_GROWTH=1 -sNODERAWFS=1 -sENVIRONMENT=node)
endif()
