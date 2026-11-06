# Minimal dlib: only the sources behind dlib::find_min_global.
#
# The engine uses one dlib facility, find_min_global (adsim/solver/optimisation.h
# and the isotherm fitter). It is header-only apart from
#   global_optimization/global_function_search.cpp   (the MaxLIPO + trust-region search)
#   threads/thread_pool_extension.cpp                (find_max_global's thread_pool overload)
#   threads/threads_kernel_shared.cpp                (thread pool support)
#   threads/threads_kernel_1.cpp                     (Windows threads; empty elsewhere)
#   threads/threads_kernel_2.cpp                     (POSIX threads; empty elsewhere)
#   test_for_odr_violations.cpp                      (the symbols dlib's headers reference
#                                                     to catch mismatched build settings;
#                                                     MSVC keeps those references, GCC drops them)
# so this file builds a static `dlib` target from those six files instead of
# add_subdirectory(external/dlib/dlib), which compiles all ~50 dlib sources.
# The target matches what dlib's own in-project build gave adsim: the
# DLIB_NO_GUI_SUPPORT switch is public, the include root is external/dlib, and
# on UNIX (including Emscripten) the thread library is linked privately.
#
# If a parent project already has a `dlib` target, that one is used instead.
#
# To go back to the full dlib build, replace include(cmake/dlib_minimal.cmake)
# in the root CMakeLists.txt with add_subdirectory(external/dlib/dlib).

if(TARGET dlib)
    if(NOT TARGET dlib::dlib)
        add_library(dlib::dlib ALIAS dlib)
    endif()
    return()
endif()

set(_adsim_dlib_dir ${ADSIM_SOURCE_DIR}/external/dlib/dlib)
add_library(dlib STATIC
    ${_adsim_dlib_dir}/global_optimization/global_function_search.cpp
    ${_adsim_dlib_dir}/threads/thread_pool_extension.cpp
    ${_adsim_dlib_dir}/threads/threads_kernel_shared.cpp
    ${_adsim_dlib_dir}/threads/threads_kernel_1.cpp
    ${_adsim_dlib_dir}/threads/threads_kernel_2.cpp
    ${_adsim_dlib_dir}/test_for_odr_violations.cpp)
add_library(dlib::dlib ALIAS dlib)

target_include_directories(dlib
    INTERFACE $<BUILD_INTERFACE:${ADSIM_SOURCE_DIR}/external/dlib>
              $<INSTALL_INTERFACE:include>
    PRIVATE ${ADSIM_SOURCE_DIR}/external/dlib)
# dlib's in-project build passes its configuration switches to every consumer
# (there is no generated config.h); DLIB_NO_GUI_SUPPORT is the only one set.
target_compile_options(dlib PUBLIC -DDLIB_NO_GUI_SUPPORT)
target_compile_features(dlib PUBLIC cxx_std_14)
if(CMAKE_COMPILER_IS_GNUCXX)
    target_compile_options(dlib PUBLIC $<$<COMPILE_LANGUAGE:CXX>:-Wreturn-type>)
endif()
if(MSVC)
    target_compile_options(dlib PUBLIC $<$<COMPILE_LANGUAGE:CXX>:/bigobj>)
endif()
if(UNIX)
    set(CMAKE_THREAD_PREFER_PTHREAD ON)
    find_package(Threads REQUIRED)
    target_link_libraries(dlib PRIVATE ${CMAKE_THREAD_LIBS_INIT})
endif()

unset(_adsim_dlib_dir)
