# Guard: every file under the bundled data directories (examples/decks/,
# data/adsorbents/) must have a shell-safe name.
#
# Front ends that consume this repository copy these files as build assets. The
# ImGui desktop app, for example, uses hello_imgui, which copies each asset with
# its own POST_BUILD `cmake -E copy_if_different` command without VERBATIM, so a
# file named "AP3-60 Activated Carbon (Schell 2012).json" breaks its link step
# on Linux/macOS with `/bin/sh: 1: Syntax error: "(" unexpected`. The Python
# wheel also ships the files as package data.
#
# Allowed in every path component: letters, digits and . _ - , = + @
# (no spaces, no ( ) [ ] & ; ' " $ ` ! * ? < > | \ # ~ { } or non-ASCII).
# Name new files in lowercase snake_case, like the adsorbent library keys
# (data/adsorbents/<key>.json).
#
# Run as a test (ctest: DataFileNamesAreShellSafe):
#   cmake "-DASSET_DIRS=<dir1>;<dir2>" -P cmake/check_data_names.cmake
# or include() it and call adsim_check_data_names(<dir>).

# Grandfathered names that are shell-safe enough in practice (unquoted [...] is
# only a glob, and matches nothing here) and are referenced by published links,
# tests and validation records, so renaming them would break those references.
# Entries are relative to the checked directory.
# Do NOT add to this list; give new files a safe name instead.
set(ADSIM_DATA_NAME_ALLOWLIST
    "Ward_2022_Q=50[ccm].json"
    "Ward_2022_Q=100[ccm].json"
    "Ward_2022_Q=200[ccm].json"
)

function(adsim_check_data_names data_dir)
    file(GLOB_RECURSE _files LIST_DIRECTORIES false RELATIVE "${data_dir}" "${data_dir}/*")
    set(_bad "")
    foreach(_rel IN LISTS _files)
        list(FIND ADSIM_DATA_NAME_ALLOWLIST "${_rel}" _allowed)
        if(NOT _allowed EQUAL -1)
            continue()
        endif()
        # A ';' in a file name splits the glob result into pieces that do not
        # exist on disk, so a missing piece is also an unsafe name.
        if(NOT _rel MATCHES "^[A-Za-z0-9._,=+@/-]+$" OR NOT EXISTS "${data_dir}/${_rel}")
            list(APPEND _bad "${_rel}")
        endif()
    endforeach()
    if(_bad)
        list(JOIN _bad "\n    " _bad_text)
        message(FATAL_ERROR
            "Unsafe data file name(s) under ${data_dir}:\n    ${_bad_text}\n"
            "File names may only contain letters, digits and . _ - , = + @ "
            "(no spaces or shell characters such as parentheses). Rename them in "
            "lowercase snake_case (e.g. ap3_60_schell2012.json) and update their "
            "references; see cmake/check_data_names.cmake.")
    endif()
endfunction()

if(CMAKE_SCRIPT_MODE_FILE AND DEFINED ASSET_DIRS)
    foreach(_dir IN LISTS ASSET_DIRS)
        if(NOT IS_DIRECTORY "${_dir}")
            message(FATAL_ERROR "Data directory not found: ${_dir}")
        endif()
        adsim_check_data_names("${_dir}")
        message(STATUS "All file names under ${_dir} are shell-safe.")
    endforeach()
endif()
