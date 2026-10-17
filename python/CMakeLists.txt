# skarstrom._core: the pybind11 extension module behind the Python package.
# Built through scikit-build-core (pyproject.toml at the repository root), which
# runs the root CMakeLists.txt with ADSIM_BUILD_PYTHON=ON and installs only the
# "python" component into the wheel: the compiled module, the licence texts and
# the bundled data (example decks and the adsorbent library). No C++ sources or
# headers are installed.

find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_core MODULE src/bindings.cpp)
target_link_libraries(_core PRIVATE adsim)
target_compile_definitions(_core PRIVATE EIGEN_MPL2_ONLY)

install(TARGETS _core
    LIBRARY DESTINATION skarstrom COMPONENT python
    RUNTIME DESTINATION skarstrom COMPONENT python)

# Licence texts shipped inside the package (the wheel is a binary
# distribution, so the project licence and the third-party notices apply).
install(FILES
    ${ADSIM_SOURCE_DIR}/LICENSE
    ${ADSIM_SOURCE_DIR}/NOTICE
    ${ADSIM_SOURCE_DIR}/THIRD_PARTY_NOTICES.md
    DESTINATION skarstrom/licenses COMPONENT python)

# Package data: every example deck (skarstrom.examples, `skarstrom examples`)
# and the curated adsorbent library (skarstrom.adsorbent_library()).
install(DIRECTORY ${ADSIM_SOURCE_DIR}/examples/decks/
    DESTINATION skarstrom/data/decks COMPONENT python
    FILES_MATCHING PATTERN "*.json")
install(DIRECTORY ${ADSIM_SOURCE_DIR}/data/adsorbents/
    DESTINATION skarstrom/data/adsorbents COMPONENT python
    FILES_MATCHING PATTERN "*.json")
