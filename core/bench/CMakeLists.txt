# Solver benchmark runner (see bench_solver.cpp). Not registered with CTest.
add_executable(bench_solver bench_solver.cpp)
target_link_libraries(bench_solver PRIVATE adsim)
