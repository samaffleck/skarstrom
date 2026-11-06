### Changed — faster builds: only the dlib sources the engine uses
- The build compiles 6 of dlib's ~50 sources: the ones behind
  `dlib::find_min_global` (the global optimiser, its thread pool and dlib's
  build-consistency check), not the
  whole library. A cold dlib build needs about a seventh of the CPU time it did.
  The `dlib::dlib` target, the installed `Adsim::dlib` export and the
  installed dlib headers are unchanged (`cmake/dlib_minimal.cmake`).
