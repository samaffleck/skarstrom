### Changed — faster builds: only the dlib sources the engine uses
- The build compiles 5 of dlib's ~50 sources: the ones behind
  `dlib::find_min_global` (the global optimiser and its thread pool), not the
  whole library. A cold build of dlib takes about a third of the time it did.
  The `dlib::dlib` target, the installed `Adsim::dlib` export and the
  installed dlib headers are unchanged (`cmake/dlib_minimal.cmake`).
