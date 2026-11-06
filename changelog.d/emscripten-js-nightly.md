### Changed — shorter pull-request CI: the JS-exceptions Emscripten build runs nightly
- The "Engine with Emscripten" CI job builds the engine with native
  WebAssembly exceptions (the default) on every pull request. The second build,
  with JavaScript-based exceptions (`ADSIM_EMSCRIPTEN_EXCEPTIONS=js`), now runs
  only on pushes to main, the nightly run and manual runs.
