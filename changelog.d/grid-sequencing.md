### Added — grid sequencing benchmark
- A slow-tier test (`GridSequencing`) compares three ways of reaching a
  fine-grid cyclic steady state on six example decks: the fine grid from
  cold, a coarse grid to CSS followed by a warm-started fine grid, and three
  levels (N/4, N/2, N cells). It asserts that every path reaches the same CSS
  within the CSS reference bound and that the warm-started fine level needs
  fewer cycles than the cold run; timings are printed, not asserted.
- `docs/cyclic-steady-state.md` tabulates cycles and wall time per path and
  recommends when sequencing pays: 1.4-2.1x on the slower VSA decks, nothing
  on a fast TSA deck or the PSA deck; a grid-refinement study should always
  warm-start each level from the previous level's CSS.
