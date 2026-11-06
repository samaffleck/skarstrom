### Changed — hint for varying a blend fraction
- The error for layer mass fractions that do not sum to 1 now also says how to
  vary one solid's fraction in a study: change the others with it (for two
  solids, x and 1 - x together at each point). A mass fraction is a number, not
  an expression. `docs/model-conventions.md` gives a Python loop and the CLI
  form, and says that a sensitivity study cannot vary a blend fraction yet.
