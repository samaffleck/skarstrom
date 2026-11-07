# Adsorbent library

- `library.json` is the **source** catalogue, curated by hand: materials, sources,
  and one record per entry with its citation and per-parameter provenance.
- `<key>.json` is the serialised adsorbent (physics) of the entry `<key>`.
- `library.merged.json` is the **artefact**: the validated, merged catalogue
  (every `library.json` record plus its `recommended` flag and its adsorbent
  inlined as `adsorbent`), in the order of `docs/adsorbent-library.md`.
  Tools that embed the library read this file; they need no validator of their own.

`library.merged.json` and `docs/adsorbent-library.md` are generated, never edited:

```bash
python scripts/docs/gen_adsorbent_library.py          # validate, then write both
python scripts/docs/gen_adsorbent_library.py --check  # validate; fail if either is out of date (CI)
```

The rules for adding an entry (every parameter cited, one recommended entry per
material) are in `AGENTS.md` and `CONTRIBUTING.md`.
