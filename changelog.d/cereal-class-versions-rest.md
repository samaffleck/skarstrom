### Fixed — class versions of the fluid, connections and CSS acceleration
- `adsim::Fluid`, `adsim::Connection` and `adsim::AccelerateCyclicSteadyState`
  now declare their `CEREAL_CLASS_VERSION` in their headers too, so every
  serialised engine class does. Saved decks do not change; the
  `CerealClassVersions.*` test now pins these three classes.
