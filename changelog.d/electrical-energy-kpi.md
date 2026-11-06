### Added — electrical energy KPI
- The specific energy KPI is now also reported under its explicit name,
  `electrical_energy` (compressor and vacuum-pump work over the cycle per kg
  of net product, J/kg): a `KPIRecord` field and Python attribute, a key of
  `KPIRecord.as_dict()` / `kpi_table()` / the CLI's `--json` KPIs, a
  sensitivity output (`Component.electrical_energy`), an appended
  `electrical energy [J/kg]` column of the CLI's `kpis.csv`, and the
  `ELECTRICAL_ENERGY` enumerator of `KPIType` and `KPIVariableType` (the
  same enumerator as `ENERGY_CONSUMPTION`). The KPI's display label is now
  "Electrical Energy", and `skarstrom run` prints it as
  `electrical energy [J/kg]`.

### Changed — `energy_consumption` is documented as the electrical energy
- `energy_consumption` keeps its name, value and meaning in KPI records,
  sensitivity outputs (still the default), optimiser objectives, the cost
  estimate and the existing `kpis.csv` columns; it is documented as the
  specific electrical energy. Step energy sources (heaters, coolers) are
  still not included. Results do not change.
