### Changed — energy KPI documentation
- The performance metrics page lists the three energy figures (electrical
  energy, thermal energy, cooling duty) in one table and shows how to convert
  J/kg to MJ/kg, GJ/t and kWh/t, and says that the thermal energy and cooling
  duty count the energy sources only: a heating medium's duty is not logged
  yet. The TSA tutorial reads the heater's thermal
  energy (6.07 MJ per kg of CO₂) and the cooling duty from the KPIs instead
  of working them out by hand, and changes the cooler power through the
  parameter registry.
- A closed-bed test checks that the logged heat duty of a step equals the
  analytic power × time, for heating and for cooling. Results do not change.
