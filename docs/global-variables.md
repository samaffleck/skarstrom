# Global Variables, Expressions and Units

Almost every number in a deck (step times, pressures, flows, temperatures, inlet and initial mole fractions, isotherm and kinetic parameters, layer and wall properties, solver tolerances) is stored as an **expression**, not just a value. You can type a plain number, a number with a unit such as `3.5 [bar]`, or a formula that refers to your own named **global variables**, such as `P_high / 2`. This page describes what an expression may contain, how units are converted and when expressions are evaluated.

---

## Global variables in a deck

A deck holds a list of global (user) variables. Each one has:

| Field | Meaning |
|---|---|
| **Name** | The identifier you use in expressions, e.g. `P_high`. Names are case-sensitive and must not clash with the built-in names listed below. |
| **Type** | The physical quantity (Pressure, Temperature, Time, Velocity, Molar Flow Rate, …, or Dimensionless / Custom). It fixes the SI base unit the value is stored in, shown next to the value. |
| **Expression** | A number or a formula, optionally ending in a unit, e.g. `3.5 [bar]` or `2*P_low`. |
| **Value** | The evaluated value in SI base units. Evaluation fails (and the deck is not runnable) if the expression cannot be evaluated. |
| **Description** | Free text for your own notes. |

Variables are saved with the deck. From Python a variable is the parameter `variable[<name>]`, and any other numeric parameter can be set to a formula that uses it:

```python
sim["variable[t_feed]"] = "40 [s]"
sim["step[Feed].time"] = "t_feed"          # both steps now follow t_feed
sim["step[Rinse].time"] = "t_feed"
sim["endpoint[3].pressure"] = "P_high / 2"
sim["endpoint[1].y[CarbonDioxide]"] = "y_feed"       # an inlet composition
sim["endpoint[1].y[Nitrogen]"] = "1 - y_feed"
```

Use a variable's name in any numeric field, for example `P_high` in a pressure profile's final value, or `t_feed` in two step times so that they always stay equal.

## What an expression can contain

Expressions are evaluated by [muparser](https://beltoforion.de/en/muparser/):

- numbers, including scientific notation (`1.5e-3`);
- the operators `+ - * / ^`, parentheses and comparisons; the ternary `cond ? a : b` is available;
- functions: `sin cos tan cot asin acos atan atan2 sinh cosh tanh coth asinh acosh atanh exp log ln log2 log10 sqrt abs sign rint min max sum avg` (`log` and `ln` are both the natural logarithm; use `^` for powers, there is no `pow` function);
- the constants **`R`** = 8.31446261815324 J mol⁻¹ K⁻¹ and **`PI`** = 3.14159265358979;
- the state variables **`T`**, **`P`**, **`c`**, **`y`** and **`v`** (see [State variables](#state-variables-t-p-c-y-v) for when they mean anything);
- the names of your global variables.

A global variable cannot be called `T`, `P`, `c`, `y`, `v`, `R` or `PI`, nor `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`, `sinh`, `cosh`, `tanh`, `exp`, `log`, `log2`, `log10`, `sqrt`, `abs`, `min`, `max`, `pow` or `sign`. The engine does not reject all of them, so it is up to you: a variable called `T` would replace the state variable in every expression. Avoid the other function names too.

## Units

### Inline units

End an expression with a unit in square brackets to enter it in that unit:

```text
3.5 [bar]          ->  350000 Pa
25 [degC]          ->  298.15 K
100 [ccm]          ->  1.6667e-6 m³/s
15 [%]             ->  0.15
(2 + 3) [bar]      ->  500000 Pa
```

The rules, from `core/src/math/unit_parser.cpp`:

- The bracket must be the **last thing** in the expression. `2 [bar] + 1` is an error.
- The unit applies to the **whole** result: the formula is evaluated first and the result is converted from that unit to SI (`base = value × factor + offset`, the offset being for °C and °F).
- The symbol must be one the parser knows; an unknown symbol is an error. Accepted symbols include, by quantity:
  - pressure `Pa kPa MPa bar atm psi Torr`; temperature `K degC °C C degF °F F`; length `m cm mm um µm nm A in ft`; time `s ms us µs min h hr d`; mass `kg g mg t tonne lb`;
  - energy `J kJ MJ cal kcal Btu` and per mole `J/mol kJ/mol …`; power `W kW MW GW mW`; amount `mol mmol kmol`; volume `m^3 m³ L l mL ml ft^3`;
  - volumetric flow `m^3/s m^3/h L/s L/min ccm mL/min SCCM ft^3/min`; mass flow `kg/s g/s kg/h t/h …`; molar flow `mol/s mol/min mol/h mmol/s kmol/h Nm^3/h Nm^3/s SCFM …` (normal and standard volumes use 0 °C and 1 atm, 22.414 L/mol);
  - velocity `m/s cm/s mm/s m/min km/h ft/s mph`; inverse time `1/s 1/min 1/h`; inverse pressure `1/Pa 1/kPa 1/bar 1/atm …`;
  - density `kg/m^3 g/L g/cm^3 …`; concentration `mol/m^3 mol/L mmol/L …`; loading `mol/kg mmol/g mol/g …`;
  - heat capacity `J/(kg·K) kJ/(kg·K) cal/(g·K)`; conductivity `W/(m·K) mW/(m·K)`; molar mass `kg/mol g/mol kg/kmol`; specific energy `J/kg kJ/kg MJ/kg kWh/kg kWh/tonne`; diffusivity `m^2/s cm^2/s mm^2/s`; area `m^2 cm^2 mm^2`; pressure gradient `Pa/m kPa/m bar/m …`;
  - fractions `- mol/mol % ppm`.

> **Variables are already in SI.** A variable's value is stored in its SI base unit, so a unit after a formula that uses variables rescales them again: with `a = 5 [bar]` (500 000 Pa), `a/2 [kPa]` gives 250 000 × 1000 = 2.5 × 10⁸ Pa. Write `a/2` instead.

### Plain numbers are SI

A plain number with no unit is taken to be in the SI base unit of the field. Front ends that offer a unit selector store a value typed in another unit as `number [unit]`, which the engine converts on evaluation. A formula is stored exactly as written, so give it an inline unit or write it in SI.

### Type checking

If a variable's **Type** is a specific quantity and its expression ends in a unit of a different quantity (say Type *Pressure* and `300 [K]`), the deck is marked not runnable with the message *inline unit is Temperature but variable category is Pressure*. A variable of type *Dimensionless* (the default for a new row) or *Custom* accepts any unit and keeps its type.

## When expressions are evaluated

The engine evaluates every expression in the deck before each run (`UserData::EvaluateAllExpressions`). The Python package also evaluates a parameter as soon as you set it, so `sim.get(path)` returns the new value straight away.

### Global variables first, in list order

Global variables are evaluated **from the top of the list down, in one pass**, before anything else:

- A variable may use any variable listed **above** it.
- A variable that uses one listed **below** it gets that variable's value from the previous evaluation, so it lags one edit behind. For example, with `A = 2*B` listed above a new variable `B = 3 [bar]`, `A` evaluates to 0 the first time and only becomes 600 000 Pa after the next edit.
- Circular references (`A = B + 1`, `B = A + 1`) are not detected; they produce values that change with every evaluation.
- A reference to an undefined name is an error and the deck is not runnable.

So **keep the list in dependency order**: define base quantities first and derived ones after them. This is the engine's current behaviour (`ResolveVariableExpressions` in `core/src/objects/user_data.cpp`).

### Then every other field

After the variables, all other expressions are evaluated once, before the run starts. Their values are fixed for the whole run.

### State variables (T, P, c, y, v)

The names `T`, `P`, `c`, `y` and `v` stand for the local temperature [K], pressure [Pa], total concentration [mol/m³], component mole fraction [–] and velocity [m/s]. They are only re-evaluated with the local state in **two places**:

| Field | Re-evaluated | Values supplied |
|---|---|---|
| **Constant** mass-transfer coefficient (Adsorbent → rate model) | in every cell, at every residual evaluation | `T` = gas temperature of the cell, `P` = cell pressure, `y` = mole fraction of the component the coefficient belongs to; `c` = 0 and `v` = 0 |
| **Arrhenius** mass-transfer coefficient: pre-exponential factor and activation energy | same | same |

So `0.05*(P/1e5)^0.5` is a valid pressure-dependent constant-LDF coefficient, but `c` and `v` evaluate to zero there. (`MassTransferCoefficientConstant` and `MassTransferCoefficientArrhenius` in `core/src/equations/mass_transfer_model.cpp`.)

**Everywhere else**, including inside global variables, a state variable is evaluated **once** with placeholder values, `T` = 298.15 K, `P` = 10⁵ Pa, `c` = 40 mol/m³, `y` = 0 and `v` = 0, and the result is used as a constant for the whole run. For example a variable `A = T*2` is 596.3, whatever the bed temperature. Use the dedicated models (temperature-dependent isotherms, the Arrhenius coefficient, transport correlations, inlet temperature profiles) for real state dependence.

## Global variables in studies

The parameter registry (`sim.parameters()`, `skarstrom params`) lists global variables as parameters (`variable[<name>]`), and so does the engine's optimiser. Setting a parameter to a number overwrites its expression, so varying a variable moves every field that refers to it, while setting a field that refers to a variable to a number cuts that link. To change several inputs together in a [parametric study or a fit](parameterisation.md), make them all refer to one variable and vary the variable. A field written as a formula is not *tunable* (`tunable_only=True` leaves it out), because setting a number would replace the formula.
