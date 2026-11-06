## Wall Energy Balance

The column wall is a separate thermal body with its own temperature $T_w(z,t)$. It stores heat, conducts it along the column, exchanges heat with the packed bed on the inside and with the surroundings on the outside. The wall model is set on **Column Properties → Wall**, and the heat transfer coefficients on each layer's **Internal HTC** and **External HTC** tabs. The implementation is `UpdateWallTemperature` in `core/src/equations/adsorbent_equations.cpp` and the wall data in `core/include/adsim/objects/wall.h`.

### Governing equation

Every term is written per unit volume of wall. For an annular wall of inner diameter $d_i$ and outer diameter $d_o$, a slice of length $\Delta z$ has volume $\tfrac{\pi}{4}(d_o^2 - d_i^2)\thinspace \Delta z$, inner area $\pi d_i \Delta z$ and outer area $\pi d_o \Delta z$, so with local thermal equilibrium in the bed:

$$
\rho_w C_{p,w}\frac{\partial T_w}{\partial t}
= k_w\frac{\partial^2 T_w}{\partial z^2}
{}+ \frac{4\thinspace d_i\thinspace h_i}{d_o^2 - d_i^2}\left(T - T_w\right)
{}- \frac{4\thinspace d_o\thinspace h_o}{d_o^2 - d_i^2}\left(T_w - T_{amb}\right)
$$

With local thermal non-equilibrium, the bed side of the exchange is split between the gas and the solid in proportion to the bed voidage:

$$
\frac{4\thinspace d_i\thinspace h_i}{d_o^2 - d_i^2}\Big[\varepsilon_b\left(T - T_w\right) + \left(1-\varepsilon_b\right)\left(T_s - T_w\right)\Big]
$$

The same exchange appears with the opposite sign in the bed's [energy balance](energy-balance.md) as $\tfrac{4 h_i}{d_i}(T - T_w)$ per unit bed volume. Both terms describe the same heat flow through the same area, so the heat the wall absorbs is exactly the heat the bed loses.

| Symbol           | Description                                             | Unit           | Where it is set |
|------------------|---------------------------------------------------------|----------------|-----------------|
| $T_w$            | Wall temperature                                        | K              | initial value: **Cycle Design → Initial Conditions → Wall temperature** |
| $T$, $T_s$       | Gas (and, in LTNE, solid) temperature of the bed        | K              | |
| $T_{amb}$        | Ambient temperature                                     | K              | **Column Properties → Geometry → Ambient temperature** |
| $\rho_w$, $C_{p,w}$, $k_w$ | Wall density, heat capacity, axial thermal conductivity | kg m⁻³, J kg⁻¹ K⁻¹, W m⁻¹ K⁻¹ | **Wall → Material** (a preset, or *Custom*) |
| $t_w$            | Wall thickness                                          | m              | **Wall → Thickness** (or calculated from a design pressure: [Wall Design](wall-design.md)) |
| $d_i$            | Inner diameter: the column diameter                     | m              | **Geometry** |
| $d_o$            | Outer diameter, $d_i + 2\thinspace t_w$                          | m              | |
| $h_i$            | Bed-to-wall (internal) heat transfer coefficient        | W m⁻² K⁻¹      | layer internal heat transfer model |
| $h_o$            | Wall-to-ambient (external) heat transfer coefficient    | W m⁻² K⁻¹      | layer external heat transfer model |
| $\varepsilon_b$  | Bed voidage                                             | –              | layer properties |

Notes:

- For a square or rectangular column, $d_i$ is the hydraulic diameter of the cross-section and the wall is treated as an annulus of thickness $t_w$ around it.
- $h_i$ and $h_o$ are constants or come from the [transport correlations](transport-correlations.md#internal-heat-transfer-coefficient) (Dixon, Leva, Holman), evaluated in every cell of the layer. In a [void region](getting-started.md) they are zero, so the wall there exchanges no heat with the gas or the surroundings and only conducts along the column.
- Setting $h_o = 0$ makes the column adiabatic towards the surroundings; setting $h_i = 0$ decouples the wall from the bed altogether.
- Axial conduction is discretised with central differences on the same cells as the bed.
- In an [isothermal run](energy-balance.md#isothermal-mode) the wall is held at the isothermal temperature and this equation is not solved.

### Boundary conditions at the column ends

The wall has two ends, at $z = 0$ and $z = L$. **Column Properties → Wall → Wall temperature boundary condition** selects what happens there (`wall_temperature_boundary_condition.cpp`):

- **Zero flux** (the default): no heat is conducted out of the wall ends, $\partial T_w/\partial z = 0$ at $z = 0$ and $z = L$.
- **Fixed**: the wall temperature at both ends is held at a set value $T_{w,bc}$ (implemented by setting the mean of the end cell and its ghost cell to $T_{w,bc}$). Use it when the column flanges are held at a known temperature by a large mass of pipework or a temperature-controlled jacket.

$$
\text{Zero flux: } \left.\frac{\partial T_w}{\partial z}\right|_{z=0,\thinspace L} = 0
\qquad\qquad
\text{Fixed: } \left.T_w\right|_{z=0,\thinspace L} = T_{w,bc}
$$

This condition only concerns the ends of the wall. Heat loss through the side of the column is governed by $h_o$ and $T_{amb}$, or, during a step that has one, by a [heating medium](#heating-medium-jacket).

### Heating medium (jacket)

A temperature swing or direct air capture column is often heated and cooled by a fluid in a jacket around it (steam, hot oil or cooling water) rather than left to the surroundings. Each step can have such a **heating medium** on the wall's outer face. While the step runs, the outer-face term of the wall balance becomes

$$
-\frac{4\thinspace d_o\thinspace h_m}{d_o^2 - d_i^2}\left(T_w - T_m(t)\right)
$$

so the medium's heat transfer coefficient $h_m$ replaces every layer's $h_o$, and the medium temperature $T_m$ replaces $T_{amb}$, over the whole length of the column (one jacket per bed; a [void region](getting-started.md) still exchanges no heat). Heat reaches the bed through the wall, so the wall's heat capacity and $h_i$ set how fast the bed follows the medium.

| Symbol  | Description | Unit | Parameter path |
|---------|-------------|------|----------------|
| –       | Whether the step has a medium (off by default) | – | `step[ID].medium.enabled` |
| $T_m$   | Medium temperature: a constant, or a ramp over the step | K | `step[ID].medium.temperature` (ramp: `.temperature_linear_start`, `_end`, `_time`, ...) |
| $h_m$   | Medium-to-wall heat transfer coefficient on the outer surface | W m⁻² K⁻¹ | `step[ID].medium.heat_transfer_coefficient` |

- $T_m$ takes the same profiles as an [inlet temperature](boundary-conditions.md#inlet-temperature-profile): constant, linear, exponential, log-linear or double exponential, with $t$ the time since the step started. The medium temperature is set, not solved: there is no energy balance on the medium, so the jacket behaves as if its flow were large enough to hold $T_m$.
- A step without a medium (the default, and every deck saved before media existed) keeps the ambient exchange. A medium is set step by step, so a heating step and a cooling step can use different fluids, and a step left without one returns to $h_o$ and $T_{amb}$.
- The logged *External Heat Transfer Coefficient* shows $h_m$ during a step with a medium.
- All three values are registry parameters, so they can be set from the command line (`--set "step[Heat].medium.enabled=true"`), swept, optimised and used in [sensitivity studies](sensitivity.md). In Python, `sim.steps["Heat"].medium.temperature = 393.15`.
- $h_m = 0$ makes the step adiabatic towards the outside, whatever the layers' $h_o$.

#### Medium duty

The solver logs, for every step, the heat the medium passes into the wall and the heat it takes back out, integrated over the step:

$$
Q_{m,\text{in}} = \int_{\text{step}} \sum_i \max\left(0,\thinspace q_i\right) dt,
\qquad
Q_{m,\text{out}} = \int_{\text{step}} \sum_i \max\left(0,\thinspace -q_i\right) dt,
\qquad
q_i = \pi\thinspace d_o\thinspace \Delta z_i\thinspace h_m\left(T_m(t) - T_{w,i}\right)
$$

$q_i$ is the outer-face flow of cell $i$, the same term the wall balance uses. The integral uses the trapezoidal rule over every logged point of the step, which means every accepted solver step plus any fill-in points. It is accumulated as the log is written, so thinning the stored rows does not change it. Both figures are zero for a step without a medium and in an isothermal run.

- In C++ they are `StepHeatDuty::medium_heat_in` and `medium_heat_out` [J], on each step's log (`EquationLog::GetLog(cycle, step).heat_duty`). They sit next to the energy sources' `heating_power` and `cooling_power`, so the two contributions can be read separately.
- The thermal-energy and cooling-duty KPIs (total and per step) add $Q_{m,\text{in}}$ to the heat supplied and $Q_{m,\text{out}}$ to the heat removed. `CalculateStepThermalDuty` also returns the medium's share on its own (`medium_heating`, `medium_cooling`).
- For a closed bed with insulated wall ends, $Q_{m,\text{in}} - Q_{m,\text{out}}$ equals the change in bed and wall internal energy over the step. The test `HeatingMediumDutyTest.LoggedDutyClosesTheEnergyBalance` checks this; the measured mismatch is below $10^{-4}$ of the change.

### Wall materials

Choosing a preset **Material** fills in the density, heat capacity and thermal conductivity, which are then read-only; *Custom* lets you type all three. The presets are reference values at about 20 °C (the source comments in `wall.h` cite the ASM Handbook, ASME B31.3 and MatWeb):

| Material | Density [kg m⁻³] | Heat capacity [J kg⁻¹ K⁻¹] | Thermal conductivity [W m⁻¹ K⁻¹] |
|---|---|---|---|
| SS316L | 7990 | 500 | 16.3 |
| SS304 | 7900 | 500 | 16.2 |
| SS304L | 7900 | 500 | 16.2 |
| Carbon Steel (A106) | 7850 | 490 | 51.9 |
| Carbon Steel (A53) | 7850 | 490 | 51.9 |
| Hastelloy C-276 | 8890 | 427 | 10.2 |
| Inconel 625 | 8440 | 410 | 9.8 |
| Aluminium 6061 | 2700 | 896 | 167 |
| Titanium Gr.2 | 4510 | 520 | 16.4 |

The wall's heat capacity is often not small. In the Perez (2019) laboratory column, for example, the 3.1 mm steel wall (about 0.8 kg) can store about 400 J K⁻¹, against about 150 J K⁻¹ for the 162 g of zeolite it contains.

### Minimum wall thickness (ASME VIII-1)

The engine also contains a preliminary pressure-vessel sizing check, `BedWall::CalculateMinThickness`, for cylindrical columns. It applies the ASME BPV Code Section VIII Division 1 formula UG-27 for the circumferential stress:

$$
t_{min} = \frac{P\thinspace R}{S\thinspace E - 0.6\thinspace P}
$$

where $R = d_i/2$, $P$ is the design pressure (gauge), $E$ the weld joint efficiency (1.0 for seamless or fully radiographed, 0.85 for spot, 0.70 for none) and $S$ the maximum allowable stress of the material at the design temperature. $S$ is interpolated linearly in temperature from tabulated ASME Section II Part D values for each preset material and held at the end value outside the table; *Custom* has no table, so no thickness is computed. For **vacuum service** the check uses a full atmosphere (101 325 Pa) as the design pressure in the same formula, a conservative simplification of the UG-28 external-pressure procedure. The result is flagged invalid if $S E \le 0.6 P$.

The design pressure, design temperature, weld efficiency and vacuum flag are stored in the deck with the wall. The check is available through the C++ API. It is a preliminary sizing aid, not a substitute for a code calculation.
