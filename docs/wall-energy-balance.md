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

This condition only concerns the ends of the wall. Heat loss through the side of the column is always governed by $h_o$ and $T_{amb}$.

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
