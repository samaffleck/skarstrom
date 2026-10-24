## Wall Design and Operating Limits

The column wall's thickness can be calculated from a design code instead of typed in, and every run is checked against the wall's design and against fluidisation of the bed. The wall design is set on **Column Properties → Wall → Pressure design**; the calculation is in `core/include/adsim/objects/wall_design.h` and the run checks in `core/include/adsim/solver/operating_limits.h`.

!!! warning "Preliminary sizing"
    The calculation is for preliminary sizing. It does not replace a design to the edition of the code in force, checked by a qualified engineer: heads, nozzles, supports, stiffening rings, fatigue from pressure cycling and local loads are not considered.

### Design inputs

| Input | Description | Unit | Registry path |
|---|---|---|---|
| Design code | ASME BPVC VIII-1 (pressure vessel shell) or ASME B31.3 (process piping, for tube and pipe columns) | | `wall.design_code` |
| Design pressure $P$ | Maximum internal pressure, gauge | Pa | `wall.design_pressure` |
| Design temperature | The allowable stress (and, for vacuum, the modulus) is taken here | K | `wall.design_temperature` |
| Joint efficiency $E$ | 1.0 seamless or fully radiographed, 0.85 spot, 0.70 none | – | `wall.weld_joint_efficiency` |
| Corrosion allowance $c$ | Added to the thickness the pressure needs | m | `wall.corrosion_allowance` |
| Vacuum service | Also design for full vacuum (101 325 Pa outside) | | `wall.vacuum_service` |

The material is the wall's **Material**; the allowable stress $S$ comes from the ASME II-D tables for the listed materials, so a *Custom* wall cannot be sized. The formulas are for cylindrical columns. The design is **in use** when the design pressure is positive or vacuum service is set; only then are the design checks below made.

**Calculate** sets **Thickness** to the required value. Python: `sim.wall_design().required_thickness`.

### Internal pressure

**ASME BPVC VIII-1, UG-27(c)(1)** (circumferential stress), with $R$ the inside radius in the corroded condition, $R = d_i/2 + c$:

$$
t = \frac{P R}{S E - 0.6 P} + c
$$

Above $P = 0.385\,S E$ (a thick shell) the rule of **Appendix 1-2(a)(1)** is used instead:

$$
t = R\left(\sqrt{Z} - 1\right) + c, \qquad Z = \frac{S E + P}{S E - P}
$$

**ASME B31.3, 304.1.2** eq. (3b), with the inside diameter $d$, $Y = 0.4$ and the weld strength reduction factor $W = 1$:

$$
t = \frac{P\,(d + 2c)}{2\left(S E W - P(1 - Y)\right)} + c
$$

With $Y = 0.4$ this equals UG-27 for a thin wall; B31.3 has no thick-wall formula, so above $P/(SE) = 0.385$ or $t \ge D/6$ (304.1.2(b)) the calculation stops. B31.3's own Table A-1 stresses equal the ASME II-D values used here for the austenitic steels and are higher for carbon steel, so the result is conservative for carbon steel. Above 427 °C (carbon steel) or 510 °C (other materials) B31.3 requires $W < 1$, which is not applied; the result notes it.

### External pressure (vacuum service)

A wall under vacuum fails by buckling, not by hoop stress, so 1 atm as an internal pressure in UG-27 is far too thin. The allowable external pressure $P_a$ is calculated by the formula method of **ASME Code Case 2286** (the same rules as **VIII-2 4.4.5.1**), which gives UG-28's result without its charts. With the corroded thickness $t$, the outside diameter $D_o$ and the column length $L$:

$$
M_x = \frac{L}{\sqrt{R_o t}}, \qquad F_{he} = 1.6\,C_h\,E_y\,\frac{t}{D_o}
$$

$$
C_h = \begin{cases}
0.55\,t/D_o & M_x \ge 2\,(D_o/t)^{0.94} \\
1.12\,M_x^{-1.058} & 13 < M_x < 2\,(D_o/t)^{0.94} \\
0.92/(M_x - 0.579) & 1.5 < M_x \le 13 \\
1 & M_x \le 1.5
\end{cases}
$$

The elastic buckling stress $F_{he}$ becomes the predicted buckling stress $F_{ic}$: $F_{ic} = S_y$ when $F_{he}/S_y \ge 2.439$, $0.7\,S_y\,(F_{he}/S_y)^{0.4}$ when $0.552 < F_{he}/S_y < 2.439$, and $F_{he}$ otherwise. The design factor is $FS = 2.0$ for $F_{ic} \le 0.55\,S_y$, $2.407 - 0.741\,F_{ic}/S_y$ up to $S_y$, and 1.667 at $S_y$. Then

$$
P_a = 2\,\frac{F_{ic}}{FS}\,\frac{t}{D_o}
$$

and the required thickness is the smallest with $P_a \ge 101\,325$ Pa (plus $c$). With vacuum service the thickness is the larger of the internal and the external requirement.

- $E_y$ is the elastic modulus at the design temperature, approximately as in ASME II-D Tables TM-1, TM-2, TM-4 and TM-5.
- $S_y$ is a lower bound on the yield strength at temperature from the allowable stress: II-D limits $S$ to $\tfrac{2}{3}S_y$ (0.9 $S_y$ for austenitic steels and nickel alloys), so $S_y \ge 1.5\,S$ ($S/0.9$). A lower $S_y$ gives a lower $P_a$, so this is conservative.
- $L$ is the packed length of the column (the sum of the layer lengths). The heads (UG-28 adds a third of each head's depth) and any stiffening rings are not modelled.

For an SS316L column at 20 °C: 100 mm × 1 m needs 0.79 mm for full vacuum, 300 mm × 1 m needs 1.50 mm, and 300 mm × 2 m needs 2.01 mm (UG-27 at 1 atm gives 0.04 mm and 0.13 mm).

### Operating-limit warnings

At the end of a run that completes, the solver checks the **last logged cycle** step by step and logs one warning per step and limit (Python: `sim.operating_limits()`; CLI: `warnings` in `skarstrom run --json`):

| Warning | When | Checked |
|---|---|---|
| Design pressure | The highest pressure in the column (packed and void cells) less 101 325 Pa exceeds the design pressure | design in use |
| Vacuum | The pressure falls below atmospheric and the wall is not in vacuum service | design in use |
| Design temperature | The wall temperature exceeds the design temperature | design in use |
| Fluidisation | The approach to fluidisation exceeds 1 | always |

Each warning names the step and cycle; the fluidisation warning also gives the height above the bottom of the bed and the time into the step at which the approach peaks. Before a run, the solver also warns when the wall is thinner than its design requires (for example after the diameter was increased) or the requirement cannot be evaluated.

The **approach to fluidisation** is the upward pressure gradient over the weight of the packing less the gas,

$$
\frac{-\partial P/\partial z}{\left(\rho_b - \rho_g\right) g}
$$

(the **Approach To Fluidisation** result variable). Above 1, the pressure drop of an upward flow exceeds the weight of the bed and an unrestrained bed lifts: in feed steps fed from the bottom, and often in pressurisation from the bottom or blowdown through the top, where the flow is fastest. It is evaluated on the faces between packed cells. Downward flow presses the bed down and gives negative values. If the bed is held down (a hold-down plate or spring), the warning can be ignored; otherwise lower the flow, slow the pressure change, or reverse the step's flow direction.

Values are those at the logged times, so a peak between two logged points is not seen. A step's first logged point is the state the step before ended in: an extreme a step reaches only there gives no warning, because the step before reports it.
