## Momentum Balance

The momentum balance relates the axial pressure gradient to the **superficial** gas velocity $u$ (volumetric flow per unit of column cross-section). Inertia of the gas and wall friction are neglected, as is usual for packed-bed adsorption columns, so the momentum balance reduces to a packed-bed pressure-drop law. Every layer selects its own law under **Layer Properties → Momentum model**:

| Momentum model | Viscous term | Inertial term | Default for |
|----------------|--------------|---------------|-------------|
| **Ergun** | ✔ (constant 150) | ✔ (constant 1.75) | new layers and new decks |
| **Darcy (Kozeny–Carman)** | ✔ (constant 180) | — | decks saved before Ergun was available |
| **Darcy (Blake–Kozeny)** | ✔ (constant 150) | — | — (choose it to match papers that use it) |

The builder and registry names are `ergun`, `kozeny_carman` and `blake_kozeny` (`layer[...].momentum_model` = 1, 0 and 2). All three are special cases of

$$
-\frac{\partial P}{\partial z} = A\thinspace u + B\thinspace \lvert u\rvert\thinspace u
$$

where $A$ [Pa s m⁻²] is the viscous (laminar) resistance and $B$ [Pa s² m⁻³] the inertial (turbulent) one.

### Darcy's law with the Kozeny–Carman permeability

$$
u = -\frac{K}{\mu}\frac{\partial P}{\partial z},
\qquad
K = \frac{d_p^2\thinspace \epsilon_b^3}{180\thinspace (1-\epsilon_b)^2}
$$

i.e. $A = \mu / K$ and $B = 0$. The pressure drop is linear in the velocity. This is accurate at low particle Reynolds number $Re_p = \rho\thinspace u\thinspace d_p / \mu$ (creeping flow) and is the model that every deck created before the Ergun option existed was solved with — such decks still load as Darcy (Kozeny–Carman) and give exactly the same results.

### Darcy's law with the Blake–Kozeny permeability

$$
u = -\frac{K}{\mu}\frac{\partial P}{\partial z},
\qquad
K = \frac{d_p^2\thinspace \epsilon_b^3}{150\thinspace (1-\epsilon_b)^2}
$$

i.e. $A = 150\thinspace \mu\thinspace (1-\epsilon_b)^2 / (\epsilon_b^3\thinspace d_p^2)$ and $B = 0$: the viscous (Blake–Kozeny) term of the Ergun equation below **without** its inertial term. This is the pressure-drop law of Haghpanah *et al.* (2013) and Estupiñan Perez *et al.* (2019, SI Eq. S5), who write it in the interstitial velocity $v = u / \epsilon_b$ and the particle radius $r_p = d_p / 2$ as

$$
-\frac{\partial P}{\partial z} = \frac{150}{4}\thinspace \frac{1}{r_p^2}\left(\frac{1-\epsilon_b}{\epsilon_b}\right)^2 \mu\thinspace v .
$$

It is solved exactly like the Kozeny–Carman law (the linear Darcy face residuals everywhere, including at boundaries), and its viscous resistance is 150/180 of the Kozeny–Carman one at every flow rate (for an isothermal ideal gas, $P_{in}^2 - P_{out}^2$ is exactly 150/180 of the Kozeny–Carman value; for a small pressure drop, so is the pressure drop itself). It equals the Ergun pressure drop in the creeping-flow limit and falls below it as the Reynolds number grows.

### Ergun equation

$$
-\frac{\partial P}{\partial z} =
\frac{150\thinspace \mu\thinspace (1-\epsilon_b)^2}{\epsilon_b^3\thinspace d_p^2}\thinspace u
{}+ \frac{1.75\thinspace \rho\thinspace (1-\epsilon_b)}{\epsilon_b^3\thinspace d_p}\thinspace \lvert u\rvert\thinspace u
$$

The first (Blake–Kozeny) term dominates at low Reynolds number; the second (Burke–Plummer) term adds the kinetic-energy losses that become important at the higher velocities of fast cycles and large particles (Ergun, *Chem. Eng. Prog.* 48 (1952) 89–94). It is the pressure-drop model used in most of the PSA literature and is therefore the default for newly created layers. (Haghpanah *et al.* (2013) and the papers built on their model keep only the viscous term: that is the Darcy (Blake–Kozeny) option above.) $\rho = \sum_i c_i M_i$ is the local gas density, so the inertial resistance follows the pressure and composition of the gas.

At low Reynolds number the inertial term vanishes and Ergun reduces to Darcy's law with the constant 150 instead of 180, so in the creeping-flow limit the Darcy (Kozeny–Carman) pressure drop is 180/150 = 1.2 times the Ergun one. Switching a layer between them is therefore a (small) physical change, not just a numerical one.

### Solving for the velocity

The velocity is an algebraic unknown on every cell face. For a given driving pressure difference $g$ across a face, the quadratic $g = A u + B\lvert u\rvert u$ has the unique physical root

$$
u = \frac{2\thinspace g}{A + \sqrt{A^2 + 4\thinspace B\thinspace \lvert g\rvert}}
$$

i.e. $u = \mathrm{sign}(g)\left(-A + \sqrt{A^2 + 4B\lvert g\rvert}\right)/(2B)$ written without the cancellation that form suffers at small $g$; it equals Darcy's $u = g/A$ exactly when $B = 0$.

The solver does not, however, impose the root in this explicit form. At high Reynolds number it behaves like $\sqrt{\lvert g\rvert/B}$ everywhere except within a fraction of a pascal of $g = 0$, and Newton iterations on such a square-root characteristic fail at the abrupt boundary switches of a PSA cycle (a closed bed end, a valve opening). Instead each face carries the momentum balance itself as its residual,

$$
r_u = u - \frac{g - B\thinspace \lvert u\rvert\thinspace u}{A} = 0
$$

which is a smooth polynomial in the unknown $u$ (slope $(A + 2B\lvert u\rvert)/A \ge 1$, continuous through flow reversal) and linear in the pressures; its solution is exactly the root above, and for $B = 0$ it is the Darcy residual. Darcy layers (Kozeny–Carman and Blake–Kozeny, $B = 0$) keep their original explicit residual, so Kozeny–Carman results are unchanged.

**Pressure and closed boundaries.** The boundary face of a **pressure** or **closed** boundary condition is the exception. These conditions fix that face's pressure drop: a pressure condition sets the ghost-cell pressure so that the face sits at the boundary pressure $P_{bc}$, giving $g = 2\thinspace (P_{bc} - P_1)$ across the face (with the sign of the flow direction), and a closed end has $g = 0$. There the face imposes the explicit root of that $g$,

$$
r_u = u - \frac{2\thinspace g}{A + \sqrt{A^2 + 4\thinspace B\thinspace \lvert g\rvert}} = 0,
$$

with $g$ taken from the boundary value $P_{bc}$ and the first bed cell, not from the ghost pressure. At a step switch the boundary-face velocity must jump (a 15 bar feed face being closed, a vent opening). With the implicit residual, the slope in $u$ changes by an order of magnitude between the old and the new velocity (from about 33 to 1 for the feed face of the 8-step biogas cycle, whose inertial term is about 16 times the viscous one), and IDA's corrector, a Newton iteration whose Jacobian is evaluated before the jump, failed to converge at every step size. The explicit root written this way is linear in $u$ with slope 1 and depends only on the first-cell pressure, which does not jump, so the corrector resolves the jump in one iteration. It has exactly the same solution, so results only change within the solver tolerance.

**Interior faces and layer interfaces.** A face between cells $i$ and $i+1$ is treated as the two half-cells between the cell centres acting in series: each contributes its own viscous resistance $\mu\thinspace \Delta z / 2K$ and inertial resistance $\beta\thinspace \rho\thinspace \Delta z / 2$ (with $\beta = 1.75(1-\epsilon_b)/(\epsilon_b^3 d_p)$ for Ergun and $0$ for both Darcy laws), evaluated with that cell's own voidage, particle diameter, viscosity and density:

$$
A = \tfrac12\left(\frac{\mu_i\thinspace \Delta z_i}{K_i} + \frac{\mu_{i+1}\thinspace \Delta z_{i+1}}{K_{i+1}}\right),\qquad
B = \tfrac12\left(\beta_i\thinspace \rho_i\thinspace \Delta z_i + \beta_{i+1}\thinspace \rho_{i+1}\thinspace \Delta z_{i+1}\right)
$$

For two Darcy cells this is exactly the harmonic-mean permeability used before, so a Darcy layer next to an Ergun layer, and an interface between two layers with different packings, are handled consistently.

**Void regions.** A void region (dead volume) is well mixed and has no packing resistance. The face between a void and the bed uses the bed-side cell's resistance over its full width — $\mu\thinspace \Delta z / K$ plus, for an Ergun layer, $\beta\thinspace \rho\thinspace \Delta z$ — in series with the optional void **resistance coefficient** (a linear Pa s m⁻¹ term, applied only at the void–bed face).

**Flow boundaries.** On a flow boundary the streams' molar flow is matched against the flow through the boundary face, which obeys the same law between the ghost cell and the first bed cell. The streams' flow may itself follow the face pressure (a valve or a pump), so the face velocity $u$ and the ghost pressure are determined together by

$$
C_f(u)\thinspace A_c\thinspace u = Fig(P_{face}(u), tig),
\qquad
P_{ghost}(u) = P_1 \pm \left(A\thinspace u + B\thinspace \lvert u\rvert\thinspace u\right)
$$

where $F$ is the streams' total molar flow, $A_c$ the cross-section and $P_1$ the first bed cell's pressure. For a Darcy layer both relations are linear and are solved as they stand. For an Ergun layer the face velocity can jump by a large factor at a step switch (a feed step opening a vent valve on a bed that the previous step pressurised above the valve's downstream pressure asks for a surge of several m/s from a standstill), across which the implicit residual's slope in $u$ changes by an order of magnitude and IDA's corrector failed. The flow boundary of an Ergun layer therefore solves the scalar face problem above for $u^\*$ inside the residual (a bracketed root find; the function increases monotonically with $u$) and imposes

$$
r_u = u - u^\*,
\qquad
r_P = P_{ghost} - P_{ghost}(u^\*),
$$

whose Jacobian block is the identity: $u^\*$ depends only on the first bed cell's state, which is continuous through the switch, the ghost temperature and the time, so the corrector resolves the jump in one iteration. The solution is the same, so results change only within the solver tolerance.

| Symbol        | Description                              | Unit           |
|---------------|------------------------------------------|----------------|
| $P$           | Pressure                                 | Pa             |
| $u$           | Superficial velocity                     | m s⁻¹          |
| $\mu$         | Gas dynamic viscosity                    | Pa·s           |
| $\rho$        | Gas density                              | kg m⁻³         |
| $K$           | Permeability of packed bed               | m²             |
| $A$, $B$      | Viscous / inertial resistance            | Pa s m⁻², Pa s² m⁻³ |
| $\beta$       | Ergun inertial coefficient               | m⁻¹            |
| $z$           | Axial coordinate                         | m              |
| $d_p$         | Particle diameter                        | m              |
| $\epsilon_b$  | Bed void fraction                        | -              |

### Approach to fluidisation

The **Approach To Fluidisation** result variable is the upward pressure gradient divided by the weight of the packing less the gas, $(-\partial P/\partial z)/((\rho_b - \rho_g) g)$. Above 1 an unrestrained bed in upward flow lifts; the solver warns at the end of a run when any step of the last cycle passes it (see [Wall Design and Operating Limits](wall-design.md#operating-limit-warnings)).
