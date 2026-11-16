## Overall Mass Balance

The overall mass balance ensures mass conservation throughout the packed bed. It is the mole balance of the gas divided by the concentration $C = P/(RT)$:

$$
\epsilon_t\left(\frac{1}{P}\frac{\partial P}{\partial t} - \frac{1}{T}\frac{\partial T}{\partial t}\right) + \frac{T}{P}\frac{\partial}{\partial z}\left(\frac{u P}{T}\right) + \frac{RT}{P}\left[\sum_{k}\rho_{b,k}\sum_{i=1}^{N}\frac{\partial q_{i,k}}{\partial t} - \sum_{i=1}^{N} S_i\right] = 0
$$

Multiplied by $C$, this is the conservative form $\epsilon_t\,\partial C/\partial t + \partial(uC)/\partial z + \sum_k \rho_{b,k}\sum_i \partial q_{i,k}/\partial t = \sum_i S_i$, which is also the form the manufactured-solution verification uses (`adsim/equations/manufactured_sources.h`). The solver carries the pressure $P$ as the unknown of this balance.

| Symbol        | Description                              | Unit           |
|---------------|------------------------------------------|----------------|
| $P$           | Pressure                                  | Pa        |
| $T$           | Gas temperature (the fluid temperature, also under local thermal non-equilibrium) | K        |
| $\epsilon_t$  | Total void fraction                         | -          |
| $u$           | Superficial velocity                      | m s⁻¹          |
| $\rho_{b,k}$  | Bulk density of adsorbent $k$ (mass of solid $k$ per bed volume) | kg m⁻³ |
| $q_{i,k}$     | Loading of component $i$ on adsorbent $k$  | mol kg⁻¹    |
| $S_i$         | Reaction source of component $i$          | mol m⁻³ s⁻¹ |
| $N$           | Number of gas components                  | -              |
| $t$           | Time                                      | s              |
| $z$           | Axial coordinate                          | m              |
| $R$           | Ideal gas constant                        | J mol⁻¹ K⁻¹ |

The sum over $k$ runs over the adsorbents of the cell's layer; a layer of one solid has a single term, $\rho_b\,\partial q_i/\partial t$. The bulk densities come from the layer's mass fractions, $\rho_{b,k} = (1-\epsilon_b)\,\rho_{mix}\,w_k$ (see [Mixed Adsorbent Layers](mixed-layers.md)). A component whose isotherm on adsorbent $k$ is inert has $q_{i,k} = 0$ and adds nothing. With a dual-rate linear driving force, $q_{i,k}$ is the total loading of the two sites, so the slow site is not counted twice.

**Reaction source.** Each reaction $r$ lives on an adsorbent $k$ and gives a rate $r_r$ per kilogram of that solid (see [Reactions](reactions.md)). Its effectiveness factor is $\eta_r$ and its stoichiometric coefficient for component $i$ is $\nu_{i,r}$. The source of component $i$ is

$$
S_i = \sum_k \sum_{r \in k} \nu_{i,r}\,\rho_{b,k}\,\eta_r\,r_r
$$

It is zero in the cells where the carrying solid is absent and in void cells, and it is zero everywhere in a deck without reactions. A reaction that changes the number of moles ($\sum_i \nu_{i,r} \ne 0$) drives the velocity field in the same way as adsorption.

## Component Mass Balance

The component mass balance accounts for axial dispersion, convection, adsorption/desorption according to the rate model selected, and reactions. Axial dispersion can be constant or calculated from different empirical equations. How $\rho_b$, $\epsilon_t$ and $\epsilon_b$ follow from the deck's parameters, and which settings to check when comparing with another code, is on [Matching Another Model](model-conventions.md#densities-and-voidages).

$$
\epsilon_t \left(\frac{\partial y_i}{\partial t} + \frac{y_i}{P}\frac{\partial P}{\partial t} - \frac{y_i}{T}\frac{\partial T}{\partial t}\right) + \frac{T}{P} \frac{\partial}{\partial z}\left(\frac{u y_i P}{T}\right) - \frac{T}{P}\frac{\partial}{\partial z}\left(\epsilon_b \frac{D_{ax} P}{T} \frac{\partial y_i}{\partial z}\right) + \frac{RT}{P} \left[\sum_k \rho_{b,k} \frac{\partial q_{i,k}}{\partial t} - S_i\right] = 0
$$

Multiplied by $C$, it is $\epsilon_t\,\partial (C y_i)/\partial t + \partial(u C y_i)/\partial z - \partial(\epsilon_b D_{ax} C\,\partial y_i/\partial z)/\partial z + \sum_k \rho_{b,k}\,\partial q_{i,k}/\partial t = S_i$. Inside a layer $\epsilon_b$ is constant, so it can be taken out of the derivative.

| Symbol        | Description                              | Unit           |
|---------------|------------------------------------------|----------------|
| $y_i$         | Mole fraction of component $i$ in gas phase | -             |
| $P$           | Pressure                                  | Pa        |
| $T$           | Gas temperature                           | K        |
| $\epsilon_t$  | Total void fraction                       | -          |
| $\epsilon_b$  | Bed void fraction (the open frontal area of a structured layer) | -          |
| $u$           | Superficial velocity                     | m s⁻¹          |
| $D_{ax}$      | Axial dispersion coefficient (one value for all components; see [below](#dispersion-with-a-coefficient-per-component) for one per component) | m² s⁻¹          |
| $\rho_{b,k}$  | Bulk density of adsorbent $k$            | kg m⁻³ |
| $q_{i,k}$     | Loading of component $i$ on adsorbent $k$ | mol kg⁻¹    |
| $S_i$         | Reaction source of component $i$         | mol m⁻³ s⁻¹ |
| $t$           | Time                                     | s              |
| $z$           | Axial coordinate                         | m              |
| $R$           | Ideal gas constant                        | J mol⁻¹ K⁻¹ |

The solver carries a component balance for **every** component, in addition to the overall balance. It does not impose $\sum_i y_i = 1$: the sum follows from the balances, to within the solver tolerance (see [Numerics](numerics.md#advection-schemes)).

**Void regions.** A [void region](getting-started.md) is one well-mixed cell with $\epsilon_t = \epsilon_b = 1$, no adsorbent and no reactions. Its dispersion coefficient is zero, so no dispersive flux crosses the faces of a void cell: gas moves between a void and the bed by convection only.

### Dispersion with a coefficient per component

!!! note "New in 0.15"
    Per-component dispersion (`component_diffusivity` = `mixture_averaged`).

With one $D_{ax}$ for all components, the dispersive fluxes $J_i = -\epsilon_b D_{ax} C\,\partial y_i/\partial z$ sum to zero, because $\sum_i y_i = 1$. That is why the overall mass balance has no dispersion term: a dispersive flux moves components relative to the molar-average velocity $u$ and carries no net moles. The molecular diffusivity is different for each component, so the Ruthven and Taylor–Aris models can give each component its own $D_{ax,i}$ ([Transport Correlations](transport-correlations.md#diffusivity-of-each-component)). Fickian fluxes with different coefficients do not sum to zero. The component balance therefore takes the self-consistent effective binary flux of Ramshaw (1990), written for the molar-average velocity:

$$
J_i = -\epsilon_b C\left(D_{ax,i}\frac{\partial y_i}{\partial z} - y_i\sum_{j=1}^{N} D_{ax,j}\frac{\partial y_j}{\partial z}\right), \qquad \sum_{i=1}^{N} J_i = 0
$$

The sum runs over all $N$ components, $i$ included. The second term subtracts from every component its share $y_i$ of the net drift $\sum_j J_j^0$ of the uncorrected fluxes. It adds the same velocity to every component, so it treats the components alike and keeps their velocities relative to each other, which are all that the multicomponent (Maxwell–Stefan) diffusion equations determine. The overall balance, the momentum balance and the meaning of $u$ stay as above. With equal coefficients the correction vanishes, and $J_i$ is the flux of the one-coefficient model. [Numerics](numerics.md#spatial-discretisation-cell-centred-finite-volumes) gives the discrete form, which sums to zero to round-off through every face.

Ramshaw derives the correction relative to the mass-average velocity (his eq. 16, with weights $\rho_i/\rho$). The momentum balance here gives the molar-average velocity, so the weights are the mole fractions. The two frames differ by a net molar flux of order $1/Pe$ against the convective flux, below 1 % at the Péclet numbers of adsorption columns. The same correction is used in combustion codes as the "correction velocity" of mixture-averaged transport (Coffee and Heimerl 1981).

- J. D. Ramshaw, "Self-consistent effective binary diffusion in multicomponent gas mixtures", J. Non-Equilib. Thermodyn. 15 (1990) 295–300, doi:10.1515/jnet.1990.15.3.295.
- T. P. Coffee and J. M. Heimerl, "Transport algorithms for premixed, laminar steady-state flames", Combust. Flame 43 (1981) 273–289, doi:10.1016/0010-2180(81)90027-4.

### Inlet and outlet conditions

At an **inlet** the boundary condition is of the Danckwerts type: the total flux of each component entering the bed (convection plus dispersion) equals the feed's,

$$
u\thinspace C\thinspace y_{i,in} = u\thinspace C\thinspace y_i\rvert_{0^+} - \epsilon_b D_{ax} C \frac{\partial y_i}{\partial z}\bigg\rvert_{0^+}
$$

so the gas in the first cell can differ from the feed, and the moles of each component the bed takes in are exactly the moles the inlet stream supplies (the amounts the stream logs, the performance metrics and the cyclic-steady-state mass balance count). At an **outlet** or a **closed** end the composition has zero gradient. The energy balance treats heat conduction at the inlet the same way. [Numerics](numerics.md#spatial-discretisation-cell-centred-finite-volumes) gives the discrete form.
