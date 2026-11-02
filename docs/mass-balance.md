## Overall Mass Balance

The overall mass balance ensures mass conservation throughout the packed bed.

$$
\epsilon_t(\frac{1}{P}\frac{\partial P}{\partial t} - \frac{1}{T}\frac{\partial T}{\partial t}) + \frac{T}{P}\frac{\partial}{\partial z}(\frac{u P}{T}) + \frac{RT}{P}\rho_b\sum_{i=1}^{N}\frac{\partial{q_i}}{\partial{t}} = 0
$$

| Symbol        | Description                              | Unit           |
|---------------|------------------------------------------|----------------|
| $P$           | Pressure                                  | Pa        |
| $T$           | Temperature                               | K        |
| $\epsilon_t$  | Total void fraction                         | -          |
| $u$           | Superficial velocity                      | m s⁻¹          |
| $q_i$         | Component i mass source                   | mol kg⁻¹ s⁻¹    |
| $t$           | Time                                      | s              |
| $z$           | Axial coordinate                          | m              |
| $R$           | Ideal gas constant                        | J mol⁻¹ K⁻¹ |

## Component Mass Balance

The component mass balance accounts for axial dispersion, convection, and adsorption/desorption accordnig to the rate model selected. Axial dispersion can be constant or calculated from different empirical equations. 

$$
\epsilon_t (\frac{\partial y_i}{\partial t} + \frac{y_i}{P}\frac{\partial P}{\partial t} - \frac{y_i}{T}\frac{\partial T}{\partial t}) + \frac{T}{P} \frac{\partial}{\partial z}(\frac{u y_i P}{T}) - \epsilon_b \frac{T}{P}\frac{\partial}{\partial z}(\frac{D_{ax} P}{T} \frac{\partial y_i}{\partial z}) + \frac{RT}{P} \rho_b \frac{\partial q_i}{\partial t} = 0
$$

| Symbol        | Description                              | Unit           |
|---------------|------------------------------------------|----------------|
| $y_i$         | Mole fraction of component $i$ in gas phase | -             |
| $P$           | Pressure                                  | Pa        |
| $T$           | Temperature                               | K        |
| $\epsilon_b$  | Bed void fraction                        | -          |
| $u$           | Superficial velocity                     | m s⁻¹          |
| $D_{ax}$         | Axial disperion coefficient              | m² s⁻¹          |
| $q_i$         | Component i mass source                   | mol kg⁻¹ s⁻¹    |
| $t$           | Time                                     | s              |
| $z$           | Axial coordinate                         | m              |
| $R$           | Ideal gas constant                        | J mol⁻¹ K⁻¹ |

### Inlet and outlet conditions

At an **inlet** the boundary condition is of the Danckwerts type: the total flux of each component entering the bed (convection plus dispersion) equals the feed's,

$$
u\thinspace C\thinspace y_{i,in} = u\thinspace C\thinspace y_i\rvert_{0^+} - \epsilon_b D_{ax} C \frac{\partial y_i}{\partial z}\bigg\rvert_{0^+}
$$

so the gas in the first cell can differ from the feed, and the moles of each component the bed takes in are exactly the moles the inlet stream supplies (the amounts the stream logs, the performance metrics and the cyclic-steady-state mass balance count). At an **outlet** or a **closed** end the composition has zero gradient. The energy balance treats heat conduction at the inlet the same way. [Numerics](numerics.md#spatial-discretisation-cell-centred-finite-volumes) gives the discrete form.
