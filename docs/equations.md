# Governing Equations

Skarstrom is a one-dimensional adsorption process simulator that uses a finite volume discretisation scheme to solve the governing equations. The governing equations are composed of coupled non-linear partial differential equations (PDEs), ordinary differential equations (ODEs) in each cell, and algebraic equations.

## Partial Differential Equations (PDEs)

The simulation is based on five conservation equations:

1. [Overall Mass Balance](mass-balance.md) - Ensures mass conservation throughout the packed bed, including adsorption and reaction sources
2. [Momentum Balance](momentum-balance.md) - Relates the pressure gradient to the velocity with the Ergun equation or Darcy's law (Kozeny–Carman or Blake–Kozeny permeability) in a packed layer, and the laminar Hagen–Poiseuille law in a structured (monolith or laminate) layer. It has no time derivative, so the face velocities are algebraic unknowns
3. [Component Mass Balance](mass-balance.md) - Accounts for axial dispersion, convection, adsorption and reactions
4. [Energy Balance](energy-balance.md) - Thermal balance between gas and adsorbent phases
5. [Wall Energy Balance](wall-energy-balance.md) - Heat transfer between wall and surroundings

The energy balance can be solved assuming local thermal equilibrium (a single lumped gas/solid temperature) or local thermal non-equilibrium (separate gas and solid temperatures) — see the [energy balance](energy-balance.md).

## Ordinary Differential Equations in each cell

- [Kinetic Rate Model](kinetic-models.md) - The rate of change of each loading, $\partial q_{i,k}/\partial t$ (linear driving force and its variants). Each loading is a differential unknown in every cell of a layer that contains the adsorbent

## Algebraic Equations

Several algebraic relationships close the system:

- [Equilibrium Isotherm](isotherms.md) - Adsorption equilibrium relationship
- Mass-transfer coefficients of the [kinetic models](kinetic-models.md)
- [Reactions](reactions.md) - Reaction rates, which give the reaction sources of the mass and energy balances
- [Transport Correlations](transport-correlations.md) - Axial dispersion, thermal conductivity and heat-transfer coefficients
- [Equation of State](equation-of-state.md) - Ideal gas law and gas-mixture properties, including molecular diffusivities

## Numerical Solution

The PDEs are discretised using a finite volume method in space. The resulting system of coupled differential-algebraic equations (DAEs) is integrated in time with a variable-order (BDF), variable-step implicit solver (SUNDIALS IDA). At each step a non-linear system is solved by a Newton iteration. The [Numerics](numerics.md) page describes the discretisation, the integrator settings, the Jacobian and cyclic steady state in detail; [Running the Solver](solver-settings.md) covers the settings you can change.
