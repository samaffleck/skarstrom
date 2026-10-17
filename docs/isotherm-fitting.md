# Isotherm Fitting

Skarstrom fits [isotherm models](isotherms.md) to measured equilibrium points with the engine's own isotherm code, so the fitted parameters mean exactly what they will mean in a simulation. It needs no column, cycle or deck: it is the quick way from a table of uptake data to model parameters. From Python the fitter is `skarstrom.fit_isotherm`; in C++ it is `adsim::isofit::Fit` (`adsim/solver/isotherm_fit.h`). The executed notebook `examples/notebooks/isotherm_fit.ipynb` walks through single- and multi-temperature fits.

## Data

Each point is a pressure $p$ [Pa], a loading $q$ [mol/kg] and a temperature $T$ [K], optionally with a mole fraction $y$ (the component's partial pressure is then $p\,y$; pure-component data use $y = 1$). Everything is in SI base units: convert first (1 bar = 10⁵ Pa, 1 mmol/g = 1 mol/kg; 22.414 cm³(STP)/g = 1 mmol/g at 0 °C and 1 atm).

```python
import numpy as np
import skarstrom

p = np.array([0.1, 0.25, 0.5, 1.0, 2.0, 4.0]) * 1e5          # Pa
q = np.array([1.06, 1.94, 2.76, 3.47, 4.01, 4.36])            # mol/kg
fit = skarstrom.fit_isotherm("single_site_langmuir_p", p, q, temperature=298.15)
fit["params"], fit["std_errors"], fit["r2"], fit["rmse"]
```

To fit the temperature dependence of a model, pass data at several temperatures (one temperature per point):

```python
fit = skarstrom.fit_isotherm("dual_site_langmuir_p", p_all, q_all, temperature=T_all)
```

## Models

`skarstrom.isotherm_kinds()` lists the kinds the fitter offers: Henry; single- and dual-site Langmuir (pressure and concentration forms); extended single- and dual-site Langmuir; Freundlich; Sips; dual-site Sips; Toth; and the temperature-dependent Sips 2, extended Sips 2 and Toth 2. `skarstrom.isotherm_parameters(kind)` gives a kind's parameter names, which are the names used in deck files and in parameter paths (`adsorbent[<id>].isotherm[<component>].<name>`). The [Isotherms](isotherms.md) page has the definitions.

`skarstrom.evaluate_isotherm(kind, params, pressure, temperature, mole_fraction=1.0)` evaluates any kind with given parameters (unspecified parameters keep the engine default), for plotting or checking a fit.

## Starting values and bounds

The fitter derives a starting value and search bounds for every parameter from the data. You can override them:

```python
fit = skarstrom.fit_isotherm(
    "sips", p, q, temperature=303.15,
    fixed={"v": 1.0},                    # hold parameters at given values
    initial={"q_sat": 5.0},              # starting guesses
    bounds={"b": (1e-7, 1e-3)},          # search bounds
)
```

- **Bounds** span several decades around what the data can resolve. Positive scale parameters, such as capacities and affinities, are searched on a log scale; exponents and energies on a linear scale.
- **Arrhenius pairs.** For terms of the form $b\,e^{Q/RT}$ (Henry $k$/$Q$, Langmuir $b_k$/$q_k$, Sips 2 $B_1$/$B_2$), the pre-exponential is searched as its value at the reference temperature, $b\,e^{Q/(R T_{ref})}$, and its bounds apply to that value. $T_{ref}$ is the mean data temperature (returned as `t_ref`). This separates the affinity from the heat and makes the search much better conditioned.
- **Single temperature.** With data at one temperature the temperature-dependence parameters ($Q$, $q_k$, $B_2$, $A_2$, $C_2$, $x$, $\Delta H_0$, $\alpha$) cannot be determined, so they are fixed at 0, which makes the affinity at the data temperature equal to the fitted pre-exponential. If you know the heat from elsewhere, pass it in `fixed`.
- The reference temperatures `T_ref` (Sips 2) and `T_0` (Toth 2) are fixed at the mean data temperature.

## The fit

The fitter minimises the sum of squared loading residuals

$$
\mathrm{SSE} = \sum_{i=1}^{n} \left(q^{*}(p_i, T_i) - q_i\right)^2
$$

over the free parameters. It first runs a global search over the bounded box (dlib's MaxLIPO with a trust region; `max_global_evals` sets its budget, by default 100 + 80 per free parameter), then refines the best few candidates with a bounded Levenberg–Marquardt least-squares search. The model $q^{*}$ is evaluated by the same C++ isotherm functions the simulator uses.

The result dictionary holds:

- `params` and `std_errors` — each parameter with its ± 1σ standard error. This is the usual linearised estimate $\mathrm{Cov} = s^2 (J^\top J)^{-1}$, where $s^2 = \mathrm{SSE}/(n-p)$ and $J$ is the Jacobian of the model with respect to the free parameters at the optimum. Treat it as indicative: it assumes independent errors of equal variance and a model that is close to linear near the optimum. It is `nan` when the parameters are too strongly correlated to separate.
- `r2` $= 1 - \mathrm{SSE}/\sum_i (q_i - \bar q)^2$, `rmse` $= \sqrt{\mathrm{SSE}/n}$ in mol/kg, `sse`, the number of points `n` and the degrees of freedom `dof` $= n - p$.
- `q_fit`, the model at each input point, and `t_ref`.
- `warning`, for example a parameter ending on a bound. If that happens, widen the bound and fit again.

## Using the result in a deck

Write the fitted values into an adsorbent of a deck through its parameter paths (the component's isotherm must already be of the fitted kind):

```python
sim = skarstrom.Simulation("my_deck.json")
for name, value in fit["params"].items():
    sim[f"adsorbent[1].isotherm[CarbonDioxide].{name}"] = value
sim.save("my_deck_fitted.json")
```

Only the isotherm is fitted. The **heat of adsorption** used by the [energy balance](energy-balance.md) is a separate adsorbent property (`adsorbent[<id>].heat_of_adsorption[<component>]`); set it from the fitted energy or from a calorimetric value.
