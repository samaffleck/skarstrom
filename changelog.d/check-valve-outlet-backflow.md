### Added — check valve on valve streams, and a warning for backflow through an outlet
- Valve streams take an opt-in check valve (`endpoint[ID].valve.check_valve`, Python
  `sk.ValveStream(..., check_valve=True)`, deck `"Check_Valve"`). An outlet then passes gas only out
  of the bed and an inlet only into it. The closure is smoothed over 10 Pa, so the boundary residual
  stays differentiable as the flow decays through zero. It applies to unconnected streams only;
  connected valves (pressure equalisation, light-product pressurisation) still carry flow both ways.
  Off by default: existing decks keep their results.
- A new end-of-run warning, `outlet_backflow` (`settings.warn_outlet_backflow`, on by default), is
  given when gas flows back into the bed through an outlet with no connection to another step. It
  fires when more than 0.1 % of the cycle feed flows back in through one outlet in one step of the
  last logged cycle. That gas carries the bed-end composition, not the stream's. The warning names
  the step, the endpoint, the moles and the fraction of the feed, and suggests the check valve.

### Fixed — the docs said an outlet's `y` and temperature enter on backflow
- `docs/python.md` said that gas flowing back in through an outlet takes the outlet's `y` and
  temperature. It takes the bed-end composition and temperature (an outlet is zero-gradient), as
  `docs/boundary-conditions.md` says.
