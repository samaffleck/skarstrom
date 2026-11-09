### Fixed — light-product pressurisation draws from a donor below the receiver
- A connected flow stream with no device (for example the light-product stream of a producing step)
  replayed its partner's logged flow whatever the pressures. A donor bed that had fallen below the
  receiving bed went on supplying it, so gas flowed through the connection from the lower to the
  higher pressure. The replayed flow is now scaled by a smooth weight of the pressure difference at
  the same instant: in full above 10 Pa in the flow's direction, none at or below zero (#241).
