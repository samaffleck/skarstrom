### Fixed — light-product pressurisation draws from a donor below the receiver
- Gas through a connection now flows only from the higher to the lower pressure (#241). Before, a
  producing step's light-product stream (a flow stream with no device) replayed what the
  pressurisation step took in, whatever the pressures, so a producing bed that had fallen below the
  receiving bed went on supplying it. Now a connected pressure boundary caps its imposed pressure at
  its partner's logged pressure at the same instant (an inlet cannot rise above the bed that feeds
  it, an outlet cannot fall below the bed it feeds), so a connected ramp may not reach its target.
  A flow stream with no device scales the flow it replays by the pressure difference in its
  direction, so a producing bed stops giving light product when it falls below the bed it feeds.
  Both close smoothly over 1 kPa, and the cap sits 1 kPa inside the partner's pressure, so once the
  cycles settle the moles through the connection balance. Results of light-product decks change
  slightly (the receiving bed now ends 1 kPa below the producing bed).

### Added — connection pressure warning
- A new end-of-run warning, `connection_pressure` (`settings.warn_connection_pressure`, on by
  default), reports gas that enters a bed through a connection while that bed is the higher (more
  than 0.1 % of the cycle feed in one step), and a connected pressure profile that the cap held
  more than 2 kPa from its own value.
