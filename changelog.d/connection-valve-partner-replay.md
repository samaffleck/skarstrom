### Fixed — pressure equalisation through a valve made gas
- A flow stream with no device that is connected to a valve stream now replays the valve's flow as it
  is, in both directions. Since connections pass gas only from the higher to the lower pressure, such
  a stream also scaled the valve's flow by the pressure difference. The valve already sets its flow
  from that difference, so after the two beds equalised the valve drew gas back from a receiving bed
  a few Pa above it, and the receiving bed gave none of it: the connection made gas (2 % of the
  equalisation moles on `Vysyaraju2024_PSA_exp1`, whose pressurisation uptake fell from 0.1681 to
  0.1666 mol). The moles through such a connection balance again, and that deck's results return to
  their values before the change.
