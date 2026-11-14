### Added — look up a step's endpoint summary by endpoint ID
- `Step::FindSummaryForEndpoint(size_t endpoint_id)` finds a stream's or pressure boundary's cycle summary by its stable endpoint ID, which does not change when a stream is renamed. The name overload is unchanged.
