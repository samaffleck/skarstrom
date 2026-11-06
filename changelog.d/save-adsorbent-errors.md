### Fixed — saving an adsorbent file fails silently
- `SaveAdsorbentToJSON` now throws `std::runtime_error` naming the file and
  the reason when it cannot write it (missing or read-only directory, full
  disk, a rename that keeps failing), instead of returning without a file.
  The temporary file is removed and the target is left as it was. The
  rename is retried for about 5 s, as on Windows it fails while another
  process has the target open.
