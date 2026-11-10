### Changed — CI: Windows all-tiers job has a longer time limit
- The C++ job on Windows now has a 240-minute limit (Linux and macOS keep 150). The all-tiers run on Windows took more than 100 minutes of tests after a 47-minute cold build, and the old 150-minute limit cancelled it before it finished. Pull-request fast-tier runs are not affected.
