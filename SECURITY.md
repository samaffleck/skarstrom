# Security Policy

## Supported versions

Security fixes are made on the latest release and on `main`. Skarstrom is
pre-1.0; older minor versions are not maintained.

## Reporting a vulnerability

Please **do not open a public issue** for a security problem. Report it
privately through GitHub's
[private vulnerability reporting](https://github.com/samaffleck/skarstrom/security/advisories/new)
(Security tab → "Report a vulnerability"), or by email to sam@skarstrom.co.uk.

Include what you found, how to reproduce it (a deck file or code snippet helps)
and the version or commit you tested. You should get an acknowledgement within
a week. Once a fix is available we will publish an advisory and credit you,
unless you prefer otherwise.

## Scope

Skarstrom is a simulation library that reads deck files (JSON). Relevant
issues include memory-safety bugs or crashes reachable from a crafted deck or
from the Python API, and problems in the build or release pipeline (for
example the published wheels). Inaccurate simulation results are bugs, not
security issues: please report them as normal issues.
