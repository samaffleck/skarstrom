"""``python -m skarstrom`` runs the command-line tool."""

import sys

from .cli import main

sys.exit(main())
