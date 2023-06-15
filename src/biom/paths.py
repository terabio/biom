from pathlib import Path

from . import __version__

CACHE = Path.home() / ".local" / "biom" / __version__
