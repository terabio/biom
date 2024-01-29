from pathlib import Path

from . import __version__

CACHE = Path.home() / ".cache" / "biom" / __version__
