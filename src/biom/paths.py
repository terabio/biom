from pathlib import Path

from . import __version__

CACHE = Path.home().joinpath(".local", "biom", __version__)
