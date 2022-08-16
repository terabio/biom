from . import __version__
from pathlib import Path

CACHE = Path.home().joinpath(".local", "biomodules", __version__)
