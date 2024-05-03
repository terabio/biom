import importlib.metadata
from pathlib import Path

CACHE = Path.home() / ".cache" / "biom" / importlib.metadata.version('biom')
