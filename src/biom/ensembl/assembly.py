from pathlib import Path
from typing import Dict

from biom.paths import CACHE

from . import data, fetch


class EnsemblAssembly:
    def __init__(self, name: str, organism: str, version: int = 110):
        self.name: str = name
        self.organism: str = organism
        self.version: int = version

        cache = CACHE.joinpath("ensembl", self.name, str(self.version))
        self._cache: Dict[str, Path] = {}

        self._cache['transcripts'] = cache.joinpath("transcripts-info.tsv.gz")
        if self._cache['transcripts'].exists():
            self.transcripts = data.TranscriptsInfo(self._cache['transcripts'])
            self.genes = data.GenesInfo(self._cache['transcripts'])

        self._cache['contigs'] = cache.joinpath("Ensembl2UCSC.txt")
        if self._cache['contigs'].exists():
            self.contigs = data.ContigsInfo(self._cache['contigs'])

    def fetch(self, what: str, force: bool = False) -> 'EnsemblAssembly':
        assert what in {"transcripts"}

        fetch.transcripts(self.organism, self.version, self._cache['transcripts'], force)

        # Reload the class
        return EnsemblAssembly(self.name, self.organism, version=self.version)
