from pathlib import Path
from typing import Optional

import pandas as pd


class ContigsInfo():
    def __init__(self, cache: Path):
        df = pd.read_csv(cache, sep='\t', header=None, names=['from', 'to'])

        # Build indexes
        self._ensembl2ucsc = df.set_index('from')['to'].to_dict()
        self._ucsc2ensembl = df.set_index('to')['from'].to_dict()

    def ensembl2ucsc(self, contig: str) -> Optional[str]:
        return self._ensembl2ucsc.get(contig, None)

    def ucsc2ensembl(self, contig: str) -> Optional[str]:
        return self._ucsc2ensembl.get(contig, None)
