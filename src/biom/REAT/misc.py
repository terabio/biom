from typing import Optional

import numpy as np
from dataclasses import dataclass

DTYPES_SITES = {
    'contig': str, 'refnuc': str, 'prednuc': str, 'trstrand': str,
    'pos': np.uint32, "A": np.uint32, "C": np.uint32, "G": np.uint32, "T": np.uint32
}

DTYPES_ROIS = {
    "contig": str, "start": np.uint32, "end": np.uint32, "strand": str, "name": str, "trstrand": str,
    "coverage": np.uint32, "nucmasked": np.uint32, "heterozygous": np.uint32, "#A": np.uint32, "A->A": np.uint32,
    "A->C": np.uint32, "A->G": np.uint32, "A->T": np.uint32, "#C": np.uint32, "C->A": np.uint32, "C->C": np.uint32,
    "C->G": np.uint32, "C->T": np.uint32, "#G": np.uint32, "G->A": np.uint32, "G->C": np.uint32, "G->G": np.uint32,
    "G->T": np.uint32, "#T": np.uint32, "T->A": np.uint32, "T->C": np.uint32, "T->G": np.uint32, "T->T": np.uint32,
}


@dataclass(frozen=True)
class SitesThreshold:
    minsamples: Optional[int] = None
    minfreq: Optional[float] = None
    minedits: Optional[int] = None
    mincoverage: Optional[int] = None


@dataclass(frozen=True)
class EIslandThreshold:
    dbscan_eps: int
    min_sites: int

    min_strong_sites: int
    strongthr: SitesThreshold
    weakthr: SitesThreshold

    fold_offset: int = 500
