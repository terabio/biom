from dataclasses import dataclass
from pathlib import Path
from typing import Generic, List, Set, TypeVar

from pysam import AlignmentFile

T = TypeVar('T')


@dataclass()
class Stranded(Generic[T]):
    fwd: T
    rev: T


def fetch_contigs(inbam: List[Path]) -> Set[str]:
    contigs = set()
    for b in inbam:
        for contig in AlignmentFile(b, 'rb').references:
            contigs.add(contig)
    return contigs
