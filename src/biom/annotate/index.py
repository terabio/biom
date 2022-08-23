import pickle
from functools import lru_cache
from pathlib import Path
from typing import List, Tuple

from HTSeq import GenomicArrayOfSets, GenomicInterval
from pybedtools import BedTool


def onthefly(bed: Path) -> GenomicArrayOfSets:
    index = GenomicArrayOfSets('auto', stranded=False)
    for interval in BedTool(bed):
        key = GenomicInterval(interval.chrom, interval.start, interval.end)
        index[key] += (interval.name, interval.strand)
    return index


def bedfile(file: Path, saveto: Path):
    index = GenomicArrayOfSets('auto', stranded=False)
    for interval in BedTool(file):
        key = GenomicInterval(interval.chrom, interval.start, interval.end)
        index[key] += (interval.name, interval.strand)
    with open(saveto, 'wb') as stream:
        pickle.dump(index, stream, protocol=pickle.HIGHEST_PROTOCOL)


def transcriptome(files: List[Tuple[str, Path]], saveto: Path):
    index = GenomicArrayOfSets('auto', stranded=False)
    for feature, path in files:
        for interval in BedTool(path):
            key = GenomicInterval(interval.chrom, interval.start, interval.end)
            index[key] += (feature, interval.name, interval.strand)
    with open(saveto, 'wb') as stream:
        pickle.dump(index, stream, protocol=pickle.HIGHEST_PROTOCOL)


@lru_cache(maxsize=1_000)
def load(file: Path) -> GenomicArrayOfSets:
    with open(file, 'rb') as stream:
        index = pickle.load(stream)
    return index
