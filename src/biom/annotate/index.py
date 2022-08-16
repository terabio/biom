import pickle
from collections import defaultdict
from functools import lru_cache
from typing import Dict, List, Set, Tuple

from HTSeq import GenomicArrayOfSets, GenomicInterval
from pybedtools import BedTool

from pathlib import Path


def onthefly(bed: Path) -> GenomicArrayOfSets:
    index = GenomicArrayOfSets('auto', stranded=False)
    for interval in BedTool(bed):
        key = GenomicInterval(interval.chrom, interval.start, interval.end)
        index[key] += (interval.geneid, interval.strand)
    return index


def bedfile(file: Path, saveto: Path):
    index = GenomicArrayOfSets('auto', stranded=False)
    for interval in BedTool(file):
        key = GenomicInterval(interval.chrom, interval.start, interval.end)
        index[key] += (interval.geneid, interval.strand)
    with open(saveto, 'wb') as stream:
        pickle.dump(index, stream, protocol=pickle.HIGHEST_PROTOCOL)


def transcriptome(files: List[Tuple[str, Path]], saveto: Path):
    index = GenomicArrayOfSets('auto', stranded=False)
    for feature, path in files:
        for interval in BedTool(path):
            key = GenomicInterval(interval.chrom, interval.start, interval.end)
            index[key] += (feature, interval.geneid, interval.strand)
    with open(saveto, 'wb') as stream:
        pickle.dump(index, stream, protocol=pickle.HIGHEST_PROTOCOL)


@lru_cache(maxsize=1_000)
def load(file: Path) -> GenomicArrayOfSets:
    with open(file, 'rb') as stream:
        index = pickle.load(stream)
    return index
