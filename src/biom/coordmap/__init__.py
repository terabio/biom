from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from typing import Dict, List

from pybedtools import BedTool, Interval


@lru_cache(maxsize=None)
def trexons(CDS: Path) -> Dict[str, List[Interval]]:
    exons = defaultdict(list)
    for interval in BedTool(CDS):
        exons[interval.name].append(interval)
    exons = {k: sorted(v, key=lambda x: x.start) for k, v in exons.items()}
    for _, trexons in exons.items():
        # exons must be sorted & without overlaps
        assert all(exb.end <= exa.start for exb, exa in zip(trexons[:-1], trexons[1:]))
        # same strand for all exons
        strand = trexons[0].strand
        assert all(x.strand == strand for x in trexons[1:]) and strand != "."
    return exons


def genome2transcriptome(exons: List[Interval], intervals: List[Interval]) -> List[Interval]:
    assert len(exons) >= 1 and len(intervals) >= 1
    assert all(exons[0].chrom == x.chrom for x in exons)
    assert all(exons[0].chrom == x.chrom for x in intervals)
    assert all(curex.end <= nexex.start for curex, nexex in zip(exons[:-1], exons[1:]))
    assert all(curi.end <= nexi.start for curi, nexi in zip(intervals[:-1], intervals[1:]))

    mapped, exind, trcoord = [], 0, 0
    for interval in intervals:
        # Skip exons if needed
        while interval.start >= exons[exind].end:
            trcoord += exons[exind].length
            exind += 1
            continue

        assert exons[exind].start <= interval.start < interval.end <= exons[exind].end
        curstart = trcoord + (interval.start - exons[exind].start)
        curend = trcoord + (interval.end - exons[exind].start)
        assert trcoord <= curstart < curend <= trcoord + exons[exind].length
        mapped.append(Interval("transcript", curstart, curend, strand=interval.strand))
    return mapped


def transcriptome2genome(exons: List[Interval], intervals: List[Interval]) -> List[Interval]:
    assert len(exons) >= 1 and len(intervals) >= 1
    assert all(exons[0].chrom == x.chrom for x in exons)
    assert all(intervals[0].chrom == x.chrom for x in intervals)
    assert all(curex.end <= nexex.start for curex, nexex in zip(exons[:-1], exons[1:]))
    assert all(curi.end <= nexi.start for curi, nexi in zip(intervals[:-1], intervals[1:]))

    mapped, exind, trcoord = [], 0, 0
    for interval in intervals:
        # Skip exons if needed
        while interval.start >= trcoord + exons[exind].length:
            trcoord += exons[exind].length
            exind += 1
            continue

        assert trcoord <= interval.start < interval.end <= trcoord + exons[exind].length
        curstart = (interval.start - trcoord) + exons[exind].start
        curend = (interval.end - trcoord) + exons[exind].start
        mapped.append(Interval(exons[exind].chrom, curstart, curend, strand=interval.strand))
    return mapped

# def map_genomic_to_transcriptomic(trindex: str, sites: List[int]) -> List[int]:
#     assert all(before <= after for before, after in zip(sites[:-1], sites[1:]))
#     exons = _trexons()[trindex]
#
#     exind, sind, trcoordinate = 0, 0, 0
#     mapped = []
#     while True:
#         if sind == len(sites) or exind == len(exons):
#             break
#         # map if needed
#         elif exons[exind].start <= sites[sind] < exons[exind].end:
#             mapped.append(trcoordinate + sites[sind] - exons[exind].start)
#             sind += 1
#         elif sites[sind] >= exons[exind].end:
#             trcoordinate += exons[exind].length
#             exind += 1
#         else:
#             assert False, "Unreachable!"
#     assert len(mapped) == len(sites)
#     # Flip if needed
#     if exons[0].strand == "-":
#         trlength = sum(x.length for x in exons)
#         mapped = [(trlength - 1) - x for x in mapped]
#     return mapped
