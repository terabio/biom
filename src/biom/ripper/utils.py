import pickle
from collections import defaultdict
from pathlib import Path
from typing import List

from HTSeq import GenomicArrayOfSets, GenomicInterval

from .cluster import Peak, PeaksCluster


def subtract(clusters: List[PeaksCluster], index: GenomicArrayOfSets) -> List[PeaksCluster]:
    result = []
    for c in clusters:
        key = GenomicInterval(c.contig, c.start, c.end)
        overlapex = False
        for iv, overlaps in index[key].steps():
            if iv.length > 0 and len(overlaps) > 0:
                overlapex = True
                break
        if not overlapex:
            result.append(c)
    return result


def savepe(peaks: List[Peak], file: Path):
    with open(file, 'wb') as stream:
        pickle.dump(peaks, stream, protocol=pickle.DEFAULT_PROTOCOL)


def loadpe(file: Path) -> List[Peak]:
    with open(file, 'rb') as stream:
        peaks = pickle.load(stream)
    return peaks


def savecl(clusters: List[PeaksCluster], file: Path):
    with open(file, 'wb') as stream:
        pickle.dump(clusters, stream, protocol=pickle.DEFAULT_PROTOCOL)


def loadcl(file: Path) -> List[PeaksCluster]:
    with open(file, 'rb') as stream:
        clusters = pickle.load(stream)
    return clusters


def cluster(peaks: List[Peak], maxdist: int) -> List[PeaksCluster]:
    contigs = defaultdict(list)
    for i in peaks:
        contigs[(i.contig, i.strand)].append(i)

    results = []
    for peaks in contigs.values():
        if len(peaks) == 1:
            results.append(PeaksCluster.from_peaks(peaks))
            continue

        peaks = sorted(peaks, key=lambda x: x.start)
        cluster = []
        for p in peaks:
            if len(cluster) == 0 or abs(cluster[-1].end - p.start) < maxdist:
                cluster.append(p)
            else:
                results.append(PeaksCluster.from_peaks(cluster))
                cluster = [p]

        if cluster:
            results.append(PeaksCluster.from_peaks(cluster))
    return results
