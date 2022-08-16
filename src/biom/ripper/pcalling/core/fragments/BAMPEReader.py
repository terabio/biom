from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

from pysam import AlignmentFile, AlignedSegment


@dataclass
class BundledFragments:
    lmates: List[AlignedSegment]
    rmates: List[AlignedSegment]


class BAMPEReader:
    def __init__(self, filename: Path, contig: str, inflags: int, exflags: int, minmapq: int):
        self.sf = AlignmentFile(filename, 'rb')
        self.contig = contig
        self.inflags = inflags
        self.exflags = exflags
        self.minmapq = minmapq
        self.cache = defaultdict(lambda *args: BundledFragments([], []))

    def makepairs(self, fragments: BundledFragments) -> List[Tuple[AlignedSegment, AlignedSegment]]:
        pairs = []
        missed_lmates = []
        for lmate in fragments.lmates:
            mate_found = False
            for ind, rmate in enumerate(fragments.rmates):
                if lmate.next_reference_id == rmate.reference_id and \
                        lmate.next_reference_start == rmate.reference_start and \
                        rmate.next_reference_id == lmate.reference_id and \
                        rmate.next_reference_start == lmate.reference_start:
                    mate_found = True
                    break

            if not mate_found:
                missed_lmates.append(lmate)
            else:
                fragments.rmates.pop(ind)
                pairs.append((lmate, rmate))
        fragments.lmates = missed_lmates
        return pairs

    def __iter__(self):
        if self.contig not in self.sf.references:
            return

        iterator = self.sf.fetch(self.contig)
        for ind, segment in enumerate(iterator):  # type: (int, AlignedSegment)
            # some required flags were not set OR some excluded flags were set
            if not segment.is_paired or \
                    segment.flag & self.inflags != self.inflags or \
                    segment.flag & self.exflags != 0 or \
                    segment.mapq < self.minmapq:
                continue

            cached = self.cache[segment.query_name]
            if segment.is_read1:
                cached.lmates.append(segment)
            else:
                assert segment.is_read2
                cached.rmates.append(segment)

            if len(cached.lmates) >= 1 and len(cached.rmates) >= 1:
                bundle = self.makepairs(cached)
                if bundle:
                    yield bundle

            if ind % 500_000 == 0:
                for k, v in list(self.cache.items()):
                    if len(v.lmates) == 0 and len(v.rmates) == 0:
                        self.cache.pop(k)

        for cached in self.cache.values():
            bundle = self.makepairs(cached)
            if bundle:
                yield bundle
