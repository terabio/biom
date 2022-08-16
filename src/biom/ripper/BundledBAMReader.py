import warnings
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Tuple

from pysam import AlignmentFile, AlignedSegment


@dataclass
class _BundledFragments:
    lmates: List[AlignedSegment]
    rmates: List[AlignedSegment]


class BundledBAMReader:
    def __init__(self, filename: str, threads: int = 1):
        self.sf = AlignmentFile(filename, 'r', threads=threads)
        self.unpaired = 0
        self.singletons: List[AlignedSegment] = []
        self.cache = defaultdict(lambda *args: _BundledFragments([], []))

    def makepairs(self, fragments: _BundledFragments) -> List[Tuple[AlignedSegment, AlignedSegment]]:
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
        for ind, segment in enumerate(self.sf):
            if not segment.is_paired:
                self.unpaired += 1
                continue
            if segment.is_paired:
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
            else:
                self.singletons.extend(cached.lmates)
                self.singletons.extend(cached.rmates)

        if self.unpaired > 1:
            warnings.warn(f"{self.unpaired} reads with missing mates encountered")
