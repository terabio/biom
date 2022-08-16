from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

from pybedtools import Interval

from ..coordmap import genome2transcriptome, transcriptome2genome


@dataclass()
class Peak:
    contig: str
    start: int
    end: int
    strand: str
    foldchange: float
    CPM: float

    repeats: Dict[str, float]
    transcripts: Dict[str, float]
    locations: Dict[str, float]

    @property
    def length(self):
        return self.end - self.start

    def __repr__(self):
        return f"{self.contig}:{self.start}-{self.end}({self.strand})"

    def __hash__(self):
        return hash(self.contig) + hash(self.start) + hash(self.end) + hash(self.strand)

    def __eq__(self, other):
        return isinstance(other, Peak) and other.contig == self.contig and other.start == self.start and \
               other.end == self.end and other.strand == self.strand

    def to_transcript(self, exons: List[Interval]) -> 'Peak':
        mapped = genome2transcriptome(exons, [Interval(self.contig, self.start, self.end, strand=self.strand)])
        return Peak(mapped[0].chrom, mapped[0].start, mapped[0].end, self.strand, self.foldchange, self.CPM,
                    self.repeats, self.transcripts, self.locations)

    def to_genome(self, exons: List[Interval]) -> 'Peak':
        mapped = transcriptome2genome(exons, [Interval(self.contig, self.start, self.end, strand=self.strand)])
        return Peak(mapped[0].chrom, mapped[0].start, mapped[0].end, self.strand, self.foldchange, self.CPM,
                    self.repeats, self.transcripts, self.locations)


@dataclass()
class PeaksCluster:
    contig: str
    start: int
    end: int
    strand: str

    foldchange: float
    CPM: float
    peaks: Tuple[Peak, ...]

    label: Optional[str] = None

    @property
    def length(self):
        return self.end - self.start

    def transcripts(self) -> Dict[str, Set[str]]:
        transcripts = defaultdict(set)
        for p in self.peaks:
            for k, v in p.transcripts.items():
                transcripts[k] |= v
        return transcripts

    def locations(self) -> Dict[str, float]:
        location = defaultdict(float)
        for p in self.peaks:
            for k, v in p.locations.items():
                location[k] += v * p.length
        total = sum(location.values())
        return {k: v / total for k, v in location.items()}

    def repeats(self) -> Dict[str, float]:
        repeats = defaultdict(float)
        for p in self.peaks:
            for k, v in p.repeats.items():
                repeats[k] += v * p.length
        total = sum(repeats.values())
        return {k: v / total for k, v in repeats.items()}

    @staticmethod
    def from_peaks(peaks: List[Peak]) -> 'PeaksCluster':
        assert len(peaks) >= 1
        assert all(x.contig == peaks[0].contig for x in peaks)
        assert all(x.strand == peaks[0].strand for x in peaks)
        peaks = sorted(peaks, key=lambda x: x.start)
        start, end = peaks[0].start, max(x.end for x in peaks)
        top = max(peaks, key=lambda x: x.foldchange)
        return PeaksCluster(top.contig, start, end, top.strand, top.foldchange, top.CPM, tuple(peaks))

    def __repr__(self):
        return f"{self.contig}:{self.start}-{self.end}({self.strand})"
