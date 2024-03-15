from biom.annotation.resolve import Resolution, AnnotationWeights
from biom.annotation.source import Source
from biom.gindex import Index, overlap


class Annotator[A]:
    def __init__(
            self,
            index: Index[A],
            source: Source,
            resolution: Resolution[list[overlap.Overlap[A]], AnnotationWeights[A]]
    ):
        self.index: Index[A] = index
        self.source: Source = source
        self.resolution: Resolution[list[overlap.Overlap[A]], AnnotationWeights[A]] = resolution
        self.processed: int = 0
        self.results: AnnotationWeights[A] = {}

    def annotate(self, contig: str, start: int, end: int) -> AnnotationWeights[A]:
        for strand, blocks in self.source.fetch(contig, start, end):
            overlaps = [self.index.overlap(contig, strand, rng=rng) for rng in blocks]
            for k, v in self.resolution(overlaps).items():
                self.results[k] = self.results.get(k, 0) + v
        self.processed += 1
        return self.results
