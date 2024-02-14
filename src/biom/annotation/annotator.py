from biom.annotation.resolve import Resolution, AnnotationWeights
from biom.annotation.source import Source
from biom.gindex import Index, overlap


class Annotator[A]:
    index: Index[A]
    source: Source
    resolver: Resolution[list[overlap.Overlap[A]], AnnotationWeights[A]]

    def __init__(
            self,
            index: Index[A],
            source: Source,
            resolution: Resolution[list[overlap.Overlap[A]], AnnotationWeights[A]]
    ):
        self.index = index
        self.source = source
        self.resolution = resolution

    def annotate(self, contig: str, start: int, end: int) -> AnnotationWeights[A]:
        results = {}
        for strand, blocks in self.source.fetch(contig, start, end):
            overlaps = [self.index.overlap(contig, strand, rng=rng) for rng in blocks]
            for k, v in self.resolution(overlaps).items():
                results[k] = results.get(k, 0) + v
        return results
