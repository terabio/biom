from biom.annotation.resolution import Resolution, AnnotationWeights
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
            strategy: Resolution[list[overlap.Overlap[A]], AnnotationWeights[A]]
    ):
        self.index = index
        self.source = source
        self.strategy = strategy

    def annotate(self, contig: str, strand: str, start: int, end: int) -> AnnotationWeights[A]:
        # results = {}
        # for blocks in self.source.fetch(contig, start, end):
        #     pass
        raise NotImplementedError("Not implemented yet")

    # def annotate(self, contig: str, strand: str, blocks: Iterable[Range]) -> AnnotationWeights:
    #     overlaps = [self.index.overlap(contig, strand, bl.start, bl.end) for bl in blocks]
    #     return self.strategy(overlaps)
