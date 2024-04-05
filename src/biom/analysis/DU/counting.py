from collections import defaultdict
from collections.abc import Iterable
from typing import Callable, cast

from attrs import define
from intervaltree import IntervalTree

from biom import algo, io
from biom.analysis import counting
from biom.analysis.nfcore import rnaseq
from biom.ds import gindex
from biom.primitives import Direction, Range, strdeductor, bioproj
from .group import Interval, Group

BamReaderFn = Callable[[rnaseq.Experiment], counting.Source]
StatsFn = Callable[[counting.Source, counting.Counter, counting.AnnotationWeights["AnnoKey | None"]], dict]


@define(slots=True, frozen=True, eq=True, hash=True, order=True)
class AnnoKey:
    """
    A class representing a counting key.

    Attributes:
        group: The identifier of the group.
        region: The identifier of the region.
    """
    group: str
    region: str


@define(slots=True, frozen=True)
class CountingResult:
    """
    A class representing the result of counting reads in a given experiment.

    Attributes:
        project: The identifier of the project.
        experiment: The identifier of the experiment.
        stats: A dictionary containing statistics about the counting process.
        weights: The weights assigned to each annotation category.
    """
    project: str
    experiment: str
    stats: dict
    weights: counting.AnnotationWeights[AnnoKey]


@define(slots=True, frozen=True)
class IntervalsPartition:
    """
    A class representing a partition of intervals.

    Attributes:
        contig: The identifier of the contig.
        strand: The strand of the contig.
        intervals: The intervals in the partition.
        keys: The keys corresponding to the intervals.
    """
    contig: str
    strand: Direction
    intervals: tuple[Interval, ...]
    keys: tuple[AnnoKey, ...]

    def __attrs_post_init__(self):
        if len(self.intervals) != len(self.keys):
            raise ValueError("Intervals and keys must have the same length")

    @staticmethod
    def from_groups(groups: Iterable[Group], maxdist: int = 1024) -> list["IntervalsPartition"]:
        """
        Partition intervals into groups based on contig, strand, and their distance from each other.

        :param groups: An iterable of Group objects.
        :param maxdist: The maximum distance between intervals for them to be considered part of the same group.
        :return: A list of interval partitions.
        """
        # Group by contig and strand
        partition = defaultdict(list)
        for group in groups:
            for region in group.regions:
                for interval in region.intervals:
                    partition[(interval.contig, interval.strand)].append((interval, group.ind, region.ind))

        # Sort by start position and make final partitions
        results = []
        for (contig, strand), part in partition.items():
            part.sort(key=lambda x: x[0].start)
            partitions = algo.misc.group_within(
                part,
                distance=lambda x, y: min(abs(y[0].start - x[0].start), abs(y[0].end - x[0].end)),
                maxdist=maxdist
            )

            for p in partitions:
                intervals, group_inds, region_inds = zip(*p)
                keys = tuple(AnnoKey(gind, rind) for gind, rind in zip(group_inds, region_inds))

                results.append(IntervalsPartition(contig, strand, intervals, keys))
        return results

    def __len__(self) -> int:
        return len(self.intervals)


def default_resolver(groups: Iterable[Group]) -> counting.Resolution[
    list[counting.Overlap[AnnoKey]], counting.AnnotationWeights[AnnoKey | None]
]:
    """
    Returns a function that resolves annotation weights.

    :param groups: A list of groups for the annotation.
    :return: A function that takes a dict of annotation weights and returns a dictionary of resolved annotation weights.
    """
    group2refind: dict[str, str] = {g.ind: g.reference.ind for g in groups}

    def fn(
            weights: counting.AnnotationWeights[AnnoKey | None]
    ) -> counting.AnnotationWeights[AnnoKey | None]:
        # If there is only one category and its None -> return it as a final result
        if len(weights) == 1 and None in weights:
            return weights
        # Otherwise, remove None category completely
        elif None in weights:
            weights.pop(None)

        # Please the type checker
        pure_weights = cast(counting.AnnotationWeights[AnnoKey], weights)

        # Partition hits by group ids
        groups: defaultdict[str, set[str]] = defaultdict(set)
        for key in pure_weights:
            groups[key.group].add(key.region)

        # Cleanup the annotation within groups
        for gind, rinds in groups.items():
            # Ignore reference region if a read also belongs to at least one region of interest
            refind = group2refind[gind]
            if len(rinds) > 1 and refind in rinds:
                weights.pop(AnnoKey(gind, refind))

        return weights

    resolution = counting.resolve.Binary[AnnoKey]().apply(fn)
    return resolution


def default_reader(experiment: rnaseq.Experiment) -> counting.Source:
    """
    Returns a source object for reading from a BAM file.

    :param experiment: An experiment object.
    :return: A reads source object for reading from a BAM file.
    """
    stranding = strdeductor.from_bioproj(experiment.source)

    runs = list(experiment.runs)
    layout = runs[0].layout
    for r in runs[1:]:
        if r.layout != layout:
            raise ValueError(
                f"All runs must have the same layout. Found different layouts: {layout} and {r.layout} in {experiment}"
            )

    match layout:
        case bioproj.SeqLayout.Single:
            reader = io.bam.Reader(experiment.bam, inflags=0, exflags=2564, minmapq=0, statistics=True)
            return counting.source.SingleEndBAM(reader, stranding)
        case bioproj.SeqLayout.Paired:
            reader = io.bam.Reader(experiment.bam, inflags=3, exflags=2564, minmapq=0, statistics=True)
            bundler = io.bam.PEReadsBundler(reader)
            return counting.source.PairedEndBAM(bundler, stranding)
        case _:
            raise ValueError(f"Unsupported layout: {layout}")


def default_stats(
        source: counting.Source,
        counter: counting.Counter,
        weights: counting.AnnotationWeights[AnnoKey | None]
) -> dict:
    """
    Returns statistics about the counting process.

    :param source: A source object for reading from a BAM file.
    :param counter: A counter object for counting reads.
    :param weights: The weights assigned to each annotation category.
    :return: A dictionary containing statistics about the counting process.
    """
    outside_annotation = weights.pop(None, 0)
    total_counts = sum(weights.values())

    stats = {
        "Total counts": total_counts,
        "Fragments beyond annotation": outside_annotation,
        "Fragment inside annotation": counter.processed - outside_annotation,
    }

    if isinstance(source, counting.source.PairedEndBAM):
        bampe = source.obj
        if bampe.reader.statistics:
            assert bampe.reader.discarded is not None
            assert bampe.reader.consumed is not None

            stats.update({
                "Discarded reads": bampe.reader.discarded.total,
                "Processed reads": bampe.reader.consumed.total,
                "Unpaired reads": bampe.unpaired,
            })
    elif isinstance(source, counting.source.SingleEndBAM):
        bamse = source.obj
        if bamse.statistics:
            assert bamse.discarded is not None
            assert bamse.consumed is not None
            stats.update({
                "Discarded reads": bamse.discarded.total,
                "Processed reads": bamse.consumed.total,
            })
    else:
        raise ValueError(f"Unsupported source type: {source}")
    return stats


def run(
        projects: Iterable[rnaseq.Project],
        partition: IntervalsPartition,
        resolution: counting.Resolution,
        bam_reader_fn: BamReaderFn = default_reader,
        stats_fn: StatsFn = default_stats
) -> Iterable[CountingResult]:
    if len(partition) == 0:
        raise ValueError("Partition must have at least one interval")

    # Build the index
    skeleton = {(partition.contig, Direction.fwd): IntervalTree(), (partition.contig, Direction.rev): IntervalTree()}
    assert all(i.contig == partition.contig for i in partition.intervals)
    assert all(i.strand == partition.strand for i in partition.intervals)

    start, end = partition.intervals[0].start, partition.intervals[0].end
    for interval, key in zip(partition.intervals, partition.keys):
        skeleton[(partition.contig, partition.strand)].addi(interval.start, interval.end, data=key)
        start = min(start, interval.start)
        end = max(end, interval.end)
    index: gindex.GenomicIndex[AnnoKey] = gindex.GenomicIndex(skeleton)
    index_range = Range(start, end)

    # Count reads
    results = []
    for project in projects:
        for experiment in project.experiments:
            source = bam_reader_fn(experiment)
            annotator = counting.Counter(index, source, resolution)
            weights = annotator.annotate(partition.contig, index_range.start, index_range.end)

            stats = stats_fn(source, annotator, weights) if stats_fn is not None else {}

            # Remove None category if it is present -> we are not interested in reads outside the annotation
            if None in weights:
                weights.pop(None)

            results.append(CountingResult(project.ind, experiment.ind, stats, weights))

    return results
