import copy
import time
from collections import defaultdict
from collections.abc import Iterable
from typing import TypeVar

from attrs import define
from attrs import field
from intervaltree import IntervalTree
from joblib import Parallel, delayed

from biom.ds.gindex import Overlap, GenomicIndex
from biom.primitives import Interval, Orientation, Range
from .partition import Partition
from .reads_counter import MultiReadsCounter
from ..resolve import Counts, Resolution
from ..source import Source

_T = TypeVar('_T')


def run(
        tag: str,
        source: Source,
        partition: Partition[_T],
        resolution: Resolution[list[Overlap[_T]], Counts[_T]],
) -> tuple[str, Counts[_T], dict]:
    if len(partition) == 0:
        raise ValueError("Partition must have at least one interval")

    launched_at = time.time()

    # Build the index
    skeleton = {
        (partition.contig, Orientation.fwd): IntervalTree(),
        (partition.contig, Orientation.rev): IntervalTree(),
        (partition.contig, Orientation.dual): IntervalTree(),
    }
    assert all(i.contig == partition.contig for i in partition.intervals)

    start, end = partition.intervals[0].rng.start, partition.intervals[0].rng.end
    for interval, data in zip(partition.intervals, partition.data):
        skeleton[(partition.contig, interval.orient)].addi(
            interval.rng.start, interval.rng.end, data=data
        )

        start = min(start, interval.rng.start)
        end = max(end, interval.rng.end)

    index: GenomicIndex[_T] = GenomicIndex(skeleton)
    index_range = Range(start, end)

    # Count reads
    counts: defaultdict[_T, float] = defaultdict(float)
    for blocks in source.fetch(partition.contig, index_range.start, index_range.end):
        overlaps = [index.overlap(partition.contig, blocks.orientation, rng=rng) for rng in blocks.blocks]
        for k, v in resolution(overlaps).items():
            counts[k] += v

    finished_at = time.time()

    stats = source.stats()

    stats['Time(s)'] = finished_at - launched_at
    stats['Partition'] = f"{partition.contig}:{start}-{end}"

    return tag, dict(counts), source.stats()


@define(slots=True)
class JoblibMultiReadsCounter(MultiReadsCounter[_T]):
    """
    A class representing a reads counter that uses joblib to parallelize the counting process.
    """
    _sources: dict[str, Source] = field(alias="sources")
    resolution: Resolution[list[Overlap[_T]], Counts[_T]]
    parallel: Parallel
    _counts: dict[str, Counts[_T]] = field(factory=dict, init=False)
    _stats: dict = field(factory=dict, init=False)

    def count(self, data: Iterable[_T], intervals: Iterable[Interval]):
        """
        Count the reads in the given intervals.

        :param data: The data to count.
        :param intervals: Intervals associated with the data.
        """
        _data, _intervals = list(data), list(intervals)
        if len(_data) != len(_intervals):
            raise ValueError(
                f"Data and intervals must have the same length, got {len(_data)} and {len(_intervals)} respectively."
            )

        # Group data into local partitions to support the parallel counting
        partitions = Partition.from_intervals(_data, _intervals)

        # Prepare the workload and run the counting process
        workloads = []
        for partition in partitions:
            for tag, source in self._sources.items():
                workloads.append((tag, copy.deepcopy(source), partition, self.resolution))

        results = self.parallel(delayed(run)(*workload) for workload in workloads)

        # Merge counts & stats
        for tag, counts, stats in results:
            _counts = self._counts.get(tag, {})
            for k, v in counts.items():
                _counts[k] = _counts.get(k, 0) + v

            _stats = self._stats.get(tag, {})
            for k, v in stats.items():
                _stats[k] = _stats.get(k, 0) + v

    def counts(self) -> dict[str, Counts[_T]]:
        return self._counts

    def sources(self) -> dict[str, Source]:
        return self._sources

    def stats(self) -> dict:
        return self._stats

    def reset(self):
        self._counts.clear()
        self._stats.clear()
