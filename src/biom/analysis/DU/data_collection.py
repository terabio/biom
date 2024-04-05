from collections.abc import Iterable

import numpy as np
import numpy.typing as npt
from attrs import define, field
from joblib import Parallel, delayed

from biom.analysis.counting import Resolution
from biom.analysis.nfcore import rnaseq
from . import counting
from .group import Group

# (Project ind -> Experiment ind -> Group ind -> Region ind) => Counts
CountsCollection = dict[
    tuple[str, str, str, str],
    int
]

# (Project ind -> Experiment ind -> Group ind -> Region ind) => list[npt.NDArray[np.float32]]
# Order of the list is the same as the order of the intervals in the region.
CoverageCollection = dict[
    tuple[str, str, str, str],
    list[npt.NDArray[np.float32]]
]


@define(slots=True)
class DataCollection:
    """
    A class representing a collection of data for differential usage analysis.

    This class holds the groups for differential usage analysis, each of which includes a reference region and
    regions of interest. It also holds the bio-projects containing the experiments to be analyzed. The counts
    and coverage attributes are pandas DataFrames that hold the counts of reads mapped to each region in each
    group in each experiment, and the coverage of reads mapped to each interval in each region in each group
    in each experiment, respectively.

    :param groups: The groups for differential usage analysis.
    :param projects: The nf-core/rnaseq containing the experiments to be analyzed.
    :param counts: The total number of reads mapped to each region in each group in each experiment.
    :param coverage: The coverage of each interval in each region in each group in each experiment.
    """
    groups: tuple[Group, ...]
    projects: tuple[rnaseq.Project, ...]
    counts: CountsCollection = field(factory=dict, init=False)
    coverage: CoverageCollection = field(factory=dict, init=False)

    def count_reads(
            self, *,
            pool: Parallel | None = None,
            cache: Parallel | None = None,
            resolution: Resolution | None = None,
            reader: counting.BamReaderFn = counting.default_reader,
            stats: counting.StatsFn = counting.default_stats
    ) -> Iterable[dict]:
        """
        Count the number of reads mapped to each region in each group in each experiment.

        :param pool: The joblib.Parallel pool to use for parallelization. If None, the function will run in serial mode.
        :param cache: The joblib.Memory cache to use for memoization. If None, no memoization will be used.
        :param resolution: Resolution is mapping between overlaps of a read and the intervals to final counts (see biom.analysis.counting).
        :param reader: A function that takes an experiment and returns a `biom.analysis.counting.Source`.
        :param stats: A function that runs after the counting is done. Should take a source, a counter, and a dictionary
            of annotation weights and returns a dictionary containing arbitrary statistics of the current run.
        :return: Collection of statistics for each counting run.
        """

        # Partition intervals into chunks to efficiently parallelize counting
        partitions = counting.IntervalsPartition.from_groups(self.groups)

        # Resolution is mapping between overlaps of a read and the intervals to final counts
        if resolution is None:
            resolution = counting.default_resolver(self.groups)

        # Cache is used to memoize the results of counting and avoid recomputing them on subsequent runs
        if cache is None:
            job = counting.run
        else:
            job = cache(counting.run)

        # Do the counting
        results: list[Iterable[counting.CountingResult]]
        if pool is None:
            results = [
                job(self.projects, partition, resolution, reader, stats) for partition in partitions
            ]
        else:
            results = pool(
                delayed(job)(self.projects, partition, resolution, reader, stats) for partition in partitions
            )

        # Flatten the results
        counts, statistics = {}, []
        for chunk in results:
            for result in chunk:
                statistics.append(result.stats)
                for key, weight in result.weights.items():
                    newkey = (result.project, result.experiment, key.group, key.region)
                    assert newkey not in counts
                    counts[newkey] = weight
        return statistics

    def fetch_coverage(self, *, pool: Parallel | None = None, cache: Parallel | None = None) -> Parallel:
        """
        Count the number of reads mapped to each region in each group in each experiment.

        :param pool: The joblib.Parallel pool to use for parallelization. If None, the function will run in serial mode.
        :param cache: The joblib.Memory cache to use for memoization. If None, no memoization will be used.
        """

        pass
