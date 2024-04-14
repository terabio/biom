from abc import ABC, abstractmethod
from typing import TypeVar, Generic

from biom.ds.gindex import Counts

_T = TypeVar("_T")


class ReadCounter(ABC, Generic[_T]):
    """
    An abstract base class representing a read counter.
    This class defines the interface for a read counter, which processes reads and accumulates counts.
    """

    @abstractmethod
    def process_all(self) -> Counts[_T]:
        """
        Process all reads in the source and return the counts.

        :return: The counts of all reads in the source.
        """
        pass

    @abstractmethod
    def stats(self) -> dict:
        """
        Return statistics about the counting process. The nature of the statistics is implementation-specific.

        :return: A dictionary of statistics about the counting process.
        """
        pass

    @abstractmethod
    def counts(self) -> Counts[_T]:
        """
        Return the counts accumulated during the counting process.

        :return: The counts accumulated during the counting process.
        """
        pass

    @abstractmethod
    def reset(self):
        """
        Reset the counter, removing all accumulated counts.
        """
        pass


class IndexedReadCounter(ReadCounter[_T]):
    """
    A class representing a read counter that can process reads within a specific region.
    This class extends the ReadCounter class with the ability to process reads within a specific region.
    """

    @abstractmethod
    def process_within(self, contig: str, start: int, end: int) -> Counts[_T]:
        """
        Process reads within a specific region and return the counts.

        :param contig: The contig of the region.
        :param start: The start position of the region.
        :param end: The end position of the region.
        :return: The counts of reads within the specified region.
        """
        pass
