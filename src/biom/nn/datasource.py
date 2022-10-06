from abc import ABC, abstractmethod
from collections import defaultdict
from pathlib import Path

import numpy as np
import pyBigWig
from pybedtools import BedTool, Interval


class DataSource(ABC):
    @abstractmethod
    def fetch(self, interval: Interval) -> np.ndarray:
        pass


class BED(DataSource):
    def __init__(self, bed: Path, presorted_and_merged: bool = False):
        bed = BedTool(bed)
        if not presorted_and_merged:
            bed = bed.sort().merge().sort()
        # Keep bed file for logging purposes
        self.bedfn = bed.saveas()

        # Create index
        index = defaultdict(list)
        for inter in bed:
            index[inter.chrom].append(inter)
        self.index = {
            chrom: (np.asarray([x.start for x in intervals]), np.asarray([x.end for x in intervals]))
            for chrom, intervals in index.items()
        }

    def fetch(self, interval: Interval) -> np.ndarray:
        result = np.zeros(interval.length, dtype=np.float32)

        if interval.chrom not in self.index:
            return result

        starts, ends = self.index[interval.chrom]
        ind = np.searchsorted(starts, interval.start)
        if ind > 0:
            # check to avoid missing hits with the end of the previous interval
            interstrt, interend = max(starts[ind - 1], interval.start), min(ends[ind - 1], interval.end)
            if interstrt < interend:
                result[interstrt - interval.start: interend - interval.start] = 1

        while ind < len(starts):
            interstrt, interend = max(starts[ind], interval.start), min(ends[ind], interval.end)
            if interstrt >= interend:
                break
            result[interstrt - interval.start:  interend - interval.start] = 1
            ind += 1
        return result


class BigWig(DataSource):
    def __init__(self, path: Path):
        self.path = path
        self.bw = None

    def fetch(self, interval: Interval) -> np.ndarray:
        if self.bw is None:
            self.bw = pyBigWig.open(self.path.as_posix())
        values = self.bw.values(interval.chrom, interval.start, interval.end, numpy=True)
        values[np.isnan(values)] = 0
        return values


# class NPArray(DataSource):
#     def __init__(self, folder: Path):
#         self.files = {}
#         for contig in folder.iterdir():
#             self.files[f"{contig.name}.np"] = contig.as_posix()
#         self.ndarrays = {}
#
#     @staticmethod
#     def cache(bigwig: Path, saveto: Path):
#         bigwig: pyBigWig.pyBigWig = pyBigWig.open(bigwig)
#
#         for contig, length in bigwig.chroms():
#             file = saveto.joinpath(f"{contig}.np")
#             values = bigwig.values(contig, 0, length, numpy=True)
#             values[np.isnan(values)] = 0
#             np.save(file.as_posix(), values)
#
#     def fetch(self, interval: Interval) -> np.ndarray:
#         assert interval.chrom in self.files, "No data available for the give contig"
#         file = self.files[interval.chrom]
#         if file not in self.ndarrays:
#             self.ndarrays[file] = np.memmap(file, dtype=np.int32, mode='r')
#
#         array = self.ndarrays[file][interval.start: interval.end]
#         return np.asarray(array)
