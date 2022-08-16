from dataclasses import dataclass
from pathlib import Path
from typing import Callable, List, Optional, Tuple

import numpy as np
import numpy.typing as npt
from pysam import AlignedSegment

from .BAMPEReader import BAMPEReader
from ..utils import Stranded
from .strdeductors import StrandDeductor


@dataclass(frozen=True)
class AlignedBlocks:
    # Transcription strand (+, -, .)
    trstrand: str
    # Start for each aligned block
    start: npt.NDArray[np.int32]
    # End for each aligned block
    end: npt.NDArray[np.int32]
    # Start index for each fragment.
    # I.e. start[index[i]: index[i + 1]] - starts for aligned block in the i-th read.
    records: npt.NDArray[np.int32]

    @staticmethod
    def from_tuples(trstrand: str, reads: List[List[Tuple[int, int]]]) -> 'AlignedBlocks':
        assert trstrand == "+" or trstrand == "-"
        starts, ends, index, i = [], [], [], 0
        for read in reads:
            index.append(i)
            for start, end in read:
                starts.append(start)
                ends.append(end)
                i += 1
        index.append(i)
        return AlignedBlocks(
            trstrand,
            np.asarray(starts, dtype=np.int32),
            np.asarray(ends, dtype=np.int32),
            np.asarray(index, dtype=np.int32)
        )

    def fragments(self):
        return self.records.size


@dataclass()
class AlignedBlocksBuilder:
    trstrand: str
    # State
    blockind: int
    start: List[int]
    end: List[int]
    index: List[int]

    def __init__(self, trstrand: str):
        self.trstrand = trstrand
        self.blockind = 0
        self.start, self.end, self.index = [], [], []

    def add(self, lmate, rmate):
        self.index.append(self.blockind)

        blocks = []
        for segment in lmate, rmate:
            # detect aligned blocks by parsing CIGAR string
            refpos = segment.reference_start
            for op, oplen in segment.cigartuples:
                # Consume reference: N, D
                if op == 2 or op == 3:
                    refpos += oplen
                # Matched: M, X, =
                elif op == 0 or op == 7 or op == 8:
                    blocks.append((refpos, refpos + oplen))
                    refpos += oplen
                # else => consume query or do nothing

        blocks = sorted(blocks)
        assert len(blocks) >= 1
        # Remove redundant blocks and store results
        curstart, curend = blocks[0]
        for blstart, blend in blocks[1:]:
            if curstart <= blstart <= curend:
                curend = max(blend, curend)
            else:
                self.start.append(curstart)
                self.end.append(curend)
                self.blockind += 1
                curstart, curend = blstart, blend

        # Store the final block
        self.start.append(curstart)
        self.end.append(curend)
        self.blockind += 1

    def finalize(self) -> Optional[AlignedBlocks]:
        # Nothing to report
        if self.blockind == 0:
            return None
        # Final index
        self.index.append(self.blockind)

        # Sort blocks by a start position and convert them to np arrays
        argsort = sorted(range(len(self.index) - 1), key=lambda ind: self.start[self.index[ind]])
        npstart, npend, npind = [np.empty(x, dtype=np.int32) for x in (len(self.start), len(self.end), len(self.index))]

        blockind, trind = 0, 0
        for ind in argsort:
            assert ind < len(self.index)
            npind[trind] = blockind

            _start, _end = self.index[ind], self.index[ind + 1]
            blstart, blend = self.start[_start: _end], self.end[_start: _end]
            size = _end - _start
            npstart[blockind: blockind + size] = blstart
            npend[blockind: blockind + size] = blend

            blockind += size
            trind += 1
        npind[trind] = blockind
        assert npind[-1] == npstart.size == npend.size and trind + 1 == npind.size
        return AlignedBlocks(self.trstrand, npstart, npend, npind)


def _oncontig(
        reader: BAMPEReader, strdeductor: StrandDeductor
) -> Tuple[Optional[AlignedBlocks], Optional[AlignedBlocks]]:
    forward, reverse = AlignedBlocksBuilder("+"), AlignedBlocksBuilder("-")

    for bundle in reader:
        for lmate, rmate in bundle:  # type: (AlignedSegment, AlignedSegment)
            strand = strdeductor(lmate, rmate)
            assert strand == "+" or strand == "-"
            if strand == "+":
                forward.add(lmate, rmate)
            else:
                reverse.add(lmate, rmate)

    return forward.finalize(), reverse.finalize()


def loadfrom(
        files: List[Path], strdeductor: StrandDeductor, contig: str, inflags: int,
        exflags: int, minmapq: int
) -> Tuple[Stranded[List[AlignedBlocks]], int]:
    assert files

    forward, reverse = [], []
    contiglens = []
    for file in files:
        reader = BAMPEReader(file, contig, inflags, exflags, minmapq)
        fwdblocks, revblocks = _oncontig(reader, strdeductor)
        if fwdblocks:
            forward.append(fwdblocks)
        if revblocks:
            reverse.append(revblocks)

        try:
            contiglens.append(reader.sf.get_reference_length(contig))
        except KeyError:
            continue

    assert contiglens, files
    contiglen = contiglens[0]
    assert all(contiglen == x for x in contiglens), f"Contradictory contig length for '{contig}': {contiglens}"
    return Stranded(forward, reverse), contiglen
