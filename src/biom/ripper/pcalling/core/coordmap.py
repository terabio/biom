import unittest
from copy import deepcopy
from typing import Optional, List, Tuple

import numba
import numpy as np

from .fragments import AlignedBlocks
from .pileup import Pileup
from .pileup.merge import _simplify
from .transcriptome import TranscriptInfo


@numba.jit(cache=True, nopython=True, nogil=True)
def _to_transcript(excoord, tr_exstart, tr_exend, block_starts, block_ends, block_records, hint):
    # Pre-allocated memory for mapped coordinates
    maxlen = block_starts.size
    mapped_starts, mapped_ends = np.empty(maxlen, dtype=np.int32), np.empty(maxlen, dtype=np.int32)
    mapped_records = np.empty(block_records.size, dtype=np.int32)

    # Pointers for the pre-allocated array
    pblock, precord = 0, 0

    # Iterate over all requested records
    for ind in range(hint, block_records.size - 1):
        # Select matched start/end coordintates
        _start, _end = block_records[ind], block_records[ind + 1]
        starts, ends = block_starts[_start: _end], block_ends[_start: _end]

        # Before the transcript
        if ends[-1] <= tr_exstart[0]:
            continue
        # After the transcript => break (blocks are sorted)
        if starts[0] >= tr_exend[-1]:
            break

        # Store the index pointer
        mapped_records[precord] = pblock
        pblock_before_mapping = pblock

        # Exon coordinates
        exind, blind = 0, 0
        exstart, exend = tr_exstart[exind], tr_exend[exind]
        # Current block
        blstart, blend = starts[blind], ends[blind]
        # Lattest mapped position in the transcriptome coordinates
        prevend = None
        while True:
            # Block to the left of the current exon -> pick next block
            if blend <= exstart:
                blind += 1
                if blind == starts.size:
                    break
                blstart, blend = starts[blind], ends[blind]
            # Exon to the left of the current block -> pick next exon
            elif exend <= blstart:
                exind += 1
                if exind == tr_exstart.size:
                    break
                exstart, exend = tr_exstart[exind], tr_exend[exind]
            # There is an intersection
            else:
                interst, internd = max(blstart, exstart), min(blend, exend)
                assert internd - interst > 0
                # Transform to the transcriptom coordinates
                interst = interst - exstart + excoord[exind]
                internd = internd - exstart + excoord[exind]

                # The first mapped block
                if prevend is None:
                    mapped_starts[pblock] = interst
                    prevend = internd
                # Current block mapped just after the previous one
                elif prevend == interst:
                    prevend = internd
                # Current block mapped with a gap relative to the previous one
                else:
                    mapped_ends[pblock] = prevend  # Finish the previous record
                    pblock += 1

                    mapped_starts[pblock] = interst  # Start a new record
                    prevend = internd

                # What to consider next, exon or a block?
                if blend < exend:
                    blind += 1
                    if blind == starts.size:
                        break
                    blstart, blend = starts[blind], ends[blind]
                else:
                    assert blend >= exend

                    exind += 1
                    if exind == tr_exstart.size:
                        break
                    exstart, exend = tr_exstart[exind], tr_exend[exind]

        # Finish the latest record
        if prevend is not None:
            mapped_ends[pblock] = prevend
            pblock += 1

        # If we mapped at least something -> advance the records index
        if pblock != pblock_before_mapping:
            precord += 1

    # Store the end index for the latest record
    mapped_records[precord] = pblock

    # Chop off unused memory
    mapped_records = mapped_records[: precord + 1]
    mapped_starts, mapped_ends = mapped_starts[:pblock], mapped_ends[:pblock]
    return mapped_records, mapped_starts, mapped_ends


def to_transcript(tr: TranscriptInfo, blocks: AlignedBlocks, hint: int = 0) -> Optional[AlignedBlocks]:
    assert 0 <= hint < blocks.records.size
    # Calculate coordinates for each exon in the transcript
    excoord = np.empty(tr.exstart.size, dtype=np.int32)
    excoord[0] = 0
    np.cumsum(tr.exend[:-1] - tr.exstart[:-1], out=excoord[1:])

    mapped_records, mapped_starts, mapped_ends = _to_transcript(
        excoord, tr.exstart, tr.exend, blocks.start, blocks.end, blocks.records, hint
    )
    if mapped_starts.size == 0:
        return None
    assert mapped_records[-1] == len(mapped_starts) and mapped_starts.size == mapped_ends.size

    return AlignedBlocks(
        trstrand=tr.trstrand,
        start=mapped_starts,
        end=mapped_ends,
        records=mapped_records
    )


@numba.jit(cache=True, nopython=True, nogil=True)
def _to_genome(contiglen, tr_exstart, tr_exend, pileup_values, pileup_interend):
    maxbreaks = (pileup_values.size + 1) + (2 * tr_exstart.size)

    mappedend, mappedval = np.empty(maxbreaks, dtype=np.int32), np.empty(maxbreaks, dtype=np.float32)

    # Before the transcript
    pos = 0
    if tr_exstart[0] != 0 and pileup_values[0] != 0:
        mappedend[pos] = tr_exstart[0]
        mappedval[pos] = 0
        pos += 1

    # Transcript coordinates
    trpos, gnpos, curexon = 0, tr_exstart[0], 0
    for i in range(pileup_values.size):
        trend, trval = pileup_interend[i], pileup_values[i]
        assert trend > trpos

        # Move in the transcriptom coordinates until the given interval is fully consumed
        while trpos != trend:
            left = trend - trpos
            # Completely inside the current exon
            if gnpos + left < tr_exend[curexon]:
                mappedend[pos] = gnpos + left
                mappedval[pos] = trval
                pos += 1

                gnpos += left
                trpos += left
            # Interval ends on the corner of the current exon
            elif gnpos + left == tr_exend[curexon]:
                mappedend[pos] = gnpos + left
                mappedval[pos] = trval
                pos += 1

                gnpos += left
                trpos += left

                # add intron if it is not the last exon
                is_last_exon = curexon == tr_exstart.size - 1
                if not is_last_exon:
                    curexon += 1
                    assert tr_exstart[curexon] > gnpos
                    gnpos = tr_exstart[curexon]

                    mappedend[pos] = gnpos
                    mappedval[pos] = 0
                    pos += 1
                else:
                    is_last_interval = i == pileup_values.size - 1
                    assert is_last_interval
            # Interval is interrupted by the intron
            else:
                is_last_exon = curexon == tr_exstart.size - 1
                assert not is_last_exon and gnpos < tr_exend[curexon]
                # consume exon before intron
                tail = tr_exend[curexon] - gnpos
                mappedend[pos] = tr_exend[curexon]
                mappedval[pos] = trval
                pos += 1

                # count consumed interval
                trpos += tail

                # inton
                curexon += 1
                mappedend[pos] = tr_exstart[curexon]
                mappedval[pos] = 0
                pos += 1

                # count consumed intron
                gnpos = tr_exstart[curexon]
    # All transcript must have been processed
    assert trpos == (tr_exend - tr_exstart).sum()

    # After the transcript
    if gnpos != contiglen:
        if mappedval[pos - 1] == 0:
            mappedend[pos - 1] = contiglen
        else:
            assert pos < maxbreaks
            mappedend[pos] = contiglen
            mappedval[pos] = 0
            pos += 1
    return mappedend[:pos], mappedval[:pos]


def to_genome(contig: str, contiglen: np.int32, tr: TranscriptInfo, pileup: Pileup) -> Pileup:
    mappedend, mappedval = _to_genome(contiglen, tr.exstart, tr.exend, pileup.values, pileup.interend)
    mappedend, mappedval = _simplify(mappedend, mappedval)
    return Pileup(contig, mappedend, mappedval)


class PileupUnitTests(unittest.TestCase):
    def _test_to_genome(self, contiglen: int, exons: List[Tuple[int, int]],
                        pileup: List[Tuple[int, int]], expected: List[Tuple[int, int]]):
        exstart, exend = zip(*exons)
        trinfo = TranscriptInfo(
            ".", contiglen, np.asarray(exstart, dtype=np.int32), np.asarray(exend, dtype=np.int32), "TMP"
        )

        pileup = Pileup.from_tuples(".", pileup)

        beforetr = deepcopy(trinfo)

        result = to_genome(trinfo, pileup)

        expected = Pileup.from_tuples(".", expected)
        np.testing.assert_array_equal(result.interend, expected.interend)
        np.testing.assert_array_equal(result.values, expected.values)

        # Make sure that arguments were not corrupted
        np.testing.assert_array_equal(trinfo.exstart, beforetr.exstart)
        np.testing.assert_array_equal(trinfo.exend, beforetr.exend)

    def test_to_genome(self):
        # setup 1
        contiglen = 18
        exons = [(1, 4), (5, 9), (11, 15), (16, 18)]

        trpileip = [(13, 1)]
        expected = [(1, 0), (4, 1), (5, 0), (9, 1), (11, 0), (15, 1), (16, 0), (18, 1)]
        self._test_to_genome(contiglen, exons, trpileip, expected)

        trpileip = [(1, 0), (2, 1), (8, 3), (13, 9)]
        expected = [(2, 0), (3, 1), (4, 3), (5, 0), (9, 3), (11, 0), (12, 3), (15, 9), (16, 0), (18, 9)]
        self._test_to_genome(contiglen, exons, trpileip, expected)

        trpileip = [(3, 1), (7, 2), (11, 3), (13, 4)]
        expected = [(1, 0), (4, 1), (5, 0), (9, 2), (11, 0), (15, 3), (16, 0), (18, 4)]
        self._test_to_genome(contiglen, exons, trpileip, expected)

        trpileip = [(5, 0), (8, 10), (13, 0)]
        expected = [(7, 0), (9, 10), (11, 0), (12, 10), (18, 0)]
        self._test_to_genome(contiglen, exons, trpileip, expected)

        trpileip = [(1, 10), (12, 0), (13, 10)]
        expected = [(1, 0), (2, 10), (17, 0), (18, 10)]
        self._test_to_genome(contiglen, exons, trpileip, expected)

        trpileip = [(1, 10), (3, 0), (4, 1), (6, 0), (8, 2), (13, 0)]
        expected = [(1, 0), (2, 10), (5, 0), (6, 1), (8, 0), (9, 2), (11, 0), (12, 2), (18, 0)]
        self._test_to_genome(contiglen, exons, trpileip, expected)

        # setup 2
        contiglen = 9
        exons = [(0, 1), (2, 3), (4, 5), (6, 7)]

        trpileip = [(3, 2), (4, 0)]
        expected = [(1, 2), (2, 0), (3, 2), (4, 0), (5, 2), (9, 0)]
        self._test_to_genome(contiglen, exons, trpileip, expected)

        trpileip = [(3, 0), (4, 1)]
        expected = [(6, 0), (7, 1), (9, 0)]
        self._test_to_genome(contiglen, exons, trpileip, expected)

        trpileip = [(2, 0), (3, 1), (4, 0)]
        expected = [(4, 0), (5, 1), (9, 0)]
        self._test_to_genome(contiglen, exons, trpileip, expected)

        # setup 3
        contiglen = 10
        exons = [(2, 7)]

        trpileip = [(1, 10), (4, 0), (5, 10)]
        expected = [(2, 0), (3, 10), (6, 0), (7, 10), (10, 0)]
        self._test_to_genome(contiglen, exons, trpileip, expected)

        trpileip = [(1, 0), (2, 1), (3, 0), (4, 5), (5, 0)]
        expected = [(3, 0), (4, 1), (5, 0), (6, 5), (10, 0)]
        self._test_to_genome(contiglen, exons, trpileip, expected)

        trpileip = [(5, 0)]
        expected = [(10, 0)]
        self._test_to_genome(contiglen, exons, trpileip, expected)

        trpileip = [(5, 10)]
        expected = [(2, 0), (7, 10), (10, 0)]
        self._test_to_genome(contiglen, exons, trpileip, expected)

    def _test_to_transcript(self, contiglen: int, exons: List[Tuple[int, int]],
                            blocks: List[List[Tuple[int, int]]], expected: List[List[Tuple[int, int]]]):
        exstart, exend = zip(*exons)
        trinfo = TranscriptInfo(
            ".", contiglen, np.asarray(exstart, dtype=np.int32), np.asarray(exend, dtype=np.int32), "TMP"
        )

        blocks = AlignedBlocks.from_tuples("", contiglen, blocks)

        beforebl, beforetr = deepcopy(blocks), deepcopy(trinfo)

        result = to_transcript(trinfo, blocks)

        trlen = sum(x[1] - x[0] for x in exons)
        expected = AlignedBlocks.from_tuples(trinfo.ensemblid, trlen, reads=expected)

        self.assertEqual(result.contig, expected.contig)
        self.assertEqual(result.contiglen, expected.contiglen)
        np.testing.assert_array_equal(result.records, expected.records)
        np.testing.assert_array_equal(result.start, expected.start)
        np.testing.assert_array_equal(result.end, expected.end)

        np.testing.assert_array_equal(trinfo.exstart, beforetr.exstart)
        np.testing.assert_array_equal(trinfo.exend, beforetr.exend)
        np.testing.assert_array_equal(beforebl.start, blocks.start)
        np.testing.assert_array_equal(beforebl.end, blocks.end)
        np.testing.assert_array_equal(beforebl.records, blocks.records)

    def test_to_transcript(self):
        # Case 1
        contiglen = 28
        exons = [(5, 8), (10, 11), (13, 15), (16, 17), (21, 26)]

        reads = [
            [(6, 8), (10, 11), (13, 14)],
            [(5, 8), (10, 11), (13, 15), (16, 17), (21, 26)],
            [(5, 6)],
            [(22, 25)],
            [(6, 8), (13, 14), (16, 17), (24, 26)],
            [(5, 6), (7, 8), (9, 11), (12, 17), (20, 23)],
            [(3, 27)],
            [(20, 21), (22, 23), (24, 25), (26, 27)],
            [(3, 12), (14, 18), (20, 21)],
            # Must not align
            [(9, 10)], [(12, 13), (15, 16), (17, 21)], [(26, 28)],
            [(0, 5), (18, 20), (27, 28)]
        ]
        expected = [
            [(1, 5)],
            [(0, 12)],
            [(0, 1)],
            [(8, 11)],
            [(1, 3), (4, 5), (6, 7), (10, 12)],
            [(0, 1), (2, 9)],
            [(0, 12)],
            [(8, 9), (10, 11)],
            [(0, 4), (5, 7)]
        ]
        self._test_to_transcript(contiglen, exons, reads, expected)

        # Case 2
        contiglen = 11
        exons = [(0, 2)]
        reads = [
            [(0, 1)],
            [(0, 2)],
            [(0, 2), (3, 4)],
            [(0, 1), (2, 3)],
            [(1, 2), (3, 4), (6, 8)],
            # Must not align
            [(3, 7), (9, 11)], [(3, 4), (5, 6), (7, 8), (9, 10)], [(2, 11)]
        ]
        expected = [
            [(0, 1)], [(0, 2)], [(0, 2)], [(0, 1)], [(1, 2)]
        ]
        self._test_to_transcript(contiglen, exons, reads, expected)

        # Case 3
        contiglen = 19
        exons = [(0, 1), (2, 5), (9, 12), (17, 19)]
        reads = [
            [(0, 2), (3, 4), (18, 19)],
            [(0, 1), (17, 19)],
            [(0, 2), (5, 9), (12, 18)],
            [(0, 13), (14, 15), (16, 17), (18, 19)],
            [(0, 6), (7, 8), (9, 18)]
        ]
        expected = [
            [(0, 1), (2, 3), (8, 9)], [(0, 1), (7, 9)], [(0, 1), (7, 8)], [(0, 7), (8, 9)],
            [(0, 8)]
        ]
        self._test_to_transcript(contiglen, exons, reads, expected)
