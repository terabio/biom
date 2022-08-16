from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import numpy.typing as npt
from joblib import Parallel, delayed
from pybedtools import BedTool, Interval

from modules import ensembl


@dataclass
class TranscriptInfo:
    # Exon coordinates for the selected transcript
    exstart: npt.NDArray[np.int32]
    exend: npt.NDArray[np.int32]
    # Ensembl Transcript Ids
    ensemblid: str
    # Transcript strand
    trstrand: str


@dataclass
class TranscriptomeIndex:
    @dataclass(frozen=True)
    class Range:
        start: int
        end: int

        def __post_init__(self):
            assert self.end > self.start

        def overlap(self, other: 'Range') -> int:
            return min(self.end, other.end) - max(self.start, other.start)

        def __len__(self):
            return self.end - self.start

    EXONCOORDS = List[Range]

    contig: str
    # "Simplified" exons for each transcript
    exstart: npt.NDArray[np.int32]
    exend: npt.NDArray[np.int32]
    # Exons for i-th transcripts: exons [trindex[i]: trindex[i+1]]
    records: npt.NDArray[np.int32]
    # Ensembl Transcript Ids
    ensemblid: List[str]
    # Transcription strands
    strands: List[str]

    @staticmethod
    def simplify(exons: List[EXONCOORDS], trids: List[str], *args) -> Tuple[List[EXONCOORDS], List[str], Any]:
        assert len(exons) == len(trids)
        # I can +- safely merge transcripts where:
        # (1): last exons differ in 5` length (pick the longest)
        # (2): first exon differ in 3` length (pick the longest)
        # (3): intermediate exon differ in 3` length (pick the shortest)
        if len(exons) == 1:
            return exons, trids, args

        results = ([], [])
        while exons:
            # 1 - pick seed, longest transcript
            seed = max(range(len(exons)), key=lambda x: exons[x][-1].end - exons[x][0].start)
            seedex = exons[seed]
            # 2 - select all "collapsing" candidates
            collapsed = []
            for ind in range(len(exons)):
                candidate = exons[ind]

                # fast skipping
                if ind == seed:
                    continue
                # More exons => skip
                if len(candidate) > len(seedex):
                    continue
                # Seed must be the longest transcript
                assert seedex[-1].end - seedex[0].start >= candidate[-1].end - candidate[0].start

                # Fetch window where seed exons contain this candidate
                left = 0
                while left < len(seedex) and candidate[0].overlap(seedex[left]) < 0:
                    left += 1
                if left == len(seedex):
                    continue

                right = len(seedex) - 1
                while right >= 0 and candidate[-1].overlap(seedex[right]) < 0:
                    right -= 1
                if right < 0:
                    continue

                window = seedex[left: right + 1]
                assert len(window) >= 1

                # seed can be longer in first/last exon
                if window[0].end != candidate[0].end or window[-1].start != candidate[-1].start:
                    continue

                # other exons MUST match
                window, candidate = window[1: -1], candidate[1: -1]
                if len(window) != len(candidate):
                    continue
                if any(w.start != c.start or w.end != c.end for w, c in zip(window, candidate)):
                    continue

                # Compatible
                collapsed.append(ind)

            results[0].append(seedex)
            results[1].append(trids[seed] + ';'.join([trids[cind] for cind in collapsed]))

            removed = {seed, *collapsed}
            exons = [x for ind, x in enumerate(exons) if ind not in removed]
            trids = [x for ind, x in enumerate(trids) if ind not in removed]
        return results[0], results[1], args

    @staticmethod
    def build(exons: Path, transcripts: Optional[Set[str]] = None) -> Dict[str, 'TranscriptomeIndex']:
        # Gather transcripts per gene/contig: contig -> (gene, strand) -> (trid) -> exons
        models = defaultdict(lambda *args: defaultdict(list))

        for record in BedTool(exons):  # type: Interval
            if transcripts is not None and record.name not in transcripts:
                continue

            gene = ensembl.transcripts.id2gene(record.name)
            gene = gene if gene else record.name

            models[(record.chrom, gene, record.strand)][record.name].append(
                TranscriptomeIndex.Range(record.start, record.end)
            )

        beforeop = 0
        jobs = []
        for (contig, _, strand), transcripts in models.items():
            exons, trids = [], []
            for trid, tridexons in transcripts.items():
                trids.append(trid)
                tridexons = sorted(tridexons, key=lambda x: x.start)
                exons.append(tridexons)
            beforeop += len(trids)
            jobs.append((exons, trids, contig, strand))
        del models

        transcripts = Parallel(n_jobs=-1)(delayed(TranscriptomeIndex.simplify)(*j) for j in jobs)
        # transcripts = [TranscriptomeIndex.simplify(*j) for j in jobs]
        del jobs

        # Group by contig
        afterop = 0
        grouped = defaultdict(list)
        for exons, trids, (contig, strand) in transcripts:
            afterop += len(trids)
            grouped[contig].append((exons, trids, strand))
        del transcripts
        print(f"After optimization: {afterop / beforeop * 100: .2f}% transcripts retained")

        # Build dense representation, ie TranscriptomeIndex
        results = {}
        for contig, transcripts in grouped.items():

            merged = ([], [], [])
            for exons, trids, strand in transcripts:
                merged[0].extend(exons)
                merged[1].extend(trids)
                merged[2].extend([strand for _ in exons])
            order = sorted(range(len(merged[0])), key=lambda x: merged[0][x][0].start)
            exons = [merged[0][ind] for ind in order]
            trids = [merged[1][ind] for ind in order]
            trstrands = [merged[2][ind] for ind in order]
            assert len(exons) == len(trids) == len(trstrands)

            exstart, exend, trindex, ensemblid, strands = [], [], [], [], []
            curind = 0
            for blocks, trid, strand in zip(exons, trids, trstrands):
                # Fill meta
                trindex.append(curind)
                ensemblid.append(trid)
                strands.append(strand)
                # Double check exons, must be sorted
                assert all(curexon.end < nxtexon.start for curexon, nxtexon in zip(blocks[:-1], blocks[1:]))
                # Fill exons coordinates
                for x in blocks:
                    exstart.append(x.start)
                    exend.append(x.end)
                    curind += 1
            # All transcripts must be sorted by the start position
            assert all(exstart[trindex[i]] <= exstart[trindex[i + 1]] for i in range(len(trindex) - 1))

            # Final index for consistency
            assert curind == len(exstart) == len(exend)
            trindex.append(curind)

            assert contig not in results
            results[contig] = TranscriptomeIndex(
                contig,
                exstart=np.asarray(exstart, dtype=np.int32),
                exend=np.asarray(exend, dtype=np.int32),
                records=np.asarray(trindex, dtype=np.int32),
                ensemblid=ensemblid,
                strands=strands
            )
        assert sum(len(x.ensemblid) for x in results.values()) == afterop
        return results
