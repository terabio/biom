from dataclasses import dataclass
from pathlib import Path
from typing import Any, List, Optional, Tuple

import numpy as np

from .. import coordmap, fragments, pileup
from ..config import PeakCallingConfig
from ..transcriptome import TranscriptInfo, TranscriptomeIndex
from ..utils import Stranded


@dataclass(frozen=True)
class Workload:
    contig: str
    bamfiles: List[Path]

    params: PeakCallingConfig.ProcessingParams
    # Transcriptome index for the given contig
    trid: Optional[TranscriptomeIndex]
    # Additional tags to identify this job later (I don't trust joblib order guarantee)
    tags: Any


@dataclass(frozen=True)
class Results:
    contig: str
    contiglen: np.int32
    fragments: int

    genomic: Stranded[pileup.Pileup]
    # list of (transcript, pileup in transcriptomic coordinates)
    transcriptomic: Stranded[List[Tuple[TranscriptInfo, pileup.Pileup]]]

    tags: Any


def _transcriptome(
        blocks: List[fragments.AlignedBlocks], trid: Optional[TranscriptomeIndex], strand: str, extensions: List[int]
) -> List[Tuple[TranscriptInfo, pileup.Pileup]]:
    if trid is None:
        return []

    # All transcripts must be sorted
    assert all(trid.exstart[trid.records[i]] <= trid.exend[trid.records[i + 1]] for i in range(trid.records.size - 2))
    # All blocks must be sorted
    for b in blocks:
        assert all(b.start[b.records[i]] <= b.start[b.records[i + 1]] for i in range(b.records.size - 2))

    results = []

    # Hinter pointer for the blocks that must overlap with the next transcript
    hinters = [0 for _ in blocks]

    for i in range(trid.records.size - 1):
        # Skip irrelevant transcripts
        if trid.strands[i] != strand:
            continue

        _st, _nd = trid.records[i], trid.records[i + 1]
        exstart, exend = trid.exstart[_st: _nd], trid.exend[_st: _nd]

        # Ignore transcript with a single exon = genomic pileup is better in this case
        # if exstart.size == 1:
        #     continue

        trinfo = TranscriptInfo(exstart, exend, trid.ensemblid[i], trid.strands[i])

        # Forward hinters if possible
        for ind, h in enumerate(hinters):
            bl = blocks[ind]
            while h < bl.records.size - 2:
                # End of the current block
                end = bl.end[bl.records[h + 1] - 1]
                if end > trinfo.exstart[0]:
                    break
                h += 1

            hinters[ind] = h

        mapped = []
        for h, block in zip(hinters, blocks):
            bl = coordmap.to_transcript(trinfo, block, h)
            if bl is not None:
                mapped.append(bl)

        # No reads overlapping the given transcript
        if len(mapped) == 0:
            continue

        trsize = np.sum(trinfo.exend - trinfo.exstart, dtype=np.int32)  # type: np.int32
        trpileup = [pileup.calculate(trid.contig, trsize, mapped, ex) for ex in extensions]
        trpileup = pileup.merge.by_max(trpileup)

        results.append((trinfo, trpileup))

    return results


def _genome(contig: str, contiglen: np.int32,
            blocks: List[fragments.AlignedBlocks], extensions: List[int]) -> pileup.Pileup:
    if len(blocks) == 0:
        return pileup.Pileup.constant(contig, contiglen, np.float32(0))

    pileups = [pileup.calculate(contig, contiglen, blocks, ext) for ext in extensions]
    return pileup.merge.by_max(pileups)


def run(workload: Workload) -> Results:
    extsize = workload.params.extsize[workload.contig]
    assert extsize and all(x >= 0 for x in extsize), f"Invalid extsize({extsize}) for contig {workload.contig}"

    # Load fragments
    blocks, contiglen = fragments.loadfrom(
        workload.bamfiles, fragments.strdeductors.get(workload.params.stranding), workload.contig,
        workload.params.inflags, workload.params.exflags, workload.params.minmapq
    )
    contiglen = np.int32(contiglen)
    total_fragments = sum(x.fragments() for x in blocks.fwd) + sum(x.fragments() for x in blocks.rev)

    # Calculate forward/reverse pileups for transcriptome and genome
    transcriptomic = Stranded(
        fwd=_transcriptome(blocks.fwd, workload.trid, '+', extsize),
        rev=_transcriptome(blocks.rev, workload.trid, '-', extsize)
    )
    genomic = Stranded(
        fwd=_genome(workload.contig, contiglen, blocks.fwd, extsize),
        rev=_genome(workload.contig, contiglen, blocks.rev, extsize)
    )
    return Results(
        contig=workload.contig,
        contiglen=contiglen,
        fragments=total_fragments,
        genomic=genomic,
        transcriptomic=transcriptomic,
        tags=workload.tags
    )
