from dataclasses import dataclass
from pathlib import Path
from typing import Any, List

import numpy as np

from .. import fragments, pileup
from ..config import PeakCallingConfig
from ..utils import Stranded


@dataclass(frozen=True)
class Workload:
    contig: str
    bamfiles: List[Path]

    params: PeakCallingConfig.ProcessingParams
    # Additional tags to identify this job later (I don't trust joblib order guarantee)
    tags: Any


@dataclass(frozen=True)
class Results:
    contig: str
    contiglen: np.int32
    fragments: int
    genomic: Stranded[pileup.Pileup]
    tags: Any


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

    genomic = Stranded(
        fwd=_genome(workload.contig, contiglen, blocks.fwd, extsize),
        rev=_genome(workload.contig, contiglen, blocks.rev, extsize)
    )
    return Results(
        contig=workload.contig,
        contiglen=contiglen,
        fragments=total_fragments,
        genomic=genomic,
        tags=workload.tags
    )
