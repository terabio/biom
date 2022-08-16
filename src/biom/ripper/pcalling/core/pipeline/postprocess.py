import logging
from dataclasses import dataclass
from typing import Any, List, Optional, Tuple

import numpy as np

from .pileup import Results
from .. import coordmap, pileup
from ..transcriptome import TranscriptInfo
from ..utils import Stranded


@dataclass(frozen=True)
class Workload:
    pileup: Results

    # Baseline values
    gmbaseline: np.float32
    trbaseline: np.float32
    minfragments: np.float32
    # Scaling coefficient
    scale: Optional[np.float32]


@dataclass(frozen=True)
class Result:
    contig: str
    contiglen: np.int32

    pileup: Stranded[pileup.Pileup]
    tags: Any


def calcscaling(trtfragments: int, cntfragments: int) -> Tuple[str, np.float32]:
    print(f"Treatment fragments: {trtfragments}, Control fragments: {cntfragments}")

    if trtfragments > cntfragments:
        scale = cntfragments / trtfragments
        scalelib = "treatment"
        logging.info(f"Treatment ({trtfragments}) is large than control ({cntfragments}), scaling multiplier: {scale}")
    else:
        scale = trtfragments / cntfragments
        scalelib = "control"
        logging.info(f"Control ({cntfragments}) is large than treatment ({trtfragments}), scaling multiplier: {scale}")
    return scalelib, np.float32(scale)


def _apply_baseline(p: Results, gmbaseline: np.float32, trbaseline: np.float32) -> Stranded[pileup.Pileup]:
    def _apply(trpileups: List[Tuple[TranscriptInfo, pileup.Pileup]], genome: pileup.Pileup) -> pileup.Pileup:
        # Apply baseline and map back to the genome
        mapped = []
        for trinfo, trpileup in trpileups:
            trpileup = pileup.merge.by_max([trpileup], trbaseline)
            mapped.append(coordmap.to_genome(p.contig, p.contiglen, trinfo, trpileup))
        trpileups = pileup.merge.by_max(mapped) if mapped else None

        if trpileups:
            return pileup.merge.by_max([trpileups, genome], baseline=gmbaseline)
        else:
            return pileup.merge.by_max([genome], baseline=gmbaseline)

    return Stranded(
        fwd=_apply(p.transcriptomic.fwd, p.genomic.fwd),
        rev=_apply(p.transcriptomic.rev, p.genomic.rev),
    )


def _normalize(p: Stranded[pileup.Pileup], scale: Optional[np.float32]) -> Stranded[pileup.Pileup]:
    if scale:
        p.fwd.values *= scale
        p.rev.values *= scale
    return p


def run(workload: Workload) -> Result:
    pileup = _apply_baseline(workload.pileup, workload.gmbaseline, workload.trbaseline)
    # Discard low-covered regions
    if workload.minfragments > 0:
        for x in pileup.fwd, pileup.rev:
            mask = x.values < workload.minfragments
            x.values[mask] = np.float32(0)

    pileup = _normalize(pileup, workload.scale)
    return Result(
        contig=workload.pileup.contig,
        contiglen=workload.pileup.contiglen,
        pileup=pileup,
        tags=workload.pileup.tags
    )
