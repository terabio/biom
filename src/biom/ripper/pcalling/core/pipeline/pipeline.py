import logging
from collections import defaultdict
from dataclasses import dataclass
from typing import List

import numpy as np
from joblib import Parallel, delayed

from . import pileup, postprocess
from ..config import PeakCallingConfig
from ..pileup import Pileup
from ..utils import fetch_contigs

TREATMENT = "treatment"
CONTROL = "control"


@dataclass(frozen=True)
class Results:
    contig: str
    contiglen: int
    trstrand: str

    trtpileup: Pileup
    cntpileup: Pileup

    def __post_init__(self):
        assert self.trtpileup.id == self.cntpileup.id == self.contig


def run(config: PeakCallingConfig, pool: Parallel) -> List[Results]:
    # Parse BAM contigs
    contigs = fetch_contigs(config.treatment + config.control)

    # Build & run pileup workloads
    contigs = tuple(contigs)
    workloads = []
    prconfigs = {
        CONTROL: config.process,
        # Disable extension for treatment
        TREATMENT: PeakCallingConfig.ProcessingParams(
            config.process.stranding, config.process.scaling,
            defaultdict(lambda *args: [0]), config.process.threads, config.process.backend,
            config.process.inflags, config.process.exflags, config.process.minmapq
        )
    }
    for contig in contigs:
        for tag, files in {TREATMENT: config.treatment, CONTROL: config.control}.items():
            workloads.append(pileup.Workload(
                contig=contig, bamfiles=files, params=prconfigs[tag], tags=tag
            ))
    results: List[pileup.Results] = pool(
        delayed(pileup.run)(w) for w in workloads
    )

    # Calculate baseline values
    trtfragments = sum(x.fragments for x in results if x.tags == TREATMENT)
    cntfragments = sum(x.fragments for x in results if x.tags == CONTROL)

    gmbaseline = cntfragments / config.geffsize
    print(f"Treatment fragments: {trtfragments}, Control fragments: {cntfragments}")
    print(f"Treatment scaling: {config.process.scaling.treatment}, "
          f"Control scaling: {config.process.scaling.control}")
    logging.info(f"Genome baseline: {gmbaseline}")

    # Build & run postprocess workloads
    workloads = []
    baselines = {TREATMENT: (0., config.callp.mintrtfrag), CONTROL: (gmbaseline, 0.)}
    for r in results:
        gmbaseline, minfragments = baselines[r.tags]

        if r.tags == TREATMENT:
            scale = config.process.scaling.treatment
        else:
            assert r.tags == CONTROL
            scale = config.process.scaling.control

        workloads.append(postprocess.Workload(
            pileup=r, gmbaseline=gmbaseline, scale=scale, minfragments=np.float32(minfragments)
        ))
    results: List[postprocess.Result] = pool(
        delayed(postprocess.run)(w) for w in workloads
    )
    # results: List[postprocess.Result] = [postprocess.run(w) for w in workloads]

    # Regroup results
    regrouped = defaultdict(dict)
    for r in results:
        # Forward
        key = (r.contig, r.contiglen, '+')
        assert key not in regrouped or r.tags not in regrouped[key]
        regrouped[key][r.tags] = r.pileup.fwd
        # Reverse
        key = (r.contig, r.contiglen, '-')
        assert key not in regrouped or r.tags not in regrouped[key]
        regrouped[key][r.tags] = r.pileup.rev

    # Final result
    results = []
    for (contig, contiglen, trstrand), pileups in regrouped.items():
        if CONTROL not in pileups or TREATMENT not in pileups:
            nonein = "control" if CONTROL not in pileups else "treatment"
            logging.warning(f"No {nonein} fragments for contig {contig}, strand {trstrand}")
            continue
        results.append(Results(
            contig=contig, contiglen=contiglen, trstrand=trstrand,
            trtpileup=pileups[TREATMENT], cntpileup=pileups[CONTROL]
        ))
    return results
