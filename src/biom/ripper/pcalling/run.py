import copy

from joblib import Parallel, delayed

from . import core
from .core import pipeline
from .core.config import PeakCallingConfig
from itertools import chain


def run(config: PeakCallingConfig):
    with Parallel(n_jobs=config.process.threads, backend=config.process.backend) as pool:
        pileups = pipeline.run(config, pool)

        # Convert to tracks and save pileups
        if config.saveto.pileup:
            for key, title in (lambda x: x.trtpileup, f"{config.saveto.title}.trt"), \
                              (lambda x: x.cntpileup, f"{config.saveto.title}.cnt"):
                tracks = [core.functors.Result.from_pileup(key(x), x.contiglen, x.trstrand) for x in pileups]
                core.io.tobigwig(tracks, config.saveto.pileup, title)
                del tracks

        # Calculate fold enrichment
        # fe = [core.functors.foldenrichment.calculate(w) for w in pileups]
        fe = pool(delayed(core.functors.foldenrichment.calculate)(w) for w in pileups)
        if config.saveto.enrichment:
            core.io.tobigwig(fe, config.saveto.enrichment, config.saveto.title)

        if config.saveto.pvpeaks is None and config.saveto.fdrpeaks is None:
            return

        # Calculate p-values
        # pvalues = [core.functors.pvalues.calculate(w) for w in pileups]
        pvalues = pool(delayed(core.functors.pvalues.calculate)(w) for w in pileups)
        pvalues, pcounts = zip(*pvalues)

        # core.io.tobigwig(pvalues, config.saveto.enrichment, f"{config.saveto.title}.pvalue")

        # Calculate q-values
        pqtable = core.functors.qvalues.make_pqtable(pcounts)
        # qvalues = [core.functors.qvalues.apply_pqtable(w, pqtable) for w in pvalues]
        qvalues = pool(delayed(core.functors.qvalues.apply_pqtable)(w, pqtable) for w in pvalues)
        # core.io.tobigwig(pvalues, config.saveto.enrichment, f"{config.saveto.title}.qvalue")

        # Call peaks using various cutoffs
        workload = []
        if config.saveto.pvpeaks and config.callp.pvcutoff is not None:
            callp = copy.deepcopy(config.callp)
            callp.qvcutoff = None
            workload.append((callp, config.saveto.pvpeaks))
        if config.saveto.fdrpeaks and config.callp.qvcutoff is not None:
            callp = copy.deepcopy(config.callp)
            callp.pvcutoff = None
            workload.append((callp, config.saveto.fdrpeaks))

        for callp, saveto in workload:
            workload = core.functors.callpeaks.PeakCalingWorkload.build(pvalues, qvalues, fe, callp)
            # peaks = [core.functors.callpeaks.calculate(w) for w in workload]
            peaks = pool(delayed(core.functors.callpeaks.calculate)(w) for w in workload)
            peaks = list(chain(*peaks))

            core.io.tobed(peaks, saveto.joinpath(f"{config.saveto.title}.narrowPeak"))
