import pickle
from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Dict, Set, Tuple, Optional, List
from pathlib import Path

from scipy.stats import fisher_exact

from .misc import SitesThreshold


@dataclass(frozen=True)
class SiteStat:
    # __slots__ = ("reference", "edits")
    reference: int
    edits: int

    def freq(self) -> float:
        if self.coverage() == 0:
            return 0
        else:
            return self.edits / self.coverage()

    def coverage(self) -> int:
        return self.edits + self.reference


@dataclass(frozen=True)
class EditingSite:
    # __slots__ = ("contig", "pos", "trstrand", "genotyped")
    contig: str
    pos: int
    trstrand: str
    genotyped: Dict[Any, SiteStat]

    def allfreq(self) -> Tuple[float, ...]:
        return tuple(gtype.freq() for gtype in self.genotyped.values())

    def alledits(self) -> Tuple[int, ...]:
        return tuple(gtype.edits for gtype in self.genotyped.values())

    def allcoverage(self) -> Tuple[int, ...]:
        return tuple(gtype.coverage() for gtype in self.genotyped.values())

    def eratio(self) -> float:
        raise Exception
        koed, koref = sum(x.edits for x in self.genotyped['p150-KO']), sum(
            x.reference for x in self.genotyped['p150-KO'])
        wted, wtref = sum(x.edits for x in self.genotyped['WT']), sum(x.reference for x in self.genotyped['WT'])

        ko, wt = koed / (koed + koref), wted / (wted + wtref)
        return ko / (wt + ko)

    def passes(self, t: SitesThreshold) -> bool:
        if t.minsamples:
            ok = 0
            for x in self.genotyped.values():
                if t.minfreq and x.freq() < t.minfreq:
                    continue
                if t.mincoverage and x.coverage() < t.mincoverage:
                    continue
                if t.minedits and x.edits < t.minedits:
                    continue
                ok += 1
            return ok >= t.minsamples
        return True

    def __hash__(self):
        return hash(self.contig) + hash(self.pos) + hash(self.trstrand)


def esites(cache: Path, db: Set[Tuple[str, int, str]], t: SitesThreshold) -> List[EditingSite]:
    results = []
    with open(cache, 'rb') as stream:
        while True:
            try:
                site: EditingSite = pickle.load(stream)
                if (site.contig, site.pos, site.trstrand) in db or site.passes(t):
                    results.append(site)
            except EOFError:
                break
    return results


def fisher(x: SiteStat, y: SiteStat, pseudocounts: int = 1) -> float:
    # # Calculate p-values / q-values
    # def _job(workload) -> float:
    #     # pvalue = ttest_ind(*workload).pvalue
    #     pvalue = fisher_exact(*workload)[-1]
    #     # pvalue = boschloo_exact(*workload).pvalue
    #     return pvalue

    # workload = []
    # for island in clusters:
    #     # ko = list(zip(*[x.genotyped["p150-KO"] for x in cluster]))
    #     # wt = list(zip(*[x.genotyped["WT"] for x in cluster]))
    #     # assert len(wt) == len(ko) == 3
    #     #
    #     # eiko = [sum(site.edits for site in sample) / sum(site.edits + site.reference for site in sample) for sample in ko]
    #     # eiwt = [sum(site.edits for site in sample) / sum(site.edits + site.reference for site in sample) for sample in wt]
    #     # workload.append((eiko, eiwt))
    #
    #     matrix = np.zeros((2, 2), dtype=np.int32)
    #     for ind, genotype in enumerate(["WT", "p150-KO"]):
    #         for site in island:
    #             for stats in site.genotyped[genotype]:
    #                 matrix[ind, 0] += stats.reference
    #                 matrix[ind, 1] += stats.edits
    #     workload.append((matrix,))
    #
    # pvalues = Parallel(n_jobs=-1)(delayed(_job)(w) for w in workload)
    # qvalues = fdrcorrection(pvalues, alpha=FDR_THRESHOLDS)[1]

    table = [
        [x.reference + pseudocounts, x.edits + pseudocounts],
        [y.reference + pseudocounts, y.edits + pseudocounts]
    ]

    return fisher_exact(table)[-1]
