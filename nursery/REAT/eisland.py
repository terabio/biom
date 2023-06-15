import pickle
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .misc import EIslandThreshold
from .site import EditingSite, SiteStat


@dataclass(frozen=True)
class EditingIslandStats:
    key: Any
    # Median EI
    trtEI: float
    # Min mean coverage
    trtcov: float
    # Median EI
    cntEI: float
    # Min mean coverage
    cntcov: float
    # Media FC
    log2fc: float
    fisher_pvalue: float

    @property
    def direction(self):
        return "Treatment ↑" if self.trtEI > self.cntEI else "Treatment ↓"

    def to_dict(self) -> Dict[str, float]:
        return {
            "log2fc": self.log2fc, "p-value": self.fisher_pvalue,
            "trtcov": self.trtcov, "cntcov": self.cntcov,
            "trtEI": self.trtEI, "cntEI": self.cntEI, "direction": self.direction
        }


@dataclass(frozen=True)
class EditingIslandLocation:
    biotype: str
    location: str
    repeat: str
    weight: float

    def to_dict(self) -> Dict[str, float]:
        return {
            "biotype": self.biotype, "location": self.location, "repeat": self.repeat, "weight": self.weight
        }


@dataclass()
class EditingIsland:
    contig: str
    start: int
    end: int
    trstrand: str
    sitescnt: int
    sites: Tuple[EditingSite, ...]
    stats: Optional[List[EditingIslandStats]] = None
    location: Optional[List[EditingIslandLocation]] = None

    def hard_assigned_location(self):
        assert self.location is not None and len(self.location) > 0
        biotype, location, repeat = defaultdict(float), defaultdict(float), defaultdict(float)
        for x in self.location:
            biotype[x.biotype] += x.weight
            location[x.location] += x.weight
            repeat[x.repeat] += x.weight
        biotype = max(biotype.items(), key=lambda x: x[1])[0]
        location = max(location.items(), key=lambda x: x[1])[0]
        repeat = max(repeat.items(), key=lambda x: x[1])[0]
        return EditingIslandLocation(biotype, location, repeat, 1.0)

    # Stats
    # pvalue: float
    # qvalue: float

    # Folding
    # fold: str
    @staticmethod
    def from_sites(sites: List[EditingSite]) -> 'EditingIsland':
        assert all(x.contig == sites[0].contig and x.trstrand == sites[0].trstrand for x in sites)
        sites = sorted(sites, key=lambda x: x.pos)
        return EditingIsland(
            sites[0].contig, sites[0].pos, sites[-1].pos + 1, sites[0].trstrand,
            len(sites), sites=tuple(sites)
        )

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def sitefolds(self) -> List[str]:
        raise Exception

    def passses(self, t: EIslandThreshold) -> bool:
        strong, total = 0, 0
        for site in self.sites:
            if site.passes(t.strongthr):
                strong += 1
                total += 1
            elif site.passes(t.weakthr):
                total += 1
        return strong >= t.min_strong_sites and total >= t.min_sites
        # total = 0
        # for site in self.sites:
        #     if site.passes(t.weakthr):
        #         total += 1
        # return total >= t.min_sites

    def subsample(self, predicate, t: EIslandThreshold) -> Optional['EditingIsland']:
        sites = []
        for s in self.sites:
            genotyped = {k: v for k, v in s.genotyped.items() if predicate(k)}
            if genotyped:
                site = EditingSite(s.contig, s.pos, s.trstrand, genotyped)
                # if site.passes(t.weakthr):
                #     sites.append(site)
                sites.append(site)
        if sites:
            island = EditingIsland.from_sites(sites)
            # if island.passses(t):
            #     return island
            return island
        return None

    def pooled(self, key) -> Dict[Any, SiteStat]:
        reference, edits = defaultdict(int), defaultdict(int)
        for site in self.sites:
            for origin, stat in site.genotyped.items():
                origin = key(origin)
                reference[origin] += stat.reference
                edits[origin] += stat.edits
        return {k: SiteStat(reference[k], edits[k]) for k in reference}

        # mapped = []
        # for site in cluster:
        #     assert start < site.pos < end
        #     pos = site.pos - start
        #     mapped.append(structure[pos])
        # return mapped, (contig, start, end, structure)

    # fold: Tuple[str, int, int, str]
    # # Localization (biotype + location + repeat -> weight)
    # localization: Tuple[Tuple[Tuple[str, str, str], float], ...]
    #
    # def pulledei(self, genotype) -> float:
    #     ref, edits = 0, 0
    #     samples = [x.genotyped[genotype] for x in self.sites]
    #     samples = list(zip(*samples))
    #     assert len(samples) == 3
    #     for sample in samples:
    #         # reference = sum(site.reference for site in sample)
    #         # edits = sum(site.edits for site in sample)
    #         # ei[genotype].append(edits / (reference + edits))
    #         ref += sum(site.reference for site in sample)
    #         edits += sum(site.edits for site in sample)
    #     return edits / (ref + edits)
    #
    # def log2fc(self, nominator: str = "p150-KO", denominator: str = "WT"):
    #     return np.log2(self.pulledei(nominator) / (self.pulledei(denominator) + 1e-12))
    #
    # # def ratio(self) -> float:
    # #     Median EI
    # #     ei = {}
    # #     for genotype in "WT", "p150-KO":
    # #         samples = [x.genotyped[genotype] for x in self.sites]
    # #         samples = list(zip(*samples))
    # #         assert len(samples) == 3
    # #         values = []
    # #         for s in samples:
    # #             values.append(sum(site.edits for site in s) / sum(site.edits + site.reference for site in s))
    # #         ei[genotype] = values
    # #     wt = np.median(ei["WT"])
    # #     ko = np.median(ei["p150-KO"])
    # #
    # #     Pooled weighted EI
    # #     ei = {}
    # #     for genotype in "WT", "p150-KO":
    # #         ei[genotype] = {"ref": 0, "edits": 0}
    # #
    # #         samples = [x.genotyped[genotype] for x in self.sites]
    # #         samples = list(zip(*samples))
    # #         assert len(samples) == 3
    # #         for sample in samples:
    # #             # reference = sum(site.reference for site in sample)
    # #             # edits = sum(site.edits for site in sample)
    # #             # ei[genotype].append(edits / (reference + edits))
    # #             ei[genotype]['ref'] += sum(site.reference for site in sample)
    # #             ei[genotype]['edits'] += sum(site.edits for site in sample)
    # #     wt, ko = ei.pop("WT"), ei.pop("p150-KO")
    # #
    # #     wt, ko = np.mean(wt), np.mean(ko)
    # #     wt, ko = wt['edits'] / (wt['ref'] + wt['edits']), \
    # #              ko['edits'] / (ko['ref'] + ko['edits'])
    # #     wt, ko = self.pulledei("WT"), self.pulledei("p150-KO")
    # #     return ko / (wt + ko)
    #
    # def direction(self):
    #     wt, ko = self.pulledei("WT"), self.pulledei("p150-KO")
    #     direction = "p150-KO ↑" if ko > wt else "p150-KO ↓"
    #     if self.qvalue <= eislands.FDR_THR and abs(self.log2fc()) >= eislands.LOG2FC_THR:
    #         return direction
    #     else:
    #         return "Not significant"
    #
    # def biotype(self) -> str:
    #     biotype = None
    #     for (contig, start, end), handloc in utils.handfilter.LOCALIZATION.items():
    #         if contig == self.contig and min(end, self.end) - max(start, self.start) > 0:
    #             biotype = handloc
    #     if biotype:
    #         return biotype
    #
    #     weights = defaultdict(float)
    #     for (biotype, _, _), w in self.localization:
    #         weights[biotype] += w
    #     return max(weights.items(), key=lambda x: x[1])[0]
    #
    # def location(self) -> str:
    #     loc = None
    #     for (contig, start, end), handloc in utils.handfilter.LOCALIZATION.items():
    #         if contig == self.contig and min(end, self.end) - max(start, self.start) > 0:
    #             loc = handloc
    #     if loc:
    #         return loc
    #
    #     weights = defaultdict(float)
    #     for (_, loc, _), w in self.localization:
    #         weights[loc] += w
    #     return max(weights.items(), key=lambda x: x[1])[0]
    #
    # def repeat(self) -> str:
    #     weights = defaultdict(float)
    #     for (_, _, repeat), w in self.localization:
    #         weights[repeat] += w
    #     return max(weights.items(), key=lambda x: x[1])[0]
    #


# class seeds:
#     @staticmethod
#     @lru_cache(maxsize=1)
#     def fetch(fdr: Optional[float] = FDR_THR, log2fc: Optional[float] = LOG2FC_THR) -> List[EditingIsland]:
#         with open(utils.path.EDITING_ISLANDS_PKL, 'rb') as stream:
#             islands = pickle.load(stream)
#         if fdr:
#             islands = [x for x in islands if x.qvalue <= fdr]
#         if log2fc:
#             islands = [x for x in islands if abs(x.log2fc()) >= log2fc]
#         return islands


def save(cache: Path, islands: List[EditingIsland]):
    with open(cache, 'wb') as stream:
        pickle.dump(islands, stream, protocol=pickle.HIGHEST_PROTOCOL)


def fetch(cache: Path) -> List[EditingIsland]:
    with open(cache, 'rb') as stream:
        islands = pickle.load(stream)
    return islands


def fold(islands: List[EditingIsland]) -> List[EditingIsland]:
    raise NotImplementedError
    # # Predict RNA folding around editing islands
    # def _job(island: REAT.eisland.EditingIsland):
    #     start, end = island.start - c.eislands.FOLD_OFFSET, island.end + c.eislands.FOLD_OFFSET
    #     assert start >= 0 and end > start
    #
    #     # Fetch the sequence
    #     sequence = utils.features.sequence(island.contig, start, end)
    #     assert set(sequence) == {"A", "C", "G", "T"}
    #
    #     # Respect the strand
    #     if island.trstrand == "-":
    #         sequence = str(Seq.Seq(sequence).reverse_complement())
    #
    #     # Predict RNA structure
    #     structure, mfe = RNA.fold(sequence)
    #     assert mfe < 0, f"MFE: {mfe}"
    #     # Flip back to the forward orientation if needed
    #     if island.trstrand == "-":
    #         structure = structure[::-1]
    #
    #     offset = island.start - start, island.end - start
    #     structure = structure[offset[0]: offset[1]]
    #     assert len(structure) == island.length
    #     return structure
    #
    # folding = Parallel(n_jobs=-1, verbose=100)(delayed(_job)(i) for i in islands)
    #
    # for f, i in zip(folding, islands):
    #     assert len(f) == i.length
    #     i.fold = f
    # REAT.eisland.save(path.EDITING_ISLANDS, islands)
