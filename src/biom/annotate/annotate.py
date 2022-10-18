from HTSeq import GenomicArrayOfSets, GenomicInterval
from collections import defaultdict
from typing import Dict, Optional, Set, Tuple

from .. import ensembl, repmasker


def byrepeats(classification: repmasker.RepmaskerClassification,
              repeats: GenomicArrayOfSets, key: GenomicInterval, how: str = "name") -> Dict[str, float]:
    assert how in {"name", "class", "family"}
    ind = {"name": 0, "family": 1, "class": 2}[how]
    weights = defaultdict(int)
    if key.chrom in repeats.chrom_vectors:
        for iv, overlaps in repeats[key].steps():
            if iv.length > 0:
                if len(overlaps) == 1:
                    for name, _ in overlaps:
                        weights[classification.classify(name)[ind]] += iv.length
                else:
                    # Same overlap -> fractional weights
                    for name, _ in overlaps:
                        weights[classification.classify(name)[ind]] += iv.length / len(overlaps)
    covered = sum(weights.values())

    if covered < key.length:
        weights['Repeat-free'] = key.length - covered
    covered = sum(weights.values())
    assert abs(covered - key.length) < 1e-3, f"{key}: {covered}"

    weights = {k: v / key.length for k, v in weights.items()}
    assert abs(sum(weights.values()) - 1) <= 1e-6
    return weights


def bytranscriptome(
        assembly: ensembl.EnsemblAssembly,
        trindex: GenomicArrayOfSets, key: GenomicInterval, strand: str, expressed: Optional[Set[str]] = None,
        priority: Tuple[str, ...] = ("cds", "utr3", "utr5", "exon", "intron", "upstream-1k", "downstream+10k",),
) -> Tuple[Dict[str, float], Dict[str, Set[str]]]:
    assert strand in ("+", "-", ".")
    priority = {k: ind for ind, k in enumerate(priority)}

    weights = defaultdict(int)
    matched_transcripts = defaultdict(set)
    if key.chrom in trindex.chrom_vectors:
        for iv, transcripts in trindex[key].steps():
            if iv.length > 0:
                options = defaultdict(set)
                for feature, trid, trstrand in transcripts:
                    if trstrand == strand or strand == ".":
                        biotype = assembly.transcripts.id2biotype(trid)
                        assert biotype is not None, f"Unknown transcript ID: {trid}"
                        if (biotype == "protein_coding" or biotype == "lncRNA") and \
                                (expressed is None or trid in expressed):
                            options[feature].add(trid)
                        # These RNAs can be part of intronic sequences and not transcribed independently
                        elif biotype == "miRNA" or biotype == "scRNA" or biotype == "snoRNA" or biotype == "ncRNA":
                            options[feature].add(trid)
                        else:
                            options[feature].add(trid)
                if options:
                    feature, trids = min(options.items(), key=lambda x: priority[x[0]])
                    weights[feature] += iv.length
                    matched_transcripts[feature] |= trids

    covered = sum(weights.values())
    if covered != key.length:
        weights['intergenic'] = key.length - covered
    weights = {k: v / key.length for k, v in weights.items()}
    assert abs(sum(weights.values()) - 1) <= 1e-6
    return weights, matched_transcripts
