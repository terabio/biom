from collections import defaultdict
from pathlib import Path

from src.biom import gindex

first = gindex.from_bed(
    Path("/home/alnfedorov/projects/Z-DoG/refactor/A-to-I editing/ld/results/Z22 RIP-seq.HSV-1.bed"),
    datafn=lambda _: 'HSV-1'
)
second = gindex.from_bed(
    Path("/home/alnfedorov/projects/Z-DoG/refactor/A-to-I editing/ld/results/Z22 RIP-seq.mock.bed"),
    datafn=lambda _: 'mock'
)

merged = gindex.merge(first, second)

for x in first, second, merged:
    cnts = defaultdict(int)
    for val in x.itrees.values():
        for v in val.all_intervals:
            cnts[v.data] += 1
    print(cnts)

    print(sum(len(val) for val in x.itrees.values()))
