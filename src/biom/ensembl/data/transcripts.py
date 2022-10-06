from pathlib import Path
from typing import Optional, Tuple

import pandas as pd

DTYPES = {
    "Gene stable ID version": str,
    "Transcript stable ID": str,
    "Transcript stable ID version": str,
    "Transcript support level (TSL)": str,
    "Transcript name": str,
    "Transcript type": str,
    "GENCODE basic annotation": str,
    "Strand": int,
    "Transcript start (bp)": int,
    "Transcript end (bp)": int,
    "Chromosome/scaffold name": str,
}


class TranscriptsInfo:
    def __init__(self, cache: Path):
        self.cache = cache

        # Prepare caches
        df = pd.read_csv(self.cache, sep="\t", dtype=DTYPES, usecols=list(DTYPES.keys()))
        assert (df['Strand'] != 0).all()
        df['Strand'] = df['Strand'].apply(lambda x: "+" if x > 0 else "-")

        idcolumns = ['Transcript stable ID', 'Transcript stable ID version']

        # Transcript id -> transcript coordinates
        self._id2coords = {}
        cols = ['Chromosome/scaffold name', 'Transcript start (bp)', 'Transcript end (bp)', 'Strand']
        for k in idcolumns:
            self._id2coords.update({
                tid: (contig, start, end, strand)
                for tid, contig, start, end, strand in zip(*[df[x] for x in [k, *cols]])
            })

        # Transcript id -> transcript biotype
        self._id2biotype = {}
        for k in idcolumns:
            self._id2biotype.update({k: v for k, v in zip(df[k].values, df["Transcript type"].values)})

        # Transcript id -> host gene
        self._id2gene = {}
        for k in idcolumns:
            self._id2gene.update({k: v for k, v in zip(df[k].values, df['Gene stable ID version'].values)})

    def _getter(self, mapping, id):
        if id not in mapping:
            id = id.split(".")[0]
        return mapping.get(id, None)

    def id2coords(self, id: str) -> Optional[Tuple[str, int, int, str]]:
        return self._getter(self._id2coords, id)

    def id2biotype(self, id: str) -> Optional[str]:
        return self._getter(self._id2biotype, id)

    def id2gene(self, id: str) -> Optional[str]:
        return self._getter(self._id2gene, id)

# @lru_cache(maxsize=None)
# def incomplete(gtf: Path) -> Set[str]:
#     transcripts = set()
#     with gzip.open(gtf, 'rt') as stream:
#         for line in stream:
#             line = line.strip()
#             if len(line) == 0 or line.startswith("#"):
#                 continue
#             line = line.split("\t")
#             if line[2] != 'transcript':
#                 continue
#             fields = line[-1].split(";")
#             tid = None
#             for field in fields:
#                 if len(field) == 0:
#                     continue
#                 k, v = field.split(maxsplit=1)
#                 if k == 'transcript_id':
#                     assert v[0] == v[-1] == '"'
#                     tid = v[1:-1]
#                 elif k == 'tag' and v in {'"cds_end_NF"', '"cds_start_NF"', '"mRNA_end_NF"', '"mRNA_start_NF"'}:
#                     assert tid is not None
#                     # print(field)
#                     transcripts.add(tid)
#     return transcripts
