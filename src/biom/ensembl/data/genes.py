from collections import defaultdict
from pathlib import Path
from typing import Optional, Tuple

import pandas as pd

DTYPES = {
    "Gene stable ID": str,
    "Gene stable ID version": str,
    "Transcript stable ID": str,
    "Transcript stable ID version": str,
    "Transcript support level (TSL)": str,
    "Gene name": str,
    "Transcript name": str,
    "Gene type": str,
    "Transcript type": str,
    "Gene Synonym": str,
    "Strand": int,
    "Transcript start (bp)": int,
    "Transcript end (bp)": int,
    "GENCODE basic annotation": str,
    "Chromosome/scaffold name": str,
}


class GenesInfo:
    def __init__(self, cache: Path):
        self.cache = cache

        # Prepare caches
        df = pd.read_csv(self.cache, sep="\t", dtype=DTYPES, usecols=list(DTYPES.keys()))
        assert (df['Strand'] != 0).all()
        df['Strand'] = df['Strand'].apply(lambda x: "+" if x > 0 else "-")

        idcolumns = ['Gene stable ID', 'Gene stable ID version']

        # Gene id -> gene coordinates
        self._id2coords = {}
        cols = ['Chromosome/scaffold name', 'Transcript start (bp)', 'Transcript end (bp)', 'Strand']
        for col in idcolumns:
            self._id2coords.update({
                tid: (contig, start, end, strand)
                for tid, contig, start, end, strand in zip(*[df[x] for x in [col, *cols]])
            })

        # Gene id -> gene biotype
        self._id2biotype = {}
        for col in idcolumns:
            self._id2biotype.update({k: v for k, v in zip(df[col].values, df["Gene type"].values)})

        # Gene id -> gene name
        self._id2name = {}
        for col in idcolumns:
            self._id2name.update({k: v for k, v in zip(df[col].values, df["Gene name"].values)})
        for k, v in self._id2name.items():
            if pd.isna(v):
                self._id2name[k] = None

        # Gene id -> all transcripts
        self._id2transcripts = defaultdict(set)
        for genecol in idcolumns:
            for tidcol in ['Transcript stable ID', 'Transcript stable ID version']:
                for gid, tid in zip(df[genecol], df[tidcol]):
                    self._id2transcripts[gid].add(tid)
        self._id2transcripts = {k: tuple(v) for k, v in self._id2transcripts.items()}

        # Gene name -> gene id
        self._name2id = {
            k: v for k, v in zip(df['Gene name'].values, df['Gene stable ID version'].values) if not pd.isna(k)
        }

        # Synonyms -> gene id
        self._synonym2id = self._name2id.copy()
        self._synonym2id.update(
            {k: v for k, v in zip(df['Gene Synonym'].values, df['Gene stable ID version'].values) if not pd.isna(k)}
        )
        # self._synonym2id.update({
        #     k: v
        #     for k, v in zip(df['NCBI gene (formerly Entrezgene) accession'].values, df['Gene stable ID version'].values)
        # })
        self._synonym2id.update({k.lower(): v for k, v in self._synonym2id.items()})
        self._synonym2id.update({k.upper(): v for k, v in self._synonym2id.items()})

    def _getter(self, mapping, id):
        if id not in mapping:
            id = id.split(".")[0]
        return mapping.get(id, None)

    def id2coords(self, id: str) -> Optional[Tuple[str, int, int, str]]:
        return self._getter(self._id2coords, id)

    def id2biotype(self, id: str) -> Optional[str]:
        return self._getter(self._id2biotype, id)

    def id2name(self, id: str) -> Optional[str]:
        return self._getter(self._id2name, id)

    def id2transcripts(self, id: str) -> Optional[Tuple[str, ...]]:
        return self._getter(self._id2transcripts, id)

    def name2id(self, name: str) -> Optional[str]:
        return self._getter(self._name2id, name)

    def synonym2id(self, synonym: str) -> Optional[str]:
        gid = self._getter(self._synonym2id, synonym)
        if gid is None:
            gid = self._getter(self._synonym2id, synonym.lower())
        if gid is None:
            gid = self._getter(self._synonym2id, synonym.upper())
        return gid
