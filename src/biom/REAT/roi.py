from pathlib import Path

import numpy as np
import pandas as pd

from .misc import DTYPES_ROIS
from .. import repmasker


def fetch_repmasker(path: Path, groupby: str = "family") -> pd.DataFrame:
    assert groupby in {"family", "name", "cls"}

    revcomp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    matches = [f"{n}->{n}" for n in "ACGT"]
    mismatches = [f"{n1}->{n2}" for n1 in "ACGT" for n2 in "ACGT" if n1 != n2]

    df = pd.read_csv(path, dtype=DTYPES_ROIS)
    print(path)
    for col in matches + mismatches:
        assert (df[col] >= 0).all()
        df[col] = df[col].astype(np.float32, copy=False)

    df['family'], df['cls'] = zip(*df['name'].apply(lambda x: repmasker.classify(x)[1:]))

    # Calculate EI
    ei = []
    for group, subdf in df.groupby(groupby):
        record = {"File": path.name, "Group": group}
        for fromn in "ACGT":
            for ton in "ACGT":
                # + strand
                key = f"{fromn}->{ton}"
                coverage = [f"{fromn}->{x}" for x in "ACGT"]
                pluscov = subdf.loc[subdf['trstrand'] == '+', coverage].sum().sum()
                plusmismatches = subdf.loc[subdf['trstrand'] == '+', key].sum()
                # - strand
                key = f"{revcomp[fromn]}->{revcomp[ton]}"
                coverage = [f"{revcomp[fromn]}->{x}" for x in "ACGT"]
                mnscov = subdf.loc[subdf['trstrand'] == '-', coverage].sum().sum()
                mnsmismatches = subdf.loc[subdf['trstrand'] == '-', key].sum()

                record[f'{fromn}->{ton}'] = (plusmismatches + mnsmismatches) / ((pluscov + mnscov) + 1e-6)
        ei.append(record)
    return pd.DataFrame(ei)


def fetchASDA(path: Path) -> pd.DataFrame:
    revcomp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    matches = [f"{n}->{n}" for n in "ACGT"]
    mismatches = [f"{n1}->{n2}" for n1 in "ACGT" for n2 in "ACGT" if n1 != n2]

    df = pd.read_csv(path, dtype=DTYPES_ROIS)
    print(path)
    for col in matches + mismatches:
        assert (df[col] >= 0).all()
        df[col] = df[col].astype(np.float32, copy=False)

    # Calculate EI
    ei = []
    # drop incorrect regions
    df = df[df['trstrand'] == df['strand']]

    for (group, strand), subdf in df.groupby(['name', 'trstrand']):
        record = {"File": path.name, "Group": group, "Strand": strand}
        for fromn in "ACGT":
            for ton in "ACGT":
                # + strand
                key = f"{fromn}->{ton}"
                coverage = [f"{fromn}->{x}" for x in "ACGT"]
                pluscov = subdf.loc[subdf['trstrand'] == '+', coverage].sum().sum()
                plusmismatches = subdf.loc[subdf['trstrand'] == '+', key].sum()
                # - strand
                key = f"{revcomp[fromn]}->{revcomp[ton]}"
                coverage = [f"{revcomp[fromn]}->{x}" for x in "ACGT"]
                mnscov = subdf.loc[subdf['trstrand'] == '-', coverage].sum().sum()
                mnsmismatches = subdf.loc[subdf['trstrand'] == '-', key].sum()

                record[f'{fromn}->{ton}'] = (plusmismatches + mnsmismatches) / ((pluscov + mnscov) + 1e-6)
        ei.append(record)
    return pd.DataFrame(ei)

# @staticmethod
# def fetch(folder: Path, groupby: str = "family") -> pd.DataFrame:
#     assert folder.is_dir()
#
#     dfs = Parallel(n_jobs=-1)(
#         delayed(repeats._fetch)(file, groupby) for file in folder.iterdir() if file.is_file())
#     return pd.concat(dfs)
