from pathlib import Path
from typing import Callable, List

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from .misc import DTYPES_ROIS


def summarize(path: Path, groupby: Callable[[str], str]) -> pd.DataFrame:
    revcomp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    matches = [f"{n}->{n}" for n in "ACGT"]
    mismatches = [f"{n1}->{n2}" for n1 in "ACGT" for n2 in "ACGT" if n1 != n2]

    df = pd.read_csv(path, dtype=DTYPES_ROIS)
    for col in matches + mismatches:
        assert (df[col] >= 0).all()
        df[col] = df[col].astype(np.float32, copy=False)

    df['group'] = df['name'].apply(groupby)

    # Calculate EI
    ei = []
    for group, subdf in df.groupby('group'):
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


def fetch(files: List[Path], groupby: Callable[[str], str]) -> pd.DataFrame:
    assert all(x.is_file() for x in files)

    dfs = Parallel(n_jobs=-1)(delayed(summarize)(file, groupby) for file in files)
    return pd.concat(dfs)
