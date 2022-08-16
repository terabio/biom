from pathlib import Path
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
from scipy import stats

from .config import Scaling


def median_of_ratios(bins: Path, treatment: List[str], control: List[str]) -> Scaling:
    bins = np.load(bins.as_posix())
    df = pd.DataFrame(data=bins['matrix'], columns=bins['labels'])
    df['treatment'] = df[treatment].sum(axis=1)
    df['control'] = df[control].sum(axis=1)
    df = df[['treatment', 'control']]

    df = df[(df > 0).all(axis=1)].copy()

    df['GeometricMean'] = df.apply(stats.gmean, axis=1)

    trtscale = (df['treatment'] / df['GeometricMean']).median()
    cntscale = (df['control'] / df['GeometricMean']).median()
    return Scaling(np.float32(1 / trtscale), np.float32(1 / cntscale))
