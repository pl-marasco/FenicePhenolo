# -*- coding: utf-8 -*-


import pandas as pd
from scipy.signal import savgol_filter


def sv(pxldrl, smp, medspan):
    if smp != 0:  # TODO Check the smp value meanong
        # Savinsky Golet filter
        ts = savgol_filter(pxldrl.ts, medspan, smp, mode='nearest')
        # TODO automatic selection of savgol window
        return pd.Series(ts, pxldrl.ts.index)
    else:
        return pxldrl.ts.rolling(pxldrl.medspan // 2 * 2, center=True).mean()
