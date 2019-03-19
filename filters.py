# -*- coding: utf-8 -*-


import pandas as pd
from scipy.signal import savgol_filter


def sv(pxldrl, param):
    if param.smp != 0:  # TODO Check the smp value meanong
        # Savinsky Golet filter
        pxldrl.ps = savgol_filter(pxldrl.ts_d, param.medspan, param.smp, mode='nearest')
        # TODO automatic selection of savgol window
        return pd.Series(pxldrl.ps, pxldrl.ts_d.index)
    else:
        ps = pxldrl.ts_d.rolling(pxldrl.medspan // 2 * 2, win_type='boxcar', center=True).mean()