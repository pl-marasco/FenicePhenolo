# -*- coding: utf-8 -*-


import pandas as pd
from scipy.signal import savgol_filter


def sv(ts_d, smp, medspan):
    if smp != 0:  # TODO Check the smp value meanong
        # Savinsky Golet filter
        # TODO automatic selection of savgol window
        return pd.Series(savgol_filter(ts_d.values, medspan, smp, mode='nearest'), ts_d.index)
    else:
        ps = ts_d.rolling(medspan // 2 * 2, center=True).mean()
