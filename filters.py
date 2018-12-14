import pandas as pd
from scipy.signal import savgol_filter


def sv(pxdrl, param):
    if param.smp != 0:  # TODO Check the smp value meanong
        # Savinsky Golet filter
        pxdrl.ps = savgol_filter(pxdrl.ts_d, param.medspan, param.smp, mode='nearest')
        # TODO automatic selection of savgol window
        return pd.Series(pxdrl.ps, pxdrl.ts_d.index)
    else:
        ps = pxdrl.ts_d.rolling(pxdrl.medspan // 2 * 2, win_type='boxcar', center=True).mean()