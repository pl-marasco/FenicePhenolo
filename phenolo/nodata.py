# -*- coding: utf-8 -*-


def climate_fx(ts, **kwargs):
    singleinterp = True

    tsm = ts.mask(ts > 250)

    # interpolate single values
    if singleinterp:
        ts_s = tsm.where(~(tsm.shift(-1).notnull() & tsm.shift(1).notnull()), tsm.interpolate(method='linear'))
    else:
        ts_s = tsm

    tsm = ts_s.mask(ts == 253, 0)

    if ts_s.isnull().sum():
        # Rough climatic indices without nan included
        count = tsm.groupby([ts.index.month, ts.index.day]).count()
        cl = tsm.groupby([ts.index.month, ts.index.day]).median()
        clm = cl.mask(count < count.max() * 0.2)

        if clm.isnull().sum():
            clm = clm.mask(clm.isnull(), tsm.groupby([ts.index.month]).min())
        if clm.isnull().sum():
            clm = clm.mask(clm.isnull(), tsm.groupby([ts.index.quarter]).min())
        if clm.isnull().sum():
            clm = clm.mask(clm.isnull(), tsm.min())

        nans = ts.where(tsm.isnull()).dropna()

        for ith in nans.index:
            tsm.loc[ith] = clm.loc[ith.month, ith.day]

    return tsm
