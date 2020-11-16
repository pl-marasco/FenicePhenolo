# -*- coding: utf-8 -*-
import phenolo.chronos as chronos
import pandas as pd


def climate_fx(ts, byr_date, eyr_date, **kwargs):

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
        agg = tsm.groupby([ts.index.month, ts.index.day])
        count = agg.count()
        cl = agg.median()
        clm = cl.mask(count < count.max() * 0.2)

        if clm.isnull().sum():
            clm = clm.mask(clm.isnull(), tsm.groupby([ts.index.month]).min())
            if clm.isnull().sum():
                clm = clm.mask(clm.isnull(), tsm.groupby([ts.index.quarter]).min())
                if clm.isnull().sum():
                    clm = clm.mask(clm.isnull(), tsm.min())

        dataindex = chronos.create(f'1/1/{byr_date}', f'31/12/{eyr_date}', 's10') # TODO make it flexible !!!
        values = pd.concat([clm] * dataindex.year.unique().size)
        climate = pd.Series(values.values, index=dataindex)
        tsm.where(tsm.notnull(), climate, inplace=True)

        # nans = ts.where(tsm.isnull()).dropna()
        #
        # for ith in nans.index:
        #     tsm.loc[ith] = clm.loc[ith.month, ith.day]

    return tsm
