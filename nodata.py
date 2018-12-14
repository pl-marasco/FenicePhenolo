# -*- coding: utf-8 -*-
# !/usr/bin/env python

import numpy as np
# import pandas as pd
# import xarray as xr
# from scipy.constants import year
# from datetime import datetime as dt
# import statsmodels.api as stm


def climate_fx(ts, **kwargs):

    tscp = ts.copy()

    param = kwargs.pop('settings', '')
    cloud = param.cloud  # 252 cloud
    snow = param.cloud  # 253 snow
    # 254 sea (not treated)
    nodata = param.mask  # 251/255 invalids or data error

    singleinterp = True

    # TODO would be correct eliminate possible trend before the groupby
    # Rought climatic indices without nan included
    tsc = tscp.mask(tscp > 250)
    clm = tsc.groupby([tscp.index.month, tscp.index.day]).median()
    clm_m = tsc.groupby([tscp.index.month]).median()
    clm_q = tsc.groupby([tscp.index.quarter]).median()
    clm_min = tsc.min()

    if singleinterp:
        # interpolate single values
        df_int = tsc.where(~ (tsc.shift(-1).notnull() & tsc.shift(1).notnull()),
                           tsc.interpolate(method='linear'))

    nans = tscp[tsc.isna()]

    for index, value in nans.items():
        agg_vls = clm.loc[index.month, index.day]
        if value in [251, 252, 255]:
            if 250 < agg_vls < 252.5 or agg_vls >= 254.5:
                if 250 < clm_m.loc[index.month] < 252.5 or clm_m.loc[index.month] > 254.5:
                    if 250 < clm_q.loc[index.quarter] < 252.5:
                        nans.loc[index] = clm_min
                    else:
                        nans.loc[index] = clm_q.loc[index.quarter]
                else:
                    if clm_m is None:
                        clm_m = tsc.groupby([tscp.index.month]).median()
                    nans.loc[index] = clm_m.loc[index.month]
            elif 252.5 <= agg_vls < 254.5:
                nans.loc[index] = 0
            else:
                nans.loc[index] = clm.loc[index.month, index.day]

        elif value == 253:
            if agg_vls > 250:
                nans.loc[index] = clm_min
            else:
                # nans.loc[index] = clm.loc[index.month, index.day] # TODO check alternatives
                nans.loc[index] = np.NaN

    tscp.update(nans)

    tscp[(tscp == 253) | (tscp == 252)] = np.NaN
    tscp = tscp.interpolate()

    return tscp

#
# def fix(param, cube, **kwargs):
#
#     type = kwargs.pop('type', '')
#
#     if type == 'SCE':
#
#         crd_x, crd_y, crd_t = _coord_names(cube)
#
#         if crd_x is not None and crd_y is not None:
#             cube = cube.load()
#             df_stk = cube.stack(z=(crd_y, crd_x)).to_pandas()
#             c3 = sce(df_stk)
#             return c3
#         else:
#             ts = cube.to_series()
#             ts_cln = pixeldrillcleaner(param, ts)
#             cube.values = ts_cln
#             return cube
#
#     else:
#         df_stk = cube.stack(z=('lat', 'lon')).to_pandas()
#         df_msk = df_stk.mask(df_stk > 250).interpolate()
#         res = stm.tsa.seasonal_decompose(df_msk, freq='', model='additive')
#
#
# def _df_expander(aggr, time):
#
#     res = pd.concat([aggr] * len(time.year.unique()), ignore_index=True)
#     res.set_axis(time, axis='index', inplace=True)
#     return res
#
#
# def sce(df):
#     # 251 invalid or no data or error
#     # 252 cloud
#     # 253 snow
#     # 254 sea (not treated)
#     # 255 invalid or no data or error
#
#     # mean without nan included
#     df_msk = df.mask(df > 250)
#
#     # interpolate single values
#     df_int = df_msk.where(~ (df_msk.shift(-1).notnull() &
#                           df_msk.shift(1).notnull()),
#                           df_msk.interpolate(method='time'))
#
#     # calculate time index
#     yrs = df_msk.index.year.unique()
#     time = pd.date_range('{}-01-01'.format(yrs[0]), '{}-12-21'.format(yrs[-1]))
#     time = time[time.day.isin([1, 11, 21])]
#
#     # climate indices
#     # day/month
#     d_aggr = df_msk.groupby([df_msk.index.month, df_msk.index.day]).median()
#     d_aggr = _df_expander(d_aggr, time)
#
#     # month aggr
#     mth_aggr = df_msk.groupby([df_msk.index.month]).median()
#     mth_aggr = _df_expander(pd.concat([mth_aggr] * 3).sort_index(), time)
#
#     # quarter aggr
#     quarter_aggr = df_msk.groupby([df_msk.index.quarter]).median()
#     mth_aggr = _df_expander(pd.concat([quarter_aggr] * 9).sort_index(), time)
#
#     # year aggr
#     year_aggr = df_msk.groupby([df_msk.index.year]).median()
#     year_aggr = pd.concat([year_aggr] * 36).sort_index()
#     year_aggr.set_axis(time, axis='index', inplace=True)
#
#     ts_min = df_msk.min()  #TODO add this as final substitution
#
#     df_int = df_int.where(df_int.notnull(), d_aggr)
#     df_int = df_int.where(df_int.notnull(), mth_aggr)
#     df_int = df_int.where(df_int.notnull(), quarter_aggr)
#     df_int = df_int.where(df_int.notnull(), year_aggr)
#
#     df_int = df_int.swapaxes('rows', 'columns').to_panel()
#     ds_int = df_int.to_xarray()
#
#     renamed = ds_int.rename({'items': 'time', 'major_axis': 'lat', 'minor_axis': 'lon'})
#
#     return renamed


#
# def _coord_names(data):
#     dims = [i.lower() for i in data.dims]
#     crd_x, crd_y, crd_t = [None] * 3
#
#     if 'lat' in dims or 'lon' in dims:
#         crd_x = data.dims[dims.index('lon')]
#         crd_y = data.dims[dims.index('lat')]
#         crd_t = data.dims[dims.index('time')]
#
#     elif 'e' in dims or 'n' in dims:
#         crd_x = data.dims[dims.index('e')]
#         crd_y = data.dims[dims.index('n')]
#         crd_t = data.dims[dims.index('time')]
#
#     else:
#         crd_t = data.dims[dims.index('time')]
#
#     return crd_x, crd_y, crd_t