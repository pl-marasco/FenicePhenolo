# -*- coding: utf-8 -*-
# !/usr/bin/env python

import nodata; import metrics; import outlier
import logging
import numpy as np
from seasonal import fit_seasons


logger = logging.getLogger(__name__)


def preingest(tuple, **kwargs):

    param = kwargs.pop('settings', '')
    tse = None

    # no data removing
    if param.sensor == 'spot':
        tuple.ts = nodata.climate_fx(tuple.ts, settings=param)

    # scaling
    if param.scale is not None:
        tuple.ts = metrics.scale(tuple.ts, settinngs=param)

    # off set
    if param.offset is not None:
        tuple.ts = metrics.offset(tuple.ts, settings=param)

    # rescaling to 0-100
    if param.min is not None and param.max is not None:
        tuple.ts = metrics.rescale(tuple.ts, settings=param)

    return tuple


def cleaning(tuple, **kwargs):
    param = kwargs.pop('settings', '')

    ts_clean = None

    # Filter outlier
    try:
        tuple.ts_t = outlier.season(tuple.ts, 1, param.yr_dys * param.outmax, 3)  # TODO make variable dek
    except (RuntimeError, Exception, ValueError):
        logger.debug('Error in filter outlayer')
        # print('Outlier research failed, in position:{0}'.format(coord))
        # err_table['outlier'] = 1
        # return err_table

    try:
        if tuple.ts_t is not None:
            # ts_cleaned = outlier.fillerSeason(ts_clean)
            tuple.ts_t = tuple.ts_t.interpolate()
    except (RuntimeError, Exception, ValueError):
        # print('Cleaning failed, in position:{0}'.format(coord))
        # err_table['cleaning'] = 1
        # return err_table
        pass

    # Estimate Season length
    try:
        tuple.seasons, tuple.trend = fit_seasons(ts_cleaned)
    except (RuntimeError, Exception, ValueError):
        # print("Error in Season and Trend estimation, in position:{0}".format(coord))
        # err_table['s&t_ext'] = 1
        # return err_table

    # Calculate season length and expected number of season
    try:
        season_lng = len(seasons) * param['dys_mlt']
        expSeason = int((ts_cleaned.index.max() - ts_cleaned.index.min()) / pd.Timedelta(season_lng, unit='d'))
    except (RuntimeError, Exception, ValueError):
        print("Error! Season conversion to days failed, in position:{0}".format(coord))
        err_table['s2d'] = 1
        return err_table

    if param['medspan'] == 0:
        medspan = season_lng / 7
    else:
        medspan = param['medspan']

    # Interpolate data to daily sample
    ts_d = ts_cleaned.resample('D').asfreq().interpolate(method='linear').fillna(0)

    try:
        if param['smp'] != 0:
            # Savinsky Golet filter
            ps = savgol_filter(ts_d, param['medspan'], param['smp'], mode='nearest')
            # TODO automatic selection of savgol window
            ps = pd.Series(ps, ts_d.index)
        else:
            ps = ts_d.rolling(medspan // 2 * 2, win_type='boxcar', center=True).mean()
    except (RuntimeError, Exception, ValueError):
        print('Error! Savinsky Golet filter problem, in position:{0}'.format(coord))
        err_table['savgol'] = 1
        return err_table

    # Valley detection
    # Detrending to catch better points
    vtrend = pd.Series(fit_trend(ps), index=ps.index)
    vdetr = ps - vtrend
    try:
        if 200.0 < season_lng < 400.0:
            mpd_val = int(season_lng * 2/3)
        elif season_lng < 200:
            mpd_val = int(season_lng * 1/3)
        else:
            mpd_val = int(season_lng * (param['tr'] - param['tr'] * 1 / 3) / 100)

        ind = dp.detect_peaks(vdetr, mph=vdetr.mean(),
                              mpd=mpd_val,
                              valley=True,

                              edge='both',
                              kpsh=False)
        oversmp = False

        if not ind.any():
            ind = dp.detect_peaks(vdetr, mph=-20,
                                  mpd=60,
                                  valley=True)
            oversmp = True

    except ValueError as e:
        print('Error in valley detection, in position:{0}, error {1}'.format(coord, e))
        err_table['vally'] = 1
        return err_table

    # Valley point time series conversion
    pks = pd.Series()
    for i in ind:
        pks[ps.index[i]] = ps.iloc[i]

    pks0 = pd.Series()
    for i in ind:
        pks0[vdetr.index[i]] = -vdetr.iloc[i]



























# #region
# def _filter_outlayer(self, cube):
#
#     # Filter outlier
#     try:
#         ts_clean = outlier.season(ts_dek_resc, 1, param['dys_x_yr'] * param['outMax'], 3)
#  # TODO make variable dek
#     except (RuntimeError, Exception, ValueError):
#         print('Outlier research failed, in position:{0}'.format(coord))
#         err_table['outlier'] = 1
#         return err_table
#
#     try:
#         if ts_clean is not None:
#             # ts_cleaned = outlier.fillerSeason(ts_clean)
#             ts_cleaned = ts_clean.interpolate()
#         else:
#             # No season strategy
#             summedYrs = ts_dek_resc.groupby(ts_dek_resc.index.year).sum()
#             season_lng = 9999
#
#             ts_table = {'sd': pd.Series(0, index=np.unique(ts_dek.index.year)),
#                         'ed': pd.Series(0, index=np.unique(ts_dek.index.year)),
#                         'sl': pd.Series(np.nan, index=np.unique(ts_dek.index.year)),
#                         'spi': pd.Series(np.nan, index=np.unique(ts_dek.index.year)),
#                         'si': summedYrs,
#                         'cf': pd.Series(0, index=np.unique(ts_dek.index.year)),
#                         'yr': pd.Series(np.unique(ts_dek.index.year), index=np.unique(ts_dek.index.year))}
#             return pd.DataFrame(ts_table), err_table, season_lng
#     except (RuntimeError, Exception, ValueError):
#         print('Cleaning failed, in position:{0}'.format(coord))
#         err_table['cleaning'] = 1
#         return err_table
#
# def _season_lenght(self):
#     # Estimate Season length
#     try:
#         seasons, trend = fit_seasons(ts_cleaned)
#     except (RuntimeError, Exception, ValueError):
#         print("Error in Season and Trend estimation, in position:{0}".format(coord))
#         err_table['s&t_ext'] = 1
#         return err_table
#     # Calculate season length and expected number of season
#     try:
#         season_lng = len(seasons) * param['dys_mlt']
#         expSeason = int((ts_cleaned.index.max() - ts_cleaned.index.min()) / pd.Timedelta(season_lng, unit='d'))
#     except (RuntimeError, Exception, ValueError):
#         print("Error! Season conversion to days failed, in position:{0}".format(coord))
#         err_table['s2d'] = 1
#         return err_table
#
#     if param['medspan'] == 0:
#         medspan = season_lng / 7
#     else:
#         medspan = param['medspan']
#
# def _interpolate_to_daily(self):
#     # Interpolate data to daily sample
#     ts_d = ts_cleaned.resample('D').asfreq().interpolate(method='linear').fillna(0)
#
#     try:
#         if param['smp'] != 0:
#             # Savinsky Golet filter
#             ps = savgol_filter(ts_d, param['medspan'], param['smp'], mode='nearest')
#             # TODO automatic selection of savgol window
#             ps = pd.Series(ps, ts_d.index)
#         else:
#             ps = ts_d.rolling(medspan // 2 * 2, win_type='boxcar', center=True).mean()
#     except (RuntimeError, Exception, ValueError):
#         print('Error! Savinsky Golet filter problem, in position:{0}'.format(coord))
#         err_table['savgol'] = 1
#         return err_table
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#         #region prosecuzione
#             # if hasattr(self, 'scale'):
#             #     cube = cube * self.scale
#             #
#             # if hasattr(self, 'offset'):
#             #     cube = cube * self.offset
#             #
#             # if hasattr(self, 'mask'):
#             #     cube = self._no_data(cube)
#             #
#             # if hasattr(self, 'min') and hasattr(self, 'max'):
#             #     cube = self._rescale(cube)
#             # else:
#             #     self.min = cube.min()
#             #     self.max = cube.max()
#             #     cube = self._rescale(cube)
#             #
#             # return cube
#         #endregion
#
