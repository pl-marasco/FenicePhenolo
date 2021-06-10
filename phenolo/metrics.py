# -*- coding: utf-8 -*-

import logging
import sys
import numpy as np
import pandas as pd
import scipy
from numba import jit, vectorize, float64
from phenolo import peaks, atoms
import xarray as xr

logger = logging.getLogger(__name__)
np.warnings.filterwarnings('ignore')


@jit(nopython=True, cache=True)
def rescale(ts, min, max):
    # Rescale values to  0-100
    #np.multiply(np.divide(np.subtract(ts, min), np.subtract(max, min)), 100)
    return ((ts - min) / (max - min)) * 100


@jit(nopython=True, cache=True)
def offset(ts,  offset, **kwargs):
    # add the offset
    try:
        return ts - offset
    except (RuntimeError,  Exception,  ValueError):
        print('Error in rescaling,  in position')
        logger.debug('Error in rescaling')
        sys.exit()


@jit(nopython=True, cache=True)
def scale(ts, scale):
    # scale according to the metadata
    # np.multiply(ts, scale)
    return ts * scale


def to_timeseries(values,  index):
    if len(values) != len(index):
        logger.debug('Lenght of the time series is different than the index provided')
        return ValueError

    return pd.Series(values,  index=index)


def scipy_valley_detection(pxldrl, tr):
    # Valley detection
    # Detrending to catch better points

    vtrend = pd.Series(pxldrl.trend_d,  index=pxldrl.ps.index)
    vdetr = pxldrl.ps - vtrend

    if 200.0 < pxldrl.season_lng < 400.0:
        mpd_val = int(pxldrl.season_lng * 2 / 3)
    elif pxldrl.season_lng < 200:
        mpd_val = int(pxldrl.season_lng * 1 / 3)
    else:
        mpd_val = int(pxldrl.season_lng * (tr - tr * 1 / 3) / 100)

    n_vdetr = np.negative(vdetr)

    ind, property = scipy.signal.find_peaks(n_vdetr, height=vdetr.mean(), distance=mpd_val)

    if not ind.any():
        ind, property = scipy.signal.find_peaks(n_vdetr, height=-20, distance=60)

    # Valley point time series conversion
    pks = pxldrl.ps.iloc[ind]

    return pks


def valley_detection(ps, trend_d, season_lng):
    # Valley detection
    # Detrending to catch better points

    # vtrend = pd.Series(trend_d,  index=ps.index)
    # vdetr = ps - vtrend
    # mph = vdetr.mean()

    vdetr = subtract(ps.values, trend_d.values)
    mph = vdetr.mean()

    if 200.0 < season_lng < 400.0:
        mpd_val = int(season_lng * 2 / 3)
    elif season_lng < 200:
        mpd_val = int(season_lng * 1 / 3)
    else:
        pass #TODO figureout what was tr
        #mpd_val = int(season_lng * (tr - tr * 1 / 3) / 100)

    ind = peaks.detect_peaks(vdetr,  mph=mph,
                             mpd=mpd_val, 
                             valley=True, 
                             edge='both', 
                             kpsh=False)
    if not ind.any():
        ind = peaks.detect_peaks(vdetr,  mph=-20,  mpd=60,  valley=True)

    # Valley point time series conversion
    pks = ps.iloc[ind]

    return pks


@jit(nopython=True, cache=True)
def subtract(val1, val2):
    return np.subtract(val1, val2)


def cycle_metrics(pks, ps, position):
    """
    Create an array of cycles with all the attributes populated

    :param pxldrl: a pixel drill object
    :return: an array of single cycles
    """

    sincys = []

    for i in range(pks.size - 1):

        # Minimum minimum time series
        sincy = atoms.SingularCycle(ps,  pks.index[i],  pks.index[i + 1])

        # avoid unusual results
        if sincy.ref_yr not in range(pks.index[i].year - 1, pks.index[i + 1].year + 1):
            logger.info(f'Warning! sbc not in a valid range,  in position:{position}')
            sincy.warn = 1  # 'sbc not in a valid range'
            continue

        sincys.append(sincy)

    return sincys


def attr_statistic(objects,  stat_type,  attribute):
    """
    Calculate a specific atrtibute stat_type over an array of objects

    :param objects: list of objects over withch must be calculated the statistics
    :param stat_type:
    :param attribute: atrtibute to be analysed
    :return:
    """

    value = None

    try:
        value = stat_type(filter(lambda x: x is not None,  [getattr(i,  attribute) for i in objects]))
    except ValueError:
        logger.debug('Statistic calculation has been unsuccessful.')
        ValueError('Statistic calculation has been unsuccessful.')

    try:
        value_d = pd.to_datetime(value,  unit='s')
    except ValueError:
        logger.debug('Date conversion of the stat_type calculation has been unsuccessful.')
        raise ValueError('Date conversion of the stat_type calculation has been unsuccessful.')
    return value_d


def __buffer_ext(sd, ed, mas, mms_b):
    """
    Add a buffer before and after the single cycle

    :param sincy:
    :return: buffered sincy
    """

    b_sd,  b_ed = None,  None
    if sd - mas >= mms_b.index[0]:
        b_sd = sd - mas
    if ed + mas <= mms_b.index[-1]:
        b_ed = ed + mas
    if sd or ed:
        return mms_b.loc[b_sd:b_ed]
    else:
        return mms_b

#@profile
def __back(smoothed, cbcd, sd, ed, delta_shift):
    """
    Calculate the curve shifted positively and truncated according to the delta and the starting date

    :param sincy: pandas time series
    :param delta_shift: pandas timedelta
    :return: pandas ts
    # """
    # shifted = smoothed.loc[:cbcd].shift(delta_shift.days, freq='d')
    # return shifted.loc[sd:].dropna()
    #shifted = smoothed.shift(delta_shift.days, freq='d')
    shifted = smoothed.shift(delta_shift.days)
    return shifted.loc[sd:ed]

#@profile
def __forward(smoothed, cbcd, sd, ed, delta_shift):
    """
    Calculate the curve shifted negatively and truncated according to the delta and the starting date

    :param sincy: pandas time series
    :param delta_shift: pandas timedelta
    :return: pandas ts
    """
    # shifted = smoothed.loc[cbcd:].shift(-delta_shift.days,  freq='d')
    # return shifted.loc[:ed].dropna()
    #shifted = smoothed.shift(-delta_shift.days,  freq='d')
    shifted = smoothed.shift(-delta_shift.days)
    return shifted.loc[sd:ed]


def __xarray_mean(data, width):
    return xr.DataArray(data, dims='x').rolling(x=width, min_periods=1, center=True).mean().to_pandas()


def phen_metrics(pxldrl,  param):
    """
    Calculate the Phenology parameter

    :param pxldrl: provide a pixel drill object from the module atoms
    :param param: provide a parameter object
    :return: list of sincy objects with added values
    """

    phen = []
    for sincy in pxldrl.sincys:

        if sincy.err:
            continue

        # specific mas
        sincy.mas = __mas(sincy.mml,  param.mavmet,  sincy.csdd)

        if sincy.mas.days < 0:
            sincy.mas = pd.to_timedelta(param.mavspan,  unit='D').days

        # buffer extractor
        # try:
        #     # sincy.buffered = __buffer_ext(sincy.sd, sincy.ed, sincy.mas, sincy.mms_b)
        #
        # except (RuntimeError,  Exception,  ValueError):
        #     logger.debug(f'Warning! Buffered curve not properly created,  in position:{pxldrl.position}')
        #     sincy.warn = 2  # 'Buffered curve'
        #     continue

        try:
            # sincy.smth_crv = sincy.buffered.rolling(sincy.mas.days, center=True).mean(numeric_only=True)
            #sincy.smth_crv = __xarray_mean(sincy.mms_b, sincy.mas.days)
            sincy.smth_crv = pd.Series(moving_average(sincy.mms_b.values, sincy.mas.days, center=True), index=sincy.mms_b.index)
        except (RuntimeError,  Exception,  ValueError):
            logger.debug(f'Warning! Smoothed curve calculation went wrong,  in position:{pxldrl.position}')
            sincy.warn = 3  # 'Smoothed curve'
            continue

        # sincy.smoothed = sincy.smth_crv.loc[sincy.sd - sincy.td:sincy.ed + sincy.td]
        sincy.smoothed = sincy.smth_crv

        # shift of the smoothed curve
        delta_shift = (sincy.mas / 2)

        # calculate the back curve
        sincy.back = __back(sincy.smoothed, sincy.cbcd, sincy.sd, sincy.ed, delta_shift)

        # calculate the forward curve
        sincy.forward = __forward(sincy.smoothed, sincy.cbcd, sincy.sd, sincy.ed,  delta_shift)

        # research the starting point of the season (SB)
        try:

            #sincy.intcpt_bk = __intercept(sincy.mms, sincy.back)

            sincy.intcpt_bk = _intercept(sincy.mms.values, sincy.back.values)
            sincy.sb = sincy.mms.iloc[sincy.intcpt_bk[0]]
            if sincy.sb.index > sincy.max_idx:
                raise Exception
        except (RuntimeError,  Exception,  ValueError):
            logger.debug(f'Warning! Start date not found in position {pxldrl.position} '
                         f'for the cycle starting in{sincy.sd}')
            sincy.warn = 4  # 'Start date'
            continue

        # research the end point of the season (SE)
        try:
            # sincy.intcpt_fw = __intercept(sincy.mms, sincy.forward)

            sincy.intcpt_fw = _intercept(sincy.mms.values, sincy.forward.values)
            sincy.se = sincy.mms.iloc[sincy.intcpt_fw[-1]]
            if sincy.se.index < sincy.max_idx:
                raise Exception
        except (RuntimeError,  Exception,  ValueError):
            logger.debug(f'Warning! End date not found in position {pxldrl.position} '
                         f'for the cycle starting in{sincy.sd}')
            sincy.warn = 5  # 'End date'
            continue

        # Season slope (SLOPE)
        try:
            sincy.sslp = ((sincy.se.values - sincy.sb.values) / (sincy.se.index.asi8 - sincy.sb.index.asi8)) * 100
        except ValueError:
            logger.debug(f'Warning! Error in slope calculation in pixel:{pxldrl.position}'
                         f'for the cycle starting in {sincy.sd}')
            pxldrl.warn = 6  # 'Slope'
            continue

        try:
            # Day of start in posix
            sincy.sbd = sincy.sb.index.asi8

            # Day of ends in posix
            sincy.sed = sincy.se.index.asi8

            # Season Length in days as ns
            sincy.sl = (sincy.se.index - sincy.sb.index).days.values

            # Season Integral [VX]
            # sincy.season = sincy.mms.loc[sincy.sb.index[0]:sincy.se.index[0]]  #TODO se.index == 0
            sincy.season = sincy.mms.iloc[sincy.intcpt_bk[0][0]:sincy.intcpt_fw[-1][0]]
            sincy.si = sincy.season.sum()

            # Season permanent
            sincy.sp = pd.Series(_min_min_line(sincy.season.index.asi8.copy(), sincy.season.values.copy()), sincy.season.index)
            # sincy.sp = sincy.season.copy()
            # sincy.sp[1:-1] = np.NaN
            # sincy.sp.interpolate(inplace=True)

            # Season permanent Integral [OX]
            #sincy.sp = sincy.sp.where(sincy.season.values - sincy.sp.values >= 0, sincy.season)
            sincy.sp = pd.Series(_sub(sincy.season.values, sincy.sp.values), index=sincy.sp.index)
            sincy.spi = sincy.sp.sum()

            # Cyclic fraction [VOX]
            sincy.cf = sincy.si - sincy.spi

            # Seasonal exceeding integral
            sincy.sei = sincy.stb - sincy.si

            # Active fraction
            # sincy.af = sincy.mms.iloc[sincy.intcpt_bk[0]:sincy.max_i] - sincy.sp[:sincy.max_i]
            # sincy.af = sincy.mms.loc[sincy.sb.index[0]:sincy.max_idx] - sincy.sp[:sincy.max_idx]
            # sincy.afi = sincy.af.sum()
            sincy.af = np.NaN
            sincy.afi = np.NaN

            # Reference yr
            sincy.ref_yr = (sincy.sb.index + pd.Timedelta(days=sincy.sl[0] * 0.66)).year

        except ValueError:
            sincy.sbd, sincy.sed, sincy.sl, sincy.sp, sincy.spi, \
                sincy.si, sincy.cf, sincy.af, sincy.afi, sincy.ref_yr, sincy.sei = [np.NaN] * 11
            sincy.warn = 100
            continue

        phen.append(sincy)

    return phen


def __mas(tsl,  mavmet,  sdd):
    """
    Calculate the mas over the single cycle.

    :param tsl: pandas time serie lenght
    :param mavmet: strenght of the equation [normally ~ 1.5-2]
    :param sdd: standard deviation expressed in yrs
    :return: mas ( moving avarage yearly)
    """
    mas = tsl - (2 * mavmet * sdd)
    if mas.days % 2 == 0:
        mas = mas + pd.to_timedelta(1,  unit='D')
    return mas
    # TODO to be reviewed


# def attribute_extractor(pxldrl, attribute, param):
#     try:
#         # values = list(
#         #     map(lambda phency:
#         #         {'index': phency.ref_yr.values[0],
#         #          'value': getattr(phency,  attribute)},  pxldrl.phen))
#         # if len(values) == 0:
#         #     raise Exception
#         #
#         # return pd.DataFrame(values).groupby('index').sum(numeric_only=True).reindex(param.dim_unq_val).squeeze()
#
#         values = {}
#         for phency in pxldrl.phen:
#             values[phency.ref_yr.values[0]] = getattr(phency,  attribute)
#
#         if len(values) == 0:
#             raise Exception
#
#         return pd.Series(values).groupby(level=0).sum(numeric_only=True).reindex(index=param.dim_unq_val, copy=False).values
#
#     except (RuntimeError,  Exception):
#         raise RuntimeError('Impossible to extract the attribute requested')
#
#
# def attribute_extractor_se(pxldrl, attribute, param):
#     try:
#         # values = list(
#         #     map(lambda phency:
#         #         {'index': phency.ref_yr.values[0],
#         #          'value': getattr(phency,  attribute)},  pxldrl.phen))
#         # if not values:
#         #     raise Exception
#         # return pd.DataFrame(values).groupby('index').min(numeric_only=True).reindex(param.dim_unq_val).squeeze()
#
#         values = {}
#         for phency in pxldrl.phen:
#             values[phency.ref_yr.values[0]] = getattr(phency,  attribute)
#
#         if len(values) == 0:
#             raise Exception
#
#         return pd.Series(values).groupby(level=0).min(numeric_only=True).reindex(index=param.dim_unq_val, copy=False).values
#
#     except (RuntimeError,  Exception):
#         raise RuntimeError('Impossible to extract the attribute requested')

def attribute_extractor(pxldrl, param):
    # try:
    #     values = dict((phency.ref_yr.values[0], (phency.stb,
    #                                            phency.mpi,
    #                                            phency.sbd[0],
    #                                            phency.sed[0],
    #                                            phency.sl[0],
    #                                            phency.spi,
    #                                            phency.si,
    #                                            phency.cf,
    #                                            phency.afi,
    #                                            phency.sei,
    #                                            phency.warn)) for phency in pxldrl.phen)
    #     df = pd.DataFrame.from_dict(values, orient='index', columns=['stb', 'mpi', 'sbd', 'sed',
    #                                                                'sl', 'spi', 'si', 'cf', 'afi', 'sei', 'warn']) \
    #         .groupby(level=0)
    #
    #     df_sum = df['sl', 'spi', 'si', 'cf', 'afi', 'sei', 'warn'].sum(numeric_only=True).reindex(
    #         index=param.dim_unq_val)
    #     df_min = df['stb', 'mpi', 'sbd', 'sed'].min(numeric_only=True).reindex(index=param.dim_unq_val)
    #     return df_min.join(df_sum, how="outer")
    #
    #     #todo add a way to keep warning in a non summed way

    try:
        values = [[phency.ref_yr.values[0],
                  phency.stb,
                  phency.mpi,
                  phency.sbd[0],
                  phency.sed[0],
                  phency.sl[0],
                  phency.spi,
                  phency.si,
                  phency.cf,
                  phency.afi,
                  phency.sei,
                  phency.warn] for phency in pxldrl.phen]
        df = pd.DataFrame(values, columns=['index', 'stb', 'mpi', 'sbd', 'sed', 'sl', 'spi',
                                           'si', 'cf', 'afi', 'sei', 'warn']).set_index('index').groupby(level=0)

        df_sum = df['sl', 'spi', 'si', 'cf', 'afi', 'sei', 'warn'].sum(numeric_only=True).reindex(
            index=param.dim_unq_val)
        df_min = df['stb', 'mpi', 'sbd', 'sed'].min(numeric_only=True).reindex(index=param.dim_unq_val)
        return df_min.join(df_sum, how="outer")

    except (RuntimeError,  Exception):
        raise RuntimeError('Impossible to extract the attribute requested')


def __intercept(reference_ts, shifted_ts):
    """
    Calculcate the intercept point
    :param s: Pandas TS
    :return:
    """
    s = (reference_ts - shifted_ts).values
    diff = np.diff(np.sign(s))
    return np.argwhere((diff != 0) & np.isfinite(diff))


@jit(nopython=True, cache=True)
def _min_min_line(x, y):
    m = (y[-1] - y[0]) / (x[-1] - x[0])
    c = y[0] - m * x[0]
    y[1:-1] = m*x[1:-1]+c
    return y


@jit(nopython=True, cache=True)
def _intercept(mms, shifted):
    delta = mms - shifted
    diff = np.diff(np.sign(delta))
    return np.argwhere((diff != 0) & np.isfinite(diff))


@jit(nopython=True, cache=True)
def _sub(val1, val2):
    # np.where(np.subtract(val1, val2) >= 0, val1, val2)
    return np.where(val1 - val2 >= 0, val2, val1)


@jit(nopython=True, cache=True)
def moving_average(array, window, center=False):
    if center:
        ret = np.cumsum(array)
        ret[window:] = ret[window:] - ret[:-window]
        ma = ret[window - 1:] / window
        n = np.empty(window//2); n.fill(np.nan)
        return np.concatenate((n.ravel(), ma.ravel(), n.ravel()))
    else:
        ret = np.cumsum(array)
        ret[window:] = ret[window:] - ret[:-window]
        ma = ret[window - 1:] / window
        n = np.empty(window-1); n.fill(np.nan)
        return np.concatenate((n.ravel(), ma.ravel(),n.ravel()))