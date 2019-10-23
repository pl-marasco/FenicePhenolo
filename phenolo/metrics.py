# -*- coding: utf-8 -*-

import logging
import sys

import numpy as np
import pandas as pd

from phenolo import peaks

logger = logging.getLogger(__name__)
np.warnings.filterwarnings('ignore')


def rescale(ts,  **kwargs):
    # Rescale values to  0-100
    param = kwargs.pop('settings',  '')
    try:
        return ((ts - param.min) / (param.max - param.min)) * 100
    except (RuntimeError,  Exception,  ValueError):
        print('Error in rescaling,  in position')
        logger.debug('Error in rescaling')
        sys.exit()


def offset(ts,  **kwargs):
    # add the offset
    param = kwargs.pop('param',  '')
    try:
        return ts - param.offset
    except (RuntimeError,  Exception,  ValueError):
        print('Error in rescaling,  in position')
        logger.debug('Error in rescaling')
        sys.exit()


def scale(ts,  **kwargs):
    # scale according to the metadata
    param = kwargs.pop('param',  '')
    try:
        return ts * param.scale
    except (RuntimeError,  Exception,  ValueError):
        print('Error in rescaling,  in position')
        logger.debug('Error in rescaling')
        sys.exit()


def to_timeseries(values,  index):
    if len(values) != len(index):
        logger.debug('Lenght of the time series is different than the index provided')
        return ValueError

    return pd.Series(values,  index=index)


def valley_detection(pxldrl,  param):
    # Valley detection
    # Detrending to catch better points

    vtrend = pd.Series(pxldrl.trend_d,  index=pxldrl.ps.index)
    vdetr = pxldrl.ps - vtrend

    if 200.0 < pxldrl.season_lng < 400.0:
        mpd_val = int(pxldrl.season_lng * 2 / 3)
    elif pxldrl.season_lng < 200:
        mpd_val = int(pxldrl.season_lng * 1 / 3)
    else:
        mpd_val = int(pxldrl.season_lng * (param.tr - param.tr * 1 / 3) / 100)

    ind = peaks.detect_peaks(vdetr,  mph=vdetr.mean(), 
                             mpd=mpd_val, 
                             valley=True, 
                             edge='both', 
                             kpsh=False)
    if not ind.any():
        ind = peaks.detect_peaks(vdetr,  mph=-20,  mpd=60,  valley=True)

    # Valley point time series conversion
    pks = pxldrl.ps.iloc[ind]

    return pks


def cycle_metrics(pxldrl):
    """
    Create an array of cycles with all the attributes populated

    :param pxldrl: a pixel drill object
    :return: an array of single cycles
    """

    sincys = []
    from phenolo import atoms

    for i in range(len(pxldrl.pks) - 1):

        # Minimum minimum time series
        sincy = atoms.SingularCycle(pxldrl.ps,  pxldrl.pks.index[i],  pxldrl.pks.index[i + 1])

        # avoid unusual results
        if sincy.ref_yr not in range(pxldrl.pks.index[i].year - 1,  pxldrl.pks.index[i + 1].year + 1):
            logger.info(f'Warning! sbc not in a valid range,  in position:{pxldrl.position}')
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


def __buffer_ext(sincy):
    """
    Add a buffer before and after the single cycle

    :param sincy:
    :return: buffered sincy
    """

    sd,  ed = None,  None
    if sincy.sd - sincy.mas >= sincy.mms_b.index[0]:
        sd = sincy.sd - sincy.mas
    if sincy.ed + sincy.mas <= sincy.mms_b.index[-1]:
        ed = sincy.ed + sincy.mas
    if sd or ed:
        return sincy.mms_b.loc[sd:ed]
    else:
        return sincy.mms_b


def __back(sincy,  delta_shift):
    """
    Calculate the curve shifted positively and truncated according to the delta and the starting date

    :param sincy: pandas time series
    :param delta_shift: pandas timedelta
    :return: pandas ts
    """
    shifted = sincy.smoothed.loc[:sincy.cbcd].shift(1,  freq=delta_shift)
    truncated = shifted.loc[sincy.sd:].dropna()
    return truncated


def __forward(sincy,  delta_shift):
    """
    Calculate the curve shifted negatively and truncated according to the delta and the starting date

    :param sincy: pandas time series
    :param delta_shift: pandas timedelta
    :return: pandas ts
    """
    shifted = sincy.smoothed.loc[sincy.cbcd:].shift(1,  freq=-delta_shift)
    truncated = shifted.loc[:sincy.ed].dropna()
    return truncated


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
            sincy.mas = pd.to_timedelta(param.mavspan,  unit='D')

        # buffer extractor
        try:
            sincy.buffered = __buffer_ext(sincy)
        except (RuntimeError,  Exception,  ValueError):
            logger.debug(f'Warning! Buffered curve not properly created,  in position:{pxldrl.position}')
            sincy.warn = 2  # 'Buffered curve'
            continue

        try:
            sincy.smth_crv = sincy.buffered.rolling(sincy.mas.days,  win_type='boxcar',  center=True) \
                .mean(numeric_only=True)
        except (RuntimeError,  Exception,  ValueError):
            logger.debug(f'Warning! Smoothed curve calculation went wrong,  in position:{pxldrl.position}')
            sincy.warn = 3  # 'Smoothed curve'
            continue

        sincy.smoothed = sincy.smth_crv.loc[sincy.sd - sincy.td:sincy.ed + sincy.td]

        # shift of the smoothed curve
        delta_shift = pd.Timedelta(days=int(sincy.mas.days / 2))

        # calculate the back curve
        sincy.back = __back(sincy,  delta_shift)

        # calculate the forward curve
        sincy.forward = __forward(sincy,  delta_shift)

        # research the starting point of the season (SB)
        try:
            sincy.intcpt_bk = __intercept((sincy.mms - sincy.back).values)
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
            sincy.intcpt_fw = __intercept((sincy.mms - sincy.forward).values)
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
            sincy.sslp = ((sincy.se.values - sincy.sb.values) / (sincy.se.index - sincy.sb.index).days) * 1e2
        except ValueError:
            logger.debug(f'Warning! Error in slope calculation in pixel:{pxldrl.position} '
                         f'for the cycle starting in {sincy.sd}')
            pxldrl.warn = 6  # 'Slope'
            continue

        try:
            # Day of start in posix
            sincy.sbd = sincy.sb.index

            # Day of ends in posix
            sincy.sed = sincy.se.index

            # Season Length in days as ns
            sincy.sl = (sincy.se.index - sincy.sb.index).to_pytimedelta()[0]

            # Season permanent
            sincy.sp = sincy.sb.append(sincy.se).resample('D').asfreq().interpolate(method='linear')
            # Season permanent Integral [OX]
            sincy.spi = sincy.sp.sum()

            # Season Integral [VX]
            sincy.si = sincy.mms.loc[sincy.sb.index[0]:sincy.se.index[0]].sum()

            # Cyclic fraction [VOX]
            sincy.cf = sincy.si - sincy.spi

            # Seasonal exceeding integral
            sincy.sei = sincy.stb - sincy.si

            # Active fraction
            sincy.af = sincy.mms.loc[sincy.sb.index[0]:sincy.max_idx] - sincy.sp[:sincy.max_idx]
            sincy.afi = sincy.af.sum()

            # Reference yr
            sincy.ref_yr = (sincy.sb.index + sincy.sl * 2 / 3).year

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
    return tsl - (2 * mavmet * sdd)
    # TODO to be reviewed


def attribute_extractor(pxldrl, dim_val_yrs, attribute):
    try:
        values = list(
            map(lambda phency:
                {'index': phency.ref_yr.values[0], 
                 'value': getattr(phency,  attribute)},  pxldrl.phen))
        if len(values) == 0:
            raise Exception

        return pd.DataFrame(values, index=dim_val_yrs).groupby('index').sum(numeric_only=True)

    except (RuntimeError,  Exception):
        raise RuntimeError('Impossible to extract the attribute requested')


def attribute_extractor_se(pxldrl,  attribute, param):
    try:
        values = list(
            map(lambda phency:
                dict([('index', phency.ref_yr.values[0]), (pxldrl.position[1],
                                                           getattr(phency,  attribute))]),  pxldrl.phen))
        if not values:
            raise Exception

        if attribute in ['stb', 'mpi',  'spi', 'si', 'cf', 'afi']:
            out = pd.DataFrame(values).groupby('index').sum(numeric_only=True).reindex(param.dim_val_yrs)
        else:  #['sbd', 'sed', 'warn']:
            out = pd.DataFrame(values).groupby('index').min(numeric_only=True).reindex(param.dim_val_yrs)
        return out

    except (RuntimeError,  Exception):
        raise RuntimeError('Impossible to extract the attribute requested')


def __intercept(s):
    """
    Calculcate the intercept point
    :param s: Pandas TS
    :return:
    """
    return np.argwhere((np.diff(np.sign(s)) != 0) & np.isfinite(np.diff(np.sign(s))))
