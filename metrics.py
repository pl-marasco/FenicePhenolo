# -*- coding: utf-8 -*-

import logging
import sys
import detect_peacks as dp
import pandas as pd; import numpy as np

logger = logging.getLogger(__name__)
np.warnings.filterwarnings('ignore')


def rescale(ts, **kwargs):
    # Rescale values to  0-100
    param = kwargs.pop('settings', '')
    try:
        return ((ts - param.min) / (param.max - param.min)) * 100
    except (RuntimeError, Exception, ValueError):
        print('Error in rescaling, in position')
        logger.debug('Error in rescaling')
        sys.exit()


def offset(ts, **kwargs):
    # add the offset
    param = kwargs.pop('param', '')
    try:
        return ts - param.offset
    except (RuntimeError, Exception, ValueError):
        print('Error in rescaling, in position')
        logger.debug('Error in rescaling')
        sys.exit()


def scale(ts, **kwargs):
    # scale according to the metadata
    param = kwargs.pop('param', '')
    try:
        return ts * param.scale
    except (RuntimeError, Exception, ValueError):
        print('Error in rescaling, in position')
        logger.debug('Error in rescaling')
        sys.exit()


def to_timeseries(values, index):
    if len(values) != len(index):
        logger.debug('Lenght of the time series is different than the index provided')
        return ValueError

    return pd.Series(values, index=index)


def valley_detection(pxldrl, param):

    # Valley detection
    # Detrending to catch better points

    vtrend = pd.Series(pxldrl.trend_d, index=pxldrl.ps.index)
    vdetr = pxldrl.ps - vtrend

    if 200.0 < pxldrl.season_lng < 400.0:
        mpd_val = int(pxldrl.season_lng * 2 / 3)
    elif pxldrl.season_lng < 200:
        mpd_val = int(pxldrl.season_lng * 1 / 3)
    else:
        mpd_val = int(pxldrl.season_lng * (param.tr - param.tr * 1 / 3) / 100)

    ind = dp.detect_peaks(vdetr, mph=vdetr.mean(),
                          mpd=mpd_val,
                          valley=True,
                          edge='both',
                          kpsh=False)
    # oversmp = False

    if not ind.any():
        ind = dp.detect_peaks(vdetr, mph=-20, mpd=60, valley=True)
        # oversmp = True

    # Valley point time series conversion
    pks = pd.Series()
    for i in ind:
        pks[pxldrl.ps.index[i]] = pxldrl.ps.iloc[i]

    # Points detrended
    pks0 = pd.Series()
    for i in ind:
        pks0[vdetr.index[i]] = -vdetr.iloc[i]

    return pks


def cycle_metrics(pxldrl):
    """
    Create an array of cycles with all the attributes populated

    :param pxldrl: a pixel drill object
    :return: an array of single cycles
    """

    sincys = []
    import atoms

    for i in range(len(pxldrl.pks) - 1):

        # Minimum minimum time series
        sincy = atoms.SingularCycle(pxldrl.ps, pxldrl.pks.index[i], pxldrl.pks.index[i + 1])

        # avoid unusual results
        if sincy.ref_yr not in range(pxldrl.pks.index[i].year - 1, pxldrl.pks.index[i + 1].year + 1):
            logger.info(f'Warning! sbc not in a valid range, in position:{pxldrl.position}')
            pxldrl.error = True
            pxldrl.errtyp = 'sbc not in a valid range'
            continue

        sincys.append(sincy)

    return sincys


def attr_statistic(objects, stat_type, attribute):
    """
    Calculate a specific atrtibute stat_type over an array of objects

    :param objects: list of objects over withch must be calculated the statistics
    :param stat_type:
    :param attribute: atrtibute to be analysed
    :return:
    """

    value = None

    try:
        value = stat_type(filter(lambda x: x is not None, [getattr(i, attribute) for i in objects]))
    except ValueError:
        logger.debug('Statistic calculation has been unsuccessful.')
        ValueError('Statistic calculation has been unsuccessful.')

    try:
        value_d = pd.to_datetime(value, unit='s')
    except ValueError:
        logger.debug('Date conversion of the stat_type calculation has been unsuccessful.')
        raise ValueError('Date conversion of the stat_type calculation has been unsuccessful.')
    return value_d


def phen_metrics(pxldrl, param):
    """
    Calculate the Phenology paramter

    :param pxldrl: provide a pixel drill object from the module atoms
    :param param: provide a paramter object
    :return: list of sincy objects with added values
    """

    phen = []
    for sincy in pxldrl.sincys:

        if sincy.err is True:
            continue

        # specific mas
        sincy.mas = _mas(sincy.mml, param.mavmet, sincy.csdd)
        # TODO verify the correctness of the standard deviation

        if sincy.mas.days < 0:
            continue

        # Maximum point and date
        # sincy.max = sincy.mms[sincy.mms == sincy.mms.max()]
        sincy.max = sincy.mms.loc[[sincy.mms.idxmax()]]

        # TODO verifying buffer use and indexing
        try:
            if sincy.sd-sincy.mas >= sincy.mms_b.index[0]:
                sd = sincy.sd-sincy.mas
            else:
                sd = sincy.mms_b.index[0]
            if sincy.ed-sincy.mas <= sincy.mms_b.index[-1]:
                ed = sincy.ed+sincy.mas
            else:
                ed = sincy.mms_b.index[-1]

            sincy.buffered = sincy.mms_b.loc[sd:ed]

        except (RuntimeError, Exception, ValueError):
            logger.debug(f'Warning! Buffered curv not properly created, in position:{pxldrl.position}')
            pxldrl.error = True
            pxldrl.errtyp = 'Bff crv'
            continue

        try:
            sincy.smth_crv = sincy.buffered.rolling(sincy.mas.days, win_type='boxcar', center=True)\
                .mean(numeric_only=True)
        except (RuntimeError, Exception, ValueError):
            logger.debug(f'Warning! Smoothed curv calculation went wrong, in position:{pxldrl.position}')
            pxldrl.error = True
            pxldrl.errtyp = 'Smth crv'
            continue

        sincy.smoothed = sincy.smth_crv.loc[sincy.sd - sincy.td:sincy.ed + sincy.td]

        # baricenter for the smoothed one
        try:
            posix_t_smth = sincy.smoothed.index.values.astype(np.int64) // 10 ** 9
            sincy.unx_sbc_Y_smth = (posix_t_smth * sincy.smoothed).sum() / posix_t_smth.sum()
        except (RuntimeError, ValueError, Exception):
            logger.debug(f'Warning! Baricenter not found in position {pxldrl.position} '
                         f'for the cycle starting in{sincy.sd}')
            pxldrl.errtyp = 'Baricenter'
            continue

        # shift of the smoothed curve
        delta = pd.Timedelta(days=int(sincy.mas.days / 2))

        # calculate the back curve
        sincy.back = sincy.smoothed.loc[:sincy.cbcd].shift(1, freq=delta).loc[sincy.sd:].dropna()
        # calculate the forward curve
        sincy.forward = sincy.smoothed.loc[sincy.cbcd:].shift(1, freq=-delta).loc[:sincy.ed].dropna()

        sincy.sbd, sincy.sed, sincy.sbd_ts, sincy.sbd_ts = 4 * [None]

        # research the starting point of the season (SBD)
        try:
            sincy.intcpt_bk = intercept((sincy.mms - sincy.back).values)
            sincy.sbd = (sincy.mms.iloc[sincy.intcpt_bk[0]])
            if sincy.sbd.index > sincy.max.index:
                raise Exception

        except (RuntimeError, Exception, ValueError):
            logger.debug(f'Warning! Start date not found in position {pxldrl.position} '
                         f'for the cycle starting in{sincy.sd}')
            sincy.sbd = None
            pxldrl.errtyp = 'Start date'

        # research the end point of the season (SED)
        try:
            sincy.intcpt_fw = intercept((sincy.mms - sincy.forward).values)
            sincy.sed = (sincy.mms.iloc[sincy.intcpt_fw[-1]])
            if sincy.sed.index < sincy.max.index:
                raise Exception

        except (RuntimeError, Exception, ValueError):
            logger.debug(f'Warning! End date not found in position {pxldrl.position} '
                         f'for the cycle starting in{sincy.sd}')
            sincy.sed = None
            pxldrl.errtyp = 'End date'

        if sincy.sed is None or sincy.sbd is None:
            sincy.sl = np.NaN
            sincy.sp = np.NaN
            sincy.spi = np.NaN
            sincy.si = np.NaN
            sincy.cf = np.NaN
            sincy.af = np.NaN
            sincy.afi = np.NaN
            sincy.ref_yr = np.NaN
            continue
        else:
            # Season slope (SLOPE)
            try:
                sincy.sslp = ((sincy.sed.values - sincy.sbd.values) / (sincy.sed.index - sincy.sbd.index).days)*1e2
            except ValueError:
                logger.debug(f'Warning! Error in slope calculation in pixel:{pxldrl.position} '
                             f'for the cycle starting in {sincy.sd}')
                pxldrl.errtyp = 'Slope'
                continue

        try:
            # week of start
            sincy.sb = sincy.sbd.index.dayofyear

            # Week of ends
            sincy.se = sincy.sed.index.dayofyear

            # Season Lenght
            sincy.sl = (sincy.sed.index - sincy.sbd.index).to_pytimedelta()[0]

            # Season permanet
            sincy.sp = sincy.sbd.append(sincy.sed).resample('D').asfreq().interpolate(method='linear')
            # da pulire
            sincy.spi = sincy.sp.sum()

            # Season Integral
            sincy.si = sincy.mms.loc[sincy.sbd.index[0]:sincy.sed.index[0]].sum()

            # Cyclic fraction
            sincy.cf = sincy.si - sincy.spi

            # Active fraction
            sincy.af = sincy.mms.loc[sincy.sbd.index[0]:sincy.max.index[0]] - sincy.sp[:sincy.max.index[0]]
            sincy.afi = sincy.af.sum()

            # reference yr
            sincy.ref_yr = (sincy.sbd.index + sincy.sl * 2 / 3).year

        except ValueError:
            sincy.sb = np.NaN
            sincy.se = np.NaN
            sincy.sl = np.NaN
            sincy.sp = np.NaN
            sincy.spi = np.NaN
            sincy.si = np.NaN
            sincy.cf = np.NaN
            sincy.af = np.NaN
            sincy.afi = np.NaN
            sincy.ref_yr = np.NaN
            continue

        phen.append(sincy)

    return phen


def _mas(tsl, mavmet, sdd):
    """
    Calculate the mas over the single cycle.

    :param tsl: pandas time serie lenght
    :param mavmet: strenght of the equation [normally ~ 1.5-2]
    :param sdd: standard deviation expressed in yrs
    :return: mas ( moving avarage yearly)
    """
    return tsl - (2 * mavmet * sdd)
    # TODO to be reviewed


def attribute_extractor(pxldrl, attribute):
    try:
        values = list(
            map(lambda phency:
                {'index': phency.ref_yr.values[0],
                 'value': getattr(phency, attribute)}, pxldrl.phen))

        return pd.DataFrame(values).groupby('index').sum(numeric_only=True).squeeze()

    except RuntimeError:
        raise RuntimeError('Impossible to extract the attribute requested')


def attribute_extractor_se(pxldrl, attribute):
    try:
        values = list(
            map(lambda phency:
                {'index': phency.ref_yr.values[0],
                 'value': getattr(phency, attribute)}, pxldrl.phen))

        return pd.DataFrame(values).groupby('index').min(numeric_only=True).squeeze()

    except RuntimeError:
        raise RuntimeError('Impossible to extract the attribute requested')


def intercept(s):
    return np.argwhere((np.diff(np.sign(s)) != 0) & np.isfinite(np.diff(np.sign(s))))
