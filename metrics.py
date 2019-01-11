import logging
import sys
import detect_peacks as dp
import pandas as pd; import numpy as np

logger = logging.getLogger(__name__)


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


def valley_detection(pxdrl, param):

    # Valley detection
    # Detrending to catch better points

    vtrend = pd.Series(pxdrl.trend_d, index=pxdrl.ps.index)
    vdetr = pxdrl.ps - vtrend

    if 200.0 < pxdrl.season_lng < 400.0:
        mpd_val = int(pxdrl.season_lng * 2 / 3)
    elif pxdrl.season_lng < 200:
        mpd_val = int(pxdrl.season_lng * 1 / 3)
    else:
        mpd_val = int(pxdrl.season_lng * (param.tr - param.tr * 1 / 3) / 100)

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

    # Valley point time series conversion
    pks = pd.Series()
    for i in ind:
        pks[pxdrl.ps.index[i]] = pxdrl.ps.iloc[i]

    # Points detrended
    pks0 = pd.Series()
    for i in ind:
        pks0[vdetr.index[i]] = -vdetr.iloc[i]

    return pks


def cycle_metrics(pxdrl, param):
    """
    Create an array of cycles with all the atrtibutes populated

    :param pxdrl: a pixel drill object
    :param param:
    :return: an array of single cycles
    """

    sincys = []

    for i in range(len(pxdrl.pks) - 1):

        import atoms

        # Minimum minimum time series
        sincy = atoms.SingularCycle(pxdrl.ps, pxdrl.pks.index[i], pxdrl.pks.index[i + 1])

        # avoid unusual results
        if sincy.ref_yr not in range(pxdrl.pks.index[i].year - 1, pxdrl.pks.index[i + 1].year + 1):
            logger.info(f'Warning! sbc not in a valid range, in position:{pxdrl.position}')
            pxdrl.error = True
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


def phen_metrics(pxdrl, param):
    """
    Calculate the Phenology paramter

    :param pxdrl: provide a pixel drill object from the module atoms
    :param param: provide a paramter object
    :return: list of sincy objects with added values
    """

    phen = []
    for sincy in pxdrl.sincys:

        # specific mas
        sincy.mas = _mas(sincy.mml, param.mavmet, sincy.csdd)
        # TODO verify the correctness of the standard devation

        if sincy.mas.days < 0:
            continue

        sincy.mmc = sincy.mms.copy()

        try:
            if sincy.sd-sincy.mas >= sincy.mms_b.index[0]:
                sd = sincy.sd-sincy.mas
            else:
                sd = sincy.mms_b.index[0]
            if sincy.ed-sincy.mas <= sincy.mms_b.index[-1]:
                ed = sincy.ed+sincy.mas
            else:
                ed = sincy.mms_b.index[-1]

            sincy.buffered = sincy.mms_b[sd:ed]
        except (RuntimeError, Exception, ValueError):
            logger.debug(f'Warning! Buffered curv not properly created, in position:{pxdrl.position}')
            pxdrl.error = True
            continue

        try:
            sincy.smth_crv = sincy.buffered.rolling(sincy.mas.days, win_type='boxcar', center=True).mean()
        except (RuntimeError, Exception, ValueError):
            logger.debug(f'Warning! Smoothed curv calculation went wrong, in position:{pxdrl.position}')
            pxdrl.error = True
            continue

        sincy.smoothed = sincy.smth_crv[sincy.sd - sincy.td:sincy.ed + sincy.td]

        # baricenter for the smoothed one
        try:
            posix_t_smth = sincy.smoothed.index.values.astype(np.int64) // 10 ** 9
            sincy.unx_sbc_Y_smth = (posix_t_smth * sincy.smoothed).sum() / posix_t_smth.sum()
        except (RuntimeError, ValueError, Exception):
            logger.debug(f'Warning! Baricenter not found in position {pxdrl.position} '
                         f'for the cycle starting in{sincy.sd}')
            continue

        # smoothed -= unx_sbc_Y_smth - unx_sbc_Y

        sincy.back = sincy.smoothed[:sincy.cbcd] \
                          .shift(1, freq=pd.Timedelta(days=int(sincy.mas.days / 2)))[sincy.sd:].dropna()

        sincy.forward = sincy.smoothed[sincy.cbcd:] \
                             .shift(1, freq=pd.Timedelta(days=-int(sincy.mas.days / 2)))[:sincy.ed].dropna()

        # sincy.mmc = sincy.mmc.reindex(pd.date_range(sincy.forward.index[0].date(),
        #                             sincy.back.index[len(sincy.forward) - 1].date()))

        sincy.sbd, sincy.sed, sincy.sbd_ts, sincy.sbd_ts = 4 * [None]

        # research the starting point of the season (SBD)
        try:
            sincy.ge_sbd = sincy.mmc.ge(sincy.back)
            change_sbd = sincy.ge_sbd.rolling(window=2, min_periods=2)\
                              .apply(lambda x: np.array_equal(x, [False, True]), raw=True)
            sincy.sbd = change_sbd[change_sbd == 1][:1] #[-1:]
            # TODO [proposal] add the possibility to select the interesection point
            sincy.sbd = pd.Series(sincy.mmc.loc[sincy.sbd.index], sincy.sbd.index)

        except (RuntimeError, Exception, ValueError):
            logger.debug(f'Warning! Start date not found in position {pxdrl.position} '
                         f'for the cycle starting in{sincy.sd}')
            continue

        # research the end point of the season (SED)
        try:
            sincy.ge_sed = sincy.mmc.ge(sincy.forward)
            change_sed = sincy.ge_sed.rolling(window=2, min_periods=2)\
                                     .apply(lambda x: np.array_equal(x, [True, False]), raw=True)
            # TODO [proposal] add the possibility to select the interesection point
            sincy.sed = change_sed[change_sed == 1][-1:]
            sincy.sed = pd.Series(sincy.mmc.loc[sincy.sed.index], sincy.sed.index)

        except (RuntimeError, Exception, ValueError):
            logger.debug(f'Warning! End date not found in position {pxdrl.position} '
                         f'for the cycle starting in{sincy.sd}')
            continue

        sincy.max = sincy.mmc[sincy.mmc == sincy.mmc.max()]
        sincy.max_date = sincy.mmc.idxmax()

        if sincy.sed.empty or sincy.sbd.empty:
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
                sincy.sslp = (sincy.sed.values[0] - sincy.sbd.values[0]) / \
                               (sincy.sed.index[0] - sincy.sbd.index[0]).days
            except ValueError:
                logger.debug(f'Warning! Error in slope calculation in pixel:{pxdrl.position} '
                             f'for the cycle starting in {sincy.sd}')
                continue

        try:
            # Season Lenght
            sincy.sl = (sincy.sed.index - sincy.sbd.index).to_pytimedelta()[0]

            # Season permanet
            sincy.sp = sincy.sbd.append(sincy.sed).resample('D').asfreq().interpolate(method='linear')
            # da pulire
            sincy.spi = sincy.sp.sum()

            # Season Integral
            sincy.si = sincy.mmc[sincy.sbd.index[0]:sincy.sed.index[0]].sum()

            # Cyclic fraction
            sincy.cf = sincy.si - sincy.spi

            # Active fraction
            sincy.af = sincy.mmc.loc[sincy.sbd.index[0]:sincy.max_date] - sincy.sp[:sincy.max_date]
            sincy.afi = sincy.af.sum()

            # reference yr
            sincy.ref_yr = (sincy.sbd.index + sincy.sl * 2 / 3).year

        except ValueError:
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


def _mas(ts, mavmet, sdd):
    """
    Calculate the mas over the single cycle.

    :param ts: pandas time serie
    :param mavmet: strenght of the equation [normally ~ 1.5-2]
    :param sdd: standard deviation expressed in yrs
    :return: mas ( moving avarage yearly)
    """
    return ts - (2 * mavmet * sdd)
    # TODO to be reviewed


def attribute_extractor(pxdrl, attribute):
    try:
        index = []
        for phency in pxdrl.phen:
            value = getattr(phency, attribute)
            index.append(pd.Series(value, phency.ref_yr))
        return pd.concat(index, axis=0)
    except RuntimeError:
        raise RuntimeError('Impossible to extract the attribute requested')
