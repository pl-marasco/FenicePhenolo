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
        sincy = atoms.SingularCycle(pxdrl.ps[pxdrl.pks.index[i]:
                                    pxdrl.pks.index[i + 1]])

        # avoid unusual results
        if sincy.ref_yr not in range(pxdrl.pks.index[i].year - 1, pxdrl.pks.index[i + 1].year + 1):
            logger.info('Warning! sbc not in a valid range, in position:{0}'.format(pxdrl.position))
            pxdrl.error = True
            continue

        sincys.append(sincy)

    return sincys


def attr_statistic(objects, statistic, attribute):
    """
    Calculate a specific atrtibute statistic over an array of objects

    :param objects: list of objects over withch must be calculated the statistics
    :param attribute: atrtibute to be analysed
    :return:
    """

    value = None

    try:
        value = statistic(getattr(i, attribute) for i in objects)
    except ValueError:
        ValueError('Statistic calculation had been unsuccessful.')

    try:
        value_d = pd.to_datetime(value, unit='s')
    except ValueError:
        raise ValueError('Date conversion of the statistic calculation had been unsuccessful.')

    return value_d


def phen_metrics(pxdrl, param):

    """

    :param pxdrl:
    :param param:
    :return:
    """

    try:

        phen = []

        for i in range(len(pxdrl.seasons)):

            sincy = pxdrl.sincys[i]  # singular cycle

            # index = ts_table.iloc[i]['sbc'] # forse riferimento temporale

            # specific mas
            sincy.mas = sincy.mml - (2 * param.mavmet * sincy.csdd)
            # TODO verify the correctness of the standard devation
            # mas = (ts_table.iloc[i]['ed'] - ts_table.iloc[i]['sd']) - (2 * param['mavmet'] * ts_table.iloc[i]['ssd'])

            if sincy.mas.days < 0:
                continue

            sincy.mmc = sincy.mms.copy()

            sincy.timedelta = sincy.mml * (1/3)

            # this can cause the interseason problem
            sincy.buffered = sincy.mms[sincy.sd - sincy.timedelta: sincy.ed + sincy.timedelta]

            try:
                sincy.smth_crv = sincy.buffered.rolling(sincy.mas.days, win_type='boxcar', center=True).mean()
            except (RuntimeError, Exception, ValueError):
                logger.info('Warning! Smoothed curv calculation went wrong, in position:{0}'.format(pxdrl.position))
                pxdrl.error = True
                continue

            sincy.smoothed = sincy.smth_crv[sincy.sd - sincy.timedelta:sincy.ed + sincy.timedelta]

            # baricenter for the smoothed one
            posix_t_smth = sincy.smoothed.index.values.astype(np.int64) // 10 ** 9
            sincy.unx_sbc_Y_smth = (posix_t_smth * sincy.smoothed).sum() / posix_t_smth.sum()

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
                sincy.sbd = change_sbd[change_sbd == 1][-1:]
                sincy.sbd = pd.Series(sincy.mmc.loc[sincy.sbd.index], sincy.sbd.index)

            except (RuntimeError, Exception, ValueError):
                continue

            # research the end point of the season (SED)
            try:
                sincy.ge_sed = sincy.mmc.ge(sincy.forward)
                change_sed = sincy.ge_sed.rolling(window=2, min_periods=2)\
                                         .apply(lambda x: np.array_equal(x, [True, False]),raw=True)
                sincy.sed = change_sed[change_sed == 1][-1:]
                sincy.sed = pd.Series(sincy.mmc.loc[sincy.sed.index], sincy.sed.index)

            except (RuntimeError, Exception, ValueError):
                continue

            sincy.max_date = sincy.mmc.idxmax()

            if sincy.sed is None or sincy.sbd is None:
                continue
            else:
                # Season slope (SLOPE)
                sincy.sslp = (sincy.sed.values[0] - sincy.sbd.values[0]) / \
                               (sincy.sed.index[0] - sincy.sbd.index[0]).days
                if not (abs(sincy.sslp) < 0.15):
                    continue

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
            sincy.ref_yr = sincy.sbd.index + sincy.sl * 2 / 3

            phen.append(sincy)

        return phen

    except (RuntimeError, Exception, ValueError):
        return
