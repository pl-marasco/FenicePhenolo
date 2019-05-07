# -*- coding: utf-8 -*-

import logging

from phenolo import chronos, filters, metrics, nodata, outlier
from seasonal import fit_seasons

logger = logging.getLogger(__name__)


def phenolo(pxldrl, **kwargs):
    param = kwargs.pop('settings', '')

    # no data removing
    try:
        if param.sensor_typ == 'spot':
            pxldrl.ts = nodata.climate_fx(pxldrl.ts_raw, settings=param)
    except(RuntimeError, ValueError, Exception):
        logger.info(f'Nodata removal error in position:{pxldrl.position}')
        pxldrl.error = True
        pxldrl.errtyp = 'No data'
        return pxldrl

    # scaling
    try:
        if param.scale is not None:
            pxldrl.ts = metrics.scale(pxldrl.ts, settinngs=param)
    except(RuntimeError, ValueError, Exception):
        logger.info(f'Scaling error in position:{pxldrl.position}')
        pxldrl.error = True
        pxldrl.errtyp = 'Scaling'
        return pxldrl

    # off set
    try:
        if param.offset is not None:
            pxldrl.ts = metrics.offset(pxldrl.ts, settings=param)
    except(RuntimeError, ValueError, Exception):
        logger.info(f'Off set error in position:{pxldrl.position}')
        pxldrl.error = True
        pxldrl.errtyp = 'Off set'
        return pxldrl

    # scaling to 0-100
    try:
        if param.min is not None and param.max is not None:
            pxldrl.ts_resc = metrics.rescale(pxldrl.ts, settings=param)
    except(RuntimeError, ValueError, Exception):
        logger.info(f'Scaling error in position:{pxldrl.position}')
        pxldrl.error = True
        pxldrl.errtyp = '0-100 Scaling'
        return pxldrl

    # Filter outlier
    try:
        if pxldrl.ts_resc.isnull().sum() > 0:
            pxldrl.ts_resc.fillna(method='bfill', inplace=True)

        pxldrl.ts_filtered = outlier.doubleMAD(pxldrl.ts_resc, param.mad_pwr)

    except (RuntimeError, ValueError, Exception):
        logger.info(f'Error in filtering outlier in position:{pxldrl.position}')
        pxldrl.error = True
        pxldrl.errtyp = 'outlier filtering'
        return pxldrl

    try:
        if pxldrl.ts_filtered is not None:
            pxldrl.ts_cleaned = pxldrl.ts_filtered.interpolate()
            # TODO make possible to use onother type of interpol
            if len(pxldrl.ts_cleaned[pxldrl.ts_cleaned.isna()]) > 0:
                pxldrl.ts_cleaned = pxldrl.ts_cleaned.fillna(method='bfill')

    except (RuntimeError, ValueError, Exception):
        logger.info(f'Error in interpolating outlier in position:{pxldrl.position}')
        pxldrl.error = True
        pxldrl.errtyp = 'gap filling'
        return pxldrl

    # Estimate Season length
    try:
        pxldrl.seasons, pxldrl.trend = fit_seasons(pxldrl.ts_cleaned)
        if pxldrl.seasons is not None and pxldrl.trend is not None:
            pxldrl.trend_ts = metrics.to_timeseries(pxldrl.trend, pxldrl.ts_cleaned.index)
        else:
            raise Exception
        # TODO add the no season option
        # pxldrl.season_lng
    except (RuntimeError, Exception, ValueError):
        logger.info(f'Error in estimate season and trend in position:{pxldrl.position}')
        pxldrl.error = True
        pxldrl.errtyp = 'Season and trend estimation'
        pxldrl.season_lng = 0
        return pxldrl

    # Calculate season length and expected number of season
    try:
        pxldrl.season_lng = len(pxldrl.seasons) * param.yr_dys
        pxldrl.expSeason = chronos.season_ext(pxldrl)
    except(RuntimeError, Exception, ValueError):
        logger.info(f'Error! Season conversion to days failed, in position:{pxldrl.position}')
        pxldrl.error = True
        pxldrl.errtyp = 'To daily conversion'
        return pxldrl

    # medspan loading
    try:
        pxldrl.medspan = chronos.medspan(pxldrl.season_lng, param)
    except(RuntimeError, Exception, ValueError):
        logger.info(f'Error! Medspan calculation:{pxldrl.position}')
        pxldrl.error = True
        pxldrl.errtyp = 'madspan error'
        return pxldrl

    # Interpolate data to daily pxldrl
    try:
        pxldrl.ts_d = chronos.time_resample(pxldrl.ts_cleaned)
        pxldrl.trend_d = chronos.time_resample(pxldrl.trend_ts)
    except(RuntimeError, Exception, ValueError):
        logger.info(f'Error! Conversion to days failed, in position:{pxldrl.position}')
        pxldrl.error = True
        pxldrl.errtyp = 'Trend conversion to daily'
        return pxldrl

    # Svainsky Golet
    try:
        pxldrl.ps = filters.sv(pxldrl, param)
    except (RuntimeError, Exception, ValueError):
        logger.info(f'Error! Savinsky Golet filter problem, in position:{pxldrl.position}')
        pxldrl.error = True
        pxldrl.errtyp = 'Savinsky Golet'
        return pxldrl

    # TODO controllare da qui in avanti il processo di analisi
    # Valley detection
    try:
        pxldrl.pks = metrics.valley_detection(pxldrl, param)
    except(RuntimeError, Exception, ValueError):
        logger.info(f'Error in valley detection in position:{pxldrl.position}')
        pxldrl.error = True
        pxldrl.errtyp = 'Valley detection'
        return pxldrl

    # Cycle with matrics
    try:
        pxldrl.sincys = metrics.cycle_metrics(pxldrl)
    except(RuntimeError, Exception, ValueError):
        logger.info(f'Error in season detection in position:{pxldrl.position}')
        pxldrl.error = True
        pxldrl.errtyp = 'Season detection'
        return pxldrl

    try:
        import statistics
        pxldrl.msdd = metrics.attr_statistic(pxldrl.sincys, statistics.median, 'csd')
    except(RuntimeError, Exception, ValueError):
        logger.info(f'Error in season mean calculation in position:{pxldrl.position}')
        pxldrl.error = True
        pxldrl.errtyp = 'Season mean'
        return pxldrl

    # Season metrics
    try:
        pxldrl.phen = metrics.phen_metrics(pxldrl, param)
    except(RuntimeError, Exception, ValueError):
        logger.info(f'Error in intercept detection in position:{pxldrl.position}')
        pxldrl.error = True
        pxldrl.errtyp = 'Season metrics'
        return pxldrl

    # General statistic agregation
    try:
        pxldrl.sb = metrics.attribute_extractor_se(pxldrl, 'sb')
        pxldrl.se = metrics.attribute_extractor_se(pxldrl, 'se')
        pxldrl.sl = metrics.attribute_extractor(pxldrl, 'sl')
        pxldrl.spi = metrics.attribute_extractor(pxldrl, 'spi')
        pxldrl.si = metrics.attribute_extractor(pxldrl, 'si')
        pxldrl.cf = metrics.attribute_extractor(pxldrl, 'cf')
        pxldrl.afi = metrics.attribute_extractor(pxldrl, 'afi')
    except(RuntimeError, Exception, ValueError):
        logger.info(f'Statistical aggregation:{pxldrl.position}')
        pxldrl.error = True
        pxldrl.errtyp = 'Statistical aggregation'
        return pxldrl

    logger.debug(f'Pixel {pxldrl.position[0]}-{pxldrl.position[1]} processed')

    return pxldrl
