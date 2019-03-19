# -*- coding: utf-8 -*-

import nodata; import metrics; import outlier; import chronos; import filters
import logging
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
        return pxldrl

    # scaling
    try:
        if param.scale is not None:
            pxldrl.ts = metrics.scale(pxldrl.ts, settinngs=param)
    except(RuntimeError, ValueError, Exception):
        logger.info(f'Scaling error in position:{pxldrl.position}')
        pxldrl.error = True
        return pxldrl

    # off set
    try:
        if param.offset is not None:
            pxldrl.ts = metrics.offset(pxldrl.ts, settings=param)
    except(RuntimeError, ValueError, Exception):
        logger.info(f'Off set error in position:{pxldrl.position}')
        pxldrl.error = True
        return pxldrl

    # rescaling to 0-100
    try:
        if param.min is not None and param.max is not None:
            pxldrl.ts_resc = metrics.rescale(pxldrl.ts, settings=param)
    except(RuntimeError, ValueError, Exception):
        logger.info(f'Scaling error in position:{pxldrl.position}')
        pxldrl.error = True
        return pxldrl

    # Filter outlier
    try:
        pxldrl.ts_filtered = outlier.madseason(pxldrl.ts_resc, param.yr_dek, param.yr_dys * param.outmax, param.mad_pwr)  # TODO make variable dek
    except (RuntimeError, ValueError, Exception):
        logger.info(f'Error in filtering outlayer in position:{pxldrl.position}')
        pxldrl.error = True
        return pxldrl

    try:
        if pxldrl.ts_filtered is not None:
            pxldrl.ts_cleaned = pxldrl.ts_filtered.interpolate()
            # TODO make possible to use onother type of interpol
            if len(pxldrl.ts_cleaned[pxldrl.ts_cleaned.isna()]) > 0:
                pxldrl.ts_cleaned = pxldrl.ts_cleaned.fillna(method='bfill')

    except (RuntimeError, ValueError, Exception):
        logger.info(f'Error in interpolating outlayer in position:{pxldrl.position}')
        pxldrl.error = True
        return pxldrl

    # Estimate Season length
    try:
        pxldrl.seasons, pxldrl.trend = fit_seasons(pxldrl.ts_cleaned)
        # pxldrl.season_ts = metrics.to_timeseries(seasons, pxldrl.ps.index)
        pxldrl.trend_ts = metrics.to_timeseries(pxldrl.trend, pxldrl.ts_cleaned.index)
        # TODO add the no season option
         # pxldrl.season_lng
    except (RuntimeError, Exception, ValueError):
        logger.info(f'Error in estimate seaason and trend in position:{pxldrl.position}')
        pxldrl.error = True
        return pxldrl

    # Calculate season length and expected number of season
    try:
        pxldrl.season_lng = len(pxldrl.seasons) * param.yr_dys
        pxldrl.expSeason = chronos.season_ext(pxldrl)
    except(RuntimeError, Exception, ValueError):
        logger.info(f'Error! Season conversion to days failed, in position:{pxldrl.position}')
        pxldrl.error = True
        return pxldrl

    # medspan loading
    try:
        pxldrl.medspan = chronos.medspan(pxldrl.season_lng, param)
    except(RuntimeError, Exception, ValueError):
        logger.info(f'Error! Medspan calculation:{pxldrl.position}')
        pxldrl.error = True
        return pxldrl

    # Interpolate data to daily pxldrl
    try:
        pxldrl.ts_d = chronos.time_resample(pxldrl.ts_cleaned)
        pxldrl.trend_d = chronos.time_resample(pxldrl.trend_ts)
    except(RuntimeError, Exception, ValueError):
        logger.info(f'Error! Conversion to days failed, in position:{pxldrl.position}')
        pxldrl.error = True
        return pxldrl

    # Svainsky Golet
    try:
        pxldrl.ps = filters.sv(pxldrl, param)
    except (RuntimeError, Exception, ValueError):
        logger.info(f'Error! Savinsky Golet filter problem, in position:{pxldrl.position}')
        pxldrl.error = True
        return pxldrl

    # TODO controllare da qui in avanti il processo di analisi
    # Valley detection
    try:
        pxldrl.pks = metrics.valley_detection(pxldrl, param)
    except(RuntimeError, Exception, ValueError):
        logger.info(f'Error in valley detection in position:{pxldrl.position}')
        pxldrl.error = True
        return pxldrl

    # Cycle with matrics
    try:
        pxldrl.sincys = metrics.cycle_metrics(pxldrl, param)
    except(RuntimeError, Exception, ValueError):
        logger.info(f'Error in season detection in position:{pxldrl.position}')
        pxldrl.error = True
        return pxldrl

    try:
        import statistics
        pxldrl.msdd = metrics.attr_statistic(pxldrl.sincys, statistics.median, 'csd')
    except(RuntimeError, Exception, ValueError):
        logger.info(f'Error in season mean calculation in position:{pxldrl.position}')
        pxldrl.error = True
        return pxldrl

    # Season metrics
    try:
        pxldrl.phen = metrics.phen_metrics(pxldrl, param)
    except(RuntimeError, Exception, ValueError):
        logger.info(f'Error in intercept detection in position:{pxldrl.position}')
        pxldrl.error = True
        return pxldrl

    # General statistic agregation
    pxldrl.sl = metrics.attribute_extractor(pxldrl, 'sl')
    pxldrl.spi = metrics.attribute_extractor(pxldrl, 'spi')
    pxldrl.si = metrics.attribute_extractor(pxldrl, 'si')
    pxldrl.cf = metrics.attribute_extractor(pxldrl, 'cf')
    pxldrl.afi = metrics.attribute_extractor(pxldrl, 'afi')

    logger.debug(f'Pixel {pxldrl.position[0]}-{pxldrl.position[1]} processed')
    import sys
    sys.getsizeof(pxldrl)
    logger.debug(f'Pixel size {sys.getsizeof(pxldrl)}')

    return pxldrl
