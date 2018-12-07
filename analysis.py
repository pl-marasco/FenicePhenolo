# -*- coding: utf-8 -*-
# !/usr/bin/env python

import nodata; import metrics; import outlier; import chronos; import filters
import logging
from seasonal import fit_seasons


logger = logging.getLogger(__name__)


def phenolo(pxdrl, **kwargs):

    param = kwargs.pop('settings', '')

    # no data removing
    try:
        if param.sensor_typ == 'spot':
            pxdrl.ts = nodata.climate_fx(pxdrl.ts, settings=param)
    except(RuntimeError, ValueError, Exception):
        logger.info('Non data removal error in position:{0}'.format(pxdrl.position))
        pxdrl.error = True
        return pxdrl

    # scaling
    try:
        if param.scale is not None:
            pxdrl.ts = metrics.scale(pxdrl.ts, settinngs=param)
    except(RuntimeError, ValueError, Exception):
        logger.info('Scaling error in position:{0}'.format(pxdrl.position))
        pxdrl.error = True
        return pxdrl

    # off set
    try:
        if param.offset is not None:
            pxdrl.ts = metrics.offset(pxdrl.ts, settings=param)
    except(RuntimeError, ValueError, Exception):
        logger.info('off set error in position:{0}'.format(pxdrl.position))
        pxdrl.error = True
        return pxdrl

    # rescaling to 0-100
    try:
        if param.min is not None and param.max is not None:
            pxdrl.ts = metrics.rescale(pxdrl.ts, settings=param)
    except(RuntimeError, ValueError, Exception):
        logger.info('Scaling error in position:{0}'.format(pxdrl.position))
        pxdrl.error = True
        return pxdrl

    # Filter outlier
    try:
        pxdrl.ts_filtered = outlier.madseason(pxdrl.ts, param.yr_dek, param.yr_dys * param.outmax, 4)  # TODO make variable dek
    except (RuntimeError, ValueError, Exception):
        logger.info('Error in filtering outlayer in position:{0}'.format(pxdrl.position))
        pxdrl.error = True
        return pxdrl

    try:
        if pxdrl.ts_filtered is not None:
            pxdrl.ts_cleaned = pxdrl.ts_filtered.interpolate()  # TODO make possible to use onother type of interpol
    except (RuntimeError, ValueError, Exception):
        logger.info('Error in interpolating outlayer in position:{0}'.format(pxdrl.position))
        pxdrl.error = True
        return pxdrl

    # Estimate Season length
    try:
        pxdrl.seasons, pxdrl.trend = fit_seasons(pxdrl.ts_cleaned)
        # pxdrl.season_ts = metrics.to_timeseries(seasons, pxdrl.ps.index)
        pxdrl.trend_ts = metrics.to_timeseries(pxdrl.trend, pxdrl.ts_cleaned.index)

    except (RuntimeError, Exception, ValueError):
        logger.info('Error in estimate seaason and trend in position:{0}'.format(pxdrl.position))
        pxdrl.error = True
        return pxdrl

    # Calculate season length and expected number of season
    try:
        pxdrl.season_lng = len(pxdrl.seasons) * param.yr_dys
        pxdrl.expSeason = chronos.season_ext(pxdrl)
    except(RuntimeError, Exception, ValueError):
        logger.info('Error! Season conversion to days failed, in position:{0}'.format(pxdrl.position))
        pxdrl.error = True
        return pxdrl

    # medspan loading
    try:
        pxdrl.medspan = chronos.medspan(pxdrl.season_lng, param)
    except(RuntimeError, Exception, ValueError):
        logger.info('Error! medspan calculation:{0}'.format(pxdrl.position))
        pxdrl.error = True
        return pxdrl

    # Interpolate data to daily pxdrl
    try:
        pxdrl.ts_d = chronos.time_resample(pxdrl.ts_cleaned)
        pxdrl.trend_d = chronos.time_resample(pxdrl.trend_ts)
    except(RuntimeError, Exception, ValueError):
        logger.info('Error! Conversion to days failed, in position:{0}'.format(pxdrl.position))
        pxdrl.error = True
        return pxdrl

    # Svainsky Golet
    try:
        pxdrl.ps = filters.sv(pxdrl, param)
    except (RuntimeError, Exception, ValueError):
        logger.info('Error! Savinsky Golet filter problem, in position:{0}'.format(pxdrl.position))
        pxdrl.error = True
        return pxdrl

    # TODO controllare da qui in avanti il processo di analisi
    # Valley detection
    try:
        metrics.valley_detection(pxdrl, param)
    except(RuntimeError, Exception, ValueError):
        logger.info('Error in valley detection pixel position:{0}'.format(pxdrl.position))
        pxdrl.error = True
        return pxdrl

    # Season metrics
    try:
        metrics.seasons_metrics(pxdrl, param)
    except(RuntimeError, Exception, ValueError):
        logger.info('Error in season detection for pixel position:{0}'.format(pxdrl.position))
        pxdrl.error = True
        return pxdrl

    # Intercept
    try:
        metrics.phen_metrics(pxdrl, param)
    except(RuntimeError, Exception, ValueError):
        logger.info('Error in intercept detection in pixel position:{0}'.format(pxdrl.position))
        pxdrl.error = True
        return pxdrl

    return pxdrl
