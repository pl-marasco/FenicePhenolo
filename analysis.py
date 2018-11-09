# -*- coding: utf-8 -*-
# !/usr/bin/env python

import nodata; import metrics; import outlier; import chronos; import filters
import logging
from seasonal import fit_seasons


logger = logging.getLogger(__name__)


def phenolo(core, **kwargs):

    param = kwargs.pop('settings', '')

    # no data removing
    try:
        if param.sensor_typ == 'spot':
            core.ts = nodata.climate_fx(core.ts, settings=param)
    except(RuntimeError, ValueError, Exception):
        logger.info('Non data removal error in position:{0}'.format(core.position))
        core.error = True
        return core

    # scaling
    try:
        if param.scale is not None:
            core.ts = metrics.scale(core.ts, settinngs=param)
    except(RuntimeError, ValueError, Exception):
        logger.info('Scaling error in position:{0}'.format(core.position))
        core.error = True
        return core

    # off set
    try:
        if param.offset is not None:
            core.ts = metrics.offset(core.ts, settings=param)
    except(RuntimeError, ValueError, Exception):
        logger.info('off set error in position:{0}'.format(core.position))
        core.error = True
        return core

    # rescaling to 0-100
    try:
        if param.min is not None and param.max is not None:
            core.ts = metrics.rescale(core.ts, settings=param)
    except(RuntimeError, ValueError, Exception):
        logger.info('Scaling error in position:{0}'.format(core.position))
        core.error = True
        return core

    # Filter outlier
    try:
        core.ts_filtered = outlier.season(core.ts, 1, param.yr_dys * param.outmax, 3)  # TODO make variable dek
    except (RuntimeError, ValueError, Exception):
        logger.info('Error in filtering outlayer in position:{0}'.format(core.position))
        core.error = True
        return core

    try:
        if core.ts_filtered is not None:
            core.interpolated = core.ts_filtered.interpolate()  # TODO make possible to use onother type of interpolation
    except (RuntimeError, ValueError, Exception):
        logger.info('Error in interpolating outlayer in position:{0}'.format(core.position))
        core.error = True
        return core

    # Estimate Season length
    try:
        core.seasons, core.trend = fit_seasons(core.ts_cleaned)
    except (RuntimeError, Exception, ValueError):
        logger.info('Error in estimate seaason and trend in position:{0}'.format(core.position))
        core.error = True
        return core

    # Calculate season length and expected number of season
    try:
        core.season_lng = len(core.seasons) * param.yr_dys
        core.expSeason = chronos.season_ext(core)
    except(RuntimeError, Exception, ValueError):
        logger.info('Error! Season conversion to days failed, in position:{0}'.format(core.position))
        core.error = True
        return core

    # medspan loading
    try:
        chronos.medspan(core, param)
    except(RuntimeError, Exception, ValueError):
        logger.info('Error! medspan calculation:{0}'.format(core.position))
        core.error = True
        return core

    # Interpolate data to daily core
    try:
        chronos.time_resample(core)
    except(RuntimeError, Exception, ValueError):
        logger.info('Error! Conversion to days failed, in position:{0}'.format(core.position))
        core.error = True
        return core

    # Svainsky Golet
    try:
        core.ps = filters.sv(core, param)
    except (RuntimeError, Exception, ValueError):
        logger.info('Error! Savinsky Golet filter problem, in position:{0}'.format(core.position))
        core.error = True
        return core

    # TODO controllare da qui in avanti il processo di analisi
    # Valley detection
    try:
        metrics.valley_detection(core, param)
    except(RuntimeError, Exception, ValueError):
        logger.info('Error in valley detection pixel position:{0}'.format(core.position))
        core.error = True
        return core

    # Season metrics
    try:
        metrics.seasons_metrics(core, param)
    except(RuntimeError, Exception, ValueError):
        logger.info('Error in season detection for pixel position:{0}'.format(core.position))
        core.error = True
        return core

    # Intercept
    try:
        metrics.intercepts(core, param)
    except(RuntimeError, Exception, ValueError):
        logger.info('Error in intercept detection in pixel position:{0}'.format(core.position))
        core.error = True
        return core

    return core
