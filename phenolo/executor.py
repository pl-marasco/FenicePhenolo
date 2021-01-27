# -*- coding: utf-8 -*-

import logging
import numpy as np
import pandas as pd
import xarray as xr
import copy

from phenolo import atoms, analysis

logger = logging.getLogger(__name__)


def _cache_da(cube, param):
    years = param.dim_unq_val

    stb = np.empty((years.size, cube[param.row_nm].size, cube[param.col_nm].size), dtype=np.float)
    mpi = np.empty((years.size, cube[param.row_nm].size, cube[param.col_nm].size), dtype=np.float)
    sbd = np.empty((years.size, cube[param.row_nm].size, cube[param.col_nm].size), dtype=np.float)
    sed = np.empty((years.size, cube[param.row_nm].size, cube[param.col_nm].size), dtype=np.float)
    sl = np.empty((years.size, cube[param.row_nm].size, cube[param.col_nm].size), dtype=np.float)
    spi = np.empty((years.size, cube[param.row_nm].size, cube[param.col_nm].size), dtype=np.float)
    si = np.empty((years.size, cube[param.row_nm].size, cube[param.col_nm].size), dtype=np.float)
    cf = np.empty((years.size, cube[param.row_nm].size, cube[param.col_nm].size), dtype=np.float)
    afi = np.empty((years.size, cube[param.row_nm].size, cube[param.col_nm].size), dtype=np.float)
    warn = np.empty((years.size, cube[param.row_nm].size, cube[param.col_nm].size), dtype=np.float)
    n_seasons = np.empty((years.size, cube[param.row_nm].size, cube[param.col_nm].size), dtype=np.int64)
    err = np.empty((years.size, cube[param.row_nm].size, cube[param.col_nm].size), dtype=np.int64)

    return xr.Dataset({'Standing_Biomass': (['time', 'lat', 'lon'], stb),
                       'Min_min_PermanentIntegral': (['time', 'lat', 'lon'], mpi),
                       'Season_Start_date': (['time', 'lat', 'lon'], sbd),
                       'Season_End_date': (['time', 'lat', 'lon'], sed),
                       'Season_Lenght': (['time', 'lat', 'lon'], sl),
                       'Seasonal_Permanent_Integral': (['time', 'lat', 'lon'], spi),
                       'Season_Integral': (['time', 'lat', 'lon'], si),
                       'Cyclic_Fraction': (['time', 'lat', 'lon'], cf),
                       'Active_Fraction_Integral': (['time', 'lat', 'lon'], afi),
                       'Cycle_Warning': (['time', 'lat', 'lon'], warn),
                       'Number_of_Seasons': (['time', 'lat', 'lon'], n_seasons),
                       'Pixel_Critical_Error': (['time', 'lat', 'lon'], err)},
                        coords={param.dim_nm: years,
                                param.row_nm: cube[param.row_nm],
                                param.col_nm: cube[param.col_nm]})


def _process(cube, **kwargs):

    param = kwargs.pop('param', '')

    cache = _cache_da(cube, param)

    for i_row in range(0, cube[param.row_nm].size):
        for i_col in range(0, cube[param.col_nm].size):
            pxldrl = atoms.PixelDrill(cube.isel(
                                                dict([(param.col_nm, i_col),
                                                      (param.row_nm, i_row)])).to_series().astype(float),
                                      [i_row, i_col])

            try:
                pxldrl = analysis.phenolo(pxldrl, settings=param)
            except Exception as e:
                cache['Pixel_Critical_Error'][:, i_row, i_col] = 1
                cache['Number_of_Seasons'][:, i_row, i_col] = -1
                logger.debug(f'Error: {_error_decoder(pxldrl.errtyp)} in position:{pxldrl.position}')

            try:
                if pxldrl.error:
                    cache['Pixel_Critical_Error'][:, i_row, i_col] = 2
                    cache['Number_of_Seasons'][:, i_row, i_col] = -1
                    logger.debug(f'Error: {_error_decoder(pxldrl.errtyp)} in position:{pxldrl.position}')
                else:
                    try:
                        cache['Standing_Biomass'][:, i_row, i_col] = copy.deepcopy(pxldrl.stb)
                        cache['Min_min_PermanentIntegral'][:, i_row, i_col] = copy.deepcopy(pxldrl.mpi)
                        cache['Season_Start_date'][:, i_row, i_col] = copy.deepcopy(pxldrl.sbd)
                        cache['Season_End_date'][:, i_row, i_col] = copy.deepcopy(pxldrl.sed)
                        cache['Season_Lenght'][:, i_row, i_col] = copy.deepcopy(pxldrl.sl)
                        cache['Seasonal_Permanent_Integral'][:, i_row, i_col] = copy.deepcopy(pxldrl.spi)
                        cache['Season_Integral'][:, i_row, i_col] = copy.deepcopy(pxldrl.si)
                        cache['Cyclic_Fraction'][:, i_row, i_col] = copy.deepcopy(pxldrl.cf)
                        cache['Active_Fraction_Integral'][:, i_row, i_col] = copy.deepcopy(pxldrl.afi)
                        cache['Cycle_Warning'][:, i_row, i_col] = copy.deepcopy(pxldrl.warn)

                        if not pd.isnull(pxldrl.season_lng):
                            if pxldrl.season_lng <= 365.0:
                                cache['Number_of_Seasons'][:, i_row, i_col] = int(365 / copy.deepcopy(pxldrl.season_lng))
                            else:
                                cache['Number_of_Seasons'][:, i_row, i_col] = int(copy.deepcopy(pxldrl.season_lng))
                                # TODO add more seasons sub division
                        else:
                            cache['Number_of_Seasons'][:, i_row, i_col] = copy.deepcopy(pxldrl.season_lng)

                        del pxldrl

                        # if param.ovr_scratch:
                        #     try:
                        #         import phenolo.output as output
                        #         output.scratch_dump(pxldrl, param)
                        #     except Exception as e:
                        #         raise Exception

                    except Exception as e:
                        logger.error(f'Error in a worker during the filling of type {e}')
                        cache['Pixel_Critical_Error'][:, i_row, i_col] = 3
                        cache['Number_of_Seasons'][:, i_row, i_col] = -1
                        pass
            except Exception as e:
                print(e)

    cache['Standing_Biomass'] = cache['Standing_Biomass'].round(2)
    cache['Min_min_PermanentIntegral'] = cache['Min_min_PermanentIntegral'].round(2)
    cache['Seasonal_Permanent_Integral'] = cache['Seasonal_Permanent_Integral'].round(2)
    cache['Season_Integral'] = cache['Season_Integral'].round(2)
    cache['Cyclic_Fraction'] = cache['Cyclic_Fraction'].round(2)
    cache['Active_Fraction_Integral'] = cache['Active_Fraction_Integral'].round(2)

    return cache


def _error_decoder(err):
    err_cod = {1: 'No data', 2: 'Scaling', 3: 'Off set', 4: '0-100 Scaling', 5: 'outlier filtering',
               6: 'gap filling', 7: 'Season and trend estimation', 8: 'To daily conversion', 9: 'madspan error',
               10: 'Trend conversion to daily', 11: 'Savinsky Golet', 12: 'Valley detection',
               13: 'Season detection', 14: 'Season mean', 15: 'Season metrics', 17: 'Statistical aggregation'}

    return err_cod[err]


def analyse(cube, client, param, template):
    """

    :param cube:
    :param client:
    :param param:
    :param action:
    :param out:
    :return:
    """
    try:
        mapped = xr.map_blocks(_process, cube, kwargs={'param': param}, template=template)

        return mapped

    except Exception as ex:

        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        print(message)

        logger.debug(f'Critical error in the main loop')
