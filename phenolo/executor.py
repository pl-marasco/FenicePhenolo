# -*- coding: utf-8 -*-

import logging

import numpy as np
import pandas as pd
from dask.distributed import as_completed

from phenolo import atoms
import time

logger = logging.getLogger(__name__)


class Processor(object):
    def __init__(self):
        pass


def process(px, **kwargs):
    """
    Wrapper for a function pass as "action" over a Pandas time series.

    :param px: pandas time series position {int}
    :param kwargs: **{'data': xarray cube,
                      'action': function to be apply,
                      'param': param object
                      'row': row position in the cube as {int}
    :return: Obj{pxdrl}
    """
    cube = kwargs.pop('data', '')
    action = kwargs.pop('action', '')
    param = kwargs.pop('param', '')
    row = kwargs.pop('row', '')

    pxldrl = atoms.PixelDrill(cube.isel(dict([(param.col_nm, px)])).to_series().astype(float), [row, px])

    return action(pxldrl, settings=param)


def _pre_feeder(nxt_row, param):
    return _pxl_lst(nxt_row, param)


def _pxl_lst(row, param):
    """
    Map pixels that must analyzed
    :param row: pd.Dataframe 2D
    :param param: param Obj
    :return: array of int representing pxl row position
    """
    reduced = row.reduce(np.percentile, dim=param.dim_nm, q=param.qt)
    med = reduced.where(((reduced > param.min_th) & (reduced < param.max_th)))
    finite = med.reduce(np.isfinite)
    y_lst = np.argwhere(finite.values).flatten()
    return y_lst


def _cache_def(dim_val, col_val):
    """

    :param dim_val:
    :param col_val:
    :return:
    """
    attr = ['sb', 'se', 'sl', 'spi', 'si', 'cf', 'afi']
    cache = {name: pd.DataFrame(index=dim_val, columns=col_val) for name in attr}
    cache['sl'] = pd.DataFrame(pd.Timedelta(0, unit='D'), index=dim_val, columns=col_val)
    cache['season'] = pd.Series(0, index=col_val)
    cache['err'] = pd.Series(0, index=col_val)
    return cache


def _filler(key, pxldrl, att, col):
    """
    Fill the dictionary with the passes key and values
    :param key: specific key to be filled
    :param pxldrl:
    :param att:
    :param col:
    :return:
    """
    key.iloc[:, col] = getattr(pxldrl, att)[:]
    return


def analyse(cube, client, param, action, out):
    """

    :param cube:
    :param client:
    :param param:
    :param action:
    :param out:
    :return:
    """
    s_param = client.scatter(param, broadcast=True)

    try:
        nxt_row, nxt_y_lst, nxt_cache = [None] * 3

        dim_val = pd.to_datetime(param.dim_val).year.unique()
        col_val = range(0, len(param.col_val))

        nxt = True

        for rowi in range(len(param.row_val)):
            if rowi == 0:
                row = cube.isel(dict([(param.row_nm, rowi)])).compute()
                y_lst = _pxl_lst(row, param)
                cache = _cache_def(dim_val, col_val)
            elif rowi == len(param.row_val)-1:
                row = cube.isel(dict([(param.row_nm, rowi)])).compute()
                y_lst = _pxl_lst(row, param)
                cache = _cache_def(dim_val, col_val)
                nxt = False
            else:
                row = nxt_row
                cache = nxt_cache
                y_lst = nxt_y_lst

            if y_lst.any():
                s_row = client.scatter(row, broadcast=True)
            else:
                continue

            futures = client.map(process, y_lst, **{'data': s_row, 'row': rowi, 'param': s_param, 'action': action})
            seq = as_completed(futures, with_results=True)

            if nxt:
                nxt_row = cube.isel(dict([(param.row_nm, rowi + 1)])).compute()
                preload = client.submit(_pre_feeder, **{'nxt_row': nxt_row,
                                                        'param': s_param})
                seq.add(preload)
                cleaner = client.submit(_cache_def, **{'dim_val': dim_val,
                                                       'col_val': col_val},
                                        priority=10)
                seq.add(cleaner)

            for future in seq:
                result = future[1]
                if isinstance(result, atoms.PixelDrill):
                    pxldrl = result
                    col = pxldrl.position[1]
                    if pxldrl.error:
                        cache['err'].iloc[col] = 1
                        cache['season'].iloc[col] = 0
                        logger.debug(f'Error: {pxldrl.errtyp} in position:{pxldrl.position}')
                        print(f'Error: {pxldrl.errtyp} in position:{pxldrl.position}')
                    else:
                        try:
                            for key in cache:
                                if key is not 'season' and key is not 'err':
                                    _filler(cache[key], pxldrl, key, col)

                            if pxldrl.season_len:
                                if pxldrl.season_len <= 365.0:
                                    cache['season'].iloc[col] = int(365 / pxldrl.season_len)
                                else:
                                    cache['season'].iloc[col] = int(pxldrl.season_len)

                        except (RuntimeError, Exception, ValueError):
                            continue
                elif isinstance(result, dict):
                    nxt_cache = result
                else:
                    nxt_y_lst = result

            out.sb[:, rowi, :] = cache['sb'].values
            out.se[:, rowi, :] = cache['se'].values
            out.sl[:, rowi, :] = cache['sl'].values
            out.spi[:, rowi, :] = cache['spi'].values
            out.si[:, rowi, :] = cache['si'].values
            out.cf[:, rowi, :] = cache['cf'].values
            out.n_seasons[rowi] = cache['season'].values
            out.err[rowi] = cache['err'].values

            # try:
            #     out.root.sync()
            # except (RuntimeError, Exception, ValueError):
            #     logger.debug(f'Error in the sync')

            logger.debug(f'Row {rowi} processed')

        return out

    except Exception as ex:

        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        print(message)

        logger.debug(f'Critical error in the main loop, latest position row {rowi}, col {col}, error type {message}')
