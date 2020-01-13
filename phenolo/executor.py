# -*- coding: utf-8 -*-

import logging
import time

import numpy as np
import pandas as pd
from dask.distributed import as_completed

from phenolo import atoms, analysis

logger = logging.getLogger(__name__)

import copy
import gc

# import line_profiler
# import atexit
# profile = line_profiler.LineProfiler()
# atexit.register(profile.print_stats)


class Processor(object):
    def __init__(self):
        pass


def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print(f'\r{prefix} |{bar}| {percent}{suffix}', end='')
    # Print New Line on Complete
    if iteration == total:
        print()


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
    param = kwargs.pop('param', '')
    row = kwargs.pop('row', '')

    pxldrl = atoms.PixelDrill(cube.isel(dict([(param.col_nm, px)])).to_series().astype(float), [row, px])

    return analysis.phenolo(pxldrl, settings=param)


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


def _cache_def(indices, dim_val, col_val):
    """
    :param indices: indices that needed to be in cache {list}
    :param dim_val: list of years {pd.series}
    :param col_val: list of couliumns {int}
    :return: cache {pd.dataframe}
    """

    cache = {name: pd.DataFrame(index=dim_val, columns=col_val) for name in indices}
    # cache['sl'] = pd.DataFrame(pd.Timedelta(0, unit='D'), index=dim_val, columns=col_val)
    cache['sl'] = pd.DataFrame(index=dim_val, columns=col_val)
    cache['season'] = pd.Series(0, index=col_val)
    cache['err'] = pd.Series(0, index=col_val)
    return cache


def _cache_cleaner(cache, indices, dim_val, col_val):

    for name in indices:
        cache[name] = cache[name][0:0].reindex(dim_val)
    cache['sl'] = pd.DataFrame(pd.Timedelta(0, unit='D'), index=dim_val, columns=col_val)
    cache['season'] = cache['season'][0:0].reindex(col_val)
    cache['err'] = cache['err'][0:0].reindex(col_val)
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
    try:
        key[col] = copy.deepcopy(getattr(pxldrl, att))
    except:
        pass
        # print(f'{att} | {pxldrl.position}')
    return


def _error_decoder(err):
    err_cod = {1: 'No data', 2: 'Scaling', 3: 'Off set', 4: '0-100 Scaling', 5: 'outlier filtering',
               6: 'gap filling', 7: 'Season and trend estimation', 8: 'To daily conversion', 9: 'madspan error',
               10: 'Trend conversion to daily', 11: 'Savinsky Golet', 12: 'Valley detection',
               13: 'Season detection', 14: 'Season mean', 15: 'Season metrics', 17: 'Statistical aggregation'}

    return err_cod[err]


def analyse(cube, client, param, out):
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

        indices = ['stb', 'mpi', 'sbd', 'sed', 'spi', 'si', 'cf', 'afi', 'warn']

        cache = _cache_def(indices, dim_val, col_val)
        prg_bar = 0

        for chunk in np.array_split(range(0, len(param.row_val)),2):

            chunked = cube.isel(dict([(param.row_nm, slice(chunk[0], chunk[-1]+1))])).compute()

            for rowi in range(0, chunked.sizes[param.row_nm]):
                row = chunked.isel(dict([(param.row_nm, rowi)]))
                y_lst = _pxl_lst(row, param)
                cache = _cache_def(indices, dim_val, col_val)

                if y_lst.any():
                    s_row = client.scatter(row, broadcast=True)
                    del row
                else:
                    print_progress_bar(rowi, len(param.row_val))
                    logger.debug(f'Row {rowi} processed')
                    continue

                futures = client.map(process, y_lst, **{'data': s_row, 'row': rowi, 'param': s_param})

                for future, pxldrl in as_completed(futures, with_results=True):
                    col = pxldrl.position[1]

                    if param.ovr_scratch:
                        try:
                            import phenolo.output as output
                            output.scratch_dump(pxldrl, param)
                        except Exception:
                            raise Exception

                    if pxldrl.error:
                        cache['err'].iloc[col] = 1
                        cache['season'].iloc[col] = np.NaN
                        logger.debug(f'Error: {_error_decoder(pxldrl.errtyp)} in position:{pxldrl.position}')
                    else:
                        try:
                            for key in cache:
                                if key is not 'season' and key is not 'err':
                                    _filler(cache[key], pxldrl, key, col)

                            if pxldrl.season_lng:
                                if pxldrl.season_lng <= 365.0:
                                    cache['season'].iloc[col] = int(365 / pxldrl.season_lng)
                                else:
                                    cache['season'].iloc[col] = int(pxldrl.season_lng)
                        except (RuntimeError, Exception, ValueError):
                            continue
                    del pxldrl, future

                abs_row = chunk[rowi]
                out.stb[:, abs_row, :] = cache['stb'].values
                out.mpi[:, abs_row, :] = cache['mpi'].values

                out.sbd[:, abs_row, :] = cache['sbd'].values
                out.sed[:, abs_row, :] = cache['sed'].values
                out.sl[:, abs_row, :] = cache['sl'].values
                out.spi[:, abs_row, :] = cache['spi'].values
                out.si[:, abs_row, :] = cache['si'].values
                out.cf[:, abs_row, :] = cache['cf'].values
                out.afi[:, abs_row, :] = cache['afi'].values

                out.warn[:, abs_row, :] = cache['warn'].values

                out.n_seasons[abs_row] = cache['season'].values
                out.err[abs_row] = cache['err'].values

                del cache
                client.cancel(s_row)
                client.cancel(futures)

                prg_bar += 1
                print_progress_bar(prg_bar, len(param.row_val))

                logger.debug(f'Row {abs_row} has been processed')

        return out

    except Exception as ex:

        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        print(message)

        logger.debug(f'Critical error in the main loop, latest position, row {abs_row}, col {col}, error type {message}')
