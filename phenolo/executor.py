# -*- coding: utf-8 -*-

import logging

import numpy as np
import pandas as pd
from dask.distributed import as_completed

from phenolo import atoms

logger = logging.getLogger(__name__)

import line_profiler
import atexit
profile = line_profiler.LineProfiler()
atexit.register(profile.print_stats)


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


def _cache_concat(cache, indices, pxldrl):
    for name in indices:
        cache[name] = pd.concat([cache[name], getattr(pxldrl, name)])
    cache['sl'] = pd.concat([cache['sl'], getattr(pxldrl, 'sl')])
    cache['season'] = pd.concat([cache['season'], getattr(pxldrl, 'season')])
    cache['err'] = pd.concat([cache['err'], getattr(pxldrl, 'err')])


def _cache_def(indices, dim_val):
    """
    :param indices: indices that needed to be in cache {list}
    :param dim_val: list of years {pd.series}
    :param col_val: list of couliumns {int}
    :return: cache {pd.dataframe}
    """

    cache = {name: pd.DataFrame(index=dim_val) for name in indices}
    cache['season'] = pd.Series()
    cache['err'] = pd.Series()
    return cache


def _cache_reindex(cache, indices, index):
    for name in indices:
        cache[name].reindex(index, axis=1)
    cache['season'] = cache['season'].reindex(index).to_frame().traspose()
    cache['err'] = cache['error'].reindex(index).to_frame().traspose()


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
        key[col] = getattr(pxldrl, att)[:]
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


def _nodata_filler(out, abs_row):

    out.stb[abs_row, :, :] = np.nan
    out.mpi[abs_row, :, :] = np.nan
    out.sbd[abs_row, :, :] = np.nan
    out.sed[abs_row, :, :] = np.nan
    out.sl[abs_row, :, :] = np.nan
    out.spi[abs_row, :, :] = np.nan
    out.si[abs_row, :, :] = np.nan
    out.cf[abs_row, :, :] = np.nan
    out.afi[abs_row, :, :] = np.nan
    out.warn[abs_row, :, :] = np.nan
    out.n_seasons[abs_row, :] = np.nan
    out.err[abs_row, :] = np.nan


def analyse(cube, client, param, action, out):
    """

    :param cube:
    :param client:
    :param param:
    :param action:
    :param out:
    :return:
    """

    try:
        param.dim_val_yrs = pd.to_datetime(param.dim_val).year.unique()
        col_int = [*range(0, len(param.col_val))]

        s_param = client.scatter(param, broadcast=True)

        indices = ['stb', 'mpi', 'sbd', 'sed', 'sl', 'spi', 'si', 'cf', 'afi', 'warn']
        prg_bar = 0

        for chunk in np.array_split(range(0, len(param.row_val)), 3):

            chunked = cube.isel(dict([(param.row_nm, slice(chunk[0], chunk[-1]+1))])).compute()

            for rowi in range(0, chunked.sizes[param.row_nm]):
                row = chunked.isel(dict([(param.row_nm, rowi)]))
                abs_row = chunk[rowi]
                y_lst = _pxl_lst(row, param)

                if y_lst.any():
                    s_row = client.scatter(row, broadcast=True)
                    del row
                    cache = _cache_def(indices, param.dim_val_yrs)
                else:
                    _nodata_filler(out, abs_row)
                    print_progress_bar(rowi, len(param.row_val))
                    logger.debug(f'Row {rowi} processed')
                    continue

                futures = client.map(process, y_lst, **{'data': s_row, 'row': rowi, 'param': s_param, 'action': action})

                for future, result in as_completed(futures, with_results=True):
                    pxldrl = result
                    col = pxldrl.position[1]

                    if param.ovr_scratch:
                        try:
                            import phenolo.output as output
                            output.scratch_dump(pxldrl, param)
                        except Exception:
                            raise Exception

                    if pxldrl.error:
                        cache['err'].append(pd.Series([1], index=[col]))
                        cache['season'].append(pd.Series([np.NaN], index=[col]))
                        logger.debug(f'Error: {_error_decoder(pxldrl.errtyp)} in position:{pxldrl.position}')
                    else:
                        try:
                            _cache_concat(cache, indices, pxldrl)
                            if pxldrl.season_lng:
                                if pxldrl.season_lng <= 365.0:
                                    cache['season'].append(pd.Series([int(365 / pxldrl.season_lng)], index=[col]))
                                else:
                                    cache['season'].iloc[col] = int(pxldrl.season_lng)
                        except (RuntimeError, Exception, ValueError):
                            continue

                _cache_reindex(cache, indices, col_int)

                client.cancel(s_row)
                client.cancel(futures)

                out.stb[abs_row, :, :] = cache['stb'].transpose().values
                out.mpi[abs_row, :, :] = cache['mpi'].transpose().values

                out.sbd[abs_row, :, :] = cache['sbd'].transpose().values
                out.sed[abs_row, :, :] = cache['sed'].transpose().values
                out.sl[abs_row, :, :] = cache['sl'].transpose().values
                out.spi[abs_row, :, :] = cache['spi'].transpose().values
                out.si[abs_row, :, :] = cache['si'].transpose().values
                out.cf[abs_row, :, :] = cache['cf'].transpose().values
                out.afi[abs_row, :, :] = cache['afi'].transpose().values

                out.warn[abs_row, :, :] = cache['warn'].transpose().values

                out.n_seasons[abs_row, :] = cache['season'].values
                out.err[abs_row, :] = cache['err'].values

                prg_bar += 1
                print_progress_bar(prg_bar, len(param.row_val))

                logger.debug(f'Row {abs_row} has been processed')
        return out

    except Exception as ex:

        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        print(message)

        logger.debug(f'Critical error in the main loop, latest position, row {abs_row}, col {col}, error type {message}')
