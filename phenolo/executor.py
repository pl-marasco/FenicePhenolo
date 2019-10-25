# -*- coding: utf-8 -*-

import logging

import numpy as np
import pandas as pd
from dask.distributed import as_completed
import copy

from phenolo import atoms

logger = logging.getLogger(__name__)


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
    for key, name in cache.items():
        if name not in ['season', 'err']:
            cache[name] = pd.concat([cache[name], getattr(pxldrl, name)], axis=1, copy=False)
        #cache[name] = pd.concat([cache[name], getattr(pxldrl, name)], axis=1, copy=False)


def _cache_def(indices, dim_val):
    """
    :param indices: indices that needed to be in cache {list}
    :param dim_val: list of years {pd.series}
    :return: cache {pd.dataframe}
    """

    cache = {name: pd.DataFrame(index=dim_val) for name in indices}
    cache['season'] = pd.Series()
    cache['err'] = pd.Series()
    return cache


def _cache_reindex(cache, indices, index):
    for name in indices:
        cache[name] = cache[name].reindex(index, axis=1, copy=False)
    cache['season'] = cache['season'].reindex(index, copy=False).to_frame().transpose()
    if not cache['err'].empty:
        cache['err'] = cache['error'].reindex(index, copy=False).to_frame().transpose()


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
    out.sl[abs_row, :, :] = np.nan  # todo verify if is correct in float instead of int
    out.spi[abs_row, :, :] = np.nan
    out.si[abs_row, :, :] = np.nan
    out.cf[abs_row, :, :] = np.nan
    out.afi[abs_row, :, :] = np.nan
    out.warn[abs_row, :, :] = np.nan
    out.n_seasons[abs_row, :] = 0
    out.err[abs_row, :] = 0


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
        param.yrs_index = pd.to_datetime(param.dim_val).year.unique()
        col_int = [*range(0, len(param.col_val))]

        s_param = client.scatter(param, broadcast=True)

        indices = ['stb', 'mpi', 'sbd', 'sed', 'sl', 'spi', 'si', 'cf', 'afi', 'warn']
        prg_bar = 0

        for chunk in np.array_split(range(0, len(param.row_val)), 3):

            chunked = cube.isel(dict([(param.row_nm, slice(chunk[0], chunk[-1] + 1))])).compute()

            for rowi in range(0, chunked.sizes[param.row_nm]):
                row = chunked.isel(dict([(param.row_nm, rowi)]))
                abs_row = chunk[rowi]
                y_lst = _pxl_lst(row, param)

                if y_lst.any():
                    s_row = client.scatter(row, broadcast=True)
                    del row
                    cache = _cache_def(indices, param.yrs_index)
                else:
                    _nodata_filler(out, abs_row)
                    print_progress_bar(rowi, len(param.row_val))
                    logger.debug(f'Row {rowi} processed')
                    continue

                futures = client.map(process, y_lst, **{'data': s_row, 'row': rowi, 'param': s_param, 'action': action})

                for future, pxldrl in as_completed(futures, with_results=True):

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
                                    # TODO can be faster to append values to a list and then concatenate
                                else:
                                    cache['season'].append(pd.Series([int(pxldrl.season_lng)], index=[col]))
                        except (RuntimeError, Exception, ValueError):
                            continue

                _cache_reindex(cache, indices, col_int)

                out.stb[abs_row, :, :] = cache['stb'].transpose().to_numpy(copy=True)
                out.mpi[abs_row, :, :] = cache['mpi'].transpose().to_numpy(copy=True)

                out.sbd[abs_row, :, :] = cache['sbd'].transpose().to_numpy(copy=True)
                out.sed[abs_row, :, :] = cache['sed'].transpose().to_numpy(copy=True)
                out.sl[abs_row, :, :] = (cache['sl'] / np.timedelta64(1, 'D')).transpose().to_numpy(copy=True)
                out.spi[abs_row, :, :] = cache['spi'].transpose().to_numpy(copy=True)
                out.si[abs_row, :, :] = cache['si'].transpose().to_numpy(copy=True)
                out.cf[abs_row, :, :] = cache['cf'].transpose().to_numpy(copy=True)
                out.afi[abs_row, :, :] = cache['afi'].transpose().to_numpy(copy=True)

                out.warn[abs_row, :, :] = cache['warn'].transpose().to_numpy(copy=True)
                out.n_seasons[abs_row, :] = cache['season'].to_numpy(copy=True)

                if not cache['err'].empty:
                    out.err[abs_row, :] = cache['err'].to_numpy(copy=True)
                else:
                    out.err[abs_row, :] = 0

                # client.cancel(s_row)
                # client.cancel(futures)

                prg_bar += 1
                print_progress_bar(prg_bar, len(param.row_val))
                logger.debug(f'Row {abs_row} has been processed')

                del cache

        return out

    except Exception as ex:

        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        print(message)

        logger.debug(
            f'Critical error in the main loop, latest position, row {abs_row}, col {col}, error type {message}')
