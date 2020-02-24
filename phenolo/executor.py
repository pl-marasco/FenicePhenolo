# -*- coding: utf-8 -*-

import logging
import time

import numpy as np
import pandas as pd
from dask.distributed import as_completed

from phenolo import atoms, analysis

from multiprocessing import Process, shared_memory, JoinableQueue, Event
import queue

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

    pxldrl = atoms.PixelDrill(cube.isel(dict([(param.col_nm, px), (param.row_nm, row)])).to_series().astype(float),
                              [row, px])

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


def _cache_def(indices, dim_sz, col_sz):
    """
    :param indices: indices that needed to be in cache {list}
    :param dim_sz: list of years {pd.series}
    :param col_sz: list of couliumns {int}
    :return: cache {pd.dataframe}
    """
    cache = {}
    shared = []

    mock = np.empty((dim_sz, col_sz), dtype=np.float64)
    mock_int = np.empty((dim_sz, col_sz), dtype=np.int64)

    for name in indices:
        shm = shared_memory.SharedMemory(create=True, size=mock.nbytes, name=name)
        shared.append(shm)
        d = {name: np.ndarray((dim_sz, col_sz), dtype=np.float64, buffer=shm.buf)}
        d[name][:] = np.NaN
        cache.update(d)

    # season length
    shm_sl = shared_memory.SharedMemory(create=True, size=mock.nbytes, name='sl')
    shared.append(shm_sl)
    sl = np.ndarray((dim_sz, col_sz), dtype=np.float64, buffer=shm_sl.buf)  # non mi torna il formato
    sl[:] = np.NaN
    cache.update({'sl': sl})

    # season
    shm_season = shared_memory.SharedMemory(create=True, size=mock_int.nbytes, name='season')
    shared.append(shm_season)
    season = np.ndarray(col_sz, dtype=np.int64, buffer=shm_season.buf)
    season[:] = 0
    cache.update({'season': season})

    # update
    shm_err = shared_memory.SharedMemory(create=True, size=mock_int.nbytes, name='err')
    shared.append(shm_err)
    update = np.ndarray(col_sz, dtype=np.int64, buffer=shm_err.buf)
    update[:] = 0
    cache.update({'err': update})

    return cache, shared


def _cache_shared(indices, dim_sz, col_sz):
    """
    :param indices: indices that needed to be in cache {list}
    :param dim_sz: list of years {pd.series}
    :param col_sz: list of couliumns {int}
    :return: cache {pd.dataframe}
    """
    cache = {}
    shared = []

    for name in indices:
        shm = shared_memory.SharedMemory(name=name)
        shared.append(shm)
        d = {name: np.ndarray((dim_sz, col_sz), dtype=np.float64, buffer=shm.buf)}
        cache.update(d)

    # season length
    shm_sl = shared_memory.SharedMemory(name='sl')
    shared.append(shm_sl)
    sl = np.ndarray((dim_sz, col_sz), dtype=np.float64, buffer=shm_sl.buf)  # non mi torna il formato
    cache.update({'sl': sl})

    # season
    shm_season = shared_memory.SharedMemory(name='season')
    shared.append(shm_season)
    season = np.ndarray(col_sz, dtype=np.int64, buffer=shm_season.buf)
    cache.update({'season': season})

    # update
    shm_err = shared_memory.SharedMemory(name='err')
    shared.append(shm_err)
    update = np.ndarray(col_sz, dtype=np.int64, buffer=shm_err.buf)
    cache.update({'err': update})

    return cache, shared


def _cache_shared_unlink(shared):
    for shm in shared:
        shm.unlink()


def _cache_shared_close(shared):
    for shm in shared:
        shm.close()


def _cache_cleaner(cache, indices, dim_val, col_val):

    for name in indices:
        cache[name] = cache[name][0:0].reindex(dim_val)
    cache['sl'] = pd.DataFrame(pd.Timedelta(0, unit='D'), index=dim_val, columns=col_val)
    cache['season'] = cache['season'][0:0].reindex(col_val)
    cache['err'] = cache['err'][0:0].reindex(col_val)
    return cache


def _pxl_filler(r_queue, param, signal):

    cache, shared = _cache_shared(param.indices, param.dim_sz, param.col_sz)

    while True:
        try:
            pxldrl = r_queue.get(False)
            try:
                col = pxldrl.position[1]

                if param.ovr_scratch:
                    try:
                        import phenolo.output as output
                        output.scratch_dump(pxldrl, param)
                    except Exception as e:
                        raise Exception

                if pxldrl.error:
                    cache['err'][col] = 1
                    cache['season'][col] = -1
                    logger.debug(f'Error: {_error_decoder(pxldrl.errtyp)} in position:{pxldrl.position}')
                else:
                    try:
                        for key in cache:
                            if key is not 'season' and key is not 'err':
                                _filler(cache[key], pxldrl, key, col)

                        if pxldrl.season_lng:
                            if pxldrl.season_lng <= 365.0:
                                cache['season'][col] = int(365 / pxldrl.season_lng)
                            else:
                                cache['season'][col] = int(pxldrl.season_lng)
                                # TODO add more seasons sub division

                    except Exception as e:
                        logger.error('Error in a worker during the filling')
                        pass
            except Exception as e:
                print(e)
            finally:
                r_queue.task_done()
        except queue.Empty:
            if signal.is_set():
                try:
                    _cache_shared_close(shared)
                except Exception as e:
                    logger.error('Error during shared memory closure')
                finally:
                    break
            else:
                continue


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
        key[:, col] = getattr(pxldrl, att) #to_numpy(copy=False)
    except Exception as e:
        pass
        # TO DO proper manage the exception
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
    try:
        param.dim_unq_val = pd.to_datetime(param.dim_val).year.unique()
        param.col_sz = len(param.col_val)
        rl_col_val = range(0, param.col_sz)
        param.dim_sz = len(param.dim_unq_val)

        param.indices = ['stb', 'mpi', 'sbd', 'sed', 'spi', 'si', 'cf', 'afi', 'warn']

        s_param = client.scatter(param, broadcast=True)

        # cache = _cache_def(indices, len(param.dim_unq_val), len(col_val))
        prg_bar = 0
        n_chunks = 2

        for chunk in np.array_split(range(0, len(param.row_val)), n_chunks):

            chunked = cube.isel(dict([(param.row_nm, slice(chunk[0], chunk[-1]+1))])).compute()

            # # -->
            # chunked.chunk({param.row_nm: 100, param.col_nm: 100, param.dim_nm: -1})
            # # <--

            chnk_scat = client.scatter(chunked, broadcast=True)

            # -->
            quantile = chunked.quantile(0.2, param.dim_nm)
            where = quantile.where(((quantile > param.min_th) & (quantile < param.max_th)))
            isfinite = where.reduce(np.isfinite)
            # <--

            pxl_lst = np.argwhere(isfinite.values)

            for rowi in range(0, chunked.sizes[param.row_nm]):

                y_lst = pxl_lst[pxl_lst[:, 0] == rowi, :][:, 1]

                if not y_lst.any():
                    print_progress_bar(rowi, len(param.row_val))
                    logger.debug(f'Row {rowi} processed')
                    continue

                cache, shared = _cache_def(param.indices, param.dim_sz, param.col_sz)
                stop_event = Event()
                r_queue = JoinableQueue()

                processes = [Process(name=str(n), target=_pxl_filler, args=(r_queue, param, stop_event)) for n in range(5)]
                for worker in processes:
                    worker.start()

                futures = client.map(process, y_lst, **{'data': chnk_scat, 'row': rowi, 'param': s_param})

                for future, pxldrl in as_completed(futures, with_results=True):
                    r_queue.put(pxldrl)

                r_queue.join()
                stop_event.set()

                for worker in processes:
                    worker.join()

                abs_row = chunk[rowi]
                out.stb[:, abs_row, :] = cache['stb']
                out.mpi[:, abs_row, :] = cache['mpi']

                out.sbd[:, abs_row, :] = cache['sbd']
                out.sed[:, abs_row, :] = cache['sed']
                out.sl[:, abs_row, :] = cache['sl']
                out.spi[:, abs_row, :] = cache['spi']
                out.si[:, abs_row, :] = cache['si']
                out.cf[:, abs_row, :] = cache['cf']
                out.afi[:, abs_row, :] = cache['afi']

                out.warn[:, abs_row, :] = cache['warn']

                out.n_seasons[abs_row] = cache['season']
                out.err[abs_row] = cache['err']

                _cache_shared_close(shared)
                _cache_shared_unlink(shared)

                # client.cancel(s_row)
                # client.cancel(futures)

                prg_bar += 1
                print_progress_bar(prg_bar, len(param.row_val))

                logger.debug(f'Row {abs_row} has been processed')

        return out

    except Exception as ex:

        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        print(message)

        logger.debug(f'Critical error in the main loop, latest position, row {abs_row}, error type {message}')
