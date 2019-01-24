import pandas as pd
from dask.distributed import Client, as_completed
import atoms
import logging
import numpy as np

logger = logging.getLogger(__name__)


class Processor(object):
    def __init__(self):
        pass


def transducer(px, **kwargs):

    cube = kwargs.pop('data', '')
    action = kwargs.pop('action', '')
    param = kwargs.pop('param', '')

    row, col = px

    ts = cube.isel(dict([(param.col_nm, col)])).to_series().astype(float)
    pxldrl = atoms.PixelDrill(ts, px)
    return action(pxldrl, settings=param)


def analyse(cube, param, action, out):

    localproc = param.processes
    n_workers = param.n_workers
    threads_per_worker = param.threads_per_worker

    if ~localproc and n_workers and threads_per_worker:
        client = Client(processes=localproc, n_workers=n_workers, threads_per_worker=threads_per_worker)
    elif localproc and n_workers:
        client = Client(n_workers=n_workers, threads_per_worker=threads_per_worker)
    else:
        client = Client()

    try:
        for rowi in range(len(param.row_val)):
            row = cube.isel(dict([(param.row_nm, rowi)])).persist()

            dim_val = pd.to_datetime(param.dim_val).year.unique() # <-- pd.to_datetime(pd.to_datetime(param.dim_val).year.unique(), format='%Y')
            col_val = range(0, len(param.col_val))
            t_sl = pd.DataFrame(index=dim_val, columns=col_val)
            t_spi = pd.DataFrame(index=dim_val, columns=col_val)
            t_si = pd.DataFrame(index=dim_val, columns=col_val)
            t_cf = pd.DataFrame(index=dim_val, columns=col_val)

            px_list = [item for item in param.pixel_list if item[0] == rowi]

            s_row = client.scatter(row, broadcast=True)
            s_param = client.scatter(param, broadcast=True)

            futures = client.map(transducer, px_list, **{'data': s_row,
                                                         'param': s_param,
                                                         'action': action})

            for future, result in as_completed(futures, with_results=True):
                    pxldrl = result
                    row, col = pxldrl.position
                    t_sl.iloc[:, col] = pxldrl.sl
                    t_spi.iloc[:, col] = pxldrl.spi
                    t_si.iloc[:, col] = pxldrl.si
                    t_cf.iloc[:, col] = pxldrl.cf

                    if pxldrl.error:
                        print(pxldrl.position)

            # client.restart()

            out.sl[rowi] = np.expand_dims(t_sl.transpose().values, axis=0)
            out.spi[rowi] = np.expand_dims(t_spi.transpose().values, axis=0)
            out.si[rowi] = np.expand_dims(t_si.transpose().values, axis=0)
            out.cf[rowi] = np.expand_dims(t_cf.transpose().values, axis=0)

            client.cancel(s_row)
            client.cancel(s_param)
            client.cancel(futures)

            del t_sl, t_cf, t_si, t_spi

        return out

    except (RuntimeError, Exception, ValueError):
        logger.debug(f'Critical error in the main loop, latest position row {rowi}, col {col}')
