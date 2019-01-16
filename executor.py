import pandas as pd
from dask.distributed import Client, as_completed
import atoms
import logging
import numpy as np

logger = logging.getLogger(__name__)


class Processor(object):
    def __init__(self):
        pass


def transducer(px, **kwarg):

    cube = kwarg.pop('data', '')
    action = kwarg.pop('action', '')
    param = kwarg.pop('param', '')

    row, col = px
    ts = cube.isel(dict([(param.col_nm, col)])).to_series().astype(float)
    pxdrl = atoms.PixelDrill(ts, px)
    return action(pxdrl, settings=param)


def analyse(cube, param, action, out):

    client = Client(processes=False, n_workers=1, threads_per_worker=1)  #

    try:
        for rowi in range(len(param.row_val)):
            row = cube.isel(dict([(param.row_nm, rowi)])).persist()

            dim_val = pd.to_datetime(param.dim_val).year.unique()
            col_val = range(0, len(param.col_val))
            t_sl = pd.DataFrame(index=dim_val, columns=col_val)
            t_spi = pd.DataFrame(index=dim_val, columns=col_val)
            t_si = pd.DataFrame(index=dim_val, columns=col_val)
            t_cf = pd.DataFrame(index=dim_val, columns=col_val)

            px_list = [item for item in param.pixel_list if item[0] == rowi]
            # px_list = np.where(param.pixel_list[:, 0] == rowi)

            # px_list = pd.DataFrame(param.pixel_list).loc[param.pixel_list[:, 0] == rowi].values.tolist()

            s_row = client.scatter(row, broadcast=True)
            s_param = client.scatter(param, broadcast=True)

            futures = client.map(transducer, px_list, **{'data': s_row,
                                                         'param': s_param,
                                                         'action': action})

            for future, result in as_completed(futures, with_results=True):
                pxdrl = result
                row, col = pxdrl.position
                t_sl.iloc[:, col] = pxdrl.sl
                t_spi.iloc[:, col] = pxdrl.spi
                t_si.iloc[:, col] = pxdrl.si
                t_cf.iloc[:, col] = pxdrl.cf
                future.cancel()

            #  TODO create a faster approach to instantiate the netcdf result file
            for column in t_sl:

                out.sl[rowi, column, :] = t_sl.iloc[:, column].values
                out.spi[rowi, column, :] = t_spi.iloc[:, column].values
                out.si[rowi, column, :] = t_si.iloc[:, column].values
                out.cf[rowi, column, :] = t_cf.iloc[:, column].values

        return out

    except (RuntimeError, Exception, ValueError):
        logger.debug(f'Critical error in the main loop, latest position row {rowi}, col {col}')
