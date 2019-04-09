# -*- coding: utf-8 -*-
import sys

import gc
import pandas as pd
from dask.distributed import Client, as_completed
import atoms
import logging
import numpy as np

logger = logging.getLogger(__name__)


class Processor(object):
    def __init__(self):
        pass


def process(px, **kwargs):

    cube = kwargs.pop('data', '')
    action = kwargs.pop('action', '')
    param = kwargs.pop('param', '')
    row = kwargs.pop('row', '')

    pxldrl = atoms.PixelDrill(cube.isel(dict([(param.col_nm, px)])).to_series().astype(float), [row, px])

    return action(pxldrl, settings=param)


def analyse(cube, client, param, action, out):

    s_param = client.scatter(param, broadcast=True)

    try:
        for rowi in range(len(param.row_val)):

            trashold = 250
            q_trashold = 0.2  # TODO change the percentage to be a parameter
            row = cube.isel(dict([(param.row_nm, rowi)])).compute()
            reduced = row.reduce(np.percentile, dim=param.dim_nm, q=q_trashold)
            med = reduced.where(reduced < trashold)
            finite = med.reduce(np.isfinite)
            y_lst = np.argwhere(finite.values).flatten()

            # med = ~np.less(row.quantile(q_trashold, dim=param.dim_nm), trashold)
            # masked = np.ma.MaskedArray(row, np.resize(med, (row.sizes[param.dim_nm], len(med))))
            # px_list = list(map(lambda x: [rowi, x], np.argwhere(np.isfinite(med.values)).flatten()))

            if y_lst.any():
                s_row = client.scatter(row, broadcast=True)
            else:
                del row; continue

            dim_val = pd.to_datetime(param.dim_val).year.unique() # <-- pd.to_datetime(pd.to_datetime(param.dim_val).year.unique(), format='%Y')
            col_val = range(0, len(param.col_val))
            t_sl = pd.DataFrame(index=dim_val, columns=col_val)
            t_spi = pd.DataFrame(index=dim_val, columns=col_val)
            t_si = pd.DataFrame(index=dim_val, columns=col_val)
            t_cf = pd.DataFrame(index=dim_val, columns=col_val)

            t_season = pd.Series(0, index=col_val)
            t_err = pd.Series(0, index=col_val)

            futures = client.map(process, y_lst, **{'data': s_row, 'row': rowi, 'param': s_param, 'action': action})

            for future, pxldrl in as_completed(futures, with_results=True):
                col = pxldrl.position[1]
                if pxldrl.error:
                    # t_err.iloc[col] = 1
                    logger.debug(f'Error: {pxldrl.errtyp} in position:{pxldrl.position}')
                    print(f'Error: {pxldrl.errtyp} in position:{pxldrl.position}')
                else:
                    t_sl.iloc[:, col] = pxldrl.sl[:]
                    t_spi.iloc[:, col] = pxldrl.spi[:]
                    t_si.iloc[:, col] = pxldrl.si[:]
                    t_cf.iloc[:, col] = pxldrl.cf[:]
                    if pxldrl.season_lng:
                        if pxldrl.season_lng < 365:
                            t_season[col] = int(365/pxldrl.season_lng)
                        else:
                            t_season[col] = int(pxldrl.season_lng)

                # client.cancel(future)
                # del future, pxldrl

            out.sl[rowi] = np.expand_dims(t_sl.transpose().values, axis=0)
            out.spi[rowi] = np.expand_dims(t_spi.transpose().values, axis=0)
            out.si[rowi] = np.expand_dims(t_si.transpose().values, axis=0)
            out.cf[rowi] = np.expand_dims(t_cf.transpose().values, axis=0)
            out.n_season = np.expand_dims(t_season.transpose().values, axis=0)
            out.err[rowi] = np.expand_dims(t_err.transpose().values, axis=0)

            try:
                out.root.sync()
            except (RuntimeError, Exception, ValueError):
                logger.debug(f'Error in the sync')

            # client.cancel(s_row)
            # client.cancel(futures)
            #
            # del futures, t_sl, t_cf, t_si, t_spi, row, s_row
            # gc.collect()
        return out

    except Exception as ex:

        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        print(message)

        logger.debug(f'Critical error in the main loop, latest position row {rowi}, col {col}, error type {message}')
