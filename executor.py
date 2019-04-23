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

            row = cube.isel(dict([(param.row_nm, rowi)])).compute()

            reduced = row.reduce(np.percentile, dim=param.dim_nm, q=param.qt)
            med = reduced.where(((reduced > param.min_th) & (reduced < param.max_th)))
            finite = med.reduce(np.isfinite)
            y_lst = np.argwhere(finite.values).flatten()

            if y_lst.any():
                s_row = client.scatter(row, broadcast=True)
            else:
                del row; continue

            dim_val = pd.to_datetime(param.dim_val).year.unique()
            col_val = range(0, len(param.col_val))
            t_sb = pd.DataFrame(index=dim_val, columns=col_val)
            t_se = pd.DataFrame(index=dim_val, columns=col_val)
            t_sl = pd.DataFrame(pd.Timedelta(0, unit='D'), index=dim_val, columns=col_val)
            t_spi = pd.DataFrame(index=dim_val, columns=col_val)
            t_si = pd.DataFrame(index=dim_val, columns=col_val)
            t_cf = pd.DataFrame(index=dim_val, columns=col_val)
            t_afi = pd.DataFrame(index=dim_val, columns=col_val)

            t_season = pd.Series(0, index=col_val)
            t_err = pd.Series(0, index=col_val)

            futures = client.map(process, y_lst, **{'data': s_row, 'row': rowi, 'param': s_param, 'action': action})

            for future, pxldrl in as_completed(futures, with_results=True):
                col = pxldrl.position[1]
                if pxldrl.error:
                    t_err.iloc[col] = 1
                    t_season.iloc[col] = 0
                    logger.debug(f'Error: {pxldrl.errtyp} in position:{pxldrl.position}')
                    print(f'Error: {pxldrl.errtyp} in position:{pxldrl.position}')
                else:
                    try:
                        t_sb.iloc[:, col] = pxldrl.sb[:]
                        t_se.iloc[:, col] = pxldrl.se[:]
                        t_sl.iloc[:, col] = pxldrl.sl[:]
                        t_spi.iloc[:, col] = pxldrl.spi[:]
                        t_si.iloc[:, col] = pxldrl.si[:]
                        t_cf.iloc[:, col] = pxldrl.cf[:]
                        t_afi.iloc[:, col] = pxldrl.afi[:]
                        if pxldrl.season_lng:
                            if pxldrl.season_lng <= 365.0:
                                t_season.iloc[col] = int(365/pxldrl.season_lng)
                            else:
                                t_season.iloc[col] = int(pxldrl.season_lng)
                    except (RuntimeError, Exception, ValueError):
                        continue

                # client.cancel(future)
                # del future, pxldrl

            out.sb[:, rowi, :] = t_sb
            out.se[:, rowi, :] = t_se
            out.sl[:, rowi, :] = t_sl
            out.spi[:, rowi, :] = t_spi
            out.si[:, rowi, :] = t_si
            out.cf[:, rowi, :] = t_cf
            out.n_seasons[rowi] = t_season.values
            out.err[rowi] = t_err.values

            try:
                out.root.sync()
            except (RuntimeError, Exception, ValueError):
                logger.debug(f'Error in the sync')

            # logger.debug(f'Row {rowi} processed')

            # client.cancel(s_row)

            # client.cancel(futures)
            # del futures, t_sl, t_cf, t_si, t_spi, row, s_row
            # gc.collect()
        return out

    except Exception as ex:

        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        print(message)

        logger.debug(f'Critical error in the main loop, latest position row {rowi}, col {col}, error type {message}')
