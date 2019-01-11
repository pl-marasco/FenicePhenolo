import pandas as pd
from dask.distributed import Client, as_completed
import atoms
import logging

logger = logging.getLogger(__name__)


class Processor(object):
    def __init__(self, param):
        self.client = Client() #(processes=False, n_workers=10, threads_per_worker=2)  #
        self.param = param

    def _transducer(self, px, **kwargs):
        cube = kwargs.pop('data', '')
        action = kwargs.pop('action', '')

        row, col = px
        ts = cube.isel(dict([(self.param.col_nm, col)])).to_series().astype(float)
        pxdrl = atoms.PixelDrill(ts, px)
        return action(pxdrl, settings=self.param)

    def analyse(self, cube, out, action):

        try:
            for rowi in range(len(self.param.row_val)):
                row = cube.isel(dict([(self.param.row_nm, rowi)])).persist()

                dim_val = pd.to_datetime(self.param.dim_val).year.unique()
                col_val = range(0, len(self.param.col_val))
                t_sl = pd.DataFrame(index=dim_val, columns=col_val)
                t_spi = pd.DataFrame(index=dim_val, columns=col_val)
                t_si = pd.DataFrame(index=dim_val, columns=col_val)
                t_cf = pd.DataFrame(index=dim_val, columns=col_val)

                px_list = [item for item in self.param.pixel_list if item[0] == rowi]

                s_row = self.client.scatter(row, broadcast=True)
                futures = self.client.map(self._transducer, px_list, **{'data': s_row,
                                                                        'param': self.param,
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
                    # if column == 0:
                    #     continue

                    out.sl[rowi, column, :] = t_sl.iloc[:, column].values
                    out.spi[rowi, column, :] = t_spi.iloc[:, column].values
                    out.si[rowi, column, :] = t_si.iloc[:, column].values
                    out.cf[rowi, column, :] = t_cf.iloc[:, column].values

            return out

        except (RuntimeError, Exception, ValueError):
            logger.debug(f'Critical error in the main loop, latest position row {rowi}, col {col}')
