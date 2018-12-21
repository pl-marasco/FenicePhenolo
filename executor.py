import numpy as np
from dask.distributed import Client, as_completed
import atoms
import logging

logger = logging.getLogger(__name__)


class Processor(object):
    def __init__(self, param):
        self.client = Client(processes=False)
        self.param = param

    def _transducer(self, px, **kwargs):
        cube = kwargs.pop('data', '')
        action = kwargs.pop('action', '')

        row, col = px
        ts = cube.isel(dict([(self.param.col_nm, col)])).to_series().astype(float)
        pxdrl = atoms.PixelDrill(ts, px)
        return action(pxdrl, settings=self.param)

    def analyse(self, cube, scratch, action):

        try:
            for rowi in range(len(self.param.row_val)):
                row = cube.isel(dict([(self.param.row_nm, rowi)])).persist()
                vrow = np.empty((1, len(self.param.col_val), len(self.param.dim_val)))
                px_list = [item for item in self.param.pixel_list if item[0] == rowi]

                s_row = self.client.scatter(row)
                futures = self.client.map(self._transducer, px_list, **{'data': s_row,
                                                                        'param': self.param,
                                                                        'action': action})

                for future, result in as_completed(futures, with_results=True):
                    # print(result)
                    ts, pixels = result
                    row, col = pixels

            return scratch

        except (RuntimeError, Exception, ValueError):
            logger.debug('Error during in the process of climat substitution')