# -*- coding: utf-8 -*-
# !/usr/bin/env python

import nodata
import sys, os, logging
from netCDF4 import Dataset
import numpy as np
from dask.distributed import Client, as_completed

logger = logging.getLogger(__name__)


class PreProcessor(object):
    def __init__(self):
        pass


def spot_no_data(param, cube, scratch):

    client = Client()

    for rowi in range(len(param.row_val)):

        row = cube.isel(dict([(param.row_nm, rowi)])).persist()
        vrow = np.empty((1, len(param.col_val), len(param.dim_val)))
        px_list = [item for item in param.pixel_list if item[0] == rowi]

        s_row = client.scatter(row)
        futures = client.map(nodata.ext_fix, px_list, **{'data': s_row, 'param': param})

        for future, result in as_completed(futures, with_results=True):
            # print(result)
            ts, pixels = result
            row, col = pixels
            vrow[:, col, :] = ts.values

        scratch.container[rowi, :, :] = vrow

    return scratch

        # try:
        #     fixed_cube = nodata.fix(self, self.cube, type='SCE')
        #     return fixed_cube
        # except (RuntimeError, Exception, ValueError):
        #     logger.debug('Error during the process of fixing')
        #     sys.exit(1)


def rescale(self, cube):
    # Rescale values to  0-100
    try:
        return ((cube - self.min) / (self.max - self.min)) * 100
    except (RuntimeError, Exception, ValueError):
        print('Error in rescaling, in position')
        logger.debug('Error in rescaling')
        sys.exit()



#region
    # def _filter_outlayer(self, cube):
    #
    #     # Filter outlier
    #     try:
    #         ts_clean = outlier.season(ts_dek_resc, 1, param['dys_x_yr'] * param['outMax'], 3)
    #  # TODO make variable dek
    #     except (RuntimeError, Exception, ValueError):
    #         print('Outlier research failed, in position:{0}'.format(coord))
    #         err_table['outlier'] = 1
    #         return err_table
    #
    #     try:
    #         if ts_clean is not None:
    #             # ts_cleaned = outlier.fillerSeason(ts_clean)
    #             ts_cleaned = ts_clean.interpolate()
    #         else:
    #             # No season strategy
    #             summedYrs = ts_dek_resc.groupby(ts_dek_resc.index.year).sum()
    #             season_lng = 9999
    #
    #             ts_table = {'sd': pd.Series(0, index=np.unique(ts_dek.index.year)),
    #                         'ed': pd.Series(0, index=np.unique(ts_dek.index.year)),
    #                         'sl': pd.Series(np.nan, index=np.unique(ts_dek.index.year)),
    #                         'spi': pd.Series(np.nan, index=np.unique(ts_dek.index.year)),
    #                         'si': summedYrs,
    #                         'cf': pd.Series(0, index=np.unique(ts_dek.index.year)),
    #                         'yr': pd.Series(np.unique(ts_dek.index.year), index=np.unique(ts_dek.index.year))}
    #             return pd.DataFrame(ts_table), err_table, season_lng
    #     except (RuntimeError, Exception, ValueError):
    #         print('Cleaning failed, in position:{0}'.format(coord))
    #         err_table['cleaning'] = 1
    #         return err_table
    #
    # def _season_lenght(self):
    #     # Estimate Season length
    #     try:
    #         seasons, trend = fit_seasons(ts_cleaned)
    #     except (RuntimeError, Exception, ValueError):
    #         print("Error in Season and Trend estimation, in position:{0}".format(coord))
    #         err_table['s&t_ext'] = 1
    #         return err_table
    #     # Calculate season length and expected number of season
    #     try:
    #         season_lng = len(seasons) * param['dys_mlt']
    #         expSeason = int((ts_cleaned.index.max() - ts_cleaned.index.min()) / pd.Timedelta(season_lng, unit='d'))
    #     except (RuntimeError, Exception, ValueError):
    #         print("Error! Season conversion to days failed, in position:{0}".format(coord))
    #         err_table['s2d'] = 1
    #         return err_table
    #
    #     if param['medspan'] == 0:
    #         medspan = season_lng / 7
    #     else:
    #         medspan = param['medspan']
    #
    # def _interpolate_to_daily(self):
    #     # Interpolate data to daily sample
    #     ts_d = ts_cleaned.resample('D').asfreq().interpolate(method='linear').fillna(0)
    #
    #     try:
    #         if param['smp'] != 0:
    #             # Savinsky Golet filter
    #             ps = savgol_filter(ts_d, param['medspan'], param['smp'], mode='nearest')
    #             # TODO automatic selection of savgol window
    #             ps = pd.Series(ps, ts_d.index)
    #         else:
    #             ps = ts_d.rolling(medspan // 2 * 2, win_type='boxcar', center=True).mean()
    #     except (RuntimeError, Exception, ValueError):
    #         print('Error! Savinsky Golet filter problem, in position:{0}'.format(coord))
    #         err_table['savgol'] = 1
    #         return err_table
#end region

        #region prosecuzione
            # if hasattr(self, 'scale'):
            #     cube = cube * self.scale
            #
            # if hasattr(self, 'offset'):
            #     cube = cube * self.offset
            #
            # if hasattr(self, 'mask'):
            #     cube = self._no_data(cube)
            #
            # if hasattr(self, 'min') and hasattr(self, 'max'):
            #     cube = self._rescale(cube)
            # else:
            #     self.min = cube.min()
            #     self.max = cube.max()
            #     cube = self._rescale(cube)
            #
            # return cube
        #endregion

