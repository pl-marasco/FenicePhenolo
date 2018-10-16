# -*- coding: utf-8 -*-
# !/usr/bin/env python

import argparse, sys, importlib, time, os, logging
from datetime import datetime
import configurator
import reader
import preprocessing as pp
import phenolo


logger = logging.getLogger(__name__)


def main(prmts):

    start_time = time.process_time()
    _log_info(logging.getLogger('paramters'), prmts)

    # read
    cube = reader.ingest(prmts)

    prmts.x_nm, prmts.y_nm, prmts.t_nm = reader._coord_names(cube)

    if len(cube.shape) is not 1:

        ppo = pp.PreProcessor(prmts)

        pp_cube = ppo.analyse(cube)

        # ph = phenolo.analyse(pp_cube)
        #
        # ph.to_netcdf(prmts.outFilePth)

    else:
        phen = phenolo.single_px(cube.to_series(), prmts)

    end = time.process_time() - start_time

    logger.info('Process ended @{} after a total time of {}'.format(datetime.now(), end))

    # cube_msk = cube.where(~cube.isin([251, 254, 255]))

    # Calculate a list of yrs used to accumulate data

    # yrs = np.unique(cube_msk.get_index('time').year)
    # yrs_date = list(map(lambda yr_val: datetime(yr_val, 1, 1), yrs))
    # yrs_ux = date2num(yrs_date, units="seconds since 1970-01-01 00:00:00.0", calendar="gregorian")

    # dys_mlt, dys_x_yr = dek.day_calc(args.dek)

    # if singleton is None:
    #
    #     # Memory dump
    #     ds, sl, spi, si, cf, sd, ed, sns = output.create(args.out, dtst, yrs)
    #
    #     sd_pnl = pd.Panel(np.NaN, items=yrs, major_axis=range(0, xsize), minor_axis=range(0, ysize))  # start season
    #     ed_pnl = pd.Panel(np.NaN, items=yrs, major_axis=range(0, xsize), minor_axis=range(0, ysize))  # end season
    #     sl_pnl = pd.Panel(np.NaN, items=yrs, major_axis=range(0, xsize), minor_axis=range(0, ysize))  # season lngth
    #     spi_pnl = pd.Panel(np.NaN, items=yrs, major_axis=range(0, xsize), minor_axis=range(0, ysize))  # ss prnt intgr
    #     si_pnl = pd.Panel(np.NaN, items=yrs, major_axis=range(0, xsize), minor_axis=range(0, ysize))  # season intgr
    #     cf_pnl = pd.Panel(np.NaN, items=yrs, major_axis=range(0, xsize), minor_axis=range(0, ysize))  # cyclic fract
    #
    #     seasons_pnl = pd.Panel(np.NaN, items=[0, ], major_axis=range(0, xsize), minor_axis=range(0, ysize))  # seasons
    #
    #     # Error pannel
    #     errIndex = ['read', 'fix', 'rescaling', 'outlier', 'cleaning', 's&t_ext', 's2d', 'savgol', 'valley']
    #     err_pnl = pd.Panel(np.NaN, items=errIndex, major_axis=range(0, xsize), minor_axis=range(0, ysize))
    #
    #     # Pixel cicle
    #     for x in range(0, xsize):
    #         for y in range(0, ysize):
    #
    #             coord = ((x, x + 1), (y, y + 1))
    #
    #             ts_table, err_table, season = phen.analyzes(dtst, prmts, coord)
    #             if ts_table is None and season is None and err_table is not None:
    #                 err_pnl.loc[:, x, y] = err_table.values
    #             elif ts_table is not None:
    #                 # Write results into panels
    #                 try:
    #
    #                     # Populate results arrays
    #                     sd_pnl.loc[:, x, y] = ts_table.loc[:, 'sd'].astype(np.int64) // 10 ** 9
    #                     ed_pnl.loc[:, x, y] = ts_table.loc[:, 'ed'].astype(np.int64) // 10 ** 9
    #
    #                     sl_pnl.loc[:, x, y] = ts_table.groupby(['yr']).sum().loc[:, 'sl']
    #                     spi_pnl.loc[:, x, y] = ts_table.groupby(['yr']).sum().loc[:, 'spi']
    #                     si_pnl.loc[:, x, y] = ts_table.groupby(['yr']).sum().loc[:, 'si']
    #                     cf_pnl.loc[:, x, y] = ts_table.groupby(['yr']).sum().loc[:, 'cf']
    #
    #                     err_pnl.loc[:, x, y] = err_table.values
    #
    #                     seasons_pnl.loc[0, x, y] = season
    #
    #                 except (RuntimeError, Exception):
    #                     print('Panel error: x = {0} y= {1}'.format(x, y))
    #             else:
    #                 continue
    #
    #             print('Pixel in x: {0} y: {1} done'.format(x, y))
    #
    #     sl[:, :, :] = sl_pnl.values
    #     spi[:, :, :] = spi_pnl.values
    #     si[:, :, :] = si_pnl.values
    #     cf[:, :, :] = cf_pnl.values
    #
    #     # err[:, :, :] = err_pnl.values
    #
    #     sns[:, :] = seasons_pnl[0, :, :].values
    #
    #     # Close the source
    #     ds.close()
    #
    #     # Close the output
    #     dtst.close()
    #
    # else:
    #     x, y = singleton[0], singleton[1]
    #
    #     window = ((x, x + 1), (y, y + 1))
    #     ts_table = phen.analyzes(cube_msk, prmts, window)
    #
    # endregion


def _log_info(logger, prmts):

    logger.debug('-------------------- start values --------------------')
    for key, value in prmts.__dict__.items():
        logger.debug('{} = {}'.format(key, value))
    logger.debug('--------------------  end values  --------------------')


if __name__ == '__main__':

    assert sys.version_info[:2] >= (3, 6), "You need at minimum python 3.6 to execute this script"

    modules = ['numpy', 'pandas', 'xarray', 'rasterio', 'netCDF4', 'scipy', 'pyhdf', 'seasonal']

    for module in modules:
        assert importlib.util.find_spec(module), "You need {0} module".format(module)

    # Options
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-c', '--conf', help='Configuration file position', required=True)
        parser.add_argument('-p', '--plot', action='store_true', help='Activate the plotting mode', required=False)
        parser.add_argument('-l', '--log', nargs='?', const=2, type=int, help='Activate the log', required=False)
        args = parser.parse_args()

        if (not args.conf) or (not os.path.isfile(os.path.normpath(args.conf))):
            parser.error("the -c argument for the config file is missing or the given path doesn't exist!!!")

        if args.log:
            logpth = os.path.join(os.path.dirname(args.conf), 'log.txt')
            log = logging.basicConfig(filename=logpth,
                                      level=args.log*10,
                                      filemode='w')
        logger = logging.getLogger(__name__)
        logger.info('*** Phenolo 2.0 ***')
        logger.info('Process started @ {}'.format(datetime.now()))

        prmts = configurator.ProjectParameters(path=args.conf, type='ini')  # TODO  type 'ini' must be flexible

        main(prmts)

    except KeyboardInterrupt:
        print('Process killed')
        raise
    except Exception:
        print('Exception occurred')
        raise


