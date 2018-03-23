#!/usr/bin/env python
import argparse
import numpy as np
import pandas as pd
import dek
import rasterio as rs
import phen
from netCDF4 import Dataset, date2num
from datetime import datetime
import output


if __name__ == '__main__':

    entire_process = datetime.now()

    # parameters definition
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('-im', help='Input time series image file', required=True)
    parser.add_argument('-out', help='Output parameter imafe file', required=True)
    parser.add_argument('-tmd', help='Define the time domain, start-end | format: [dd/mm/yyyy-dd/mm/yyyy]',
                        required=True)
    parser.add_argument('-dek', help='Define the dekad time | format: s10, s15, s30', default='s10')
    parser.add_argument('-msk', help='Define values to be masked | format [0, 251, 252 ...]',
                        default='[0.0,251.0,252.0,253.0,254.0,255.0]')
    parser.add_argument('-ovrlp', type=int, help='Detect peaks that are at least separated by minimum peak distance '
                                                 '| expressed in % of the extimated season lenght',
                        default=75)
    parser.add_argument('-mavspan', type=int, help='Maximal moving avarage span in days', default=180)
    parser.add_argument('-mavmet', type=float, help='Power of the equation of not growing season', default=1.5)
    parser.add_argument('-rng', help='Data range | format [0,250]', default='[0,250]')
    parser.add_argument('-medspan', type=int, help='Lenght of Savitzky-Golay window', default=51)
    parser.add_argument('-smp', type=int, help='Order of Savitzky-Golay polynomial', default=4)
    parser.add_argument('-outmax', type=int, help='Maximum window multiplication value to calculate outlayer',
                        default=4)

    args = parser.parse_args()

    # General variables
    inpth = args.im

    # Convert mask values to float
    msk = None
    try:
        msk = list(map(float, args.msk.translate({ord(c): None for c in '[]'}).split(',')))
    except ValueError:
        print('Error: mask values miss written')
        exit()

    # Convert min and max values
    min_v, max_v = None, None
    try:
        min_v = float(args.rng.translate({ord(c): None for c in '[]'}).replace(" ", "").split(',')[0])
        max_v = float(args.rng.translate({ord(c): None for c in '[]'}).replace(" ", "").split(',')[1])
    except ValueError:
        print('Error in min max values')
        exit()

    # Force oddity
    medspan = None
    try:
        medspan = int(np.ceil(args.medspan) // 2 * 2 + 1)
    except ValueError:
        print('Medspan force oddity error')
        exit()

    # Dek dates
    dates = None
    try:
        dates = dek.create(args.tmd, args.dek)
    except ValueError:
        print('Error in dek creation')
        exit()

    # Image loader
    dtst, xsize, ysize, bdsIdx = None, None, None, None

    try:
        with rs.Env(CPL_DEBUG=True):
                dtst = rs.open(inpth, 'r')
                # Get image metadata info
                xsize = dtst.height
                ysize = dtst.width
                bdsIdx = dtst.indexes
    except IOError as e:
        print('\nError!\n\n{0}'.format(e))
        exit()

    # Check the alignment between dek and layers
    if dtst.count != dates.size:
        print('Error in dekad file, number of decads it\'s different from number of bands \n'
              'Number of dekads: {0} \nNumber of layers in the stack: {1}'.format(dates.size, dtst.count))
        exit()

    # Calculate a list of yrs used to accumulate data
    yrs = np.unique(dates.year)
    yrs_date = list(map(lambda yr_val: datetime(yr_val, 1, 1), yrs))
    yrs_ux = date2num(yrs_date, units="seconds since 1970-01-01 00:00:00.0", calendar="gregorian")

    dys_mlt, dys_x_yr = dek.day_calc(args.dek)

    # define the dictionary
    prmts = {
        'dektyp': args.dek,
        'dates': dates,
        'bdsIdx': bdsIdx,
        'yrs': yrs,
        'msk': msk,
        'tr': args.ovrlp,
        'mavspan': args.mavspan,
        'mavmet': args.mavmet,
        'min_v': min_v,
        'max_v': max_v,
        'medspan': medspan,
        'smp': args.smp,
        'outMax': args.outmax,
        'dys_mlt': dys_mlt,
        'dys_x_yr': dys_x_yr}

    #singleton = [775, 3227] #,None 3227,775 330 500
    #2619, 717
    singleton = [0, 0]

    if singleton is None:

        # Memory dump
        ds, sl, spi, si, cf, sd, ed, sns = output.create(args.out, dtst, yrs)

        sd_pnl = pd.Panel(np.NaN, items=yrs, major_axis=range(0, xsize), minor_axis=range(0, ysize))  # start season
        ed_pnl = pd.Panel(np.NaN, items=yrs, major_axis=range(0, xsize), minor_axis=range(0, ysize))  # end season
        sl_pnl = pd.Panel(np.NaN, items=yrs, major_axis=range(0, xsize), minor_axis=range(0, ysize))  # season lngth
        spi_pnl = pd.Panel(np.NaN, items=yrs, major_axis=range(0, xsize), minor_axis=range(0, ysize))  # ss prnt intgr
        si_pnl = pd.Panel(np.NaN, items=yrs, major_axis=range(0, xsize), minor_axis=range(0, ysize))  # season intgr
        cf_pnl = pd.Panel(np.NaN, items=yrs, major_axis=range(0, xsize), minor_axis=range(0, ysize))  # cyclic fract

        seasons_pnl = pd.Panel(np.NaN, items=[0,], major_axis=range(0, xsize), minor_axis=range(0, ysize))  # seasons

        # Error pannel
        errIndex = ['read', 'fix', 'rescaling', 'outlier', 'cleaning', 's&t_ext', 's2d', 'savgol', 'valley']
        err_pnl = pd.Panel(np.NaN, items=errIndex, major_axis=range(0, xsize), minor_axis=range(0, ysize))

        # Pixel cicle
        for x in range(0, xsize):
            for y in range(0, ysize):

                coord = ((x, x + 1), (y, y + 1))

                ts_table, err_table, season = phen.analyzes(dtst, prmts, coord)
                if ts_table is None and season is None and err_table is not None:
                    err_pnl.loc[:, x, y] = err_table.values
                elif ts_table is not None:
                    # Write results into panels
                    try:

                        # Populate results arrays
                        sd_pnl.loc[:, x, y] = ts_table.loc[:, 'sd'].astype(np.int64) // 10 ** 9
                        ed_pnl.loc[:, x, y] = ts_table.loc[:, 'ed'].astype(np.int64) // 10 ** 9

                        sl_pnl.loc[:, x, y] = ts_table.groupby(['yr']).sum().loc[:, 'sl']
                        spi_pnl.loc[:, x, y] = ts_table.groupby(['yr']).sum().loc[:, 'spi']
                        si_pnl.loc[:, x, y] = ts_table.groupby(['yr']).sum().loc[:, 'si']
                        cf_pnl.loc[:, x, y] = ts_table.groupby(['yr']).sum().loc[:, 'cf']

                        err_pnl.loc[:, x, y] = err_table.values

                        seasons_pnl.loc[0, x, y] = season

                    except (RuntimeError, Exception):
                        print('Panel error: x = {0} y= {1}'.format(x, y))
                else:
                    continue

                print('Pixel in x: {0} y: {1} done'.format(x, y))

        sl[:, :, :] = sl_pnl.values
        spi[:, :, :] = spi_pnl.values
        si[:, :, :] = si_pnl.values
        cf[:, :, :] = cf_pnl.values

        # err[:, :, :] = err_pnl.values

        sns[:, :] = seasons_pnl[0, :, :].values

        # Close the source
        ds.close()

        # Close the output
        dtst.close()

    else:
        x, y = singleton[0], singleton[1]

        window = ((x, x + 1), (y, y + 1))
        ts_table = phen.analyzes(dtst, prmts, window)

