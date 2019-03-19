# -*- coding: utf-8 -*-

import pandas as pd
import output
import os
from netCDF4 import Dataset
import numpy as np


def main():
    print('entrato!')

    import output

    pth = r'c:\Temp\temp'
    dates = pd.date_range('1998/04/01', '2013/12/21')
    dates = dates[dates.day.isin([1, 11, 21])]

    a = output.scratcher(pth, 50, 50, dates, name='test_1')

    carrier = a.variables['data']

    carrier[10, 10, 10] = 10

    a.close()

    # pth = os.path.join(pth, 'scratch.nc')
    # root = Dataset(pth, 'w', format='NETCDF4')
    #
    # xi, yi = 15680, 400
    #
    # x = root.createDimension('x', xi)
    # y = root.createDimension('y', yi)
    # time = root.createDimension('time', len(dates))
    # xv = root.createVariable('x', 'i8', ('x',))
    # yv = root.createVariable('y', 'i8', ('y',))
    # timev = root.createVariable('time', 'f8', ('time',))
    # random = root.createVariable('random', 'f8', ('x', 'y', 'time'))
    # xv[:] = np.arange(xi)
    # yv[:] = np.arange(yi)
    # timev[:] = dates
    #
    # try:
    #     for x in range(xi):
    #         if x == 100:
    #             continue
    #         random[x, :, :] = np.random.random((1, yi, len(dates)))
    #         if x in range(0, 15680, 1000):
    #             print(x)
    # except (RuntimeError, Exception):
    #     print(x)
    #
    # print('done!')
    # root.close()

    # time_idx = (pd.Timestamp('1998/04/01')-pd.Timestamp('1970/01/01'))//pd.Timedelta('1s')
    # random[1, time_idx, 1] = 2


if __name__ is '__main__':
    print('start')
    main()
    # timev[2].data.item()
    # 8.931168e+17
    # time_idx = pd.Timestamp('1998/04/01')
    # time_idx
    # Timestamp('1998-04-01 00:00:00')
    # timev[2].data.item()
    # 8.931168e+17
    # pd.to_datetime(timev[2].data.item())
    # Timestamp('1998-04-21 00:00:00')
    # time_idx = pd.Timestamp('1998/04/01 00:00:00')
    # time_idx
    # Timestamp('1998-04-01 00:00:00')
    # time_idx = (pd.Timestamp('1998/04/01 00:00:00')-pd.Timestamp('1970/01/01 00:00:00'))//pd.Timedelta('1s')
    # time_idx
    # 891388800
    # time_idx = pd.Timestamp('1998/04/01 00:00:00')
    # time_idx.asm8
    # numpy.datetime64('1998-04-01T00:00:00.000000000')
    # time_idx.astype('f8')