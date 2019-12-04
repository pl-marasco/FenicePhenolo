# -*- coding: utf-8 -*-

import os
from datetime import datetime

import numpy as np
import pandas as pd
from netCDF4 import Dataset, date2num


def create(path, orig_ds, yrs_in):
    """
    Create a netCDF file to be used as memory dump for the pixeldrill analysis.

    :param param: configuration parameters object
    :param kwargs: name of the object
    """

    # calculate the center of pixels
    delta_x = ((orig_ds.bounds.right - orig_ds.bounds.left) / orig_ds.width) / 2
    delta_y = ((orig_ds.bounds.top - orig_ds.bounds.bottom) / orig_ds.height) / 2

    # Create an array of coordinates
    lon_val = np.linspace(orig_ds.bounds.left + delta_x, orig_ds.bounds.right + delta_x, orig_ds.width, endpoint=False)
    lat_val = np.linspace(orig_ds.bounds.top - delta_y, orig_ds.bounds.bottom - delta_y, orig_ds.height, endpoint=False)

    # Create the list of dates
    yrs_out = list(map(lambda x_val: datetime(x_val, 1, 1), yrs_in))

    # Create the new file
    root_ds = Dataset(path, 'w', format="NETCDF4", clobber=True)

    # Definition of the convention
    root_ds.Conventions = 'CF-1.6'

    root_ds.geospatial_lat_min = np.float32(lat_val.min())
    root_ds.geospatial_lat_max = np.float32(lat_val.max())
    root_ds.geospatial_lon_min = np.float32(lon_val.min())
    root_ds.geospatial_lon_max = np.float32(lon_val.max())
    root_ds.geospatial_lat_units = 'degrees_north'
    root_ds.geospatial_lon_units = 'degrees_east'

    # Dimension
    w = orig_ds.width
    h = orig_ds.height
    root_ds.createDimension('time', None)
    root_ds.createDimension('latitude', h)
    root_ds.createDimension('longitude', w)

    # Coordinates
    time_var = root_ds.createVariable('time', 'i4', ('time',))
    lat_var = root_ds.createVariable('latitude', 'f4', ('latitude',))
    lat_var.actual_range = np.float64([lat_val[-1], lat_val[0]])
    lon_var = root_ds.createVariable('longitude', 'f4', ('longitude',))
    lon_var.actual_range = np.float64([lon_val[-1], lon_val[0]])

    # Coordinates attribute
    lat_var.long_name = "Latitude"
    lat_var.standard_name = "latitude"
    lat_var.units = "degrees_north"
    lat_var.coordinates = 'latitude'
    lat_var.axis = 'Y'
    lat_var.valid_min = np.float32(lat_val.min())
    lat_var.valid_max = np.float32(lat_val.max())

    lon_var.long_name = "Longitude"
    lon_var.standard_name = "longitude"
    lon_var.units = "degrees_east"
    lon_var.coordinates = 'longitude'
    lon_var.axis = 'X'
    lon_var.valid_min = np.float32(lon_val.min())
    lon_var.valid_max = np.float32(lon_val.max())

    time_var.long_name = 'time'
    time_var.units = "seconds since 1970-01-01 00:00:00.0"
    time_var.calendar = "gregorian"
    time_var.axis = 'T'

    # Assign the values to the dimensions
    lon_var[:] = lon_val
    lat_var[:] = lat_val
    time_var[:] = date2num(yrs_out, units=time_var.units, calendar=time_var.calendar)

    # Create the variables
    sbw_int = root_ds.createVariable('Start_week_of_the_season', 'i4', ('time', 'latitude', 'longitude'),
                                     zlib=True)
    sew_int = root_ds.createVariable('End_week_of_the_season', 'i4', ('time', 'latitude', 'longitude'),
                                     zlib=True)
    sl_int = root_ds.createVariable('Season_lenght', 'f4', ('time', 'latitude', 'longitude'),
                                    least_significant_digit=2,
                                    zlib=True)
    spi_int = root_ds.createVariable('Season_permanent_integral', 'f4', ('time', 'latitude', 'longitude'),
                                     least_significant_digit=2,
                                     zlib=True)
    si_int = root_ds.createVariable('Season_integral', 'f4', ('time', 'latitude', 'longitude'),
                                    least_significant_digit=2,
                                    zlib=True)
    cf_int = root_ds.createVariable('Cyclic_fraction', 'f4', ('time', 'latitude', 'longitude'),
                                    least_significant_digit=2,
                                    zlib=True)
    sns_int = root_ds.createVariable('season lenght', 'i4', ('latitude', 'longitude'),
                                     least_significant_digit=2,
                                     zlib=True)

    # Attributes
    root_ds.title = 'Phen'
    root_ds.description = 'Phenolo results'
    root_ds.srs = 'GEOGCS["WGS 84",DATUM["WGS_1984",' \
                  'SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],' \
                  'AUTHORITY["EPSG","6326"]],' \
                  'PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],' \
                  'UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],' \
                  'AUTHORITY["EPSG","4326"]]'

    sl_int.coordinates = 'time latitude longitude '
    spi_int.coordinates = 'time latitude longitude'
    si_int.coordinates = 'time latitude longitude'
    cf_int.coordinates = 'time latitude longitude'

    # output the file and the variables
    return root_ds, sl_int, spi_int, si_int, cf_int, sbw_int, sew_int, sns_int


class OutputCointainer(object):
    """
    Create a netCDF file to be used as memory dump for the pixeldrill analysis.

    :param param: configuration parameters object
    :param kwargs: name of the object
    """

    def __init__(self, cube, param, **kwargs):
        pth = os.path.join(param.outFilePth, '.'.join((kwargs.pop('name', ''), 'nc')))
        self.root = Dataset(pth, 'w', format='NETCDF4')

        row = self.root.createDimension(param.row_nm, len(param.row_val))
        col = self.root.createDimension(param.col_nm, len(param.col_val))
        dim = self.root.createDimension(param.dim_nm, len(self._yrs_reducer(param.dim_val)))

        self.row_v = self.root.createVariable(param.row_nm, 'f8', (param.row_nm,))
        self.col_v = self.root.createVariable(param.col_nm, 'f8', (param.col_nm,))
        self.dim_v = self.root.createVariable(param.dim_nm, 'f8', (param.dim_nm,))

        self.stb = self.root.createVariable('StandingBiomass', 'f8', (param.dim_nm, param.row_nm, param.col_nm), zlib=True, complevel=4)
        self.mpi = self.root.createVariable('MinimumminimumPermanentIntegral', 'f8', (param.dim_nm, param.row_nm, param.col_nm), zlib=True, complevel=4)

        self.sbd = self.root.createVariable('Startdate', 'f8', (param.dim_nm, param.row_nm, param.col_nm), zlib=True, complevel=4)
        self.sed = self.root.createVariable('Enddate', 'f8', (param.dim_nm, param.row_nm, param.col_nm), zlib=True, complevel=4)

        self.sl = self.root.createVariable('SeasonLenght', 'i8', (param.dim_nm, param.row_nm, param.col_nm), zlib=True, complevel=4)
        self.spi = self.root.createVariable('SeasonalPermanentIntegral', 'f8', (param.dim_nm, param.row_nm, param.col_nm), zlib=True, complevel=4)
        self.si = self.root.createVariable('SeasonIntegral', 'f8', (param.dim_nm, param.row_nm, param.col_nm), zlib=True, complevel=4)
        self.cf = self.root.createVariable('CyclicFraction', 'f8', (param.dim_nm, param.row_nm, param.col_nm), zlib=True, complevel=4)

        self.afi = self.root.createVariable('ActiveFractionIntegral', 'f8', (param.dim_nm, param.row_nm, param.col_nm), zlib=True, complevel=4)

        self.warn = self.root.createVariable('CycleWarning', 'f8', (param.dim_nm, param.row_nm, param.col_nm), zlib=True, complevel=4)

        self.n_seasons = self.root.createVariable('NumberOfSeasons', 'i8', (param.row_nm, param.col_nm), zlib=True, complevel=4)
        self.err = self.root.createVariable('PixelCriticalError', 'i8', (param.row_nm, param.col_nm), zlib=True, complevel=4)

        self.row_v[:] = param.row_val
        self.col_v[:] = param.col_val
        self.dim_v[:] = pd.to_datetime(param.dim_val).year.unique().tolist()
        # ^^^ pd.to_datetime(pd.to_datetime(param.dim_val).year.unique(), format='%Y') ^^^

        # General attributes.
        self.spatial_ref = '''GEOGCS["WGS 84",
    DATUM["WGS_1984",
        SPHEROID["WGS 84",6378137,298.257223563,
            AUTHORITY["EPSG","7030"]],
        AUTHORITY["EPSG","6326"]],
    PRIMEM["Greenwich",0,
        AUTHORITY["EPSG","8901"]],
    UNIT["degree",0.01745329251994328,
        AUTHORITY["EPSG","9122"]],
    AUTHORITY["EPSG","4326"]]'''

    @staticmethod
    def _yrs_reducer(dim_val):
        return pd.DatetimeIndex(dim_val).year.unique()

    def close(self):
        self.root.close()


def scratch_dump(pxldrl, param):
    import os
    import pickle
    """
    Create a feather file for this parameters:
        - ts 
        - ts_resc 
        - ts_cleaned
        - ts_filtered 
        - ts_d
        - trend_d
        - ps
        - pks 
    and a pickle file for:
        - sincy 
        - phen

    :param pxldrl: a pixel drill as by phenolo 
    :return: N/A
    """
    lst = ['ts', 'ts_resc', 'ts_cleaned', 'ts_filtered', 'ts_d', 'trend_d', 'ps', 'pks']
    #
    # for nm in lst:
    #     file_name = nm + '.' + str(pxldrl.position[0]) + '_' + str(pxldrl.position[1]) + '.pickle'
    #
    #     getattr(pxldrl, nm).to_pickle(file_path)
    #
    # lst1 = ['sincy', 'phen']
    #
    # for nm in lst:
    #     file_name = nm + '.' + str(pxldrl.position[0]) + '_' + str(pxldrl.position[1]) + '.pickle'
    #     with open(file_name, 'wb') as handle:
    #         pickle.dump(getattr(pxldrl, nm), handle, protocol=pickle.HIGHEST_PROTOCOL)

    file_name = 'pxldrl_' + '.' + str(pxldrl.position[0]) + '_' + str(pxldrl.position[1]) + '.pickle'
    file_path = os.path.join(param.scratch_pth, file_name)
    with open(file_path, 'wb') as handle:
        pickle.dump(pxldrl, handle, protocol=pickle.HIGHEST_PROTOCOL)
