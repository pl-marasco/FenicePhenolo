# -*- coding: utf-8 -*-
# !/usr/bin/env python

from netCDF4 import Dataset, date2num
from datetime import datetime
import numpy as np
import os
import pandas as pd


def create(path, orig_ds, yrs_in):

    # calculate the center of pixels
    delta_x = ((orig_ds.bounds.right - orig_ds.bounds.left) / orig_ds.width)/2
    delta_y = ((orig_ds.bounds.top - orig_ds.bounds.bottom) / orig_ds.height)/2

    # Create an array of coordinates
    lon_val = np.linspace(orig_ds.bounds.left+delta_x, orig_ds.bounds.right+delta_x, orig_ds.width, endpoint=False)
    lat_val = np.linspace(orig_ds.bounds.top-delta_y, orig_ds.bounds.bottom-delta_y, orig_ds.height, endpoint=False)

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
    sd_int = root_ds.createVariable('Start day of the season', 'i4', ('time', 'latitude', 'longitude'),
                                    zlib=True)
    ed_int = root_ds.createVariable('End date of the season', 'i4', ('time', 'latitude', 'longitude'),
                                    zlib=True)
    sl_int = root_ds.createVariable('Season lenght', 'f4', ('time', 'latitude', 'longitude'),
                                    least_significant_digit=2,
                                    zlib=True)
    spi_int = root_ds.createVariable('Season permanent integral', 'f4', ('time', 'latitude', 'longitude'),
                                     least_significant_digit=2,
                                     zlib=True)
    si_int = root_ds.createVariable('Season integral', 'f4', ('time', 'latitude', 'longitude'),
                                    least_significant_digit=2,
                                    zlib=True)
    cf_int = root_ds.createVariable('Cyclic fraction', 'f4', ('time', 'latitude', 'longitude'),
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
    return root_ds, sl_int, spi_int, si_int, cf_int, sd_int, ed_int, sns_int


def scratch(pth, dim_x, dim_y, dim_t, timeidx):

    pth = os.path.join(pth, 'scratch.nc')
    root = Dataset(pth, "a", format="NETCDF4")

    x = root.createDimension('x', dim_x)
    y = root.createDimension('y', dim_y)
    time = root.createDimension('time', dim_t)

    xv = root.createVariable('x', 'i32', ('x',))
    yv = root.createVariable('y', 'i32', ('y',))
    timev = root.createVariable('time', 'f32', ('time',))
    random = root.createVariable('random', 'f32', ('y', 'time', 'x'))

    xv[:] = np.arange(0, len(dim_x))
    yv[:] = np.arange(0, len(dim_y))
    time[:] = pd.DatetimeIndex(timeidx).astype(np.int64)/1e6





