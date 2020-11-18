# -*- coding: utf-8 -*-

import fnmatch
import glob
import logging
import os
import sys
import time

import numpy as np
import numpy.ma as ma
import pandas as pd
import xarray as xr
from pyhdf.SD import *

logger = logging.getLogger(__name__)


class Reader(object):
    def __init__(self):
        return


def _get_rasterio(prmts, dim):

    from phenolo import chronos as dk

    try:
        # dataset = xr.open_rasterio(prmts.inFilePth, chunks={'x': 500, 'y': 500})
        dataset = xr.open_rasterio(prmts.inFilePth)
    except IOError:
        logger.debug('Error reading file through RasterIO')
        sys.exit(1)

    if prmts.start_time is not pd.NaT and prmts.end_time is not pd.NaT:
        logger.debug('Dates comes from the setting file')
        time_dom = dk.create(prmts.start_time, prmts.end_time, prmts.dek)
    else:
        if 'band_names' in dataset.attrs:
            import re
            time_dom = pd.to_datetime(re.findall(r'\d\d\d\d\d\d\d\d', dataset.attrs['band_names']))
        else:
            logger.debug('Bands name doesn\'t contains dates; please add to the setting file under obs_start obs_end')
            raise sys.exit(1)

    if dict(dataset.sizes)['band'] == time_dom.size:
        dt = dataset.assign_coords(band=time_dom)
    else:
        logger.debug('Ups! we have a problem with the time size. Doesn\'t match the data lenght')
        raise sys.exit(1)

    dt = dt.rename(x='lon', y='lat', band='time')

    return _slice_cube(dt, dim)


def _get_gs_zarr(prmts, dim):
    import gcsfs

    fs = gcsfs.GCSFileSystem(project='pheno', token='cloud')
    bucket = prmts.inFilePth.replace('gs://','')
    gcsmap = gcsfs.mapping.GCSMap(bucket, gcs=fs, check=True, create=False)
    dataset = xr.open_zarr(gcsmap, mask_and_scale=False)

    return _slice_cube(dataset, dim)


def _get_netcdf(prmts, dim):
    dataset = xr.open_dataset(prmts.inFilePth,
                              mask_and_scale=False,
                              decode_times=True)

    return _slice_cube(dataset, dim)


def _get_multi_netcdf(path, dim, prmts):
    dataset = xr.open_mfdataset(path,
                                mask_and_scale=False,
                                decode_times=prmts.decode)
    return _slice_cube(dataset, dim)


def _get_hls(path):
    hdf = SD(path, SDC.READ)

    # dataset_dic = hdf.datasets()
    attrib_dic = hdf.attributes()

    if 'SPACECRAFT_NAME' in attrib_dic and ('Sentinel-2A' in attrib_dic['SPACECRAFT_NAME'] or
                                            'NONE' in attrib_dic['SPACECRAFT_NAME']):
        nir_band_nm = 'B8A'
        r_band_nm = 'B04'
        qa_band_nm = 'QA'
        scene_name = attrib_dic['TILE_ID']
        n_rows = int(attrib_dic['NROWS'])
        n_cols = int(attrib_dic['NCOLS'])
        sp_res = int(attrib_dic['SPATIAL_RESOLUTION'])
        ul_x = int(attrib_dic['ULX'])
        ul_y = int(attrib_dic['ULY'])
        lr_x = int(ul_x + n_cols * sp_res)
        lr_y = int(ul_y - n_rows * sp_res)
        s_time = attrib_dic['SENSING_TIME']
        if "+" in s_time:
            s_time = s_time.split(' + ')
            s_time = pd.to_datetime(s_time[0]) + pd.to_datetime(pd.Series([s_time[0], s_time[1]])).diff().mean() / 2

    elif 'SENSOR' in attrib_dic and 'OLI_TIRS' in attrib_dic['SENSOR']:
        nir_band_nm = 'band05'
        r_band_nm = 'band04'
        qa_band_nm = 'QA'
        if ";" in attrib_dic['LANDSAT_SCENE_ID']:
            scene_name = attrib_dic['LANDSAT_SCENE_ID'].split(';')[0]
        else:
            scene_name = attrib_dic['LANDSAT_SCENE_ID']

        n_rows = attrib_dic['NROWS']
        n_cols = attrib_dic['NCOLS']
        sp_res = attrib_dic['SPATIAL_RESOLUTION']
        ul_x = int(attrib_dic['ULX'])
        ul_y = int(attrib_dic['ULY'])
        lr_x = int(ul_x + n_cols * sp_res)
        lr_y = int(ul_y - n_rows * sp_res)

        s_time = attrib_dic['SENSING_TIME']
        if ";" in s_time:
            s_time = s_time.split(';')
            s_time = pd.to_datetime(s_time[0]) + pd.to_datetime(pd.Series([s_time[0], s_time[1]])).diff().mean() / 2
    else:
        return

    # NIR band
    nir_band = hdf.select(nir_band_nm)
    nir_add_offset, nir_scale_factor = _scale(nir_band)

    # Red band
    r_band = hdf.select(r_band_nm)
    r_add_offset, r_scale_factor = _scale(r_band)

    # QA mask application
    q_band = hdf.select(qa_band_nm)

    # x_arr, y_arr = [None] * 2

    # X/Y array for the DataArray
    x_arr = np.linspace(ul_x, lr_x, n_cols, endpoint=False)
    y_arr = np.linspace(ul_y, lr_y, n_rows, endpoint=False)

    # Nir_band
    nir = nir_band.get()
    nir = (ma.array(nir, mask=np.isin(nir, np.array([-1000]))) - nir_add_offset) * nir_scale_factor

    # Red_band
    r = r_band.get()
    r = (ma.array(r, mask=np.isin(r, np.array([-1000]))) - r_add_offset) * r_scale_factor

    # Q_mask
    q_mask = np.invert(np.isin(q_band.get(), np.array([128, 192, 64, 4, 68, 132, 196, 0])))

    r = ma.array(r, mask=q_mask)
    nir = ma.array(nir, mask=q_mask)

    # NDVI calculation and 3d dimension expantion
    ndvi = np.expand_dims(((nir - r) / (nir + r)), 2)

    # NDVI Outlayer removal
    ndvi_msk = ma.masked_outside(ndvi, -1, 1)

    data_array = xr.DataArray(ndvi_msk, coords=[y_arr, x_arr, [pd.to_datetime(s_time)]],
                              dims=['N', 'E', 'Time'], name=scene_name)

    # region old script
    # if _coord['E'] is not None:
    #     # Tests area
    #     dat = xr.DataArray(ndvi_msk, coords=[Y_arr, X_arr, [pd.to_datetime(TIME)]],
    #                        dims=['N', 'E', 'Time'], name=scene_name)
    #     # attrs=attrib_dic)
    #
    #     dat = dat.sel(E=_coord['E'], N=_coord['N'])
    # else:
    #     dat = xr.DataArray(ndvi_msk, coords=[Y_arr, X_arr, [pd.to_datetime(TIME)]],
    #                        dims=['N', 'E', 'Time'], name=scene_name)
    #     # attrs=attrib_dic)

    # if isinstance(dim['x'], slice) is False and _coord['E'] is not None:
    #     x = int((_coord['E'] - UL_X) / SP_RES)
    #     y = int((UL_Y - _coord['N']) / SP_RES)
    #
    #     # X/Y array for the DataArray
    #     X_arr = [_coord['E']]
    #     Y_arr = [_coord['N']]
    #
    #     # Nir band
    #     nir = np.array(nir_band[y, x])
    #     nir = (ma.array(nir, mask=np.isin(nir, np.array([-1000]))) - add_offset) * scale_factor
    #     # Red band
    #     r = np.array(r_band[y, x])
    #     r = (ma.array(r, mask=np.isin(r, np.array([-1000]))) - add_offset) * scale_factor
    #     # Q_mask
    #     q_band = np.array(q_band[x, y])
    #     q_mask = np.invert(np.isin(q_band, np.array([128, 192, 64, 4, 68, 132, 196, 0])))
    #
    #     # NDVI calculation and 3d dimension expantion
    #
    #     ndvi = np.atleast_3d(((nir - r) / (nir + r)))
    #
    #     # NDVI Outlayer removal
    #     ndvi_msk = ma.masked_outside(ndvi, -1, 1)
    #
    #     # Tests area
    #     dat = xr.DataArray(ndvi_msk, coords=[Y_arr, X_arr, [pd.to_datetime(TIME)]],
    #                        dims=['N', 'E', 'Time'], name=scene_name)
    #     # attrs=attrib_dic)
    #
    #
    # else:
    #     # X/Y array for the DataArray
    #     X_arr = np.linspace(UL_X, LR_X, NCOLS, endpoint=False)
    #     Y_arr = np.linspace(UL_Y, LR_Y, NROWS, endpoint=False)
    #
    #     # Nir_band
    #     nir = nir_band.get()
    #     nir = (ma.array(nir, mask=np.isin(nir, np.array([-1000]))) - add_offset) * scale_factor
    #
    #     # Red_band
    #     r = r_band.get()
    #     r = (ma.array(r, mask=np.isin(r, np.array([-1000]))) - add_offset) * scale_factor
    #
    #     # Q_mask
    #     q_mask = np.invert(np.isin(q_band.get(), np.array([128, 192, 64, 4, 68, 132, 196, 0])))
    #
    #     r = ma.array(r, mask=q_mask)
    #     nir = ma.array(nir, mask=q_mask)
    #
    #     # NDVI calculation and 3d dimension expantion
    #     ndvi = np.expand_dims(((nir - r) / (nir + r)), 2)
    #
    #     # NDVI Outlayer removal
    #     ndvi_msk = ma.masked_outside(ndvi, -1, 1)
    #
    #     if _coord['E'] is not None:
    #         # Tests area
    #         dat = xr.DataArray(ndvi_msk, coords=[Y_arr, X_arr, [pd.to_datetime(TIME)]],
    #                            dims=['N', 'E', 'Time'], name=scene_name)
    #         # attrs=attrib_dic)
    #
    #         dat = dat.sel(E=_coord['E'], N=_coord['N'])
    #     else:
    #         dat = xr.DataArray(ndvi_msk, coords=[Y_arr, X_arr, [pd.to_datetime(TIME)]],
    #                            dims=['N', 'E', 'Time'], name=scene_name)
    #         # attrs=attrib_dic)
    #
    # return dat
    # end region

    return data_array


def _get_multi_hdf(f_list, dim):
    hls_lst = [_get_hls(file) for file in f_list]

    data_array = xr.concat(hls_lst, dim='Time', **{'decode_cfbool': False})

    cube = _slice_cube(data_array, dim)

    return cube


def _get_hdf():
    # return xr.concat(HLS_lst, dim='Time')
    return


def _coord_names(data):
    dims = [i.lower() for i in data.dims]
    crd_x, crd_y, crd_t = [None] * 3

    if 'lat' in dims or 'lon' in dims:
        crd_x = data.dims[dims.index('lon')]
        crd_y = data.dims[dims.index('lat')]
        crd_t = data.dims[dims.index('time')]

    elif 'e' in dims or 'n' in dims:
        crd_x = data.dims[dims.index('e')]
        crd_y = data.dims[dims.index('n')]
        crd_t = data.dims[dims.index('time')]

    return crd_x, crd_y, crd_t


def _slice_cube(dataset, dim):
    if hasattr(dataset, 'NDVI'):
        data = dataset.NDVI
    else:
        data = dataset

    crd_x, crd_y, crd_t = _coord_names(data)

    logger.debug(f'Coordinates names are: {crd_x},{crd_y},{crd_t}')

    # no slice
    if dim['time'] == slice(None) and dim['x'] == slice(None) and dim['y'] == slice(None):
        return data
    # single pixel
    elif dim['x'].start == dim['x'].stop and dim['y'].start == dim['y'].stop:
        # Absolute values
        try:
            dc = dict([(crd_x, dim['x'].start), (crd_y, dim['y'].start)])
            if isinstance(dim['x'].start, int) and isinstance(dim['y'].start, int):
                _coord_range_check(dim, data, crd_x, crd_y, crd_t)
                itrm_cube = data[dc]
            # Relative values (prj)
            elif isinstance(dim['x'].start, float) and isinstance(dim['y'].start, float):
                _coord_range_check(dim, data, crd_x, crd_y, crd_t)
                itrm_cube = data.sel(dc, method='nearest')
            else:
                raise ValueError

            if dim['time'] == slice(pd.NaT, pd.NaT, None) or dim['time'] == slice(None, None, None):
                return itrm_cube
            else:
                return itrm_cube.sel(dict([(crd_t, dim['time'])]))

        except ValueError:
            logger.debug('Coordinates slice has created an error')
            raise sys.exit(1)
    # Area
    else:
        try:
            # Absolute
            if isinstance(dim['x'].start, int) or isinstance(dim['x'].stop, int) and \
                    isinstance(dim['y'].start, int) or isinstance(dim['y'].stop, int):
                dc = dict([(crd_x, dim['x']), (crd_y, dim['y'])])
                _coord_range_check(dim, data, crd_x, crd_y, crd_t)
                itrm_cube = data[dc]
                if dim['time'] == slice(pd.NaT, pd.NaT, None):
                    return itrm_cube
                else:
                    return itrm_cube.sel(dict([(crd_t, dim['time'])]))
            # Relative (prj)
            elif isinstance(dim['x'].start, float) or isinstance(dim['x'].stop, float) and \
                    isinstance(dim['y'].start, float) or isinstance(dim['y'].stop, float):
                _coord_range_check(dim, data, crd_x, crd_y, crd_t)
                if dim['time'] == slice(pd.NaT, pd.NaT, None):
                    return data.sel(dict([(crd_x, dim['x']), (crd_y, dim['y'])]))
                else:
                    return data.sel(dict([(crd_t, dim['time']), (crd_x, dim['x']), (crd_y, dim['y'])]))
                # TODO split time slice
        except ValueError:
            logger.debug('Coordinates slice has created an error')
            raise ValueError


def _coord_range_check(dim, data, crd_x, crd_y, crd_t):
    sizes = dict(data.sizes)
    # single pixel
    if dim['x'].start == dim['x'].stop and dim['y'].start == dim['y'].stop:
        # Absolute
        if isinstance(dim['x'].start, int) and isinstance(dim['y'].start, int):
            if not (dim['x'].start <= sizes[crd_x]):
                logger.info('Coordinates {} out of range'.format(crd_x))
                raise ValueError
            if not (dim['y'].start <= sizes[crd_y]):
                logger.info('Coordinates {} out of range'.format(crd_y))
                raise ValueError

        # Relative
        elif isinstance(dim['x'].start, float) and isinstance(dim['y'].start, float):
            if not (data.coords[crd_x].values.min() <= dim['x'].start <= data.coords[crd_x].values.max()) or \
                    not (data.coords[crd_x].values.min() <= dim['x'].stop <= data.coords[crd_x].values.max()):
                logger.info('Coordinates {} out of range'.format(crd_x))
                raise ValueError

            if not (data.coords[crd_y].values.min() <= dim['y'].start <= data.coords[crd_y].values.max()) or \
                    not (data.coords[crd_y].values.min() <= dim['y'].stop <= data.coords[crd_y].values.max()):
                logger.info('Coordinates {} out of range'.format(crd_y))
                raise ValueError

    # multi pixel
    else:
        # Absolute
        if isinstance(dim['x'].start, int) or isinstance(dim['x'].stop, int) and \
                isinstance(dim['y'].start, int) or isinstance(dim['y'].stop, int):
            if not (dim['x'].start <= sizes[crd_x]) or not (dim['x'].stop <= sizes[crd_x]):
                logger.info('Coordinates {} out of range'.format(crd_x))
                raise ValueError
            if not (dim['y'].start <= sizes[crd_y]) or not (dim['y'].stop <= sizes[crd_y]):
                logger.info('Coordinates {} out of range'.format(crd_y))
                raise ValueError

        # Relative (prj)
        elif isinstance(dim['x'].start, float) or isinstance(dim['x'].stop, float) and \
                isinstance(dim['y'].start, float) or isinstance(dim['y'].stop, float):

            if not (data.coords[crd_x].values.min() <= dim['x'].start <= data.coords[crd_x].values.max()) or \
                    not (data.coords[crd_x].values.min() <= dim['x'].stop <= data.coords[crd_x].values.max()):
                logger.info('Coordinates {} out of range'.format(crd_x))
                raise ValueError

            if not (data.coords[crd_y].values.min() <= dim['y'].start <= data.coords[crd_y].values.max()) or \
                    not (data.coords[crd_y].values.min() <= dim['y'].stop <= data.coords[crd_y].values.max()):
                logger.info('Coordinates {} out of range'.format(crd_y))
                raise ValueError

        else:
            logger.debug('Problem in the coordinate, int and float mixed')
            raise ValueError

    # time check
    if dim['time'].start is not None:
        # todo add the possibility to build the time dimension if unaviable, this should avoid the possibility to subsample
        tm_rng = pd.date_range(data.coords[crd_t].min().values, data.coords[crd_t].max().values).floor('D')

        if dim['time'].start is not None:
            if dim['time'].start.floor('D') not in tm_rng:
                logger.info('Time slice Start out of range')
                raise ValueError

        if dim['time'].stop is not None:
            if dim['time'].stop.floor('D') not in tm_rng:
                logger.info('Time slice Start out of range')
                raise ValueError


def _coord(coord):
    try:
        if coord is not '':
            if ':' in coord:
                arr = coord.translate({ord(c): None for c in '[]'}).split(':')
                # convert string value and null
                arr = list(map(lambda x: None if x == "" else float(x) if '.' in x else int(x), arr))
                return slice(arr[0], arr[1])
            else:
                arr = coord.translate({ord(c): None for c in '[]'})
                arr = float(arr) if ',' in arr or '.' in arr else int(arr)
                return arr
        else:
            return slice(None, None, None)

    except RuntimeError:
        raise Exception('Transforming coordinates went wrong')


def _scale(band_obj):
    add_offset = None
    scale_factor = None

    for key, value in band_obj.attributes().items():
        if key == 'add_offset':
            add_offset = np.float32(value)
        if key == 'scale_factor':
            scale_factor = np.float32(value)

    return add_offset, scale_factor


def _get_slicers(prmts):
    x_slice, y_slice, time_slice = [slice(None)] * 3

    if hasattr(prmts, 'exm_start') or hasattr(prmts, 'exm_end'):

        tm_sl = [None] * 2

        if type(prmts.exm_start) is not type(pd.NaT):
            tm_sl[0] = prmts.exm_start

        if type(prmts.exm_end) is not type(pd.NaT):
            tm_sl[1] = prmts.exm_end
            # pd._libs.tslibs.nattype.NaTType

        time_slice = slice(tm_sl[0], tm_sl[1])

    if hasattr(prmts, 'ext'):

        extent = prmts.ext

        if extent is not None:
            if len(extent) in [2, 4]:
                if len(extent) == 2:
                    logger.debug('Two coordinates found')
                    x_slice = slice(extent[0], extent[0])
                    y_slice = slice(extent[1], extent[1])

                elif len(extent) == 4:
                    logger.debug('Four cooridnates found')
                    if type(extent[1]) is not int:
                        if not (extent[1] > extent[3]) or not (extent[0] < extent[1]):
                            logger.debug('Coordinates doesn\'t respect the descending order, please check')
                            raise sys.exit(0)
                    x_slice = slice(extent[0], extent[2])
                    y_slice = slice(extent[1], extent[3])
            else:
                logger.debug('Error in the slice coordinates')
                raise Exception('Error in the coordinates subsample')
        else:
            logger.debug('No coordinate slice selected')
            x_slice, y_slice = [slice(None)] * 2

    logger.debug(f'Time slicer:\
        {time_slice.start}{"->" + str(time_slice.stop) if time_slice.stop is not time_slice.start else ""}\n'
            f'\t\t\t X slicer:    {x_slice.start}{" -> " + str(x_slice.stop)if x_slice.stop != x_slice.start else ""}\n'
            f'\t\t\t Y slicer:    {y_slice.start}{" -> " + str(y_slice.stop)if y_slice.stop != y_slice.start else ""}')

    return {'time': time_slice, 'x': x_slice, 'y': y_slice}


def _dasker(dataset, dim_bloks, col_bloks, row_bloks):
    crd_y, crd_x, crd_t = _coord_names(dataset)
    return dataset.chunk({crd_t: dim_bloks, crd_x: row_bloks, crd_y: col_bloks})


def ingest(prmts):
    start = time.time()

    dim = _get_slicers(prmts)

    try:
        cube = None
        # TODO adapt to GCFS file
        if os.path.isfile(prmts.inFilePth) or 'gs://' in prmts.inFilePth:
            if fnmatch.fnmatch(prmts.inFilePth, '*.nc'):
                cube = _get_netcdf(prmts, dim)
            elif fnmatch.fnmatch(prmts.inFilePth, '*.hdf'):
                cube = _get_hdf(prmts, dim)
            elif fnmatch.fnmatch(prmts.inFilePth, '*.vrt'):
                cube = _get_rasterio(prmts, dim)
            elif fnmatch.fnmatch(prmts.inFilePth, '*.img'):
                cube = _get_rasterio(prmts, dim)
            elif fnmatch.fnmatch(prmts.inFilePth, '*.zarr'):
                cube = _get_gs_zarr(prmts, dim)
            else:
                cube = _get_rasterio(prmts.inFilePth, dim)
        elif os.path.isdir(os.path.dirname(prmts.inFilePth)):
            if '*.nc' in prmts.inFilePth:
                cube = _get_multi_netcdf(prmts.inFilePth, dim, prmts)
            elif '*.hdf' in prmts.inFilePth:
                cube = _get_multi_hdf(glob.glob(os.path.join(prmts.inFilePth, '*.hdf')), dim)
            else:
                cube = _get_rasterio(prmts.inFilePth, dim)
        else:
            logger.info('File or directory not found')
            raise FileNotFoundError

        if prmts.ext is None:
            deltatime = time.time() - start
            logger.info('Loading data required:{}'.format(deltatime))
            prmts.add_dims(cube)

            if cube.chunks is None:
                col_blocks = prmts.col_val.size
                row_blocks = 1 # prmts.row_val.size
                dim_blocks = prmts.dim_val.size
                cube = _dasker(cube, dim_blocks, col_blocks, row_blocks)

            return cube

        else:
            prmts.add_dims(cube)
            deltatime = time.time() - start
            logger.info('Loading data required:{}'.format(deltatime))
            return cube

    except Exception as ex:
        raise ex