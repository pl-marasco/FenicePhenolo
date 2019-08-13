# -*- coding: utf-8 -*-

import configparser as cp
import logging
import os
import sys

import numpy as np
import pandas as pd

from phenolo import chronos

logger = logging.getLogger(__name__)


class ProjectParameters(object):

    def __init__(self, **kwargs):
        """
        Create an object containg all the parameters for the run.
        The beset approach is to create an ini fille as below descripted:

        # Configuration for Phenolo 2.0

        [GENERAL_SETTINGS]
        # input time series image file
        InFile  =
        # output in netCDF
        OutFile =
        # Scratch files (#True retain the scratch files, False overwrite over single file)
        Retain_scratch = True
        ScratchPath =
        Sensor_type = Right now the only option is Spot but if not specyfy only the specific cleaning for
                      that particular sensor isn't applayed

        [RUN_PARAMETERS_INPUT]
        # time span in format dd/mm/yyyy,dd/mm/yyyy
        obs_start             =
        obs_end               =
        # analysis window
        exm_start             =
        exm_end               =
        area                  =
        # Single point (E, N) , Area (Top left E, N ; Bottom right E, N)
        extent                =
        # temporal range (decad) of input values (s10, s15, s30, ...)
        dek                   = s10
        # data range of values (comma separated 2 values max)
        rng                   = 0, 250
        scale                 =
        offset                =
        # data values to be used during masking (comma separated values)
        msk                   = 251, 255
        cloud                 = 252
        snow                  = 253
        sea                   = 254

        [RUN_PARAMETERS_FILTER]
        mad_power = 1.5

        [RUN_PARAMETERS_SEGMENTATION]
        # Detect peaks that are at least separated by the minimum peak distance, expressed in % of the estimated season length
        ovrlp = 75
        # Maximal moving avarage span in days
        mavspan = 180
        # Power of equation of not growing season
        mavmet = 1.5

        [RUN_PARAMETERS_SMOOTH]
        # length of Savitzky-Golay window
        medspan               = 51
        # order of Savitzky-Golay polynomial
        smp                   = 3
        #Maximum window multiplication value to calculate outlayer
        outmax                = 4

        :param kwargs:
        """

        if 'type' in kwargs and 'path' in kwargs:
            if kwargs['type'] == 'ini':
                logger.debug('Options are structured in .ini file')
                path_config = kwargs['path']

                if not os.path.isfile(path_config):
                    logger.debug("Config file don't exists or is in the wrong location!{}".format(path_config))
                    print("Config file don't exists or is in the wrong location!! I quit. Check: " + str(path_config))
                    sys.exit(0)

                config = cp.ConfigParser()
                config.read(path_config)

                # [GENERAL_SETTINGS]
                section = 'GENERAL_SETTINGS'
                self.inFilePth = self.__read(config, section, 'in_file')

                if self.__read(config, section, 'out_file') != '':
                    read_out_pth = self.__read(config, section, 'out_file')

                    root, file = os.path.split(read_out_pth)
                    if file == '':
                        file = os.path.split(self.inFilePth)[1]
                        file_nm = os.path.splitext(file)[0]
                    else:
                        file_nm = os.path.splitext(file)[0]

                    self.outFilePth = root
                    self.outName = file_nm
                else:
                    root, file_nm = os.path.split(self.inFilePth)
                    self.outFilePth = root
                    self.outName = file_nm

                if self.__read(config, section, 'retain_scratch').lower() == 'true':
                    self.ovr_scratch = True
                else:
                    self.ovr_scratch = False
                self.scratch_pth = self.__read(config, section, 'scratch_path')

                self.sensor_typ = self.__read(config, section, 'sensor_type').lower()

                self.decode = self.__read(config, section, 'data_decode').lower()

                # [INFRASTRUCTURE_PARAMETERS]
                section = 'INFRASTRUCTURE_PARAMETERS'
                if self.__read(config, section, 'processes').lower() == 'true':
                    self.processes = True
                else:
                    self.processes = False

                if self.__read(config, section, 'n_workers') is not '':
                    self.n_workers = self.__read(config, section, 'n_workers', type='int')
                else:
                    self.n_workers = None

                if self.__read(config, section, 'threads_per_worker') is not '':
                    self.threads_per_worker = self.__read(config, section, 'threads_per_worker', type='int')
                else:
                    self.threads_per_worker = None

                # [RUN_PARAMETERS_INPUT]
                # Time dimension
                section = 'RUN_PARAMETERS_INPUT'
                self.start_time = self.__read(config, section, 'obs_start', type='time')
                self.end_time = self.__read(config, section, 'obs_end', type='time')
                self.exm_start = self.__read(config, section, 'exm_start', type='time')
                self.exm_end = self.__read(config, section, 'exm_end', type='time')
                self.exm_end = self.__read(config, section, 'exm_end', type='time')

                # area
                self.area = self.__read(config, section, 'area')
                # extent
                self.typExt = 'Full'
                ext = self.__read(config, section, 'extent', type='coord')

                if ext is not None and (len(ext) in [1, 3] or len(ext) > 4):
                    logger.debug('Error: Founded [{}] coordinates'.format(len(ext)))
                    raise exit(1)

                self.ext = ext

                # dekad
                dek = self.__read(config, section, 'dek')
                if dek in ['s5', 's10', 's15', 's30']:
                    self.dek = dek
                    self.yr_dys = chronos.day_calc(self.dek)[0]
                    self.yr_dek = chronos.day_calc(self.dek)[1]

                else:
                    print("Dekad type unrecognised, please check: " + str(dek))
                    sys.exit(0)

                # range, scale, offset
                rng = self.__read(config, section, "rng", type='list')  # lan= can be used to define the minimum lenght
                try:
                    if rng is not None:
                        self.min = float(rng[0])
                        self.max = float(rng[1])
                    else:
                        self.min = None
                        self.max = None
                except ValueError:
                    print('Error in range min, max values')
                    sys.exit(0)

                self.scale = self.__read(config, section, "scale", type='float')
                self.offset = self.__read(config, section, "offset", type='float')

                # validity range

                self.qt = self.__read(config, section, "threshold", type='float')
                if self.qt is None:
                    self.qt = 0.25
                self.min_th = self.__read(config, section, "min_th", type='float')
                if self.min_th is None:
                    self.min_th = self.min
                self.max_th = self.__read(config, section, "max_th", type='float')
                if self.max_th is None:
                    self.max_th = self.max

                # mask, cloud, snow, sea
                self.mask = self.__read(config, section, "msk", type='list')
                self.cloud = self.__read(config, section, "cloud", type='int')
                self.snow = self.__read(config, section, "snow", type='int')
                self.sea = self.__read(config, section, "sea", type='int')

                # [RUN_PARAMETERS_FILTER]
                section = 'RUN_PARAMETERS_FILTER'
                self.mad_pwr = self.__read(config, section, "mad_power", type='float')

                # [RUN_PARAMETERS_SEGMENTATION]
                section = 'RUN_PARAMETERS_SEGMENTATION'
                self.ovrlp = self.__read(config, section, "ovrlp", type='int')
                self.mavspan = self.__read(config, section, "mavspan", type='int')
                self.mavmet = self.__read(config, section, "mavmet", type='float')

                # [RUN_PARAMETERS_SMOOTH]
                section = 'RUN_PARAMETERS_SMOOTH'
                self.medspan = self.__read(config, section, "medspan", type='int')
                try:
                    self.medspan = int(np.ceil(self.medspan) // 2 * 2 + 1)
                except ValueError:
                    print('Medspan force oddity error')
                    sys.exit(0)

                self.smp = self.__read(config, section, "smp", type='int')
                self.outmax = self.__read(config, section, "outmax", type='int')

                self.row_nm, self.col_nm, self.dim_nm, = [None] * 3
                self.row_val, self.col_val, self.dim_val = [None] * 3
                self.pixel_list = None

                return
            elif kwargs['type'] == 'CopernicusNetCDF':
                pass
                # path_config = kwargs['path']
                # TODO create the attribute reading to populate min max ...

            elif kwargs['type'] == 'Datacube':
                pass
                # path_config = kwargs['path']
                # TODO create the attribute reading to populate min max ...

        else:
            logger.debug('Default parameters loaded')
            self.ovrlp = 75
            self.mavspan = 180
            self.mavmet = 1.5
            self.medspan = 51
            self.smp = 4
            self.outmax = 5

    @staticmethod
    def __read(config, section, parameter, **kwargs):

        if 'type' in kwargs:
            if kwargs['type'] == 'float':
                try:
                    return config.getfloat(section, parameter)
                except (Exception, TypeError):
                    return None

            elif kwargs['type'] == 'int':
                try:
                    return config.getint(section, parameter)
                except(Exception, TypeError):
                    return None

            elif kwargs['type'] == 'list':
                if 'len' in kwargs:
                    min_length = kwargs['len']
                else:
                    min_length = 0

                values = config.get(section, parameter).split(',')

                if values == ['']:
                    return None
                elif len(values) < min_length:
                    print("Config file does not contain list of values!! I quit. Check: ")
                    sys.exit(0)
                else:
                    if '.' in values[0]:
                        readed_list = list(map(float, values))
                    else:
                        readed_list = list(map(int, values))
                return readed_list

            elif kwargs['type'] == 'coord':
                if config.get(section, parameter) is not '':
                    values = config.get(section, parameter)
                    charachter = {' ': '', ';': ',', '-': ',', ':': ','}
                    for i, j in charachter.items():
                        values = values.replace(i, j)
                    values = values.split(',')
                    if '.' in values[0]:
                        readed_list = list(map(float, values))
                    else:
                        readed_list = list(map(int, values))
                    return readed_list
                else:
                    return None

            elif kwargs['type'] == 'time':
                time = config.get(section, parameter)
                if time is not '':
                    return pd.to_datetime(config.get(section, parameter), format='%d/%m/%Y')
                else:
                    return pd.to_datetime(config.get(section, parameter))
            else:
                return config.get(section, parameter)
        else:
            return config.get(section, parameter)

    @staticmethod
    def __coord_names(data):

        dims = [i.lower() for i in data.dims]
        crd_row, crd_col, crd_dim = [None] * 3

        if 'lat' in dims or 'lon' in dims:
            crd_col = data.dims[dims.index('lon')]
            crd_row = data.dims[dims.index('lat')]
            crd_dim = data.dims[dims.index('time')]

        elif 'e' in dims or 'n' in dims:
            crd_col = data.dims[dims.index('e')]
            crd_row = data.dims[dims.index('n')]
            crd_dim = data.dims[dims.index('time')]

        elif 'time' in dims:
            crd_dim = data.dims[dims.index('time')]

        return crd_row, crd_col, crd_dim

    def add_dims(self, data):
        self.row_nm, self.col_nm, self.dim_nm = self.__coord_names(data)

        if self.row_nm is not None and self.col_nm is not None:
            self.row_val = data.coords[self.row_nm].data
            self.col_val = data.coords[self.col_nm].data
            self.dim_val = data.coords[self.dim_nm].data
        else:
            self.dim_val = data.coords[self.dim_nm].data

    def add_px_list(self, cube):
        # TODO add other sensors or structures
        #  Verify the stability over the time period of water classification

        # Create a list of pixels to be analyzed
        med = cube.quantile(0.80, dim=self.dim_nm)
        # med = cube.median(dim=self.dim_nm)

        mask = np.append(self.sea, self.mask)
        masked = np.ma.MaskedArray(np.ma.MaskedArray(med, np.in1d(med, mask)))
        self.pixel_list = np.argwhere(masked)
