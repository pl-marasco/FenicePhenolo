from netCDF4 import Dataset
import os


class ScratchFile(object):
    def __init__(self, param, **kwargs):
        """
        Create a netCDF file to be used as memory dump for the pixeldrill analysis.

        :param param: configuration parameters object
        :param kwargs: name of the object
        """

        pth = os.path.join(param.scratch_pth, '.'.join((kwargs.pop('name', 'scratch'), 'nc')))
        self.root = Dataset(pth, 'w', format='NETCDF4')

        row = self.root.createDimension(param.row_nm, None)
        col = self.root.createDimension(param.col_nm, None)
        dim = self.root.createDimension(param.dim_nm, len(param.dim_val))
        self.row_v = self.root.createVariable(param.row_nm, 'f8', (param.row_nm,))
        self.col_v = self.root.createVariable(param.col_nm, 'f8', (param.col_nm,))
        self.dim_v = self.root.createVariable(param.dim_nm, 'f8', (param.dim_nm,))
        self.container = self.root.createVariable('container', 'f8', (param.row_nm, param.col_nm, param.dim_nm))

        self.row_v[:] = param.row_val
        self.col_v[:] = param.col_val
        self.dim_v[:] = param.dim_val
