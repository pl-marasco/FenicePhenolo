# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import logging

logger = logging.getLogger(__name__)


class PixelDrill(object):
    """
    Pixel drill:

    a single pixel drill rappresent the minimum unit of analye
    """
    def __init__(self, ts, px):
        self.ts_raw = ts
        self.position = px
        self.tst = None
        self.ts_filtered = None
        self.ts_interpolated = None
        self.season_ts = None
        self.season_len = None
        self.expSeason = None
        self.trend_ts = None
        self.medspan = None
        self.ts_d = None
        self.ts_sv = None
        self.mpd_val = None
        self.pks = None
        self.seasons = None
        self.sincys = []
        self.error = None
        self.phen = []


class SingularCycle(object):
    def __init__(self, ts, sd, ed):
        """
        Rappresent a singular cycle defined as the curve between two minima

        Attributes:
            mms: minimum minimum time series
            sd: Start date - MBD
            ed: End date - MED
            mml: Cycle lenght in days
            sb: Standing biomas
            mpf: permanent fration
            mpi: permanent fration
            vox: Values between two minima substracted the permanet fraction
            voxi: integral
            cbc: cycle baricenter / ex season baricenter
            csd: cycle deviation standard / Season deviation standard
            cdsdcycle: deviation standard in days /Season deviation standard in days
            ref_yr: reference yr

        :param mms: Time series as pandas.Series object
        """

        self.sd = sd  # Start date - MBD
        self.ed = ed  # End date - MED
        self.mml = self._time_delta(self.sd, self.ed)  # Cycle lenght in days
        self.mms = ts[sd:ed]  # minimum minimum time series
        self.td = self.mml*2/3  # time delta
        self.mms_b = ts[sd-self.td:ed+self.td]
        self.sb = self._integral(self.mms)  # Standing biomas
        self.mpf = self._min_min_line(self.mms)  # permanent fration
        self.mpi = self._integral(self.mpf)  # permanent fration integral
        self.vox = self._difference(self.mms, self.mpf)  # Values between two minima substracted the permanet fraction
        self.voxi = self._integral(self.vox)  # integral of vox
        self.cbc = self._baricenter()  # cycle baricenter / ex season baricenter
        self.cbcd = self._to_gregorian_date(self.cbc)
        self.csd = self._cycle_deviation_standard()  # cycle deviation standard / Season deviation standard
        self.csdd = self._to_gregorian(self.csd)  # cycle deviation standard in days /Season deviation standard in days
        self.ref_yr = self.cbcd.year  # reference yr

        self.sfs = None
        self.mas = None
        self.unx_sbc = None

    @staticmethod
    def _time_delta(sd, ed):
        """Minimum minimum lenght"""
        return ed - sd

    @staticmethod
    def _integral(ts):
        """Return the integral of a time series"""
        return ts.sum()
        # if s > 0:
        #     return s
        # else:
        #     raise ValueError('Cycle integral is negative')

    @staticmethod
    def _min_min_line(ts):
        """Interpolated line between two min and give back a time series"""
        pf = ts.copy()
        pf.iloc[1:-1] = np.nan
        return pf.interpolate()

    @staticmethod
    def _difference(crv_1, crv_2):
        """Return the differences between two time series"""
        if crv_2.sum() > 0:
            out = crv_1 - crv_2
        else:
            out = crv_1 + crv_2
        return out

    @staticmethod
    def _to_gregorian_date(value):
        """Convert to pandas date format"""
        try:
            if value is not None:
                return pd.to_datetime(value, unit='s')
            else:
                raise ValueError('date value is null')
        except (RuntimeError, Exception, ValueError):
            logger.debug('Warning! Datetime conversion went wrong')

    @staticmethod
    def _to_gregorian(value):
        try:
            return pd.Timedelta(value, unit='s')
        except (RuntimeError, Exception, ValueError):
            logger.debug('Warning! Datetime conversion went wrong')

    def _baricenter(self):
        """Barycenter"""
        cbc = 0
        try:
            index = self.vox.index
            self.posix_time = [np.int64(i.timestamp()) for i in index]
            cbc = (self.posix_time * self.vox).sum() / self.vox.sum()
        except(RuntimeError, Exception, ValueError):
            logger.debug('Warning! Baricenter calculation whent wrong in reference')
        if cbc > 0:
            return cbc
        else:
            logger.debug('Warning! Baricenter has a negative value')
            raise Exception()

    def _cycle_deviation_standard(self):
        """Season deviation standard"""
        try:
            sd = np.sqrt((np.square(self.posix_time) * self.vox).sum()
                         / self.vox.sum() - np.square(self.cbc))
            if not np.isnan(sd):
                return sd
            else:
                logger.debug('Cycle standard deviation is Nan')
                raise ValueError('Cycle standard deviation is Nan')
        except ValueError:
            logger.debug('Warning! Season deviation standard failed')
