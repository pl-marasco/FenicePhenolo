# -*- coding: utf-8 -*-
# !/usr/bin/env python

import numpy as np
import pandas as pd
from IPython.core.magic_arguments import kwds
from seasonal import fit_seasons, adjust_seasons, fit_trend
from statsmodels.robust import mad
from scipy.signal import savgol_filter


class MAD(object):
    """
    Based on:

    http://eurekastatistics.com/using-the-median-absolute-deviation-to-find-outliers/

    Leys, C., et al.,
    Detecting outliers: Do not use standard deviation around the mean, use absolute deviation around the median,
    Journal of Experimental Social Psychology, Volume 49, Issue 4, July 2013, pp. 764-766. *

    Rousseeuw, P.J. and Croux C. (1993)
    Alternatives to the Median Absolute Deviation,
    Journal of the American Statistical Association, December 1993, pp. 1273-1283.
    """

    def __init__(self, mad_pwr=2.575, **kwargs):

        self.median = None
        self.mad_pwr = mad_pwr
        self.c = kwargs.pop('c',  1.4826)
        self.min_w = kwargs.pop('min_w',  1)
        self.max_w = kwargs.pop('max_w', None)
        self.left = None
        self.right = None
        self.MAD = None

    def singular(self, residual):
        # TODO not rechecked
        self.median = np.median(residual)

        residual[np.abs(residual - self.median) / mad(residual, c=self.c) > self.mad_pwr] = np.nan

        return residual

    def double(self, x):
        # TODO create an alternative in case the window is too short

        if (self.min_w and self.max_w is not None) and len(x) > self.max_w:  # limit the lenght of the window to the imposed maxW value
            x = x[-self.max_w:]

        m = np.median(x)  # calculate the median inside the window
        abs_dev = np.abs(x.copy() - m)  # absolute deviation

        self.left = np.median(abs_dev[x <= m])  # median of the leftish part
        self.right = np.median(abs_dev[x >= m])  # median of the rigtish part

        # TODO implement different strategy
        if self.left == 0:
            self.left = np.NaN
        if self.right == 0:
            self.right = np.NaN

        return

    def dbl_frommedian(self, x):

        x = x.astype(float)
        self.median = np.median(x)
        self.double(x)

        madMap = x.copy()
        madMap[x < self.median] = self.left
        madMap[x > self.median] = self.right

        self.MAD = np.abs(x - self.median) / madMap

        self.MAD[x == self.median] = 0

        # if self.MAD[-1] > self.sigma:
        #     return np.NaN
        # else:
        #     return x[-1]

        x[self.MAD >= self.mad_pwr] = np.NaN

        return x


def madseason(t, minW, maxW, mad_pwr):

    # seasonal decomposition by Season module
    seasons, trend = fit_seasons(t)

    # adjusted season
    adjusted = adjust_seasons(t, seasons=seasons)

    if adjusted is not None:
        # Residuals
        residual = adjusted - trend
        # rs = residual.copy()

        # Seasons
        seasons = t - adjusted

        # Trend time series
        trend = pd.Series(trend, index=adjusted.index)

        # Cleaner
        cleaned = dbl_mad_clnr(residual, minW, maxW, mad_pwr)

        # Reconstructed time series
        timeseries = trend + seasons + cleaned

        return timeseries

    else:
        return None


def dbl_mad_clnr(residual, min_w=1, max_w=108, mad_pwr=2.575):

    mado = MAD(mad_pwr)

    # TODO implement a moving window that manage propeerly the first period
    # double_mad = lambda x: mado.dbl_frommedian(x)
    # return residual.expanding(min_w).apply(double_mad, raw=True)

    return mado.dbl_frommedian(residual)


def filler(ts):

    tsC = ts.copy()

    filled = tsC.interpolate().fillna(0)

    # Find the trend inside the series
    trend = fit_trend(filled, kind='spline')

    trend = trend - trend.min()

    trendSeries = pd.Series(trend, tsC.index)

    detrended = tsC - trendSeries

    stdseason = detrended.groupby([detrended.index.month, detrended.index.day]).median()

    # Create fake yrs
    reshape = np.tile(stdseason, 3)
    reindex = np.tile(stdseason.index, 3)
    t = pd.Series(reshape, reindex)

    # TODO decide which filter it's the best one
    # Smooth by Savinsky Golet
    tSV = savgol_filter(t, 5, 2)  # TODO change to parametric

    # # Smooth by boxcar
    # tSV = t.rolling(5, win_type='bartlett', center=True).mean()  #parzen

    tsv = tSV[stdseason.count(): 2 * stdseason.count()]
    ps = pd.Series(tsv, stdseason.index)

    nanlist = tsC[tsC.isnull()]

    for index, value in nanlist.iteritems():

        nanlist.loc[index] = stdseason.loc[index.month, index.day]

    tsC.update(nanlist)

    return tsC


def fillerSeason(ts):

    tsC = ts.copy()

    filled = tsC.interpolate().fillna(0)

    # seasonal decomposition by Season module
    seasons, trend = fit_seasons(filled)

    # adjusted season
    adjusted = adjust_seasons(filled, seasons=seasons)

    if adjusted is not None:
        # Residuals
        residual = adjusted - trend
        rs = residual.copy()
        rsTS = pd.Series(rs, index=adjusted.index)

        # Seasons
        seasonsTS = filled - adjusted

        # Trend time series
        trendTS = pd.Series(trend, index=adjusted.index)

    stdseason = rsTS.groupby([rsTS.index.month, rsTS.index.day]).median()

    # Create fake yrs
    reshape = np.tile(stdseason, 3)
    reindex = np.tile(stdseason.index, 3)
    t = pd.Series(reshape, reindex)

    # TODO decide which filter it's the best one
    # Smooth by Savinsky Golet
    #tSV = savgol_filter(t, 5, 2)  # TODO change to parametric

    # # Smooth by boxcar
    tSV = t.rolling(5, win_type='bartlett', center=True).mean()  #parzen

    tsv = tSV[stdseason.count(): 2 * stdseason.count()]
    ps = pd.Series(tsv, stdseason.index)

    nanlist = tsC[tsC.isnull()]

    for index, value in nanlist.iteritems():

        nanlist.loc[index] = stdseason.loc[index.month, index.day] + \
                             trendTS.loc[index] + \
                             seasonsTS.loc[index]

    tsC.update(nanlist)

    return tsC


# def phenolo1_cleaner(s, range=[0, 250], t=5):
#
#     minVal = range[0]
#     maxVal = range[1]
#
#     for i in range(1, s.size-1):
#         vb = (s.iloc[i-1]-minVal)/(maxVal-minVal)
#         va = (s.iloc[i+1]-minVal)/(maxVal-minVal)
#         out = (s.iloc[i]-minVal)/(maxVal-minVal)
#
#         if out < vb and (out-vb)*(out-va) > t**2:
#             s[s.index[i]] = np.nan
#
#     s = s.interpolate()
#
#     return s
