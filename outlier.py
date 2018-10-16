# -*- coding: utf-8 -*-
# !/usr/bin/env python

import numpy as np
import pandas as pd
from seasonal import fit_seasons, adjust_seasons, fit_trend
from statsmodels.robust import mad
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter


def cleanOLD(rs, sigma=2.575, method='polynomial', order=3):

    # Initial standard deviation
    std = rs.std()

    # Initial mean
    mean = rs.mean()

    # Initial upper and lower bounds
    LB = mean - sigma * std

    UB = mean + sigma * std

    # check
    rs[rs > UB] = np.nan
    rs[rs < LB] = np.nan

    # clean
    rs = rs.interpolate(method=method, order=order)

    return rs


def old_clean(s, range=[0, 250], t=5):

    minVal = range[0]
    maxVal = range[1]

    for i in range(1, s.size-1):
        vb = (s.iloc[i-1]-minVal)/(maxVal-minVal)
        va = (s.iloc[i+1]-minVal)/(maxVal-minVal)
        out = (s.iloc[i]-minVal)/(maxVal-minVal)

        if out < vb and (out-vb)*(out-va) > t**2:
            s[s.index[i]] = np.nan

    s = s.interpolate()

    return s


def season(t, minW, maxW, sigma):

    # seasonal decomposition by Season module
    seasons, trend = fit_seasons(t)

    # adjusted season
    adjusted = adjust_seasons(t, seasons=seasons)

    if adjusted is not None:

        # Residuals
        residual = adjusted - trend
        rs = residual.copy()

        # Seasons
        seasonsTS = t - adjusted

        # Trend time series
        trendTS = pd.Series(trend, index=adjusted.index)

        # Cleaner
        # cleanedRS = cleanOLD(residual, sigma, method='linear')
        cleanedRS = madCleaner(residual, minW, maxW, sigma)
        #cleanedRS = MAD(residual, sigma)

        # Reconstructed time series
        timeseries = trendTS + seasonsTS + cleanedRS

        return timeseries

    else:

        return None


def MAD(residual, sigma=2.575):

    median = np.median(residual)

    residual[np.abs(residual - median) / mad(residual, c=0.5) > sigma] = np.nan

    return residual


def madCleaner(residual, minW=1, maxW=108, sigma=2.575):

    doubleMad = lambda x: doubleMADsFromMedian(x, sigma, maxW)

    return residual.expanding(minW).apply(doubleMad)


def doubleMAD(x, maxW):

    if len(x) > maxW:
        x = x[-maxW:]

    m = np.median(x)
    absDev = np.abs(x - m)
    leftMAD = np.median(absDev[x <= m])
    rightMAD = np.median(absDev[x >= m])

    if leftMAD == 0:
        leftMAD = np.NaN
    if rightMAD == 0:
        rightMAD = np.NaN

    return [leftMAD, rightMAD]


def doubleMADsFromMedian(x, sigma, maxW):

    x = x.astype(float)

    twoSidedMad = doubleMAD(x, maxW)

    m = np.median(x)

    madMap = x.copy()
    madMap[x < m] = twoSidedMad[0]
    madMap[x > m] = twoSidedMad[1]

    madDistance = np.abs(x-m)/madMap

    madDistance[x == m] = 0

    if madDistance[-1] > sigma:
        return np.NaN
    else:
        return x[-1]


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
