# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from scipy.stats import median_absolute_deviation as MAD
from numba import jit

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


@jit(nopython=True, cache=True)
def mad_segments(x):
    m = np.nanmedian(x)  # calculate the median inside the window
    abs_dev = np.abs(x - m)  # absolute deviation

    lower = np.nanmedian(np.where(np.less_equal(x, m), abs_dev, np.NaN))  # median of the lower part
    higher = np.nanmedian(np.where(np.greater_equal(x, m), abs_dev, np.NaN))  # median of the higher part

    return lower, higher


@jit(nopython=True, cache=True)
def dblMAD(x, mad_pwr=2.575):

    median = np.median(x)
    if median == 0:
        return x
    lower, higher = mad_segments(x)
    madMap = np.where(x < median, lower, higher)
    MAD = np.divide(np.abs(np.subtract(x, median)), madMap)
    MAD_cld = np.where(x == median, 0, MAD)
    fx = np.where(MAD_cld >= mad_pwr, np.NaN, x)
    return fx


@jit(nopython=True, cache=True)
def median(values):
    return np.median(values)


def npMAD(x, mad_pwr):

    median = np.nanmedian(x)
    x_left = np.where(x <= median, x, np.NaN)
    x_right = np.where(x >= median, x, np.NaN)
    left_MAD = MAD(x_left, scale=1, nan_policy='omit')
    right_MAD = MAD(x_right, scale=1, nan_policy='omit')
    MAD_val = np.where(x <= median, left_MAD, right_MAD)
    MAD_dst = np.abs(x-median)/MAD_val
    MAD_dst = np.where(x == median, 0, MAD_dst)
    fx = np.where(MAD_dst >= mad_pwr, np.NaN, x)
    return fx


def doubleMAD(ts, mad_pwr=2.575):

    r = dblMAD(ts.values, mad_pwr)
    return pd.Series(r, ts.index)

#TODO rivedere doubleMAD