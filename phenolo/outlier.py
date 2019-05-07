# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

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


def mad_segments(x):
    m = np.nanmedian(x)  # calculate the median inside the window
    abs_dev = np.abs(x - m)  # absolute deviation

    lower = np.nanmedian(np.where(np.less_equal(x, m), abs_dev, np.NaN))  # median of the lower part
    higher = np.nanmedian(np.where(np.greater_equal(x, m), abs_dev, np.NaN))  # median of the higher part

    return lower, higher


def dblMAD(x, mad_pwr=2.575):
    median = np.median(x)

    lower, higher = mad_segments(x)

    madMap = np.where(x < median, lower, higher)

    MAD = np.divide(np.abs(np.subtract(x, median)), madMap)

    MAD_cld = np.where(x == median, 0, MAD)

    fx = np.where(MAD_cld >= mad_pwr, np.NaN, x)

    return fx


def doubleMAD(ts, mad_pwr=2.575):
    if ts.median() == 0:
        return ts
    else:
        r = dblMAD(ts.values, mad_pwr)
        return pd.Series(r, ts.index)
