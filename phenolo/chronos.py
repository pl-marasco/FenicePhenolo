# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np


def create(dts, dte, dektyp):
    ms = pd.date_range(dts, dte)

    if dektyp.lower() == 's5':
        dd = ms[ms.day.isin([1, 5, 10, 20, 25])]
    elif dektyp.lower() == 's10':
        dd = ms[ms.day.isin([1, 11, 21])]
    # Monthly creator
    elif dektyp.lower() == 's30':
        dd = pd.date_range(dts, dte, freq='MS')
    # Gimm model
    elif dektyp.lower() == 's15':
        dd = pd.date_range(dts, dte, freq='SMS')
    # ProbaV
    elif dektyp.lower() == 'p5':
        dd = ms[ms.day.isin([1, 6, 11, 21, 26])]
    else:
        dd = ms
    return dd


def day_calc(dekstr):
    dys_mlt = None
    dek_xyr = None

    if dekstr == 's5':
        dys_mlt = 5
        dek_xyr = 73
    if dekstr == 's10':
        dys_mlt = 10.13888888888889
        dek_xyr = 36
    elif dekstr == 's30':  # TODO to be tested
        dys_mlt = 30.41666666666667
        dek_xyr = 12
    elif dekstr == 's15':  # TODO to be tested
        dys_mlt = 24.33333333333333
        dek_xyr = 24
    elif dekstr == 'p5':  # TODO to be tested
        dys_mlt = 5.069444444444444
        dek_xyr = 72

    return dys_mlt, dek_xyr


def season_ext(ts, season_lng):
    return int(
        (ts.index.max() - ts.index.min()) / pd.Timedelta(season_lng, unit='d'))


def medspan(season_lng, param):
    if param.medspan == 0:
        medspan = season_lng / 7
    else:
        medspan = param.medspan
    return medspan


def time_resample(ts):
    ts = ts.astype(np.float64)
    return ts.asfreq('D').interpolate(method='linear').fillna(0)
