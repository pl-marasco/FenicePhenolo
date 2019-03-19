# -*- coding: utf-8 -*-

import pandas as pd


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
    return dys_mlt, dek_xyr


def season_ext(pxldrl):
    return int((pxldrl.ts_cleaned.index.max() - pxldrl.ts_cleaned.index.min()) / pd.Timedelta(pxldrl.season_lng, unit='d'))


def medspan(season_lng, param):

    if param.medspan == 0:
        medspan = season_lng / 7
    else:
        medspan = param.medspan
    return medspan


def time_resample(ts):
    return ts.resample('D').asfreq().interpolate(method='linear').fillna(0)
