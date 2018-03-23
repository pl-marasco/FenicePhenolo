#!/usr/bin/env python

import pandas as pd
from datetime import datetime


def create(timedomain, dektyp):
    dds = timedomain.translate({ord(c): None for c in '[]'}).split('-')
    dts = datetime.strptime(dds[0], '%d/%m/%Y').date()
    dte = datetime.strptime(dds[1], '%d/%m/%Y').date()

    if dektyp.lower() == 's10':

        date_11 = pd.date_range(dts, dte, freq='SMS-11')
        date_21 = pd.date_range(dts, dte, freq='SMS-21')
        dd = date_11.union(date_21)

    # Monthly creator #TODO to be tested
    elif dektyp.lower() == 's30':
        dd = pd.date_range(dts, dte, freq='MS')

    # Gimm model #TODO to be tested
    elif dektyp.lower() == 's15':
        dd = pd.date_range(dts, dte, freq='SMS')

    # Modis dek creator #TODO to be tested
    else:
        d1 = pd.date_range(dts, dte, freq='MS')
        d15 = pd.date_range(dts, dte, freq='SMS')
        dd = d15.drop(d1)

    return dd


def day_calc(dekstr):
    dys_mlt = None
    dek_xyr = None

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
