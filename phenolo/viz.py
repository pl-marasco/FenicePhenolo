# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import pandas as pd


def plot(pxldrl, param):
    fig, axes = plt.subplots(3, 1, figsize=(20, 10))

    pxldrl.ts_raw.plot(ax=axes[0], style='k', title='Original Time Series for pixel in {0}'.format(pxldrl.position))

    pxldrl.ts_resc.plot(ax=axes[1], style='y', title='TS without masked values')

    # pxldrl.ts_cleaned.plot(ax=axes[1], style='y', title='TS interplated')

    gaps = pxldrl.ts_filtered[~(pxldrl.ts_filtered.shift(-1).notnull() & pxldrl.ts_filtered.shift(1).notnull())][1:-1]

    if len(gaps.values) > 0:
        gaps.plot(ax=axes[1], style='ro', title='TS without strong outlayers')

    pxldrl.ps.plot(ax=axes[2], style='g', title='TS smoothed with Savinsky Golet and braking points')

    if len(pxldrl.pks > 0):
        pxldrl.pks.plot(ax=axes[2], style='ro')

    plt.tight_layout()

    import math
    col = 3
    rows = math.ceil(len(pxldrl.phen) / col)

    fig, axes = plt.subplots(rows, col, figsize=(22, 10))
    plt.tight_layout()

    sb_list = pd.Series()
    se_list = pd.Series()

    for i in range(0, len(pxldrl.phen)):
        phency = pxldrl.phen[i]
        plt.subplot(rows, col, i + 1)

        if phency.buffered is not None:
            phency.mms.plot(style='b', title='{}'.format(phency.ref_yr.values[0]))
        if phency.back is not None:
            phency.back.plot(style='-', color='gold')
        if phency.forward is not None:
            phency.forward.plot(style='-', color='olive')
        if phency.sbd is not None:
            sbd = pd.to_datetime(phency.sbd)
            sb = phency.mms.loc[sbd]
            sb.plot(style='r>', )
            sb_list = sb_list.append(sb)
        if phency.sed is not None:
            sed = pd.to_datetime(phency.sed)
            se = phency.mms.loc[sed]
            se.plot(style='r<', )
            se_list = se_list.append(se)
        if phency.max_idx is not None:
            phency.mms.loc[[phency.max_idx]].plot(style='rD')

    fig, axes = plt.subplots(6, 1, figsize=(22, 10))
    plt.tight_layout()

    plt.subplot(6, 1, 1)
    plt.tight_layout()
    if len(sb_list) is not None:
        sb_list.plot(style='r>', title='Season beginning')

    plt.subplot(6, 1, 2)
    plt.tight_layout()
    if len(se_list) is not None:
        se_list.plot(style='r<', title='Season end')

    plt.subplot(6, 1, 3)
    plt.tight_layout()
    if pxldrl.sl is not None:
        sl = pd.Series(pxldrl.sl, index=param.dim_unq_val)
        sl.plot.bar(style='r', title='Season Lenght')

    plt.subplot(6, 1, 4)
    plt.tight_layout()
    if pxldrl.spi is not None:
        spi = pd.Series(pxldrl.spi, index=param.dim_unq_val)
        spi.plot.bar(style='r', title='Season permanet')

    plt.subplot(6, 1, 5)
    plt.tight_layout()
    if pxldrl.si is not None:
        si = pd.Series(pxldrl.si, index=param.dim_unq_val)
        si.plot.bar(style='r', title='Season Integral')

    plt.subplot(6, 1, 6)
    plt.tight_layout()
    if pxldrl.cf is not None:
        cf = pd.Series(pxldrl.cf, index=param.dim_unq_val)
        cf.plot.bar(style='r', title='Cyclic fraction')

    # plt.subplot(7, 1, 7)
    # if pxldrl.afi is not None:
    #     afi = pd.Series(pxldrl.afi, index=param.dim_unq_val)
    #     afi.plot.bar(style='r', title='Active fraction')

    plt.show(block=True)
