# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt


def plot(pxldrl):

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

    for i in range(0, len(pxldrl.phen)):
        phency = pxldrl.phen[i]
        plt.subplot(rows, col, i+1)

        if phency.buffered is not None:
            phency.mms.plot(style='b', title='{}'.format(phency.ref_yr.values[0]))
        if phency.back is not None:
            phency.back.plot(style='-', color='gold')
        if phency.forward is not None:
            phency.forward.plot(style='-', color='olive')
        if phency.sbd is not None:
            phency.sbd.plot(style='r>',)
        if phency.sed is not None:
            phency.sed.plot(style='r<',)
        if phency.max is not None:
            phency.max.plot(style='rD')

    fig, axes = plt.subplots(7, 1, figsize=(22, 10))
    plt.tight_layout()

    plt.subplot(7, 1, 1)
    plt.tight_layout()
    if pxldrl.sbw is not None:
        pxldrl.sbw.plot(style='r', title='Season beginning week')

    plt.subplot(7, 1, 2)
    plt.tight_layout()
    if pxldrl.sew is not None:
        pxldrl.sew.plot(style='r', title='Season end week')

    plt.subplot(7, 1, 3)
    plt.tight_layout()
    if pxldrl.sl is not None:
        pxldrl.sl.plot(style='r', title='Season Lenght')

    plt.subplot(7, 1, 4)
    plt.tight_layout()
    if pxldrl.spi is not None:
        pxldrl.spi.plot(style='r', title='Season permanet')

    plt.subplot(7, 1, 5)
    plt.tight_layout()
    if pxldrl.si is not None:
        pxldrl.si.plot(style='r', title='Season Integral')

    plt.subplot(7, 1, 6)
    plt.tight_layout()
    if pxldrl.cf is not None:
        pxldrl.cf.plot(style='r', title='Cyclic fraction')

    plt.subplot(7, 1, 7)
    if pxldrl.afi is not None:
        pxldrl.afi.plot(style='r', title='Active fraction')

    plt.show(block=True)



